#!/usr/bin/env python3
"""Build the experimental Flint API reference, then validate local HTML links."""
import argparse
from html.parser import HTMLParser
import json
import shutil
from pathlib import Path
from urllib.parse import unquote, urlsplit

from flint.tools.gendoc import generate_docs
from render_site import GUIDE_PAGES, render_site
from render_content import SCIENCE_SOURCES, science_text


DEFAULT_SOURCES = (
    'src',
    'config_src/drivers/solo_driver',
    'config_src/external',
    'config_src/infra/FMS2',
    'config_src/memory/dynamic_symmetric',
)
REPOSITORY = Path(__file__).resolve().parents[1]
DEPENDENCY_SOURCES = {'GSW': 'pkg/GSW-Fortran', 'CVMix': 'pkg/CVMix-src'}


class HtmlLinks(HTMLParser):
    def __init__(self):
        super().__init__()
        self.ids = set()
        self.duplicate_ids = set()
        self.links = []

    def handle_starttag(self, tag, attrs):
        attrs = dict(attrs)
        if 'id' in attrs:
            if attrs['id'] in self.ids:
                self.duplicate_ids.add(attrs['id'])
            self.ids.add(attrs['id'])
        if tag == 'a':
            href = attrs.get('href', attrs.get('xlink:href'))
            if href is not None:
                self.links.append(href)
        if tag in ('img', 'script') and attrs.get('src'):
            self.links.append(attrs['src'])
        if tag == 'link' and attrs.get('href'):
            self.links.append(attrs['href'])


def check_links(directory):
    """Check relative page links and fragments in the rendered site."""
    pages = {}
    for path in sorted(directory.glob('*.html')):
        page = HtmlLinks()
        page.feed(path.read_text(encoding='utf-8'))
        pages[path.name] = page
    count = cross_file = 0
    errors = []
    for filename, page in pages.items():
        for identifier in sorted(page.duplicate_ids):
            errors.append(f'{filename}: duplicate HTML id {identifier!r}')
        for href in page.links:
            url = urlsplit(href)
            if url.scheme or url.netloc:
                continue
            target = unquote(url.path) or filename
            fragment = unquote(url.fragment)
            count += 1
            cross_file += target != filename
            if (not (directory / target).is_file()
                    or (fragment and (target not in pages or fragment not in pages[target].ids))):
                errors.append(f'{filename}: {href}')
    if errors:
        raise ValueError('HTML validation errors:\n' + '\n'.join(errors))

    return count, cross_file


def generate_reference(source, rst_dir, reports, *, graph_depth=1):
    """Generate the API with Flint, then stage MOM6's maintained chapters."""
    excluded_lib_roots = {library: REPOSITORY / path
                          for library, path in DEPENDENCY_SOURCES.items()}
    sources = [source] if isinstance(source, Path) else list(source)

    result = generate_docs(
        sources, rst_dir, excludes=excluded_lib_roots, link_format='html',
        call_graphs='svg', graph_depth=graph_depth, markup='doxygen',
        include_private=True, graph_fontsize=12, reports=reports,
        index_file='api-index.rst', index_title='MOM6 API reference')

    shutil.copyfile(Path(__file__).with_name('api-home.rst'), rst_dir / 'index.rst')

    for name in GUIDE_PAGES:
        shutil.copyfile(Path(__file__).with_name(name + '.rst'), rst_dir / (name + '.rst'))

    for name, path in SCIENCE_SOURCES.items():
        text = science_text(path, Path(__file__).parent, 'images/')
        (rst_dir / (name + '.rst')).write_text(text, encoding='utf-8')

    shutil.copytree(Path(__file__).parent / 'images', rst_dir / 'images')

    return result.summary


def build(source, output, *, graph_depth=1):
    if graph_depth < 0:
        raise ValueError('Graph depth must be nonnegative')

    # Generate output directories
    output.mkdir(parents=True, exist_ok=False)
    rst_dir, html_dir, reports = (
        output / name for name in ('rst', 'html', 'reports')
    )
    for directory in (rst_dir, html_dir, reports):
        directory.mkdir()

    # Generate the reST and call graph figures
    summary = generate_reference(source, rst_dir, reports,
                                 graph_depth=graph_depth)

    # Generate the HTML from the reST output
    with (reports / 'render.log').open('w', encoding='utf-8') as warnings:
        render_site(rst_dir, html_dir, warnings)

    # Validate HTML links (broken anchors, missing pages, etc)
    links, cross_file = check_links(html_dir)

    # Update the documentation statistics for the workflow summary
    summary.update(html_pages=len(list(html_dir.glob('*.html'))),
                   checked_links=links, cross_file_links=cross_file)

    # Disable GitHub Page's Jekyll build process with a sentinel file.
    (html_dir / '.nojekyll').touch()

    # Print the documentation process in the workflow summary.
    (reports / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n',
                                          encoding='utf-8')
    print(json.dumps(summary, indent=2))

    # Return the summary data to any callers
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source',
        type=Path,
        action='append',
        help=('Source directory; repeat to combine roots '
              '(replaces the MOM6 defaults)')
    )
    parser.add_argument('--output',
        type=Path,
        default=Path('docs/_build/flint-api')
    )
    parser.add_argument('--graph-depth',
        type=int,
        default=1,
        help='Caller/callee graph depth (default: 1)'
    )
    args = parser.parse_args()

    build(args.source or [Path(root) for root in DEFAULT_SOURCES], args.output,
          graph_depth=args.graph_depth)


if __name__ == '__main__':
    main()
