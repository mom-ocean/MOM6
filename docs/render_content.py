"""Render API and manual content with Docutils and standalone BibTeX support."""
from html import escape
import contextlib
import os
from pathlib import Path
import re
import shutil
from urllib.parse import urlsplit

from docutils import nodes
from docutils.core import publish_doctree, publish_from_doctree
from docutils.parsers.rst import Parser
from docutils.parsers.rst.states import Inliner
from docutils.writers.html5_polyglot import Writer
from pybtex.database import BibliographyData, parse_file
from pybtex.plugin import find_plugin


def science_sources(docs):
    """Locate maintained science chapters alongside their implementation sources."""
    pages = {}
    for path in sorted((docs.parent / 'src').rglob('*.rst')):
        if path.stem in pages:
            raise ValueError(f'Duplicate science page {path.stem}: {pages[path.stem]} and {path}')
        pages[path.stem] = path
    return pages


def science_text(path, docs, image_prefix):
    """Adjust figure paths when staging a chapter away from its source directory."""
    prefix = os.path.relpath(docs / 'images', path.parent) + '/'
    text = path.read_text(encoding='utf-8')
    return re.sub(r'(?m)^(\s*\.\. (?:figure|image):: )' + re.escape(prefix),
                  lambda match: match[1] + image_prefix, text)


MATHJAX_URL = 'https://cdn.jsdelivr.net/npm/mathjax@3.2.2/es5/tex-svg-full.js'


GUIDE_PAGES = {
    'about': 'About this documentation',
    'equations': 'Equations',
    'discrete_space': 'Spatial Discretization',
    'discrete_time': 'Time Discretization',
    'tracers': 'Tracers in MOM6',
    'grids': 'Grids',
    'parameterizations': 'Parameterizations',
    'parameterizations_vertical': 'Vertical Parameterizations',
    'parameterizations_lateral': 'Lateral Parameterizations',
    'other_physics': 'Other Physics',
    'working_with_MOM6': 'Working with MOM6',
    'code_organization': 'Organization of the code',
    'forcing': 'Forcing',
    'parallel': 'Parallel Implementation',
    'testing': 'Testing of MOM6',
    'zzbibliography': 'Bibliography',
}
DOCS = Path(__file__).resolve().parent
SCIENCE_SOURCES = science_sources(DOCS)
SCIENCE_PAGES = {name: path.read_text(encoding='utf-8').splitlines()[2]
                 for name, path in SCIENCE_SOURCES.items()}


class GuideInliner(Inliner):
    """Capture manual references for resolution after all chapters are parsed."""
    # Docutils interpolates these patterns from vars(type(self)), rather than
    # looking up inherited attributes when initializing an Inliner subclass.
    non_whitespace_after = Inliner.non_whitespace_after
    non_whitespace_escape_before = Inliner.non_whitespace_escape_before
    non_unescaped_whitespace_escape_before = Inliner.non_unescaped_whitespace_escape_before
    simplename = Inliner.simplename
    uri_end = Inliner.uri_end
    uric = Inliner.uric
    emailc = Inliner.emailc

    def interpreted(self, rawsource, text, role, lineno):
        if role in ('cite', 'ref', 'doc'):
            node = nodes.inline(rawsource, text)
            node['manual_role'] = role
            node['manual_text'] = text
            node.line = lineno
            return [node], []
        return super().interpreted(rawsource, text, role, lineno)


def chapter_destination(name, titles):
    name = name.removeprefix('/').removesuffix('.rst')
    if name in titles:
        return name + '.html', titles[name]
    raise ValueError(f'Unknown manual chapter: {name}')


def prepare_source(source, titles):
    """Express chapter navigation and bibliography directives as Docutils nodes."""
    def toctree(match):
        lines = [line.strip() for line in match[0].splitlines()[1:] if line.strip()]
        depth = next((int(line.split(':')[-1]) for line in lines if line.startswith(':maxdepth:')), 1)
        result = f'.. container:: manual-toctree manual-depth-{depth}\n\n'
        for line in lines:
            if line.startswith(':'):
                if not line.startswith(':maxdepth:'):
                    raise ValueError(f'Unsupported chapter navigation option: {line}')
                continue
            url, title = chapter_destination(line, titles)
            result += f'   * `{title} <{url}>`__\n'
        return result + '\n'
    source = re.sub(r'(?m)^\.\. toctree::[^\n]*\n(?:(?:[ \t]+[^\n]*|)\n)*', toctree, source)
    return re.sub(r'(?m)^\.\. bibliography::[^\n]*\n(?:(?:[ \t]+[^\n]*|)\n)*',
                  '.. container:: manual-bibliography\n\n   .. Filled from the cited BibTeX entries.\n\n', source)


def bibliography(keys):
    data = BibliographyData()
    for filename in ('ocean.bib', 'references.bib', 'zotero.bib'):
        for key, entry in parse_file(DOCS / filename).entries.items():
            if key in data.entries:
                raise ValueError(f'Duplicate bibliography key: {key}')
            data.add_entry(key, entry)
    missing = keys - set(data.entries)
    if missing:
        raise ValueError(f'Unknown citation keys: {sorted(missing)}')
    style = find_plugin('pybtex.style.formatting', 'plain')()
    entries = list(style.format_bibliography(data, sorted(keys)))
    destinations, html = {}, ['<dl class="citation">']
    anchors = set()
    for entry in entries:
        target = 'cite-' + nodes.make_id(entry.key)
        if target in anchors:
            raise ValueError(f'Duplicate citation anchor: {target}')
        anchors.add(target)
        destinations[entry.key] = (target, entry.label, entry.text.render_as('text'))
        html.append(f'<dt id="{target}">[{escape(entry.label)}]</dt>'
                    f'<dd>{entry.text.render_as("html")}</dd>')
    html.append('</dl>')
    return destinations, '\n'.join(html)


def section_links(document, depth):
    result = nodes.bullet_list()
    if depth <= 0:
        return result
    for section in document.children:
        if not isinstance(section, nodes.section):
            continue
        item = nodes.list_item()
        title = section.next_node(nodes.title).astext()
        item += nodes.paragraph('', '', nodes.reference('', title, refuri='#' + section['ids'][0]))
        children = section_links(section, depth - 1)
        if children.children:
            item += children
        result += item
    return result


def render_guides(rst_dir, warning_stream, html_dir=None, *, api_reader=None):
    """Resolve chapters and optional API pages together before writing HTML.

    Supplying an API reader includes the generated pages in the same reference
    and citation index. The reader adds presentation transforms such as graphs.
    """
    titles = {name: title for name, title in {**GUIDE_PAGES, **SCIENCE_PAGES}.items()
               if (rst_dir / (name + '.rst')).is_file()}
    manual_pages = set(titles)
    if api_reader is not None:
        for path in rst_dir.glob('*.rst'):
            if path.stem not in titles:
                titles[path.stem] = path.stem
    if not titles:
        return {}
    documents, targets, keys = {}, {}, set()
    settings = {'warning_stream': warning_stream, 'halt_level': 2,
                'stylesheet_path': [], 'embed_stylesheet': False,
                'math_output': 'MathJax ' + MATHJAX_URL,
                'output_encoding': 'unicode', 'initial_header_level': 2,
                'syntax_highlight': 'short'}
    for name in titles:
        path = rst_dir / (name + '.rst')
        source = prepare_source(path.read_text(encoding='utf-8'), titles)
        source_path = str(path.resolve())
        with contextlib.chdir(rst_dir):
            document = publish_doctree(source, source_path=source_path,
                                      reader=api_reader() if api_reader and name not in manual_pages else None,
                                      parser=Parser(inliner=GuideInliner()), settings_overrides=settings)
        documents[name] = document
        for label, target in document.nameids.items():
            if not document.nametypes.get(label) or target is None:
                continue
            if label in targets:
                raise ValueError(f'Duplicate manual target: {label}')
            element = document.ids[target]
            heading = next(iter(element.findall(nodes.title)), None)
            title = (heading.astext() if heading is not None else
                     element.astext() if isinstance(element, nodes.rubric) else titles[name])
            targets[label] = (name + '.html#' + target, title)
        for node in document.findall(nodes.inline):
            if node.get('manual_role') == 'cite':
                keys.update(key.strip() for key in node['manual_text'].split(','))
    citations, bibliography_html = bibliography(keys)
    if keys and 'zzbibliography' not in documents:
        raise ValueError('Citations require the bibliography chapter')
    bodies = {}
    for name, document in documents.items():
        for node in list(document.findall(nodes.inline)):
            role = node.get('manual_role')
            if role is None:
                continue
            text = node['manual_text']
            if role == 'cite':
                links = []
                for key in text.split(','):
                    target, label, tooltip = citations[key.strip()]
                    if links:
                        links.append(nodes.Text(', '))
                    links.append(nodes.reference('', '[' + label + ']',
                                                 refuri='zzbibliography.html#' + target,
                                                 reftitle=tooltip, classes=['bibtex-citation']))
                node.replace_self(links)
                continue
            explicit = re.fullmatch(r'(.+?)\s*<([^>]+)>', text)
            label, target = explicit.groups() if explicit else (None, text)
            if role == 'doc':
                url, title = chapter_destination(target, titles)
            else:
                destination = targets.get(nodes.fully_normalize_name(target))
                if destination is None:
                    raise ValueError(f'{name}:{node.line}: Unknown manual reference {target}')
                url, title = destination
            node.replace_self(nodes.reference('', label or title, refuri=url))
        for container in list(document.findall(nodes.container)):
            if 'manual-bibliography' in container['classes']:
                container.replace_self(nodes.raw('', bibliography_html, format='html'))
            elif 'manual-toctree' in container['classes']:
                depth = int(next(c.removeprefix('manual-depth-') for c in container['classes']
                                 if c.startswith('manual-depth-')))
                for link in list(container.findall(nodes.reference)):
                    page = link['refuri'].removesuffix('.html')
                    if page in documents and depth > 1:
                        outline = section_links(documents[page], depth - 1)
                        for child in outline.findall(nodes.reference):
                            child['refuri'] = page + '.html' + child['refuri']
                        if outline.children:
                            link.parent.parent += outline
        for image in document.findall(nodes.image):
            uri = image['uri']
            if urlsplit(uri).scheme or urlsplit(uri).netloc:
                image['loading'] = 'link'
                continue
            if name not in manual_pages and uri.endswith('.svg'):
                continue
            if name not in manual_pages and not uri.startswith('images/'):
                uri = 'images/' + uri
            if not uri.startswith('images/'):
                raise ValueError(f'{name}: Unsupported manual image path: {uri}')
            original = rst_dir / uri
            if not original.is_file():
                raise FileNotFoundError(original)
            image['uri'] = '_images/' + uri.removeprefix('images/')
            image['loading'] = 'link'
            if html_dir is not None:
                destination = html_dir / image['uri']
                destination.parent.mkdir(parents=True, exist_ok=True)
                shutil.copyfile(original, destination)
        # Use the ordinary HTML5 writer without global role/directive registrations.
        writer = Writer()
        with contextlib.chdir(rst_dir):
            publish_from_doctree(document, writer=writer,
                                 settings_overrides={**settings, 'image_loading': 'embed'
                                                     if name not in manual_pages else 'link'})
        bodies[name] = writer.parts['html_body']
    return bodies
