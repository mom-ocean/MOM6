"""Deployment checks use rendered files rather than trusting generated URLs."""
import contextlib
import io
import re
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from build_docs import build as build_site, check_links, generate_reference
from render_site import render_site
from render_content import GUIDE_PAGES, SCIENCE_PAGES, render_guides


MANUAL_API_FIXTURES = {
    'mom_hor_visc': 'section_horizontal_viscosity',
    'mom_meke': 'section_MEKE',
    'mom_mixed_layer_restrat': 'section_mle',
    'mom_interface_filter': 'section_interface_filter',
    'mom_lateral_mixing_coeffs': 'section_Resolution_Function',
    'mom_tidal_forcing': 'section_tides',
    'mom_load_love_numbers': 'section_Love_numbers',
    'mom_self_attr_load': 'section_SAL',
    'mom_spherical_harmonics': 'section_spherical_harmonics',
    'mom_unit_scaling': None,
}


def build(*args, **kwargs):
    """Supply the manual's API destinations for tiny synthetic Fortran projects.

    These tests exercise the real manual as well as synthetic API pages. The
    manual's referenced science sections must be local, even in test fixtures.
    Full-project validation uses build_docs.build directly.
    """
    def generate_with_manual_fixtures(source, rst_dir, reports, **options):
        summary = generate_reference(source, rst_dir, reports, **options)
        for page, label in MANUAL_API_FIXTURES.items():
            path = rst_dir / (page + '.rst')
            assert not path.exists(), f'Fixture would overwrite generated page: {page}'
            target = f'.. _{label}:\n\n' if label else ''
            path.write_text(target + page + '\n' + '=' * len(page) + '\n\n'
                            'Local API destination for the manual fixture.\n')
        return summary
    with patch('build_docs.generate_reference', generate_with_manual_fixtures):
        return build_site(*args, **kwargs)


class ApiDocumentationTests(unittest.TestCase):
    def test_link_validation_rejects_duplicate_ids_within_each_page(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            first = root / 'first.html'
            first.write_text('<h2 id="example">First</h2><a href="#example">Link</a>')
            (root / 'second.html').write_text('<h2 id="example">Second</h2>')
            self.assertEqual(check_links(root), (1, 0))
            first.write_text(first.read_text() + '<h2 id="example">Duplicate</h2>')
            with self.assertRaisesRegex(ValueError, "first.html: duplicate HTML id 'example'"):
                check_links(root)
            # Duplicate IDs are invalid even without any links targeting them.
            first.write_text('<div id="unused"></div><span id="unused"></span>')
            with self.assertRaisesRegex(ValueError, "first.html: duplicate HTML id 'unused'"):
                check_links(root)

    def test_science_links_resolve_direct_section_labels(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / 'first.rst').write_text(
                'First\n=====\n\n:ref:`section_MEKE`\n\n'
                ':doc:`units <mom_unit_scaling>`\n')
            (root / 'mom_meke.rst').write_text(
                'MEKE\n====\n\n.. _section_MEKE:\n\nScience\n-------\n\nLocal science.\n')
            (root / 'mom_unit_scaling.rst').write_text('Units\n=====\n\nLocal units.\n')
            titles = {'first': 'First', 'mom_meke': 'MEKE', 'mom_unit_scaling': 'Units'}
            with patch('render_content.GUIDE_PAGES', titles), patch('render_content.SCIENCE_PAGES', {}):
                pages = render_guides(root, io.StringIO())
                self.assertIn('href="mom_meke.html#section-meke"', pages['first'])
                self.assertIn('href="mom_unit_scaling.html"', pages['first'])
                self.assertNotIn('readthedocs', pages['first'])
                # References follow the label when its owning page changes.
                (root / 'mom_meke.rst').rename(root / 'other.rst')
                titles['other'] = titles.pop('mom_meke')
                pages = render_guides(root, io.StringIO())
                self.assertIn('href="other.html#section-meke"', pages['first'])
                (root / 'first.rst').write_text(
                    'First\n=====\n\n:doc:`missing`\n')
                with self.assertRaisesRegex(ValueError, 'Unknown manual chapter'):
                    render_guides(root, io.StringIO())

    def test_api_comment_references_share_manual_targets_and_bibliography(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / 'source'
            source.mkdir()
            (source / 'demo.F90').write_text(
                'module demo\ncontains\n!> Follow \\cite Adcroft2019 and \\ref science.\n'
                'subroutine run()\nend subroutine run\n'
                '!> \\namespace demo\n!! \\section science Scientific details\n'
                '!!\n!! \\f[ x = y^2 \\f]\nend module demo\n')
            with contextlib.redirect_stdout(io.StringIO()):
                build([source], root / 'site', graph_depth=0)
            page = (root / 'site/html/demo.html').read_text()
            self.assertIn('href="demo.html#science"', page)
            self.assertIn('href="zzbibliography.html#cite-adcroft2019"', page)
            self.assertIn('src="https://cdn.jsdelivr.net/npm/mathjax@3.2.2/es5/tex-svg-full.js"', page)
            self.assertIn('Scientific details', page)
            self.assertNotIn('\\namespace', page)

    def test_manual_references_citations_and_code_are_resolved_without_sphinx(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / 'first.rst').write_text(
                'First\n=====\n\n.. toctree::\n   :maxdepth: 2\n\n   second\n\n'
                'See :ref:`the result <result>` and :cite:`Adcroft2019,reichl2018`.\n\n'
                '.. code-block:: text\n\n   :cite:`not-a-citation`\n')
            (root / 'second.rst').write_text(
                'Second\n======\n\nIntroduction.\n\n.. _result:\n\n'
                'Result\n------\n\nA result.\n')
            (root / 'zzbibliography.rst').write_text(
                'Bibliography\n============\n\n.. bibliography:: ocean.bib\n   :style: plain\n')
            titles = {'first': 'First', 'second': 'Second', 'zzbibliography': 'Bibliography'}
            with patch('render_content.GUIDE_PAGES', titles), patch('render_content.SCIENCE_PAGES', {}):
                warnings = io.StringIO()
                pages = render_guides(root, warnings)
                self.assertEqual(warnings.getvalue(), '')
                self.assertIn('href="second.html#result"', pages['first'])
                self.assertIn('the result</a>', pages['first'])
                self.assertIn(':cite:`not-a-citation`', pages['first'])
                self.assertIn('href="zzbibliography.html#cite-adcroft2019"', pages['first'])
                self.assertIn('id="cite-adcroft2019"', pages['zzbibliography'])
                for name, body in pages.items():
                    (root / (name + '.html')).write_text(body)
                check_links(root)
                (root / 'first.rst').write_text('First\n=====\n\n:ref:`missing-target`\n')
                with self.assertRaisesRegex(ValueError, 'Unknown manual reference missing-target'):
                    render_guides(root, io.StringIO())
                (root / 'first.rst').write_text('First\n=====\n\n:cite:`missing-paper`\n')
                with self.assertRaisesRegex(ValueError, 'Unknown citation keys'):
                    render_guides(root, io.StringIO())

    def test_renderer_uses_title_case_without_changing_urls(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            rst = root / 'rst'
            rst.mkdir()
            title = 'Ocean_API module reference'
            (rst / 'ocean_api.rst').write_text(
                '=' * len(title) + '\n' + title + '\n' + '=' * len(title) + '\n\n'
                '.. container:: api-procedure\n\n'
                '   .. _flint-ocean-0api-1advance:\n\n'
                '   subroutine Advance(State, dt)\n'
                '     Update the state.\n')
            warnings = io.StringIO()
            render_site(rst, root / 'html', warnings)
            page = (root / 'html/ocean_api.html').read_text()
            self.assertEqual(warnings.getvalue(), '')
            self.assertIn('href="ocean_api.html" aria-current="page">Ocean_API</a>', page)
            self.assertIn('<span>Ocean_API</span>', page)
            self.assertIn('api-procedure', page)
            self.assertIn('subroutine Advance(State, dt)', page)
            check_links_directory = root / 'html'
            (check_links_directory / 'index.html').write_text('<p>Index</p>')
            (check_links_directory / 'api-index.html').write_text('<p>API reference</p>')
            check_links(check_links_directory)

    def test_dependency_sources_are_excluded_by_resolved_ownership(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / 'src'
            source.mkdir()
            gsw = root / 'pkg/GSW-Fortran'
            cvmix = root / 'pkg/CVMix-src'
            gsw.mkdir(parents=True)
            cvmix.mkdir()
            # Deliberately invalid: these implementation files must not be parsed.
            (gsw / 'library.F90').write_text('not Fortran\n')
            (cvmix / 'library.F90').write_text('not Fortran\n')
            (source / 'unrelated_name.F90').symlink_to(gsw / 'library.F90')
            (source / 'client.F90').write_text(
                'module client\nuse gsw_mod_toolbox, only: gsw_rho\n'
                'use cvmix_kpp, only: cvmix_init_kpp\ncontains\n'
                'subroutine run()\ncall gsw_rho()\ncall cvmix_init_kpp()\n'
                'end subroutine run\nend module client\n')
            output = root / 'site'
            with patch('build_docs.REPOSITORY', root), contextlib.redirect_stdout(io.StringIO()):
                summary = build([source, cvmix], output)
            self.assertEqual(summary['source_files'], 1)
            self.assertEqual(summary['module_pages'], 1)
            self.assertEqual(summary['unresolved_references'], 0)
            self.assertEqual(summary['external_libraries'], {})
            self.assertEqual(summary['excluded_dependency_files'], {'GSW': 1, 'CVMix': 1})
            self.assertEqual(summary['external_references'], 2)
            self.assertIn('gsw_mod_toolbox::gsw_rho',
                          (output / 'reports/external.txt').read_text())
            self.assertEqual({path.name for path in (output / 'rst').glob('*.rst')},
                             {'client.rst', 'index.rst', 'api-index.rst'}
                              | {name + '.rst' for name in {**GUIDE_PAGES, **SCIENCE_PAGES, **MANUAL_API_FIXTURES}})
            self.assertIn('unrelated_name.F90', (output / 'reports/excluded-dependencies.txt').read_text())
            (gsw / 'library.F90').unlink()
            with patch('build_docs.REPOSITORY', root), contextlib.redirect_stdout(io.StringIO()):
                missing_summary = build([source, cvmix], root / 'without-submodule')
            self.assertEqual(missing_summary, summary)

    def test_missing_submodule_target_has_actionable_error(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / 'src'
            source.mkdir()
            (source / 'gsw_gibbs.f90').symlink_to('../pkg/GSW-Fortran/gsw_gibbs.f90')
            with self.assertRaisesRegex(FileNotFoundError, 'Source symlink has a missing target'):
                build(source, root / 'site')

    def test_cross_module_build_and_reports(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / 'src'
            source.mkdir()
            infra = root / 'config_src' / 'infra' / 'FMS2'
            infra.mkdir(parents=True)
            (infra / 'library.F90').write_text(
                'module library\nprivate :: secret\ncontains\nsubroutine helper()\n'
                'call secret()\nend subroutine helper\n'
                'subroutine secret()\nend subroutine secret\nend module library\n')
            (source / 'client.F90').write_text(
                'module client\nuse library, only: work => helper\n'
                'use netcdf, only: nf90_close\ncontains\n'
                'subroutine run()\ncall work()\ncall nf90_close()\ncall unavailable()\n'
                'end subroutine run\nend module client\n')
            output = root / 'site'
            with contextlib.redirect_stdout(io.StringIO()):
                summary = build([source, infra], output)
            self.assertEqual(summary['source_files'], 2)
            self.assertEqual(summary['module_pages'], 2)
            self.assertEqual(summary['html_pages'],
                             4 + len(GUIDE_PAGES) + len(SCIENCE_PAGES) + len(MANUAL_API_FIXTURES))
            self.assertEqual(summary['call_graphs'], 4)
            self.assertEqual(summary['unresolved_references'], 1)
            self.assertEqual(summary['external_references'], 1)
            self.assertEqual(summary['external_libraries'], {})
            self.assertGreater(summary['cross_file_links'], 0)
            self.assertIn('library.html#', (output / 'html/client.html').read_text())
            self.assertIn('unavailable', (output / 'reports/unresolved.txt').read_text())
            self.assertEqual((output / 'reports/render.log').read_text(), '')
            self.assertTrue((output / 'rst/index.rst').exists())
            home = (output / 'html/index.html').read_text()
            self.assertIn('MOM6 documentation', home)
            self.assertIn('href="api-index.html"', home)
            self.assertIn('href="code_organization.html"', home)
            self.assertIn('href="about.html"', home)
            about = (output / 'html/about.html').read_text()
            self.assertIn('<title>About this documentation — MOM6</title>', about)
            self.assertIn('href="about.html" aria-current="page"', about)
            self.assertIn('Flint', about)
            self.assertIn(f'<span class="count">{2 + len(MANUAL_API_FIXTURES)}</span>', about)
            chapter = (output / 'html/code_organization.html').read_text()
            self.assertIn('<title>Organization of the code — MOM6</title>', chapter)
            self.assertIn('href="code_organization.html" aria-current="page"', chapter)
            self.assertIn('configuration directories</a>', chapter)
            self.assertIn('FMS_cap', chapter)
            self.assertNotIn(':ref:', chapter)
            self.assertIn(f'<span class="count">{2 + len(MANUAL_API_FIXTURES)}</span>', chapter)
            self.assertIn('href="equations.html"', home)
            equations = (output / 'html/equations.html').read_text()
            self.assertIn('href="Notation.html"', equations)
            governing = (output / 'html/Governing_Equations.html').read_text()
            self.assertIn('src="https://cdn.jsdelivr.net/npm/mathjax@3.2.2/es5/tex-svg-full.js"', governing)
            self.assertIn('class="math"', governing)
            indexing = (output / 'html/Horizontal_Indexing.html').read_text()
            self.assertIn('_images/Arakawa_C_grid.png', indexing)
            self.assertIn('<ul', indexing)
            self.assertNotIn('HI%%', indexing)
            self.assertIn('href="zzbibliography.html#', equations)
            working = (output / 'html/working_with_MOM6.html').read_text()
            self.assertIn('href="code_organization.html"', working)
            lateral = (output / 'html/parameterizations_lateral.html').read_text()
            self.assertNotIn(':cite:', lateral)
            self.assertNotIn(':ref:', lateral)
            self.assertIn('mom_hor_visc.html#section-horizontal-viscosity', lateral)
            self.assertNotIn('https://mom6.readthedocs.io/en/main/api/', lateral)
            module_index = (output / 'html/api-index.html').read_text()
            self.assertIn('href="client.html"', module_index)
            self.assertNotIn('aria-current="page">api-index</a>', module_index)
            page = (output / 'html/client.html').read_text()
            library = (output / 'html/library.html').read_text()
            self.assertIn('External calls', page)
            self.assertNotIn('<dt>Calls into</dt>', page)
            self.assertNotIn('<dt>Called by</dt>', library)
            text_views = re.findall(
                r'<details class="call-graph-text"><summary>Text view</summary>(.*?)</details>',
                library, flags=re.DOTALL)
            self.assertEqual(len(text_views), 3)
            self.assertTrue(any('library::secret</a>' in view for view in text_views))
            self.assertTrue(any('client::run</a>' in view for view in text_views))
            self.assertIn('netcdf', page)
            self.assertNotIn('unavailable', page)
            self.assertIn('netcdf::nf90_close', (output / 'reports/external.txt').read_text())
            self.assertIn('aria-label="Reference navigation"', page)
            self.assertIn('aria-current="page">client</a>', page)
            self.assertIn('assets/api.css', page)
            self.assertTrue((output / 'html/assets/api.css').is_file())
            self.assertTrue((output / 'html/assets/api.js').is_file())
            self.assertIn('<svg', page)
            self.assertIn('class="call-graph" tabindex="0"', page)
            self.assertIn('<details class="call-graph-panel"><summary>Callees call graph</summary>', page)
            self.assertIn('font-size="12.00"', page)
            self.assertIn('Callees', page)
            self.assertIn('Callers', library)
            self.assertIn('<details class="call-graph-panel"><summary>Callers call graph</summary>', library)
            self.assertEqual(library.count('<details class="call-graph-panel">'), library.count('<svg'))
            self.assertIn('[private/internal]', library)
            ids = re.findall(r'\bid="([^"]+)"', library)
            self.assertEqual(len(ids), len(set(ids)))
            self.assertEqual(len(list((output / 'rst').glob('*.dot'))), 4)
            self.assertEqual(len(list((output / 'rst').glob('*.svg'))), 4)
            with contextlib.redirect_stdout(io.StringIO()):
                graphless = build([source, infra], root / 'without-graphs', graph_depth=0)
            self.assertEqual(graphless['call_graphs'], 0)
            self.assertIn('<dt>Calls into</dt>', (root / 'without-graphs/html/client.html').read_text())

    def test_alternative_modules_cannot_overwrite_pages(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            sources = [root / 'solo', root / 'ice_solo']
            for source in sources:
                source.mkdir()
                (source / 'fluxes.F90').write_text('module fluxes\nend module fluxes\n')
            with self.assertRaisesRegex(ValueError, 'Duplicate module fluxes'):
                build(sources, root / 'site')

    def test_broken_page_and_fragment_links_fail(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / 'index.html').write_text('<a href="missing.html">Missing page</a>')
            with self.assertRaisesRegex(ValueError, 'missing.html'):
                check_links(root)
            (root / 'missing.html').write_text('<p id="present">Text</p>')
            (root / 'index.html').write_text('<a href="missing.html#absent">Missing anchor</a>')
            with self.assertRaisesRegex(ValueError, 'missing.html#absent'):
                check_links(root)
            (root / 'index.html').write_text(
                '<svg><a xlink:href="missing.html#absent">Graph node</a></svg>')
            with self.assertRaisesRegex(ValueError, 'missing.html#absent'):
                check_links(root)


if __name__ == '__main__':
    unittest.main()
