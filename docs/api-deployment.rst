Documentation build and deployment
==================================

See ``README.rst`` for dependencies, local build commands, and authoring
instructions, and ``documentation-checks.rst`` for the source-documentation
checking policy.

Source selection
----------------

The default API build scans ``src`` and these configuration directories:

* ``config_src/drivers/solo_driver``
* ``config_src/external`` (all subdirectories)
* ``config_src/infra/FMS2``
* ``config_src/memory/dynamic_symmetric`` (also supplies include headers)

Repeat ``--source DIRECTORY`` to replace these defaults. Alternative drivers
can define the same module names; duplicate names fail the build. Partial
selections must include modules referenced by the maintained manual.

GSW and CVMix implementations are excluded, including symlinks into their
dependency trees. MOM6's wrappers remain documented. Exclusions are configured
in ``DEPENDENCY_SOURCES`` in ``build_docs.py``.
The documentation build does not require submodule checkouts. Other broken
source symlinks fail the build.

Rendering and validation
------------------------

``build_docs.py`` calls Flint's ``generate_docs()`` to discover and parse sources,
generate module reST, the API index, and call graphs, and write diagnostic reports. It stages
the maintained manual and invokes the HTML renderer. ``render_content.py`` lists the guide
pages and discovers science chapters in ``src/**/*.rst``.
The template, stylesheet, and navigation script live in ``api_theme``.

Named references and citations resolve across API and manual pages. Pybtex
formats the existing BibTeX bibliographies. Unsupported comment markup,
renderer warnings, missing targets or citations, duplicate anchors, and broken
HTML links fail the build. Scientific/API references have no automatic
external-site fallback.

Pages containing equations load MathJax from a version-pinned jsDelivr URL,
defined by ``MATHJAX_URL`` in ``render_content.py``. MathJax requires no
build-time installation; viewing rendered equations requires internet access.

Graphviz generates caller/callee diagrams with one hop by default. Use
``--graph-depth 2`` to expand them or ``--graph-depth 0`` to disable them.
Diagrams show possible calls, not execution order. External and unresolved
calls are omitted from diagrams; explicit imports identify external references
without guessed documentation URLs.

Build output contains ``rst``, ``html``, and ``reports`` directories. Reports
include ``summary.json``, ``excluded-dependencies.txt``, ``external.txt``, and
``unresolved.txt``. Unresolved calls remain in diagnostics rather than published
call lists and should be investigated.

Update ``requirements-api.txt`` explicitly when adopting Flint changes, and
validate using the installed pin rather than only a local Flint checkout.

GitHub Actions and Pages
------------------------

``.github/workflows/docs-and-style.yml`` checks source whitespace and
documentation, runs builder tests, and builds and validates the complete site.
Linux verification handles compilation and numerical tests separately.

The workflow retains source-check diagnostics as ``docstring-check-reports``
and site output as ``flint-api-reference``. Download the latter to preview PR
builds locally. Required-check configuration should use
``Check source style and documentation`` for the source checks.

The deployment job's condition selects the repository and branches eligible
for Pages publication; pull requests do not deploy. To enable publication,
select **GitHub Actions** in **Settings > Pages** and allow the selected branches
in the ``github-pages`` environment. Update the workflow's repository and branch
condition when adopting it in another repository. Generated files are uploaded
as artifacts rather than committed to a branch.

Read the Docs
-------------

``.readthedocs.yml`` installs the same pinned dependencies and Graphviz, runs
``docs/build_docs.py``, and copies HTML to ``$READTHEDOCS_OUTPUT/html``.
GitHub Pages and Read the Docs therefore publish the same site structure.
Version activation and hosting settings are managed on Read the Docs.
Redirects from old ``api/generated/...`` URLs require separate configuration.
