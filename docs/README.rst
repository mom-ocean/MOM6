Building MOM6 documentation
===========================

The documentation site combines a Flint-generated Fortran API reference with
the user and science guide. Docutils renders reStructuredText, Pybtex formats
citations, MathJax renders equations, and Graphviz draws call diagrams.
MathJax loads from a pinned CDN URL when viewing pages containing equations,
so equation rendering requires internet access.
Doxygen, Sphinx, Pandoc, and TeX are not build dependencies.

Quick start
-----------

Install Python 3.13 or newer and Graphviz (the ``dot`` executable). From the
repository root::

    python3 -m venv .venv-api
    .venv-api/bin/python -m pip install -r docs/requirements-api.txt
    .venv-api/bin/python -m unittest discover -s docs -p 'test_docs.py'
    .venv-api/bin/python docs/build_docs.py --output docs/_build/flint-api

Open ``docs/_build/flint-api/html/index.html``. The output directory must be
new; choose a different ``--output`` directory for subsequent builds.

For a faster content and link check, add ``--graph-depth 0``. To rerender
existing generated reST after a theme change::

    .venv-api/bin/python docs/render_site.py \
        docs/_build/flint-api/rst docs/_build/flint-api/html

The Makefile provides ``html`` and ``check`` shortcuts. For example::

    make -C docs check PYTHON="$PWD/.venv-api/bin/python"
    make -C docs html PYTHON="$PWD/.venv-api/bin/python" \
        OUTPUT=docs/_build/preview GRAPH_DEPTH=0

``OUTPUT`` is relative to the repository root unless it is absolute.

Authoring
---------

* Keep the existing Doxygen-style Fortran comments. Flint converts supported
  markup directly to reST; the Doxygen program is not involved. Unsupported
  commands fail generation with the owning source file in the error.
* Edit the user-guide chapters in ``docs/*.rst`` and the 44 detailed scientific
  chapters in ``src/**/*.rst``. Images live in ``docs/images``; bibliography
  entries live in ``ocean.bib``, ``references.bib``, and ``zotero.bib``.
* Edit ``api-home.rst`` for the landing page and ``api_theme`` for presentation.
* Use explicit targets, ``:ref:``, ``:cite:``, and standard reST mathematics.
  Literal blocks preserve ASCII diagrams and their spacing.

See ``api-deployment.rst`` for source selection, dependencies, rendering,
validation, and deployment.
``documentation-checks.rst`` describes the separate API documentation-presence
check run by source verification.

For running MOM6 within CESM, see the `NCAR quick start guide
<https://github.com/ESCOMP/MOM_interface/wiki/Quick-Start>`_.

Hosting and history
-------------------

GitHub Actions builds the site and publishes eligible branches through Pages.
Read the Docs runs the same builder through the custom jobs in
``.readthedocs.yml`` and publishes its HTML. Both hosts use the same pinned
dependencies.

The previous documentation pipeline was developed with contributions from
Angus Gibson, Rob Cermak, and Marshall Ward. Its implementation and build
instructions remain available in Git history.
