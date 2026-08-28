import os.path
from lxml import etree as ET
from sphinx.errors import ExtensionError


def set_doxygen_xml(app):
    """Load all doxygen XML files from the app config variable
    `app.config.doxygen_xml` which should be a path to a directory
    containing doxygen xml output. If the configured path is relative,
    it is resolved against `app.confdir` rather than the current working
    directory -- Sphinx may have any cwd by the time builder-inited fires,
    and in particular RTD runs sphinx-build from the repo root.
    """
    doxygen_xml = app.config.doxygen_xml
    if not os.path.isabs(doxygen_xml):
        doxygen_xml = os.path.join(app.confdir, doxygen_xml)

    err = ExtensionError(
        '[autodoc_doxygen] No doxygen '
        'xml output found in doxygen_xml="%s"' % doxygen_xml)

    if not os.path.isdir(doxygen_xml):
        raise err

    files = [os.path.join(doxygen_xml, f)
             for f in os.listdir(doxygen_xml)
             if f.lower().endswith('.xml') and not f.startswith('._')]
    if len(files) == 0:
        raise err

    setup.DOXYGEN_ROOT = ET.ElementTree(ET.Element('root')).getroot()
    for file in files:
        root = ET.parse(file).getroot()
        for node in root:
            setup.DOXYGEN_ROOT.append(node)


def get_doxygen_root():
    """Get the root element of the doxygen XML document.
    """
    if not hasattr(setup, 'DOXYGEN_ROOT'):
        setup.DOXYGEN_ROOT = ET.Element("root")  # dummy
    return setup.DOXYGEN_ROOT


def get_doxygen_id_index():
    """Return a dict mapping every ``@id`` in the merged doxygen tree to
    the element that owns it. Built lazily on first use and memoized
    on the :func:`setup` function object.

    Profiling a serial build at full MOM6 input (XML_PROGRAMLISTING=YES,
    109 MB merged tree) showed ``xmlutils.visit_ref`` burning 250 s of
    self time -- 27% of total wall clock -- in a single ``findall('.//*
    [@id=X]')`` call that linearly scanned the entire merged tree once
    per ``<ref>`` in prose. This index turns that scan into an O(1)
    dict lookup. Same shape of fix as the scanNode `//` -> `.//` patch
    in commit 8a217135e.
    """
    if not hasattr(setup, 'DOXYGEN_ID_INDEX'):
        root = get_doxygen_root()
        index = {}
        for el in root.iter():
            eid = el.get('id')
            if eid is not None:
                index[eid] = el
        setup.DOXYGEN_ID_INDEX = index
    return setup.DOXYGEN_ID_INDEX


def setup(app):
    import sphinx
    from .autodoc import (
        DoxygenMethodDocumenter,
        DoxygenTypeDocumenter,
        DoxygenModuleDocumenter,
    )
    from .autosummary import DoxygenAutosummary, DoxygenAutoEnum
    from .autosummary.generate import process_generate_options
    from .autodoxysource import AutoDoxySourceDirective

    app.connect("builder-inited", set_doxygen_xml)
    app.connect("builder-inited", process_generate_options)

    app.setup_extension('sphinx.ext.autodoc')
    app.setup_extension('sphinx.ext.autosummary')

    app.add_autodocumenter(DoxygenModuleDocumenter)
    app.add_autodocumenter(DoxygenMethodDocumenter)
    app.add_autodocumenter(DoxygenTypeDocumenter)

    app.add_config_value("doxygen_xml", "", 'env')
    # Used in autodoc_doxygen/autosummary/generate.py
    app.add_config_value('autosummary_toctree', '', 'html')

    app.add_directive('autodoxysummary', DoxygenAutosummary)
    app.add_directive('autodoxyenum', DoxygenAutoEnum)
    app.add_directive('autodoxysource', AutoDoxySourceDirective)

    app.add_css_file('autodoxysource.css')

    return {'version': sphinx.__display_version__, 'parallel_read_safe': True}
