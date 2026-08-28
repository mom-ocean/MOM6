# Makefile for Sphinx documentation

# You can set these variables from the command line.
SPHINXOPTS       ?=
SPHINXBUILD       = sphinx-build
PAPER            ?=
BUILDDIR         ?= _build
DOXYGEN_BIN      ?= doxygen
DOXYGEN_CONF     ?=
FINDDOXYGENBIN    = $(DOXYGEN_BIN)
DOXYGEN_RELEASE  ?= Release_1_8_13

# User-friendly check for sphinx-build
ifeq ($(shell which $(SPHINXBUILD) >/dev/null 2>&1; echo $$?), 1)
$(error The '$(SPHINXBUILD)' command was not found. Make sure you have Sphinx installed, then set the SPHINXBUILD environment variable to point to the full path of the '$(SPHINXBUILD)' executable. Alternatively you can add the directory with the executable to your PATH. If you don\'t have Sphinx installed, grab it from http://sphinx-doc.org/)
endif

# Set DOXYGENBIN variable if the executable is found.  If it is not found, then we set it to a makefile rule that will attempt to build it.
ifeq ($(shell which $(FINDDOXYGENBIN) >/dev/null 2>&1; echo $$?), 0)
DOXYGENBIN = $(shell which $(FINDDOXYGENBIN))
else
DOXYGENBIN = doxygen/bin/doxygen
endif

# Internal variables.
PAPEROPT_a4     = -D latex_paper_size=a4
PAPEROPT_letter = -D latex_paper_size=letter
ALLSPHINXOPTS   = -d $(BUILDDIR)/doctrees $(PAPEROPT_$(PAPER)) $(SPHINXOPTS)
# the i18n builder cannot share the environment and doctrees with the others
I18NSPHINXOPTS  = $(PAPEROPT_$(PAPER)) $(SPHINXOPTS) .

.PHONY: help clean html dirhtml singlehtml pickle json htmlhelp qthelp devhelp epub latex latexpdf text man changes linkcheck doctest gettext nortd nortd_latex all

help:
	@echo "Please use \`make <target>' where <target> is one of"
	@echo "  html       to make standalone HTML files"
	@echo "  dirhtml    to make HTML files named index.html in directories"
	@echo "  singlehtml to make a single large HTML file"
	@echo "  pickle     to make pickle files"
	@echo "  json       to make JSON files"
	@echo "  htmlhelp   to make HTML files and a HTML help project"
	@echo "  qthelp     to make HTML files and a qthelp project"
	@echo "  devhelp    to make HTML files and a Devhelp project"
	@echo "  epub       to make an epub"
	@echo "  latex      to make LaTeX files, you can set PAPER=a4 or PAPER=letter"
	@echo "  latexpdf   to make LaTeX files and run them through pdflatex"
	@echo "  latexpdfja to make LaTeX files and run them through platex/dvipdfmx"
	@echo "  text       to make text files"
	@echo "  man        to make manual pages"
	@echo "  texinfo    to make Texinfo files"
	@echo "  info       to make Texinfo files and run them through makeinfo"
	@echo "  gettext    to make PO message catalogs"
	@echo "  changes    to make an overview of all changed/added/deprecated items"
	@echo "  xml        to make Docutils-native XML files"
	@echo "  pseudoxml  to make pseudoxml-XML files for display purposes"
	@echo "  linkcheck  to check all external links for integrity"
	@echo "  doctest    to run all doctests embedded in the documentation (if enabled)"
	@echo "  nortd      to run the doxygen html version of the APIs manual"
	@echo "  nortd_latex to run the doxygen pdf version of the APIs manual"

clean:
	rm -rf $(BUILDDIR)/* api/generated/* MOM6.tags xml

$(BUILDDIR):
	mkdir -p $@

html: $(DOXYGENBIN) $(BUILDDIR)
	$(SPHINXBUILD) -M html . $(BUILDDIR) $(ALLSPHINXOPTS) -j auto
# Optional equation post-processing (was a patch in jr3cermak/sphinx@v3.2.1mom6.4
# applied to sphinx/cmd/build.py); now invoked from the Makefile so we can run
# against stock upstream Sphinx. The nortd target below has the same hook.
ifeq ($(UPDATEHTMLEQS), Y)
	@echo "Post processing equations."
	@python3 ./postProcessEquations.py -d $(BUILDDIR) -p html -b sphinx -s index.html $(UPDATEHTMLEQSVERBOSE)
endif
	@echo
	@echo "Build finished. The HTML pages are in $(BUILDDIR)/html."

dirhtml:
	$(SPHINXBUILD) -M dirhtml . $(BUILDDIR) $(ALLSPHINXOPTS)
	@echo
	@echo "Build finished. The HTML pages are in $(BUILDDIR)/dirhtml."

singlehtml:
	$(SPHINXBUILD) -M singlehtml . $(BUILDDIR) $(ALLSPHINXOPTS)
	@echo
	@echo "Build finished. The HTML page is in $(BUILDDIR)/singlehtml."

pickle:
	$(SPHINXBUILD) -M pickle . $(BUILDDIR) $(ALLSPHINXOPTS)
	@echo
	@echo "Build finished; now you can process the pickle files."

json:
	$(SPHINXBUILD) -M json . $(BUILDDIR) $(ALLSPHINXOPTS)
	@echo
	@echo "Build finished; now you can process the JSON files."

htmlhelp:
	$(SPHINXBUILD) -M htmlhelp . $(BUILDDIR) $(ALLSPHINXOPTS)
	@echo
	@echo "Build finished; now you can run HTML Help Workshop with the" \
	      ".hhp project file in $(BUILDDIR)/htmlhelp."

qthelp:
	$(SPHINXBUILD) -M qthelp . $(BUILDDIR) $(ALLSPHINXOPTS)
	@echo
	@echo "Build finished; now you can run "qcollectiongenerator" with the" \
	      ".qhcp project file in $(BUILDDIR)/qthelp, like this:"
	@echo "# qcollectiongenerator $(BUILDDIR)/qthelp/MOM6.qhcp"
	@echo "To view the help file:"
	@echo "# assistant -collectionFile $(BUILDDIR)/qthelp/MOM6.qhc"

devhelp:
	$(SPHINXBUILD) -M devhelp . $(BUILDDIR) $(ALLSPHINXOPTS)
	@echo
	@echo "Build finished."
	@echo "To view the help file:"
	@echo "# mkdir -p $$HOME/.local/share/devhelp/MOM6"
	@echo "# ln -s $(BUILDDIR)/devhelp $$HOME/.local/share/devhelp/MOM6"
	@echo "# devhelp"

epub:
	$(SPHINXBUILD) -M epub . $(BUILDDIR) $(ALLSPHINXOPTS)
	@echo
	@echo "Build finished. The epub file is in $(BUILDDIR)/epub."

latex: $(DOXYGENBIN) $(BUILDDIR)
	$(SPHINXBUILD) -M latex . $(BUILDDIR) $(ALLSPHINXOPTS)
	@echo
	@echo "Build finished; the LaTeX files are in $(BUILDDIR)/latex."
	@echo "Run \`make' in that directory to run these through (pdf)latex" \
	      "(use \`make latexpdf' here to do that automatically)."

latexpdf: $(DOXYGENBIN) $(BUILDDIR)
	$(SPHINXBUILD) -M latex . $(BUILDDIR) $(ALLSPHINXOPTS)
	@echo "Running LaTeX files through pdflatex..."
	$(MAKE) -C $(BUILDDIR)/latex LATEXMKOPTS="-f -silent" all-pdf
	@echo "pdflatex finished; the PDF files are in $(BUILDDIR)/latex."

latexpdfja:
	$(SPHINXBUILD) -M latex . $(BUILDDIR) $(ALLSPHINXOPTS)
	@echo "Running LaTeX files through platex and dvipdfmx..."
	$(MAKE) -C $(BUILDDIR)/latex all-pdf-ja
	@echo "pdflatex finished; the PDF files are in $(BUILDDIR)/latex."

text:
	$(SPHINXBUILD) -M text . $(BUILDDIR) $(ALLSPHINXOPTS)
	@echo
	@echo "Build finished. The text files are in $(BUILDDIR)/text."

man:
	$(SPHINXBUILD) -M man . $(BUILDDIR) $(ALLSPHINXOPTS)
	@echo
	@echo "Build finished. The manual pages are in $(BUILDDIR)/man."

texinfo:
	$(SPHINXBUILD) -M texinfo . $(BUILDDIR) $(ALLSPHINXOPTS)
	@echo
	@echo "Build finished. The Texinfo files are in $(BUILDDIR)/texinfo."
	@echo "Run \`make' in that directory to run these through makeinfo" \
	      "(use \`make info' here to do that automatically)."

info:
	$(SPHINXBUILD) -M texinfo . $(BUILDDIR) $(ALLSPHINXOPTS)
	@echo "Running Texinfo files through makeinfo..."
	make -C $(BUILDDIR)/texinfo info
	@echo "makeinfo finished; the Info files are in $(BUILDDIR)/texinfo."

gettext:
	$(SPHINXBUILD) -b gettext $(I18NSPHINXOPTS) $(BUILDDIR)/locale
	@echo
	@echo "Build finished. The message catalogs are in $(BUILDDIR)/locale."

changes:
	$(SPHINXBUILD) -M changes . $(BUILDDIR) $(ALLSPHINXOPTS)
	@echo
	@echo "The overview file is in $(BUILDDIR)/changes."

linkcheck:
	$(SPHINXBUILD) -M linkcheck . $(BUILDDIR) $(ALLSPHINXOPTS)
	@echo
	@echo "Link check complete; look for any errors in the above output " \
	      "or in $(BUILDDIR)/linkcheck/output.txt."

doctest:
	$(SPHINXBUILD) -M doctest . $(BUILDDIR) $(ALLSPHINXOPTS)
	@echo "Testing of doctests in the sources finished, look at the " \
	      "results in $(BUILDDIR)/doctest/output.txt."

xml:
	$(SPHINXBUILD) -M xml . $(BUILDDIR) $(ALLSPHINXOPTS)
	@echo
	@echo "Build finished. The XML files are in $(BUILDDIR)/xml."

pseudoxml:
	$(SPHINXBUILD) -M pseudoxml . $(BUILDDIR) $(ALLSPHINXOPTS)
	@echo
	@echo "Build finished. The pseudo-XML files are in $(BUILDDIR)/pseudoxml."

doxygen:
	git clone -b $(DOXYGEN_RELEASE) --depth 1 https://github.com/doxygen/doxygen

doxygen/bin/doxygen: doxygen
	(cd $(<); cmake -DCMAKE_POLICY_VERSION_MINIMUM=3.5 -G "Unix Makefiles" .)
	(cd $(<); $(MAKE))

nortd: $(DOXYGENBIN) $(BUILDDIR)
ifdef DOXYGEN_CONF
	$< $(DOXYGEN_CONF)
else
	$< Doxyfile_nortd
endif
ifeq ($(UPDATEHTMLEQS), Y)
	@echo "Post processing equations."
	@python3 ./postProcessEquations.py -d $(BUILDDIR) -p APIs -b doxygen -s index.html $(UPDATEHTMLEQSVERBOSE)
endif
	@echo "Build finished. The html files are in $(BUILDDIR)/APIs."

nortd_latex: $(DOXYGENBIN) $(BUILDDIR)
ifdef DOXYGEN_CONF
	$< $(DOXYGEN_CONF)
else
	$< Doxyfile_nortd_latex
endif
	# Make sure we are in nonstopmode for generating the pdf by patching the Makefile
	# Older doxygen versions do not have LATEX_CMD
	-grep -q "pdflatex refman" $(BUILDDIR)/doxygen_latex/Makefile && sed -i.bak 's/pdflatex refman/pdflatex -interaction=nonstopmode refman/' $(BUILDDIR)/doxygen_latex/Makefile
	make -C $(BUILDDIR)/doxygen_latex LATEX_CMD="pdflatex -interaction=nonstopmode"

Clean: clean
	rm -rf $(BUILDDIR) xml MOM6.tags

CLEAN: Clean
	rm -rf doxygen
