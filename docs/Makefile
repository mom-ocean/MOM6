# Standalone Flint/Docutils documentation build.
PYTHON ?= python3
OUTPUT ?= docs/_build/flint-api
GRAPH_DEPTH ?= 1
ROOT := $(abspath $(dir $(lastword $(MAKEFILE_LIST)))/..)

.PHONY: help html check

help:
	@echo "make html   Generate and validate the documentation site"
	@echo "make check  Run builder/renderer tests (requires Graphviz)"
	@echo "Set PYTHON, OUTPUT (a new directory), and GRAPH_DEPTH as needed."
	@echo "OUTPUT is relative to the repository root; GRAPH_DEPTH=0 skips diagrams."

html:
	cd "$(ROOT)" && "$(PYTHON)" docs/build_docs.py --output "$(OUTPUT)" --graph-depth "$(GRAPH_DEPTH)"

check:
	cd "$(ROOT)" && "$(PYTHON)" -m unittest discover -s docs -p 'test_docs.py'
