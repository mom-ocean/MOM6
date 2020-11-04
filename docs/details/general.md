# General troublehshooting (FAQ)

## Latex

### doxygen
ERROR: `! Dimension too large`

There may be a pdf that is too large.  The default `\maxdimen`
is 16384 pts.  To check the size of various `pdf` documents is to use
`pdfinfo`.  In Ubuntu, the package `poppler-utils` contains
this tool.

The class graph for `MOM__memory_8h__dep__incl.pdf` gets very large
(22607 pts) and swamps `pdflatex` with a `DOT_GRAPH_MAX_NODES` set at
the maximum value of 10000.  The default value, 50, works.  A value 
of 75 also worked in a recent test.

## PDF: Contents

The generated PDF shows a Contents heading with no entries or is a
blank page.  This is indicative of excessive failures with latex.
Latex generally requires multiple runs to properly create the table
of contents, index, references, equations and pages.

## Python

Be very careful of string processing in python when handling
latex escaped math.  If a latex string for `\theta` is read in as
`\\theta`.  Mishandling the string can cause the `\\t` to become a
tab (`\t`) and cause mahem with the processing.

Reference: [string literals](https://docs.python.org/3/reference/lexical_analysis.html?#literals)
