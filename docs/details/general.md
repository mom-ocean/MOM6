# General troublehshooting (FAQ)

## Tags and References

### equation-

Avoid using the prefix `equation` for tags and references.  Actual equations within
the document are promoted with an `equation-` prefix.  While using them is not a
problem.  Just be aware that casual use of them may unexpectedly conflict with an
equation reference.

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
latex escaped math in equations.  If a latex string for `\theta` is read in as
`\\theta`.  Mishandling the string can cause the `\\t` to become a
tab (`\t`) and cause mahem with the processing.

Reference: [string literals](https://docs.python.org/3/reference/lexical_analysis.html?#literals)

# Debugging

## python3 virtual enviroment

Setup a virtual environment for processing:

```bash
python3 -m venv venv/mom6Doc
source venv/mom6Doc/bin/activate
# cd to the docs directory within the MOM6 repo
pip3 install -r requirements.txt
```

The `deactivate` command allows you to exit from the virtual environment.

NOTE: RTD will not upgrade the sphinx module if `#egg=` is specified in the `requirements.txt` file.

### debugging

A useful commnad line tool for debugging sphinx and extensions is the python debugger.
Add the following line to stop to any portion of the python code to create a break
point.

```python
import pdb; pdb.set_trace()
```

Run `make html` without redirection to a log file.

For only processing .dox files and some specific F90 files, edit and use the
`Doxyfile_rtd_dox` file.  This limits the document processing to fewer files and
allows for rapid testing.

`make clean; make html DOXYGEN_CONF=Doxyfile_rtd_dox UPDATEHTMLEQS=Y`

## Example execution

The following example assumes a virtual environment as setup above using `mom6Doc`.
A similar environment is possible using anaconda.

```
$ source venv/mom6Doc/bin/activate
(mom6Doc) $ cd docs
(mom6Doc) $ make clean
(mom6Doc) $ make html >& _build/html_log.txt
(mom6Doc) $ make latexpdf >& _build/latex_log.txt
```

The last command may appear to hang.  On error, latex will request input from the keyboard.
Press `R` and enter.  This will keep latex running to completion or stop after 100 errors
are reached.

Once the documentation is built, you can use a web browser to look around in the `_build`
directory.

## Local web server

Python provides a way to quickly stand up a private web server for checking documentation. It requires knowledge of
the IP address if you are using a remote server, otherwise `localhost` should work.

You can start the server on any port. Port 8080 is shown here as an example.
```bash
python3 -m http.server 8080
```

After starting the server, you can browse to the known IP using `http://IP:8080/` or if you are on the same
machine use `http://localhost:8080/`.
