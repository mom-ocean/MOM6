# Generated documentation

We use [Doxygen](http://www.doxygen.org/) for in-code documentation of APIs (i.e. arguments to subroutines/functions and members of types).
The guide for using doxygen in MOM6 is hosted on the [MOM6 developer's wiki](https://github.com/NOAA-GFDL/MOM6/wiki/Doxygen).

You can preview documentation by running doxygen or sphinx.

## Sphinx based documentation

The full documentation can be generated locally with
```bash
make html
```
which will generate html in `docs/_build/html/`. Start at `docs/_build/html/index.html`.

## Doxygen generated HTML

The doxygen generated HTML can be obtained locally (and slightly more quickly) with
```bash
make nortd SPHINXBUILD=false
```
which will generate html in `docs/APIs/`. Start at `docs/APIs/index.html`. If doxygen is not already available this will install a local copy of doxygen.

## Dependencies

If you do not have doxygen, to build a local version of the doxygen processor you will need the following packages:
- cmake
- g++ (or a c++ compiler)
- flex
- bison
- graphviz

(e.g. `apt-get install cmake g++ flex bison graphviz`)

If you are building the full generated sphinx documentation you will need the following packages in addition to those for doxygen above:
- libxml2-dev
- libxslt-dev

(.e.g `apt-get install libxml2-dev libxslt-dev`)

Before running sphinx (`make html`) you will need to issue:
```bash
pip install -r requirements.txt
```

## Credits

The sphinx documentation of MOM6 is made possible by modifications by [Angus Gibson](https://github.com/angus-g) to two packages, [sphinx-fortran](https://github.com/angus-g/sphinx-fortran) and [autodoc_doxygen](https://github.com/angus-g/sphinxcontrib-autodoc_doxygen).
