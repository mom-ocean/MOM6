# Generated documentation

We use [Doxygen](http://www.doxygen.org/) for in-code documentation of APIs (i.e. arguments to subroutines/functions and members of types).
The guide for using doxygen in MOM6 is hosted on the [MOM6 developer's wiki](https://github.com/NOAA-GFDL/MOM6/wiki/Doxygen).

The full documentation can be generated locally with
```bash
cd docs/
make html
```
which will generate html in `docs/_build/html/`. Start at `docs/_build/html/index.html`.

The generated documentation can be obtained locally with
```bash
cd docs/
doxygen Doxyfile_nortd
```
which will generate html in `docs/APIs/`. Start at `docs/APIs/index.html`.
