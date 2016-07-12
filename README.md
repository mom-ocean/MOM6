NOAA-GFDL MOM6
==============

The Modular Ocean Model, MOM6.

This branch (gh-pages) is used for storing the source code documentation which can be viewed at http://NOAA-GFDL.github.io/MOM6/.

The files on this branch are mostly generated and we will re-write the history of this branch to keep the repository size under control.

The generated API documentation can be viewed http://NOAA-GFDL.github.io/MOM6/APIs/.

Updating the generated API documentation
========================================

The workflow for re-generating the API documentation with doxygen (http://www.doxygen.org/) is:
```bash
git checkout dev/master
git clone https://github.com/doxygen/doxygen
(cd doxygen/; cmake -G "Unix Makefiles" .)
(cd doxygen/; make -j 4)
./doxygen/bin/doxygen .doxygen
git checkout gh-pages
git reset --hard gh-pages-static
rm -rf APIs
mv html APIs
git add APIs
git commit -m "Generated API documentation with doxygen"
```
