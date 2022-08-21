===============
MOM6 Test Suite
===============

This directory contains test configurations used to evaluate submissions to the
MOM6 codebase.  The tests are designed to run either locally or in a CI
environment.


Usage
=====

``make -j``
   Build the FMS library and test executables.

``make -j test``
   Run the test suite, defined in the ``tc`` directores.

``make clean.build``
	Delete only the MOM6 test executables.

``make clean``
   Delete the MOM6 test executables and dependency builds (FMS).


Configuration
=============

The test suite includes many configuration flags and variables which can be
configured at either the command line, or can be stored in a ``config.mk``
file.

Several of the following may require configuration for particular systems.

``MPIRUN`` (*default:* ``mpirun``)
   Name of the MPI launcher.  Often this is ``mpirun`` or ``mpiexec`` but may
   all need to run through a scheduler, e.g. ``srun`` if using Slurm.

``FRAMEWORK`` (*default:* ``fms1``)
   Select either the legacy FMS framework (``fms1``) or an FMS2 I/O compatible
   version (``fms2``).

``DO_REPRO_TESTS`` (*default:* *none*)
   Set to ``true`` to test the REPRO build and confirm equivalence of DEBUG and
   REPRO builds.

   For compilers with aggressive optimization, DEBUG and REPRO may not produce
   identical results and this test should not be used.

``DO_REGRESSION_TESTS`` (*default:* *none*)
   Set to ``true`` to compare output with a defined target branch, set by
   ``MOM_TARGET_LOCAL_BRANCH``.  (NOTE: This defaults to ``dev/gfdl``).

``DO_COVERAGE`` (*default:* *none*)
   Set to ``true`` to enable code coverage.  Currently only configured for
   ``gcov``.

``REQUIRE_COVERAGE_UPLOAD`` (*default:* *none*)
   Set to ``true`` if failure to upload the coverage report to codecov.io
   should result in an error.  This should only be enabled if codecov.io has
   already been configured for the user, or by a supporting CI.

``DO_PROFILE`` (*default:* *none*)
   Set to ``true`` to enable performance profile monitoring.  Models are
   compiled using ``OPT_FCFLAGS`` (see below) and performance of various
   functions are reported and compared to the target branch.

   Results from these tests should only be considered if the platform has been
   configure for benchmarking.


Build configuration
-------------------

Compilation is controlled with the following variables.  Defaults are chosen
for the widest compatibility across platforms.  Users should modify these to
reflect their own needs.

``FCFLAGS_DEBUG`` (*default:* ``-g -O0``)
   The "DEBUG" build, for rapid compilation and debugging.

``FCFLAGS_REPRO`` (*default:* ``-g -O2``)
   The "REPRO" build, for reproducible production runs.

``FCFLAGS_OPT`` (*default:* ``-g -O3``)
   The "OPT" build, for aggressive optimization and profiling.

``FCFLAGS_COVERAGE`` (*default:* ``-g -O0 -fbacktrace --coverage``)
   Flags used for producing code coverage reports.  Defaults are for gcc,
   although ``--coverage`` is relatively common across compilers.

``FCFLAGS_INIT`` (*default:* *none*)
   A placeholder flag for aggressive initialization testing.  This is appended
   to existing flags (usually ``FCFLAGS_DEBUG``).

``FCFLAGS_FMS`` (*default:* ``FCFLAGS_DEBUG``)
   Compiler flags used for the supporting FMS library.  In most cases, it is
   sufficient to use ``FCFLAGS_DEBUG``.

``LDFLAGS_COVERAGE`` (*default:* ``--coverage``)
   Linker flags to enable coverage.

``LDFLAGS_USER`` (*default:* *none*)
   A placeholder for supplemental linker flags, such as an external library not
   configured by autoconf.

The following flags are passed as environment variables to other Makefiles.

``FC``, ``MPIFC``
   The Fortran compiler and its MPI wrapper.

``CC``, ``MPICC``
   The C compiler and its MPI wrapper.  This is primarily used by FMS, but may
   be used in some MOM6 autoconf tests.

If unset, these will be configured by autoconf or from the user's environment
variables.

Additional settings for particular tasks are explained below.


Example ``config.mk``
---------------------

An example config.mk file configured for GFortran is shown below.::

   DO_REPRO_TESTS = true
   DO_REGRESSION_TESTS = true
   DO_COVERAGE = true
   DO_PROFILE = true

   FCFLAGS_DEBUG = -g -O0 -Wextra -Wno-compare-reals -fbacktrace -ffpe-trap=invalid,zero,overflow -fcheck=bounds
   FCFLAGS_REPRO = -g -O2 -fbacktrace
   FCFLAGS_OPT = -g -O3 -mavx -mfma
   FCFLAGS_INIT = -finit-real=snan -finit-integer=2147483647 -finit-derived
   FCFLAGS_COVERAGE = --coverage

The file follows Makefile syntax, so quotations are generally not required and
spaces are permitted between assignment operators (``=``).


Builds
======

Run ``make`` to build the test executables.::

   $ make

This will fetch external dependencies, compile the FMS framework library, and
compile the executables used in the test suite.

The following executables will be created.

``build/symmetric/MOM6``
   Use symmetric grids for model fields, using DEBUG flags.

   A symmetric grid is one where each finite-volume cell has grid points along
   all faces.  Often this results in a redundant row of points along each side
   of a regular domain.

   This is the recommended production configuration, and is the reference build
   for all tests in the suite.

``build/asymmetric/MOM6``
   Use asymmetric grids for model fields.

   Asymmetric grids eliminate a redundant fields along western and southern
   boundaries, which reduces the total number of points.  They also ensure
   that center, face, and vertex field arrays are the same size.

   The disadvantages are greater computational complexity along these
   boundaries.  They also do not support open boundary conditions.

   Asymmetric grids were traditionally used in many legacy ocean models.

``build/repro/MOM6``
   Optimized build for doing reproducible runs, based REPRO flags.

   This is only built if ``DO_REPRO_TESTS`` is set to ``true``.

``build/target/MOM6``
   A reference build for regression testing.

   The reference branch is set by ``MOM_TARGET_LOCAL_BRANCH``.  This would
   generally be configured by a CI to a pull request's target branch.  This is
   only built if ``DO_REGRESSION_TESTS`` is set to ``true``.

``build/openmp/MOM6``
   A DEBUG build with OpenMP enabled.


Tests
=====

The ``test`` rule will run all of the tests.::

   $ make test

Tests are based on configurations which are designed to give identical output.
When the output differs, the test reports a failure.


Test groups
-----------

The tests are gathered into the following groups.

``test.grid``
   Compare symmetric and nonsymmetric grids.

``test.regression``
   Compare the current codebase to a target branch (e.g. ``dev/gfdl``).

``test.layout``
   Compare a serial (one domain) and a parallel (two domain) simulation.

``test.restart``
   Compare a single run to two runs separated by a restart.

``test.repro``
   Compare the unoptimized (DEBUG) and optimized (REPRO) builds.

``test.nan``
   Enable NaN-initialization of allocated (heap) arrays.

   This relies on internal features of glibc and may not work on other
   platforms.

``test.dim``
   Enable dimension rescaling tests.

Each tests uses the symmetric build for its reference state.

These rules can be used to run individual groups of tests.::

   $ make test.grid


Test experiments
----------------

For each group, we test each of the following configurations, which represent
idealizations of various production experiments.

``tc0``
   Unit testing of various model components, based on ``unit_tests``

``tc1``
   A low-resolution version of the ``benchmark`` configuration

   ``tc1.a``
      Use the un-split mode with Runge-Kutta 3 time integration

   ``tc1.b``
      Use the un-split mode with Runge-Kutta 2 time integration

``tc2``
   An ALE configuration based on tc1 with tides

   ``tc2.a``
      Use sigma, PPM_H4 and no tides

``tc3``
   An open-boundary condition (OBC) test based on ``circle_obcs``

``tc4``
   Sponges and initialization using I/O


Test procedure
--------------

The test suite checks for numerical consistency of the model output across
different model configurations when subjected to relevant numerical and
mathematical transformations, such as grid layout or dimensional rescaling.  If
the model state is unchanged after each transformation, then the test is
reported as passing.  Any discrepancy in the model state causes the test to
fail.

Model state is currently defined by the ``ocean.stats`` output file, which
reports the total energy (per unit mass) at machine precision alongside similar
global metrics at lower precision, such as mass or mean sea level.

Diagnostics are based on the MOM checksum function, which includes the mean,
minimum, and maximum values, alongside a bitcount checksum, in the physical
domain, which are saved in the ``chksum_diag`` output file.


Regression testing
==================

When ``DO_REGRESSION_TESTS`` is enabled, the Makefile will check out a second
copy of the repository from a specified URL and branch given by
``MOM_TARGET_URL`` and ``MOM_TARGET_BRANCH``, respectively.  The code is
checked out into the ``TARGET_CODEBASE`` directory.

The default settings, with resolved values as comments, are shown below.::

   MOM_TARGET_SLUG = NOAA-GFDL/MOM6
   MOM_TARGET_URL = https://github.com/$(MOM_TARGET_SLUG)
                 #= https://github.com/NOAA-GFDL/MOM6
   MOM_TARGET_LOCAL_BRANCH = dev/gfdl
   MOM_TARGET_BRANCH = origin/$(MOM_TARGET_LOCAL_BRANCH)
                    #= origin/dev/gfdl
   TARGET_CODEBASE = $(BUILD)/target_codebase

These default values can be configured to target a particular development
branch.

Currently the target can only be specified by branch name, rather than hash.

New diagnostics do not report as a fail, and are not tracked by any CIs, but
the test will report a warning to the user.


Code coverage
=============

Code coverage reports the lines of code which have been tested, and can be used
to determine if a particular section is untested.

To enable code coverage, set ``DO_COVERAGE`` to ``true``.

Reports are stored in the build directories.  There is one report per source
file, and each ends in the ``.gcov`` suffix.  Two sets of coverage reports are
generated.

``build/cov``
   Test suite code coverage

``build/unit``
   Unit test code coverage

To upload the tests to codecov.io, use the following rules.::

   $ make report.cov             # Test suite
   $ make report.cov.unit        # Unit test

Note that any uploads will require a valid CodeCov token.  If uploading through
the CI, this can be set up through your GitHub account.

Pull request coverage reports for the CI can be checked at
https://codecov.io/gh/NOAA-GFDL/MOM6


CI configuration
================

Whenever code is pushed to GitHub or a pull request (PR) is created, the test
suite is run.

When the tests are run on the CI, the following variables are re-defined:

- ``DO_REPRO_TESTS`` is set to ``true`` for all tests.

- ``DO_REGRESSION_TESTS`` is set to ``true`` for a PR submission, and is unset for
  code pushes.

- ``DO_COVERAGE`` is set to ``true``.

   - For pull requests, ``REQUIRE_COVERAGE_UPLOAD`` is set to ``true``.

- ``MOM_TARGET_SLUG`` is set to the URL stub of the model to be built.

  For submissions to NOAA-GFDL, this will be set to ``NOAA-GFDL/MOM6`` and the
  reference URL will be ``https://github.com/NOAA-GFDL/MOM6``.

- ``MOM_TARGET_LOCAL_BRANCH``

  For a code push, this is set to the name of the active branch at GitHub.  For
  a PR, this is the name of the branch which is receiving the PR.
