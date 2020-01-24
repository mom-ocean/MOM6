# .testing

This directory contains the Makefile and test configurations used to evaluate
submissions to the MOM6 codebase.  The tests are designed to run either locally
or in a Travis-CI.


## Overview

This section gives a very brief overview of the test suite and how to use it.

To build and run the model tests
```
make
make test
```

Regression testing is disabled on default.  To include regression tests:
```
make DO_REGRESSION_TESTS=true
make test DO_REGRESSION_TESTS=true
```

On platforms other than Gaea, a MKMF build template may be required.  To
specify the path to the template:
```
make MKMF_TEMPLATE=/path/to/template.mk
```

To run individual classes of tests, use the subclass name:
```
make test.grids
make test.layouts
make DO_REGRESSION_TESTS=true test.regressions
```

To test an individual test configuration (TC):
```
make tc0.grid
```

The rest of the document describes the test suite in more detail, including
names and descriptions of the test classes and configurations.


## Testing overview

The test suite consists of many comparisons of model output for different model
configurations when subjected to relevant numerical and mathematical
transformations, such as grid layout or dimensional rescaling, for which the
model output should be invariant.  If the model state is unchanged after each
transformation, then the test is reported as passing.  Any discrepancy in the
model state causes the test to fail.

Model state is currently defined by the `ocean.stats` output file, which
reports the total energy (per unit mass) at machine precision alongside similar
global metrics, such as mass or mean sea level, at lower precision.

Checksums for every available diagnostic are also compared and the Makefile
will report any differences, but such differences are not yet considered a fail
condition.


## Building the executables

Run `make` to build the test executables.
```
make
```
This will fetch the MKMF build toolchain, fetch and compile the FMS framework
library, and compile the executables used in the test suite.  The default
configuration uses the symmetric grid in the debug-compile mode, with
optimizations disabled and stronger quality controls.  The following
executables will be created:

- `build/symmetric/MOM6`: Symmetric grid configuration (extended grids along
  western and/or southern boundaries).  This is the default configuration.

- `build/asymmetric/MOM6`: Non-symmetric grid (equal-sized grids)

- `build/repro/MOM6`: Optimized reproducible mode

- (optional) `build/target/MOM6`: A reference build for regression testing

The `target` build is only created when the `DO_REGRESSION_TESTS` flag is set
to `true`:
```
make DO_REGRESSION_TESTS=true
```
When set, the build will check out a second copy of the repository from a
specified URL and branch given by `MOM_TARGET_URL` and `MOM_TARGET_BRANCH`,
respectively.  The code is checked out into the `TARGET_CODEBASE` directory.

The current default settings are
```
MOM_TARGET_SLUG = NOAA-GFDL/MOM6
MOM_TARGET_URL = https://github.com/$(MOM_TARGET_SLUG)
#              = https://github.com/NOAA-GFDL/MOM6
MOM_TARGET_LOCAL_BRANCH = dev/gfdl
MOM_TARGET_BRANCH = origin/$(MOM_TARGET_LOCAL_BRANCH)
#                 = origin/dev/gfdl
TARGET_CODEBASE = $(BUILD)/target_codebase
```
These default values can be configured to target a particular development
branch.


#### MKMF template

The MKMF build toolchain requires a template file when building the model.  The
default template, `ncrc-gnu.mk`, is part of the MKMF repository, but has been
specifically configured for use on NOAA's Gaea computer, and other institutes
will require their own template files.

The template can be specified as a Make flag.
```
make MKMF_TEMPLATE=/path/to/template.mk
```
The `linux-ubuntu-xenial-gnu.mk` template is provided in the `.testing`
directory, and is intended for Travis-CI builds, but may also be a good
reference point for other Linux distributions.

In the future, this step may be replaced with a more generalized build system,
such as CMake or automake.


## Tests

Using `test` will run through the test suite.
```
make test
```
This will run through the following tests:

- `test.regressions`: Regression tests relative to a code state (when enabled)
- `test.grids`: Symmetric vs nonsymmetric grids
- `test.layouts`: Domain decomposition, based on parallelization
- `test.restarts`: Resubmission by restarts
- `test.repros`: Optimized (REPRO) and unoptimized (DEBUG) compilation
- `test.nans`: NaN initialization of allocated arrays
- `test.dims`: Dimensional scaling (length, time, thickness, depth)

To enable the regression tests, use `DO_REGRESSION_TEST=true`.
```
make test DO_REGRESSION_TESTS=true
```

Each test can also be run individually.  For example, the following command
will only run the grid tests.
```
make test.grids
```

Each configuration is tested relative to the `symmetric` build, and reports a
fail if the answers differ from this build.


## Test configurations

The following test configurations (TCs) are supported:

- tc0: Unit testing of various model components, based on `unit_tests`
- tc1: A low-resolution version of the `benchmark` configuration
  - tc1.a: Use the un-split mode with Runge-Kutta 3 time integration
  - tc1.b: Use the un-split mode with Runge-Kutta 2 time integration
- tc2: An ALE configuration based on tc1 with tides
  - tc2.a: Use sigma, PPM_H4 and no tides
- tc3: An open-boundary condition (OBC) test based on `circle_obcs`


## Code coverage

Code coverage reports the lines of code which have been tested, and can
explicitly demonstrate when a particular operation is untested.

Coverage is measured using `gcov` and is reported for TCs using the `symmetric`
executable.

Coverage reporting is optionally sent to the `codecov.io` site.
```
https://codecov.io/gh/NOAA-GFDL/MOM6
```
This is disabled on default, but can be enabled by the `REPORT_COVERAGE` flag.
```
make test REPORT_COVERAGE=true
```
Note that any uploads will require a valid token generated by CodeCov.


## Running on Travis

Whenever code is pushed to GitHub or a pull request (PR) is created, the test
suite is triggered and the code changes are tested.

When the tests are run on Travis, the following variables are re-defined:

- `DO_REGRESSION_TESTS` is set to `true` for a PR submission, and is unset for
  code pushes.

- `MOM_TARGET_SLUG` is set to `TRAVIS_REPO_SLUG`, the URL stub of the model to
  be built.

  For submissions to NOAA-GFDL, this will be set to `NOAA-GFDL/MOM6` and the
  reference URL will be `https://github.com/NOAA-GFDL/MOM6`.

- `MOM_TARGET_LOCAL_BRANCH` is set to `TRAVIS_BRANCH`.

  For a code push, this is set to the name of the active branch at GitHub.  For
  a PR, this is the name of the branch which is receiving the PR.

- `REPORT_COVERAGE` is set to `true`.

## Running under slurm

By default the executables are invoked using `mpirun`. Under slurm you might need to
use `srun` (such as on GFDL's gaea HPC):
```
make MPIRUN=srun test
```
