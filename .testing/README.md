# .testing

This directory contains the Makefile and test configurations used to evaluate
submissions to the MOM6 codebase.  The tests are designed to run either locally
or in a CI environment such as Travis.


## Overview

This section gives a very brief overview of the test suite and how to use it.

To build and run the model tests:
```
make -j
make -j test
```
For new users, the default configuration should be suitable for most platforms.
If not, then the following options may need to be configured.

`MPIRUN` (*default:* `mpirun`)

  Name of the MPI launcher.  Often this is `mpirun` or `mpiexec` but may all
  need to run through a scheduler, e.g. `srun` if using Slurm.

`DO_REGRESSION_TESTS` (*default: none*)

  Set to `true` to compare output with `dev/gfdl`.

`DO_REPRO_TESTS` (*default: none*)

  Set to `true` to compare DEBUG and REPRO builds, which typically correspond
  to unoptimized and optimized builds.  See TODO for more information.

These settings can either be specified at the command line, as shown below
```
make DO_REGRESSION_TESTS=true
make test DO_REGRESSION_TESTS=true
```
or saved in a configuration file, `config.mk`.

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
See "Tests" and "Test Configurations" for the complete list of tests.

The rest of the document describes the test suite in more detail, including
names and descriptions of the test classes and configurations.


## Testing overview

The test suite checks for numerical consistency of the model output across
different model configurations when subjected to relevant numerical and
mathematical transformations, such as grid layout or dimensional rescaling.  If
the model state is unchanged after each transformation, then the test is
reported as passing.  Any discrepancy in the model state causes the test to
fail.

Model state is currently defined by the `ocean.stats` output file, which
reports the total energy (per unit mass) at machine precision alongside similar
global metrics at lower precision, such as mass or mean sea level.

Diagnostics are based on the MOM checksum function, which includes the mean,
minimum, and maximum values, alongside a bitcount checksum, in the physical
domain, which are saved in the `chksum_diag` output file.


## Build configuration

The test suite defines a DEBUG and a REPRO build, which resemble targets used
at GFDL.  The DEBUG build is intended for detecting potential errors and
troubleshooting, while the REPRO build has typically been optimized for
production runs.

Ideally, the DEBUG and REPRO runs will produce identical results, although this
is often not the case for many compilers and platforms.  The `DO_REPRO_TEST`
flag is used to test DEBUG/REPRO equivalency.

The following options are provided to configure your compiler flags.

`FCFLAGS_DEBUG` (*default:* `-g -O0`)

  Specify the flags used in the DEBUG build.  These are the flags used for all
  tests excepting the REPRO builds.  They are also used to build the FMS
  library.

  These should be used to enable settings favorable to debugging, such as no
  optimizations, backtraces, range checking, and warnings.

  For more aggressive debugging flags which cannot be used with FMS, see
  `FCFLAGS_INIT`.

`FCFLAGS_REPRO:` (*default:* `-g -O2`)

  Specify the optimized reproducible run, typically used in production runs.

  Ideally, this should consist of optimization flags which improve peformance
  but do not change model output.  In practice, this is difficult to achieve,
  and should only used in certain environments.

`FCFLAGS_INIT` (*default: none*)

  This flag was historically used to specify variable initialization, such as
  nonzero integers or floating point values, and is still generally used for
  this purpose.

  As implemented, it is used for all MOM6 builds.  It is not used for FMS
  builds, so can also act as a debugging flag independent of FMS.

`FCFLAGS_COVERAGE` (*default: none*)

  This flag is used to define a build which supports some sort of code
  coverage, often one which is handled by the CI.

  For many compilers, this is set to `--coverage`, and is applied to both the
  compiler (`FCFLAGS`) and linker (`LDFLAGS`).

Example values used by GFDL and Travis for the GFortran compiler are shown
below.
```
FCFLAGS_DEBUG="-g -O0 -Wextra -Wno-compare-reals -fbacktrace -ffpe-trap=invalid,zero,overflow -fcheck=bounds"
FCFLAGS_REPRO="-g -O2 -fbacktrace"
FCFLAGS_INIT="-finit-real=snan -finit-integer=2147483647 -finit-derived"
FCFLAGS_COVERAGE="--coverage"
```

Note that the default values for these flags are very minimal, in order to
ensure compatibility over the largest possible range of compilers.

Like all configuration variables, these can be specified in a `config.mk` file.


## Building the executables

Run `make` to build the test executables.
```
make
```
This will fetch the MKMF build toolchain, fetch and compile the FMS framework
library, and compile the executables used in the test suite.  The default
configuration uses the symmetric grid in the debug-compile mode, with
optimizations disabled and stronger quality controls.  The following
executables will be created.

- `build/symmetric/MOM6`: Symmetric grid configuration (i.e. extended grids
  along western and/or southern boundaries for selected fields).  This is the
  default configuration.

- `build/asymmetric/MOM6`: Non-symmetric grid (equal-sized grids)

- `build/repro/MOM6`: Optimized reproducible mode

- `build/target/MOM6`: A reference build for regression testing

- `build/openmp/MOM6`: OpenMP-enabled build

The `target` and `repro` builds are only created when their respective tests
are set to `true`.


### Regression testing

When regression tests are enabled, the Makefile will check out a second copy of
the repository from a specified URL and branch given by `MOM_TARGET_URL` and
`MOM_TARGET_BRANCH`, respectively.  The code is checked out into the
`TARGET_CODEBASE` directory.

The default settings, with resolved values as comments, are shown below.
```
MOM_TARGET_SLUG = NOAA-GFDL/MOM6
MOM_TARGET_URL = https://github.com/$(MOM_TARGET_SLUG)
              #= https://github.com/NOAA-GFDL/MOM6
MOM_TARGET_LOCAL_BRANCH = dev/gfdl
MOM_TARGET_BRANCH = origin/$(MOM_TARGET_LOCAL_BRANCH)
                 #= origin/dev/gfdl
TARGET_CODEBASE = $(BUILD)/target_codebase
```
These default values can be configured to target a particular development
branch.

Currently the target can only be specifed by branch name, rather than hash.

New diagnostics do not report as a fail, and are not tracked by any CIs, but
the test will report a warning to the user.


## Tests

Using `test` will run through the full test suite.
```
make test
```
The tests are gathered into the following groups.

- `test.regressions`: Regression tests relative to a code state (when enabled)
- `test.grids`: Symmetric vs nonsymmetric grids
- `test.layouts`: Domain decomposition, based on parallelization
- `test.restarts`: Resubmission by restarts
- `test.repros`: Optimized (REPRO) and unoptimized (DEBUG) compilation
- `test.nans`: NaN initialization of allocated arrays
- `test.dims`: Dimensional scaling (length, time, thickness, depth)

Each group of tests can also be run individually, such as in the following
example.
```
make test.grids
```

Each configuration is tested relative to the `symmetric` build, and reports a
fail if the answers differ from this build.


## Test configurations

The following model test configurations (TCs) are supported, and are based on
configurations in the MOM6-examples repository.

- `tc0`: Unit testing of various model components, based on `unit_tests`
- `tc1`: A low-resolution version of the `benchmark` configuration
  - `tc1.a`: Use the un-split mode with Runge-Kutta 3 time integration
  - `tc1.b`: Use the un-split mode with Runge-Kutta 2 time integration
- `tc2`: An ALE configuration based on tc1 with tides
  - `tc2.a`: Use sigma, PPM_H4 and no tides
- `tc3`: An open-boundary condition (OBC) test based on `circle_obcs`
- `tc4`: Sponges and initialization using I/O


## Code coverage

Code coverage reports the lines of code which have been tested, and can be used
to determine if a particular section is untested.

Coverage is measured using `gcov` and is reported for TCs using the `symmetric`
executable.

Coverage reporting is optionally uploaded to the `codecov.io` site.
```
https://codecov.io/gh/NOAA-GFDL/MOM6
```
This is disabled on default, but can be enabled by the `REPORT_COVERAGE` flag.
```
make test REPORT_COVERAGE=true
```
Note that any uploads will require a valid CodeCov token.


## Running on Travis

Whenever code is pushed to GitHub or a pull request (PR) is created, the test
suite is triggered and the code changes are tested.

When the tests are run on Travis, the following variables are re-defined:

- `DO_REPRO_TESTS` is set to `true` for all tests.

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
