# CI script pipeline-ci-tool.sh

pipeline-ci-tool.sh contains functions corresponding to each job within the gitlab CI pipeline for MOM6 at GFDL, specifically on the gaea HPC.
Each function can be run by a parser function so that the functions can be invoked from the command line or a shell.
Some functions take arguments.
Encapsulating the job commands in a function allows us to develop/debug the pipeline by issuing the same, relatively short, commands at the command line.

pipeline-ci-tool.sh relies on three environment variables to execute. They are mandatory.
  - JOB_DIR is a scratch location that will be created and populated
  - CI_PROJECT_DIR is normally set by gitlab and will point to the working directory where MOM6 is cloned
  - CI_COMMIT_SHA is the commit of MOM6 to be tested

To use pipeline-ci-tool.sh interactively from an existing MOM6 clone, you could use
  `JOB_DIR=tmp CI_PROJECT_DIR=. CI_COMMIT_SHA=`git rev-parse HEAD` .gitlab/pipeline-ci-tool.sh ...`
This will use the HEAD commit in the current working dir and setup an independent test suite under tmp/.

## Usage
  `pipeline-ci-tool.sh FUNCTION [-x|+x] [-n|+n] [ARG1] [ARG2] [...]`
  `pipeline-ci-tool.sh FUNCTION [-x|+x] [-n|+n] [ARG1] [ARG2] [[-x|+x] [-n|+n] FUNCTION [ARG1] [ARG2] [...]] [...]`

FUNCTION can be one of
  - `create-job-dir` : Create a "job directory" using the environment variable JOB_DIR. This is a where all the compilation and running takes place.
  - `clean-job-dir` : Not used by .gitlab-ci.yml but useful for resetting an interactive session.
  - `copy-test-space LABEL` : Within $JOB_DIR, clones MOM6-examples to tmp-MOM6-examples-LABEL to use as a workspace for tests
  - `mrs-compile TARGET` : Invokes tools/MRS/Makefile.build to build MODE_VENDER. VENDER can be gnu, intel, or pgi. MODE can be repro, debug, static, etc.
  - `nolibs-ocean-only-compile VENDER` : Compiles the "no libraries" executables. These are not used elsewhere in the CI but check we have no namespace problems. VENDER can be gnu, intel, or pgi.
  - `run-suite VENDER CODE` : runs subsets of the MOM6-examples according to CODE using the VENDER executables. CODE is a string of the characters S (symmetric), N (non-symmetric), L (layout), D (debug), or R (restart), and if present executes the corresponding tests.
  - `check-stats VENDER CODE` : check the stats files for the corresponding VENDOR/CODE resulting from run-suite
  - `check-params VENDER CODE` : check the parameter documentation files for the corresponding VENDOR/CODE resulting from run-suite
  - `check-diags VENDER CODE` : check the available diagnostics files for the corresponding VENDOR/CODE resulting from run-suite

Options:
  - `-x` : shows commands as they are executed. `+x` turns back to silent executions. You can precede each function as needed so that only commands from selected functions are shown.
  - `-n` : for many function, disables all functionality and simply prints the banner that each sections was reached. `+n` turns the functions back on.

## Correspondance to jobs in .gitlab-ci.yml

The .gitlab-ci.yml jobs names and pipeline-ci-tool.sh commands are:

  clone:
    `pipeline-ci-tool.sh create-job-dir`

  work-space:pgi:
    `pipeline-ci-tool.sh copy-test-space pgi`

  work-space:intel:
    `pipeline-ci-tool.sh copy-test-space intel`

  work-space:gnu:
    `pipeline-ci-tool.sh copy-test-space gnu`

  work-space:gnu-restarts:
    `pipeline-ci-tool.sh copy-test-space gnu-rst`

  compile:pgi:repro:
    `pipeline-ci-tool.sh mrs-compile repro_pgi`

  compile:intel:repro:
    `pipeline-ci-tool.sh mrs-compile repro_intel`

  compile:gnu:repro:
    `pipeline-ci-tool.sh mrs-compile repro_gnu mrs-compile static_gnu`

  compile:gnu:debug:
    `pipeline-ci-tool.sh mrs-compile debug_gnu`

  compile:gnu:ocean-only-nolibs:
    `pipeline-ci-tool.sh nolibs-ocean-only-compile gnu`

  compile:gnu:ice-ocean-nolibs:
    `pipeline-ci-tool.sh nolibs-ocean-ice-compile gnu`

  run:pgi:
    `pipeline-ci-tool.sh run-suite pgi SNL`

  run:intel:
    `pipeline-ci-tool.sh run-suite intel SNL`

  run:gnu:
    `pipeline-ci-tool.sh run-suite gnu SNLD`

  run:gnu-restarts:
    `pipeline-ci-tool.sh run-suite gnu R`

  t:pgi:symmetric:
    `pipeline-ci-tool.sh check-stats pgi S`

  t:pgi:non-symmetric:
    `pipeline-ci-tool.sh check-stats pgi N`

  t:pgi:layout:
    `pipeline-ci-tool.sh check-stats pgi L`

  t:pgi:params:
    `pipeline-ci-tool.sh check-params pgi S`

  t:intel:symmetric:
    `pipeline-ci-tool.sh check-stats intel S`

  t:intel:non-symmetric:
    `pipeline-ci-tool.sh check-stats intel N`

  t:intel:layout:
    `pipeline-ci-tool.sh check-stats intel L`

  t:intel:params:
    `pipeline-ci-tool.sh check-params intel S`

  t:gnu:symmetric:
    `pipeline-ci-tool.sh check-stats gnu S`

  t:gnu:non-symmetric:
    `pipeline-ci-tool.sh check-stats gnu N`

  t:gnu:layout:
    `pipeline-ci-tool.sh check-stats gnu L`

  t:gnu:static:
    `pipeline-ci-tool.sh check-stats gnu T`

  t:gnu:symmetric-debug:
    `pipeline-ci-tool.sh check-stats gnu D`

  t:gnu:restart:
    `pipeline-ci-tool.sh check-stats gnu R`

  t:gnu:params:
    `pipeline-ci-tool.sh check-params gnu S`

  t:gnu:diags:
    `pipeline-ci-tool.sh check-diags gnu S`

### Duplicating the pipeline interactively

You can run a sequence of commands as follows. The setup and compile phases of the CI pipeline can be summarized with
```
pipeline-ci-tool.sh create-job-dir copy-test-space pgi copy-test-space intel copy-test-space gnu copy-test-space gnu-rst mrs-compile repro_pgi mrs-compile repro_intel mrs-compile repro_gnu mrs-compile static_gnu mrs-compile debug_gnu nolibs-ocean-only-compile gnu nolibs-ocean-ice-compile gnu
```

The run stage (works on compute nodes only) can be summarized with:
```
pipeline-ci-tool.sh run-suite pgi SNL run-suite intel SNL run-suite gnu SNLDT run-suite gnu R
```

The test stage is summarized by:
```
pipeline-ci-tool.sh check-stats pgi S check-stats pgi N check-stats pgi L check-params pgi S check-stats intel S check-stats intel N check-stats intel L check-params intel S check-stats gnu S check-stats gnu N check-stats gnu L check-stats gnu T check-stats gnu D check-stats gnu R check-params gnu S check-diags gnu S
```
