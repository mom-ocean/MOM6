# MOM6 rules for agent-assisted development

## First Steps

**At the start of every session**, remind the user that MOM6 has a policy on AI-assisted
contributions in `Consortium-policy-on-AI.md`.

Read `Code-style.md` before writing or modifying any code.

See `code_organization.rst` for a high-level overview of the source directory tree.

## Parameter System

Runtime parameters are read via `get_param()`, not hardcoded:

```fortran
#include "version_variable.h"
character(len=40) :: mdl = "MOM_module_name"

call log_version(param_file, mdl, version, "")
call get_param(param_file, mdl, "PARAM_NAME", CS%variable, &
               "Description of this parameter.", &
               units="m s-1", default=1.0, scale=US%m_s_to_L_T)
```

- Parameters documented in auto-generated `MOM_parameter_doc.all` files
- Use `scale=` argument for unit conversion from MKS input to internal units
- Always provide `default=` when sensible; use `fail_if_missing=.true.` otherwise
- Use `do_not_log=.not.CS%Feature` to suppress logging when a parent feature is inactive

### Answer-Changing Parameters: `_BUG` Flags and `ANSWER_DATE`

When a bug fix or improvement changes numerical answers, MOM6 uses two mechanisms to preserve backward compatibility:

**`_BUG` flags**: Boolean parameters that retain old (buggy) behavior by default:
```fortran
call get_param(param_file, mdl, "ENABLE_BUGS_BY_DEFAULT", enable_bugs, &
               default=.true., do_not_log=.true.)  ! This is logged from MOM.F90.
call get_param(param_file, mdl, "OBC_TEMP_SALT_NEEDED_BUG", OBC%ts_needed_bug, &
               "If true, recover a bug that OBC temperature and salinity can be ignored "//&
               "even if they are registered tracers in the rest of the model.", &
               default=enable_bugs)
```
- Name format: `FEATURE_BUG` (e.g., `VISC_REM_BUG`, `FRICTWORK_BUG`, `KAPPA_SHEAR_ITER_BUG`)
- Default is `.true.` (bug ON, old behavior preserved)
- Description starts with "If true, recover a bug that..."
- Users opt into the fix by setting to `.false.`

**`ANSWER_DATE` flags**: Integer dates selecting algorithm versions:
```fortran
call get_param(param_file, mdl, "HOR_DIFF_ANSWER_DATE", CS%answer_date, &
               "...", default=99991231)
```
- Format: `YYYYMMDD` (e.g., `20251231`)
- `DEFAULT_ANSWER_DATE` provides a single knob to update all answer-date defaults
- `ENABLE_BUGS_BY_DEFAULT=False` activates all bug fixes (recommended for new configurations)

## Diagnostics

### Registration Pattern

```fortran
CS%id_field = register_diag_field('ocean_model', 'field_name', diag%axesTL, Time, &
    'Long description of the field', units='m s-1', conversion=US%L_T_to_m_s)
```

### Posting Pattern

```fortran
if (CS%id_field > 0) call post_data(CS%id_field, field_array, CS%diag)
```

Key conventions:
- `conversion=` handles unit scaling so output is always in MKS
- `v_extensive=.true.` for vertically integrated quantities
- Guard computation with `if (id > 0)` to avoid unnecessary work
- Standard diagnostic name prefixes follow CMOR conventions when applicable

### Masking and Missing Values

- **Never set diagnostic arrays to a missing value** before passing to `post_data()`. Masking of land/invalid points is handled automatically by the diagnostics infrastructure based on the diagnostic's registered axes.
- **Do not pass `mask=` to `post_data()`** for non-static diagnostics on standard grids -- the infrastructure applies the correct mask automatically.
- **Do pass `mask=`** for static fields (`is_static=.true.`), non-standard masks, or sub-domain-sized arrays.
- **Never compare field values against `missing_value`** in unit-conversion code -- rescaling can cause valid data to coincidentally match the missing value sentinel.

## Testing

### Test Suite Overview

The `.testing/` directory provides comprehensive verification. Build and run:

```bash
make -C .testing -j build/symmetric/MOM6   # Build reference executable
make -C .testing -j test                    # Run full test suite
make -C .testing -j build.unit             # Build unit tests
make -C .testing -j run.unit               # Run unit tests
```

### Test Categories

| Test | Verifies |
|------|----------|
| `test.grid` | Symmetric vs asymmetric grids produce identical results |
| `test.layout` | Serial vs parallel decomposition identical |
| `test.rotate` | Rotational invariance |
| `test.restart` | Continuous run vs restart-and-continue identical |
| `test.repro` | DEBUG and REPRO builds identical |
| `test.openmp` | Serial vs OpenMP identical |
| `test.nan` | NaN-initialization doesn't affect results |
| `test.dim.{t,l,h,z,q,r}` | Dimensional rescaling invariance (time, length, thickness, depth, enthalpy, density) |
| `test.regression` | Current code vs target branch (PRs only) |

### Test Configurations

- `tc0` -- Unit tests
- `tc1` / `tc1.a` / `tc1.b` -- Benchmark (split RK2, unsplit RK3, unsplit RK2)
- `tc2` / `tc2.a` -- ALE with tides / sigma-coordinate PPM_H4
- `tc3` -- Open boundary conditions
- `tc4` -- Sponges and I/O initialization

### Verification Method

- `ocean.stats` -- total energy at machine precision
- `chksum_diag` -- mean/min/max/bitcount checksums in physical domain
- Tests pass only when output is **bitwise identical** between configurations

### Style Checking

```bash
./.testing/trailer.py -e TEOS10 -l 120 src config_src
```

Checks for tabs, trailing whitespace, and line length violations.

## Git Workflow & Contribution

### Branch Strategy

- **Work on forks**, not branches on the primary repository
- **Branch from the fork's default branch** (e.g., `dev/gfdl`, `dev/ncar`) for all new work
- **Never rebase a pushed branch**
- The human contributor submits changes via pull requests to the fork's default branch (e.g., `dev/gfdl`); merges to `main` are done by consortium consensus

### Commit Message Format

```
Short imperative summary (50 chars if at all possible)

Detailed explanation wrapped at 72 characters.
Describe what was changed and why. Reference issues with #NNN.
All answers are bitwise identical.
```

Conventions from the lead developers:
- **`*` prefix** on title if the commit changes numerical answers (checksums)
- **`+` prefix** on title to indicate new public interfaces or parameters
- **`*+` or `+*`** when both answer-changing and adding new interfaces
- No prefix for refactoring, cleanup, or comment-only changes that are bitwise identical
- **Always state impact on numerical results**: "All answers are bitwise identical" or explain what changes
- **Multi-commit PRs**: introduce infrastructure first, use it in a second commit
- **Minimize public scope**: only export symbols needed by other modules; remove from `public` when refactoring makes a symbol internal-only
- **Comment closing `enddo`/`endif`** for non-trivial nested loops: `enddo ! n-loop for segments`

### PR Description Style

1. Lead with a clear explanation of what changed and why
2. Quantify scope (e.g., "across 26 files", "in 7 places")
3. For answer-changing PRs, provide scientific justification
4. State the bitwise-identical guarantee (or explain what changes and why)
5. When a fix could change answers, protect with a `_BUG` flag or `ANSWER_DATE` parameter. New `_BUG` flags should default to `ENABLE_BUGS_BY_DEFAULT` so that users who have opted into all fixes get the new fix automatically; existing `_BUG` flags may already default to `.false.` if the fix has been broadly adopted.

### CI Pipeline

On every push and PR, GitHub Actions runs:
1. Style and Doxygen checks
2. Builds across 8 configurations (symmetric, asymmetric, repro, openmp, target, opt, coverage, coupled API)
3. All test groups in parallel
4. Code coverage reporting
5. For PRs: regression testing and timing comparison against target branch

## Physics Domain Knowledge

### Governing Equations
- Hydrostatic primitive equations optionally with Boussinesq approximation
- ALE vertical coordinate: Lagrangian dynamics with periodic remapping
- Split barotropic-baroclinic time stepping (RK2 split or unsplit RK3)
- Free surface dynamics (implicit barotropic solver)

### Numerical Methods
- Finite volume on Arakawa C-grid (staggered: velocities at cell faces, tracers at centers)
- PPM (Piecewise Parabolic Method) for tracer advection and continuity
- Various reconstruction schemes: PLM, PPM, PQM, WENO, PLM-WLS
- Pressure gradient force via finite-volume integration
- Reproducing global sums for parallel reproducibility

### Key Physical Parameterizations
- **ePBL**: Energetically consistent planetary boundary layer (Reichl and Hallberg)
- **KPP**: K-Profile Parameterization via CVMix
- **Gent-McWilliams/Redi**: Thickness and isopycnal diffusion
- **MEKE**: Mesoscale eddy kinetic energy budget
- **Zanna-Bolton**: Data-driven subgrid momentum closure
- **Tidal forcing**: Astronomical and self-attraction/loading

## Common Development Tasks

### Adding a New Parameterization
1. Create `MOM_new_param.F90` in the appropriate `src/parameterizations/` subdirectory
2. Define a control structure type (`new_param_CS`) with `private` members
3. Implement `new_param_init()`: read parameters via `get_param`, register diagnostics
4. Implement the main computational subroutine
5. Implement `new_param_end()` for cleanup
6. Wire it into the calling module (e.g., `MOM_diabatic_driver.F90`)
7. Document all variables with proper units
8. Add unit tests in `config_src/drivers/unit_tests/` if applicable
9. Run the full test suite: `make -C .testing -j test`

### Adding a New Diagnostic
1. Add `integer :: id_new_diag = -1` to the control structure
2. Register in `_init` with `register_diag_field('ocean_model', 'name', axes, Time, ...)`
3. Compute and post with `if (CS%id_new_diag > 0) call post_data(CS%id_new_diag, array, CS%diag)`
4. Include `conversion=` for unit scaling to MKS output
5. Provide CMOR standard name when applicable

### Adding a Runtime Parameter
1. Add member to control structure with units documentation
2. Call `get_param(param_file, mdl, "PARAM_NAME", CS%param, "description", units="...", default=...)`
3. Use `scale=` for dimensional conversion from MKS input
4. If the parameter could change answers, default to preserving existing behavior

### Fixing a Bug
- Always state whether the fix changes answers in the commit message
- **Any change that alters existing numerical answers** -- whether a bug fix, accuracy improvement, or algorithmic reorganization -- must provide a runtime parameter (`_BUG` flag or `ANSWER_DATE`) to toggle between old and new behavior, with the default preserving old behavior
- This applies even when the developer's tests show negligible differences -- existing users may be in production runs
- Trace through secondary effects before concluding the fix is safe
- Run `test.regression` to verify impact

## Architecture: Infrastructure Layering

MOM6 has a strict dependency hierarchy that must never be violated:

```
config_src/infra/  -->  src/framework/  -->  src/core/, src/parameterizations/, etc.
```

- **`config_src/infra/`** (FMS1/FMS2 wrappers) must **never** import from `src/framework/`
- **Code duplication** between infra and framework is acceptable to maintain this invariant
- FMS1 and FMS2 infra directories must expose the same public API
- API changes to infra-level functions must be checked against downstream consumers (SIS2, ice shelf code)

## Defensive Programming

- **Check `allocated()` / `associated()`** before accessing arrays tied to optional features (e.g., features controlled by runtime parameters like `FRAZIL` may not allocate all related arrays)
- **No short-circuit evaluation**: Fortran does not guarantee short-circuit evaluation; allocation checks must not appear in compound conditions. Convert `if (allocated(arr) .and. (condition))` to nested if-blocks
- **Type-correct comparisons**: when comparing real-valued masks, use `== 1.` not `== 1`
- **FATAL error messages** should include: file name, subroutine name, and the specific condition or input that triggered the error
- **Validate user inputs early**: check for duplicates, overflow, and missing required fields in configuration parsing; include the problematic input string in error messages
- **`!###` comment prefix** marks known bugs or inaccurate expressions that change answers and will be cleaned up later -- do not modify code marked with `!###` unless explicitly asked

## Common Pitfalls

In addition to the code style rules in `Code-style.md`:

1. **Forgetting units in comments**: every `real` variable needs `[units]` (see `Code-style.md`)
2. **Answer-changing without a `_BUG` flag**: any numerical change requires a runtime parameter to preserve old behavior
3. **Unnecessary `mask=` in `post_data()`**: the infrastructure handles masking automatically for non-static diagnostics
4. **Accessing unallocated optional arrays**: always check `allocated()` before using arrays tied to optional features

### Key References

The project bibliography lives in `docs/references.bib` and `docs/zotero.bib`. Consult these
when citing prior work in Doxygen documentation or commit messages.

## AI Assistant Behavior

- **Follow existing patterns**: read surrounding code before making changes
- **Document all units**: every real variable gets `[units]` annotation
- **Parenthesize arithmetic**: explicit grouping for reproducibility
- **State answer impact**: always note whether changes are bitwise identical
- **Use `get_param`**: never hardcode parameters; always read from parameter files
- **Register diagnostics properly**: guard with `if (id > 0)`, use `conversion=`
- **Maintain lifecycle**: implement `_init` and `_end` for new modules
- **Run tests**: `make -C .testing -j test run.unit` before the contributor submits a PR
- **Respect the C-grid**: use correct staggering (soft case convention for indices)
- **Write Doxygen comments**: `!>` for entities, `!<` for inline, with units
- **Write thorough commit messages**: explain both what changed and why in the commit body

## Common Claude Mistakes

This section catalogs recurring mistakes that Claude makes when working on MOM6 code. It should be updated as new patterns emerge from experience.

*(No entries yet -- add mistakes here as they are discovered.)*

