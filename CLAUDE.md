# MOM6 rules for agent-assisted development

## Purpose of this document

This file guides [Claude Code](https://docs.anthropic.com/en/docs/claude-code) -- Anthropic's agentic
command-line tool for software development -- when working on the MOM6 codebase. Claude Code uses Claude,
Anthropic's AI assistant, to read and edit files, run shell commands, search codebases, and interact with
git and GitHub, all from the terminal. It is available via `npm install -g @anthropic-ai/claude-code`
or at https://claude.ai/download.

When Claude Code operates inside a repository, it reads this `CLAUDE.md` file automatically
to learn project-specific conventions, coding standards, and development workflows. Everything
below captures the patterns, rules, and best practices that govern MOM6 development -- drawn
from the official wiki, ReadTheDocs documentation, CI infrastructure, and the coding style
of important historical PRs. Following these guidelines ensures that AI-generated contributions
match the quality and consistency expected by the MOM6 community.

This is a **living document** intended to be constantly updated as both the MOM6 codebase
and Claude's capabilities evolve. As developers gain experience using Claude Code on MOM6,
common mistakes and new best practices should be added here. See the
"Common Claude Mistakes" section at the end of this file -- it is expected to grow over time.

## Project Overview

MOM6 (Modular Ocean Model, version 6) is a next-generation open-source ocean model developed by NOAA-GFDL. It combines the best of GOLD and MOM5 into a modern Fortran codebase solving the primitive equations for ocean dynamics on an Arakawa C-grid. Key features:

- Arbitrary Lagrangian-Eulerian (ALE) vertical coordinate
- Boussinesq and non-Boussinesq modes
- Flexible equation of state (Wright, TEOS-10, linear)
- Comprehensive parameterization library (ePBL, KPP, lateral mixing, tidal forcing)
- Coupled to SIS2 (sea ice), ice shelves, and Earth system models via FMS/NUOPC/MCT
- Dimensional unit scaling for consistency testing

### Language & Environment

- **Language**: Fortran 2003+ (free-form `.F90`, preprocessed)
- **Build systems**: Autoconf (primary), CMake (experimental), GNU Make for testing
- **Dependencies**: FMS framework, CVMix, GSW-Fortran (TEOS-10)
- **Compilers**: Must compile under GNU, Intel, and PGI
- **Testing**: Comprehensive CI via GitHub Actions and GFDL GitLab pipeline
- **Repository**: https://github.com/NOAA-GFDL/MOM6
- **Documentation**: https://mom6.readthedocs.io/en/main/
- **Examples**: https://github.com/NOAA-GFDL/MOM6-examples/wiki
- **Main branch**: `dev/gfdl` (GFDL development); `main` (inter-lab coordination)

## Code Style & Conventions

### Formatting Rules

- **Indentation**: 2 spaces (no tabs, ever)
- **Continuation lines**: minimum 4 spaces indent
- **Line length**: target 100 characters for code; absolute maximum 120 (enforced by `.testing/trailer.py`)
- **No trailing whitespace** (enforced by CI)
- **`enddo`** and **`endif`** (single words); but `end module`, `end subroutine`, `end function`, `end type`
- Space after language tokens: `if (x > 0)` not `if(x > 0)`
- No space between function name and `(`: `call fn(x)` not `call fn (x)`
- Space around assignment `=`; but no spaces in loop bounds: `do i=is,ie`
- Named arguments: `call fn(arg_name=val)` (no spaces around `=`)

### Module Structure

Every module follows this pattern:

```fortran
!> Brief module description
module MOM_module_name

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_some_module, only : specific_symbol
use MOM_other,       only : other_thing

implicit none ; private

#include <MOM_memory.h>

public :: exported_routine_1, exported_routine_2

!> Control structure for this module
type, public :: module_CS ; private
  real :: param         !< Description [units]
  integer :: id_diag = -1 !< Diagnostic ID for some_field
end type module_CS

contains

!> Initialize the module, read parameters, register diagnostics
subroutine module_init(Time, G, GV, US, param_file, diag, CS)
  ! ... call log_version, get_param, register_diag_field ...
end subroutine module_init

!> Deallocate module memory
subroutine module_end(CS)
  ! ... cleanup ...
end subroutine module_end

!> \namespace MOM_module_name
!! Extended description, references, and equations
end module MOM_module_name
```

Key rules:
- **`implicit none ; private`** on a single line in every module
- **Explicit `only` imports** on all `use` statements -- no blanket imports
- **No global/module data** -- all state lives in control structures passed as arguments
- **All arguments must have declared `intent`** (`in`, `out`, or `inout`); pointers exempt
- **Every module has `_init` and `_end` subroutines** for lifecycle management

### Naming Conventions

- **Files**: `MOM_something.F90` (module inside is `MOM_something`)
- **Variables**: `snake_case` for multi-word names (Doxygen-compatible)
- **Control structures**: `module_CS` (e.g., `energetic_PBL_CS`), always `private`
- **Diagnostic IDs**: `id_diag_name`, initialized to `-1`
- **Inverses**: prefix with `I` (e.g., `IdxCu` = `1/dxCu`, `IareaT` = `1/areaT`)
- **Grid objects**: `G` (ocean_grid_type), `GV` (verticalGrid_type), `US` (unit_scale_type)
- **Public functions**: self-documenting names; private helpers may use short names

### Loop Index Conventions (Soft Case)

This is a critical MOM6 convention for the Arakawa C-grid staggering:

- **Lowercase `i`, `j`, `k`**: cell-center / layer-center (tracer points). Example: `h(i,j,k)`, `T(i,j,k)`
- **Uppercase `I`, `J`**: cell-edge / staggered points (north-east convention). `I` = i+½. Example: `u(I,j,k)`, `v(i,J,k)`, `q(I,J,k)`
- **Uppercase `K`**: vertical interface. `K=1` is above layer `k=1`; `K` = k-½. Example: `e(i,j,K)`
- **Domain bounds**: `is, ie, js, je` (computational); `isd, ied, jsd, jed` (data/halo)
- **Grid locations**: T (tracer/h-points), Cu (u-points), Cv (v-points), Bu (q/vertex-points)

### Memory Macros

Array dimensions use preprocessor macros from `MOM_memory.h`:
- `SZI_(G)`, `SZJ_(G)`, `SZK_(GV)` for explicit-shape arrays
- `NIMEM_`, `NJMEM_`, `NKMEM_` for allocatable arrays

## Unit Documentation (Critical Convention)

MOM6 uses a dimensional annotation system for every real variable. Units are documented in square brackets at the end of comments, using a two-part notation:

```
[rescaled_dimensions ~> MKS_equivalent]
```

### Dimensional Symbols

| Symbol | Physical Dimension | MKS Unit |
|--------|-------------------|----------|
| Z | Vertical depth/distance | m |
| L | Horizontal length | m |
| T | Time | s |
| H | Layer thickness | m (Boussinesq) or kg m-2 |
| R | Density | kg m-3 |
| Q | Enthalpy | J kg-1 |
| C | Temperature | degC |
| S | Salinity | ppt |
| A | Arbitrary/generic units | a |

### Examples

```fortran
real :: velocity     !< Horizontal velocity [L T-1 ~> m s-1]
real :: pressure     !< Hydrostatic pressure [R L2 T-2 ~> Pa]
real :: thickness    !< Layer thickness [H ~> m or kg m-2]
real :: diffusivity  !< Vertical diffusivity [Z2 T-1 ~> m2 s-1]
real :: slope        !< Isopycnal slope [Z L-1 ~> nondim]
real :: efficiency   !< Mixing efficiency [nondim]
real :: field        !< A field in arbitrary units [A]
real :: Z_to_m       !< Scaling factor [m Z-1 ~> 1]
```

### Unit Annotation Rules

1. **Every real variable** must have units in `[brackets]` at the end of its comment
2. **Canonical symbol ordering**: consistent order (e.g., `H L2` not `L2 H`)
3. **Boussinesq variants first**: `[H ~> m or kg m-2]` when units differ by approximation
4. **Simplified expressions only**: write `[T2 Z-1 ~> s2 m-1]`, not `[H Z T-1 / H Z2 T-3 = T2 Z-1 ~> s2 m-1]`
5. **Exponent notation**: `m-1`, `s-2`, `kg-3` (no slashes like `1/m`)
6. **No extra spaces** inside brackets
7. **Nondimensional**: use `[nondim]`
8. **Arbitrary/generic**: use `[A]` or `[A ~> a]`, never `[arbitrary]`
9. **Scaling factors**: `[target source-1 ~> 1]`, e.g., `[Z m-1 ~> 1]`

## Arithmetic and Reproducibility

MOM6 demands bitwise reproducibility across processor counts and build modes. These rules are non-negotiable:

1. **Parenthesize all additions**: `z = (a + b) + c` not `z = a + b + c`
2. **Never use `sum()`, `prod()`, or `matmul()` intrinsics** -- operation order is undefined
3. **Pre-compute reciprocals**: `Q * G%IareaT(i,j)` not `Q / G%areaT(i,j)`
4. **Never write `B / C * D`**: use `(B * D) / C` (explicit grouping)
5. **Avoid the exponent operator `**`**: write `a * a * a` not `a**3` (compilers emit transcendental `pow()`)
6. **Avoid transcendental functions** where possible (sin, cos, log, non-integer powers are implementation-dependent)
7. **`sqrt()` is safe** (IEEE-754 exactly rounded); use MOM6's `cuberoot` for cube roots
8. **Explicit parenthesization for FP precision**: group unit-conversion factors before multiplying data
   ```fortran
   tmp = ( ( tv%C_p * GV%H_to_RZ ) * h(i,j,k) ) * tv%T(i,j,k)
   ```
9. **Vanished layer pattern**: `h + h_neglect` (not `max(h, h_neglect)`)

### Array Syntax

- **Prohibited**: `tv%S = 0.` (scalar-looking whole-array assignment)
- **Allowed**: `tv%S(:,:,:) = 0.` (explicit colon notation)
- **Prohibited**: array-syntax math on arrays that include halos (halo data may be invalid)

## Doxygen Documentation

### Comment Syntax

- `!>` for documentation comments on the following entity
- `!<` for inline documentation on the preceding entity (same line)
- `!!` for multi-line continuation (no blank lines between)
- `!>@{` and `!>@}` for grouping related declarations

### Requirements

- **All public subroutines/functions**: `!>` header describing purpose
- **All arguments**: documented with `!<` or `!>` including units
- **All type members**: documented with `!<` including units
- **All real variables**: must include physical description and units
- **Equations**: LaTeX with `\f$ ... \f$` (inline) or `\f[ ... \f]` (display)
- **Extended descriptions**: placed before `end module` using `\namespace`

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
- **Branch from `dev/gfdl`** for all new work
- **Never rebase a pushed branch**
- Submit changes via pull requests to `dev/gfdl`

### Commit Message Format

```
Short imperative summary (aim for ~50 chars)

  Detailed explanation indented by 2 spaces, wrapped at ~80 columns.
  Describe what was changed and why. Reference issues with #NNN.
  All answers are bitwise identical.
```

Conventions from the lead developers:
- **`*` prefix** on title if the commit changes numerical answers (checksums)
- **`+` prefix** on title to indicate new features/additions
- **Always state impact on numerical results**: "All answers are bitwise identical" or explain what changes
- **Multi-commit PRs**: introduce infrastructure first, use it in a second commit

### PR Description Style

1. Lead with a clear explanation of what changed and why
2. Quantify scope (e.g., "across 26 files", "in 7 places")
3. For answer-changing PRs, provide scientific justification
4. State the bitwise-identical guarantee (or explain what changes and why)
5. When a fix could change answers, protect with a runtime parameter defaulting to `.false.`

### CI Pipeline

On every push and PR, GitHub Actions runs:
1. Style and Doxygen checks
2. Builds across 8 configurations (symmetric, asymmetric, repro, openmp, target, opt, coverage, coupled API)
3. All test groups in parallel
4. Code coverage reporting
5. For PRs: regression testing and timing comparison against target branch

Additionally, GFDL's internal GitLab runs ~400 tests across 59 configurations on Gaea HPC (GNU/Intel/PGI, debug/repro).

## Source Directory Structure

```
src/
  core/                          # Dynamical core
    MOM.F90                      # Main stepping routines (~5000 lines)
    MOM_barotropic.F90           # Barotropic solver (~6700 lines)
    MOM_continuity_PPM.F90       # PPM-based continuity
    MOM_dynamics_split_RK2.F90   # Split RK2 time stepping
    MOM_grid.F90                 # Horizontal grid type
    MOM_variables.F90            # Common variable types
    MOM_verticalGrid.F90         # Vertical grid type
    MOM_PressureForce_FV.F90     # Finite-volume pressure gradient
  ALE/                           # Vertical remapping/regridding
    MOM_ALE.F90                  # ALE driver
    MOM_regridding.F90           # Vertical grid generation
    MOM_remapping.F90            # Conservative remapping
    Recon1d_*.F90                # 1D reconstruction schemes (PLM, PPM, PQM, WENO, etc.)
  diagnostics/                   # Diagnostic calculations
    MOM_diagnostics.F90          # Standard diagnostics
    MOM_diagnose_MLD.F90         # Mixed layer depth
  equation_of_state/             # EOS implementations
    MOM_EOS.F90                  # EOS wrapper
    MOM_EOS_Wright.F90           # Wright (1997) EOS
    MOM_EOS_TEOS10.F90           # TEOS-10 EOS
  framework/                     # Infrastructure
    MOM_diag_mediator.F90        # Diagnostics framework
    MOM_file_parser.F90          # Parameter file parsing
    MOM_unit_scaling.F90         # Dimensional scaling system
    MOM_domains.F90              # Domain decomposition
    MOM_restart.F90              # Restart I/O
  ice_shelf/                     # Ice shelf dynamics
  initialization/                # Grid/state initialization
  ocean_data_assim/              # Data assimilation
  parameterizations/
    lateral/                     # Lateral parameterizations
      MOM_hor_visc.F90           # Horizontal viscosity
      MOM_thickness_diffuse.F90  # Thickness diffusion (GM)
      MOM_MEKE.F90               # Mesoscale eddy kinetic energy
      MOM_Zanna_Bolton.F90       # Zanna-Bolton backscatter
    vertical/                    # Vertical mixing
      MOM_energetic_PBL.F90      # ePBL mixed layer (~4500 lines)
      MOM_CVMix_KPP.F90          # KPP via CVMix
      MOM_diabatic_driver.F90    # Diabatic processes driver
      MOM_set_diffusivity.F90    # Background diffusivity
      MOM_vert_friction.F90      # Vertical friction
  tracer/                        # Tracer transport and specific tracers
  user/                          # Idealized configuration initialization
config_src/
  drivers/
    solo_driver/                 # Ocean-only standalone (simplest; testing)
    FMS_cap/                     # GFDL coupler
    nuopc_cap/                   # NUOPC/CESM coupling
    unit_tests/                  # Unit test executables
    timing_tests/                # Performance benchmarks
  memory/
    dynamic_symmetric/           # Symmetric memory layout (default)
    dynamic_nonsymmetric/        # Asymmetric memory layout
  infra/                         # FMS1/FMS2 wrappers
  external/                      # Null hooks for optional components
pkg/
  CVMix-src/                     # Community Vertical Mixing
  GSW-Fortran/                   # TEOS-10 Gibbs Seawater
```

## Physics Domain Knowledge

### Governing Equations
- Primitive equations with hydrostatic or Boussinesq approximation
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
- **ePBL**: Energetically consistent planetary boundary layer (Hallberg)
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
- If it could change answers, consider a runtime parameter defaulting to `.false.` (old behavior)
- Trace through secondary effects before concluding the fix is safe
- Run `test.regression` to verify impact

## Common Pitfalls

1. **Forgetting units in comments**: every `real` variable needs `[units]`
2. **Unparenthesized arithmetic**: causes non-reproducible results across compilers
3. **Using `**` operator**: triggers transcendental `pow()` -- write explicit multiplications
4. **Module-level data**: no globals; pass everything through arguments
5. **Missing `only` on imports**: all `use` statements require explicit imports
6. **Array syntax with halos**: halo data is not guaranteed valid; use explicit loops
7. **Blanket `use` imports**: never `use module` without `only`
8. **Tabs in source**: CI will fail; use spaces only
9. **Trailing whitespace**: CI will fail
10. **`sum()` intrinsic**: operation order undefined across compilers

## Helpful Resources

- MOM6 documentation: https://mom6.readthedocs.io/en/main/
- MOM6 developers wiki: https://github.com/NOAA-GFDL/MOM6/wiki
- MOM6 code style guide: https://github.com/NOAA-GFDL/MOM6/wiki/Code-style-guide
- MOM6 Doxygen conventions: https://github.com/NOAA-GFDL/MOM6/wiki/Doxygen
- MOM6 examples wiki: https://github.com/NOAA-GFDL/MOM6-examples/wiki
- MOM6 repository policies: https://github.com/NOAA-GFDL/MOM6-examples/wiki/MOM6-repository-policies
- MOM6 developer workflow: https://github.com/NOAA-GFDL/MOM6/wiki/Developer-workflow
- MOM6 runtime parameters: https://github.com/NOAA-GFDL/MOM6/wiki/MOM6-run-time-parameter-system
- MOM6 forum: https://bb.cgd.ucar.edu/cesm/forums/mom6.148/
- CVMix documentation: https://github.com/CVMix/CVMix-src
- TEOS-10 (GSW): http://www.teos-10.org/
- GOTM (General Ocean Turbulence Model): https://gotm.net/

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
- **Run tests**: `make -C .testing -j test` before any PR
- **Respect the C-grid**: use correct staggering (soft case convention for indices)
- **Write Doxygen comments**: `!>` for entities, `!<` for inline, with units

## Common Claude Mistakes

This section catalogs recurring mistakes that Claude makes when working on MOM6 code.
It should be updated as new patterns emerge from experience.

*(No entries yet -- add mistakes here as they are discovered.)*
