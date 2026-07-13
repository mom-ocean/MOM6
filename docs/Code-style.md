# MOM6 code style guide

MOM6 makes use of Fortran2003 and later extensions but only when supported by all available compilers. We try to avoid the use of very modern Fortran constructs that can limit portability. To help keep the code readily understandable, this page makes recommendations about how to use or not use the various features of the modern Fortran language.

Code style is typically a personal choice, but when styles clash it can lead to discord.  These
standards have been adopted in an attempt to promote harmony and clarity.

## White space

- No tabs
- No trailing white space
- Indents should be consistent (same amount used in contiguous blocks)
- Preferred indent is 2 spaces
  - "preferred" might understate the reaction invoked by other indent amounts! :wink:
- Continuation lines should be indented at least 4 spaces, but more space can be used if it helps align lines
  - One exception is the continuation line after a simple `if` test with no `then`, which could be 2 space indent as though the `then` and `endif` were present.  I.e.,
```
        if (test) &
          var = expression
```
and
```
        if (test) then
          var = expression
        endif
```
should have the same indenting on `var = expression` so that there would be no need to change the indenting if another statement were added to the `if` block, or if the latter expression were converted to the former.
- No white space between a function name and the opening parenthesis
- White space after all language token
  - `if(a==0)` is legal fortran but bad style. Use `if (a==0)` instead.
  - `if (a == 0)` is even better, since `==` is a language token.
- Use space around the equal sign in variable assignment, but not when using a named optional argument
  - `a = b` is strongly preferred over `a=b`
    - One exception is loop indices, where `do i=is,ie` is acceptable
  - `call fn(arg1, arg_name=a)` is strongly preferred over `call fn(arg1, arg_name = a)`
- Use a space after the comma separating subroutine or function arguments in calls and declarations (e.g., `call fn(arg1, arg2, arg3)`), but usually not after the comma separating array indices for multidimensional arrays (e.g., `A(i,j,k) = B(i,j,k)`).

## Line length

Some compilers handle very long lines of code gracefully, but MOM6 needs to adhere to the Fortran standard, which is 132 characters for code, after any macro expansion.  MOM6 does use macros for some memory declarations,
so we need to build in some added space in setting MOM6 guidelines:
- The target maximum length for MOM6 code lines is 120 characters, including comments.
  - 80 character lines can be much easier to read if printed; smaller lines are encouraged where they make sense

## Local variables

- Local variable declarations appear after all the dummy argument declarations, and in the case of a function the return value. We often use `  ! Local variables` to delineate between the argument and local variable declarations.
- Local variables should preferably be descriptive multi-character names meaningful in their context, e.g. `del_rho_int` (delta rho at interface).
 - If using a highly abbreviated or short name, the declaration **MUST** be commented.
 - Units should be provided in the comments describing real variables.
 - Multi-word names should use [snake_case](https://en.wikipedia.org/wiki/Snake_case) (e.g. `delta_rho`).
   - snake_case admittedly used more characters than [camelCase](https://en.wikipedia.org/wiki/CamelCase) but unfortunately Doxygen interprets the Fortran standard too literally and throws away any attempts to use CamelCase. We briefly adopted CamelCase for new code but are systematically replacing it as we Doxygen-ize existing code.

## Block constructs

- `do` and `if` constructs should be terminated with the combined `enddo` and `endif` statements, respectively.
- All other block end statements separate the `end` token (e.g. `end program [label]`)
  - Examples: `program`, `module`, `type`, `subroutine`, `function`, `interface`, `select`

## Loop indices

### Soft case convention
- `i`,`j`,`k` are used for cell-center, layer-center references, e.g. `h(i,j,k)`, `T(i+1,j,k)`.
- `I`,`J` are used for staggered, cell-edge references, e.g. `u(I,j,k)`, `v(i,J,k)`, `q(I,J,k)`, `u(I-1,j,k)`. We use a north-east staggering convention so the `I` means i+1/2 and `I-1` means i-1/2.
- `K` is used for the interface (between layer) references, e.g. `del_t(i,j,K) = T(i,j,K+1) - T(i,j,K)`. The vertical staggering is such that interface `K=1` is above layer `k=1` so that `K` means k-1/2 and `K+1` means k+1/2.

## Global / module data

- Absolutely **NO**!
- There are a few exceptions which are strictly for debugging non-shared memory applications. Do not use these as an excuse for adding module data.

## Module use statements
- Modules may use interfaces, data-types, and constant parameters from other modules via module use statements
  - Modules may not use variables from other modules via use statements
  - All MOM variables are passed around as explicit arguments in interfaces.
- All module use statements must include the `, only` modifier

## Implicit variables
- Absolutely **NO**!
- All MOM6 modules must declare `implicit none ; private`
  - Top-level drivers (i.e., files declaring a program main()) only need `implicit none`

## Array syntax

- We **do not permit scalar-style expressions without the colon notation**, e.g.
  - `tv%S = 0.` is forbidden.
- We do allow array syntax for whole array initialization, e.g.
  - `tv%S(:,:,:) = 0.`
- We do allow array syntax for identical copies, e.g.
  - `S_tmp(:,:,:) = tv%S(:,:,:)`
- We do not allow whole array-syntax for math expressions that include halos  because halos are not guaranteed to have valid data:
  - `tmp(:,:) = 1.0 / G%areaT(:,:)` might have zeros in the halo region.
  - `call post_data(id_AT, G%areaT(:,:)*tv%(T(:,:,1))` is wrong because it can use uninitialized data in halos.

## Data flow

- All needed data is passed via arguments to subroutines and functions, or as the returned value of a function.
- All arguments must have declared intent, with the exception of pointers: i.e. `intent(in)`, `intent(out)`, `intent(inout)`.
- Opaque types are preferred, i.e. referencing members of types defined in other modules is discouraged.

## Documentation in code

- Do it when you are writing the code in the first place!
- All subroutines, functions, arguments, and elements of public types should be described in with Doxygen comments.
- All real variables should have a full physical description, including units.
- All comments should be clearly written and grammatically correct; American spelling is preferred.

## Optimization for divides

Divisions are prone to NaNs and relatively expensive. An optimizing compiler will often rearrange math which makes debugging divisions by zero harder to catch.
- Many common reciprocals are pre-computed
  - Use `Q(i,j) * G%IareaT(i,j)` instead of `Q(i,j) / G%areaT(i,j)`.
- Never write `B / C * D` which is ambiguous to humans (not the compiler)
  - Use `( B * D ) / C`
- Never double divide: `A / ( A + B / C)`
  - Use `( A * C ) / ( A * C + B)`


## Arithmetic reproducibility

Floating point operations are sensitive to the order of operations (associativity), which can not generally be guaranteed due to compiler serialization and optimization.


### Addition

Addition operations must be done in pairs.  When more than one addition is required, the order should be specified using parentheses.

- This is bad:
  - `z = a + b + c`
- This is good:
  - `z = (a + b) + c`

Ideally, the order of operation should be chosen to give the best accuracy.  For example, if `a = 1.` `b = -1.` and `c = 1.e-20`, then the order should be chosen to preserve the residual value.

- This is bad:
  - `a + (b + c) == 0.`
- This is good:
  - `(a + b) + c == 1.e-20`

Not only does this impact reproducibility, but the second choice is more accurate and avoids a potential division by zero.

All operations should be ordered, but no particular ordering is enforced.  Contributors are encouraged to consider the most accurate order of operations.  In some cases the order of sums can be chosen to give expressions that yield identical answers if the underlying horizontal coordinate is rotated by 90 or 180 degrees, which would define the preferred order of operations.


### `sum()` intrinsic

We avoid the Fortran `sum()` intrinsic since the result is dependent on the order of operations within the summation. Using explicit loops allows us to define the order of summation. So
```
a = sum(b(:))
```
should be 
```
a = 0.
do k = 1, nz
  a = a + b(k)
enddo
```
The `prod()` and `matmul()` intrinsics should also not be used.


### Global summation

Floating point operations across MPI ranks are volatile, since the order can change depending on the state of the network.  Functions such as `MPI_Reduce` will not generally be reproducible when used for floating point arithmetic.

When performing summations over MPI ranks, use the `reproducing_sum` function.
```
use MOM_coms, only: reproducing_sum
...
sum = reproducing_sum(array(:,:))
```


### Multiplication

Multiplication is also non-associative and thus not reproducible, but the impact is typically small.  Results may depend on the order of operations, most often in the least significant bit of the fractional component.

In single precision, if `a = b = 1 + 2**-23` and `c = 1.5`, then the following calculations differ:

- `(a*b)*c = 1.50000036`
- `a*(b*c) = 1.50000048`

Parentheses in multiplication operations are currently not enforced, but contributors should consider using them when applicable.


### Transcendental functions

Use of transcendental functions, such as trigonometric functions, non-integer powers, and logarithms, are often implementation-dependent and should be avoided when possible.


### Exponent operator

The exponent operator, `a**b` should be used sparingly, since compilers will often internally replace it with `pow(a, b)`, which is often computed as a transcendental function, `exp(b * log(a))`.  Even small integral power, such as `a**3`, have been known to be replaced with `pow(a, 3)`.  To maximize reproducibility, integral powers should be explicitly computed, e.g. `a3 = a * a * a`.

Square roots (`a**0.5`) should always use the `sqrt()` intrinsic.  An IEEE-754 compliant `sqrt` function must be exactly rounded.

Cube roots (`a**(1./3.)`) should be avoided, the MOM6 intrinsic `cuberoot` should be used.  This is not exactly rounded, but it is reproducible.


## Module structure

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
- Most modules have `_init` and `_end` subroutines for lifecycle management
- `! Local variables` comment separates dummy arguments from local declarations
- **Prefer `allocatable` over `pointer`** for control structure members


## Naming conventions

- **Files**: one module per file; module name must match file name (e.g., `MOM_something.F90` contains module `MOM_something`)
- **Control structures**: `module_CS` (e.g., `energetic_PBL_CS`), always `public` but opaque (with `private` contents)
- **Diagnostic IDs**: `id_diag_name`, initialized to `-1`
- **Inverses**: prefix with `I` (e.g., `IdxCu` = `1/dxCu`, `IareaT` = `1/areaT`)
- **Grid objects**: `G` (ocean_grid_type), `GV` (verticalGrid_type), `US` (unit_scale_type)
- **Public functions**: self-documenting names; private helpers may use short names


## Memory macros

Array dimensions use preprocessor macros from `MOM_memory.h`:
- `SZI_(G)`, `SZJ_(G)`, `SZK_(GV)` for explicit-shape cell-center arrays
- `SZIB_(G)`, `SZJB_(G)`, `SZKB_(GV)` for explicit-shape face/edge-point arrays
- `NIMEM_`, `NJMEM_`, `NKMEM_` for allocatable arrays


## Unit documentation

MOM6 uses a dimensional annotation system for every real variable. Units are documented in square brackets at the end of comments, using a two-part notation:

```
[rescaled_dimensions ~> MKS_equivalent]
```

### Dimensional symbols

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

### Unit annotation rules

1. **Every real variable** must have units in `[brackets]` at the end of its comment
2. **Canonical symbol ordering**: consistent order (e.g., `H L2` not `L2 H`)
3. **Boussinesq variants first**: `[H ~> m or kg m-2]` when units differ by approximation
4. **Simplified expressions only**: write `[T2 Z-1 ~> s2 m-1]`, not `[H Z T-1 / H Z2 T-3 = T2 Z-1 ~> s2 m-1]`
5. **Exponent notation**: `m-1`, `s-2`, `kg-3` (no slashes like `1/m`)
6. **No extra spaces** inside brackets
7. **Nondimensional**: use `[nondim]`
8. **Arbitrary/generic**: use `[A]` or `[A ~> a]`, never `[arbitrary]`
9. **Scaling factors**: `[target source-1 ~> 1]`, e.g., `[Z m-1 ~> 1]`


## Doxygen documentation

### Comment syntax

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

