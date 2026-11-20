Source documentation completeness
=================================

``docs-and-style.yml`` checks whitespace and source documentation with
``flint check-docs``, then tests the builder and generates/validates the published site.
These are separate checks: an undocumented variable can still be rendered, but
is a source-policy failure. Neither check rewrites Fortran source.

Local command, from the repository root::

    flint check-docs src config_src/drivers/solo_driver \
        config_src/memory/dynamic_symmetric config_src/external \
        config_src/drivers/FMS_cap \
        -I config_src/memory/dynamic_symmetric \
        --defer-include src/framework/version_variable.h \
        --exclude src/equation_of_state/TEOS10 \
        --exclude pkg/GSW-Fortran --exclude pkg/CVMix-src

Install the Flint revision from ``requirements-api.txt`` first.

Scope
-----

The roots match the former Doxygen check, including both the solo driver and
FMS cap. They intentionally differ from site generation: files are checked
individually, so alternative drivers need not share a module namespace.
TEOS10 and dependency implementations are excluded, including symlink targets.

The default check requires descriptions for parsed concrete types,
named interfaces, module variables, and type components. Private entities are included, as in
the old ``EXTRACT_PRIVATE = YES`` configuration. Ordinary procedure-local
variables are excluded. Explicit ``!>``, ``!<``, and ``!!`` comments and
declaration groups qualify; empty markers and ordinary comments do not.
Descriptions shared by a declaration on one line count for all its entities;
a comment on a continued declaration does not automatically cover other lines.

Module/submodule, program, subroutine, function, and abstract-type descriptions
are also required by the intended
policy, but enforcement is deferred while the existing omissions are addressed.
Add ``--strict`` to include those descriptions. CI currently uses relaxed mode.
This does not skip their contents: variables, types, and interfaces inside
undocumented modules and procedures remain checked. Procedure descriptions are
strict-only regardless of public/private accessibility; the requirement is
deferred, not abandoned.
Abstract-type components remain checked in relaxed mode, as do concrete types
that extend an abstract base. Regression tests verify that ``--strict`` restores
each deferred description category. There are no permanent exclusions for them.

For initial adoption, CI also defers descriptions originating in the specific
header ``src/framework/version_variable.h`` using ``--defer-include``. The
declaration of ``version`` remains subject to the intended documentation policy;
``--strict`` restores this check too. Other headers and local declarations with
macro initializers are not exempt. Included declarations are reported using
their original header path and line, rather than the next line of the caller.

Private module parameters are also deferred to ``--strict``. Both properties
must hold: public parameters and private non-parameter variables are still
checked. The checker honors module default access, named ``public``/``private``
statements, and access attributes on the declaration.

The former ``WARN_NO_PARAMDOC = NO`` setting did not mandate documentation for
every dummy argument. Accordingly, this CI check does not enable Flint's optional
``--arguments`` mode. Specific type-bound procedure and finalization-procedure
descriptions are strict-only, like other procedures. Generic binding descriptions
are checked in relaxed mode, like named interfaces. Function-result descriptions,
units syntax, and prose quality are not checked. This policy is
explicit; it does not claim exact equivalence with Doxygen's heuristics or
coverage of unsupported parser constructs.

Dummy procedures declared with ``procedure(interface)`` or ``procedure()`` are
included in argument checking, with or without ``pointer``. Their declaration
comments and locations are retained. Module-level non-pointer procedure
declarations use the strict-only procedure-description policy; procedure pointers
stored in modules remain subject to the ordinary variable policy.

Diagnostics
-----------

Missing descriptions are reported as ``path:line: missing-documentation: kind
qualified_name``. The command exits 1 for omissions and 2 for invalid inputs or
reported incomplete parses, including missing include files. Other parser
failures also terminate unsuccessfully. There is no suppressed-finding baseline
and no success status for an empty source selection.

Passing relaxed mode does not establish complete documentation. Use ``--strict``
to report the deferred categories as well.
