config_src/external
===================

Subdirectories in here provide null versions of external packages that
can be called by, or used with, MOM6 but that are not needed in all
configurations/executables.

The APIs in these modules should be consistent with the actual external
package. To build with the actual external package include it in the
search path for your build system and remove the associated null version.
