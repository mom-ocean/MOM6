.. _Domain_Decomposition:

Domain Decomposition
====================

.. _section_domain_decomp:

Domain Decomposition
--------------------

MOM6 supports serial, OpenMP, and MPI computations, with the user choosing between them at run time. All are accomplished through domain decomposition in the horizontal. All of the horizontal operations are explicit with a relatively small footprint, so the tiling is a logical choice. Some goals in the parallel design of MOM6 were:

- Don’t hard-code the number of processes.

- MPI and OpenMP share the same basic structure.

- Don’t break the serial optimizations.

- Same result as serial code for any number of processes.

- Portability - able to run on any (Unix) system.

The whole horizontal MOM6 grid is shown in `Declaration of variables <Horizontal_Indexing.html#section-memory>`__. The computations are done over the cells inside the darker line; the cells are numbered 1 to NIGLOBAL in the :math:`x`-direction and 1 to NJGLOBAL in the :math:`y`-direction. Those looking ahead to running in parallel would be wise to include factors of two and three in their choice of NIGLOBAL and NJGLOBAL when building new grids. MOM6 will run in parallel with any values of these, but the computations might not be load-balanced.

.. _section_wide_halos:

Wide Halos
----------
