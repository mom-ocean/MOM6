.. _Tracer_Advection:

Tracer Advection
================

MOM6 implements a generalised tracer advection scheme, which is a combination of the modified flux advection scheme :cite:`easter1993` with reconstructed tracer distributions. The tracer distributions may be piecewise linear (PLM) or piecewise parabolic (PPM), which may itself use either the :cite:`colella1984` (CW84) or :cite:`huynh1997` (H3) reconstruction.

.. _Flux_advection:

Flux advection
--------------

The modified flux advection scheme preserves the tracer mixing ratio in a cell across directional splitting by accounting for changes in mass changes. Fluxes are applied to alternate directions in turn, restricting the applied flux so as not to evacuate all mass out of a cell. Because of this, we need to know the stencil used during the calculation of the reconstruction. Every iteration of the splitting algorithm, cells at the edge of a processor’s data domain are invalidated. When this invalidation region extends below the halo, a group pass is required to refresh the halo. A larger stencil (such as for the CW84 reconstruction) therefore introduces more frequent updates, and may impact performance.

.. _Tracer_reconstruction:

Tracer reconstruction
---------------------

While MOM6 only carries the mean tracer concentration in a cell, a higher order reconstruction is computed for the purpose of advection. Reconstructions are also modified to ensure that monotonicity is preserved (i.e. spurious minima or maxima cannot be introduced).

The piecewise linear (PLM) reconstruction uses the monotonic modified van Leer scheme :cite:`lin1994`. One might think to use the average of the one-sided differences of mean tracer concentration within a cell to calculate the slope of the linear reconstruction, however this method guarantees neither monotonicity, nor positive definiteness. Instead, the method is locally limited to the minimum of this average slope and each of the one-sided slopes, i.e. 

.. math::

   \Delta \Phi_i = \min\left\{\left|[\Delta
   \Phi_i]_\text{avg}\right|, 2\left(\Phi_i - \Phi_i^\text{min}\right),
   2\left(\Phi_i^\text{max} - \Phi_i\right)\right\}

(where :math:`\Phi_i^\text{min}` is the minimum in the 3-point stencil).

In a PPM scheme (`PPM Advection Scheme <PPM.html#ppm>`__), for a cell with mean tracer concentration :math:`\Phi_i`, the values at the left and right interfaces, :math:`\Phi_{L,i}` and :math:`\Phi_{R,i}` must be estimated. First, an interpolation is used to calculate :math:`\Phi_{i-1/2}` and :math:`\Phi_{i+1/2}`. These values are then modified to preserve monotonicity in each cell, which introduces discontinuities between cell edges (e.g. :math:`\Phi_{R,i}` and :math:`\Phi_{L,i+1}`).

The reconstruction :math:`\Phi_i(\xi)` then satisfies three properties:

- total amount of tracer is conserved, :math:`\int_{\xi_{i-1/2}}^{\xi_{i+1/2}} \Phi_i(\xi') \,\mathrm d\xi' = \Phi_i`

- left interface value matches, :math:`\Phi(\xi_{i-1/2}) = \Phi_{L,i}`

- right interface value matches, :math:`\Phi(\xi_{i+1/2}) = \Phi_{R,i}`

There are two methods of reconstruction for a piecewise parabolic (PPM) profile. They differ in the estimate of interface values :math:`\Phi_{i+1/2}` prior to monotonicity limiting. The CW84 scheme makes use of the limited slope :math:`\Delta\Phi_i` from PLM, above. This has the effect of requiring a larger stencil for each reconstruction. On the other hand, the H3 scheme reduces the requirement of this stencil, by only examining the tracer concentrations in adjacent cells, at the same time reducing order of accuracy of the reconstruction.
