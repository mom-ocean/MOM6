.. _Tracer_Timestep:

Tracer Timestep
===============

Overview of Tracer Timestep

The MOM6 code handles advection and lateral diffusion of all tracers. For potential temperature and salinity, it also timesteps the thermodynamics and vertical mixing (column physics). Since evaporation and precipitation are handled as volume changes, the layer thicknesses need to be updated:

.. math::

   \frac{\partial h_k}{\partial t} = (P - E)_k

The full tracer equation for tracer :math:`\theta` is:

.. math::

   \frac{\partial}{\partial t} (h_k\theta_k) + \nabla_s \cdot
      (\vec{u}h_k \theta_k) = Q_k^\theta h_k + \frac{1}{h_k} \Delta \left(
      \kappa \frac{\partial \theta}{\partial z} \right) + \frac{1}{h_k}
      \nabla_s (h_k K \nabla_s \theta)

Here, the advection is on the left hand side of the equation while the right hand side contains thermodynamic processes, vertical diffusion, and horizontal diffusion. There is more than one choice for vertical diffusion; these will be described elsewhere. Also, the lateral diffusion is handled in novel ways so as to avoid introduction of new extrema and to avoid instabilities associated with rotated mixing tensors. The lateral diffusion is described in `Horizontal Diffusion <Horizontal_Diffusion.html#horizontal-diffusion>`__.
