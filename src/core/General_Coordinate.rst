.. _General_Coordinate:

Generalized vertical coordinate equations
=========================================

The ocean equations discretized by MOM6 are formulated using generalized vertical coordinates. Motivation for using generalized vertical coordinates, and a full accounting of the ocean equations written using these coordinates, can be found in Griffies, Adcroft and Hallberg (2020) :cite:`Griffies_Adcroft_Hallberg2020`. Here we provide a brief summary.

Consider a smooth function of space and time, :math:`r(x,y,z,t)`, that has a single-signed and non-zero vertical derivative known as the specific thickness

.. math::

   \begin{aligned}
   \partial z/\partial r = (\partial r/\partial z)^{-1} =  \mbox{specific thickness.}
   \end{aligned}

The specific thickness measures the inverse vertical stratification of the vertical coordinate surfaces. As so constrained, :math:`r` can uniquely prescribe a positiion in the vertical. Consequently, the ocean equations can be mapped one-to-one from geopotential vertical coordinates to generalized vertical coordinate. Upon transforming to :math:`r`-coordinates, the material time derivative of :math:`r` appears throughout the equations, playing the role of a pseudo-vertical velocity, and we make use of the following shorthand for this derivative

.. math::

   \begin{aligned}
   \dot{r} = D_{t} r.
   \end{aligned}

The Boussinesq hydrostatic ocean equations take the following form using generalized vertical coordinates (:math:`r`-coordinates)

.. _html:r-equations:

.. _eq:r-horz-momentum:

.. _eq:r-hydrostatic-equation:

.. _eq:r-non-divergence:

.. _eq:r-temperature-equation:

.. _eq:r-salinity-equation:

.. math::

   \begin{aligned}
   \\
   \rho_o \left[
    \partial_{t} \mathbf{u} + (f + \zeta) \, \hat{\mathbf{z}} \times \mathbf{u}
    + \dot{r} \, \partial_{r} \mathbf{u} \right]
    &= -\nabla_r \, (p + \rho_{o} \, K) -\rho \nabla_r \Phi + \rho_{o} \, \mathbf{\mathcal{F}}
    &\mbox{horizontal momentum}
   \\
   \rho \, \partial_{r} \Phi + \partial_{r}p
    &= 0
   &\mbox{hydrostatic}
   \\
    \partial_{t}( z_r)
   + \nabla_r \cdot ( z_r \, \mathbf{u} )
   + \partial_{r} ( z_r \, \dot{r} )
   &= 0
   &\mbox{specific thickness}
   \\
    \partial_{t} ( \theta \, z_r )
   + \nabla_r \cdot ( \theta z_r \, \mathbf{u} )
   + \partial_{r} ( \theta \, z_r \, \dot{r} )
   &= z_r \mathbf{\mathcal{N}}_\theta^\gamma
   - \partial_{r} J_\theta^{(z)}
   &\mbox{potential/Conservative temp}
   \\
   \partial_{t} ( S \, z_r)
   + \nabla_r \cdot ( S \, z_r \, \mathbf{u} )
   + \partial_{r} ( S \, z_r \, \dot{r} )
   &= z_r \mathbf{\mathcal{N}}_S^\gamma
   - \partial_{r} J_S^{(z)}
   &\mbox{salinity}
   \\
   \rho &= \rho( S, \theta, -g \rho_0 z )
   &\mbox{equation of state.}
   \end{aligned}

The time derivatives appearing in these equations are computed with the generalized vertical coordinate fixed rather than the geopotential. It is a common misconception that the horizontal velocity, :math:`\mathbf{u}`, is rotated to align with constant :math:`r` surfaces. Such is not the case. Rather, the horizontal velocity, :math:`\mathbf{u}`, is precisely the same horizontal velocity used with geopotential coordinates. However, its evolution has here been formulated using generalized vertical coordinates.

As a finite volume model, MOM6 is discretized in the vertical by integrating between surfaces of constant :math:`r`. The layer thickness is a basic term appearing in these equations, which results from integrating the specific thickness over a layer

.. math::

   \begin{aligned}
   h = \int z_r \, \mathrm{d}r.
   \end{aligned}

Correspondingly, the model variables are treated as finite volume averages over each layer, with full accounting of this finite volume approach presented in Griffies, Adcroft and Hallberg (2020) :cite:`Griffies_Adcroft_Hallberg2020`, and with the semi-discrete model ocean model equations written as follows.

.. _eq:h-horz-momentum:

.. _eq:h-hydrostatic-equation:

.. _eq:h-thickness-equation:

.. _eq:h-temperature-equation:

.. _eq:h-salinity-equation:

.. _eq:h-equation-of-state:

.. math::

   \begin{aligned}
   \rho_0
   \left[ \frac{\partial \mathbf{u}}{\partial t} + \frac{( f + \zeta )}{h} \,
   \hat{\mathbf{z}} \times h \, \mathbf{u} + \underbrace{ \dot{r} \,
   \frac{\partial \mathbf{u}}{\partial r} }
   \right]
   &= -\nabla_r \, (p + \rho_{0} \, K) -
   \rho \nabla_r \, \Phi + \mathbf{\mathcal{F}}
   &\mbox{horizontal momentum}
   \\
   \rho \, \delta_r \Phi + \delta_r p
   &= 0
   &\mbox{hydrostatic}
   \\
   \frac{\partial h}{\partial t} + \nabla_r \cdot \left( h \, \mathbf{u} \right) +
   \underbrace{ \delta_r ( z_r \dot{r} ) }
    &= 0
   &\mbox{thickness}
   \\
   \frac{\partial ( \theta \, h )}{\partial t} + \nabla_r \cdot \left( \theta h \,
   \mathbf{u} \right) + \underbrace{ \delta_r ( \theta \, z_r \dot{r} ) }
   &=
   h \mathbf{\mathcal{N}}_\theta^\gamma - \delta_r J_\theta^{(z)}
   &\mbox{potential/Conservative temp}
   \\
   \frac{\partial ( S \, h )}{\partial t} + \nabla_r \cdot \left( S \, h \,
   \mathbf{u} \right) + \underbrace{ \delta_r ( S \, z_r \dot{r} ) }
   &=
   h \mathbf{\mathcal{N}}_S^\gamma - \delta_r J_S^{(z)}
   &\mbox{salinity}
   \\
   \rho &= \rho\left( S, \theta, -g \rho_0 z(r) \right)
   &\mbox{equation of state,}
   \end{aligned}

where

.. math::

   \begin{aligned}
   \delta_{r} = \mathrm{d}r \, (\partial/\partial r)
   \end{aligned}

is the discrete vertical difference operator. The pressure gradient accelerations in the momentum equation are written in continuous-in-the-vertical form for brevity; the exact discretization is detailed in :cite:`adcroft2008` and :cite:`Griffies_Adcroft_Hallberg2020`. The :math:`1/h` and :math:`h` appearing in the horizontal momentum equation are carefully handled in the code to ensure proper cancellation even when the layer thickness goes to zero i.e., l’Hospital’s rule is respected.

The MOM6 time-stepping algorithm integrates the above layer-averaged equations forward in time allowing the vertical grid to follow the motion, i.e. :math:`\dot{r}=0`, so that the underbraced terms are dropped. This approach is generally known as a Lagrangian method, with the Lagrangian approach in MOM6 limited to the vertical direction. After each Lagrangian step, a regrid step is applied that generates a new vertical grid of the user’s choosing. The ocean state is then remapped from the old to the new grid. The physical state is not meant to change during the remap step, yet truncation errors make remapping imperfect. We employ high-order accurate reconstructions to minimize errors introduced during the remap step (:cite:`white2008`, :cite:`white2009`). The connection between time-stepping and remapping is described in section `Vertical Lagrangian method in pictures <ALE_Timestep.html#ale-timestep>`__.
