.. _Barotropic_Momentum_Equations:

Barotropic Momentum Equations
=============================

Barotropic Momentum Equations

The barotropic equations are timestepped on a relatively short timestep compared to the rest of the model. Since the timestep constraints for this are known, the barotropic timestep is computed at runtime.

The 2-d linear momentum equations with integrated continuity are:

.. math::

   \frac{\partial \eta}{\partial t} +  \nabla \cdot \left( ( D + \eta) \vec{u}_{BT}
     h_k \right) = P - E

.. math::

   \frac{\partial \vec{u}_{BT}}{\partial t} = - g \nabla \eta - f \hat{z} \times
     \vec{u}_{BT} + \vec{F}_{BT}

where

.. math::

   \vec{u}_{BT} \equiv \frac{1}{D + \eta} \int_{-D}^\eta \vec{u}dz

and :math:`\vec{F}_{BT}` is the barotropic momentum forcing from baroclinic processes. Note that explicit mass fluxes such as evaporation and precipitation change the model volume explicitly.

In the mode splitting between baroclinic and barotropic processes, it is important to include the contribution of free surface waves on the internal interface heights on the pressure gradient force, shown here as :math:`g_{Eff}`:

.. math::

   \frac{\partial p}{\partial z} = -\rho g

.. math::

   g_{Eff} = g + \frac{\partial}{\partial \eta} \left[ \frac{1}{D + \eta}
      \int_{-D}^\eta p dz \right]

The barotropic momentum equation then becomes:

.. math::

   \frac{\partial \vec{u}_{BT}}{\partial t} + f \hat{z} \times
     \vec{u}_{BT} + \frac{1}{\rho_0} \nabla g_{Eff} \eta = \mbox{Residual}

Without including the internal wave motion in the barotropic equations, one can generate instabilities (:cite:`bleck1990`, :cite:`hallberg1997a`).
