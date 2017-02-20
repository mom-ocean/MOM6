.. equations-governing:

Governing equations
===================

The Boussinesq hydrostatic equations of motion in height coordinates are

.. math::
  D_t \vec{u} + f \hat{k} \wedge \vec{u} + \nabla_z \Phi + \frac{1}{\rho_o} \nabla_z p &= \nabla \cdot \vec{\underline{\tau}} \\
  \rho \partial_z \Phi + \partial_z p &= 0 \\
  \nabla_z \cdot \vec{u} + \partial_z w &= 0 \\
  D_t \theta &= \nabla \cdot \vec{Q}_\theta \\
  D_t S &= \nabla \cdot \vec{Q}_S \\
  \rho &= \rho(S, \theta, z)
 
where notation is described in :ref:`equations-notation`. :math:`\vec{\underline{\tau}}` is the stress tensori and
:math:`\vec{Q}_\theta` and :math:`\vec{Q}_S` are fluxes of heat and salt respectively.

The total derivative is

.. math::
  D_t &\equiv \partial_t + \vec{v} \cdot \nabla \\
      &= \partial_t + \vec{u} \cdot \nabla_z + w \partial_z

The non-divergence of flow allows a total derivative to be re-written in flux form:

.. math::
  D_t \theta &= \partial_t + \nabla \cdot ( \vec{v} \theta ) \\
             &= \partial_t + \nabla_z \cdot ( \vec{u} \theta ) + \partial_z ( w \theta )

An identity allows the total derivative of velocity to be written in the vector-invariant form:

.. math::
  D_t \vec{u} &= \partial_t \vec{u} + \vec{v} \cdot \nabla \vec{u} \\
              &= \partial_t \vec{u} + \vec{u} \cdot \nabla_z \vec{u} + w \partial_z \vec{u} \\
              &= \partial_t \vec{u} + \left( \nabla \wedge \vec{u} \right) \wedge \vec{v} + \nabla \frac{1}{2} \left|\vec{u}\right|^2

The above equations of motion can thus be written as:

.. math::
  \partial_t \vec{u} + \left( f \hat{k} + \nabla \wedge \vec{u} \right) \wedge \vec{u} + \nabla K
  + \frac{\rho}{\rho_o} \nabla \Phi + \frac{1}{\rho_o} \nabla p &= \nabla \cdot \vec{\underline{\tau}} \\
  \nabla_z \cdot \vec{u} + \partial_z w &= 0 \\
  \partial_t \theta + \nabla_z \cdot ( \vec{u} \theta ) + \partial_z ( w \theta ) &= \nabla \cdot \vec{Q}_\theta \\
  \partial_t S + \nabla_z \cdot ( \vec{u} S ) + \partial_z ( w S ) &= \nabla \cdot \vec{Q}_S \\
  \rho &= \rho(S, \theta, z)

where the horizontal momentum equations and vertical hydrostatic balance equation have been written as a single three-dimensional equation.
