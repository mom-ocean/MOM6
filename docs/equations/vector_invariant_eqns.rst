.. vector-invariant-eqn:

Vector Invariant Equations
==========================

MOM6 solve the momentum equations written in vector-invariant form.

An identity allows the total derivative of velocity to be written in the vector-invariant form:

.. math::
  D_t \vec{u} &= \partial_t \vec{u} + \vec{v} \cdot \nabla \vec{u} \\
              &= \partial_t \vec{u} + \vec{u} \cdot \nabla_z \vec{u} + w \partial_z \vec{u} \\
              &= \partial_t \vec{u} + \left( \nabla \wedge \vec{u} \right) \wedge \vec{v} + \nabla \frac{1}{2} \left|\vec{u}\right|^2

The flux-form equations of motion in height coordinates can thus be written succinctly as:

.. math::
  \partial_t \vec{u} + \left( f \hat{k} + \nabla \wedge \vec{u} \right) \wedge \vec{v} + \nabla K
  + \frac{\rho}{\rho_o} \nabla \Phi + \frac{1}{\rho_o} \nabla p &= \nabla \cdot \vec{\underline{\tau}} \\
  \nabla_z \cdot \vec{u} + \partial_z w &= 0 \\
  \partial_t \theta + \nabla_z \cdot ( \vec{u} \theta ) + \partial_z ( w \theta ) &= \nabla \cdot \vec{Q}_\theta \\
  \partial_t S + \nabla_z \cdot ( \vec{u} S ) + \partial_z ( w S ) &= \nabla \cdot \vec{Q}_S \\
  \rho &= \rho(S, \theta, z)

where the horizontal momentum equations and vertical hydrostatic balance equation have been written as a single three-dimensional equation.
