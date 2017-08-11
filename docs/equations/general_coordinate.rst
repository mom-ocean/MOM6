.. general-coordinate-eqns:

General coordinate equations
============================

Transforming to a vertical coordinate :math:`r(z,x,y,t)` ...

The Boussinesq hydrostatic equations of motion in general-coordinate
:math:`r` are

.. math::
  D_t \vec{u} + f \hat{k} \wedge \vec{u} + \nabla_z \Phi + \frac{1}{\rho_o} \nabla_z p &= \nabla \cdot \vec{\underline{\tau}} \\
  \rho \partial_z \Phi + \partial_z p &= 0 \\
  \partial_t z_r + \nabla_r \cdot ( z_r \vec{u} ) + \partial_r ( z_r \dot{r} ) &= 0 \\
  \partial_t z_r \theta + \nabla_r \cdot ( z_r \vec{u} \theta ) + \partial_r ( z_r \dot{r} \theta ) &= \nabla \cdot \vec{Q}_\theta \\
  \partial_t z_r S + \nabla_r \cdot ( z_r \vec{u} S ) + \partial_r ( z_r \dot{r} S ) &= \nabla \cdot \vec{Q}_S \\
  \rho &= \rho(S, \theta, z)
