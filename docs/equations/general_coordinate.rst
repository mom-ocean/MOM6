.. general-coordinate-eqns:

General coordinate equations
============================

Transforming to a vertical coordinate :math:`r(z,x,y,t)` ...

The Boussinesq hydrostatic equations of motion in general-coordinate
:math:`r` are

.. math::
  D_t \boldsymbol{u} + f \widehat{\boldsymbol{k}} \wedge \boldsymbol{u} + \frac{\rho}{\rho_o}\boldsymbol{\nabla}_z \Phi + \frac{1}{\rho_o} \boldsymbol{\nabla}_z p &= \boldsymbol{\nabla} \cdotp \boldsymbol{\underline{\tau}} ,\\
  \rho \partial_z \Phi + \partial_z p &= 0 ,\\
  \partial_t z_r + \boldsymbol{\nabla}_r \cdotp ( z_r \boldsymbol{u} ) + \partial_r ( z_r \dot{r} ) &= 0 ,\\
  \partial_t (z_r \theta) + \boldsymbol{\nabla}_r \cdotp ( z_r \boldsymbol{u} \theta ) + \partial_r ( z_r \dot{r} \theta ) &= \boldsymbol{\nabla} \cdotp \boldsymbol{Q}_\theta ,\\
  \partial_t (z_r S) + \boldsymbol{\nabla}_r \cdotp ( z_r \boldsymbol{u} S ) + \partial_r ( z_r \dot{r} S ) &= \boldsymbol{\nabla} \cdotp \boldsymbol{Q}_S ,\\
  \rho &= \rho(S, \theta, z) .
