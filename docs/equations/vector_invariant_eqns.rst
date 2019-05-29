.. vector-invariant-eqn:

Vector Invariant Equations
==========================

MOM6 solves the momentum equations written in vector-invariant form.

A vector identity allows the total derivative of velocity to be written in the vector-invariant form:

.. math::
  D_t \boldsymbol{u} &= \partial_t \boldsymbol{u} + \boldsymbol{v} \cdotp \boldsymbol{\nabla} \boldsymbol{u} \\
              &= \partial_t \boldsymbol{u} + \boldsymbol{u} \cdotp \boldsymbol{\nabla}_z \boldsymbol{u} + w \partial_z \boldsymbol{u} \\
              &= \partial_t \boldsymbol{u} + \left( \boldsymbol{\nabla} \wedge \boldsymbol{u} \right) \wedge \boldsymbol{v} + \boldsymbol{\nabla} \underbrace{\frac{1}{2} \left|\boldsymbol{u}\right|^2}_{\equiv K} .

The flux-form equations of motion in height coordinates can thus be written succinctly as:

.. math::
  \partial_t \boldsymbol{u} + \left( f \widehat{\boldsymbol{k}} + \boldsymbol{\nabla} \wedge \boldsymbol{u} \right) \wedge \boldsymbol{v} + \boldsymbol{\nabla} K
  + \frac{\rho}{\rho_o} \boldsymbol{\nabla} \Phi + \frac{1}{\rho_o} \boldsymbol{\nabla} p &= \boldsymbol{\nabla} \cdotp \boldsymbol{\underline{\tau}} ,\\
  \boldsymbol{\nabla}_z \cdotp \boldsymbol{u} + \partial_z w &= 0 ,\\
  \partial_t \theta + \boldsymbol{\nabla}_z \cdotp ( \boldsymbol{u} \theta ) + \partial_z ( w \theta ) &= \boldsymbol{\nabla} \cdotp \boldsymbol{Q}_\theta ,\\
  \partial_t S + \boldsymbol{\nabla}_z \cdotp ( \boldsymbol{u} S ) + \partial_z ( w S ) &= \boldsymbol{\nabla} \cdotp \boldsymbol{Q}_S ,\\
  \rho &= \rho(S, \theta, z) ,

where the horizontal momentum equations and vertical hydrostatic balance equation have been written as a single three-dimensional equation.
