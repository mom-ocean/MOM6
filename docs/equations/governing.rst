.. equations-governing:

Governing equations
===================

The Boussinesq hydrostatic equations of motion in height coordinates are

.. math::
  D_t \boldsymbol{u} + f \widehat{\boldsymbol{k}} \wedge \boldsymbol{u} + \frac{\rho}{\rho_o} \boldsymbol{\nabla}_z \Phi + \frac{1}{\rho_o} \boldsymbol{\nabla}_z p &= \boldsymbol{\nabla} \cdotp \boldsymbol{\underline{\tau}} , \\
  \rho \partial_z \Phi + \partial_z p &= 0 , \\
  \boldsymbol{\nabla}_z \cdotp \boldsymbol{u} + \partial_z w &= 0 , \\
  D_t \theta &= \boldsymbol{\nabla} \cdotp \boldsymbol{Q}_\theta , \\
  D_t S &= \boldsymbol{\nabla} \cdotp \boldsymbol{Q}_S , \\
  \rho &= \rho(S, \theta, z) ,

where notation is described in :ref:`equations-notation`. :math:`\boldsymbol{\underline{\tau}}` is the stress tensori and
:math:`\boldsymbol{Q}_\theta` and :math:`\boldsymbol{Q}_S` are fluxes of heat and salt respectively.


.. :ref:`vector_invariant`
The total derivative is

.. math::
  D_t & \equiv \partial_t + \boldsymbol{v} \cdotp \boldsymbol{\nabla} \\
      &= \partial_t + \boldsymbol{u} \cdotp \boldsymbol{\nabla}_z + w \partial_z .

The non-divergence of flow allows a total derivative to be re-written in flux form:

.. math::
  D_t \theta &= \partial_t + \boldsymbol{\nabla} \cdotp ( \boldsymbol{v} \theta ) \\
             &= \partial_t + \boldsymbol{\nabla}_z \cdotp ( \boldsymbol{u} \theta ) + \partial_z ( w \theta ) .

The above equations of motion can thus be written as:

.. math::
  D_t \boldsymbol{u} + f \widehat{\boldsymbol{k}} \wedge \boldsymbol{u} + \frac{\rho}{\rho_o}\boldsymbol{\nabla}_z \Phi + \frac{1}{\rho_o} \boldsymbol{\nabla}_z p &= \boldsymbol{\nabla} \cdotp \boldsymbol{\underline{\tau}} ,\\
  \rho \partial_z \Phi + \partial_z p &= 0 ,\\
  \boldsymbol{\nabla}_z \cdotp \boldsymbol{u} + \partial_z w &= 0 ,\\
  \partial_t \theta + \boldsymbol{\nabla}_z \cdotp ( \boldsymbol{u} \theta ) + \partial_z ( w \theta ) &= \boldsymbol{\nabla} \cdotp \boldsymbol{Q}_\theta ,\\
  \partial_t S + \boldsymbol{\nabla}_z \cdotp ( \boldsymbol{u} S ) + \partial_z ( w S ) &= \nabla \cdotp \boldsymbol{Q}_S ,\\
  \rho &= \rho(S, \theta, z) .

.. toctree::
  vector_invariant_eqns
