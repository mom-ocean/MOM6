.. _equations-notation:

Notation for equations
======================

Symbols for variables
---------------------

:math:`z` refers to elevation (or height), increasing upward so that for much of the ocean :math:`z` is negative.

:math:`x` and :math:`y` are the Cartesian horizontal coordinates.

:math:`\lambda` and :math:`\phi` are the geographic coordinates on a sphere (longitude and latitude respectively).

Horizontal components of velocity are indicated by :math:`u` and :math:`v` and vertical component by :math:`w`.

:math:`p` is pressure and :math:`\Phi` is geo-potential:

.. math::
  \Phi = g z .

The thermodynamic state variables are usually salinity, :math:`S`, and potential temperature, :math:`\theta` or the absolute salinity and conservative temperature, depending on the equation of state. :math:`\rho` is in-situ density.

Vector notation
---------------

The three-dimensional velocity vector is denoted :math:`\boldsymbol{v}`

.. math::
  \boldsymbol{v} = \boldsymbol{u} + \widehat{\boldsymbol{k}} w ,

where :math:`\widehat{\boldsymbol{k}}` is the unit vector pointed in the upward vertical direction and :math:`\boldsymbol{u} = (u, v, 0)` is the horizontal
component of velocity normal to the vertical.

The gradient operator without a suffix is three dimensional:

.. math::
  \boldsymbol{\nabla} = ( \boldsymbol{\nabla}_z, \partial_z  ) .

but a suffix indicates a lateral gradient along a surface of constant property indicated by the suffix:

.. math::
  \boldsymbol{\nabla}_z = \left( \left. \partial_x \right|_z, \left. \partial_y \right|_z, 0 \right) .
