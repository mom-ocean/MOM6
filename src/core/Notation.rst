.. _Notation:

Notation for equations
======================

.. _Symbols:

Symbols for variables
---------------------

:math:`z` refers to geopotential elevation (or height), increasing upward and with :math:`z=0` defining the resting ocean surface. Much of the ocean has :math:`z < 0`.

:math:`x` and :math:`y` are the Cartesian horizontal coordinates. MOM6 uses generalized orthogonal curvilinear horizontal coordinates. However, the equations are simpler to write using Cartesian coordinates, and it is very straightforward to generalize the horizontal coordinates using the methods in Chapters 20 and 21 of :cite:`SMGbook`.

:math:`\lambda` and :math:`\phi` are the geographic coordinates on a sphere (longitude and latitude, respectively).

Horizontal components of velocity are indicated by :math:`u` and :math:`v` and vertical component by :math:`w`.

:math:`p` is the hydrostatic pressure.

:math:`\Phi` is the geopotential. In the absence of tides, the geopotential is given by :math:`\Phi = g z,` whereas more general expressions hold when including astronomical tide forcing.

The thermodynamic state variables can be salinity, :math:`S`, and potential temperature, :math:`\theta`. Alternatively, one can choose the Conservative Temperature if using the TEOS10 equation of state from :cite:`TEOS2010`.

:math:`\rho` is the in-situ density computed as a function :math:`\rho(S,\theta,p)` for non-Boussinesq ocean or :math:`\rho(S,\theta,p=-g \, \rho_o \, z)` for Boussinesq ocean. See Young (2010) :cite:`Young2010` or Section 2.4 of Vallis (2017) :cite:`GVbook` for reasoning behind the simplified pressure used in the Boussinesq equation of state.

.. _vector_notation:

Vector notation
---------------

The three-dimensional velocity vector is denoted :math:`\mathbf{v}` and it is decomposed into its horizontal and vertical components according to

.. math::

   \begin{aligned}
   \mathbf{v}
     = \mathbf{u} + \hat{\mathbf{z}} \, w
     = \hat{\mathbf{x}} \, u + \hat{\mathbf{y}} \, v + \hat{\mathbf{z}} \, w,
   \end{aligned}

where :math:`\hat{\mathbf{z}}` is the unit vector pointed in the upward vertical direction and :math:`\mathbf{u} = (u, v, 0)` is the horizontal component of velocity normal to the vertical.

The three-dimensional gradient operator is denoted :math:`\nabla`, and it is decomposed into its horizontal and vertical components according to

.. math::

   \begin{aligned}
   \nabla
      = \nabla_z + \hat{\mathbf{z}} \, \partial_z
      = \hat{\mathbf{x}} \, \partial_x + \hat{\mathbf{y}} \, \partial_y + \hat{\mathbf{z}} \, \partial_z.
   \end{aligned}
