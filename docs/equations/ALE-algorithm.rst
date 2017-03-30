ALE algorithm
=============

The semi-discrete, vertically integrated, Boussinesq hydrostatic equations of
motion in general-coordinate :math:`r` are

.. math::
  D_t \vec{u} + f \hat{k} \wedge \vec{u} + \nabla_z \Phi + \frac{1}{\rho_o} \nabla_z p &= \nabla \cdot \vec{\underline{\tau}} \\
  \rho \delta_k \Phi + \delta_k p &= 0 \\
  \partial_t h + \nabla_r \cdot ( h \vec{u} ) + \delta_k ( z_r \dot{r} ) &= 0 \\
  \partial_t h \theta + \nabla_r \cdot ( h \vec{u} \theta ) + \delta_k ( z_r \dot{r} \theta ) &= \nabla \cdot \vec{Q}_\theta \\
  \partial_t h S + \nabla_r \cdot ( h \vec{u} S ) + \delta_k ( z_r \dot{r} S ) &= \nabla \cdot \vec{Q}_S \\
  \rho &= \rho(S, \theta, z)

The Arbitrary-Lagrangian-Eulerian algorithm we use is quasi-Lagrangian in
that in the first (Lagrangian) phase, regardless of the current mesh (or coordinate
:math:`r`) we integrate the equations forward with :math:`\dot{r}=0`, i.e.:

.. math::
  D_t \vec{u} + f \hat{k} \wedge \vec{u} + \nabla_z \Phi + \frac{1}{\rho_o} \nabla_z p &= \nabla \cdot \vec{\underline{\tau}} \\
  \rho \delta_k \Phi + \delta_k p &= 0 \\
  \partial_t h + \nabla_r \cdot ( h \vec{u} ) &= 0 \\
  \partial_t h \theta + \nabla_r \cdot ( h \vec{u} \theta ) &= \nabla \cdot \vec{Q}_\theta \\
  \partial_t h S + \nabla_r \cdot ( h \vec{u} S ) &= \nabla \cdot \vec{Q}_S \\
  \rho &= \rho(S, \theta, z)

Notice that by setting :math:`\dot{r}=0` all the terms with the metric
:math:`z_r` disappeared.

After a finite amount of time, the mesh (:math:`h`) may become very distorted
or unrelated to the intended mesh. At any point in time, we can simply define
a new mesh and remap from the current mesh to the new mesh without an
explicit change in the physical state.

