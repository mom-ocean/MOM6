.. _CVMix_KPP:

The K-Profile Parameterization
==============================

The K-Profile Parameterization (KPP) of :cite:`large1994` is implemented via the Community Vertical Mixing package, `CVMix <http://cvmix.github.io/>`__, which is called directly by this module.

The formulation and implementation of KPP is described in great detail in the `CVMix manual <https://github.com/CVMix/CVMix-description/raw/master/cvmix.pdf>`__ (written by our own Steve Griffies).

.. _section_KPP_nutshell:

KPP in a nutshell
-----------------

Large et al., :cite:`large1994`, decompose the parameterized boundary layer turbulent flux of a scalar, :math:`s`, as

.. math::

   \overline{w^\prime s^\prime} = -K \partial_z s + K \gamma_s(\sigma),

where :math:`\sigma = -z/h` is a non-dimensional coordinate within the boundary layer of depth :math:`h`. :math:`K` is the eddy diffusivity and is a function of position within the boundary layer as well as a function of the surface forcing:

.. math::

   K = h w_s(\sigma) G(\sigma) .

Here, :math:`w_s` is the vertical velocity scale of the boundary layer turbulence and :math:`G(\sigma)` is a “shape function” which is described later. The last term is the “non-local transport” which involves a function :math:`\gamma_s(\sigma)` that is matched to the forcing but is not actually needed in the final implementation. Instead, the entire non-local transport term can be equivalently written

.. math::

   K \gamma_s(\sigma) = C_s G(\sigma) Q_s

where :math:`Q_s` is the surface flux of :math:`s` and :math:`C_s` is a constant. The vertical structure of the redistribution (non-local) term is solely due to the shape function, :math:`G(\sigma)`. In our implementation of KPP, we allow the shape functions used for :math:`K` and for the non-local transport to be chosen independently.

The particular shape function most widely used in the atmospheric community is

.. math::

   G(\sigma) = \sigma (1-\sigma)^2

which satisfies the boundary conditions :math:`G(0) = 0`, :math:`G(1) = 0`, :math:`G^\prime(0) = 1`, and :math:`G^\prime(1) = 0`. Large et al, 1994, alter the function so as to match interior diffusivities but we have found that this leads to inconsistencies within the formulation (see google groups thread `Extreme values of non-local transport <https://groups.google.com/forum/#!msg/CVMix-dev/i6rF-eHOtKI/Ti8BeyksrhAJ>`__). Instead, we use either the above form, or even simpler forms that use alternative upper boundary conditions.

The KPP boundary layer depth is a function of the bulk Richardson number, Rib. But to compute Rib, we need the boundary layer depth. To address this circular logic, we compute Rib for each vertical cell in a column, assuming the BL depth equals to the depth of the given grid cell. Once we have a vertical array of Rib(k), we then call the OBLdepth routine from CVMix to compute the actual OBLdepth. We optionally then “correct” the OBLdepth by cycling through once more, this time knowing the OBLdepth from the first pass. This “correction” step is not used by NCAR. It has been found in idealized MOM6 tests to not be necessary.

See also: kpp_calculate(), kpp_applynonlocaltransport()
