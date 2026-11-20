.. _BML:

Bulk Surface Mixed Layer
========================

This bulk surface mixed layer scheme was designed to be used with a purely isopycnal model. Following :cite:`niiler1977`, :cite:`oberhuber1993`, and Hallberg (:cite:`muller2003`) the TKE budget is used to construct a time-evolving homogeneous mixed layer. A buffer layer sits between the mixed layer and the interior ocean to mediate between the two.

The following processes are executed, in the order listed.

1. Undergo convective adjustment into mixed layer.

2. Apply surface heating and cooling.

3. Starting from the top, entrain whatever fluid the TKE budget permits. Penetrating shortwave radiation is also applied at this point.

4. If there is any unentrained fluid that was formerly in the mixed layer, detrain this fluid into the buffer layer. This is equivalent to the mixed layer detraining to the Monin- Obukhov depth.

5. Divide the fluid in the mixed layer evenly into CS%nkml pieces.

6. Split the buffer layer if appropriate.

Layers 1 to nkml are the mixed layer, nkml+1 to nkml+nkbl are the buffer layers. The results of this subroutine are mathematically identical if there are multiple pieces of the mixed layer with the same density or if there is just a single layer. There is no stability limit on the time step.

The key parameters for the mixed layer are found in the control structure. These include mstar, nstar, nstar2, pen_SW_frac, pen_SW_scale, and TKE_decay. For the :cite:`oberhuber1993` and :cite:`kraus1967` mixed layers, the values of these are:

.. list-table:: Model variables used in the bulk mixed layer
   :header-rows: 1

   * - Symbol
     - Value in Oberhuber (1993)
     - Value in Kraus-Turner (1967)
   * - pen_SW_frac
     - 0.42
     - 0.0
   * - pen_SW_scale
     - 15.0 m
     - 0.0 m
   * - mstar
     - 1.25
     - 1.25
   * - nstar
     - 1
     - 0.4
   * - TKE_decay
     - 2.5
     - 0.0
   * - conv_decay
     - 0.5
     - 0.0


TKE_decay is :math:`1/\kappa` in eq. 28 of :cite:`oberhuber1993`, while conv_decay is :math:`1/\mu`. Conv_decay has been eliminated in favor of the well-calibrated form for the efficiency of penetrating convection from :cite:`wang2003`.
