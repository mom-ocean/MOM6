from __future__ import unicode_literals

import numpy as np
from numpy import ma
from matplotlib import scale as mscale
from matplotlib import transforms as mtransforms
from matplotlib.ticker import Formatter, FixedLocator, MaxNLocator, AutoLocator


class VerticalSplitScale(mscale.ScaleBase):
    """
    Scales data in range -pi/2 to pi/2 (-90 to 90 degrees) using
    the system used to scale latitudes in a Mercator projection.

    """

    # The scale class must have a member ``name`` that defines the
    # string used to select the scale.  For example,
    # ``gca().set_yscale("splitscale")`` would be used to select this
    # scale.
    name = 'splitscale'

    def __init__(self, axis, **kwargs):
        """
        Any keyword arguments passed to ``set_xscale`` and
        ``set_yscale`` will be passed along to the scale's
        constructor.

        thresh: The degree above which to crop the data.
        """
        mscale.ScaleBase.__init__(self)
        #thresh = kwargs.pop("thresh", (85 / 180.0) * np.pi)
        #if thresh >= np.pi / 2.0:
        #    raise ValueError("thresh must be less than pi/2")
        zval = kwargs.pop("zval", None)
        if zval == None:
            raise Exception("zval must be specified")
        if len(zval) < 3:
            raise Exception("zval must be at least 3 long")
        zfrac = kwargs.pop("zfrac", None)
        if zfrac == None:
            zfrac = np.linspace(0., 1., len(zval))
        if len(zfrac) != len(zval):
            raise Exception("zval and zfrac must have the same length")
        #zval[0] = zval[1] + 1e4*(zval[0] - zval[1])
        #zfrac[0] = zfrac[1] + 1e4*(zfrac[0] - zfrac[1])
        #zval[-1] = zval[-2] + 1e4*(zval[-1] - zval[-2])
        #zfrac[-1] = zfrac[-2] + 1e4*(zfrac[-1] - zfrac[-2])
        self.zval = np.array(zval)
        self.zfrac = np.array(zfrac)

    def get_transform(self):
        """
        Override this method to return a new instance that does the
        actual transformation of the data.

        The VerticalSplitScaleTransform class is defined below as a
        nested class of this one.
        """
        return self.VerticalSplitScaleTransform(self.zval, self.zfrac)

    def set_default_locators_and_formatters(self, axis):
        """
        Override to set up the locators and formatters to use with the
        scale.  This is only required if the scale requires custom
        locators and formatters.  Writing custom locators and
        formatters is rather outside the scope of this example, but
        there are many helpful examples in ``ticker.py``.

        In our case, the Mercator example uses a fixed locator from
        -90 to 90 degrees and a custom formatter class to put convert
        the radians to degrees and put a degree symbol after the
        value::
        """
        #class DegreeFormatter(Formatter):
        #    def __call__(self, x, pos=None):
        #        # \u00b0 : degree symbol
        #        return "%d\u00b0" % ((x / np.pi) * 180.0)

        #axis.set_major_locator(MaxNLocator(10))
        ticks = None
        for k in range(1,len(self.zval)):
          nbins = 16 * ( self.zfrac[k] - self.zfrac[k-1] )
          newticks = MaxNLocator(nbins=nbins, steps=[1,2,2.5,5,10]).tick_values(self.zval[k-1], self.zval[k])
          if ticks==None: ticks = newticks
          else: ticks = np.append( ticks, newticks )
        ticks = [x for x in ticks if (x>+self.zval.min() and x<=self.zval.max())] # Only used tickes within range
        ticks = np.sort( ticks ) # Fix due to different python versions on PP and workstations!
        axis.set_major_locator( FixedLocator( ticks ) )
        #axis.set_major_locator(AutoLocator())
        #deg2rad = np.pi / 180.0
        #axis.set_major_locator(FixedLocator(
                #np.arange(-90, 90, 10) * deg2rad))
        #axis.set_major_formatter(DegreeFormatter())
        #axis.set_minor_formatter(DegreeFormatter())

    def limit_range_for_scale(self, vmin, vmax, minpos):
        """
        Override to limit the bounds of the axis to the domain of the
        transform.  In the case of Mercator, the bounds should be
        limited to the threshold that was passed in.  Unlike the
        autoscaling provided by the tick locators, this range limiting
        will always be adhered to, whether the axis range is set
        manually, determined automatically or changed through panning
        and zooming.
        """
        return min(vmin, self.zval[0]), max(vmax, self.zval[-1])

    class VerticalSplitScaleTransform(mtransforms.Transform):
        # There are two value members that must be defined.
        # ``input_dims`` and ``output_dims`` specify number of input
        # dimensions and output dimensions to the transformation.
        # These are used by the transformation framework to do some
        # error checking and prevent incompatible transformations from
        # being connected together.  When defining transforms for a
        # scale, which are, by definition, separable and have only one
        # dimension, these members should always be set to 1.
        input_dims = 1
        output_dims = 1
        is_separable = True

        def __init__(self, zval, zfrac):
            mtransforms.Transform.__init__(self)
            self.zval = zval
            self.zfrac = zfrac

        def transform_non_affine(self, a):
            """
            This transform takes an Nx1 ``numpy`` array and returns a
            transformed copy.  Since the range of the Mercator scale
            is limited by the user-specified threshold, the input
            array must be masked to contain only valid values.
            ``matplotlib`` will handle masked arrays and remove the
            out-of-range data from the plot.  Importantly, the
            ``transform`` method *must* return an array that is the
            same shape as the input array, since these values need to
            remain synchronized with values in the other dimension.
            """
            return np.interp(-a, -self.zval, self.zfrac)

        def inverted(self):
            """
            Override this method so matplotlib knows how to get the
            inverse transform for this transform.
            """
            return VerticalSplitScale.InvertedVerticalSplitScaleTransform(self.zval, self.zfrac)

    class InvertedVerticalSplitScaleTransform(mtransforms.Transform):
        input_dims = 1
        output_dims = 1
        is_separable = True

        def __init__(self, zval, zfrac):
            mtransforms.Transform.__init__(self)
            self.zval = zval
            self.zfrac = zfrac

        def transform_non_affine(self, a):
            #return np.arctan(np.sinh(a))
            #return 0.5*(np.array(a)-1000.)
            return np.interp(a, self.zfrac, -self.zval)

        def inverted(self):
            return VerticalSplitScale.VerticalSplitScaleTransform(self.zval, self.zfrac)

# Now that the Scale class has been defined, it must be registered so
# that ``matplotlib`` can find it.
mscale.register_scale(VerticalSplitScale)


if __name__ == '__main__':
    import matplotlib.pyplot as plt

    z = np.linspace(-6500., 0., 43)
    s = -1. * z

    plt.plot(s, z, '.-', lw=2)
    plt.axhline(-1000.)
    plt.gca().set_yscale('splitscale', zval=[0.,-1000.,-9000.])

    plt.xlabel('Depth')
    plt.ylabel('Z')
    plt.grid(True)

    plt.show()
