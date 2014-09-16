import numpy
import math

def dzIter(nk, Htot, dzTop, Huniform, fnPow, prec):
    """
    Optimizes the highest ratio dzFn so that the sum(dz)=Htot.
    """
    def dzFn(nk, Huniform, dzTop, fnPow, zFac):
        """
        Returns dz = dzTop * ( 1 + zFac \int \fn(k) dk ) and sum(dz).
        """
        dz = numpy.ones(nk) * dzTop
        k = math.ceil( Huniform/dzTop )
        dz = dzTop * numpy.cumprod( 1. + zFac * (fn( numpy.linspace(1, nk, nk).astype(numpy.float64), k, nk )**fnPow) )
        return dz, numpy.sum(dz)

    def fn(z, z0, z1):
        """
        Cosine bell function between z0 and z1 s.t. f(z<z0)=0 and f(z>z1)=0.
        """
        zStar = (z-z0)/(z1-z0) # non-dimensional coordinate 0..1
        zStar = numpy.maximum(0., zStar)
        zStar = numpy.minimum(1., zStar)
        return 0.5*(1. - numpy.cos(2.*math.pi*zStar))

    def optimizeZfac(nk, Htot, dzTop, Huniform, fnPow,  prec):
        """
        Optimizes the highest ratio dzFn() so that the sum(dz)=Htot.
        """
        it = 0
        zc0 = 0; dz, H0 = dzFn( nk, Huniform, dzTop, fnPow,  zc0 )
        zc2 = 2; dz, H2 = dzFn( nk, Huniform, dzTop, fnPow,  zc2 )
        while H2-H0 > prec/8  and it<200: # Binary search
            zc1 = (zc0 + zc2)/2; dz, H1 = dzFn( nk, Huniform, dzTop, fnPow,  zc1 )
            if Htot<H1: zc2, H2 = 1.*zc1, 1.*H1
            else: zc0, H0 = 1.*zc1, 1.*H1
            it += 1
        #print 'zFac=',zc1,'max ratio=',(dz[1:]/dz[:-1]).max()
        return dz

    def roundDz(nk, Htot, dzTop, Huniform, fnPow,  prec):
        """
        Returns dz = optimizeZfac() rounded to the precision "prec", and sum(dz).
        """
        dz = prec * numpy.rint( optimizeZfac(nk, Htot, dzTop, Huniform, fnPow,  prec)/prec ) # Round to a given precision
        return dz, numpy.sum(dz)

    dz, H = roundDz( nk, Htot, dzTop, Huniform, fnPow,  prec )
    return dz

def zFromDz(dz):
    """
    Sums dz to return zInterface and zCenter.
    """
    zInt = zeros(dz.shape[0]+1)
    zInt[1:] = -numpy.cumsum(dz)
    zCenter = 0.5*( zInt[:-1] + zInt[1:] )
    return zInt, zCenter

nk = 75         # number of levels
Htot=6500       # deepest ocean point
Huniform = 0    # upper region of uniform resolution
dzTop = 2       # thickness of top level
fnPow = 1.42865 # ???
prec = .01      # precision to round thicknesses/depths todz = dzIter(nk, Htot, dzTop, Huniform, fnPow, prec)
dz = dzIter(nk, Htot, dzTop, Huniform, fnPow,  prec)
print 'dz=',dz
print 'sum(dz)=',numpy.sum(dz)
print 'max ratio=',(dz[1:]/dz[:-1]).max()
