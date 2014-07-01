"""
A collection of useful functions...
"""
import numpy as np

def section2quadmesh(x, z, q, representation='pcm'):
  """
  Creates the appropriate quadmesh coordinates to plot a scalar q(1:nk,1:ni) at
  horizontal positions x(1:ni+1) and between interfaces at z(nk+1,ni), using
  various representations of the topography.

  Returns X(2*ni+1), Z(nk+1,2*ni+1) and Q(nk,2*ni) to be passed to pcolormesh.

  TBD: Optionally, x can be dimensioned as x(ni) in which case it will be extraplated as if it had 
  had dimensions x(ni+1).
  
  Optional argument:
  
  representation='pcm' (default) yields a step-wise visualization, appropriate for
           z-coordinate models.
  representation='plm' yields a piecewise-linear visualization more representative
           of general-coordinate (and isopycnal) models.
  representation='linear' is the aesthetically most pleasing but does not
           represent the data conservatively.

  """

  if x.ndim!=1: raise Exception('The x argument must be a vector')
  if z.ndim!=2: raise Exception('The z argument should be a 2D array')
  if q.ndim!=2: raise Exception('The z argument should be a 2D array')
  qnk, qni = q.shape
  znk, zni = z.shape
  xni = x.size
  if zni!=qni: raise Exception('The last dimension of z and q must be equal in length')
  if znk!=qnk+1: raise Exception('The first dimension of z must be 1 longer than that of q. q has %i levels'%qnk)
  if xni!=qni+1: raise Exception('The length of x must 1 longer than the last dimension of q')

  if type( z ) == np.ma.core.MaskedArray: z[z.mask] = 0
  if type( q ) == np.ma.core.MaskedArray: qmin = np.amin(q); q[q.mask] = qmin

  periodicDomain =  abs((x[-1]-x[0])-360. ) < 1e-6 # Detect if horizontal axis is a periodic domain

  if representation=='pcm':
    X = np.zeros((2*qni))
    X[::2] = x[:-1]
    X[1::2] = x[1:]
    Z = np.zeros((qnk+1,2*qni))
    Z[:,::2] = z
    Z[:,1::2] = z
    Q = np.zeros((qnk,2*qni-1))
    Q[:,::2] = q
    Q[:,1::2] = ( q[:,:-1] + q[:,1:] )/2.
  elif representation=='linear':
    X = np.zeros((2*qni+1))
    X[::2] = x
    X[1::2] = ( x[0:-1] + x[1:] )/2.
    Z = np.zeros((qnk+1,2*qni+1))
    Z[:,1::2] = z
    Z[:,2:-1:2] = ( z[:,0:-1] + z[:,1:] )/2.
    Z[:,0] = z[:,0]
    Z[:,-1] = z[:,-1]
    Q = np.zeros((qnk,2*qni))
    Q[:,::2] = q
    Q[:,1::2] = q
  elif representation=='plm':
    X = np.zeros((2*qni))
    X[::2] = x[:-1]
    X[1::2] = x[1:]
    # PLM reconstruction for Z
    dz = np.roll(z,-1,axis=1) - z # Right-sided difference
    if not periodicDomain: dz[:,-1] = 0 # Non-periodic boundary
    d2 = ( np.roll(z,-1,axis=1) - np.roll(z,1,axis=1) )/2. # Centered difference
    d2 = ( dz + np.roll(dz,1,axis=1) )/2. # Centered difference
    s = np.sign( d2 ) # Sign of centered slope
    s[dz * np.roll(dz,1,axis=1) <= 0] = 0 # Flatten extrema
    dz = np.abs(dz) # Only need magnitude from here on
    S = s * np.minimum( np.abs(d2), np.minimum( dz, np.roll(dz,1,axis=1) ) ) # PLM slope
    Z = np.zeros((qnk+1,2*qni))
    Z[:,::2] = z - S/2.
    Z[:,1::2] = z + S/2.
    Q = np.zeros((qnk,2*qni-1))
    Q[:,::2] = q
    Q[:,1::2] = ( q[:,:-1] + q[:,1:] )/2.
  else: raise Exception('Unknown representation!')

  return X, Z, Q


def rho_Wright97(S, T, P=0):
  """
  Returns the density of seawater for the given salinity, potential temperature
  and pressure.

  Units: salinity in PSU, potential temperature in degrees Celsius and pressure in Pascals.
  """
  a0 = 7.057924e-4; a1 = 3.480336e-7; a2 = -1.112733e-7
  b0 = 5.790749e8;  b1 = 3.516535e6;  b2 = -4.002714e4
  b3 = 2.084372e2;  b4 = 5.944068e5;  b5 = -9.643486e3
  c0 = 1.704853e5;  c1 = 7.904722e2;  c2 = -7.984422
  c3 = 5.140652e-2; c4 = -2.302158e2; c5 = -3.079464
  al0 = a0 + a1*T + a2*S
  p0 = b0 + b4*S + T * (b1 + T*(b2 + b3*T) + b5*S)
  Lambda = c0 + c4*S + T * (c1 + T*(c2 + c3*T) + c5*S)
  return (P + p0) / (Lambda + al0*(P + p0))


def ice9(i, j, source, xcyclic=True, tripolar=True):
  """
  An iterative (stack based) implementation of "Ice 9".

  The flood fill starts at [j,i] and treats any positive value of "source" as
  passable. Zero and negative values block flooding.

  xcyclic = True allows cyclic behavior in the last index. (default)
  tripolar = True allows a fold across the top-most edge. (default)

  Returns an array of 0's and 1's.
  """
  wetMask = 0*source
  (nj,ni) = wetMask.shape
  stack = set()
  stack.add( (j,i) )
  while stack:
    (j,i) = stack.pop()
    if wetMask[j,i] or source[j,i] <= 0: continue
    wetMask[j,i] = 1
    if i>0: stack.add( (j,i-1) )
    elif xcyclic: stack.add( (j,ni-1) )
    if i<ni-1: stack.add( (j,i+1) )
    elif xcyclic: stack.add( (0,j) )
    if j>0: stack.add( (j-1,i) )
    if j<nj-1: stack.add( (j+1,i) )
    elif tripolar: stack.add( (j,ni-1-i) ) # Tri-polar fold
  return wetMask


def maskFromDepth(depth, zCellTop):
  """
  Generates a "wet mask" for a z-coordinate model based on relative location of
  the ocean bottom to the upper interface of the cell.

  depth (2d) is positiveo
  zCellTop (scalar) is a negative position of the upper interface of the cell..
  """
  wet = 0*depth
  wet[depth>-zCellTop] = 1
  return wet

# Tests
if __name__ == '__main__':

  import matplotlib.pyplot as plt
  import numpy.matlib

  # Test data
  x=np.arange(5)
  z=np.array([[0,0.2,0.3,-.1],[1,1.5,.7,.4],[2,2,1.5,2],[3,2.3,1.5,2.1]])*-1
  q=np.matlib.rand(3,4)
  print 'x=',x
  print 'z=',z
  print 'q=',q

  X, Z, Q = section2quadmesh(x, z, q)
  print 'X=',X
  print 'Z=',Z
  print 'Q=',Q
  plt.subplot(3,1,1)
  plt.pcolormesh(X, Z, Q)

  X, Z, Q = section2quadmesh(x, z, q, representation='linear')
  print 'X=',X
  print 'Z=',Z
  print 'Q=',Q
  plt.subplot(3,1,2)
  plt.pcolormesh(X, Z, Q)

  X, Z, Q = section2quadmesh(x, z, q, representation='plm')
  print 'X=',X
  print 'Z=',Z
  print 'Q=',Q
  plt.subplot(3,1,3)
  plt.pcolormesh(X, Z, Q)

  plt.show()
