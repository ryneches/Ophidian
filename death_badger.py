#!/usr/bin/env python
# vim: set ts=4 sw=4 et:

from Numeric import *
from scipy import *
from pylab import *

R  = []
Z  = []
Z2 = []

last_i = 0

class SplineError(Exception) :
    pass

def build (R_new, Z_new, natural = 'yes'):
    """
    Takes list of points in a Cartesian plane and computes the spline
    coefficeints. See Numerical Recipies in C, 2nd edition, section
    3.3 for details.

    Defaults to natural boundary conditions.
    """
    
    if len(R_new) != len(Z_new) :
        raise SplineError, "Dimension Mismatch: R: " + \
        str(R_new) + " Z: " + str(Z_new)
        
    R = R_new
    Z = Z_new
  

    data = []
    for i in range(len(R)) :
        data.append([ float(R[i]), float(Z[i]) ])

    data.sort()

    for i in range(len(R)) :
        R[i] = data[i][0]
        Z[i] = data[i][1]

    for i in range(len(R) - 1) :
        if R[i] == R[i + 1] :
            raise SplineError, "R values must be unique. Value: " + \
                str(R[i])
    
    # assorted variables
    
    i = 1
    k = 1
    n = len(R) # number of elements
    p = 0.0
    qn = 0.0
    sig = 0.0
    un = 0.0
    u = []
       
    # first derivatives at the boundaries
    Zp1 = 1e5
    Zpn = -1e5
       
    # initialize Z2 and u arrays to the correct size 
    # (values don't matter)

    Z2 = range(n)
    u = range(n)
      
    if natural == 'yes' :
        Z2[0] = u[0] = 0.0  # first real value here is also zero, 
                            # because this is a natural spline
        
        qn = un = Z2[n - 1] = u[n - 1] = 0.0 # "natural" upper boundary
     
    else :
        # set the lower boundary condition to match Zp1
        Z2[0] = -0.5
        u[0] = ( 3.0 / ( R[1] - R[0] ) ) * \
            ( (Z[1] - Z[0] ) / ( R[1] - R[0] ) - Zp1 )
      
        qn = 0.5
        un = ( 3.0 / ( R[n - 1] - R[n - 2] ) ) * \
            ( Zpn - ( Z[n - 1] - Z[n - 2] ) / ( R[n - 1] - R[n - 2] ) )
    
    # tridiagonal decomposition (Z2 used as temporary storage)
    # We loop from the second element to the second to last 
    # element
    for i in range(1, n - 1) :
        sig = ( R[i] - R[i - 1] ) / ( R[i + 1] - R[i - 1] )
        p = sig * Z2[i - 1] + 2.0
        Z2[i] = ( sig - 1.0 ) / p
        u[i] = (6.0 * ((Z[i + 1] - Z[i]) / (R[i + 1] - R[i]) - \
            (Z[i] - Z[i - 1]) / (R[i] - R[i - 1])) / ( R[i + 1] - \
            R[i - 1]) - sig*u[i - 1]) / p
       
    for k in range( n - 2, 0, -1) :
        Z2[k] = Z2[k] * Z2[k + 1] + u[k]
   

def getZ ( r ):
    """
    find the value Z for a given r
    """
    
    # make sure r is treated as a float
    r = float(r)
       
    if r < R[0] or r > R[len(R) - 1] :
        raise SplineError, "Called value out of range. Value: " + str(r)
    
    # find the right place in the table of R values
    for i in range(len(R) - 1) :
        if r >= R[i] and r < R[i + 1] :
            break
   
    print "index :", i

    h = R[i + 1] - R[i]

    if h == 0.0 :
        raise SplineError, "Bad spline! Value: ", str(r)
    
    a = (R[i + 1] - r) / h
    b = (r - R[i]) / h
    
    return a * Z[i] + b * Z[i + 1] + \
        ((a**3 - a) * Z2[i] + (b**3 - b) * Z2[i + 1]) * (h**2) / 6.0


def getiZ ( r ):
    """
    find the value Z for a given r
    """
    
    # make sure r is treated as a float
    r = float(r)
       
    if r < R[0] or r > R[len(R) - 1] :
        raise SplineError, "Called value out of range. Value: " + str(r)
    
    # find the right place in the table of R values
    for i in range(len(R) - 1) :
        if r >= R[i] and r < R[i + 1] :
            break
    
    print "index :", i

    h = R[i + 1] - R[i]

    if h == 0.0 :
        raise SplineError, "Bad spline! Value: ", str(r)
    
    return (1/(6*h))*(-3*r**2*Z[i] + 6*r*R[i+1]*Z[i] + 3*r**2*Z[i+1] -              \
        6*r*R[i]*Z[i+1] - (r**4*Z[i+1][i])/4 + r**3*R[i+1]*Z[i+1][i] +              \
        (1/2)*r**2*(h**2 - 3*R[i+1]**2)*Z[i+1][i] - r*R[i+1]*(h**2 -                \
        R[i+1]**2)*Z2[i] + (r**4*Z2[i+1])/4 - r**3*R[i]*Z2[i+1] + r*R[i]*(h**2      \
        - R[i]**2)*Z2[i+1] + (1/2)*r**2*(-h**2 + 3*R[i]**2)*Z2[i+1])



