#!/usr/bin/env python
# vim: set ts=4 sw=4 et:

from scipy import *
from pylab import *
from ophidian import blob

def edge_points( a ) :
    """
    Find the points just inside the edge of the plasma.
    """
    l = []
    for i in range(len(a)) :
        for j in range(len(a[0])) :
            if a[i][j] > 0.0 and            \
                    ( a[i+1][j] < 0.0 or    \
                    a[i-1][j] < 0.0 or      \
                    a[i][j+1] < 0.0 or      \
                    a[i][j-1] < 0.0 ) :
                l.append([i,j])
    return l

def lcfs_points( a, i, j, r_coor, z_coor ) :
    """
    Find the interpolated points for the LCFS.
    """
    p = []
    s = blob.spline()
   
    def ss( RR ) :
        """
        Wrapper function for the spline evaluation. This is used
        to catch out-of-bounds guesses from the numerical root
        finder.

        Note that this function uses the spline objects, so beware
        of side effects.
        """
        if RR >= s.R[0] and RR <= s.R[-1] :
            try :
                ZZ = s.get_Z_for_R( RR )
            except blob.SplineError, e :
                # print e
                return 1.0
            if type(ZZ) == float :
                return ZZ
            else :
                return 1.0
        else :
            return 1.0
    
    # in the R direction
    if a[i, j+1] < 0.0 or a[i, j-1] < 0.0 :
        R = []
        for k in range( len(a[0, :]) ) :
            R.append( r_coor(k) )
        s.build( R, a[i, :] )
        p.append( [ optimize.fsolve( ss, r_coor(j) ), z_coor(i) ] )
            
    # in the Z direction
    if a[i+1, j] < 0.0 or a[i-1, j] < 0.0 :
        Z = []
        for k in range( len(a[:, 0]) ) :
            Z.append( z_coor(k) )
        s.build( Z, a[:, j] )
        p.append( [ r_coor(j), optimize.fsolve( ss, z_coor(i) ) ] )
            
    # FIXME : handle points with two adjacent roots 
    # along one axis!
    
    return p
