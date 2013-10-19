#!/usr/bin/env python
# vim: set ts=4 sw=4 et:

import ophidian.simplesnake
from scipy import *
from pylab import *

def BuildGrid(p) :
    """
    Build and return flux grid for problem p.
    
    Remember that the python array A[i,j] looks like this :
    	     
	         j
      +--------------
      |     0 1 2 3 4
      |   +----------
      | 0 | a b c d e
      | 1 | f g h i j
    i | 2 | k l m n o
      | 3 | p q r s t
      | 4 | u v w x y

    Mapping to configuration space, i maps to Z, j maps to R :
        
        z = p.z_coor( i )
        r = p.r_coor( j )

    """
    
    A = zeros( [p.NZ,p.NR], dtype=float )
    
    for i in range( p.NZ ) :
        for j in range( p.NR ) :
            if      p.r_coor( j ) > p.Rmin                              \
                and p.r_coor( j ) < p.Rmax                              \
                and p.z_coor( i ) > p.l2_s.get_Z_for_R( p.r_coor( j ) ) \
                and p.z_coor( i ) < p.l1_s.get_Z_for_R( p.r_coor( j ) ) :
                
                Xi = p.xi_geometric( p.r_coor( j ), p.z_coor( i ) )
                
                def a(r) :
                    print r
                    if r >= p.r_coor( j ) :
                        return r 
                    try :
                        xi = p.xi_s( r, Xi, p.r_coor( j ), force='no' )
                    except ophidian.simplesnake.GLimitError, e :
                        return r
                    if xi == None :
                        return r
                    else :
                        return (Xi - xi)**2
                
                Rh = optimize.fminbound( a, p.Rmin, p.Rmax, full_output=1 )[0]
                
                print [i,j], Rh
                
                A[i,j] = p.Psi_s.get_Z_for_R( Rh )

    return A
