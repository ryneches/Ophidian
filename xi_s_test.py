#!/usr/bin/env python
# vim: set ts=4 sw=4 et:

from Numeric import *
from scipy import *
from pylab import *
from ACUBE2 import *
from ACUBE import *

a = io.array_import.read_array( '../CUBE_data/cube_results_2/rerun/bin/flux.txt' )
t = exp_pressure.tokamak()

t.p0 = -2.29626365
t.p1 = 17.80667626
t.p2 = 2.66593021
t.p3 = 10.16174372
t.a_psi = -0.03166988
t.b_psi = 0.0345172
t.c_psi = 0.24642886
t.d_psi = 1.07627846
t.BtRmin = 0.268526591
t.axis = 5.76492437379
t.axis_label = 5.76492437379

s = t.find_lcfs( a, 3.9, 6.1, -1.1, 1.1 )

t.Rmin = s[0].R[0]
t.Rmax = s[0].R[-1]

T = data.CUBEdata()
T.read( '../CUBE_data/cube_results_2/rerun/bin/results.txt' )

p_s = blob.spline()
Psi_s = blob.spline()

p_s.build( T.R_raw, T.p_raw )
Psi_s.build( T.R_raw, T.psi_raw )

t.p_s   = p_s
t.Psi_s = Psi_s

RR = arange( t.Rmin, t.Rmax, 0.001 ) 

def Z ( xi, R ) :
    R0 = ( t.Rmin + t.Rmax ) / 2.0
    a0 = t.Rmax - R0

    b = ( R - R0 ) / ( a0 - xi )

    if b**2 > 1 :
        return 0

    return ( a0 - xi ) * sqrt( 1 - b**2 )

rr = arange(t.Rmin, 5.7, 0.1)

for rh in rr :
    
    xixi = []
    XiXi = []

    for i in RR :
        xixi.append(Z(t.xi( rh, i ),i) )
        XiXi.append(Z(t.xi_s( rh, i ),i))

    plot(RR,XiXi, color='blue')
    plot(RR,xixi, color='green')

show()
