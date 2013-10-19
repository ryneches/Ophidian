#!/usr/bin/env python
# vim: set ts=4 sw=4 et:
"""
Fit ACUBE imput functions to CUBE data.
"""

from ACUBE import *
from scipy import *
from pylab import *
from optparse import OptionParser

# Option parsing!

parser = OptionParser()
parser.add_option("-f", "--file", dest="filename",
                help="read CUBE output from FILE", metavar="FILE")
parser.add_option("-m", "--Rmin", dest="Rmin",
                help="value for Rmin", metavar="RMIN")
parser.add_option("-x", "--Rmax", dest="Rmax",
                help="value for Rmax", metavar="RMAX")
parser.add_option("-c", "--cutoff", dest="cuttoff",
                help="truncate data after this radial position", 
                metavar="CUTOFF")

(options, args) = parser.parse_args()

file = options.filename
Rmin = float(options.Rmin)
Rmax = float(options.Rmax)
cuttoff = float(options.cuttoff)
# Slurp up cube data!

T = data.CUBEdata()
T.read( file )

Rmin, Rmax = T.plasma_range()
axis = T.plasma_axis()

#T.slice( Rmin, cuttoff ) 

print "[Rmin:Rmax] =", [Rmin,Rmax]

T.slice( Rmin, axis )
print "sliced at Rmin =", Rmin, "axis =", axis

# Fitting!

#V = [ 1500000.0, -1000.0, 100.0, 100.0, 100.0, 100.0, 100, 100, 100, 100, 100, 100]

#def pressure( V, R ) :
#    p0, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12 = V
#    w = Rmax - Rmin
#    r = Rmax - R
#    a = p2 * r**2 + p3 * r**3 + p4 * r**4 + p5 * r**5 + p6 * r**6 + p7 * r**7 + p8 * r**8 + p9 * r**9 + p10 * r**10 + p11 * r**11 + p12 * r**12
#    b = p2 * w**2 + p3 * w**3 + p4 * w**4 + p5 * w**5 + p6 * r**6 + p7 * w**7 + p8 * w**8 + p9 * w**9 + p10 * r**10 + p11 * w**11 + p12 * w**12
#
#    return p0 * ( 1 - ( a / b ) )

#V = [ 0, 1000, 2, 120, 5, 3, 120]

#def pressure( V, R ) :
#    p0, p1, p2, p3, p4, p5, p6 = V
#    return p0 + (p1 - e**(-p2 * R - p3) ) / ( p4 + e**(-p5 * R - p6) )

V = [ 1, 1, 1, 1 ]

def pressure( V, R ) :
    p0, p1, p2, p3 = V
    return p0 + p1 * e**( p2 * R - p3 )

def p_residuals( V, y, x ) :
    return y - pressure( V, x )

W = optimize.leastsq(p_residuals, V, args=(T.p, T.R))[0]

O = [ 1, 1, 1, 1 ]
def psi ( O, R ) :
    a, b, c, d = O
    return a + b * e**( c * R - d )

def psi_residuals( O, y, x ) :
    #a, b, c, d = O
    return y - psi( O, x )

P = optimize.leastsq(psi_residuals, O, args=(T.psi, T.R))[0]

#N = [ 1, 1, 1 ]
#def dpsi( N, R ) :
#    a_psi, b_psi, c_psi = N
#    return sqrt(R - Rmin) * sqrt(Rmax - R) * ( a_psi + b_psi *  \
#        (R - Rmin) + c_psi * (R - Rmin)**2 )

N = [ 1, 1, 1, 1 ]
def dpsi( N, R ) :
    a_psi, b_psi, c_psi, d_psi = N
    return b_psi * c_psi * e**( c_psi * R - d_psi )

def dpsi_residuals( N, y, x ) :
    return y - dpsi( N, x )

s = blob.spline()
s.build( T.R, T.psi, natural = 'yes' )
T.dpsi = []
for i in T.R :
    T.dpsi.append( s.get_dZ_for_R(i) )

M = optimize.leastsq(dpsi_residuals, N, args=(T.dpsi, T.R) )[0]

print "pressure coefficients:", W
print "psi coefficients:", P
print "dpsi coefficients:", M
subplot(311)
plot(T.R, T.p)
plot(T.R, pressure(W, T.R))
#plot(T.R, pressure(V, T.R))
subplot(312)
plot(T.R, T.psi)
plot(T.R, psi(P, T.R))
subplot(313)
#plot(T.R, T.psi) 
plot(T.R, T.dpsi)
plot(T.R, dpsi(M, T.R))
show()
