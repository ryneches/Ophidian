#!/usr/bin/env python
# vim: set ts=4 sw=4 et:
"""
Plot R verses Rhat from a pickled fluxbucket.
"""

from Numeric import *
from scipy import *
from pylab import *
from optparse import OptionParser
from ACUBE import *

parser = OptionParser()
parser.add_option("-f", "--file", dest="filename",
                  help="write results to FILE", metavar="FILE")

(options, args) = parser.parse_args()

picklefile = options.filename

def tangent( x1, z1, x2, z2, x ) :
    return ( ( z2 - z1 ) / ( x2 - x1 ) ) * ( x - x1 ) + z1

# read the fluxbucket from disk...
print "Unpickling fluxbucket..."
f = open(picklefile, 'r')
tokamak = pickle.load(f)
f.close()
print "Done."

# spill off surfaces that don't have any data
cleanbucket = []
for surface in tokamak.surfaces :
    if sum(surface.spline.Z) != 0 :
        cleanbucket.append(surface)   

# plot the raw R-vs-Rhat curve
Rhs = []
roots = []
for fs in cleanbucket :
    Rhs.append(fs.Rhat)
    roots.append(fs.root)

roots.reverse()
Rhs.reverse()
RhRv = []
for i in Rhs :
    RhRv.append(i)
Rhs.reverse()

RvRhatSpline = blob.spline()
RvRhatSpline.build(Rhs + roots, Rhs + RhRv, natural = 'yes')

ynew = []
xnew = []
x = arange(min(RvRhatSpline.R), max(RvRhatSpline.R), 0.001)
for i in x :
    z = RvRhatSpline.get_Z_for_R(i)
    if type(z) == float:
        ynew.append(RvRhatSpline.get_Z_for_R(i))
        xnew.append(i)

# plot the current
def J(R) :
    return tokamak.dp(RvRhatSpline.get_Z_for_R(R))
    #return tokamak.dPsi(RvRhatSpline.get_Z_for_R(R))
    #return RvRhatSpline.get_Z_for_R(R)**3 / R - R
    #return tokamak.p(RvRhatSpline.get_Z_for_R(R))

    #return physics.constants.mu_0 * (-1/R)          \
    #    * tokamak.dp(RvRhatSpline.get_Z_for_R(R))   \
    #    / tokamak.dPsi(RvRhatSpline.get_Z_for_R(R)) \
    #    * ( R**2 - RvRhatSpline.get_Z_for_R(R)**2 )
    
    #return physics.constants.mu_0 * (1/R) * tokamak.dp(R)   \
    #    / tokamak.dPsi(RvRhatSpline.get_Z_for_R(R))         \
    #    * ( R**2 - RvRhatSpline.get_Z_for_R(R)**2 )
q = []
for i in x :
    q.append(J(i))

xlabel("R")
ylabel('dp')
plot(x, q, color="green")

show()
