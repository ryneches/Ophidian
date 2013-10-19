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

x = arange(5.8, 7.0, 0.1)

# spill off surfaces that don't have any data or have
# bogus roots
cleanbucket = []
for surface in tokamak.surfaces :
    if sum(surface.spline.Z) != 0       \
        and surface.root > tokamak.Rmin \
        and surface.root < tokamak.Rmax \
        and surface.root > 5.66 :
        cleanbucket.append(surface)   
     
print "surfaces:", len(cleanbucket)

theRs = []
theRhats = []

subplot(421)

x = arange(tokamak.Rmin, tokamak.Rmax, 0.01)
y = tokamak.l(x)
xlabel("R")
ylabel("Z")
plot(x,y)
for fs in cleanbucket :
    plot(fs.spline.R, fs.spline.Z, color="blue")

axvline(x=tokamak.axis_label, linewidth=2, color='r')

subplot(422)

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
plot(Rhs + roots, Rhs + RhRv, color="blue" )

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

plot(roots, RhRv, 'ro')

xlabel("R")
ylabel("Rhat")
plot(xnew, ynew, color="red")

# plot the current
subplot(423)

#x = arange(tokamak.Rmin, tokamak.Rmax, 0.01)
y = []
for i in x :
    y.append(tokamak.p(RvRhatSpline.get_Z_for_R(i)))

pressure = blob.spline()
pressure.build( x, y, natural='yes')

def J(R) :
    #return tokamak.dp(RvRhatSpline.get_Z_for_R(R))
    #return tokamak.dPsi(RvRhatSpline.get_Z_for_R(R))
    #return RvRhatSpline.get_Z_for_R(R)**3 / R - R
    #return tokamak.p(RvRhatSpline.get_Z_for_R(R))
    return pressure.get_dZ_for_R(R)
    
    #if tokamak.dPsi(RvRhatSpline.get_Z_for_R(R)) == 0.0 :
    #    return 0.0
    
    #return physics.constants.mu_0 * (-1/R)          \
    #    * pressure.get_dZ_for_R(R)                  \
    #    / tokamak.dPsi(RvRhatSpline.get_Z_for_R(R)) \
    #    * ( R**2 - RvRhatSpline.get_Z_for_R(R)**2 )
    
    #return physics.constants.mu_0 * (1/R) * tokamak.dp(R)   \
    #    / tokamak.dPsi(RvRhatSpline.get_Z_for_R(R))         \
    #    * ( R**2 - RvRhatSpline.get_Z_for_R(R)**2 )

q = []
for i in x :
    #q.append(J(i))
    q.append(tokamak.F(RvRhatSpline.get_Z_for_R(i)))

xlabel("R")
ylabel('J')
plot(x, q, color="green")

# dPsi
subplot(424)
x = arange(tokamak.Rmin, tokamak.Rmax, 0.01)
xlabel("R")
ylabel('\over{d\psi}{d\hat{R}')
plot(x, tokamak.dPsi(x))

# psi
subplot(425)
xlabel("Rhat")
ylabel("Psi")
plot(x, tokamak.psi(x))

# pressure
subplot(426)
xlabel("Rhat")
ylabel("Pressure")
plot(x, tokamak.p(x))

# q(R)
subplot(427)
xlabel("R")
ylabel("q(R)")
x = arange(tokamak.Rmin, tokamak.Rmax, 0.01)
y = []
if not hasattr(tokamak, 'BtRmin') :
    tokamak.BtRmin = 1.0

for i in x:
    y.append(tokamak.q(RvRhatSpline.get_Z_for_R(i)))
plot(x, y)

# dp(R)
subplot(428)
xlabel("Rhat")
ylabel("dp")
x = arange(tokamak.Rmin, tokamak.Rmax, 0.01)
plot(x, tokamak.dp(x))

show()

