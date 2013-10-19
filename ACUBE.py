#!/usr/bin/env python
# vim: set ts=4 sw=4 et:
"""
Simple analytic plotting utility for high beta tokamaks.
"""

from __future__ import division
from Numeric import *
from scipy import *
from pylab import *
from optparse import OptionParser
from ACUBE import *

# Option parsing!
parser = OptionParser()
parser.add_option("-f", "--file", dest="filename",
                help="write results to FILE", metavar="FILE")
parser.add_option("-n", "--surfaces", dest="n",
                help="Number of flux surfaces to compute", metavar="N")
parser.add_option("-r", "--resolution", dest="steps",
                help="Number of points per flux surface", metavar="STEPS")

(options, args) = parser.parse_args()

picklefile  = options.filename
n           = int(options.n)
steps       = int(options.steps)

# build and populate our machine...
tokamak = data.solution_circular()

tokamak.Rmin   = 4.0
tokamak.Rmax   = 6.0
tokamak.BtRmin = 3.0
tokamak.al     = 1.0
tokamak.bl     = 1.0
tokamak.cl     = 1.0
tokamak.a_psi  = 0.01
tokamak.b_psi  = 0.01
tokamak.c_psi  = 0.1
tokamak.p0     = 1500000.0
tokamak.p2     = 1.0
tokamak.p3     = 1.0
tokamak.p4     = 1.0
tokamak.root_err = 1e-1
tokamak.axis_err = 1e-1
tokamak.resolution = steps

# plot the last closed flux surface
x = arange(tokamak.Rmin, tokamak.Rmax, 0.01)
y = tokamak.l(x)
plot(x,y, color="black")

# find the magnetic axis

tokamak.find_axis()
axvline(x=tokamak.axis_label, linewidth=2, color='r')

rhats = arange(tokamak.Rmin, tokamak.axis_label, \
    (tokamak.axis_label - tokamak.Rmin) / (n + 1) )

#if len(rhats) == n + 1 :
#    # pop the last element off -- It'll always be Rmax, and 
#    # we don't need to plot that.
#    rhats = rhats.tolist()
#    rhats.reverse()
#    rhats.pop()
#    rhats.reverse()

#if len(rhats) == n + 2 :
#    # pop the last element off -- It'll always be Rmax, and 
#    # we don't need to plot that.
#    rhats = rhats.tolist()
#    rhats.pop()
#    rhats.reverse()
#    rhats.pop()
#    rhats.reverse()

# build the surfaces
for j in rhats :
    print "Plotting", j
    tokamak.add_surface(j, steps)

for fs in tokamak.surfaces:
    plot(fs.spline.R, fs.spline.Z, color="blue")

show()

# pickle the results
print "Pickling solution..."
f = open(picklefile, 'w')
pickle.dump(tokamak, f)
f.close()
print "Done."
