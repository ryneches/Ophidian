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

# override default functions with corrected versions
class solution_circular_corrected(data.solution_circular) :

    def __init__ (self) :
        data.solution_circular.__init__(self)
        self.Raxis = 0.0

    def p( self, R) :
        """
        Corrected pressure fuction.
        """
        # shortcuts!
        Rmin  = self.Rmin
        Rmax  = self.Rmax
        Raxis = self.Raxis
        p0    = self.p0
        p2    = self.p2
        p3    = self.p3
        p4    = self.p4

        if R < Rmin or R > Rmax :
            return 0

        w = Raxis - Rmin
        r = Rraxis - R
        a = p2 * r**2 + p3 * r**3 + p4 * r**4
        b = p2 * w**2 + p3 * w**3 + p4 * w**4

        return p0 * ( 1 - ( a / b ) )

    def xi( self, R, r ) :
        """
        Corrected xi function
        """

# build and populate our machine...
tokamak = data.solution_circular()

tokamak.Rmin  = 4.0
tokamak.Rmax  = 6.0
tokamak.Raxis = 5.8
tokamak.al    = 1.0
tokamak.bl    = 1.0
tokamak.cl    = 1.0
tokamak.a_psi = 0.01
tokamak.b_psi = 1.0
tokamak.c_psi = 0.001
tokamak.p0    = 100000000.0
tokamak.p2    = 8.0
tokamak.p3    = -16.0 / 3.0
tokamak.p4    = 1.0
tokamak.root_err = 1e-4
tokamak.axis_err = 1e-4
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
print "Done."
