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

# experimental version of tokamak object

class tokamak :
    """
    An experimental version of the basic class.
    """

    # data memebers

    def __init__ (self) :
        self.Rmin   = 0.0
        self.Rmax   = 0.0
        self.BtRmin = 0.0

        self.al, self.bl, self.cl           = 0.0, 0.0, 0.0
        self.p0, self.p1, self.p2, self.p3  = 0.0, 0.0, 0.0, 0.0
        self.a_psi, self.b_psi, self.c_psi  = 0.0, 0.0, 0.0
        self.axis                           = -1
        self.axis_label                     = -1
        self.root_err                       = 1e-4
        self.axis_err                       = 1e-4

    def l (self, R) :
        """
        The last closed flux surface, as a function of the major
        radius (R). This version is strictly circular.
        """    
        # shortcuts!
        Rmin = self.Rmin
        Rmax = self.Rmax

        if not Rmin <= R <= Rmax :
            return 0

        return sqrt(R - Rmin) * sqrt(Rmax - R)

    def dldR (self, R) :
        """
        Derivative of the last closed flux surface as a function
        of major radius.
        """
        # shortcuts!
        Rmin = self.Rmin
        Rmax = self.Rmax

        return -(sqrt(R - Rmin)/(2*sqrt(-R + Rmax))) +  \
            sqrt(-R + Rmax)/(2*sqrt(R - Rmin))

    def p (self, R) :
        """
        The pressure as a function of major radius (R).
        """
        # shortcuts!
        p0 = self.p0
        p1 = self.p1
        p2 = self.p2
        p3 = self.p3
        
        return p0 + p1 * e**( p2 * R - p3 )

    def dp (self, R) :
        """
        The gradient of pressure as a function of major radius (R).
        """
        # shortcuts!
        p0 = self.p0
        p1 = self.p1
        p2 = self.p2
        p3 = self.p3

        return p1 * p2 e**( p2 * R - p3 )

    def dpsi (self, r) :
        """
        The derivative of the poloidal flux as a function of
        flux label (r).
        """
        # shortcuts!
        Rmin  = self.Rmin
        Rmax  = self.Rmax
        a_psi = self.a_psi
        b_psi = self.b_psi
        c_psi = self.c_psi

        if r < Rmin or r > Rmax :
            return 0

        return self.l(r) * ( a_psi + b_psi * (r - Rmin) + \
            c_psi * (r - Rmin)**2 )

    def G ( self, a, b ) :
        """
        Utility function for computing xi(r,R). Not for general
        use.
        """
        # shortcuts!
        p0 = self.p0
        p1 = self.p1
        p2 = self.p2
        p3 = self.p3
        
        return -( 2 * ( a * p2 - 1 ) * e**( a * p2 ) - ( a**2 * p2**2 - \
            b**2 * p2**2 + 2 * b * p2 - 2 ) * e**( b * p2 ) ) * p1 *    \
            e**(-p3) ) / p2**2

    def xi ( self, r, R ) :
        """
        The xi function. 
        """
        # shortcuts!
        Rmin  = self.Rmin
        Rmax  = self.Rmax
        p0 = self.p0
        p1 = self.p1
        p2 = self.p2
        p3 = self.p3
        a_psi = self.a_psi
        b_psi = self.b_psi
        c_psi = self.c_psi

        def a ( x ) :
            """
            xi integrand.
            """

            return dpsi(x) / sqrt( -2 * G( r, x ) )

        return integrate.quad(a, Rmin, r, full_output=1)[0] 

    def Z ( self, r, R ) :
        """
        Real position of a flux surface point as a function of 
        flux label and major radius.
        """

        return self.l(R) - sqrt( 1 + self.dldR(R)**2 ) * self.xi( r, R )

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
print "Done."
