#!/usr/bin/env python
# vim: set ts=4 sw=4 et:
"""
Build a table of data for plotting a comparison graph.
"""

from pylab import *
from Numeric import *
from scipy import *
from optparse import OptionParser
from ACUBE import *
from ACUBE2 import *
from os import *

parser = OptionParser()
parser.add_option("-d", "--directory", dest="directory",
                  help="directory where CUBE data can be found", metavar="DIR")
parser.add_option("-n", "--surfaces", dest="n",
                  help="number of flux surfaces to use in the analytic solution", metavar="N")

(options, args) = parser.parse_args()

picklefile = options.directory+'/output.pickle'
n = int(options.n)

DataHeap = []
for i in listdir(options.directory) :
    if i.find('cube_result') >= 0 :
        DataHeap.append(i)

allplots = len(DataHeap)
plotnum = 1
for directory in DataHeap :

    directory = options.directory+'/'+directory
    print "importing data from", directory, " ..."
    
    # import CUBE data
    print "reading CUBE flux data...."
    a = io.array_import.read_array( directory+'/rerun/bin/flux.txt' )
    NX = len(a)
    NY = len(a[0])
    
    print "    [NX, NY]:", NX, NY
    
    print "reading CUBE radial functions..."
    T = data.CUBEdata()
    T.read( directory+'/rerun/bin/results.txt' )
    
    Rmin, Rmax = T.plasma_range()
    Rmin_i, Rmax_i = T.plasma_range_index()
    axis = T.plasma_axis()
    
    print "    [Rmin, Rmax] =", [Rmin, Rmax]
    print "    axis =", axis
    T.slice( Rmin, axis )
    print "    sliced at Rmin =", Rmin, "axis =", axis
    
    print "fitting CUBE data..."
    
    #### pressure fit
    V = [ 10, 10, 10, 10 ]
    
    def pressure( V, R ) :
        p0, p1, p2, p3 = V
        return p0 + p1 * e**( p2 * R - p3 )
    
    def p_residuals( V, y, x ) :
        return y - pressure( V, x )
    
    fit = optimize.leastsq(p_residuals, V, args=(T.p, T.R), full_output=1)
    if fit[2] == 1 :
        print "    pressure profile fitted to exponential of form p0 + p1 * e**( p2 * R - p3 )"
        print "    [p0, p1, p2, p3] =", fit[0]
        W = fit[0]
    else :
        print "    Error: Pressure profile fit failed."
        print fit[3]
        
    #### flux fit
    #O = [ 1, 1, 1, 1 ]
    O = [ 0.1, 0.1, 0.1, 0.1 ]
    
    def psi ( O, R ) :
        a, b, c, d = O
        return a + b * e**( c * R - d )
    
    def psi_residuals( O, y, x ) :
        #a, b, c, d = O
        return y - psi( O, x )
    
    fit = optimize.leastsq(psi_residuals, O, args=(T.psi, T.R), full_output=1)
    if fit[2] == 1 :
        print "    flux profile fitted to exponential of form a + b * e**( c * R - d )"
        print "    [a, b, c, d] =", fit[0]
        P = fit[0]
    else :
        print "    Error: Flux profile fit failed."
        print fit[3]
    
    #### [d psi]/[d R] fit
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
    
    fit = optimize.leastsq(dpsi_residuals, N, args=(T.dpsi, T.R), full_output=1 )
    if fit[2] == 1 :
        print "    [d psi]/[d R] profile fitted to exponential of form b_psi * c_psi * e**( c_psi * R - d_psi )"
        print "    [b_psi, c_psi, d_psi] =", fit[0]
        M = fit[0]
    else :
        print "    Error: [d psi]/[d R] profile fit failed."
        print fit[3]
    
    
    # build analytic solution
    print "building analytic solution..."
    t = exp_pressure.tokamak()
    t.Rmin = Rmin
    t.Rmax = Rmax
    t.p0 = W[0]
    t.p1 = W[1]
    t.p2 = W[2]
    t.p3 = W[3]
    t.a_psi = P[0]
    t.b_psi = P[1]
    t.c_psi = P[2]
    t.d_psi = P[3]
    t.BtRmin = T.Btoro[Rmin_i]
    
    print "    finding magnetic axis..."
    t.axis = t.axis_label = t.find_axis_label()
    print "    axis found at", t.axis
    
    print "    interpolating to grid...."
    b = array(t.contour( n, NX, NY ))
    c = array(a)
    
    for i in range(NX) :
        for j in range(NY) :
            if b[i][j] == 0.0 :
                c[i][j] = 0.0
    
    d = abs(c - b)
    print "    scanning for numerical clutter..."
    junk = []
    for i in range(NY) :
        if b[i][0] != 0.0 :
            junk.append(i)
    
    print "    clutter found on rows", junk
    
    for i in junk :
        for j in range(NX) :
            d[i][j] = 0.0
    
    max_flux = []
    max_error = []
    for i in a :
        max_flux.append(max(i))
    max_flux = max(max_flux)
    for i in d :
        max_error.append(max(i))
    max_error = max(max_error)
    
    print "    Maximum Flux:", max_flux
    print "    Maximum error:", max_error
    
    print "building contour plots..."
    
    subplot(allplots, 3, plotnum)
    #bone()
    title('CUBE result')
    contourf(a, 50, origin='lower', extent=(3.9,6.1,-1.1,1.1))
    colorbar()
    contour(a, 50, origin='lower', extent=(3.9,6.1,-1.1,1.1), colors='black')
    
    plotnum = plotnum + 1
    
    subplot(allplots, 3, plotnum)
    #bone()
    title('analytic result')
    contourf(b, n + 1, origin='lower', extent=(3.9,6.1,-1.1,1.1))
    #contour(b, n + 1, origin='lower', extent=(3.9,6.1,-1.1,1.1), colors='black')
    colorbar()
    
    plotnum = plotnum + 1
    
    subplot(allplots, 3, plotnum)
    #jet()
    title('normalized max error')
    contourf(d/max_flux, origin='lower', extent=(3.9,6.1,-1.1,1.1))
    colorbar()
    #contour(d, origin='lower', extent=(3.9,6.1,-1.1,1.1), colors='black')
    
    plotnum = plotnum + 1
    
show()
