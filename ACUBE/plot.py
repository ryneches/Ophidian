#!/usr/bin/env python
# vim: set ts=4 sw=4 et:
"""
Assorted utility functions for plotting.
"""


def plot_diagnostics(  Rmin = machine['Rmin'], Rmax = machine['Rmax'] ) :
    """
    Plot the constructor functions.
    """
    
    print "Global resolution :", (Rmax - Rmin) / step
    
    r_points = arange(Rmin, Rmax, step)
    Z = []
  
    # Last Closed Flux Surface
    subplot(231)
    title("Profile diagnostics")
    xlabel("R")
    ylabel("LCFS")
    plot(r_points, l(r_points))
    
    # Derivative of LCFS
    subplot(232)
    xlabel("R")
    ylabel("dldR")
    
    for R in r_points :
        Z.append(dldR(R))
    
    plot(r_points, Z)

    # d psi / d r_hat 
    subplot(233)
    xlabel("R")
    ylabel("d psi / dR_hat")
    plot(r_points, d_psi(r_points))

    # J_phi
    subplot(234)
    xlabel("R")
    ylabel("J_phi")

    Z=[]
    for R in r_points :
    	Z.append(J(R))
    
    plot(r_points, Z, color="blue")

    # Presssure
    subplot(235)
    xlabel("R")
    ylabel("Pressure")
    plot(r_points, p(r_points))

    # Flux surfaces
    subplot(236)
     
    # find the magnetic axis and snap down a vertical line
    
    print "Locating magnetic axis..."
    axis = findaxis()
    print "    axis located at", axis
    axvline(x=axis, linewidth=2, color='r')
    plotL()
    plotflux()

def plotL ( Rmin = machine['Rmin'], Rmax = machine['Rmax'] ) :
   
    r_points = arange(Rmin, Rmax, step)
    
    z = l(r_points)
    plot(r_points, z, linewidth=1, color="black")


