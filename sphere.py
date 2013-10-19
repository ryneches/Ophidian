#!/usr/bin/env python
# vim: set ts=4 sw=4 et:
"""
Simple analytic plotting utility for high beta tokamaks.
"""

from __future__ import division 
from Numeric import *
from scipy import *
from pylab import *

machine  = { 'Rmin': 4.0, 'Rmax': 6.0 }
boundary = { 'al': 1.0, 'bl': 0.2, 'cl': 0.0 }
profile  = 'D'  # D or reverse-D
pressure = { 'p0': 4000.0, 'p2': -1.0, 'p3': 0.0, 'p4': 0.0 }
flux     = { 'af': 0.01, 'bf': 1.0, 'cf': 1.0 }
mu_0     = 1.25663706e-6
#mu_0     = 1.0
step     = 0.01
Raxis     = "unknown"

Rmin = machine['Rmin']
Rmax = machine['Rmax']

def l( R, Rmin = machine['Rmin'], Rmax = machine['Rmax'],   \
    al = boundary['al'],                                    \
    bl = boundary['bl'],                                    \
    cl = boundary['cl'],                                    \
    profile = profile ) :
    """
    For a given radial position (R), compute the vertical position
    (Z) of the plasma boundary.
    
    Rmin and Rmax default to the machine parameters of UCLA's
    Electric Tokamak (4 and 6 meters). The shaping parameters
    (al, bl, and cl) default to a circular profile. The profile
    argument can be either 'D' or 'reverse-D', and defaults to
    'D'.
    
    The only required parameter is R.
    """
    if not Rmin <= R <= Rmax :
        return 0
    
    a = sqrt( R - Rmin ) * sqrt( Rmax - R )
	
    if profile == 'D' :
        return a * sqrt( al / ( R - Rmin + bl * ( Rmax - Rmin )**cl ) )
    
    if profile == 'reverse-D' :
		return a * sqrt( al / ( Rmax - R + bl * ( Rmax - Rmin )**cl ) )

def dldR( R, Rmin = machine['Rmin'], Rmax = machine['Rmax'],    \
    al = boundary['al'],                	                    \
    bl = boundary['bl'],                                        \
    cl = boundary['cl'],                                        \
    profile = profile ) :
    """ 
    For a given radial position (R), compute the derivative of the
    vertical position (Z) of the plasma boundary.
    
    Rmin and Rmax default to the machine parameters of UCLA's
    Electric Tokamak (4 and 6 meters). The shaping parameters
    (al, bl, and cl) default to a circular profile. The profile
    argument can be either 'D' or 'reverse-D', and defaults to
    'D'.
    
    The only required parameter is R.
    """
    if Rmin >= R or R >= Rmax :
        return 0
    
    a = -((al*cl*sqrt(-R + Rmax)*sqrt(R - Rmin)* (R + bl*(Rmax       \
    - Rmin) - Rmin)**(-1 - cl))/ (2*sqrt(al/(R + bl*(Rmax - Rmin)    \
    - Rmin)**cl))) + (sqrt(-R + Rmax)*sqrt(al/(R + bl*(Rmax -        \
    Rmin) - Rmin)** cl))/(2*sqrt(R - Rmin)) - (sqrt(R -              \
    Rmin)*sqrt(al/(R + bl*(Rmax - Rmin) - Rmin)** cl))/(2*sqrt(-R    \
    + Rmax))

    # don't return junk

    if isnan(a) or isinf(a) :
        return 0
    
    return a
		    
def p( R, Rmin = machine['Rmin'], Rmax = machine['Rmax'],       \
    p0 = pressure['p0'],                                        \
    p2 = pressure['p2'],                                        \
    p3 = pressure['p3'],                                        \
    p4 = pressure['p4'] ) :
    """
    For a given radial position (R), compute the pressure (p(R)).
    """
     
    if R < Rmin or R > Rmax :
        return 0

    w = Rmax - Rmin
    r = Rmax - R
    a = p2 * r**2 + p3 * r**3 + p4 * r**4
    b = p2 * w**2 + p3 * w**3 + p4 * w**4
    
    return p0 * ( 1 - ( a / b ) )

def d_psi( Rh, Rmin = machine['Rmin'], Rmax = machine['Rmax'],  \
    af = flux['af'],                                            \
    bf = flux['bf'],                                            \
    cf = flux['cf'] ) :
    """
    For a given flux label (\hat{R}), compute the derivative of
    the flux function psi(\hat{R}). Note that psi(\hat{R}) itself 
    is never calculated, only d(psi)/(d\hat{R}).
    """
    
    if Rh < Rmin or Rh > Rmax :
        return 0

    r = Rh - Rmin
    
    return l(Rh) * ( af + bf * r + cf * r**2 )

def J( R, Rmin = machine['Rmin'], Rmax = machine['Rmax'],          \
    af = flux['af'],                                               \
    bf = flux['bf'],                                               \
    cf = flux['cf'],                                               \
    al = boundary['al'],                                           \
    bl = boundary['bl'],                                           \
    cl = boundary['cl'] ) :
    """
    For a given flux label (\hat{R}), compute the derivative of
    the flux function psi(\hat{R}). Note that psi(\hat{R}) itself 
    is never calculated, only d(psi)/(d\hat{R}).
    """
    
    if R <= Rmin or R >= Rmax :
        return 0

    a = -((al*cl*(af + bf*(R - Rmin) +                         \
              cf*(R - Rmin)**2) *                              \
                      sqrt(R - Rmin) *                         \
                      sqrt(-R + Rmax) *                        \
                          (R - Rmin + bl*                      \
                  (-Rmin + Rmax))**(-1 - cl)) /                \
                  (2*sqrt(al/(R - Rmin +                       \
                  bl*(-Rmin + Rmax))**cl))) -                  \
                  ((af + bf*(R - Rmin) +                       \
                    cf*(R - Rmin)**2)*sqrt(R - Rmin) *         \
                    sqrt(al/(R - Rmin +                        \
                 bl*(-Rmin + Rmax))**cl)) /                    \
                  (2*sqrt(-R + Rmax)) +                        \
              ((af + bf*(R - Rmin) +                           \
            cf*(R - Rmin)**2)*sqrt(-R + Rmax) *                \
                           sqrt(al/(R - Rmin +                 \
                    bl*(-Rmin + Rmax))**cl)) /                 \
              (2*sqrt(R - Rmin)) +                             \
              (bf + 2*cf*(R - Rmin))*sqrt(R - Rmin) *          \
                  sqrt(-R + Rmax) *                            \
      sqrt(al/(R - Rmin + bl*(-Rmin + Rmax))**cl)
    	
    if isnan(a) or isinf(a) :
       	return 0
	
    return ( (-1.0 / R ) * a ) / mu_0

def J_profile( n = 10, err = 2e-4,                          \
    Rmin = machine['Rmin'], Rmax = machine['Rmax'] ) :
    """
    Find 2*n points in the J profile (we locate n on each side of
    the axis).

    err is the termination threshold for the shooting method.
    """

    # find the magnetic axis (redundant!)
    print "Locating magnetic axis..."
    Raxis = findaxis()

    # FIXME : debuging only!
    #Raxis = 5.8389879999999916
    #print "    axis found at", Raxis
    
    # We snap down 10 points between the magnetic axis and the 
    # high field side of the LCFS, and try to find the flux 
    # surfaces that hit them.
    
    print "Beginning J profile calculation for", n, "points..."
    surfaces = []
    points = arange( Raxis, Rmax, (Rmax - Raxis)/n )

    # For our first value of \hat{R}_c, our candidate flux
    # label, we chose Rmin. This should produce a value of 
    # R_c roughly equal to Rmax, which will be good enough to
    # get things going.
    
    Rhc = Rmin + step
    
    # Stupid initial values to make sure we get at least one 
    # iteration

    for Rt in points :
        
        print "Target point:", Rt

        # If the magnitude normalized discrepency vector is 
        # below our termination threshold, halt the problem
        # and continue to the next point.

        F = 1.0
        i = 10
 
        while abs(F) > err and i > 0 :

            # Find the last non-zero point in the curve. We use 
            # this and the point before it to construct a linear 
            # extrapolation of the location of the root. 
            # FIXME : This is dumb. We need a real value for the root.
            
            #a = fluxsurface(Rhc)
            #
            #N = len(a) - 1
            #b = a[N]
            #while b == 0 :
            #    N = N - 1
            #    print [N]
            #    b = a[N]
           
            # our stupid linear interpolation... ick.
            
            #Zc1 = a[N]
            #Rc1 = Rmin + N * step
            #Zc2 = a[N-1]
            #Rc2 = Rmin + (N-1) * step
            
            #Rc = ( Rc2 * Zc1 - Rc1 * Zc2 ) / ( Zc1 - Zc2 )
            
            Rc = HFSfluxroot( Rhc, Raxis=Raxis )
            
            # Having obtained Rc, we compute our normalized 
            # discrepency vector, F.
            
            F = ( Rt - Rc ) / ( Rmax - Raxis ) 
      
            # We use F to compute the next guess for \hat{R}_c
            
            Rhc = Rhc - F * ( Raxis - Rmin )

            if Rhc > Rmax :
                Rhc = Rmax

            if Rhc < Rmin :
                Rhc = Rmin

            # decrement our iteration limiter
            i = i - 1
            
            print "    F:", F," Rc:",Rc," Rhc:",Rhc
           
        # We've finished with the shooting method. Store our candidate
        # position.
        surfaces.append( Rhc )

    return surfaces

def xi ( R, r, Rmin = machine['Rmin'], Rmax = machine['Rmax'],	\
    al = boundary['al'],                                        \
    bl = boundary['bl'],                                        \
    cl = boundary['cl'],                                        \
    p0 = pressure['p0'],                                        \
    p2 = pressure['p2'],                                        \
    p3 = pressure['p3'],                                        \
    p4 = pressure['p4'],                                        \
    af = flux['af'],                                            \
    bf = flux['bf'],                                            \
    cf = flux['cf'],                                            \
    supress_negative = 'yes') :
    """
    In the HAC paper, xi(R, r) is the distance from the last closed
    flux surface to the surface r, perpendicular to the LCFS. R is
    the radial position from the Z axis.
    
    Well, actually, it doesn't. It returns the distance perpendicular
    to the LCFS to to the flux surface r. See fluxsurface() for the
    transform to {R,Z} coordinates.
    
    The guts of this function are going to look horrible no how I stir
    them around. Deal.
    """

    # test some bounds 

    if supress_negative == 'yes' :

        if R == r :
            return 0

        if Rmin > R or R > Rmax :
            return 0

        if Rmin > r or r > Rmax :
            return 0

    # This is horrible. Just go read the paper if you want to
    # understand this. Even then, it probably won't help.
    
    a = (( sqrt(15) *                                         \
      sqrt(-r + Rmax) *                                       \
      (af + bf*(r - Rmin) + cf*(r - Rmin)**2) *               \
      sqrt(r - Rmin) *                                        \
      sqrt(al/(r + bl*(Rmax - Rmin) - Rmin)**cl))  /          \
          sqrt( -(                                            \
    	(p0*(-r + R)**2 *                                     \
    	( -6*p3*(                                             \
    	          3*r**3 +                                    \
    	      6*r**2*R +                                      \
              4*r*R**2 +                                      \
    	      2*R**3                                          \
    	    ) + 45*p3*(r + R)**2 * Rmax -                     \
                30*p3*(r + 2*R)*Rmax**2 +                     \
    	5*p2*( 3*(r + R)**2 -                                 \
                       4*(r + 2*R)*Rmax) +                    \
                2*p4*( 5*(r + R)**2 *                         \
    	       (2*r**2 + R**2) -                              \
                       12*(3*r**3 +                           \
    	       6*r**2*R +                                     \
    	       4*r*R**2 +                                     \
    	       2*R**3) * Rmax +                               \
                       45*(r + R)**2*Rmax**2 -                \
    	       20*(r + 2*R)*Rmax**3)) *                       \
                mu_0))/                                       \
    	(( p2 +                                               \
    	   (p3 + p4*(Rmax - Rmin))*(Rmax - Rmin))*            \
                (Rmax - Rmin)**2)))

    # don't return junk

    if isnan(a) or isinf(a) :
        return 0
    
    return a

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

def fluxsurface ( r, Rmin = machine['Rmin'], Rmax = machine['Rmax'] ) :
    """
    Returns an array of points on a flux surface.
    """
    
    r_points = arange(Rmin, Rmax, step)
    
    def xi_rhat( R ) :
        return xi(R, r)
    
    z = []
    for R in r_points :

        if R >= r :
            a = l(R) - sqrt( 1 + (dldR(R))**2 ) *            \
                integrate.quad(xi_rhat, R, (R + step), full_output=1)[0]
        else :
            a = 0
        
        if a < 0 :
            z.append(0)
        else :
            z.append(a)
    
    return z

def HFSfluxroot ( r, err = 1e-4,                                    \
    Rmin = machine['Rmin'], Rmax = machine['Rmax'],                 \
    Raxis = Raxis ) :
    """
    For a given flux label r, returns the radial position R of the 
    root that falls between the axis and the last closed flux surface.
    """

    # first, check to see if r is nonzero at least *somewhere*

    if sum(fluxsurface( r )) == 0 :
        print "    Flux label", r, "has no roots."
        return Raxis

    # start from the axis

    Rc = Raxis

    def xi_rhat( Rc ) :
        return xi(Rc, r)    
        
    # the scaling factor
    n =  ( l(Raxis) - sqrt( 1 + (dldR(Raxis))**2 ) *            \
         integrate.quad(xi_rhat,                                \
                        Raxis - err,                            \
                        Raxis + err,                            \
                        full_output=1)[0] )
    
    print "scaling factor:", n
    
    i = 0
    a = 187 
    step = 1
    
    while abs(a) > err or Rc > Rmax :
         
        a = l(Rc) - sqrt( 1 + (dldR(Rc))**2 ) *                 \
            integrate.quad(xi_rhat, (Rc - err), (Rc + err), full_output=1)[0]
        
        F = ( n / step ) * a
        
        Rc0 = Rc
        Rc = Rc + F
        
        # bounce off the axis and Rmax so we don't find
        # the wrong root
        
        if Rc < Raxis or Rc > Rmax :
            print "    --> bounce ::"
            print "          F:", F, "Rc:", Rc
            step = step * 2
            Rc = Rc0

        print "   [", i, "] a: ", a, " Rc:", Rc
        
        i = i + 1

    print "    flux surface root found after", i, "iterations"
    return Rc
 
def plotfluxsurface ( r ) :
    """
    Use matplotlib to plot a flux surface with flux label r.
    """
    r_points = arange(Rmin, Rmax, step)
    
    plot(r_points,fluxsurface(r), color="blue")

def plotflux ( n = 40, Rmin = machine['Rmin'], Rmax = machine['Rmax'], \
    Raxis = Raxis ) :
    """
    Plot n flux surfaces. Defaults to 10.
    """
    
    s = ( Rmax - Rmin ) / n
    
    surfaces = arange(Rmin + s, Rmax, 1/n)

    r_points = arange(Rmin, Rmax, step)
    
    print "Plotting flux surfaces..."
    for R in surfaces :
    
        z_points = fluxsurface(R)

        if sum(z_points) == 0 :
            print "    Largest non-null flux label coodinate plotted:", R - s
            break
            
        plot(r_points, z_points, color="blue")

        
def plotfluxHFS( n = 10, Rmin = machine['Rmin'], Rmax = machine['Rmax'], \
    axis = axis ) :
    """
    Kinda like plotflux, except this will choose an even distribution of
    points between the axis and the LCFS.
    """
     
    if axis == "unknown" : 
        print "Locating magnetic axis..."
        axis = findaxis()
        print "    axis located at", axis
    
    points = arange(axis, Rmax, (Rmax - axis) / n)

       # FIXME this function doesn't do anything. 
        
def findaxis ( err=2e-5, Rmin = machine['Rmin'], Rmax = machine['Rmax'] ) :
    """
    Searches for the magnetic axis for n iterations.
    """
    
    # make sure we get at least one loop

    F    = 1
    a    = 1
    Rhc  = Rmin
    Rhc0 = Rmin
    step = 1

    # normalization factor ( kinda hokey )

    n = sum(fluxsurface(Rmin)) / ( Rmax - Rmin )

    while a > err or a == 0.0 :

        Rhc = Rhc + F
        
        s = fluxsurface(Rhc)
        
        a = sum(s)
        
        if a == 0.0 :
            
            step = step * 2
            Rhc = Rhc0
            print "    --> bounce ::"
            print "          F:", F, "R: ", Rhc

        F = n / ( n * step )
        Rhc0 = Rhc

    print "    Axis flux label at:", Rhc

    N = 1
    A = 0
    while A == 0 :
        A = s[len(s) - N]
        N = N + 1
    
    return Rmax - ( ( N + 1 ) / len(s) ) * ( Rmax - Rmin )

def rhat_of_r_guess( R, Rmin = machine['Rmin'], Rmax = machine['Rmax'] ) :
    """
    This is a parabolic approximation of the function rhat(R). We
    use this as our initial guess when solving for the real function.
    """
    return - ( Rmin * ( Rmax - 2*Rmin ) - ( Rmax + Rmin ) * R + R**2 )
