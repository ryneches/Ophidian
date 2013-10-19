#!/usr/bin/env python
# vim: set ts=4 sw=4 et:
"""
The function used for calculating flux surfaces, and friends.
"""

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

