#!/usr/bin/env python
# vim: set ts=4 sw=4 et:
"""
The functions used to construct the problem (p, l, and
dPsi/d_rhat).
"""

from Numeric import *
from scipy import *
import constants

class queue :
    """
    Simple queue class
    """

    def __init__ (self, n) :
        self.q = []
        self.n = n

    def add (self, x) :
        if len(self.q) >= self.n :
            self.q.reverse()
            self.q.pop()
            self.q.reverse()

        self.q.append(x)

class machine :
    
    # data memebers
    
    def __init__ (self) :
        self.Rmin   = 0.0
        self.Rmax   = 0.0
        self.BtRmin = 0.0

    # class methods

    def valid( self ) :
    
        if self.Rmin <= 0 or self.Rmax <= 0 :
        
            return 0

        if self.Rmin > self.Rmax :
            return 0

        return 1

class plasma(machine) :

    # data members

    def __init__ (self) :
        machine.__init__(self)
        self.al, self.bl, self.cl           = 0.0, 0.0, 0.0
        self.p0, self.p2, self.p3, self.p4  = 0.0, 0.0, 0.0, 0.0
        self.a_psi, self.b_psi, self.c_psi  = 0.0, 0.0, 0.0
        self.axis                           = -1
        self.axis_label                     = -1
        self.root_err                       = 1e-4
        self.axis_err                       = 1e-4        

    # class methods

    ########################################################
    # the constructor functions as define in the HAC paper #
    ########################################################

    def l( self, R ) :
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
        # shortcuts!
        Rmin = self.Rmin
        Rmax = self.Rmax
        al   = self.al
        bl   = self.bl
        cl   = self.cl
        
        if not Rmin <= R <= Rmax :
            return 0
        
        a = sqrt( R - Rmin ) * sqrt( Rmax - R )
	    
        return a * sqrt( al / ( R - Rmin + bl * ( Rmax - Rmin )**cl ) )
        
    def p( self, R ):
        """
        For a given radial position (R), compute the pressure (p(R)).
        """
        # shortcuts!
        Rmin = self.Rmin
        Rmax = self.Rmax
        p0   = self.p0
        p2   = self.p2
        p3   = self.p3
        p4   = self.p4
        
        if R < Rmin or R > Rmax :
            return 0
        
        w = Rmax - Rmin
        r = Rmax - R
        a = p2 * r**2 + p3 * r**3 + p4 * r**4
        b = p2 * w**2 + p3 * w**3 + p4 * w**4
        
        return p0 * ( 1 - ( a / b ) )

    def dp( self, R ) :
        """
        For a given radial position (R), compute the derivative of the
        pressure ([d/dR]p(R)).
        """
        # shortcuts!
        Rmin = self.Rmin
        Rmax = self.Rmax
        p0   = self.p0
        p2   = self.p2
        p3   = self.p3
        p4   = self.p4
        
        if R < Rmin or R > Rmax :
            return 0

        return -((p0*(-2*p2*(Rmax - R) - 3*p3*(Rmax - R)**2 -   \
            4*p4*(Rmax - R)**3))/(p2*(Rmax - Rmin)**2 +         \
            p3*(Rmax - Rmin)**3 + p4*(Rmax - Rmin)**4))

    def dPsi( self, Rh ) :
        """
        For a given flux label (\hat{R}), compute the derivative of
        the flux function psi(\hat{R}). Note that psi(\hat{R}) itself 
        is never calculated, only d(psi)/(d\hat{R}).
        """
        # shortcuts!
        Rmin  = self.Rmin
        Rmax  = self.Rmax
        a_psi = self.a_psi
        b_psi = self.b_psi
        c_psi = self.c_psi
               
        if Rh < Rmin or Rh > Rmax :
            return 0
        
        r = Rh - Rmin
        
        return self.l(Rh) * ( a_psi + b_psi * r + c_psi * r**2 )

    ##########################################################
    # extra functions derived from the constructor functions #
    ##########################################################
    
    def dldR( self, R ) :
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
        # shortcuts!
        Rmin = self.Rmin
        Rmax = self.Rmax
        al   = self.al
        bl   = self.bl
        cl   = self.cl

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

    def xi ( self, R, r ) :
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
        # shortcuts!
        Rmin  = self.Rmin
        Rmax  = self.Rmax
        al    = self.al
        bl    = self.bl
        cl    = self.cl
        p0    = self.p0
        p2    = self.p2
        p3    = self.p3
        p4    = self.p4
        a_psi = self.a_psi
        b_psi = self.b_psi
        c_psi = self.c_psi
        mu_0  = constants.mu_0
      
        # test some bounds 
        
        if R == r :
            return 0
            
        if Rmin > R or R > Rmax :
            return 0
            
        if Rmin > r or r > Rmax :
            return 0
        
        if r > R :
            return 0
        
        # This is horrible. Just go read the paper if you want to
        # understand this. Even then, it probably won't help.
        
        a = (( sqrt(15) *                                           \
            sqrt(-r + Rmax) *                                       \
            (a_psi + b_psi*(r - Rmin) + c_psi*(r - Rmin)**2) *      \
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

    def surface_point ( self, R, r, step = 0.001 ) :
        """
        Returns a point on a given flux surface.
        
        R is the position along the R axis for which you want a 
        corresponding Z value. r is the flux surface label.
        """
        def xi_rhat( R ) :
            return self.xi(R, r)
        
        if R >= r:
            a = self.l(R) - sqrt( 1 + self.dldR(R)**2) *             \
            integrate.quad(xi_rhat, R, (R + step), full_output=1)[0]
        else :
            a = 0

        if a < 0 :
            a = 0

        return a

class plasma_circular (plasma) :
    """
    Extends plasma class, overriding constructor functions with
    their circular corralaries.
    """
    
    def l( self, R ) :
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
        # shortcuts!
        Rmin = self.Rmin
        Rmax = self.Rmax
        al   = self.al
        bl   = self.bl
        cl   = self.cl
        
        if not Rmin <= R <= Rmax :
            return 0
        
        return sqrt(R - Rmin) * sqrt(Rmax - R)

    def dldR( self, R ) :
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
        # shortcuts!
        Rmin = self.Rmin
        Rmax = self.Rmax
        al   = self.al
        bl   = self.bl
        cl   = self.cl

        if Rmin >= R or R >= Rmax :
            return 0
        
        a = -(sqrt(R - Rmin)/(2*sqrt(-R + Rmax))) +     \
             sqrt(-R + Rmax)/(2*sqrt(R - Rmin))
  
        # don't return junk
        
        if isnan(a) or isinf(a) :
            return 0
    
        return a

    def psi ( self, R ) :
        """
        Psi as a function of Rhat. Sorry, but this time R is really Rhat.
        """
        Rmin = self.Rmin
        Rmax = self.Rmax
        a_psi = self.a_psi
        b_psi = self.b_psi
        c_psi = self.c_psi
        
        return (1/192.0)*(sqrt(R - Rmin)*sqrt(-R + Rmax)* (48.0*c_psi   \
            *(R - Rmin)**3.0 + 8.0*(R - Rmin)**2.0*(8.0*b_psi + c_psi   \
            *(Rmin - Rmax)) + 3.0*(16.0*a_psi + (-8.0*b_psi + 5.0*      \
            c_psi*(Rmin - Rmax))*(Rmin - Rmax))* (Rmin - Rmax) + 2.0*   \
            (R - Rmin)*(48.0*a_psi + (Rmin - Rmax)*(8.0*b_psi + 5.0     \
            *c_psi*(-Rmin + Rmax)))) + 3.0*(16.0*a_psi +  (-8.0*b_psi   \
            + 5.0*c_psi*(Rmin - Rmax))*(Rmin - Rmax))* (Rmin -          \
            Rmax)**2.0*arctan(sqrt(R - Rmin)/sqrt(-R + Rmax))) 

    def F ( self, r ) :
        """
        The toroidal flux function as a function of Rhat.
        """
        # shortcuts!
        Rmin   = self.Rmin
        Rmax   = self.Rmax
        BtRmin = self.BtRmin
        p0     = self.p0
        p2     = self.p2
        p3     = self.p3
        p4     = self.p4
        mu_0   = constants.mu_0

        if r >= Rmax or r <= Rmin :
            return 0
        
        return sqrt(2)*sqrt((BtRmin**2*Rmin**2)/2 - (p0* ( r**3 *       \
           (-20*p4*r**3 + 18*r**2*(p3 + 4*p4*Rmax) - 15*r*(p2 + 3 *     \
           Rmax * ( p3 + 2*p4*Rmax)) + 10*Rmax*(2*p2 + Rmax*(3*p3 +     \
           4*p4*Rmax))) - Rmin**3*(-20*p4*Rmin**3 + 18*Rmin**2*(p3 +    \
           4*p4*Rmax) - 15*Rmin*(p2 + 3*Rmax*(p3 + 2*p4*Rmax)) + 10 *   \
           Rmax * ( 2 * p2 + Rmax*(3*p3 + 4 * p4 * Rmax)))) * mu_0 ) /  \
           (30*(Rmin - Rmax)**2*(p2 - (Rmin - Rmax)*(p3 + p4*(-Rmin +   \
           Rmax)))))

    def B ( self, R, r ) :
        """
        The total magnetic field as a function of R and Rhat
        """
        # shortcuts!
        Rmin   = self.Rmin
        Rmax   = self.Rmax
        BtRmin = self.BtRmin
        p0     = self.p0
        p2     = self.p2
        p3     = self.p3
        p4     = self.p4
        mu_0   = constants.mu_0

        if r >= Rmax or r <= Rmin :
            return 0
      
        return sqrt(-(( p0*(r - R)**2 * (-6*p3*(2*r**3 + 4*r**2*R +       \
            6*r*R**2 + 3*R**3) + 45*p3*(r + R)**2 * Rmax - 30*p3*(2*r +   \
            R) * Rmax**2 + 5*p2*(3*(r + R)**2 - 4*(2*r + R)*Rmax) + 2 *   \
            p4 * ( 5 * (r + R)**2*(r**2 + 2*R**2) - 12*(2*r**3 +          \
            4*r**2*R + 6*r*R**2 + 3*R**3)*Rmax + 45*(r + R)**2*Rmax**2 -  \
            20*(2*r + R) * Rmax**3 )) *mu_0) / (15*R**2 * (Rmin -         \
            Rmax)**2 * (p2 - (Rmin - Rmax) * (p3 + p4*(-Rmin + Rmax)))))  \
            + ( 1 / R**2 ) * ( 2 * (( BtRmin** 2 * Rmin**2)/2 - (p0 * (   \
            r**3 * (-20*p4*r**3 + 18*r**2*(p3 + 4*p4*Rmax) - 15*r*(p2 +   \
            3*Rmax*(p3 + 2*p4*Rmax)) + 10*Rmax*(2*p2 + Rmax*(3*p3 +       \
            4*p4*Rmax))) - Rmin**3*(-20*p4*Rmin**3 + 18*Rmin**2*(p3 +     \
            4*p4*Rmax) - 15*Rmin*(p2 + 3*Rmax*(p3 + 2*p4*Rmax)) +         \
            10*Rmax*(2*p2 + Rmax*(3*p3 + 4*p4*Rmax))))*mu_0)/ (30*(Rmin   \
            - Rmax)**2*(p2 - (Rmin - Rmax)*(p3 + p4*(-Rmin + Rmax)))))))
          
    def q ( self, r ) :
        """
        Return q as a function of Rhat
        """
        # shortcuts!
        Rmin   = self.Rmin
        Rmax   = self.Rmax
        BtRmin = self.BtRmin
        al     = self.al
        bl     = self.bl
        cl     = self.cl
        p0     = self.p0
        p2     = self.p2
        p3     = self.p3
        p4     = self.p4
        a_psi  = self.a_psi
        b_psi  = self.b_psi
        c_psi  = self.c_psi
        mu_0   = constants.mu_0

        if r >= Rmax or r <= Rmin :
            return 0

        return (sqrt(2)*sqrt(r - Rmin)*sqrt(-r + Rmax)* sqrt((          \
           BtRmin**2 * Rmin**2 ) / 2 - (p0*(-r + Rmin)**2 * ( -6 * p3   \
           * (3 * r**3 + 6 * r**2 * Rmin + 4*r*Rmin**2 + 2 * Rmin**3 )  \
           + 45 * p3 * (r + Rmin)**2 * Rmax - 30 * p3 * (r + 2 * Rmin   \
           ) * Rmax**2 + 5 * p2 * ( 3 * (r + Rmin)**2 - 4 * (r + 2 *    \
           Rmin ) * Rmax) + 2 * p4 * ( 5 * (r + Rmin )**2 * (2 * r**2   \
           + Rmin**2) - 12*(3*r**3 + 6*r**2*Rmin + 4 * r * Rmin**2 + 2  \
           * Rmin**3 ) * Rmax + 45*(r + Rmin )**2 * Rmax**2 - 20*(r +   \
           2*Rmin ) * Rmax**3))*mu_0)/ ( 30 * (Rmin - Rmax)**2 * (p2 -  \
           (Rmin - Rmax)* (p3 + p4 * (-Rmin + Rmax))))))/ ( pi * r *    \
           abs((a_psi + b_psi*(r - Rmin) + c_psi*(r - Rmin)**2)*        \
           sqrt(r - Rmin)*sqrt(-r + Rmax)))
      
    def beta ( self, R, r ) :
        """
        Plasma beta as a function of R and Rhat.
        """
        # shortcuts!
        Rmin   = self.Rmin
        Rmax   = self.Rmax
        BtRmin = self.BtRmin
        p0     = self.p0
        p2     = self.p2
        p3     = self.p3
        p4     = self.p4
        mu_0   = constants.mu_0

        if r >= Rmax or r <= Rmin :
            return 0

        return (2*p0*(1 - (p2*(-r + Rmax)**2 + p3*(-r + Rmax)**3 +      \
            p4*(-r + Rmax)**4)/ (p2*(-Rmin + Rmax)**2 + p3*(-Rmin +     \
            Rmax)**3 + p4*(-Rmin + Rmax)**4))*mu_0)/ ((BtRmin**2 *      \
            Rmin**2)/r**2 - (p0*(r - Rmin)**2 * (-6*p3*(2*r**3 + 4 *    \
            r**2 * Rmin + 6*r*Rmin**2 + 3*Rmin**3) + 45*p3*(r + Rmin    \
            )**2 * Rmax - 30*p3*(2*r + Rmin)*Rmax**2 + 5*p2*(3*(r +     \
            Rmin)**2 - 4*(2*r + Rmin)*Rmax) + 2*p4*(5*(r + Rmin)**2 *   \
            (r**2 + 2*Rmin**2) - 12*(2*r**3 + 4*r**2*Rmin + 6 * r *     \
            Rmin**2 + 3*Rmin**3)*Rmax + 45*(r + Rmin)**2*Rmax**2 - 20   \
            * ( 2 * r + Rmin)*Rmax**3))*mu_0)/ (15*r**2*(Rmin -         \
            Rmax)**2 * (p2 - (Rmin - Rmax)*(p3 + p4*(-Rmin + Rmax)))))
   
    def xi ( self, R, r ) :
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
        # shortcuts!
        Rmin  = self.Rmin
        Rmax  = self.Rmax
        al    = self.al
        bl    = self.bl
        cl    = self.cl
        p0    = self.p0
        p2    = self.p2
        p3    = self.p3
        p4    = self.p4
        a_psi = self.a_psi
        b_psi = self.b_psi
        c_psi = self.c_psi
        mu_0  = constants.mu_0
      
        # test some bounds 
        
        # I've removed this contraint to allow plotting
        # below the R-axis. This is also necessary for solving
        # for the roots of the flux surface curves.
        #
        #if Rmin > R or R > Rmax :
        #    return 0
            
        if Rmin > r or r > Rmax :
            return 0
        
        if r >= R :
            return 0
       
        # This is the integrand of the xi function. Note that "x" is the
        # dummy variable within the integration (\hat{R}' in the paper by
        # Hsu et al.)
        #
        # Warning: This code was written while the author was being forced
        # to wear a dress.
        def xi_integrand(x) :
            a = (sqrt(15)*sqrt(Rmax - x)*sqrt(-Rmin + x)*(a_psi +       \
            b_psi*(-Rmin + x) + c_psi*(-Rmin + x)**2))/ sqrt(-((p0*(R   \
            - x)**2*(45*p3*Rmax*(R + x)**2 - 30*p3*Rmax**2*(2*R + x) -  \
            6*p3*(2*R**3 + 4*R**2*x + 6*R*x**2 + 3*x**3) + 5*p2*(3*(R   \
            + x)**2 - 4*Rmax*(2*R + x)) + 2*p4*(45*Rmax**2*(R + x)**2   \
            - 20*Rmax**3*(2*R + x) + 5*(R + x)**2*(R**2 + 2*x**2) -     \
            12*Rmax*(2*R**3 + 4*R**2*x + 6*R*x**2 + 3*x**3)))*mu_0)/    \
            ((Rmin - Rmax)**2*(p2 - (Rmin - Rmax)*(p3 + p4*(-Rmin +
            Rmax))))))
            
            # don't return junk
            
            if isnan(a) or isinf(a) :
                return 0
        
            return a
        
        return integrate.quad(xi_integrand, Rmin, r, full_output=1)[0]

           
    def surface_point( self, R, r ) :
        """
        Retern the Z coordinate for a given R of the Rhat'th flux surface
        """

        if R <= r :
            return 0
        
        #a = self.l(R) - sqrt( 1 + self.dldR(R)**2) * self.xi( R, r )
        #
        #if a < 0 :
        #    return 0
        #return a

        # This is the new method of computing Z(R, r) in circular 
        # cases. Note that this is *NOT* how it is described in the HAC
        # paper.
        
        Rmin    = self.Rmin
        Rmax    = self.Rmax
        R0      = ( Rmax + Rmin ) / 2.0
        a       = Rmax - R0
        Xi      = self.xi(R, r)
        
        b = ( R - R0 ) / ( a - Xi )
        
        if b**2 > 1 :
            return 0

        return ( a - Xi ) * sqrt( 1 - b**2 )

    def surface_root( self, r, err ) :
        """
        Find the high-field-side root of the Rhat'th flux surface,
        and return its radial position (R).
        """
        l       = self.l
        dldR    = self.dldR
        xi      = self.xi
        Rmin    = self.Rmin
        Rmax    = self.Rmax
        
        # initial guess
        Rg = Rgnew = Rmin + ( Rmax - Rmin ) / 2
       
        # result queue (for detectic oscilations)
        q = queue(3)
       
        a = err + 1
        b = 1.0
        c = 1.0
        while ( abs(a) > err ) :
            
            if Rgnew < Rmin + ( Rmax - Rmin ) / 4 :
                Rg = c * err + Rmin + ( Rmax - Rmin ) / 2
                c = c + 1
                a = err + 1
                print "--> Noise! Restarting..."
            
            if  q.q.__contains__(Rgnew) :
                Rg = Rg + err
                b = b / 2
                print "--> Attractor!"      
            else :
                Rg = Rgnew
                q.add(Rgnew)
             

            a = l(Rg) - sqrt( 1 + dldR(Rg)**2) * xi( Rg, r )
            
            #if abs(a) > 0.01 :
            #    Rgnew = Rg + ( a / (Rmax - Rmin) )**3
            #else :
            #    Rgnew = Rg + 10**(-abs(1/a))
            
            Rgnew = Rg + ( a / (Rmax-Rmin) )**2
            
            print "[a, R]:", a, Rg

        return Rg
