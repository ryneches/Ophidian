#!/usr/bin/env python
# vim: set ts=4 sw=4 et:

from Numeric import *
from scipy import *
from pylab import *
from ACUBE import blob
from __future__ import division

mu_0     = 1.25663706e-6

class fluxsurface :
    """
    Container class for a splined flux surface.
    """

    def __init__ (self) :
        self.Rhat = -1
        self.root = -1
        self.spline = blob.spline()
        
    def valid (self) :
        """
        Return 1 if the we contain data that could concievably be a
        viable flux surface, 0 if not.
        """
        if self.spline.valid() == 1 and Rhat != -1 :
            return 1
        else :
            return 0

class tokamak :
    """
    An experimental version of the basic class.
    """

    # data memebers

    def __init__ (self) :
        self.Rmin   = 0.0
        self.Rmax   = 0.0
        self.BtRmin = 0.0

        self.al, self.bl, self.cl                       = 0.0, 0.0, 0.0
        self.p0, self.p1, self.p2, self.p3              = 0.0, 0.0, 0.0, 0.0
        self.a_psi, self.b_psi, self.c_psi, self.d_psi  = 0.0, 0.0, 0.0, 0.0
        self.axis                                       = -1
        self.axis_label                                 = -1
        self.root_err                                   = 1e-4
        self.axis_err                                   = 1e-4

        self.surfaces = []

        # memoization values for psi_solver
        self.R_saved    = 0.0
        self.Rh_saved   = 0.0
        
        # pressure and flux splines...
        self.p_s    = blob.spline()
        self.Psi_s  = blob.spline()

        # ...and their derivatives
        self.dp_s   = blob.spline()
        self.dPsi_s = blob.spline()

        # memoizatized quantities for for G_s
        self.G_s_last_i     = 0
        self.G_s_last_j     = 0
        self.G_s_K          = []
        
        # and the R-Rhat spline for handling the branch cut in G_s
        self.rPsi_s = blob.spline()

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

        if R <= Rmin or R >= Rmax :
            return 0

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

        return p1 * p2 * e**( p2 * R - p3 )

    def psi (self, r) :
        """
        The poloidal flux as a function of flux label (r).
        """
        # shortcuts!
        Rmin  = self.Rmin
        Rmax  = self.Rmax
        a_psi = self.a_psi
        b_psi = self.b_psi
        c_psi = self.c_psi
        d_psi = self.d_psi

        if r < Rmin or r > Rmax :
            return 0

        return a_psi + b_psi * e**( c_psi * r - d_psi )
    
    def F ( self, r ) :
        """
        The F function.
        """
        # shortcuts!
        Rmin    = self.Rmin
        Rmax    = self.Rmax
        p0      = self.p0
        p1      = self.p1
        p2      = self.p2
        p3      = self.p3
        BtRmin  = self.BtRmin

        C = (0.5) * ( BtRmin * Rmin )**2

        M = (e**(p2*r)*p1*(2 + p2*r*(-2 + p2*r)) -          \
            e**(p2*Rmin)*p1*(2 + p2*Rmin*(-2 + p2*Rmin))) / \
            (e**p3*p2**2)
        
        return sqrt( 2.0 * ( C - mu_0 * M ) ) 

    def dF ( self, r ) :
        """
        Derivative of F function.
        """
        # shortcuts!
        Rmin    = self.Rmin
        Rmax    = self.Rmax
        p0      = self.p0
        p1      = self.p1
        p2      = self.p2
        p3      = self.p3
        BtRmin  = self.BtRmin

        C = (0.5) * ( BtRmin * Rmin )**2
        
        return -((0.7071067811865476*mu_0*                      \
            (e**(p2*r)*p1*(p2**2*r + p2*(-2 + p2*r)) +          \
            e**(p2*r)*p1*p2*(2 + p2*r*(-2 + p2*r))*log(e)))/    \
            e**p3/(p2**2*sqrt(C - (1/p2**2)*                    \
            ((mu_0*(e**(p2*r)*p1*(2 + p2*r*(-2 + p2*r)) -       \
            e**(p2*Rmin)*p1*(2 + p2*Rmin*(-2 + p2*Rmin))))/     \
            e**p3))))

    def J ( self, r, R ) :
        """
        On-axis current as a function of radius R and flux label r.
        """
        
        return ( self.dp(r) / self.dpsi(r) ) * ( r**2 / R - R )

    def curlB ( self, r ) :
        """
        On-axis current as a function of radius
        """
        # shortcuts!
        Rmin  = self.Rmin
        Rmax  = self.Rmax
        a_psi = self.a_psi
        b_psi = self.b_psi
        c_psi = self.c_psi
        d_psi = self.d_psi

        return -((e**(r*(c_psi) - (d_psi))*r*(b_psi)*(c_psi)**2 \
            + e**(r*(c_psi) - (d_psi))*(b_psi)*(c_psi))/r)/mu_0

    def q ( self, r ) :
        """
        The safety factor as a function of Rhat.
        """
        if r <= self.Rmin or r >= self.Rmax :
            return self.F(self.Rmin) /  \
                ( pi * self.Rmin * self.dpsi(self.Rmin) )
       
        #return ( self.F(r) * self.l(r) ) / ( pi * r * self.dpsi(r) )
    
        # NOTE: I have redefined q(r) by removing l(r) from the
        # numerator. This is because my definition of dpsi(r) no
        # longer includes l(r).
        return ( self.F(r) ) / ( pi * r * self.dpsi(r) )
    
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
        d_psi = self.d_psi
        
        if r < Rmin or r > Rmax :
            return 0

        #return self.l(r) * ( a_psi + b_psi * (r - Rmin) + \
        #    c_psi * (r - Rmin)**2 )

        return b_psi * c_psi * e**( c_psi * r - d_psi )

    def rhat_of_psi( self, psi ) :
        """
        The flux label as a function of flux.
        """

        return (log((psi - self.a_psi) / self.b_psi ) + \
            self.d_psi ) / self.c_psi

    def beta (self, r) :
        """
        Plasma beta as a function of flux label (r).
        """
        
        return ( 2 * mu_0 * self.p(r) ) / self.modB(r, r)**2

    def modB ( self, r, R ) :
        """
        Mod B as a function of R and Rhat.
        """
        
        return sqrt( ( -2 * mu_0 * self.G(R, r) ) /     \
            R**2 + ( self.F(r) / R**2 ) )

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

        if b < self.Rmin or b > self.Rmax :
            return 0

        #return (-( 2 * ( a * p2 - 1 ) * e**( a * p2 ) - ( a**2 * p2**2 -    \
        #    b**2 * p2**2 + 2 * b * p2 - 2 ) * e**( b * p2 ) ) * p1 *        \
        #    e**(-p3) ) / p2**2

        return (p1*(e**(a*p2)*(2 - 2*a*p2) + e**(b*p2)*(-2 + 2*b*p2 + \
            (a**2 - b**2)*p2**2)))/(e**p3*p2**2)

    def Rbranch_build (self, aa, RRR ) :
        """
        Builds the rPsi_s spline used to handle the branch cut in
        the G_s function.

        aa is a table of flux values across the midplane, and 
        RRR is a table of their radial coordinates.
        """
        
        psi_max = max(aa)
        for i in range(len(aa)) :
            if aa[i] >= psi_max :
                break
        i = i + 1

        psi = blob.spline()
        psi.build(RRR,aa)
        self.rPsi_s.build(aa[:i], RRR[:i])

        RR = []
        Rbranch = []
        psi_min = min(aa[:i])

        for j in range(len(aa)) :
            if j < i :
                Rbranch.append(RRR[j])
                RR.append(RRR[j])
            if j >= i and psi.get_Z_for_R(RRR[j]) >= psi_min :
                Rbranch.append(self.rPsi_s.get_Z_for_R(psi.get_Z_for_R(RRR[j])))
                RR.append(RRR[j])
        
        self.rPsi_s.build(RR, Rbranch)
        
    def G_s ( self, rh1, R, slow='no', plots='no' ) :
        """
        Utility function for computing xi_s(r,R). Not for general
        use.

        We must take care here to remember that our definition for
        dp/dr is a spline, and therefor defined piecewise. However,
        we have computed the value of the G function analytically
        assuming that it is a cubic polynomial. To properly compute
        the G function for our splined definition of dp/dr, we must
        also evaluate the integral piecewise. 
        
        The intervals of integration are like so :
        
            R --> R[i]  
                  R[i] --> R[i+2] ... 
                                      R[i+n] --> rh1
        
        At the beginning of the serise, we have a quantity resulting
        from integration with the lower boundary of R and an upper 
        boundary corresponding to the boundary of the local spline
        interval.
        
        At the end of the series, we have a quanitity resulting 
        from integration from the beginning of the local spline 
        interval and ending with rh1.
        
        Between these two quantities, each element in the series is
        the result of integration bounded by the spline intervals
        themselves. The G function is the sum of these quantities.

        In the larger problem, R is ALWAYS larger tha Rhat. However,
        this function returns results for any R and Rhat between
        Rmin and Rmax. This allows us to define the xi function
        in a way that it will return values below the Z axis. These
        results are physically meaningless, but they are useful for
        improving convergance when locating the roots of Z.
        """
        # shortcuts!
        p_s = self.p_s
        
        Rp = p_s.R
        p  = p_s.Z
        dp = p_s.Z2
        
        K = self.G_s_K
        
        if R < self.Rmin or R > self.Rmax :
            return 0
        
        #if rh1 > self.axis :
        #    rh1 = self.rPsi_s.get_Z_for_R(rh1)
        
        if rh1 < Rp[0] or rh1 > Rp[-1] or R < Rp[0] or R > Rp[-1] :
            raise blob.SplineError,                 \
            "Called value out of range. Value: "    \
            + str(R)
        
        # check our cached value for i; if it's right, use it
        if rh1 >= Rp[self.G_s_last_i] and rh1 < Rp[self.G_s_last_i] :
            
            i = self.G_s_last_i
        
        else :
            
            # if we have a cache miss, it is most likly that
            # we've steped to the next interval
            if self.G_s_last_i + 1 < len(Rp) :
                if rh1 >= Rp[self.G_s_last_i + 1] and \
                    rh1 < Rp[self.G_s_last_i + 1] :
                    i = self.G_s_last_i + 1

            # if not, then find the right place in the 
            # table of R values...
            for i in range(len(Rp) - 1) :
                if rh1 >= Rp[i] and rh1 < Rp[i + 1] :
                    break

            # ...and cache our result.
            self.G_s_last_i = i
        
        # check our cached value for j; if it's right, use it
        if R >= Rp[self.G_s_last_j] and R < Rp[self.G_s_last_j] :
            
            j = self.G_s_last_j

        else :

            # if we have a cache miss, it is most likly that
            # we've steped to the next interval
            if self.G_s_last_j + 1 < len(Rp) :
                if rh1 >= Rp[self.G_s_last_j + 1] and \
                    rh1 < Rp[self.G_s_last_j + 1] :
                    j = self.G_s_last_j + 1

            # if not, find the right place in the
            # table of R values...
            for j in range(len(Rp) - 1) :
                if R >= Rp[j] and R < Rp[j + 1] :
                    break
        
            # ...and cache our result.
            self.G_s_las_j = j

        # For testing only! Dangerously slow!
        if slow == 'yes' :
         
            GG1 = R**2 * (self.p_s.get_Z_for_R(rh1) - \
                self.p_s.get_Z_for_R(R))
        
            def foo ( rh ) :
                return self.dp_s.get_Z_for_R(rh) * rh**2
       
            return GG1 - integrate.quad(foo, R, rh1, full_output=1)[0]
        
        def g( a, b ) :
            """
            Analytic solution to the definite integral of p_s(x) * x,
            where p_s(x) is a spline interval of the pressure function,
            and the bounds of the integral are :
                
                a -> lower bound
                b -> upper bound
            
            NOTE :: a and b must be WITHIN the interval of the spline!
            That is, 
                Rp1 <= a <= Rp2
                
            and 
                
                Rp1 <= b <= Rp2
            
            This function is the result of integration by parts of the
            inner integral in the xi(r,R) function, and applying it to
            the definition of a spline found in Numerical Recipes.
            """
            h = Rp2 - Rp1
            
            AA = ( ( (Rp2 / 2) * b**2 - b**3 / 3) - \
                   ( (Rp2 / 2) * a**2 - a**3 / 3) ) / h
            
            BB = ( (b**3 / 3 - (Rp1 / 2) * b**2) - \
                   (a**3 / 3 - (Rp1 / 2) * a**2) ) / h
       
            CC = (h**2/6) * ( (1/h**3) * ( (1/4) * ( (Rp2 - b)**4*b      \
                                                  - ( Rp2 - a)**4*a )) - \
                 ( (1/20) * ( (Rp2 - b)**5 - (Rp2 - a)**5 ) ) - AA )
           
            DD = (h**2/6) * ( (1/h**3) * ( (1/4) * ( (Rp1 - b)**4*b      \
                                                  - ( Rp1 - a)**4*a )) - \
                 ( (1/20) * ( (b - Rp1)**5 - (a - Rp1)**5 ) ) - BB )
            
            return AA * p1 + BB * p2 + CC * dp1 + DD * dp2 
            
        # if the memoization table for the 
        # inner quantities doesn't exist...
        if K == [] :
           
            # ...then build it

            G = 0
            for k in range(len(Rp) - 1) :
                Rp1 = Rp[k]
                Rp2 = Rp[k + 1]
                p1  = p[k]
                p2  = p[k + 1]
                dp1 = dp[k]
                dp2 = dp[k + 1]
        
                G += g(Rp1, Rp2)
                K.append(G)

        # The residual from the analytic integration by parts.
        # It's invariant over spline intervals, so we put it
        # outside the g(a,b) function.
        LL = (R**2 - rh1**2) * p_s.get_Z_for_R(rh1)

        # If R and rh1 are in the same spline interval,
        # then we've got a very simple task.
        if i == j :

            Rp1 = Rp[i]
            Rp2 = Rp[i + 1]
            p1  = p[i]
            p2  = p[i + 1]
            dp1 = dp[i]
            dp2 = dp[i + 1]
 
            return LL + 2*g(R,rh1)

        # If R and rh1 span spline intervals, then we
        # have to compute the integral as a sum over
        # multiple intervals.
        if R < rh1 :
            
            G = 0.0            
            
            # contribution from the
            # interval containing R
            Rp1 = Rp[j]
            Rp2 = Rp[j + 1]
            p1  = p[j]
            p2  = p[j + 1]
            dp1 = dp[j]
            dp2 = dp[j + 1]
            
            G += g(R, Rp2)
            
            # use pre-computed inner quantities
            G += K[i-1] - K[j]
            
            # contribution from the
            # interval containing rh1
            Rp1 = Rp[i]
            Rp2 = Rp[i + 1]
            p1  = p[i]
            p2  = p[i + 1]
            dp1 = dp[i]
            dp2 = dp[i + 1]
            
            G += g(Rp1, rh1)
            
            return LL + 2*G
            
        else :
            
            G = 0.0
            
            # contribution from the 
            # interval containing R
            Rp1 = Rp[j]
            Rp2 = Rp[j + 1]
            p1  = p[j]
            p2  = p[j + 1]
            dp1 = dp[j]
            dp2 = dp[j + 1]
            
            G += g(R, Rp1)

            # use pre-computed inner quantities
            G += K[i] - K[j-1]
            
            # contribution from the
            # interval containing rh1
            Rp1 = Rp[i]
            Rp2 = Rp[i + 1]
            p1  = p[i]
            p2  = p[i + 1]
            dp1 = dp[i]
            dp2 = dp[i + 1]
     
            G += g(Rp2, rh1)

            return LL + 2*G

    def xi ( self, r, R, limit='none' ) :
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

        if limit != 'none' :
            if abs( r - R ) < limit :
                return 0

        if r >= R :
            return 0

        def a ( x ) :
            """
            xi integrand.
            """
            
            if x <= Rmin or x >= Rmax :
                return 0
            
            if x == R :
                return 0

            b = sqrt( -2 * mu_0 * self.G( R, x ) )

            if b == 0 :
                return 0

            c = self.dpsi(x) / b
            
            if isnan(b) or isinf(b) :
                return 0
            
            return c

        return integrate.quad(a, Rmin, r, full_output=1)[0]

    def Z ( self, r, R ) :
        """
        Real position of a flux surface point as a function of 
        flux label and major radius.
        """
        l       = self.l
        dldR    = self.dldR
        xi      = self.xi_s
        Rmin    = self.Rmin
        Rmax    = self.Rmax


        if R <= r or abs(r - R) < 2e-9 :
            return 0
            
        #a = self.l(R) - sqrt( 1 + self.dldR(R)**2 ) * self.xi( r, R )

        Xi = xi( r, R )
        R0 = ( Rmax + Rmin ) / 2.0
        a  = Rmax - R0
       
        #if ( a - Xi )**2 < 1 :
        #    return 0

        b = ( R - R0 )/( a - Xi )

        if b**2 > 1 :
            return 0
        
        #d = self.l( R + Xi * c ) - Xi * sqrt( 1 - c**2 )
        #d = self.l(R) - Xi * sqrt( 1 - c**2 )

        d = (a - Xi) * sqrt( 1 - b**2 )
        
        return d
 
    def xi_s ( self, r, R ) :
        """
        The xi function, splined. 
        """
        # shortcuts!
        Rmin  = self.Rmin
        Rmax  = self.Rmax

        if r >= R :
            return 0

        def a ( x ) :
            """
            xi integrand.
            """
            
            if x <= Rmin or x >= Rmax :
                return 0
            
            if x == R :
                return 0
            
            b = sqrt( -2 * mu_0 * self.G_s( x, R ) )

            if b == 0 :
                return 0

            c = self.dPsi_s.get_Z_for_R(x) / b
            
            if isnan(b) or isinf(b) :
                return 0
            
            return c

        return integrate.quad(a, Rmin, r, full_output=1)[0]

    def find_axis_label( self ) :
        """
        Return the location in R of the magnetic axis within
        an error of err.
        """
        
        Rmin = self.Rmin
        Rmax = self.Rmax
        xi   = self.xi
        

        def a( r ) :
    
            def b ( R ) : 
                return R - Rmax + xi( r, R, limit=2e-8 )
           
            return integrate.quad( b, r, Rmax, full_output=1 )[0]

        # NOTE: this is the Rh value for the smallest 
        # area flux surface
        return optimize.fsolve( a, Rmin )
        
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

        if r >= self.axis_label :
            return -1
        
        # Low field side root
        #def a( R ) :
        #    if R < r :
        #        return -1
        #    return R - Rmax + xi( r, R )

        def b( R ) :
            if R <= r :
                # return something large
                return 10
            #print r, R
            return R - Rmax + xi( r, R )
        
        HFr = optimize.fsolve(b,r)
        
        LFr = optimize.fsolve(b,Rmax)
        
        # WARNING: Must check against invalid Rhat values (i.e.,
        # values that do not produce non-negative Z values)!
        # For such Rhat values, if you plot a, it will not cross
        # the R-axis where xi is defined.
        return [ HFr, LFr ]

    def build_surface(self, r, n, build='yes') :
        """
        Compute a flux surface with a given \hat{R} (r) and 
        return a fluxsurface object with a resolution of n 
        points.
        """

        fs = fluxsurface()

        # don't forget to set Rhat!
        fs.Rhat = r

        # We only allow surface_root to be invoked
        # if the magnetic axis has been found. Otherwise,
        # the function will choke when we hand it surfaces
        # that don't rise above the R-axis.
        if self.axis_label != -1 and r < self.axis_label :
            # find the high field side root
            fs.root = self.surface_root(r, self.root_err)

        # note -- we plot only from \hat{R} to Rmax, so the 
        # spline is undefined from Rmin to \hat{R}
        step = (self.Rmax - r) / n
        R = arange(r - step, self.Rmax, step)

        Z = []
        for i in R :
            # numerical fuzz
            if i - r < 1e-15 :
                Z.append(0.0)
            else :
                Z.append(self.Z(r, i))

        if build == 'yes' :
            fs.spline.build(R, Z)
        else :
            fs.spline.R = R
            fs.spline.Z = Z

        return fs

    def add_surface(self, r, n) :
        """
        Build a surface of a given flux label \hat{R} and resolution 
        n and append it to the surfaces array.
        """

        self.surfaces.append(self.build_surface(r, n))

    def r_coor( self, n ) :
        """
        Return R position of a coordinate N for imported
        CUBE data.
        """
        if len(self.CUBEdata) > 0 :
            return self.Rin + (self.Rout - self.Rin) * ( float(n) / float(self.NX) )
        else :
            return 0

    def z_coor( self, n ) :
        if len(self.CUBEdata) > 0 :
            return self.Zin + (self.Zout - self.Zin) * ( float(n) / float(self.NY) )
        else :
            return 0

    def import_cube_data ( self, path ) :
        """
        Import CUBE data, perform functional fits, and set
        assorted parameters.
        """

    def contour( self, n, NX, NY ) :
        """
        Builds the countour matrix for the analytic solution and the numerical
        solution. Accepts an argument n for the number of contours.
        """
        
        d = (self.axis - self.Rmin)/n 
        Rhats = arange( self.Rmin, self.axis, d )
        
        #print "reading CUBE data...."
        #a = io.array_import.read_array( '../CUBE_data/cube_results_1/rerun/bin/flux.txt' )
        
        U = [6.1, 1.1]
        L = [3.9, -1.1]
        
        b = []
        for i in range(NX) :
            temp = []
            for j in range(NY) :
                temp.append(0.0)
            b.append(temp)
        b = array(b)
        
        for r in Rhats :
            print "      ", r
            for j in range(NX) :
                R = L[0] + (U[0] - L[0]) * ( ( float(j) ) / ( NY - 1 ) )
                Z = self.Z(r, R)
                for i in range(NY) :
                    if Z >= abs( L[1] + (U[1] - L[1]) * ( ( float(i) ) / ( NX - 1 )) ) :
                        psi = self.psi(r + d)
                        if psi > b[i][j] :
                            b[i][j] = psi
        
        #print "construction plots..."
        #subplot(121)
        #contourf(a, 50, origin='lower', extent=(3.9,6.1,-1.1,1.1))
        #contour(a, 50, origin='lower', extent=(3.9,6.1,-1.1,1.1), colors='black')
        #subplot(122)
        #contourf(b, n + 1, origin='lower', extent=(3.9,6.1,-1.1,1.1))
        #contour(b, n + 1, origin='lower', extent=(3.9,6.1,-1.1,1.1), colors='black')
        #show()
        return b

    def find_lcfs( self, a, Rin, Rout, Zin, Zout ) :
        """
        Return a spline of a splined interpolation of the 
        last closed flux surface of a CUBE solution.
        """
        
        NX = len(a)
        NY = len(a[0])
        
        def r_coor( n ) :
            return Rin + (Rout - Rin) * ( float(n) / float(NX) )
            
        def z_coor( n ) :
            return Zin + (Zout - Zin) * ( float(n) / float(NY) )
        
        R = []
        for i in range(len(a)) :
            R.append( r_coor(i) )
        
        roots_lx_top = []
        roots_ly_top = []
        roots_rx_top = []
        roots_ry_top = []
        
        roots_lx_bot = []
        roots_ly_bot = []
        roots_rx_bot = []
        roots_ry_bot = []
        
        s_top = blob.spline()
        s_bot = blob.spline()
        
        s = blob.spline()
        
        # WARNING! This function receives input by side
        # effect! [ R1, R2, s ]
        def ss( R ) :
            if R > R1 and R < R2 :
                try :
                    Z = s.get_Z_for_R( R )
                except blob.SplineError, e :
                    print e
                    return 1.0
                if type(Z) == float :
                    return Z
                else :
                    return 1.0
            else :
                return 1.0
        
        for i in range(0, len(a)) :
            if max(a[i]) >= 0.0 :
                
                s.build( R, a[i] )
                
                # set bounds to left half
                R1 = Rin
                R2 = Rin + (Rout-Rin)/2.0
                
                lr = optimize.fsolve( ss, Rin  + (0.1*(Rout-Rin)) )
                
                # set bounds to right half
                R1 = Rin + (Rout-Rin)/2.0
                R2 = Rout
                
                rr = optimize.fsolve( ss, Rout - (0.1*(Rout-Rin)) )
                
                if i > int(len(a)/2.0) :
                
                    roots_lx_top.append(lr)
                    roots_ly_top.append(z_coor(i))
                    roots_rx_top.append(rr)
                    roots_ry_top.append(z_coor(i))
                
                if i <= int(len(a)/2.0) :
                    
                    roots_lx_bot.append(lr)
                    roots_ly_bot.append(z_coor(i))
                    roots_rx_bot.append(rr)
                    roots_ry_bot.append(z_coor(i))

        roots_rx_top.reverse()
        roots_ry_top.reverse()
        roots_x_top = roots_lx_top + roots_rx_top
        roots_y_top = roots_ly_top + roots_ry_top
        
        s_top.build(roots_x_top, roots_y_top)
        
        roots_rx_bot.reverse()
        roots_ry_bot.reverse()
        roots_x_bot = roots_lx_bot + roots_rx_bot
        roots_y_bot = roots_ly_bot + roots_ry_bot
        
        s_bot.build(roots_x_bot, roots_y_bot)
        
        return (s_top, s_bot)

    def psi_solver( self, s, R1, Z1, RhGuess='none' ) :
        """
        Return the value of psi for a given point { R1, Z1 }, given
        splined plasma boundary l_s.

        We require two l_s splines; l_s1 in the usual coordinates,
        and l_s2 in coordinates rotated -90 degrees. This is to
        allow us to avoid solving for differences of quantities
        tending toward infinity.
        """
        
        def a( R ) :
            """
            Private function for calculating geometric
            values of Xi.
            """
        
            if R > self.Rmin and R < self.Rmax :
            
                try :
                    if Z1 >= 0.0 :
                        # use the top spline
                        result = s[0].get_Z_for_R(R) - Z1 + \
                            (R - R1) / s[0].get_dZ_for_R(R)
                    if Z1 < 0.0 :
                        # use the bottom spline
                        result = s[1].get_Z_for_R(R) - Z1 + \
                            (R - R1) / s[1].get_dZ_for_R(R)
                    
                except blob.SplineError, e :
                    print e
                    return 1.0
                
                else :
                    return result
            
            else :
                return 1.0
        
        def b( Rh ) :
            """
            Private function for calculating values of
            Rhat by comparing geometric values of Xi
            at (R1, Z1) to values of xi(Rhat,R).
            """
            
            if Rh < self.Rmin :
                return 1.0 
            if Rh == R1 :
                return -1.0
            if Rh >= R1 :
                return -1.0 - Rh + R1
            else :
                return Xi - self.xi_s( Rh, R1 )
       
        a0 = ( self.Rmax - self.Rmin ) / 2.0
        R0 = self.Rmin + a0
               
        #if self.R_saved != 0.0 and \
        #    (( R0 <= self.R_saved and R0 <= R1 ) or \
        #    ( R0 > self.R_saved and R0 > R1 )) :
        
            # if the R value from the last call to psi_solver
            # is in the right quadrant, use it as our guess...
            
        #    RGuess = self.R_saved
            
        #else :
            # ...otherwise, at least make sure our guess
            # is in the right quadrant
            
        #    if R1 >= R0 :
        #        RGuess = R0 + a0 / sqrt( 1 + (Z1/(R0-R1))**2 )
        #    else :
        #        RGuess = R0 - a0 / sqrt( 1 + (Z1/(R0-R1))**2 )
        
        if R1 >= R0 :
            RGuess = R0 + a0 / sqrt( 1 + (Z1/(R0-R1))**2 )
        else :
            RGuess = R0 - a0 / sqrt( 1 + (Z1/(R0-R1))**2 )
 
        R = optimize.fsolve(a, RGuess)
        
        # stash the result for next time
        self.R_saved = R
        
        # Find the location of Z on the plasma boundary
        try :
        
            if Z1 >= 0.0 :
                # use top spline
                Z = s[0].get_Z_for_R(R)
            if Z1 < 0.0 :
                # use bottom spline
                Z = s[1].get_Z_for_R(R)
        
        except blob.SplineError, e :
        
            print "Xi was not calculated correctly for point R: " + str(R1) \
                + " Z: " + str(Z1)
            print e
            return 0.0
            
        
        Xi = sqrt( ( R - R1 )**2 + ( Z - Z1 )**2 )
       
        #return Xi
        
        #if self.Rh_saved > self.Rmin and self.Rh_saved < R1 :
        
            #try the Rhat value from the last call to
            #psi_solver as our initial guess
        
        #    RhGuess = self.Rh_saved
            
        #else :
            
        #    RhGuess = self.Rmin + (R1 - self.Rmin) / 2
       
        if RhGuess == 'none' :
            
            RhGuess = self.Rmin + (self.axis + self.Rmin) / 2.0
        
        # find our local Rhat value
        Rh = optimize.fsolve(b, RhGuess)
        
        v = (b(Rh))**2
        
        if v > 0.1 :
            print "--> Failed once, trying again"
            Rh = optimize.fsolve(b, self.Rmin + 2e-4)
            v = (b(Rh))**2
            
        if v**2 > 0.1 :
            print "--> Failed twice, trying again"
            Rh = optimize.fsolve(b, self.axis - 2e-2)
            v = (b(Rh))**2
            
        if v**2 > 0.1 :
            print "--> Failed thrice, trying again"
            Rh = optimize.fsolve(b, R1 - 2e-4)
            v = (b(Rh))**2
         
        if v**2 > 0.1 :
            print "--> Failed four times, trying again"
            Rh = optimize.fsolve(b, (self.axis+R1)/2.0)
            v = (b(Rh))**2
         
        if v**2 > 0.1 :
            print "Bad point!"
            self.Rh_saved = -1
            return -1
        
        # save our Rhat value for next time
        self.Rh_saved = Rh
     
        return self.Psi_s.get_Z_for_R(Rh)
