#!/usr/bin/env python
# vim: set ts=4 sw=4 et:

import pdb
import blob
import pickle

from scipy import io, optimize, integrate
from numpy import array, isnan, isinf, sqrt
from os import access, mkdir, F_OK, R_OK

mu_0     = 1.25663706e-6

class DataImportError(Exception) :
    pass

class NoSuchDirectory(Exception) :
    pass

class GLimitError(Exception) :
    pass

class BadPointError(Exception) :
    pass

class PList :
    """
    Dumb container class for pickling a list of problems.
    """

    def __init__( self ) :
        
        self.list = []

    def add( self, p, rr, zz ) :
        
        del( p.data )
        p.data = []
        p.Psi_ss.build( rr, zz )
        self.list.append(p)
        
class Problem :

    def __init__ (self) :
        
        self.problemname = ""
        
        # radial quantities
        self.p_s        = blob.spline()
        self.dp_s       = blob.spline()
        self.Psi_s      = blob.spline()
        self.Psi_ss     = blob.spline()
        self.dPsi_s     = blob.spline()
        self.q_s        = blob.spline()
        self.F_s        = blob.spline()
        self.Bt_s       = blob.spline()
        self.Bp_s       = blob.spline()
        self.beta_s     = blob.spline()
        self.betaP_s    = blob.spline()
        self.modB_s     = blob.spline()
        
        # flux boundary curves
        self.l1_s       = blob.spline()
        self.l2_s       = blob.spline()
        
        # scalars
        self.Rmin       = 0.0
        self.Rmax       = 0.0
        self.axis       = 0.0
        self.Rin        = 0.0
        self.Rout       = 0.0
        self.Zin        = 0.0
        self.Zout       = 0.0
        self.NR         = 0
        self.NZ         = 0
        
        self.data       = []
        self.datalist   = []
       
        # memoizatized quantities for for G_s
        self.G_s_last_i     = 0
        self.G_s_last_j     = 0
        self.G_s_K          = []
        
        # mapping function vars
        self.Rh_cached      = 0.0
        self.R_cached       = 0.0

    def r_coor( self, n ) :
        return self.Rin + (self.Rout - self.Rin) * ( float(n) / float(self.NR) )
    
    def z_coor( self, n ) :
        return self.Zin + (self.Zout - self.Zin) * ( float(n) / float(self.NZ) )

    def r_n( self, r ) :
        return int(self.NR * ((r - self.Rin)/(self.Rout - self.Rin)))
    
    def z_n( self, z ) :
        return int(self.NZ * ((z - self.Zin)/(self.Zout - self.Zin)))
    
    def matrix_max( self, array ) :
        q = (0,0)
        for i in range(len(array)) :
            for j in range(len(array[0])) :
                if array[q[0]][q[1]] < array[i][j] :
                    q = (i,j)
        return q
    
    def BuildArray ( self, a ) :
        """
        Build the 2D array of points from the CUBE solution and
        the 1D array of matrix coordinates for those points.
        """
        
        b = array(a)
        for i in range(len(a)) :
            for j in range(len(a[0])) :
                if a[i][j] >= 0.0 :
                    b[i][j] = 1
                    self.datalist.append([i,j,'null'])
                else :
                    b[i][j] = 0
                        
        self.data = b
    
    def edge_points( self, a ) :
        """
        Find the points just inside the edge of the plasma.
        """
        l = []
        for i in range(len(a)) :
            for j in range(len(a[0])) :
                if a[i][j] > 0.0 and            \
                        ( a[i+1][j] < 0.0 or    \
                        a[i-1][j] < 0.0 or      \
                        a[i][j+1] < 0.0 or      \
                        a[i][j-1] < 0.0 ) :
                    l.append([i,j])
        return l
    
    def lcfs_points( self, a, i, j ) :
        """
        Find the interpolated points for the LCFS.
        """
        Rin     = self.Rin
        Rout    = self.Rout
        Zin     = self.Zin
        Zout    = self.Zout
        NR      = self.NR
        NZ      = self.NZ
        r_coor  = self.r_coor
        z_coor  = self.z_coor
        
        p = []
        s = blob.spline()
        
        # offset for guessing LCFS position
        #dR = 1.5 * ( Rout - Rin ) / NR
        #dZ = 1.5 * ( Zout - Zin ) / NZ
        dR = 0.0
        dZ = 0.0
        
        def ss( RR ) :
            """
            Wrapper function for the spline evaluation. This is used
            to catch out-of-bounds guesses from the numerical root
            finder.
            
            Note that this function uses the spline objects, so beware
            of side effects.
            """
            if RR >= s.R[0] and RR <= s.R[-1] :
                try :
                    ZZ = s.get_Z_for_R( RR )
                except blob.SplineError, e :
                    # print e
                    return 1.0
                if type(ZZ) == float :
                    return ZZ
                else :
                    return 1.0
            else :
                return 1.0
        
        def sss( RR ) :
            """
            Unprotected wrapper.
            """
            return s.get_Z_for_R( RR )
        
        # in the R direction
        if a[i, j+1] < 0.0 or a[i, j-1] < 0.0 :
            R = []
            for k in range( len(a[0, :]) ) :
                R.append( r_coor(k) )
            s.build( R, a[i, :] )
            gR = r_coor(j)
            p.append( [ optimize.fsolve( sss, gR ), z_coor(i) ] )
        
        # in the Z direction
        if a[i+1, j] < 0.0 or a[i-1, j] < 0.0 :
            Z = []
            for k in range( len(a[:, 0]) ) :
                Z.append( z_coor(k) )
            s.build( Z, a[:, j] )
            gZ = z_coor(i)
            p.append( [ r_coor(j), optimize.fsolve( sss, gZ ) ] )
            
        # FIXME : handle points with two adjacent roots 
        # along one axis!
        
        return p
       
    def Import( self, a,            \
        pR, pZ,                     \
        dpR, dpZ,                   \
        dPsiR, dPsiZ,               \
        FR, FZ,                     \
        qR, qZ,                     \
        BtR, BtZ,                   \
        BpR, BpZ,                   \
        Rin, Rout, Zin, Zout ) :
        """
        Load the problem data into the right data members.
        """
        
        self.BuildArray(a)
        
        self.NR = len(a[0])
        self.NZ = len(a)
        
        # Set the right variables
        
        self.Rin    = Rin
        self.Rout   = Rout
        self.Zin    = Zin
        self.Zout   = Zout
        
        r_coor = self.r_coor
        z_coor = self.z_coor
        
        ##########################
        # start building LCFS... #
        ##########################
        
        # ...find grid points just inside the LCFS
        edge_p = self.edge_points( a )
        
        lcfs_p = []
        for p in edge_p :
            pp = self.lcfs_points( a, p[0], p[1] )
            for q in pp :
                lcfs_p.append(q)
        
        lcfs_p.sort()
        
        # ...separate positive and negative Z values
        lcfs_p1 = []
        lcfs_p2 = []
        for p in lcfs_p :
            if p[1] >= 0.0 :
                lcfs_p1.append(p)
            else :
                lcfs_p2.append(p)
        
        lcfs_p1 = array(lcfs_p1)
        lcfs_p2 = array(lcfs_p2)
        
        # ...build the LCFS splines
        self.l1_s.build( lcfs_p1[:,0], lcfs_p1[:,1] )
        self.l2_s.build( lcfs_p2[:,0], lcfs_p2[:,1] )
        
        # find the min and max R values        
        self.Rmin = self.l1_s.R[0]
        self.Rmax = self.l2_s.R[-1]
        
        # find the grid location of the magnetic axis
        self.axis = r_coor(self.matrix_max(a)[1])
        
        # Build pressure spline
        self.p_s.build( pR, pZ )
        
        # Build pressure gradient spline
        self.dp_s.build( dpR, dpZ )
        
        # build flux derivative spline
        self.dPsi_s.build( dPsiR, dPsiZ )
        
        # build the toroidal function spline
        self.F_s.build( FR, FZ )
        
        # build the safty function spline
        self.q_s.build( qR, qZ )
        
        # build the toroidal mangnetic component spline
        self.Bt_s.build( BtR, BtZ )
        
        # build the poloidal magnetic component spline
        self.Bp_s.build( BpR, BpZ )

        # Build flux spline
        R = []
        for k in range( len(a[0, :]) ) :
           R.append( r_coor(k) )
        self.Psi_s.build( R, a[len(a)/2] )
        
    def ImportCUBEdata( self, CUBEpath, problemname ) :
        """
        Import a solution from CUBE data and initialize the
        problem information.
        """
        
        self.problemname = problemname
        
        flux_path       = CUBEpath + '/bin/flux.txt'
        radials_path    = CUBEpath + '/bin/results.txt'
        input_path      = CUBEpath + '/input.txt'
        
        # make sure we can read the CUBE result!
        if not access( flux_path, R_OK ) :
            raise DataImportError, 'Cannot read from ' + flux_path
        if not access( radials_path, R_OK ) :
            raise DataImportError, 'Cannot read from ' + radials_path
        if not access( input_path, R_OK ) :
            raise DataImportError, 'Cannot read from ' + input_path
        
        print "   ::> Importing " + CUBEpath + " ..."
        
        # read flux data
        data = io.array_import.read_array( flux_path )
        
        print "         Imported " + str(len(data)) + "x" + \
            str(len(data[0])) + " flux grid"
       
        # For various builds of CUBE, the format of results.txt 
        # is different. Pierre occasionally adds columns to his
        # output. So, we must pay attention to the column labels!
        
        f = open( radials_path, 'r' )
        labels = f.readline().split( "\t" )
        f.close()
        labels.pop() # we don't want the "\n" at the end
        
        R_index     = labels.index("R")
        p_index     = labels.index("p")
        dpsi_index  = labels.index("dpsi")
        F_index     = labels.index("F")
        q_index     = labels.index("q")
        Bt_index    = labels.index("Btoro")
        Bp_index    = labels.index("Bpolo")
        
        # read radial quantities
        a = io.array_import.read_array( radials_path )

        R       = a[1:,R_index]
        pZ      = a[1:,p_index]
        dPsiZ   = a[1:,dpsi_index]
        F       = a[1:,F_index]
        q       = a[1:,q_index]
        Bt      = a[1:,Bt_index]
        Bp      = a[1:,Bp_index]

        # some versions of CUBE don't output dp !
        if labels.__contains__("dp") :
            dpZ = a[1:, labels.index("dp")]
        else :
            p  = blob.spline()
            dp = blob.spline()
            p.build( R, pZ )
            dp.build( R, map( p.get_dZ_for_R, R ) )
            del( p )
            dpZ = dp.Z
        
        print "         Imported radial quantity arrays of length " + \
            str(len(R))
        
        # read grid coordinates
        f = open( input_path, 'r' )
        s = "\n"
        while s != ['']  :
            s = f.readline().split(',')
            if s[0] == 'lower_mesh_point' :
                Rin = float(s[1])
                Zin = float(s[2])
            if s[0] == 'upper_mesh_point' :
                Rout = float(s[1])
                Zout = float(s[2])
        f.close()
        
        print "         Grid boundaries :"
        print "            Rin  : " + str(Rin)
        print "            Rout : " + str(Rout)
        print "            Zin  : " + str(Zin)
        print "            Zout : " + str(Zout)
        
        # initialize the problem
        self.Import( data, R, pZ,   \
            R, dpZ,                 \
            R, dPsiZ,               \
            R, F,                   \
            R, q,                   \
            R, Bt,                  \
            R, Bp,                  \
            Rin, Rout, Zin, Zout )
        
        print "   :: Imported Problem_" + problemname
        
        self.make_beta()
        self.make_betaP()
        self.make_modB()
        
        print "   :: Generated beta and |B| splines"

    def make_beta( self ) :
        """
        Build spline for beta values.
        """
        def beta( r ) :
            return ( 2 * mu_0 * self.p_s.get_Z_for_R( r ) ) /   \
                ( self.Bt_s.get_Z_for_R( r )**2 +               \
                self.Bp_s.get_Z_for_R( r )**2 )
        
        self.beta_s.build( self.p_s.R, map( beta, self.p_s.R ) )
 
    def make_betaP( self ) :
        """
        Build spline for beta values.
        """
        def betaP( r ) :
            betap = ( 2 * mu_0 * self.p_s.get_Z_for_R( r ) ) /  \
                self.Bp_s.get_Z_for_R( r )**2
            if isinf(betap) :
                return 0.0
            else :
                return betap
        
        self.betaP_s.build( self.p_s.R, map( betaP, self.p_s.R ) )
    
    def make_modB( self ) :
        """
        Build spline for |B| values.
        """
        def modB( r ) :
            return sqrt( self.Bt_s.get_Z_for_R( r )**2 +    \
                self.Bp_s.get_Z_for_R( r )**2 )
    
        self.modB_s.build( self.Bt_s.R, map( modB, self.Bt_s.R ) )

    def G_s( self, rh1, R, force_null='no' ) :
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
        p_s = self.dp_s
        
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
        if rh1 >= Rp[self.G_s_last_i] and rh1 < Rp[self.G_s_last_i + 1] :
            
            i = self.G_s_last_i
        
        else :
            
            # if we have a cache miss, it is most likly that
            # we've steped to the next interval
            if self.G_s_last_i + 1 < len(Rp) :
                if rh1 >= Rp[self.G_s_last_i + 1] and \
                    rh1 < Rp[self.G_s_last_i + 2] :
                    i = self.G_s_last_i + 1
            
            # if not, then find the right place in the 
            # table of R values...
            for i in range(len(Rp) - 1) :
                if rh1 >= Rp[i] and rh1 < Rp[i + 1] :
                    break
            
            # ...and cache our result.
            self.G_s_last_i = i
        
        # check our cached value for j; if it's right, use it
        if R >= Rp[self.G_s_last_j] and R < Rp[self.G_s_last_j + 1] :
            
            j = self.G_s_last_j
        
        else :
        
            # if we have a cache miss, it is most likly that
            # we've steped to the next interval
            if self.G_s_last_j + 1 < len(Rp) :
                if rh1 >= Rp[self.G_s_last_j + 1] and \
                    rh1 < Rp[self.G_s_last_j + 2] :
                    j = self.G_s_last_j + 1
            
            # if not, find the right place in the
            # table of R values...
            for j in range(len(Rp) - 1) :
                if R >= Rp[j] and R < Rp[j + 1] :
                    break
            
            # ...and cache our result.
            self.G_s_last_j = j
            
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
        
        # The residual from the analytic integration by parts.  It's
        # invariant over spline intervals, so we put it outside the
        # g(a,b) function. If the force_null condition is true, we
        # ignore this value, using instead the observation that the G
        # function must be zero when Rh=R in order to produce the
        # singularity in the right place.
        
        if force_null == 'yes' :
        
            #if R <= rh1 :
            #    return 0.0
            
            LL = 0.0

        else :
        
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
    
    def G_limit_Rh( self, Rh ) :
        """
        Find the point where the G_s function becomes positive.
        Beyond this point, xi is undefined.
        """
        
        def a( R ) :
            return self.G_s( Rh, R, force_null='yes' )**2
        
        return optimize.fminbound( a, Rh, self.Rmax, full_output=1 )[0]
        
    def G_limit_R( self, R ) :
        """
        Find the point where the G_s function becomes positive.
        Beyond this point, xi is undefined.
        """
        
        def a( Rh ) :
            return self.G_s( Rh, R, force_null='yes' )**2
        
        return optimize.fminbound( a, self.Rmin, R, full_output=1 )[0]
        
    def xi_s( self, r, R, return_spline='no', force='yes' ) :
        """
        The xi function, splined. 
        """
        # shortcuts!
        Rmin  = self.Rmin
        Rmax  = self.Rmax
         
        def a ( x ) :
            """
            xi integrand.
            """
            
            if x <= Rmin or x >= Rmax :
                return 0
            
            if x == R :
                return 0
            
            G = self.G_s( x, R, force_null=force )
            
            if G > 0.0 :
                raise GLimitError
            
            b = ( -2.0 * mu_0 * G )**(0.5)
            
            if b == 0 :
                return 0
            
            c = self.dPsi_s.get_Z_for_R(x) / b
            
            if isnan(b) or isinf(b) :
                return 0
            
            return c
        
        if return_spline == 'no' :
        
            if r >= R or ( R - r ) < 0.0 :
                return 1.0
            
            try :
                
                xi = integrate.quad(a, Rmin, r, full_output=1)[0]
            
            except integrate.quadpack.error, e :
            
                print "       ", e
                print "        Arguments ::> Rhat =", r, "R =", R
                return 0.0
                
            except GLimitError, e :
                
                print "Xi is infinite or undefined at Rh:", r, "R:", R
                print e
                raise e
                
            return xi
             
        if return_spline == 'yes' :
            
            if R < Rmin or R > Rmax :
                return 0
            
            # prune impossible values from Rhat list
            
            rr = []

            for i in r :
                if i <= R               \
                    and (R - i) > 6e-6  \
                    and i >= Rmin       \
                    and i <= Rmax :
                    rr.append(i)
            
            if len(rr) < 2 :
                return 0
            
            xixi = []

            for i in rr :
                
                try :
                    
                    xi = integrate.quad(a, Rmin, i, full_output=1)[0]
                
                except integrate.quadpack.error, e :
                    
                    print "       ", e
                    print "        Arguments ::> Rhat =", r, "R =", R
                    return 0.0
                
                xixi.append(xi)
            
            s = blob.spline()
            s.build( rr, xixi )
        
	    return s
        
    def __xi_match__( self, R, R1, Z1 ) :
        """
        Private function used by xi_geometric for calculating geometric
        values of Xi.
        """
        
        if R > self.Rmin and R < self.Rmax :
            
            try :
                if Z1 >= 0.0 :
                    # use the top spline
                    result = self.l1_s.get_Z_for_R(R) - Z1 + \
                        (R - R1) / self.l1_s.get_dZ_for_R(R)
                if Z1 < 0.0 :
                    # use the bottom spline
                    result = self.l2_s.get_Z_for_R(R) - Z1 + \
                        (R - R1) / self.l2_s.get_dZ_for_R(R)
                    
            except blob.SplineError, e :
                print e
                return 1.0
                
            else :
                return result
            
        else :
            return 1.0
        
    def __Rh_match__(self, Rh, Xi, R1 ) :
        """
        Private function used to calculate Rhat.
        """
        
        try :   
                
            if Rh < self.Rmin :
                result = 1.0 
            if Rh == R1 :
                result = -1.0
            if Rh >= R1 :
                #result = - 0.5 - Rh + R1
                restult = sqrt(Rh)
            else :
                try :
                    result = Xi - self.xi_s( Rh, R1 )
                except GLimitError, e :
                    #result = - 0.5 - Rh + R1 
                    result = (Rh)

            return result**2

        except :
                
            print "Rh:", Rh
                
    def xi_geometric( self, R1, Z1 ) :
        """
        Calculate the value of xi geometrically.
        """
        
        a0 = ( self.Rmax - self.Rmin ) / 2.0
        R0 = self.Rmin + a0
               
        # throw an error if the point is outiside the LCFS
        
        if R1 < self.Rmin or R1 > self.Rmax :
            raise BadPointError, "Bad grid point: " + str(R1) + "," \
                + str(Z1) + " is not inside the plasma boundary."
        
        if Z1 > self.l1_s.get_Z_for_R(R1) or \
                Z1 < self.l2_s.get_Z_for_R(R1) :
            raise BadPointError, "Bad grid point: " + str(R1) + "," \
                + str(Z1) + " is not inside the plasma boundary."
         
	# R represents the position ON THE BOUNDARY from which a
        # perpendicular line segment of length xi will extend to
        # the point { R1, Z1 }
        
        if Z1 == 0.0 :
        
            # the radial case
            if R1 > R0 :
                R = self.Rmax
            else :
                R = self.Rmin
        
        else :
        
            # the general case
            
            def a ( r ) : 
                return self.__xi_match__( r, R1, Z1 )
            
            if R1 >= R0 :
                RGuess = R0 + a0 / ( 1 + (Z1/(R0-R1))**2 )**(0.5)
            else :
                RGuess = R0 - a0 / ( 1 + (Z1/(R0-R1))**2 )**(0.5)
            
            R = optimize.fsolve(a, RGuess, full_output=1)[0]
            
        self.R_cached = R
        
        # Find the location of Z on the plasma boundary
        try :
        
            if Z1 > 0.0 :
                # use top spline
                Z = self.l1_s.get_Z_for_R(R)
            if Z1 < 0.0 :
                # use bottom spline
                Z = self.l2_s.get_Z_for_R(R)
            if Z1 == 0.0 :
                # radial case
                Z = 0.0
        
        except blob.SplineError, e :
        
            print "Xi was not calculated correctly for point R: " + str(R) \
                + " Z: " + str(Z1)
            print e
            return 0.0
 
        return ( ( R - R1 )**2 + ( Z - Z1 )**2 )**(0.5)
    
    def psi_radial( self, R1 ) :
        """
        Return the value of psi at a given position along the
        R-axis.
        """
        
        print self.problemname, "-> R1:", R1
        
        if R1 < self.Rmin or R1 > self.Rmax :
            return 0.0
        
        Z1 = 0.0
        
        Xi = self.xi_geometric( R1, 0.0 )
        
        R_limit = self.G_limit_R( R1 )
        
        def a( r ) :
            return self.__Rh_match__( r, Xi, R1 )
        
        Rh = optimize.fminbound( a, self.Rmin, R_limit, full_output=1 )[0]
        
        return self.Psi_s.get_Z_for_R( Rh )

