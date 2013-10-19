#!/usr/bin/env python
# vim: set ts=4 sw=4 et:

from ACUBE import *
from Numeric import *
from scipy import *
from pylab import *
from __future__ import division

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
        mu_0    = constants.mu_0

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
        mu_0    = constants.mu_0

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
        mu_0  = constants.mu_0

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
        mu_0 = constants.mu_0
        
        return ( 2 * mu_0 * self.p(r) ) / self.modB(r, r)**2

    def modB ( self, r, R ) :
        """
        Mod B as a function of R and Rhat.
        """
        mu_0 = constants.mu_0
        
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
        mu_0 = constants.mu_0

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
        xi      = self.xi
        Rmin    = self.Rmin
        Rmax    = self.Rmax


        if R <= r or abs(r - R) < 2e-9 :
            return 0
            
        #a = self.l(R) - sqrt( 1 + self.dldR(R)**2 ) * self.xi( r, R )

        Xi = self.xi( r, R )
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
                return R - Rmax + xi( r, R )
            
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

        fs = data.fluxsurface()

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
        
        print "buidling analytic solution..."
        for r in Rhats :
            for j in range(NX) :
                R = L[0] + (U[0] - L[0]) * ( ( float(j) ) / ( NY - 1 ))
                Z = self.Z(r, R)
                for i in range(NY) :
                    if Z >= abs( L[1] + (U[1] - L[1]) * ( ( float(i) ) / ( NX - 1 )) ) :
                        psi = self.psi(r + d)
                        if psi > b[i][j] :
                            b[i][j] = psi
        
        print "construction plots..."
        #subplot(121)
        #contourf(a, 50, origin='lower', extent=(3.9,6.1,-1.1,1.1))
        #contour(a, 50, origin='lower', extent=(3.9,6.1,-1.1,1.1), colors='black')
        #subplot(122)
        contourf(b, n + 1, origin='lower', extent=(3.9,6.1,-1.1,1.1))
        contour(b, n + 1, origin='lower', extent=(3.9,6.1,-1.1,1.1), colors='black')
        show()
        return b
        
t = tokamak()
t.Rmin = 4.0964285699999996
t.Rmax = 5.9919642800000004
t.p0 = -1.09829042e+01
t.p1 = 4.39420946e-03
t.p2 = 2.98816778e+00
t.p3 = 3.83557094e+00
t.a_psi = -0.02736385
t.b_psi = 6.66980231
t.c_psi = 0.27757685
t.d_psi = 6.64166538
t.BtRmin = 0.268526591 
t.axis = 5.81517857
t.axis_label = 5.81517857

#print "Finding the magnetic axis..."
#t.axis_label = t.find_axis_label()
print "Axis label found at", t.axis_label

T = data.CUBEdata()
T.read('../CUBE_data/cube_results_1/bin/results.txt')

Rh = arange(t.Rmin, t.axis_label - 0.01, 0.1) 

for i in Rh :
    print "building", i
    t.add_surface( i, 200 )

# plot flux surfaces
subplot(331)

#axes([ 0, t.Rmax+0.1, -0.6, 0.6 ])
for fs in t.surfaces :
    plot(fs.spline.R, fs.spline.Z, color='blue')
    plot(fs.spline.R, [ i * -1.0 for i in fs.spline.Z ], color='blue')
xlabel("R")
ylabel("Z")

# plot R vs. Rhat
subplot(332)
Rhs = []
roots = []
for fs in t.surfaces :
    if sum(fs.spline.Z) != 0    \
        and fs.root[1] >= t.Rmin    \
        and fs.root[1] <= t.Rmax :
        Rhs.append(fs.Rhat)
        roots.append(fs.root[1])

roots.reverse()
Rhs.reverse()
RhRv = []
for i in Rhs :
    RhRv.append(i)
Rhs.reverse()

Xc = t.axis
Yc = RhRv[0]
a = roots[0] - t.axis
b = RhRv[0] - t.axis

#def ellipse ( x ) :
#    return sqrt( b**2 * (1 - (x - Xc)**2 / a**2 ) ) + Yc
#
#X = arange(Xc, Xc + a, a / 10)
#Y = ellipse(X)
#
#X = X.tolist()
#Y = Y.tolist()

#plot( Rhs + X + roots, Rhs + Y + RhRv, 'ro' )
plot( Rhs + roots, Rhs + RhRv, 'ro' )
RvRhatSpline = blob.spline()
#RvRhatSpline.build( Rhs + X + roots, Rhs + Y + RhRv, natural='yes' )
RvRhatSpline.build( Rhs + roots, Rhs + RhRv, natural='yes' )

x = arange(t.Rmin, t.Rmax, 0.01)
x_filt = []
y_filt = []
for i in x :
    z = RvRhatSpline.get_Z_for_R(i)
    if type(z) == float:
        x_filt.append(i)
        y_filt.append(z)

plot(x_filt, y_filt, color='red')
xlabel("R")
ylabel("Rhat")

# plot p(R)

subplot(333)
x = arange(t.Rmin, t.Rmax, 0.01)
y = []
for i in x :
    y.append(t.p(RvRhatSpline.get_Z_for_R(i)))

plot(x,y)
plot(T.R_raw, T.p_raw)

xlabel("R")
ylabel("p")

# plot psi(R)

subplot(334)
x = arange(t.Rmin, t.Rmax, 0.01)
y = []
for i in x :
    y.append(t.psi(RvRhatSpline.get_Z_for_R(i)))

plot(x,y)
plot(T.R_raw, T.psi_raw)

xlabel("R")
ylabel("psi")

# plot F(R)

subplot(335)
x = arange(t.Rmin, t.Rmax, 0.01)
y = []
for i in x :
    y.append(t.F(RvRhatSpline.get_Z_for_R(i)))

plot(x,y)
plot(T.R_raw, T.F_raw)

xlabel("R")
ylabel("F")

# plot q(R)

subplot(336)
x = arange(t.Rmin, t.Rmax, 0.01)
y = []
for i in x :
    y.append(t.q(RvRhatSpline.get_Z_for_R(i)))

plot(x,y)
plot(T.R_raw, T.q_raw)

xlabel("R")
ylabel("q")

# plot beta(R)

subplot(337)
x = arange(t.Rmin, t.Rmax, 0.01)
y = []
for i in x :
    y.append(t.beta(RvRhatSpline.get_Z_for_R(i)))

plot(x,y)
plot(T.R_raw,   \
    ( 2 * constants.mu_0 * T.p_raw ) / abs( T.Bpolo_raw**2 + T.Btoro_raw**2 ))

xlabel("R")
ylabel("beta")

# plot J(R)

subplot(338)

y = []
for i in x :
    y.append( t.J(RvRhatSpline.get_Z_for_R(i), i) )

plot( x, y, color='purple' )

plot(T.R_raw, -T.Jtoro_raw, color="green")

# plot modB(R)

subplot(339)

x = arange(t.Rmin, t.Rmax, 0.01)
for i in Rh :
    y = []
    for j in x :
        y.append(t.modB(i, j))
    plot(x, y, color='red')

xlabel("R")
ylabel("|B|")

show()

#for i in Rh :
#    print "plotting", i
#    y = []
#    for j in x :
#        y.append(t.Z(i,j))
#    plot(x,y, color='blue')
#
#show()

### plots!

subplot(211)
x = arange(t.Rmin, t.Rmax, 0.01)

y = []
for i in x :
    y.append(t.p(RvRhatSpline.get_Z_for_R(i)))

plot(x,y, color='blue', linewidth="1.5", label='pressure (fitted)')
plot(T.R_raw, T.p_raw, 'r--', linewidth='1.5', label='CUBE data')
xlabel(r'$R$')
ylabel(r'$p$')
legend(loc='upper left')

subplot(212)
y = []
for i in x :
    y.append(t.psi(RvRhatSpline.get_Z_for_R(i)))

plot(x,y, color='blue', linewidth='1.5', label='poloidal flux (fitted)')
plot(T.R_raw, T.psi_raw, 'r--', linewidth='1.5', label='CUBE data')
xlabel(r'$R$')
ylabel(r'$\psi$')
legend(loc='upper left')
show()

# flux plot
RR = arange(t.Rmin, t.axis_label, 0.113)
fs = []
for i in RR :
    fs.append(t.build_surface(i, 400))
    print "building", i

for i in fs :
    plot(i.spline.R, i.spline.Z, color='blue', linewidth=1.2)

show()

# dump table
f = open('expy_table.txt', 'w')
RRR = arange(t.Rmin, t.Rmax, 0.001)
f.write('R Rh l dl p dp psi dpsi curlB F dF J q beta\n')
for i in RRR :
    r = RvRhatSpline.get_Z_for_R(i)
    
    R       = str(i)
    Rh      = str(r)
    l       = str(t.l(i))
    dl      = '0.0'  # str(t.dl(i))
    p       = str(t.p(r))
    dp      = str(t.dp(r))
    psi     = str(t.psi(r))
    dpsi    = str(t.dpsi(r))
    curlB   = str(t.curlB(r))
    F       = str(t.F(r))
    dF      = str(t.dF(r))
    J       = str(t.J(r,i))
    q       = str(t.q(r))
    beta    = str(t.beta(r))
    
    f.write(R  + ' ' + Rh  + ' ' + l    + ' ' + dl    + ' ' + p + ' ' + \
            dp + ' ' + psi + ' ' + dpsi + ' ' + curlB + ' ' + F + ' ' + \
            dF + ' ' + J   + ' ' + q    + ' ' + beta  + '\n' )
    

f.close()
