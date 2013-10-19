#!/usr/bin/env python
# vim: set ts=4 sw=4 et:
"""
Default values for the problem and physical constants.
"""

from scipy import *
from Numeric import *
import blob
import physics

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

class solution_circular(physics.plasma_circular) :
    """
    Problem solution!
    """
    
    def __init__ (self) :
        physics.plasma_circular.__init__(self)
        self.surfaces = []
        self.resolution = -1

    def find_axis( self ) :
        """
        Searches for the magnetic axis.
        """
        err = self.axis_err
        Rmin = self.Rmin
        Rmax = self.Rmax
   
        # make sure we get at least one loop

        F    = 1
        a    = 1
        Rhc  = self.Rmin
        Rhc0 = self.Rmin
        step = 1

        # hokey normalization factor. no physical
        # significance whatsoever.

        n = sum(self.build_surface(Rmin, self.resolution, build='no').spline.Z) / ( Rmax - Rmin )

        while a > err or a == 0.0 :

            Rhc = Rhc + F

            s = self.build_surface(Rhc, self.resolution, build='no').spline.Z
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

        self.axis = Rmax - ( ( N + 1 ) / len(s) ) * ( Rmax - Rmin )
        self.axis_label = Rhc

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
        if self.axis != -1 and r < self.axis_label :
            # find the high field side root
            fs.root = self.surface_root(r, self.root_err)
        
        # note -- we plot only from \hat{R} to Rmax, so the 
        # spline is undefined from Rmin to \hat{R}
        step = (self.Rmax - r) / n
        R = arange(r - step, self.Rmax, step)
        
        Z = []
        for i in R :
            Z.append(self.surface_point(i, r))

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

class CUBEdata :
    """
    A dataset imported from CUBE
    """

    def __init__ (self) :
        self.Rmin   = 0.0
        self.Rmax   = 0.0

        self.i_min  = 0
        self.i_max  = 0

        self.R_raw      = []
        self.psi_raw    = []
        self.p_raw      = []
        self.F_raw      = []
        self.q_raw      = []
        self.Jtoro_raw  = []
        self.Btoro_raw  = []
        self.Bpolo_raw  = []
        self.dp_raw     = []
        self.dpsi_raw   = []

        self.R      = []
        self.psi    = []
        self.p      = []
        self.F      = []
        self.q      = []
        self.Jtoro  = []
        self.Btoro  = []
        self.Bpolo  = []
        self.dp     = []
        self.dpsi   = []


    def read (self, filename) :
       
        a = io.array_import.read_array( filename )
        
        self.R_raw      = a[1:,0]
        self.psi_raw    = a[1:,1]
        self.p_raw      = a[1:,2]
        self.F_raw      = a[1:,3]
        self.q_raw      = a[1:,4]
        self.Jtoro_raw  = a[1:,5]
        self.Btoro_raw  = a[1:,6]
        self.Bpolo_raw  = a[1:,7]
        self.dp_raw     = a[1:,8]
        self.dpsi_raw   = a[1:,9]

    def plasma_range ( self ) :
        """
        Return the Rmin and Rmax values of the plasma according to
        the CUBE boundary conditions that p = const (zero in all cases 
        we care about) on the LCFS.
        """
        
        i = -1
        j = 2
        
        p_edge = self.p_raw[0]
        
        for x in self.p_raw :
            if x == p_edge :
                i = i + 1
            else :
                break
        
        for x in self.p_raw[i+1:] :
            if x != p_edge :
                j = j + 1
            else :
                break

        return (self.R_raw[i], self.R_raw[i+j])

    def plasma_range_index ( self ) :
        """
        Return the Rmin and Rmax values of the plasma according to
        the CUBE boundary conditions that p = const (zero in all cases 
        we care about) on the LCFS.
        """
        
        i = -1
        j = 2
        
        p_edge = self.p_raw[0]
        
        for x in self.p_raw :
            if x == p_edge :
                i = i + 1
            else :
                break
        
        for x in self.p_raw[i+1:] :
            if x != p_edge :
                j = j + 1
            else :
                break

        return (i, i+j)

    def plasma_axis ( self ) :
        """
        Returns the position of the magnetic axis by finding
        the maximum of p(R).
        """
        
        pmax = max(self.p_raw)
        
        i = 0
        while self.p_raw[i] < pmax :
            i = i + 1
        
        return self.R_raw[i+1]

    def slice (self, Rmin, Rmax) :
        
        self.Rmin = Rmin
        self.Rmax = Rmax 

        # find index of first R value larger than Rmin
        i_min = 0
        for r in self.R_raw :
            if r >= Rmin :
                break
            i_min = i_min + 1
        
        # find index of last R value smaller than Rmax
        i_max = 0
        for r in self.R_raw :
            if r >= Rmax :
                break
            i_max = i_max + 1
        i_max = i_max - 1

        self.i_min = i_min
        self.i_max = i_max
        
        self.R      = self.R_raw[i_min:i_max+1]
        self.psi    = self.psi_raw[i_min:i_max+1]  
        self.p      = self.p_raw[i_min:i_max+1]
        self.F      = self.F_raw[i_min:i_max+1]
        self.q      = self.q_raw[i_min:i_max+1]
        self.Jtoro  = self.Jtoro_raw[i_min:i_max+1]
        self.Btoro  = self.Btoro_raw[i_min:i_max+1]
        self.Bpolo  = self.Bpolo_raw[i_min:i_max+1]
        self.dp     = self.dp_raw[i_min:i_max+1]
        self.dpsi   = self.dpsi_raw[i_min:i_max+1]

