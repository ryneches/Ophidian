#!/usr/bin/env python
# vim: set ts=4 sw=4 et:
"""
Utility functions for building objects representing functions from
computed values.
"""

class SplineError(Exception) :
    pass

class spline :
    """
    This is a "natural" spline; the boundary conditions at the beginning
    and end of the spline are such that the second derivative across
    the boundary is zero.
    """

    # data members

    def __init__ (self) :
        self.last_i = 0 
        self.R = []
        self.Z = []
        self.Z2 = []

    # class methods

    def get_R_for_Z (self, Z) :
        """
        find the value R for a given Z
        """
        
        if self.valid() == 0 :
            raise SplineError, "Spline not initialized. Value: " + str(Z) 
        # avoid typing self.foo 
        
        R  = self.R
        Z  = self.Z
        Z2 = self.Z2
        
        # find the right place in the table of R values

        
    def get_Z_for_R (self, r) :
        """
        find the value Z for a given r
        """
        
        # make sure r is treated as a float
        r = float(r)
       
        if self.valid() == 0 :
            raise SplineError,  "Spline not initialized. Value: " + str(r)

        if r < self.R[0] or r > self.R[len(self.R) - 1] :
            raise SplineError, "Called value out of range. Value: " + str(r)
        
        # avoid typing self.foo 
        
        R  = self.R
        Z  = self.Z
        Z2 = self.Z2
        
        # check to see if most recent index is correct

        if r >= R[self.last_i] and r < R[self.last_i + 1] :
            i = self.last_i

        else :
        
            # find the right place in the table of R values
            for i in range(len(R) - 1) :
                if r >= R[i] and r < R[i + 1] :
                    self.last_i = i
                    break
        
        h = R[i + 1] - R[i]

        if h == 0.0 :
            raise SplineError, "Bad spline! Value: ", str(r)
        
        a = (R[i + 1] - r) / h
        b = (r - R[i]) / h
        
        return a * Z[i] + b * Z[i + 1] + \
            ((a**3 - a) * Z2[i] + (b**3 - b) * Z2[i + 1]) * (h**2) / 6.0
        
    def get_dZ_for_R (self, r) :
        """
        find the value Z for a given r
        """
        
        # make sure r is treated as a float
        r = float(r)
       
        if self.valid() == 0 :
            raise SplineError,  "Spline not initialized. Value: " + str(r)


        if r < self.R[0] or r > self.R[len(self.R) - 1] :
            raise SplineError, "Called value out of range. Value: " + str(r)
        
        # avoid typing self.foo 
        
        R  = self.R
        Z  = self.Z
        Z2 = self.Z2
        
        # check to see if most recent index is correct

        if r >= R[self.last_i] and r < R[self.last_i + 1] :
            i = self.last_i

        else :
        
            # find the right place in the table of R values
            for i in range(len(R) - 1) :
                if r >= R[i] and r < R[i + 1] :
                    self.last_i = i
                    break
        
        h = R[i + 1] - R[i]

        if h == 0.0 :
            raise SplineError, "Bad spline! Value: ", str(r)
 
        a = (R[i + 1] - r) / h
        b = (r - R[i]) / h
        
        return - Z[i] / h + Z[i+1] / h + \
            (h/6) * ( ( -3 * a**2 - 1 )*Z2[i] + ( 3 * b**2 - 1 )*Z2[i+1] )

    def build (self, R_new, Z_new, natural = 'yes') :
        """
        Takes list of points in a Cartesian plane and computes the spline
        coefficeints. See Numerical Recipies in C, 2nd edition, section
        3.3 for details.

        Defaults to natural boundary conditions.
        """
    
        if len(R_new) != len(Z_new) :
            raise SplineError, "Dimension Mismatch: R: " + \
                str(R_new) + " Z: " + str(Z_new)
            
        self.R = R_new
        self.Z = Z_new
  
        # create some local names for class members. This is to avoid
        # having to type self.foo all the time

        R  = self.R
        Z  = self.Z
        Z2 = self.Z2

        data = []
        for i in range(len(R)) :
            data.append([ float(R[i]), float(Z[i]) ])

        data.sort()

        for i in range(len(R)) :
            R[i] = data[i][0]
            Z[i] = data[i][1]

        for i in range(len(R) - 1) :
            if R[i] == R[i + 1] :
                raise SplineError, "R values must be unique. Value: " + \
                    str(R[i])
        
        # assorted variables
        
        i = 1
        k = 1
        n = len(self.R) # number of elements
        p = 0.0
        qn = 0.0
        sig = 0.0
        un = 0.0
        u = []
       
        # first derivatives at the boundaries
        Zp1 = 1e5
        Zpn = -1e5
       
        # initialize Z2 and u arrays to the correct size 
        # (values don't matter)

        Z2 = range(n)
        u = range(n)
      
        if natural == 'yes' :
            Z2[0] = u[0] = 0.0  # first real value here is also zero, 
                                # because this is a natural spline
            
            qn = un = Z2[n - 1] = u[n - 1] = 0.0 # "natural" upper boundary
         
        else :
            # set the lower boundary condition to match Zp1
            Z2[0] = -0.5
            u[0] = ( 3.0 / ( R[1] - R[0] ) ) * \
                ( (Z[1] - Z[0] ) / ( R[1] - R[0] ) - Zp1 )
      
            qn = 0.5
            un = ( 3.0 / ( R[n - 1] - R[n - 2] ) ) * \
                ( Zpn - ( Z[n - 1] - Z[n - 2] ) / ( R[n - 1] - R[n - 2] ) )
        
        # tridiagonal decomposition (Z2 used as temporary storage)
        # We loop from the second element to the second to last 
        # element
        for i in range(1, n - 1) :
            sig = ( R[i] - R[i - 1] ) / ( R[i + 1] - R[i - 1] )
            p = sig * Z2[i - 1] + 2.0
            Z2[i] = ( sig - 1.0 ) / p
            u[i] = (6.0 * ((Z[i + 1] - Z[i]) / (R[i + 1] - R[i]) - \
                (Z[i] - Z[i - 1]) / (R[i] - R[i - 1])) / ( R[i + 1] - \
                R[i - 1]) - sig*u[i - 1]) / p
       
        for k in range( n - 2, 0, -1) :
            Z2[k] = Z2[k] * Z2[k + 1] + u[k]
   
        self.Z2 = Z2
        self.last_i = 0

    def valid (self) : 
        """
        Returns true if the spline has been constructed (Z2 is fully
        populated) and the boundary conditions match
        """

        if len( self.Z2 ) == 0 :
            return 0

        #if self.Z2[0] != 0.0 or self.Z2[ len(self.Z2) - 1 ] != 0.0 :
        #    return 0

        Rtest = self.R
        for i in range( len(Rtest) - 1 ) :
            if Rtest[i] == Rtest[i + 1] :
                print "Independant variable cannot have repeated values!"
                print "Shitting the bed."
                return "alkdfjaldsfjadsf"
                
        return 1
