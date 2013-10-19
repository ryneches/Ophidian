#!/usr/bin/env python
# vim: set ts=4 sw=4 et:

import scipy
import numpy
from ophidian import *
import pickle
import pylab

class radials :
    
    def __init__( self ) :
        
        self.Rmin   = 0.0
        self.Rmax   = 0.0
        
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
        
        self.opsi   = []
        self.Rh     = []

    def __init__( self, radials_path ) :
        
        a = scipy.io.array_import.read_array( radials_path )
        
        self.R      = a[1:,0]
        self.psi    = a[1:,1]
        self.p      = a[1:,2]
        self.F      = a[1:,3]
        self.q      = a[1:,4]
        self.Jtoro  = a[1:,5]
        self.Btoro  = a[1:,6]
        self.Bpolo  = a[1:,7]
        self.dp     = a[1:,8]
        self.dpsi   = a[1:,9]
        
        del(a)
        
        self.B      = scipy.sqrt( self.Btoro**2 + self.Bpolo**2 )
        self.beta   = (2.0 * 1.25663706e-6 * self.p) / self.B**2
        
        self.opsi   = []
        self.Rh     = []

    def build_psi( self, p ) :
        
        def a( x ) :
            return ( Xi - p.xi_s( x, R ) )**2

        s = blob.spline()
        s.build( self.R, self.psi )

        Rh = []
        for R in self.R :
            
            if R < (r.Rmin + r.Rmax)/2.0 :
                Xi = R - r.Rmin
            else :
                Xi = r.Rmax - R
            print "R:", R, "\nXi:", Xi
            
            if R >= r.Rmin and R <= r.Rmax :
                rr = scipy.optimize.fminbound( a, r.Rmin, R )
            else :
                rr = 0.0
                
            Rh.append( rr )
        
        psi = []
        for i in Rh :
            
            if i >= self.Rmin and i <= self.Rmax :
                psi.append( s.get_Z_for_R(i) )
            else :
                psi.append(0.0)
            
        self.Rh = numpy.array(Rh)
        self.opsi = numpy.array(psi)

p = ophidian.ProblemClient()
p.RunDir = 'ClientRunDir'
p.InitProblemList()

def bload ( problem ) :
    f = open( p.RunDir + '/Problem_' + problem + '/radials.pickle' )
    r = pickle.load( f )
    f.close()
    return r

def bwrite ( problem, r ) :
    f = open( p.RunDir + '/Problem_' + problem + '/radials.pickle', 'w' )
    pickle.dump( r, f )
    f.close()

def bplot( problem ) :
    r = bload( problem )
    for i in range(len(r.R)) :
        if r.R[i] < r.Rmin or r.R[i] > r.Rmax :
            r.psi[i] = 0.0
            r.opsi[i] = 0.0
    pylab.plot( r.R, r.psi, color='blue', label='CUBE' )
    pylab.plot( r.R, r.opsi, color='green', label='Ophidian' )

for problem in p.ProblemList :
    print "problem ::", problem
    #path = 'circular/' + problem
    #r = radials( path + '/bin/results.txt' )
    p.LoadProblem( problem )
    r = bload( problem )
    r.build_psi( p )
    bwrite( problem, r )
   
    #if r.opsi != [] :
        #p.LoadProblem( problem )
        #for i in range(len(r.R)) :
        #    if r.R[i] < p.Problem.Rmin :
        #        r.opsi[i] = 0.0
        #    if r.R[i] > p.Problem.Rmax :
        #        r.opsi[i] = 0.0
        
        #r.Rmin = p.Problem.Rmin
        #r.Rmax = p.Problem.Rmax
        
        #bwrite( problem, r )
        
    pylab.plot( r.R, r.psi, color='black', label=problem )
    pylab.plot( r.R, r.opsi, color='green', label=problem )
    
