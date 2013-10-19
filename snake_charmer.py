#!/usr/bin/env python
# vim: set ts=4 sw=4 et:

from ophidian import simplesnake
from scipy import *
from pylab import *

f = open( 'snake_herder.txt', 'r' )
p = simplesnake.Problem()

while True :
    problem = f.readline().strip()
    if problem == '' : 
        break
    rr = map( float, f.readline().split(' ')[1:] )
    zz = map( float, f.readline().split(' ')[1:] )
    plot(rr,zz, color='green')
    
    p.ImportCUBEdata( 'circular/'+problem, problem )
    plot(p.Psi_s.R, p.Psi_s.Z, color='black')
    
f.close()

show()
