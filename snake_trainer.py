#!/usr/bin/env python
# vim: set ts=4 sw=4 et:

from ophidian import simplesnake
from ophidian import blob
from scipy import *
import pickle

f = open( 'snake_herder_noshift.txt', 'r' )

P = simplesnake.PList()

while True :
    
    p = simplesnake.Problem()

    problem = f.readline().strip()
    if problem == '' : 
        break
    rr = map( float, f.readline().split(' ')[1:] )
    zz = map( float, f.readline().split(' ')[1:] )
    
    p.ImportCUBEdata( 'noshift/cube_' + problem + '.txt', problem )
    
    P.add( p, rr, zz )

f.close()

print "writing to disk...."
f = open( 'snake_poop_noshift.dump', 'w' )
pickle.dump( P, f )
f.close()
