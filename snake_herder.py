#!/usr/bin/env python
# vim: set ts=4 sw=4 et:

from ophidian import simplesnake
from scipy import *

f = open( 'snake_herder_noshift.txt', 'a' )

L = range( 11, 21 )

for i in L :
    p = simplesnake.Problem()
    p.ImportCUBEdata( 'noshift/cube_case_' + str(i) + '.txt', 'case_' + str(i) )
    rr = p.dPsi_s.R
    ZZ = map( p.psi_radial, rr )
    
    f.write( p.problemname + '\n' )
    f.write( 'R: ' + ' '.join([`num` for num in rr]) + '\n' )
    f.write( 'Z: ' + ' '.join([`num` for num in ZZ]) + '\n' )
    f.flush()
    del(p)
f.close()
