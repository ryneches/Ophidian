#!/usr/bin/env python
# vim: set ts=4 sw=4 et:

from ophidian import *
from pylab import *
from numpy import *
from scipy import *
from os import access, mkdir, F_OK, R_OK
import pickle
from optparse import OptionParser

parser = OptionParser()
parser.add_option( "-c", "--case", dest="case",     \
                    help="Case to run", \
                    metavar="CASE" )

(options, args) = parser.parse_args()
case = options.case

if not case :
    raise ophidian.DataImportError,     \
        "Please specify a case"

p = ophidian.ProblemClient()
p.RunDir = "ClientRunDir"
p.InitProblemList()

if p.ProblemList.__contains__( case ) :
    print "Loading case", case
    p.LoadProblem(case)
else :
    raise ophidian.DataImportError,     \
        "Case " + case + " cannot be found in " + p.RunDir

bucket_path = 'xi_bucket_' + case + ".pickle"
flux_path   = 'flux_data_' + case + ".pickle"

if access( bucket_path, F_OK ) :
    
    print "Loading interpolated xi functions..."
    f = open( bucket_path )
    bucket = pickle.load( f )
    f.close()

else :

    # if we can't find a pre-calculated xi bucket for
    # this case, build one
    print "Building interpolated xi functions..."
    rr = []
    for i in range(p.Problem.NR) :
        rr.append( p.Problem.r_coor( i ) )
    rr = array(rr)
    
    bucket = []
    for i in rr :
        bucket.append( p.xi_s( rr, i, return_spline='yes' ) )
    
    print "Writing xi_bucket.pickle...."
    f = open( bucket_path, 'w' )
    pickle.dump( bucket, f )
    f.close()

    
print "Calculating problem array..."
data = array(p.Problem.data)
for i in range(p.Problem.NR) :
    for j in range(p.Problem.NZ) :
        if p.Problem.data[i,j] == 1.0 :
            data[i,j] = p.psi_solver(   \
                p.Problem.r_coor(j),    \
                p.Problem.z_coor(i),    \
                xi_bucket=bucket )
            print [i,j], ":", data[i,j]

print "Plotting flux information for inspection..."
contour( data, 16 )
show()

print "Saving flux information..."
f = open( flux_path, 'w' )
pickle.dump( data, f )
f.close()

