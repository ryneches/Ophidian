#!/usr/bin/env python2.4
# vim: set ts=4 sw=4 et:

#from scipy import *
from SimpleXMLRPCServer import *
from ophidian import *

print "Creating server instance..."
ProblemServer = ophidian.ProblemServer()
ProblemServer.RunDir = 'RunDir'

# check the run directory for cached problems
ProblemServer.InitProblemList()

# load the first problem, if it exists, exit otherwise
if ProblemServer.ProblemList != [] :
    
    print str(len(ProblemServer.ProblemList)) + " problems found :"
    for i in ProblemServer.ProblemList :
        print "   " + i

    print "Loading CUBE data for problem " + ProblemServer.ProblemList[0]
    ProblemServer.LoadProblem( ProblemServer.ProblemList[0] )
    
    print "Starting server..."
    server = SimpleXMLRPCServer(('127.0.0.1', 8000 ))
    server.register_function(ProblemServer.GetProblem)
    server.register_function(ProblemServer.GetBlock)
    server.register_function(ProblemServer.CheckInBlock)
    
    print "Now serving problem " + ProblemServer.Problem.problemname
    server.serve_forever()
else :
    print "No problems found. Ophidian server shutting down."
