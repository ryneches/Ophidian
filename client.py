#!/usr/bin/env python2.4
# vim: set ts=4 sw=4 et:

from ophidian import *

print "Launching Ophidian client..."
ProblemClient = ophidian.ProblemClient()
ProblemClient.clientname = "Your Name Here"
ProblemClient.InitServer( "http://am.orpho.us:8765" )
ProblemClient.RunDir = 'ClientRunDir'

#print "Downloading problem data..."
#ProblemClient.GetProblem()

ProblemClient.InitRunDir()
ProblemClient.InitProblemList()

print "Starting calculation..."

ProblemClient.run()
