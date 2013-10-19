#!/usr/bin/env python2.4
# vim: set ts=4 sw=4 et:

from scipy import *
from ophidian import *

print "Loading problem data..."
ProblemServer = ophidian.ProblemServer()
ProblemServer.RunDir = 'WatchDir'
ProblemServer.InitProblemList()

for p in ProblemServer.ProblemList :
    
    print "Problem :", p
    
    ProblemServer.LoadProblem( p )
    
    Clients = {}
    
    for w in ProblemServer.doneblocks :
        
        if Clients.has_key(w.client) :
            Clients[w.client] += 1
        else :
            Clients[w.client] = 1

    print "Standings :"
    print Clients
    print "Problem status :"
    print 100 * (float(len(ProblemServer.doneblocks))       \
                / float(len(ProblemServer.workblocks) +     \
                len(ProblemServer.outblocks))), "%"

    print "Missing Blocks :"
    print len(ProblemServer.outblocks) - len(ProblemServer.doneblocks)

