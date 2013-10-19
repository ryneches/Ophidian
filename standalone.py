#!/usr/bin/env python2.4
# vim: set ts=4 sw=4 et:

from SimpleXMLRPCServer import *
from ophidian import *
from optparse import OptionParser
import time
import threading
#from scipy import io
#import pickle

parser = OptionParser()
parser.add_option( "-c", "--case", dest="case",     \
                    help="Cases to dump data from", \
                    metavar="CASE" )

(options, args) = parser.parse_args()
case = options.case

class ServerThread( threading.Thread ) :
    
    def __init__( self, case ) :
        self.case = case
        threading.Thread.__init__( self )

        
        print "starting local Ophidian server..."
        self.ProblemServer = ophidian.ProblemServer()
        self.ProblemServer.RunDir = 'RunDir'
        
        # check the run directory for cached problems
        self.ProblemServer.InitProblemList()
        
        # load the first problem, if it exists, exit otherwise
        if self.ProblemServer.ProblemList.__contains__(case) :
         
            self.ProblemServer.ProblemList = [case]
            
            print "Loading CUBE data for problem " + case
            self.ProblemServer.LoadProblem( case )
            
            print "Starting server..."
            self.server = SimpleXMLRPCServer(('127.0.0.1', 8000 ))
            self.server.register_function(self.ProblemServer.GetProblem)
            self.server.register_function(self.ProblemServer.GetBlock)
            self.server.register_function(self.ProblemServer.CheckInBlock)
        else :
            print "No problems found. Ophidian server shutting down."
    
    def run( self ) :
        print "Now serving problem " + self.ProblemServer.Problem.problemname
        self.server.serve_forever()
        
S = ServerThread( case )
S.start()

print "sleeping for 5 seconds to let server come up..."
time.sleep( 5 )

print "Launching Ophidian client..."
ProblemClient = ophidian.ProblemClient()
ProblemClient.clientname = "Your Name Here"
ProblemClient.InitServer( "http://localhost:8000" )
ProblemClient.RunDir = 'ClientRunDir'

#print "Downloading problem data..."
#ProblemClient.GetProblem()

ProblemClient.InitRunDir()
ProblemClient.InitProblemList()

print "Starting calculation..."

ProblemClient.run()
