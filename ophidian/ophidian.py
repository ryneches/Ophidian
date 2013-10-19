#!/usr/bin/env python
# vim: set ts=4 sw=4 et:

import pdb
import blob
import pickle

from SimpleXMLRPCServer import *
from scipy import io, optimize, integrate
from numpy import array, isnan, isinf
#from Numeric import array, isnam, isinf
from os import access, mkdir, F_OK, R_OK
from time import time

mu_0     = 1.25663706e-6

class WorkBlock :

    def __init__(self) :
        problemname = "none"
        blocknumber = 0
        client      = ""
        done        = False
        points      = []

class BadBlockError(Exception) :
    pass

class ClientError(Exception) :
    pass

class DataImportError(Exception) :
    pass

class NoSuchDirectory(Exception) :
    pass

class Problem :

    def __init__ (self) :
        
        self.problemname = ""
        
        self.p_s        = blob.spline()
        self.Psi_s      = blob.spline()
        self.dPsi_s     = blob.spline()
        self.l1_s       = blob.spline()
        self.l2_s       = blob.spline()
        
        self.Rmin       = 0.0
        self.Rmax       = 0.0
        self.axis       = 0.0
        self.Rin        = 0.0
        self.Rout       = 0.0
        self.Zin        = 0.0
        self.Zout       = 0.0
        self.NR         = 0
        self.NZ         = 0
        
        self.data       = []
        self.datalist   = []
        
    def r_coor( self, n ) :
        return self.Rin + (self.Rout - self.Rin) * ( float(n) / float(self.NR) )
    
    def z_coor( self, n ) :
        return self.Zin + (self.Zout - self.Zin) * ( float(n) / float(self.NZ) )

    def r_n( self, r ) :
        return int(self.NR * ((r - self.Rin)/(self.Rout - self.Rin)))
    
    def z_n( self, z ) :
        return int(self.NZ * ((z - self.Zin)/(self.Zout - self.Zin)))
    
    def matrix_max( self, array ) :
        q = (0,0)
        for i in range(len(array)) :
            for j in range(len(array[0])) :
                if array[q[0]][q[1]] < array[i][j] :
                    q = (i,j)
        return q
    
    def BuildArray ( self, a ) :
        """
        Build the 2D array of points from the CUBE solution and
        the 1D array of matrix coordinates for those points.
        """
        
        b = array(a)
        for i in range(len(a)) :
            for j in range(len(a[0])) :
                if a[i][j] >= 0.0 :
                    b[i][j] = 1
                    self.datalist.append([i,j,'null'])
                else :
                    b[i][j] = 0
                        
        self.data = b
    
    def InitStorage ( self, RunDir ) :
        """
        Write initialized problem to the run directory.
        """
        
        path = RunDir + '/' + 'Problem_' + self.problemname + '/'
        
        # if the storage directory for the problem
        # doesn't exist, create it
        
        if not access( path, F_OK ) :
            mkdir( path )
            
            # This no longer seems to work. Now done where
            # InitStorage is called.
            #f = open( path + 'Problem_' + self.problemname + '.pickle', 'w' )
            #pickle.dump( self, f )
            #f.close()
            
        return path + 'Problem_' + self.problemname + '.pickle'
    
    def edge_points( self, a ) :
        """
        Find the points just inside the edge of the plasma.
        """
        l = []
        for i in range(len(a)) :
            for j in range(len(a[0])) :
                if a[i][j] > 0.0 and            \
                        ( a[i+1][j] < 0.0 or    \
                        a[i-1][j] < 0.0 or      \
                        a[i][j+1] < 0.0 or      \
                        a[i][j-1] < 0.0 ) :
                    l.append([i,j])
        return l
    
    def lcfs_points( self, a, i, j ) :
        """
        Find the interpolated points for the LCFS.
        """
        Rin     = self.Rin
        Rout    = self.Rout
        Zin     = self.Zin
        Zout    = self.Zout
        NR      = self.NR
        NZ      = self.NZ
        r_coor  = self.r_coor
        z_coor  = self.z_coor
        
        p = []
        s = blob.spline()
        
        # offset for guessing LCFS position
        #dR = 1.5 * ( Rout - Rin ) / NR
        #dZ = 1.5 * ( Zout - Zin ) / NZ
        dR = 0.0
        dZ = 0.0
        
        def ss( RR ) :
            """
            Wrapper function for the spline evaluation. This is used
            to catch out-of-bounds guesses from the numerical root
            finder.
            
            Note that this function uses the spline objects, so beware
            of side effects.
            """
            if RR >= s.R[0] and RR <= s.R[-1] :
                try :
                    ZZ = s.get_Z_for_R( RR )
                except blob.SplineError, e :
                    # print e
                    return 1.0
                if type(ZZ) == float :
                    return ZZ
                else :
                    return 1.0
            else :
                return 1.0
        
        def sss( RR ) :
            """
            Unprotected wrapper.
            """
            return s.get_Z_for_R( RR )
        
        # in the R direction
        if a[i, j+1] < 0.0 or a[i, j-1] < 0.0 :
            R = []
            for k in range( len(a[0, :]) ) :
                R.append( r_coor(k) )
            s.build( R, a[i, :] )
            gR = r_coor(j)
            p.append( [ optimize.fsolve( sss, gR ), z_coor(i) ] )
        
        # in the Z direction
        if a[i+1, j] < 0.0 or a[i-1, j] < 0.0 :
            Z = []
            for k in range( len(a[:, 0]) ) :
                Z.append( z_coor(k) )
            s.build( Z, a[:, j] )
            gZ = z_coor(i)
            p.append( [ r_coor(j), optimize.fsolve( sss, gZ ) ] )
            
        # FIXME : handle points with two adjacent roots 
        # along one axis!
        
        return p
        
    def Import( self, a, pR, pZ, dPsiR, dPsiZ, Rin, Rout, Zin, Zout ) :
        """
        Load the problem data into the right data members.
        """
        
        self.BuildArray(a)
        
        self.NR = len(a[0])
        self.NZ = len(a)
        
        # Set the right variables
        
        self.Rin    = Rin
        self.Rout   = Rout
        self.Zin    = Zin
        self.Zout   = Zout
        
        r_coor = self.r_coor
        z_coor = self.z_coor
        
        ##########################
        # start building LCFS... #
        ##########################
        
        # ...find grid points just inside the LCFS
        edge_p = self.edge_points( a )
        
        lcfs_p = []
        for p in edge_p :
            pp = self.lcfs_points( a, p[0], p[1] )
            for q in pp :
                lcfs_p.append(q)
        
        lcfs_p.sort()
        
        # ...separate positive and negative Z values
        lcfs_p1 = []
        lcfs_p2 = []
        for p in lcfs_p :
            if p[1] >= 0.0 :
                lcfs_p1.append(p)
            else :
                lcfs_p2.append(p)
        
        lcfs_p1 = array(lcfs_p1)
        lcfs_p2 = array(lcfs_p2)
        
        # ...build the LCFS splines
        self.l1_s.build( lcfs_p1[:,0], lcfs_p1[:,1] )
        self.l2_s.build( lcfs_p2[:,0], lcfs_p2[:,1] )
        
        # find the min and max R values        
        self.Rmin = self.l1_s.R[0]
        self.Rmax = self.l2_s.R[-1]
        
        # find the grid location of the magnetic axis
        self.axis = r_coor(self.matrix_max(a)[1])
        
        # Build pressure spline
        self.p_s.build( pR, pZ )
        
        # build flux derivative spline
        self.dPsi_s.build( dPsiR, dPsiZ )
        
        # Build flux spline
        R = []
        for k in range( len(a[0, :]) ) :
           R.append( r_coor(k) )
        self.Psi_s.build( R, a[len(a)/2] )
        
class ProblemLoader :
    """
    Base class for loading and saving problem state information.
    """

    def __init__( self ) :
         
        self.Problem = ""

        # number of points in each block (default is 50)
        self.block_size = 50
        
        # Work queues
        self.outblocks  = []
        self.workblocks = []
        self.doneblocks = []
         
    def InitRunDir( self ) :
        """
        Create the run directory if it doesn't exist.
        """
        if not os.access( self.RunDir, F_OK ) :
            mkdir( self.RunDir )
    
    def InitProblemList( self ) :
        """
        Populate the ProblemList from the items found in RunDir
        matching the Problem_<problemname> format.
        """
        
        # FIXME ::> make sure this function works as advertised
        
        for d in os.listdir( self.RunDir ) :
            if os.access( self.RunDir + '/' + d, F_OK ) :
                dd = d.split('Problem_')
                if len(dd) > 1 :
                    self.ProblemList.append(dd[1])
    
    def UnserializeWorkBlocks( self ) :
        """
        Reassemble the 2D array.
        """
        for w in self.doneblocks :
            for p in w.points :
                self.Problem.data[p[0]][p[1]] = p[2]
    
    def BuildBlock( self, blocknumber, n, m ) :
        """
        Create and populate a workblock.
        """
        
        w = WorkBlock()
        w.problemname   = self.Problem.problemname
        w.blocknumber   = blocknumber
        w.points        = self.Problem.datalist[n:m]
        w.client        = ""
        
        return w
        
    def Blockify( self, n ) :
        """
        Carve up the problem space into work blocks, and append
        them to the workblock queue.
        """
        i = 0
        num = 0
        
        for i in range( 0, len(self.Problem.datalist), n ) :
            self.workblocks.append( self.BuildBlock( num, i, i+n ) )
            num = num + 1
        
    def LoadProblem( self, problemname ) :
        """
        Page up a problem from disk.
        """
        
        # Make sure there is a current problem
        if self.Problem != "" :

            # only test if the server is initialized
            if problemname == self.Problem.problemname :
                # problem already loaded
                return
        
            # save the problem's state before proceeding
            self.SaveProblemState()
         
        if self.ProblemList.__contains__( problemname ) :
            
            f = open( self.RunDir + '/Problem_' + problemname   \
                + '/' + 'Problem_' + problemname + '.pickle' )
            
            self.Problem = pickle.load( f )
            f.close()
            
            self.ResumeState( problemname )
            
        else :
            
            raise NoSuchDirectory,                              \
                "Problem isn't in problem list :" + problemname
        
    def SaveProblemState( self ) :
        """
        Save queues.
        """
        
        path = self.RunDir + '/Problem_' + self.Problem.problemname + '/'
        
        f = open( path + self.Problem.problemname + \
                  '_workblocks.pickle', 'w' )
        
        pickle.dump( self.workblocks, f )
        f.close()
        
        f = open( path + self.Problem.problemname + \
                  '_doneblocks.pickle', 'w' )
        
        pickle.dump( self.doneblocks, f )
        f.close()
        
        f = open( path + self.Problem.problemname + \
                  '_outblocks.pickle', 'w' )
        
        pickle.dump( self.outblocks, f )
        f.close()
        
    def ResumeState( self, problemname ) :
        """
        Restore queues from saved state.
        """
        
        path = self.RunDir + '/Problem_' + problemname + '/'
        
        if not access( path, F_OK ) :
            raise NoSuchDirectory, "Received a workblock that has no \
                problem directory : " + problemname
        
        workpath = path + problemname + '_workblocks.pickle' 
        donepath = path + problemname + '_doneblocks.pickle' 
        outpath  = path + problemname + '_outblocks.pickle' 
        
        # If the problem doesn't have associated queues, create
        # them.
        if not ( access( workpath, F_OK )       \
                and access( donepath, F_OK )    \
                and access( outpath, F_OK ) ) :
           
            self.workblocks = []
            self.doneblocks = []
            self.outblocks  = []
           
            # self.Blockify( self.block_size )
            
            return
        
        # If the problem's queues do exist, page them up from 
        # the disk.
        else :
            
            f = open( workpath )
            self.workblocks = pickle.load( f )
            f.close()
            
            f = open( donepath )
            self.doneblocks = pickle.load( f )
            f.close()
            
            f = open( outpath )
            self.outblocks = pickle.load( f )
            f.close()
        
class ProblemServer(ProblemLoader) :
    
    def __init__ ( self ) :
        
        ProblemLoader.__init__( self )
        
        # server information
        self.serverport = 8765
               
        # Run Directory (no trailing '/' !)
        self.RunDir         = ""
        self.ProblemList    = []
        
    def ImportCUBEdata( self, CUBEpath, problemname ) :
        """
        Import a solution from CUBE data, initialize the
        problem information, and create its storage directory.
        """
        
        self.InitRunDir()
        
        if self.Problem != "" :
            self.SaveProblemState()
            self.workblocks = []
            self.doneblocks = []
            self.outblocks  = []

        self.Problem = Problem()
        self.Problem.problemname = problemname
        
        flux_path       = CUBEpath + '/bin/flux.txt'
        radials_    = CUBEpath + '/bin/results.txt'
        input_path      = CUBEpath + '/input.txt'
        
        # make sure we can read the CUBE result!
        if not os.access( flux_path, R_OK ) :
            raise DataImportError, 'Cannot read from ' + flux_path
        if not os.access( radials_path, R_OK ) :
            raise DataImportError, 'Cannot read from ' + radials_path
        if not os.access( input_path, R_OK ) :
            raise DataImportError, 'Cannot read from ' + input_path
        
        print "   ::> Importing " + CUBEpath + " ..."
        
        # read flux data
        data = io.array_import.read_array( flux_path )
        
        print "         Imported " + str(len(data)) + "x" + \
            str(len(data[0])) + " flux grid"
        
        # read radial quantities
        a = io.array_import.read_array( radials_path )
        R = a[1:,0]
        pZ = a[1:,2]
        dPsiZ = a[1:,10]
        
        print "         Imported radial quantity arrays of length " + \
            str(len(R))
        
        # read grid coordinates
        f = open( input_path, 'r' )
        s = "\n"
        while s != ['']  :
            s = f.readline().split(',')
            if s[0] == 'lower_mesh_point' :
                Rin = float(s[1])
                Zin = float(s[2])
            if s[0] == 'upper_mesh_point' :
                Rout = float(s[1])
                Zout = float(s[2])
        f.close()
        
        print "         Grid boundaries :"
        print "            Rin  : " + str(Rin)
        print "            Rout : " + str(Rout)
        print "            Zin  : " + str(Zin)
        print "            Zout : " + str(Zout)
        
        # initialize the problem
        self.Problem.Import( data, R, pZ, R, dPsiZ, Rin, Rout, Zin, Zout )
        
        # write the problem to disk
        #self.Problem.InitStorage( self.RunDir )
        path = self.Problem.InitStorage( self.RunDir )
        f = open( path, 'w' )
        pickle.dump( self.Problem, f )
        f.close()

        # carve up the problem space into blocks
        self.Blockify( self.block_size )
        
        # write the problem queues to disk
        self.SaveProblemState()
        
        # add the problme to the ProblemList
        self.ProblemList.append( problemname )
        
        print "   :: Imported Problem_" + problemname
        
        
    def CheckInBlock( self, ws ) :
        """
        Accept a block from a client.
        """
        
        w = pickle.loads( ws )
        
        if w.problemname != self.Problem.problemname :
            if self.ProblemList.__contains__( w.problemname ) :
                
                # if we receive a block that belongs to a problem
                # that isn't loaded, load that problem, queue the
                # block, and reload the current problem
                
                current_problem = self.Problem.problemname
                
                print "   ### Received wayward block, loading problem"  \
                    + w.problemname
                
                self.LoadProblem( w.problemname )
                self.doneblocks.append(w)
                print "   <== queued wayward block", w.blocknumber,     \
                    "from", w.client
                self.SaveProblemState()
                print "   ### Reloading current problem"
                self.LoadProblem( current_problem )
                return "Thanks!"

            else :
                print "   %%% Received block for unknown problem"
        
        for i in w.points :
            if i[2] == 'null' :
                raise BadBlockError, "Unfinished work block."
        
        # At the moment, we do not remove blocks from the outblocks
        # queue. We can handle that later.
        
        #for i in range(len(self.outblocks)) :
        #    if self.outblocks[i].blocknumber == int(w['blocknumber']) :
        #        self.outblocks.pop(i)
        
        self.doneblocks.append(w)
        
        print "   <== received block", w.blocknumber, "from", w.client
        
        self.SaveProblemState()
        
        return "Thanks!"
        
    def GetBlock( self, clientname ) :
        """
        Serve a block to a client, attach it to the queue of
        outblocks.
        """
        
        if len(self.workblocks) == 0 :
            # The work queue is empty. Check for missing blocks,
            # and requeue them.

            print "   ### Work Queue empty."
            self.ReissueMissingBlocks()
            
            if len(self.workblocks) == 0 :
                # if there are no missing blocks to requeue, switch problems
                
                for p in self.ProblemList :
                    
                    self.LoadProblem( p )
                    print "   ### Loading problem", p
                    
                    if len(self.workblocks) > 0 :
                        # Found an incomplete problem. Get to work.
                        self.ProblemList.insert( 0, self.ProblemList.pop() )
                        break
                
                if len(self.workblocks) == 0 :
                    # All problems seem to be complete. Tell clients
                    # to shut down.
                    print "   ### All problems complete." 
                    print "   ### Sending shutdown message to client",  \
                        clientname
                    return "Shutdown."
            
            else :
                print "   ### ", len(self.workblocks), "blocks re-queued"
        
        w = self.workblocks.pop()
        w.client = clientname
        self.outblocks.append(w)
        
        # FIXME ::> This function will not properly handle problem
        # switching case or shutdown case
        
        print "   ==> served block", w.blocknumber, "to", clientname
        
        return pickle.dumps( w )
        
    def GetProblem( self, problemname, clientname ) :
        """
        Package up the current problem and serve it
        to the client as a pickled string.
        """
        
        path = self.RunDir + '/Problem_' + problemname + '/'
        
        if not access( path, F_OK ) :
            print "   *** Received request for unknown problem : " \
                + problemname
            return "No such problem."
        
        f = open( path + 'Problem_' + problemname + '.pickle' )
        problem = f.read()
        f.close()
    
        print '   ==> served ' + problemname + ' to ' + clientname
        
        return problem
        
    def ContainsBlock( self, queue, blocknumber ) :
        """
        Return True if queue contains a block with a given
        blocknumber.
        """
        
        for w in queue :
    
            if w.blocknumber == blocknumber :
                return True
            
        return False
            
    def ReissueMissingBlocks( self ) :
        """
        Reissue missing blocks.
        """
        
        missingblocks = []
        
        # for each issued block that has not yet been returned,
        # append it to the missing block queue if it has alredy
        # been placed there (blocks can often be issued more than
        # once, and therefor the outblocks queue may contain
        # duplicates)

        for ow in self.outblocks :
            if not self.ContainsBlock( self.doneblocks, ow.blocknumber )    \
                and not self.ContainsBlock( missingblocks, ow.blocknumber ) \
                and not self.ContainsBlock( self.workblocks, ow.blocknumber ) :
                
                missingblocks.append( ow )
    
        self.workblocks = self.workblocks + missingblocks
    
class ProblemClient(ProblemLoader) :
    """
    The Ophidian client object. All the physics is in here. 
    """
    
    def __init__(self) :
        
        ProblemLoader.__init__( self )
        
        # client information
        self.clientname     = ""
        
        # server information
        self.serverURL      = ""
        self.server         = ""
        
        # client information
        self.RunPath        = ""
        self.ProblemList    = []
        
        # memoizatized quantities for for G_s
        self.G_s_last_i     = 0
        self.G_s_last_j     = 0
        self.G_s_K          = []
        
        # current block
        self.w              = ""
        
    def InitServer( self, serverURL ) :
        """
        Initialize the XML-RPC server object.
        """
        
        self.serverURL = serverURL
        
        # initialize XML-RPC
        self.server = xmlrpclib.Server(self.serverURL)
        
    def GetProblem( self, problemname ) :
        """
        Fetch the problem data from the server or load it from disk.
        """
        
        if self.Problem != "" :
            self.Problem.InitStorage( self.RunDir )
            self.SaveProblemState()
        
        if self.ProblemList.__contains__( problemname ) :
            
            # if we already have the problem, load it from disk...
            print "   ### Loading", problemname, "from disk"
            self.LoadProblem( problemname )
            print "   ### Loaded problem", self.Problem.problemname
             
        else :
            
            # ...otherwise, request it from the server
            print "   ### Fetching", problemname, "from server"
            self.Problem = pickle.loads(    \
                self.server.GetProblem( problemname, self.clientname ))
             
            self.workblocks = []
            self.outblocks  = []
            self.doneblocks = []
            
            # initialize the problem's storage directory
            # and pickle it pickle the problem into it
            path = self.Problem.InitStorage(self.RunDir)
            f = open( path, 'w' )
            pickle.dump( self.Problem, f )
            f.close()
            
            print "   <== Downloaded problem", self.Problem.problemname
            
            self.ProblemList.append(self.Problem.problemname)
     
    def DoBlock( self, w ) :
        """
        Execute a work block.
        """
        
        # If there is no current problem, load one
        if self.Problem == "" :
            self.LoadProblem( w.problemname )
        
        # Make sure we have the right problem loaded, if not, load it
        if w.problemname != self.Problem.problemname :
            self.GetProblem( w.problemname )
        
        for k in range(len(w.points)) :
            w.points[k][2] = self.psi_solver(               \
                self.Problem.r_coor(w.points[k][1]),        \
                self.Problem.z_coor(w.points[k][0]) )
        
        return w
        
    def GetBlock( self ) :
        """
        Ask the server for a work block. 
        """
        ws = self.server.GetBlock( self.clientname )
        
        # raise an error if server reports shutdown state
        if ws == 'Shutdown.' :
            raise ClientError, "Server has invoked shutdown"
        
        w = pickle.loads( ws )
        
        # If we don't have a current problem, request it from
        # the server
        if self.Problem == '' :
            self.GetProblem(w.problemname)
        
        # If we don't already have the problem, request it from
        # the server. It will become the current problem.
        if not self.ProblemList.__contains__(w.problemname) :
            self.GetProblem(w.problemname)
        
        print "   <== Downloaded block", w.blocknumber
        
        return w
    
    def run( self ) :
        """
        Run forever!
        """
        haltfile = self.RunDir + "/halt"
        
        while True :
            
            if os.access( haltfile, F_OK ) :
                print "   ### Halt file located :", haltfile
                print "   ### Client shutting down"
                break
                
            try :
                self.w = self.GetBlock()
            except ClientError, e :
                print e
                break
                
            print "   ooo Computing block", self.w.blocknumber
            
            t1 = time()
            self.w = self.DoBlock( self.w )
            t2 = time()
            
            print "   ooo Block %s completed in %1.2f seconds" %    \
                (self.w.blocknumber, t2 - t1 )
            
            print "   ==> Uploading block", self.w.blocknumber
            
            try :
                self.server.CheckInBlock( pickle.dumps( self.w ) )
            except :
                # sometimes xmlrpclib shits itself
                pdb.set_trace()
                
        print "Shutting down..."
        
    def G_s( self, rh1, R ) :
        """
        Utility function for computing xi_s(r,R). Not for general
        use.
        
        We must take care here to remember that our definition for
        dp/dr is a spline, and therefor defined piecewise. However,
        we have computed the value of the G function analytically
        assuming that it is a cubic polynomial. To properly compute
        the G function for our splined definition of dp/dr, we must
        also evaluate the integral piecewise. 
        
        The intervals of integration are like so :
        
            R --> R[i]  
                  R[i] --> R[i+2] ... 
                                      R[i+n] --> rh1
        
        At the beginning of the serise, we have a quantity resulting
        from integration with the lower boundary of R and an upper 
        boundary corresponding to the boundary of the local spline
        interval.
        
        At the end of the series, we have a quanitity resulting 
        from integration from the beginning of the local spline 
        interval and ending with rh1.
        
        Between these two quantities, each element in the series is
        the result of integration bounded by the spline intervals
        themselves. The G function is the sum of these quantities.
        
        In the larger problem, R is ALWAYS larger tha Rhat. However,
        this function returns results for any R and Rhat between
        Rmin and Rmax. This allows us to define the xi function
        in a way that it will return values below the Z axis. These
        results are physically meaningless, but they are useful for
        improving convergance when locating the roots of Z.
        """
        # shortcuts!
        p_s = self.Problem.p_s
        
        Rp = p_s.R
        p  = p_s.Z
        dp = p_s.Z2
        
        K = self.G_s_K
        
        if R < self.Problem.Rmin or R > self.Problem.Rmax :
            return 0
        
        #if rh1 > self.axis :
        #    rh1 = self.rPsi_s.get_Z_for_R(rh1)
        
        if rh1 < Rp[0] or rh1 > Rp[-1] or R < Rp[0] or R > Rp[-1] :
            raise blob.SplineError,                 \
            "Called value out of range. Value: "    \
            + str(R)
        
        # check our cached value for i; if it's right, use it
        if rh1 >= Rp[self.G_s_last_i] and rh1 < Rp[self.G_s_last_i + 1] :
            
            i = self.G_s_last_i
        
        else :
            
            # if we have a cache miss, it is most likly that
            # we've steped to the next interval
            if self.G_s_last_i + 1 < len(Rp) :
                if rh1 >= Rp[self.G_s_last_i + 1] and \
                    rh1 < Rp[self.G_s_last_i + 2] :
                    i = self.G_s_last_i + 1
            
            # if not, then find the right place in the 
            # table of R values...
            for i in range(len(Rp) - 1) :
                if rh1 >= Rp[i] and rh1 < Rp[i + 1] :
                    break

            # ...and cache our result.
            self.G_s_last_i = i
        
        # check our cached value for j; if it's right, use it
        if R >= Rp[self.G_s_last_j] and R < Rp[self.G_s_last_j + 1] :
            
            j = self.G_s_last_j

        else :

            # if we have a cache miss, it is most likly that
            # we've steped to the next interval
            if self.G_s_last_j + 1 < len(Rp) :
                if rh1 >= Rp[self.G_s_last_j + 1] and \
                    rh1 < Rp[self.G_s_last_j + 2] :
                    j = self.G_s_last_j + 1

            # if not, find the right place in the
            # table of R values...
            for j in range(len(Rp) - 1) :
                if R >= Rp[j] and R < Rp[j + 1] :
                    break
        
            # ...and cache our result.
            self.G_s_last_j = j

        def g( a, b ) :
            """
            Analytic solution to the definite integral of p_s(x) * x,
            where p_s(x) is a spline interval of the pressure function,
            and the bounds of the integral are :
                
                a -> lower bound
                b -> upper bound
            
            NOTE :: a and b must be WITHIN the interval of the spline!
            That is, 
                Rp1 <= a <= Rp2
                
            and 
                
                Rp1 <= b <= Rp2
            
            This function is the result of integration by parts of the
            inner integral in the xi(r,R) function, and applying it to
            the definition of a spline found in Numerical Recipes.
            """
            h = Rp2 - Rp1
            
            AA = ( ( (Rp2 / 2) * b**2 - b**3 / 3) - \
                   ( (Rp2 / 2) * a**2 - a**3 / 3) ) / h
            
            BB = ( (b**3 / 3 - (Rp1 / 2) * b**2) - \
                   (a**3 / 3 - (Rp1 / 2) * a**2) ) / h
       
            CC = (h**2/6) * ( (1/h**3) * ( (1/4) * ( (Rp2 - b)**4*b      \
                                                  - ( Rp2 - a)**4*a )) - \
                 ( (1/20) * ( (Rp2 - b)**5 - (Rp2 - a)**5 ) ) - AA )
           
            DD = (h**2/6) * ( (1/h**3) * ( (1/4) * ( (Rp1 - b)**4*b      \
                                                  - ( Rp1 - a)**4*a )) - \
                 ( (1/20) * ( (b - Rp1)**5 - (a - Rp1)**5 ) ) - BB )
            
            return AA * p1 + BB * p2 + CC * dp1 + DD * dp2 
            
        # if the memoization table for the 
        # inner quantities doesn't exist...
        if K == [] :
           
            # ...then build it
            
            G = 0
            for k in range(len(Rp) - 1) :
                Rp1 = Rp[k]
                Rp2 = Rp[k + 1]
                p1  = p[k]
                p2  = p[k + 1]
                dp1 = dp[k]
                dp2 = dp[k + 1]
        
                G += g(Rp1, Rp2)
                K.append(G)
        
        # The residual from the analytic integration by parts.
        # It's invariant over spline intervals, so we put it
        # outside the g(a,b) function.
        LL = (R**2 - rh1**2) * p_s.get_Z_for_R(rh1)
        
        # If R and rh1 are in the same spline interval,
        # then we've got a very simple task.
        if i == j :
            
            Rp1 = Rp[i]
            Rp2 = Rp[i + 1]
            p1  = p[i]
            p2  = p[i + 1]
            dp1 = dp[i]
            dp2 = dp[i + 1]
 
            return LL + 2*g(R,rh1)

        # If R and rh1 span spline intervals, then we
        # have to compute the integral as a sum over
        # multiple intervals.
        if R < rh1 :
            
            G = 0.0            
            
            # contribution from the
            # interval containing R
            Rp1 = Rp[j]
            Rp2 = Rp[j + 1]
            p1  = p[j]
            p2  = p[j + 1]
            dp1 = dp[j]
            dp2 = dp[j + 1]
            
            G += g(R, Rp2)
            
            # use pre-computed inner quantities
            G += K[i-1] - K[j]
            
            # contribution from the
            # interval containing rh1
            Rp1 = Rp[i]
            Rp2 = Rp[i + 1]
            p1  = p[i]
            p2  = p[i + 1]
            dp1 = dp[i]
            dp2 = dp[i + 1]
            
            G += g(Rp1, rh1)
            
            return LL + 2*G
            
        else :
            
            G = 0.0
            
            # contribution from the 
            # interval containing R
            Rp1 = Rp[j]
            Rp2 = Rp[j + 1]
            p1  = p[j]
            p2  = p[j + 1]
            dp1 = dp[j]
            dp2 = dp[j + 1]
            
            G += g(R, Rp1)
            
            # use pre-computed inner quantities
            G += K[i] - K[j-1]
            
            # contribution from the
            # interval containing rh1
            Rp1 = Rp[i]
            Rp2 = Rp[i + 1]
            p1  = p[i]
            p2  = p[i + 1]
            dp1 = dp[i]
            dp2 = dp[i + 1]
     
            G += g(Rp2, rh1)
            
            return LL + 2*G
    
    def xi_s( self, r, R, return_spline='no' ) :
        """
        The xi function, splined. 
        """
        # shortcuts!
        Rmin  = self.Problem.Rmin
        Rmax  = self.Problem.Rmax
 
        def a ( x ) :
            """
            xi integrand.
            """
            
            if x <= Rmin or x >= Rmax :
                return 0
            
            if x == R :
                return 0
            
            G = self.G_s( x, R )
            
            if G > 0.0 :
                return 0
            
            b = ( -2.0 * mu_0 * G )**(0.5)
            
            if b == 0 :
                return 0
            
            c = self.Problem.dPsi_s.get_Z_for_R(x) / b
            
            if isnan(b) or isinf(b) :
                return 0
            
            return c
       
        if return_spline == 'no' :
        
            if r >= R or ( R - r ) < 0.0 :
                return 0        
            
            try :
                
                xi = integrate.quad(a, Rmin, r, full_output=1)[0]
            
            except integrate.quadpack.error, e :
            
                print "       ", e
                print "        Arguments ::> Rhat =", r, "R =", R
                return 0.0
        
            return xi
    
        if return_spline == 'yes' :
            
            if R < Rmin or R > Rmax :
                return 0
            
            # prune impossible values from Rhat list
            
            rr = []

            for i in r :
                if i <= R               \
                    and (R - i) > 6e-6  \
                    and i >= Rmin       \
                    and i <= Rmax :
                    rr.append(i)
            
            if len(rr) < 2 :
                return 0
            
            xixi = []

            for i in rr :
                
                try :
                    
                    xi = integrate.quad(a, Rmin, i, full_output=1)[0]
                
                except integrate.quadpack.error, e :
                    
                    print "       ", e
                    print "        Arguments ::> Rhat =", r, "R =", R
                    return 0.0
                
                xixi.append(xi)
            
            s = blob.spline()
            s.build( rr, xixi )
            
            return s

    def psi_solver( self, R1, Z1, RhGuess='none', xi_bucket='none' ) :
        """
        Return the value of psi for a given point { R1, Z1 }, given
        splined plasma boundary l_s.
        
        We require two l_s splines; l_s1 in the usual coordinates,
        and l_s2 in coordinates rotated -90 degrees. This is to
        allow us to avoid solving for differences of quantities
        tending toward infinity.
        """
        
        Rmin = self.Problem.Rmin
        Rmax = self.Problem.Rmax
        axis = self.Problem.axis
        l1_s = self.Problem.l1_s
        l2_s = self.Problem.l2_s
        
        def xi_b( rr, RR ) :
            NR = self.Problem.r_n(RR)
            if xi_bucket[NR] == 0 :
                return 0
            else :
                if rr > xi_bucket[NR].R[0] and rr < xi_bucket[NR].R[-1] :
                    return xi_bucket[NR].get_Z_for_R(rr)
                else :
                    return 0
        
        if xi_bucket == 'none' :
            xi_s = self.xi_s
        else :
            xi_s = xi_b
            
        def a( R ) :
            """
            Private function for calculating geometric
            values of Xi.
            """
        
            if R > Rmin and R < Rmax :
            
                try :
                    if Z1 >= 0.0 :
                        # use the top spline
                        result = l1_s.get_Z_for_R(R) - Z1 + \
                            (R - R1) / l1_s.get_dZ_for_R(R)
                    if Z1 < 0.0 :
                        # use the bottom spline
                        result = l2_s.get_Z_for_R(R) - Z1 + \
                            (R - R1) / l2_s.get_dZ_for_R(R)
                    
                except blob.SplineError, e :
                    print e
                    return 1.0
                
                else :
                    return result
            
            else :
                return 1.0
        
        def b( Rh ) :
            """
            Private function for calculating values of
            Rhat by comparing geometric values of Xi
            at (R1, Z1) to values of xi(Rhat,R).
            """
            try :   
                
                if Rh < Rmin :
                    return 1.0 
                if Rh == R1 :
                    return -1.0
                if Rh >= R1 :
                    return -1.0 - Rh + R1
                else :
                    return Xi - self.xi_s( Rh, R1 )
                
            except :
                
                pdb.set_trace()
                
        a0 = ( Rmax - Rmin ) / 2.0
        R0 = Rmin + a0
        
        if R1 >= R0 :
            RGuess = R0 + a0 / ( 1 + (Z1/(R0-R1))**2 )**(0.5)
        else :
            RGuess = R0 - a0 / ( 1 + (Z1/(R0-R1))**2 )**(0.5)
        
        # The radial case
        if Z1 == 0.0 :
            if R1 > R0 :
                R = Rmax
            else :
                R = Rmin
        
        else :
            R = optimize.fsolve(a, RGuess, full_output=1)[0]
        
        # stash the result for next time
        self.R_saved = R
        
        # Find the location of Z on the plasma boundary
        try :
        
            if Z1 > 0.0 :
                # use top spline
                Z = l1_s.get_Z_for_R(R)
            if Z1 < 0.0 :
                # use bottom spline
                Z = l2_s.get_Z_for_R(R)
            if Z1 == 0.0 :
                # radial case
                Z = 0.0
        
        except blob.SplineError, e :
        
            print "Xi was not calculated correctly for point R: " + str(R1) \
                + " Z: " + str(Z1)
            print e
            return 0.0
            
        Xi = ( ( R - R1 )**2 + ( Z - Z1 )**2 )**(0.5)
        
        if RhGuess == 'none' :
            
            RhGuess = Rmin + (axis + Rmin) / 2.0
        
        # find our local Rhat value
        Rh = optimize.fsolve(b, RhGuess)
        
        v = (b(Rh))**2
        
        if v > 0.1 :
            print "       ::> Failed once, trying again"
            Rh = optimize.fsolve(b, Rmin + 2e-4)
            v = (b(Rh))**2
            
        if v**2 > 0.1 :
            print "       :::> Failed twice, trying again"
            Rh = optimize.fsolve(b, axis - 2e-2)
            v = (b(Rh))**2
            
        if v**2 > 0.1 :
            print "       ::::> Failed thrice, trying again"
            Rh = optimize.fsolve(b, R1 - 2e-4)
            v = (b(Rh))**2
         
        if v**2 > 0.1 :
            print "       :::::> Failed four times, trying again"
            Rh = optimize.fsolve(b, (axis+R1)/2.0)
            v = (b(Rh))**2
         
        if v**2 > 0.1 :
            print "       ::::::> Bad point:", [R1,Z1]
            self.Rh_saved = -1
            return -1
        
        # save our Rhat value for next time
        self.Rh_saved = Rh
     
        return self.Problem.Psi_s.get_Z_for_R(Rh)    

