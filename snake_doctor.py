#!/usr/bin/env python
# vim: set ts=4 sw=4 et:

from ophidian import simplesnake, blob
from scipy import *
import pickle
from pylab import *

#f = open( 'snake_poop_noshift_trimmed.dump', 'r' )
#P = pickle.load(f)
#f.close()

P = P1

def cmp_neg( r ) :
    if r > p.Rmin and r < p.Rmax :
        return - ( p.Psi_s.get_Z_for_R( r ) - p.Psi_ss.get_Z_for_R( r ) )
    else :
        return 0.0

def cmp_pos( r ) :
    if r > p.Rmin and r < p.Rmax :
        return psi_bar * ( p.Psi_s.get_Z_for_R( r ) - p.Psi_ss.get_Z_for_R( r ) )
    else :
        return 0.0

def beta_neg( r ) :
    return - p.beta_s.get_Z_for_R( r )

def psi_neg( r ) :
    return - p.Psi_s.get_Z_for_R( r )

R_max_errors_hfs  = []
max_errors_hfs    = []
R_max_errors_lfs  = []
max_errors_lfs    = []
int_errors        = []
R_beta_max        = []
beta_max          = []
shift             = []
qbar              = []
qavg              = []
qmax              = []
theory_error_qbar = []
theory_error_qavg = []
theory_error_qmax = []

for p in P.list :
    
    a = ( p.Rmax - p.Rmin ) / 2.0
    
    psi_bar = ( p.Rmax - p.Rmin ) / integrate.quad( p.Psi_s.get_Z_for_R, p.Rmin, p.Rmax )[0]
    
    print p.problemname 
    R_max_errors_hfs.append(  optimize.fmin( cmp_neg, [p.Rmin + 0.001], full_output='yes' )[0] )
    max_errors_hfs.append(    psi_bar * cmp_pos( R_max_errors_hfs[-1] ) )
    R_max_errors_lfs.append(  optimize.fmin( cmp_neg, [p.Rmax - 0.001], full_output='yes' )[0] )
    max_errors_lfs.append(    psi_bar * cmp_pos( R_max_errors_lfs[-1] ) )
    int_errors.append(        integrate.quad( cmp_pos, p.Rmin, p.Rmax )[0] / \
                                integrate.quad( p.Psi_s.get_Z_for_R, p.Rmin, p.Rmax )[0] ) 
    R_beta_max.append(        optimize.fminbound( beta_neg, p.Rmin, p.Rmax ) )
    beta_max.append( p.beta_s.get_Z_for_R( R_beta_max[-1] ) )
    shift.append( ( a - ( p.Rmax - optimize.fminbound( psi_neg, p.Rmin, p.Rmax ) ) ) / a )
    
    # The predicted theoretical error
    qavg.append( p.q_s.Z.sum() / float( len(p.q_s.Z) ) )
    qbar.append( min(p.F_s.Z) / p.Psi_s.Z.max() )
    qmax.append( p.q_s.Z.max() )
    # theory_error.append( ( a**2 / R_max_errors_lfs[-1]**2 ) / \
    #   ( (p.Psi_s.Z.max() / p.F_s.Z.max())**2 * ( beta_max[-1] + 1 ) ) )
    theory_error_qbar.append( ( a**2 / R_max_errors_lfs[-1]**2 ) / ( qbar[-1]**2 * ( beta_max[-1] + 1 ) ) )
    theory_error_qavg.append( ( a**2 / R_max_errors_lfs[-1]**2 ) / ( qavg[-1]**2 * ( beta_max[-1] + 1 ) ) )
    theory_error_qmax.append( ( a**2 / R_max_errors_lfs[-1]**2 ) / ( qmax[-1]**2 * ( beta_max[-1] + 1 ) ) )

    plot( p.Psi_s.R, psi_bar * array( map( cmp_neg, p.Psi_s.R ) ) )
    r = optimize.fmin( cmp_neg, [p.Rmin+0.001], full_output='yes' )[0]
    plot( [r], [ psi_bar * cmp_neg( r ) ], 'ro' )
