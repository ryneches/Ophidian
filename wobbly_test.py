#!/usr/bin/env python
# vim: set ts=4 sw=4 et:

from scipy import *
from pylab import *
from ACUBE import *
from ACUBE2 import *
import psyco


def r_coor( n ) :
    return Rin + (Rout - Rin) * ( float(n) / float(NX) )

def z_coor( n ) :
    return Zin + (Zout - Zin) * ( float(n) / float(NY) )
    
def matrix_min( array ) :
    q = (0,0)
    for i in range(len(array)) :
        for j in range(len(array[0])) :
            if array[q[0]][q[1]] > array[i][j] :
                q = (i,j)
    return q

def matrix_max( array ) :
    q = (0,0)
    for i in range(len(array)) :
        for j in range(len(array[0])) :
            if array[q[0]][q[1]] < array[i][j] :
                q = (i,j)
    return q

a = io.array_import.read_array( '../CUBE_data/cube_results_5/rerun/bin/flux.txt' )
b = array(a) 

T = data.CUBEdata()
T.read( '../CUBE_data/cube_results_5/rerun/bin/results.txt' )

cube_axis = T.plasma_axis()
Rmin, Rmax = T.plasma_range()
T.slice( T.R_raw[0], cube_axis )

dp_s    = blob.spline()
dPsi_s  = blob.spline()
p_s     = blob.spline()
Psi_s   = blob.spline()
rPsi_s  = blob.spline()

dp_s.build( T.R_raw, T.dp_raw )
dPsi_s.build( T.R_raw, T.dpsi_raw )

p_s.build( T.R_raw, T.p_raw )
Psi_s.build( T.R_raw, T.psi_raw )

NX = len(a)
NY = len(a[0])
Rin = 3.9
Rout = 6.1
Zin = 1.1
Zout = -1.1

t = exp_pressure.tokamak()

# initialize the R-Rhat branch cut spline
aa = a[len(a)/2]
RRR = []
for j in range(len(aa)) :
    RRR.append(r_coor(j))

t.Rbranch_build( aa, RRR )

#plot(t.rPsi_s.R, t.rPsi_s.Z)
#show()


####### Quick and nasty test for the G_s function
#RRR = arange(1.0, 11.0 , 0.05)
#ZZZ = array(RRR)**(-1)
#zzz = -array(RRR)**(-2)
#ZZZ = array(zeros(200) + 1)
#ZZZ = sin(array(RRR))-sin(1)
#zzz = []
#for i in RRR :
#    zzz.append(cos(i)/i)

#foo = blob.spline()
#bar = blob.spline()
#foo.build(RRR,zzz)
#bar.build(RRR,zzz)

#t.Rmin = 0.0
#t.Rmax = 10.0
#t.p_s = foo
#t.dp_s = bar

#ppp = []
#ggg = []
#rrr = []

#for i in arange(1,10,0.01) :
#    ppp.append(t.p_s.get_Z_for_R(i))
#    ggg.append(t.G_s(1.0,i))
    #GGG.append(t.G_s(1.0,i, slow='yes'))
#    rrr.append(i)
    


#plot(RRR,ZZZ, color='black')
#plot(rrr,ppp, color='blue')
#plot(rrr,-array(ggg), color='red' )
#plot(rrr,GGG, color='green')
#show()

########## back to work...

t.dp_s      = dp_s
t.dPsi_s    = dPsi_s
t.p_s       = p_s
t.Psi_s     = Psi_s

t.p0 = -1.09829042e+01
t.p1 = 4.39420946e-03
t.p2 = 2.98816778e+00
t.p3 = 3.83557094e+00
t.a_psi = -0.02736385
t.b_psi = 6.66980231
t.c_psi = 0.27757685
t.d_psi = 6.64166538
t.BtRmin = 0.268526591 
#t.axis = 5.81517857
#t.axis_label = 5.81517857

t.CUBEdata = a
t.NX = len(a)
t.NY = len(a[0])
t.Rin = 3.9
t.Rout = 6.1
t.Zin = -1.1
t.Zout = 1.1

s = t.find_lcfs( a, 3.9, 6.1, -1.1, 1.1 )

for i in range(len(b)) :
    for j in range(len(b[0])) :
        b[i][j] = 0.0
    if a[i][j] < 0.0 :
        a[i][j] = 0.0

#s2 = t.find_lcfs( rot90(a, k=-1), -1.1, 1.1, 3.9, 6.1 )

t.Rmin = s[0].R[0]
t.Rmax = s[0].R[-1]

l = matrix_max(a)
psi_max = a[l[0]][l[1]]

t.axis = r_coor(l[1])

# assorted diagnostics
if 1==0 :

    # verify that the G_s function works
    
    RRR = arange(3.9, 6.1, 0.01)
    for Rh in arange(4.1, t.axis, 0.2) :
        GGG = []
        ggg = []
        #gGg = []
        print ":: plotting ", Rh
        for R in RRR :
            GGG.append(t.G_s(Rh, R))
            #gGg.append(t.G_s(Rh, R, slow='yes'))
            ggg.append(t.G(Rh, R))
        plot(RRR,GGG,color='red')
        plot(RRR,ggg,color='blue')
        #plot(RRR,gGg,color='green')
    show()

    # plot the xi(r,R) function
    RRR = arange(4.0, 6.0, 0.02)
    for Rh in arange(4.1, t.axis, 0.5) :
        xixi = []
        print "plotting xi for Rh =", Rh
        for i in RRR :
            xixi.append(t.xi_s(Rh,i))
        plot(RRR,xixi)
    show()

psyco.full()

#for i in range(0, len(a)) :
for i in range(373, 374) :
    print "row :", i
    for j in range(289, 336) :
        print [i,j]
        if a[i][j] >= 2e-6 :
            b[i][j] = t.psi_solver( s, r_coor(j), z_coor(i) )

plot(b[373])
show()

#c = abs(a-b)

#for i in range(120, len(a)-120) :
#    for j in range(len(a[0])) :
#        if 100*c[i][j]/psi_max > 20.0 :
#            print [i,j]
#            for q in arange(0.01, 2.0, 0.1) :
#                print "   ."
#                A = t.psi_solver( s, r_coor(j), z_coor(i), RhGuess=t.Rmax-q )
#                if 100*(abs(A - a[i][j]))/psi_max < 20.0 :
#                    b[i][j] = A
#                    print "   $"
#                    break
