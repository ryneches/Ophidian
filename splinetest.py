#!/usr/bin/env python
# vim: set ts=4 sw=4 et:
"""
Test the spline class.
"""

from ACUBE import *
from Numeric import *
from pylab import *

s = blob.spline()

def qdr(x) :
    #return x**2 + 2*x - 1
    return sin(x)

p_lowres = arange(0, 9.5, 0.5)
q_lowres = qdr(p_lowres)

p = arange(0, 9, 0.01)
q = qdr(p) 

s.build(p_lowres,q_lowres)

y_s = []

for i in p :
    y_s.append(s.get_Z_for_R(i))
    
subplot(211)

plot(s.R, s.Z, 'ro')
plot(p, y_s)

subplot(212)

d_s = []
for i in range(len(p)) :
    d_s.append( q[i] - y_s[i] )

plot(p, d_s)

show()
