#!/usr/bin/env python
# vim: set ts=4 sw=4 et:
"""
Default values for the problem and physical constants.
"""

machine  = { 'Rmin': 4.0, 'Rmax': 6.0 }
boundary = { 'al': 1.0, 'bl': 0.2, 'cl': 0.0 }
profile  = 'D'  # D or reverse-D
pressure = { 'p0': 1000000.0, 'p2': -1.0, 'p3': 0.0, 'p4': 0.0 }
flux     = { 'af': 0.01, 'bf': 1.0, 'cf': 1.0 }
mu_0     = 1.25663706e-6
