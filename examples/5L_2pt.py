 
    # ======
    # 5L_2pt
    # ======

# import feyntrop from parent directory
import os
import sys
path = os.path.abspath('..')
sys.path.append(path)
from py_feyntrop import *

# Feynman diagram
# current bug: mm10 would get replaced by mm1 in phase_space_point, so calling it MM10
edges = [
    ((0,6), 1, 'mm0'), ((0,5), 1, 'mm1'), ((5,6), 1, 'mm2'), 
    ((6,4), 1, 'mm3'), ((5,3), 1, 'mm4'), ((5,4), 1, 'mm5'), 
    ((4,3), 1, 'mm6'), ((4,2), 1, 'mm7'), ((3,2), 1, 'mm8'), 
    ((3,1), 1, 'mm9'), ((2,1), 1, 'MM10')]

# replace scalar product
replacement_rule = [(sp[0,0], 'pp')]

# numerically evaluate at this point
phase_space_point = [
        ('mm0', 1), ('mm1', 2), ('mm2', 3), 
        ('mm3', 4), ('mm4', 5), ('mm5', 6), 
        ('mm6', 7), ('mm7', 8), ('mm8', 9), 
        ('mm9', 10), ('MM10', 11), ('pp', 100)]

# D = D0 - 2*eps dimensions
D0 = 3

# expand up to, but not including, eps_order
eps_order = 11

# deformation tuning
Lambda = 0.02

# number of sampling points
N = int(1e7)

# epsilon expansion without prefactor (trop_res) and normalization of tropical measure (Itr)
trop_res, Itr = tropical_integration(N, D0, Lambda, eps_order, edges, replacement_rule, phase_space_point)

# expansion of gamma functions to 11th order is slow in python, please try e.g. mathematica instead
# eps_expansion(trop_res, edges, D0)
