 
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
graph = [((0,6), 1), ((0,5), 1), ((5,6), 1), ((6,4), 1), ((5,3), 1), ((5,4), 1), ((4,3), 1), ((4,2), 1), ((3,2), 1), ((3,1), 1), ((2,1), 1)]

# squared momentum 
momentum_vars = [(p_sqr[0], 100)]

# internal masses
masses_sqr = [i for i in range(1,12)]

# D = D0 - 2*eps dimensions
D0 = 3

# expand up to, but not including, eps_order
eps_order = 11

# deformation tuning
Lambda = 0.02

# number of sampling points
N = int(1e8)

# epsilon expansion without prefactor (trop_res) and normalization of tropical measure (Itr)
trop_res, Itr = tropical_integration(graph, masses_sqr, momentum_vars, D0, eps_order, Lambda, N)

# expansion of gamma functions to 11th order is slow in python, please try e.g. mathematica instead
# eps_expansion(trop_res, graph, D0)
