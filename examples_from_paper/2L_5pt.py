
    # ======
    # 2L_5pt
    # ======

# import feyntrop from parent directory
import os
import sys
path = os.path.abspath('..')
sys.path.append(path)
from py_feyntrop import *

# Feynman diagram
graph = [((0,1), 1), ((1,2), 1), ((2,6), 1), ((6,3), 1), ((3,4), 1), ((4,5), 1), ((5,0), 1), ((5,6), 1)]

# some quark mass squared
m = 1/2

# squared momenta and Mandelstam variables
momentum_vars = [(p_sqr[0], 0), (p_sqr[1], m), (p_sqr[2], m), (p_sqr[3], m), (s[0,1], 2.2), (s[0,2], 2.3), (s[0,3], 2.4), (s[1,2], 2.5), (s[1,3], 2.6), (s[2,3], 2.7)]

# internal masses
masses_sqr = [0, m, 0, m, 0, m, 0, m]

# D = D0 - 2*eps dimensions
D0 = 6

# expand up to, but not including, eps_order
eps_order = 5

# deformation tuning
Lambda = 0.28

# number of sampling points
N = int(1e8)

# epsilon expansion without prefactor (trop_res) and normalization of tropical measure (Itr)
trop_res, Itr = tropical_integration(graph, masses_sqr, momentum_vars, D0, eps_order, Lambda, N)

# epsilon expansion with prefactor
expansion = eps_expansion(trop_res, graph, D0)
print("\n" + str(expansion) + "\n")
