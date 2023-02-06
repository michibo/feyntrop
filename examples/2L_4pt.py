
    # ======
    # 2L_4pt
    # ======

# import feyntrop from parent directory
import os
import sys
path = os.path.abspath('..')
sys.path.append(path)
from py_feyntrop import *

# Feynman diagram
graph = [((0,1), 1), ((0,4), 1), ((1,5), 1), ((5,2), 1), ((5,3), 1), ((4,3), 1), ((4,2), 1)]

# muon mass squared
M = 1

# electron mass squared
m = M/200

# Mandelstam variables similar to section 4.1.2 of [1806.08241]
S = -1/7
T = -1/3
U = 2*M + 2*m - S - T 

# squared momenta and Mandelstam variables
momentum_vars = [(p_sqr[0], M), (p_sqr[1], m), (p_sqr[2], m), (s[0,1], S), (s[1,2], T), (s[0,2], U)]

# internal masses
masses_sqr = [0, M, m, m, 0, M, 0]

# D = D0 - 2*eps dimensions
D0 = 6

# expand up to, but not including, eps_order
eps_order = 5

# deformation tuning
Lambda = 1.29

# number of sampling points
N = int(1e8)

# epsilon expansion without prefactor (trop_res) and normalization of tropical measure (Itr)
trop_res, Itr = tropical_integration(graph, masses_sqr, momentum_vars, D0, eps_order, Lambda, N)

# epsilon expansion with prefactor
expansion = eps_expansion(trop_res, graph, D0)
print("\n" + str(expansion) + "\n")
