
    # ======
    # 3L_4pt
    # ======

# import feyntrop from parent directory
import os
import sys
path = os.path.abspath('..')
sys.path.append(path)
from py_feyntrop import *

# Feynman diagram
graph = [((0,1), 1), ((1,2), 1), ((2,3), 1), ((3,0), 1), ((0,2), 2), ((1,3), 2)]

# squared momenta and Mandelstam variables
momentum_vars = [(p_sqr[0], 1.1), (p_sqr[1], 1.2), (p_sqr[2], 1.3), (s[0,1] , 2.1), (s[0,2] , 2.2), (s[1,2] , 2.3)]

# internal masses
masses_sqr = [0.05, 0.06, 0.07, 0.08, 0.09, 0.10]

# D = D0 - 2*eps dimensions
D0 = 4

# expand up to, but not including, eps_order
eps_order = 7

# deformation tuning
Lambda = 1.24

# number of sampling points
N = int(1e8)

# epsilon expansion without prefactor (trop_res) and normalization of tropical measure (Itr)
trop_res, Itr = tropical_integration(graph, masses_sqr, momentum_vars, D0, eps_order, Lambda, N)

# epsilon expansion with prefactor
expansion = eps_expansion(trop_res, graph, D0)
print("\n" + str(expansion) + "\n")
