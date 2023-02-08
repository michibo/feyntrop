
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

# internal masses (indices of m_sqr correspond to graph)
masses_sqr = [(m_sqr[0], 0.05), (m_sqr[1], 0.06), (m_sqr[2], 0.07), (m_sqr[3], 0.08), (m_sqr[4], 0.09), (m_sqr[5], 0.10)]

# D = D0 - 2*eps dimensions
D0 = 4

# expand up to, but not including, eps_order
eps_order = 7

# deformation tuning
Lambda = 1.24

# number of sampling points
N = int(1e8)

# epsilon expansion without prefactor (trop_res) and normalization of tropical measure (Itr)
trop_res, Itr = tropical_integration(N, D0, Lambda, eps_order, graph, momentum_vars, masses_sqr)

# epsilon expansion with prefactor
expansion = eps_expansion(trop_res, graph, D0)
print("\n" + str(expansion) + "\n")
