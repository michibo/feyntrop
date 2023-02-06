
    # ================
    # 1L_3pt_conformal
    # ================

# import feyntrop from parent directory
import os; import sys
path = os.path.abspath('..')
sys.path.append(path)
from py_feyntrop import *

# Feynman diagram with half edge weights
graph = [((0,1), 1/2), ((1,2), 1/2), ((2,0), 1/2)]

# squared momenta
momentum_vars = [(p_sqr[0], -2), (p_sqr[1], -3), (s[0,1] , -5)]

# internal masses
masses_sqr = [0] * len(graph)

# D = D0 - 2*eps dimensions
D0 = 2

# only leading term in epsilon expansion
eps_order = 1

# deformation tuning (Euclidean regime, hence no analytic continuation)
Lambda = 0

# number of sampling points
N = int(1e8)

# epsilon expansion without prefactor (trop_res) and normalization of tropical measure (Itr)
trop_res, Itr = tropical_integration(graph, masses_sqr, momentum_vars, D0, eps_order, Lambda, N)
