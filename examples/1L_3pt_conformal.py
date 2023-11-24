
    # ================
    # 1L_3pt_conformal
    # ================

# import feyntrop from parent directory
import os
import sys
path = os.path.abspath('..')
sys.path.append(path)
from py_feyntrop import *

# Feynman diagram with half edge weights and massless internal lines
edges = [((0,1), 1/2, '0'), ((1,2), 1/2, '0'), ((2,0), 1/2, '0')]

# replace scalar products in terms of chosen kinematic variables. Here s = (p0 + p1)^2
replacement_rules = [(sp[0,0], 'pp0'), (sp[1,1], 'pp1'), (sp[0,1], '(s-pp0-pp1)/2')]

# numerically evaluate at this Euclidean point
phase_space_point = [('pp0', -2), ('pp1', -3), ('s', -5)]

# D = D0 - 2*eps dimensions
D0 = 2

# only leading term in epsilon expansion
eps_order = 1

# deformation tuning (Euclidean regime, hence no analytic continuation)
Lambda = 0

# number of sampling points
N = int(1e8)

# epsilon expansion without prefactor (trop_res) and normalization of tropical measure (Itr)
trop_res, Itr = tropical_integration(N, D0, Lambda, eps_order, edges, replacement_rules, phase_space_point)
