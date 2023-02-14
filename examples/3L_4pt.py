
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
edges = [((0,1), 1, 'mm0'), ((1,2), 1, 'mm1'), ((2,3), 1, 'mm2'), 
         ((3,0), 1, 'mm3'), ((0,2), 2, 'mm4'), ((1,3), 2, 'mm5')]

# replace scalar products in terms of chosen kinematic variables s_ij = (p_i+p_j)^2
replacement_rules = [(sp[0,0], 'pp0'), (sp[1,1], 'pp1'), (sp[2,2], 'pp2'), 
                    (sp[0,1], '(s01-pp0-pp1)/2'), (sp[0,2], '(s02-pp0-pp2)/2'), (sp[1,2], '(s12-pp1-pp2)/2')]

# numerically evaluate at this point
phase_space_point = [('pp0', 1.1), ('pp1', 1.2), ('pp2', 1.3), 
                     ('s01', 2.1), ('s02', 2.2), ('s12', 2.3), 
                     ('mm0', 0.05), ('mm1', 0.06), ('mm2', 0.07), ('mm3', 0.08), ('mm4', 0.09), ('mm5', 0.10)]

# D = D0 - 2*eps dimensions
D0 = 4

# expand up to, but not including, eps_order
eps_order = 7

# deformation tuning
Lambda = 1.24

# number of sampling points
N = int(1e8)

# epsilon expansion without prefactor (trop_res) and normalization of tropical measure (Itr)
trop_res, Itr = tropical_integration(N, D0, Lambda, eps_order, edges, replacement_rules, phase_space_point)

# epsilon expansion with prefactor
expansion = eps_expansion(trop_res, edges, D0)
print("\n" + str(expansion) + "\n")
