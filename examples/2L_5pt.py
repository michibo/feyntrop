
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
# 'mm' stands for some quark mass squared
edges = [((0,1), 1, '0'), ((1,2), 1, 'mm'), ((2,6), 1, '0'), ((6,3), 1, 'mm'), 
         ((3,4), 1, '0'), ((4,5), 1, 'mm'), ((5,0), 1, '0'), ((5,6), 1, 'mm')]

# replace scalar products in terms of chosen kinematic variables s_ij = (p_i + p_j)^2
replacement_rules = [(sp[0,0], '0'), (sp[1,1], 'mm'), (sp[2,2], 'mm'), (sp[3,3], 'mm'),
                     (sp[0,1], '(s01-mm)/2'), (sp[0,2], '(s02-mm)/2'), (sp[0,3], '(s03-mm)/2'),
                     (sp[1,2], '(s12-2*mm)/2'), (sp[1,3], '(s13-2*mm)/2'), (sp[2,3], '(s23-2*mm)/2')]

# numerically evaluate at this point
phase_space_point = [('mm',  1/2), ('s01', 2.2), ('s02', 2.3), ('s02', 2.3), 
                     ('s03', 2.4), ('s12', 2.5), ('s13', 2.6), ('s23', 2.7)]

# D = D0 - 2*eps dimensions
D0 = 6

# expand up to, but not including, eps_order
eps_order = 5

# deformation tuning
Lambda = 0.28

# number of sampling points
N = int(1e7)

# epsilon expansion without prefactor (trop_res) and normalization of tropical measure (Itr)
trop_res, Itr = tropical_integration(N, D0, Lambda, eps_order, edges, replacement_rules, phase_space_point)

# epsilon expansion with prefactor
expansion = eps_expansion(trop_res, edges, D0)
print("\n" + str(expansion) + "\n")
