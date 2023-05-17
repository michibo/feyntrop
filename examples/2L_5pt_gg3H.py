

    # ===========
    # 2L_5pt_gg3H
    # ===========

# import feyntrop from parent directory
import os
import sys
path = os.path.abspath('..')
sys.path.append(path)
from py_feyntrop import *

# Feynman diagram
# 'mm_top' stands for the top quark mass squared
edges = [
    ((0,1), 1, 'mm_top'), ((1,6), 1, 'mm_top'), ((5,6), 1, '0'),  ((6,2), 1, 'mm_top'),
    ((2,3), 1, 'mm_top'), ((3,4), 1, 'mm_top'), ((4,5), 1, 'mm_top'), ((5,0), 1, 'mm_top')
]

# replace scalar products in terms of Mandelstam variables s_ij = (p_i + p_j)^2 and the Higgs mass squared mm_higgs
replacement_rules = [
   (sp[0,0], '0'),
   (sp[1,1], '0'),
   (sp[2,2], 'mm_higgs'),
   (sp[3,3], 'mm_higgs'),
   (sp[0,1], '(5*mm_higgs - s02 - s03 - s12 - s13 - s23)/2'),
   (sp[0,2], '(-mm_higgs + s02)/2'),
   (sp[0,3], '(-mm_higgs + s03)/2'),
   (sp[1,2], '(-mm_higgs + s12)/2'),
   (sp[1,3], '(-mm_higgs + s13)/2'),
   (sp[2,3], '(-2*mm_higgs + s23)/2')
]

# numerically evaluate at this point
phase_space_point = [('mm_top', 1.8995), ('mm_higgs', 1),
                     ('s02', -4.4), ('s03', -0.5), ('s12', -0.6), ('s13', -0.7), ('s23', 1.8)
]

# D = D0 - 2*eps dimensions
D0 = 4

# expand up to, but not including, eps_order
eps_order = 5

# deformation tuning
Lambda = 0.64

# number of sampling points
N = int(1e7)

# epsilon expansion without prefactor (trop_res) and normalization of tropical measure (Itr)
trop_res, Itr = tropical_integration(N, D0, Lambda, eps_order, edges, replacement_rules, phase_space_point)

# epsilon expansion with prefactor
expansion = eps_expansion(trop_res, edges, D0)
print("\n" + str(expansion) + "\n")
