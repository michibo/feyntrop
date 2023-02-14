
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
# MM and mm stand for the muon and electron masses respectively
edges = [((0,1), 1, '0'), ((0,4), 1, 'MM'), ((1,5), 1, 'mm'), ((5,2), 1, 'mm'), 
         ((5,3), 1, '0'), ((4,3), 1, 'MM'), ((4,2), 1, '0')]

# replace scalar products in terms of chosen kinematic variables
replacement_rules = [(sp[0,0], 'MM'), (sp[1,1], 'mm'), (sp[2,2], 'mm'), (sp[0,1], '(s-MM-mm)/2'),
                     (sp[1,2], '(t-2*mm)/2'), (sp[0,2], '(MM+mm-s-t)/2')]

# numerically evaluate at this point
phase_space_point = [('MM', 1), ('mm', 1/200), ('s', -1/7), ('t', -1/3)]

# D = D0 - 2*eps dimensions
D0 = 6

# expand up to, but not including, eps_order
eps_order = 5

# deformation tuning
Lambda = 1.29

# number of sampling points
N = int(1e7)

# epsilon expansion without prefactor (trop_res) and normalization of tropical measure (Itr)
trop_res, Itr = tropical_integration(N, D0, Lambda, eps_order, edges, replacement_rules, phase_space_point)

# epsilon expansion with prefactor
expansion = eps_expansion(trop_res, edges, D0)
print("\n" + str(expansion) + "\n")
