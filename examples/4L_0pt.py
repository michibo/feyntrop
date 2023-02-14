 
    # ======
    # 4L_0pt
    # ======

# import feyntrop from parent directory
import os
import sys
path = os.path.abspath('..')
sys.path.append(path)
from py_feyntrop import *

# Feynman diagram
# 'mm' is some fermion mass
edges = [((0,1), 1, 'mm'), ((1,2), 1, 'mm'), ((2,0), 1, 'mm'), 
         ((0,5), 1, '0' ), ((1,4), 1, '0' ), ((2,3), 1, '0' ), 
         ((3,4), 1, 'mm'), ((4,5), 1, 'mm'), ((5,3), 1, 'mm')]

# no replacement rules for momenta, since this is a vacuum graph
replacement_rules = []

# numerically evaluate at this point
phase_space_point = [('mm', 1)]

# D = D0 - 2*eps dimensions
D0 = 4

# expand up to, but not including, eps_order
eps_order = 9

# deformation tuning (no momenta, hence no analytic continuation)
Lambda = 0

# number of sampling points
N = int(1e8)

# epsilon expansion without prefactor (trop_res) and normalization of tropical measure (Itr)
trop_res, Itr = tropical_integration(N, D0, Lambda, eps_order, edges, replacement_rules, phase_space_point)

# epsilon expansion with prefactor
expansion = eps_expansion(trop_res, edges, D0)
print("\n" + str(expansion) + "\n")

