 
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
graph = [((0,1), 1), ((1,2), 1), ((2,0), 1), ((0,5), 1), ((1,4), 1), ((2,3), 1), ((3,4), 1), ((4,5), 1), ((5,3), 1)]

# electron mass squared
me = 1

# squared momentum 
momentum_vars = [(p_sqr[0], 0)]

# internal masses
masses_sqr = [me, me, me, 0, 0, 0, me, me, me]

# D = D0 - 2*eps dimensions
D0 = 4

# expand up to, but not including, eps_order
eps_order = 9

# deformation tuning (no momenta, hence no analytic continuation)
Lambda = 0

# number of sampling points
N = int(1e8)

# epsilon expansion without prefactor (trop_res) and normalization of tropical measure (Itr)
trop_res, Itr = tropical_integration(graph, masses_sqr, momentum_vars, D0, eps_order, Lambda, N)

# epsilon expansion with prefactor
expansion = eps_expansion(trop_res, graph, D0)
print("\n" + str(expansion) + "\n")

