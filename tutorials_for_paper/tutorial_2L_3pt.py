
            #######################
            # Example: 2L_3pt graph
            #######################

# feyntrop module
import os; import sys
path = os.path.abspath('..')
sys.path.append(path)
from pytrop import *

# Feynman graph
graph = [((0,1), 1), ((1,3), 1), ((3,2), 1), ((2,0), 1), ((0,3), 1)]

# Kinematics
masses = [0.2] * len(graph)
momentum_vars = [
    (sp[0,0], 0),
    (sp[1,1], 0),
    (s[0,1],  1)
]

# Settings
D0 = 2
eps_order = 5
Lambda = 8.9
N = int(1e8)

# Tropical integration
trop_res, Itr = tropical_integration(graph, masses, momentum_vars, D0, eps_order, Lambda, N)

# Epsilon expansion with prefactor
eps_expansion(trop_res, graph, D0)
