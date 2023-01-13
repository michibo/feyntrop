
            #######################
            # Example: 1L_4pt graph
            #######################

    # ---------------------------------
    # Import package from parent folder
    # ---------------------------------

import os
import sys
path = os.path.abspath('..')
sys.path.append(path)
from pytrop import *

    # ---------
    # The graph
    # ---------

# Defining the graph in terms of vertices v = [0,1,2,3].
# Notation: ((v1,v2), w), where an edge goes between vertices v1 and v2, and w is the edge weight.
graph = [((0, 1), 1), ((1, 2), 1), ((2, 3), 1), ((3, 0), 1)]

    # ------------------
    # Momenta and masses
    # ------------------

# The momenta live in D = D0 - 2*eps dimensions.
D0 = 4

# Momentum variables.
# 'sp[i,i]' stands for scalar product, so (sp[0,0],1.1) means p0^2 = 1.1.
# 's[i,j]' stands for a Mandelstam variable, so (sp[0,1], 1.2) means (p0+p1)^2 = 1.2.
momentum_vars = [
    (sp[0,0], 1.1),
    (sp[1,1], 1.2),
    (sp[2,2], 1.3),
    (s[0,1] , 1.4),
    (s[0,2] , 1.5),
    (s[1,2] , 1.6)
]

# Mass for each internal edge in 'graph'.
masses = [0.1, 0.2, 0.3, 0.4]

    # ---------------------
    # Numerical integration
    # ---------------------

# Number of Monte Carlo sampling points.
N = int(1e7)

# Expand up to O(eps^epsorder).
eps_order = 3

# Contour deformation tuning parameter.
# The size of the errors is sensitive to the value of Lambda.
Lambda = 1.8

# Tropical integration.
# trop_res is the numerical value of the Feynman integral (without any prefactor).
# Itr is the tropical normalization factor.
trop_res, Itr = tropical_integration(graph, masses, momentum_vars, D0, eps_order, Lambda, N)

# Epsilon expansion with prefactor.
# Prefactor = (-1)^sum(a_e) * Gamma(omega)/(Gamma(a1)...Gamma(aE)).
eps_expansion(trop_res, graph, D0)
