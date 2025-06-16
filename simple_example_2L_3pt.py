# Authors: Michael Borinsky & Henrik Munch 2023
#
# Email: michael.borinsky@eth-its.ethz.ch
#
# This example script accompanies the publication:
#
# M. Borinsky, H. J. Munch, F. Tellander: 'Tropical Feynman integration in
# the Minkowski regime', Computer Physics Communications, 292 (2023), 108874
# arXiv:2302.08955

# Running it in the terminal with
#
# > python simple_example_2L_3pt.py
#
# integrates a 2-loop 3-point Feynman graph.
#
# The comments in this file should completely explain the python interface of
# feyntrop.

# The first step is to import the `feyntrop` module. If there is an error
# here, make sure the files `feyntrop` and `py_feyntrop.py` are in the working
# directory.

from py_feyntrop import sp, tropical_integration, eps_expansion, prepare_kinematic_data

# A Feynman graph is given as a list of edges (propagators). Each edge or
# propagator is put in as a tuple in the format:
#
# ((v,u), nu, 'msqr')
#
# where
#
# * u,v are the vertices that the edge/propagator connects
# * nu is the edge weight
# * msqr is the squared mass of the edge/propagator.
#
# The vertices are numbered starting with 0.
#
# IMPORTANT:
#
# In the Python interface, external vertices have to have smaller numbers
# than internal ones. I.e. external vertices *have to be* labeled 0,1,2,...
# and internal ones *have to be labeled* labeled n,n+1,... where n is the
# number of external vertices.
#
#
# The edge weights must be numbers > 0.
#
#
# The squared masses can either be floating point numbers (in quotes) or
# strings that represent variables that are later replaced by numbers.

# Here, we integrate the 2-loop 3-point graph from arXiv:2302.08955 (see
# the picture on page 21 in the arXiv version). The edge/propagator list is:

edges = [((0, 1), 1, 'mm'), ((1, 3), 1, 'mm'), ((3, 2), 1, 'mm'),
         ((2, 0), 1, 'mm'), ((0, 3), 1, 'mm')]

# All edge weights are 1. The external vertices are 0, 1 and 2. All
# propagators have a uniform squared mass m^2 represented by the symbol 'mm'.

# The next step is to fix kinematics. We must fix the three external momenta
# p_0, p_1, p_1 in-coming into vertices 0, 1 and 2. This is done by fixing all
# the scalar products p_u*p_v for u,v in {0,1}, as the remaining momentum
# coming into vertex 2 is fixed by momentum conservation.
#
# feyntrop uses the convention that all momenta are in-coming and that
# Euclidean kinematics give rise to a negative semi-definite matrix sp[u,v].
# (see the script examples/1L_3pt_conformal.py for a Euclidean example)
#
# In the running example, we want to set p_0^2 = p_1^2 = 0, and the remaining
# free paramter fixes is the scalar product p_0 * p_1. Due to momentum
# conservation this is equivalent to fixing p_2^2, because
# p_2^2 = (-p_0 - p_1)^2 = 2 \cdot p_0 \cdot p_1.
#
# In code, we write

replacement_rules = [(sp[0, 0], '0'), (sp[1, 1], '0'), (sp[0, 1], 'pp2/2')]

# sp[u,v] stands for the scalar product p_u * p_v. Implicitly, we define a new
# variable here, 'pp2', which represents p_2^2, which is equal to 2 p_1*p_2.

# The two symbolic variables in `edges` and `replacement_rules`, namely 'mm'
# and 'pp2', must still get explicit numerical values.
#
# Here, we choose $m^2 = 0.2$ and $p_2^2 = 1$:

phase_space_point = [('mm', 0.2), ('pp2', 1)]

# We can print the resulting $V \times V$ dimensional scalar product matrix
# $\mathcal{P}^{u,v} = p_u \cdot p_v$ (where $V = |V_\text{ext}| +
# |V_\text{int}|$) and the list of squared edge masses as follows:

P_uv, m_sqr_list = prepare_kinematic_data(edges, replacement_rules,
                                          phase_space_point)

print("P_uv-matrix:", P_uv)
print("masses^2 list:", m_sqr_list)

# ============================================================================
# OPTIONAL
#
# At Euclidean kinematic points, the matrix `P_uv` is negative semi-definite.
# At Minkowski kinematic points, it is generally indefinite. If you want to
# integrate in the Euclidean regime, you can check the signature of `P_uv`
# e.g. using the commands
#
# import numpy as np
# eig_vals, eig_vecs = np.linalg.eig(np.array(P_uv, dtype=float))
# print("Eigenvalues of P_uv:", eig_vals)
#
# ============================================================================

# Additional parameters need to be fixed for the integration:
#
# * `D0` is the integer part of the spacetime dimension $D = D_0 - 2\epsilon$.
# * We want to expand up to, but not including, `eps_order`.
# * `Lambda` denotes the deformation parameter $\lambda$.
# * `N` is the number of Monte Carlo sampling points.

D0 = 2
eps_order = 5
Lambda = 7.6
N = int(1e7)  # 1e7 is Python notation for 10^7

# Now, we are ready for the tropical integration:

print("=" * 80)
print("Starting tropical integration")

trop_res, Itr = tropical_integration(N, D0, Lambda, eps_order, edges,
                                     replacement_rules, phase_space_point)

print("Tropical integration finished")
print("=" * 80)

# The output variable
#
# * `trop_res` contains the value of the Feynman integral *without* prefactor.
# * `Itr` is the normalization factor in the tropical measure.
#
# The command also prints statistics or outputs errors if something goes
# wrong. E.g., if the input graph has a divergence.

# The following command still translates the output into the epsilon-expansion
# *with* the prefactor $\Gamma(\omega)/\Gamma(\nu_1) \cdots \Gamma(\nu_E) =
# \Gamma(2\epsilon+3)$ included, where \omega is the superficial degree of
# divergence, and $E$ is the number of edges. The output in this convention is
# equal to the usual momentum space representation.

expansion = eps_expansion(trop_res, edges, D0)
print(expansion)
