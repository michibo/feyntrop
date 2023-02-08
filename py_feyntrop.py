import feyntrop
from math import log10
from sympy import prod, gamma, series, symbols, zeros, IndexedBase
from sympy.tensor import get_indices

    #=====
    # Misc 
    #=====

# m_sqr[i] = m_i^2, p_sqr[i] = p_i^2 and s[i,j] = (p_i+p_j)^2 are used in replacement rules
# eps is used as an expansion parameter 
m_sqr = IndexedBase('m_sqr')
p_sqr = IndexedBase('p_sqr')
s = IndexedBase('s')
eps = symbols('eps')

# un-nest a list
def flatten(list):
    return [item for sublist in list for item in sublist]

# vertices of a graph
def vertices(graph):
    graph = [e[0] for e in graph] # remove edge weights --> only edges remain
    graph = flatten(graph)
    return list(set(graph)) # unique integers in graph

# find largest index in momentum_vars to get V_ext = n_particles
# last line: + 2 = 1 (zero indexing) + 1 (only provide data for V_ext-1 particles in momentum_vars)
def get_V_ext(momentum_vars):
    indices = [
        momentum_vars[i][0] . indices for i in range(len(momentum_vars))
    ]
    indices = flatten(indices)
    return max(indices) + 2

# those masses not specified in masses_sqr are set to zero
def masses_sqr_zero(masses_sqr, graph):
    E = len(graph)
    n = len(masses_sqr)
    m_indices = []
    for i in range(n):
        msqr = masses_sqr[i][0]
        m_indices . append(list( get_indices(msqr)[0] ))
    m_indices = flatten(m_indices)
    g_indices = list(range(E))
    complement = list(set(g_indices) - set(m_indices))
    return [(m_sqr[i], 0) for i in complement]

    #==================
    # Formatting result
    #==================

# round result to two significant digits
def round_res(res_err_list):
    res = res_err_list[0]
    err = res_err_list[1]
    if err == 0:
        return [str(res), str(err)]
    nzeros = abs(int( log10(err) )) # number of zeros after decimal pt.
    nround = nzeros + 2
    res = round(res, nround)
    err = round(err, nround)
    err = f'{err:.{nround}f}'
    res = f'{res:.{nround}f}'
    return [res, err]

# print each epsilon order with rounded result
def print_res(res):
    round_re    , round_im    = [] , []
    lens_re     , lens_im     = [] , []
    lens_re_err , lens_im_err = [] , []
    eps_order = len(res)
    for i in range(eps_order):
        #
        round_re . append(round_res( res[i][0] ))
        round_im . append(round_res( res[i][1] ))
        #
        lens_re     . append(len( round_re[i][0] ))
        lens_re_err . append(len( round_re[i][1] ))
        #
        lens_im     . append(len( round_im[i][0] ))
        lens_im_err . append(len( round_im[i][1] ))
        #
    max_re = str(max( lens_re ))
    max_im = str(max( lens_im ))
    max_re_err = str(max( lens_re_err ))
    max_im_err = str(max( lens_im_err ))
    for i in range(eps_order):
        re , re_err = round_re[i]
        im , im_err = round_im[i]
        eps_str = '-- eps^' + str(i) + ': '
        #
        print(f"{eps_str : <5}{'[' : ^1}{re : ^{max_re}}{' +/- ' : ^5}{re_err : ^{max_re_err}}{']' : ^1}{'  +  i * ' : ^9}{'[' : ^1}{im : ^{max_im}}{' +/- ' : ^5}{im_err : ^{max_im_err}}{']' : ^1}")

# prefactor in Symanzik representation
def prefactor(graph, D0, eps_order):
    V = len(vertices(graph))
    E = len(graph)
    nu = [e[1] for e in graph] # propagator powers
    nu_tot = sum(nu)
    L = E - V + 1 # loops
    w = nu_tot - L * (D0-2*eps)/2 # superficial degree of divergence
    num = gamma(w)
    den = prod([gamma(nu[i]) for i in range(E)])
    return num/den

# mulitplying result by prefactor and expanding in epsilon
def eps_expansion(res, graph, D0):
    eps_order = len(res)
    terms = 0
    for i in range(eps_order):
        re = res[i][0][0]
        im = res[i][1][0]
        terms += complex(re, im)*eps**i
    pf = prefactor(graph, D0, eps_order)
    terms = pf * terms
    terms = series(terms, eps, 0, eps_order)
    return terms . evalf(10)

    #=========================================
    # Inverting Mandelstams to scalar products
    #=========================================

def pi_pj(graph, momentum_vars):
    #
    # 1-pt function
    #
    if vertices(graph) == [0]:
        return [0]
    #
    # 2-pt or higher
    #
    V = len(vertices(graph)) # total number of vertices (external + internal)
    V_ext = get_V_ext(momentum_vars) # number of particles
    pi_pj_mat = zeros(V) # V x V matrix with pi*pj = 0 if i or j internal vertex
    #
    # first 1, ..., (V_ext-1) diagonals
    #
    for i in range(V_ext-1):
        pi_pj_mat[i,i] = p_sqr[i] 
    #
    # upper triangular (V_ext-1) x (V_ext-1) block
    #
    for i in range(V_ext-1):
        for j in range(i+1, V_ext-1):
            pi_pj_mat[i,j] = (s[i,j] - p_sqr[i] - p_sqr[j]) / 2
    #
    # copy to lower triangular (V_ext-1) x (V_ext-1) block
    #
    for i in range(V_ext-1):
        for j in range(i+1, V_ext-1):
            pi_pj_mat[j,i] = pi_pj_mat[i,j]
    #
    # last entry of column = (-1) * everything above 
    #
    for j in range(V_ext-1):
         colsum = 0
         for i in range(V_ext-1):
             colsum += pi_pj_mat[i,j]
         pi_pj_mat[V_ext-1,j] = (-1) * colsum
    #
    # last entry of row = (-1) * everything to the left
    #
    for i in range(V_ext):
         rowsum = 0
         for j in range(V_ext-1):
             rowsum += pi_pj_mat[i,j]
         pi_pj_mat[i,V_ext-1] = (-1) * rowsum
    #
    # subsitute phase space point
    #
    return pi_pj_mat . subs(momentum_vars) . tolist()

    #=====================
    # Tropical integration
    #=====================

# trop_int is the value of the Feynman integral (without prefactor!)
# Itr is the normalization factor in the tropical measure
# if no list of masses is given, all are assumed zero
# if a list of masses *is* given, but some masses are missing, these are assumed zero
def tropical_integration(N, D0, Lambda, eps_order, graph, momentum_vars, masses_sqr = None):
    #
    E = len(graph)
    if masses_sqr == None:
        masses_sqr = [(m_sqr[e], 0) for e in range(E)] # m_e^2 = 0 if no list is given
    masses_sqr_rest = masses_sqr_zero(masses_sqr, graph) # insert zeros for missing masses if != None
    masses_sqr = masses_sqr + masses_sqr_rest
    masses_sqr = [m_sqr[e] . subs(masses_sqr) for e in range(E)]
    #
    g = feyntrop.graph(graph)
    pipj = pi_pj(graph, momentum_vars)
    print("Prefactor: " + str(prefactor(graph, D0, eps_order)) + ".")
    trop_res, Itr  = feyntrop.integrate_graph(g, D0, pipj, masses_sqr, eps_order, Lambda, N)
    print("")
    print_res(trop_res)
    return trop_res, Itr
