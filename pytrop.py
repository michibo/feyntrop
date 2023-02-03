import feyntrop
from math import log10
from sympy import prod, gamma, series, symbols, zeros, IndexedBase

    ######
    # Misc 
    ######

# sp[i,i] = scalar product = p_i^2 and s[i,j] = Mandelstam = (p_i+p_j)^2 are used in replacement rules
# eps is used as an expansion parameter 
sp  = IndexedBase('sp')
s   = IndexedBase('s')
eps = symbols('eps')

# un-nest a list
def flatten(list):
    return [item for sublist in list for item in sublist]

# vertices of a graph
def vertices(graph):
    graph = [e[0] for e in graph] # remove edge weights --> only edges remain
    graph = flatten(graph)
    return list(set(graph)) # unique integers in graph

# find largest index in pspt to get Vext = n_particles
# last line: + 2 = 1 (zero indexing) + 1 (only provide data for Vext-1 particles in pspt)
def get_Vext(pspt):
    indices = [
        pspt[i][0] . indices for i in range(len(pspt))
    ]
    indices = flatten(indices)
    return max(indices) + 2

    ###################
    # Formatting result
    ###################

# round result to two significant digits
def round_res(res, err):
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
    eps_order = len(res)
    for i in range(eps_order):
        #
        re , re_err = res[i][0]
        im , im_err = res[i][1]
        re , re_err = round_res(re , re_err)
        im , im_err = round_res(im , im_err)
        #
        eps_str = '-- eps^' + str(i) + ': '
        pm_str = " +/- "
        #
        re_str = "[" + re + pm_str + re_err + "] + "
        im_str = "[" + im + pm_str + im_err + "] * i"
        #
        print(eps_str + re_str + im_str)

# prefactor in Symanzik representation
def prefactor(graph, D0, eps_order):
    V = len(vertices(graph))
    E = len(graph)
    nu = [e[1] for e in graph] # propagator powers
    nu_tot = sum(nu)
    L = E - V + 1 # loops
    w = nu_tot - L * (D0-2*eps)/2 # superficial degree of divergence
    num = gamma(w)
    den = [gamma(nu[i]) for i in range(E)]
    den = prod(den)
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

    ##########################################
    # Inverting Mandelstams to scalar products
    ##########################################

def set_momenta(graph, pspt):
    #
    # 1-pt function
    #
    if vertices(graph) == [0]:
        return [0]
    #
    # 2-pt or higher
    #
    V = len(vertices(graph)) # total number of vertices (external + internal)
    Vext = get_Vext(pspt)    # number of particles
    sp_mat = zeros(V)        # V x V matrix with pi*pj = 0 if i or j internal vertex
    #
    # first 1, ..., (Vext-1) diagonals
    #
    for i in range(Vext-1):
        sp_mat[i,i] = sp[i,i] 
    #
    # upper triangular (Vext-1) x (Vext-1) block
    #
    for i in range(Vext-1):
        for j in range(i+1, Vext-1):
            sp_mat[i,j] = (s[i,j] - sp[i,i] - sp[j,j]) / 2
    #
    # copy to lower triangular (Vext-1) x (Vext-1) block
    #
    for i in range(Vext-1):
        for j in range(i+1, Vext-1):
            sp_mat[j,i] = sp_mat[i,j]
    #
    # last entry of column = (-1) * everything above 
    #
    for j in range(Vext-1):
         colsum = 0
         for i in range(Vext-1):
             colsum += sp_mat[i,j]
         sp_mat[Vext-1,j] = (-1) * colsum
    #
    # last entry of row = (-1) * everything to the left
    #
    for i in range(Vext):
         rowsum = 0
         for j in range(Vext-1):
             rowsum += sp_mat[i,j]
         sp_mat[i,Vext-1] = (-1) * rowsum
    #
    # subsitute phase space point
    #
    return sp_mat . subs(pspt) . tolist()

    ######################
    # Tropical integration
    ######################

# trop_int is the value of the Feynman integral (without prefactor!)
# Itr is the normalization factor in the tropical measure
def tropical_integration(graph, masses, mom_var, D0, eps_order, Lambda, N):
    print("Prefactor: " + str(prefactor(graph, D0, eps_order)) + ".")
    g = feyntrop.graph(graph)
    pi_pj = set_momenta(graph, mom_var)
    trop_int, Itr  = feyntrop.integrate_graph(g, D0, pi_pj, masses, eps_order, Lambda, N)
    print("")
    print_res(trop_int)
    return trop_int, Itr
