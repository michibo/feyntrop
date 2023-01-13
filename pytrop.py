import feyntrop
import re
from math import  prod
from sympy import gamma, series, symbols, zeros, Matrix, IndexedBase

sp  = IndexedBase('sp')
s   = IndexedBase('s')
eps = symbols('eps')

    ######
    # MISC
    ######

# un-nest a list
def flatten(list):
    return [item for sublist in list for item in sublist]

# vertices of a graph
def vertices(edges):
    edges = [e[0] for e in edges] # remove edge weights
    edges = flatten(edges)
    return list(set(edges)) # unique integers in edges

# find largest index in pspt to get Vext = n_particles
# last line: + 2 = 1 (zero indexing) + 1 (only provide data for Vext-1 particles in pspt)
def get_Vext(pspt):
    indices = [
        pspt[i][0] . indices for i in range(len(pspt))
    ]
    indices = flatten(indices)
    return max(indices) + 2

    ################
    # FORMATTING RES
    ################

# format epsilon terms into tuples with associated error
def num_err_str(number_error):
    number, error = number_error
    return "%f +/- %f" % (number, error)

# print epsilon terms (without prefactor)
def print_res(res):
    for i, term in enumerate(res):
        real, imag = term
        string = "eps^%d : [%s] + I * [%s]" % (i, num_err_str(real), num_err_str(imag))
        string = re.sub('\[+(\d)' , '[(+1) * \\1' , string)
        string = re.sub('\[-(\d)' , '[(-1) * \\1' , string)
        print(string)

# prefactor of Symanzik representation
def prefactor(edges, D0, eps_order):
    V = len(vertices(edges))
    E = len(edges)
    nu = [e[1] for e in edges] # propagator powers
    nutot = sum(nu)
    L = E - V + 1 # loops
    sign = (-1)**nutot
    w = nutot - L * (D0-2*eps)/2 # superficial degree of divergence
    num = sign * gamma(w)
    den = [gamma(nu[i]) for i in range(E)]
    den = prod(den)
    return num/den

# mulitplying res by prefactor and expanding in epsilon
def eps_expansion(res, edges, D0):
    eps_order = len(res)
    terms = 0
    for i in range(eps_order):
        re = res[i][0][0]
        im = res[i][1][0]
        terms += complex(re, im)*eps**i
    pf = prefactor(edges, D0, eps_order)
    terms = pf * terms
    terms = series(terms, eps, 0, eps_order)
    return terms . evalf(10)

    ##########################################
    # INVERTING MANDELSTAMS TO SCALAR PRODUCTS
    ##########################################

def set_momenta(edges, pspt):
    #
    # 1-pt function
    #
    if vertices(edges) == [0]:
        return [0]
    #
    # 2-pt or higher
    #
    V = len(vertices(edges)) # total number of vertices (external + internal)
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

# trop_int is the value of the Feynman integral (without prefactor)
# Itr is the normalization factor in the tropical measure
def tropical_integration(graph, masses, mom_var, D0, eps_order, Lambda, N):
    print("\nprefactor:")
    print(prefactor(graph, D0, eps_order))
    print("")
    g = feyntrop.graph(graph)
    pi_pj = set_momenta(graph, mom_var)
    trop_int, Itr  = feyntrop.integrate_graph(g, D0, pi_pj, masses, eps_order, Lambda, N)
    print("")
    print_res(trop_int)
    return trop_int, Itr

