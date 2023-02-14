import feyntrop
from math import log10
from sympy import prod, gamma, series, symbols, zeros, IndexedBase

    #=====
    # Misc 
    #=====

# s[i,j] = (p_i+p_j)^2 is used in replacement rules
# eps is used as an expansion parameter 
sp = IndexedBase('sp')
eps = symbols('eps')

# un-nest a list
def flatten(list):
    return [item for sublist in list for item in sublist]

# vertices of an edge list
def vertices(edges):
    edges = [e[0] for e in edges]
    edges = flatten(edges)
    return list(set(edges)) 

# find largest index in momentum_vars to get V_ext = n_particles
# last line: + 2 = 1 (zero indexing) + 1 (only provide data for V_ext-1 particles in momentum_vars)
def get_V_ext(momentum_vars):
    indices = [
        momentum_vars[i][0] . indices for i in range(len(momentum_vars))
    ]
    indices = flatten(indices)
    return max(indices) + 2

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
    print("")
    for i in range(eps_order):
        re , re_err = round_re[i]
        im , im_err = round_im[i]
        eps_str = '-- eps^' + str(i) + ': '
        #
        print(f"{eps_str : <5}{'[' : ^1}{re : ^{max_re}}{' +/- ' : ^5}{re_err : ^{max_re_err}}{']' : ^1}{'  +  i * ' : ^9}{'[' : ^1}{im : ^{max_im}}{' +/- ' : ^5}{im_err : ^{max_im_err}}{']' : ^1}")

# prefactor in Symanzik representation
def prefactor(edges, D0, eps_order):
    V = len(vertices(edges))
    E = len(edges)
    nu = [e[1] for e in edges] # propagator powers
    nu_tot = sum(nu)
    L = E - V + 1 # loops
    w = nu_tot - L * (D0-2*eps)/2 # superficial degree of divergence
    num = gamma(w)
    den = prod([gamma(nu[i]) for i in range(E)])
    return num/den

# mulitplying result by prefactor and expanding in epsilon
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

    #=======================
    # Prepare kinematic data
    #=======================

# convert phase_space_point into a dictionay with it being optional whether 'val' is string or not
def to_dict(phase_space_point):
    n = len(phase_space_point)
    temp = []
    for i in range(n):
        var, val = phase_space_point[i][0], phase_space_point[i][1]
        temp . append( (var, str(val)) )
    return dict(temp)

# replace all elements from a dictionary for expressions such as 's + t' containing more than one variable
def replace_all(expr, dic):
    for i, j in dic . items():
        expr = expr . replace(i, j)
    return expr

# replace sp[i,j] in replacement_rules with the chosen values in phase_space_point
def replace_sp(replacement_rules, phase_space_point):
    n = len(replacement_rules)
    phase_space_point = to_dict(phase_space_point)
    temp = []
    for i in range(n):
        var, val  = replacement_rules[i][0], replacement_rules[i][1]
        val = replace_all(str(val), phase_space_point)
        val = eval(val)
        temp . append( (var, val) )
    return temp

# list of squared masses
def get_m_sqr(edges, phase_space_point):
    E = len(edges)
    phase_space_point = to_dict(phase_space_point)
    temp = []
    for e in range(E):
        val = edges[e][2]
        val = replace_all(str(val), phase_space_point)
        val = eval(val)
        temp . append(val)
    return temp

# V x V matrix of scalar products
def get_P_uv(edges, phase_space_point, replacement_rules):
   V = len(vertices(edges)) # total number of vertices (external and internal)
   P_uv_mat = zeros(V)
   #
       # 1-pt function: all scalar products are zero
   #
   if replacement_rules == []:
       return P_uv_mat . tolist()
   #
       # 2-pt or higher
   #
   V_ext = get_V_ext(replacement_rules) # number of external particles
   n = V_ext - 1 # to simplify for-loops
   #
   # diaonal and upper triangular part
   #
   for i in range(n):
       for j in range(i, n):
           P_uv_mat[i,j] = sp[i,j]
   #
   # copy to upper to lower triangular part
   #
   for i in range(n):
       for j in range(i+1, n):
           P_uv_mat[j,i] = P_uv_mat[i,j]
   #
   # last entry of column = (-1) * everything above 
   #
   for j in range(n):
        colsum = 0
        for i in range(n):
            colsum += P_uv_mat[i,j]
        P_uv_mat[n,j] = (-1) * colsum
   #
   # last entry of row = (-1) * everything to the left
   #
   for i in range(V_ext):
        rowsum = 0
        for j in range(n):
            rowsum += P_uv_mat[i,j]
        P_uv_mat[i,n] = (-1) * rowsum
   #
   # subsitute phase space point
   #
   rep = replace_sp(replacement_rules, phase_space_point)
   return P_uv_mat . subs(rep) . tolist()

    #=====================
    # Tropical integration
    #=====================

# returns matrix of scalar products and a list of squared masses
def prepare_kinematic_data(edges, replacement_rules, phase_space_point):
    P_uv = get_P_uv(edges, phase_space_point, replacement_rules)
    m_sqr = get_m_sqr(edges, phase_space_point)
    return P_uv, m_sqr

# trop_int is the value of the Feynman integral (without prefactor!)
# Itr is the normalization factor in the tropical measure
def tropical_integration(N, D0, Lambda, eps_order, edges, replacement_rules, phase_space_point):
    P_uv, m_sqr = prepare_kinematic_data(edges, replacement_rules, phase_space_point)
    edges = [e[:-1] for e in edges]
    graph = feyntrop.graph(edges)
    print("Prefactor: " + str(prefactor(edges, D0, eps_order)) + ".")
    trop_res, Itr  = feyntrop.integrate_graph(graph, D0, P_uv, m_sqr, eps_order, Lambda, N)
    print_res(trop_res)
    return trop_res, Itr
