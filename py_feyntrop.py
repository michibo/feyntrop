
from math import log10
from sympy import prod, gamma, series, symbols, zeros, IndexedBase

import json
import subprocess
import os

# ======
#  Misc
# ======

# s[i,j] = (p_i+p_j)^2 is used in replacement rules
# eps is used as an expansion parameter
sp = IndexedBase('sp')
eps = symbols('eps')


def flatten(list):
    """
    un-nest a list
    """
    return [item for sublist in list for item in sublist]


def vertices(edges):
    """
    vertices of an edge list
    """
    edges = [e[0] for e in edges]
    edges = flatten(edges)
    return list(set(edges))


def get_V_ext(momentum_vars):
    """
    find largest index in momentum_vars to get V_ext = n_particles
    last line: + 2 = 1 (zero indexing) + 1 (only provide data for V_ext-1 particles in momentum_vars)
    """
    indices = [
        momentum_vars[i][0].indices for i in range(len(momentum_vars))
    ]
    indices = flatten(indices)
    return max(indices) + 2

# ===================
#  Formatting result
# ===================


def round_res(res_err_list):
    """
    round result to two significant digits
    """
    res = res_err_list[0]
    err = res_err_list[1]
    if err == 0:
        return [str(res), str(err)]
    nzeros = abs(int(log10(err)))  # number of zeros after decimal pt.
    nround = nzeros + 2
    res = round(res, nround)
    err = round(err, nround)
    err = f'{err:.{nround}f}'
    res = f'{res:.{nround}f}'
    return [res, err]


def print_res(res):
    """
    print each epsilon order with rounded result
    """
    round_re, round_im = [], []
    lens_re, lens_im = [], []
    lens_re_err, lens_im_err = [], []
    eps_order = len(res)
    for i in range(eps_order):
        #
        round_re.append(round_res(res[i][0]))
        round_im.append(round_res(res[i][1]))
        #
        lens_re    .append(len(round_re[i][0]))
        lens_re_err.append(len(round_re[i][1]))
        #
        lens_im    .append(len(round_im[i][0]))
        lens_im_err.append(len(round_im[i][1]))
        #
    max_re = str(max(lens_re))
    max_im = str(max(lens_im))
    max_re_err = str(max(lens_re_err))
    max_im_err = str(max(lens_im_err))
    print("")
    for i in range(eps_order):
        re, re_err = round_re[i]
        im, im_err = round_im[i]
        eps_str = '-- eps^' + str(i) + ': '
        #
        print(f"{eps_str : <5}{'[' : ^1}{re : ^{max_re}}{' +/- ' : ^5}{re_err : ^{max_re_err}}{']' : ^1}{'  +  i * ' : ^9}{'[' : ^1}{im : ^{max_im}}{' +/- ' : ^5}{im_err : ^{max_im_err}}{']' : ^1}")


def prefactor(edges, D0, eps_order):
    """
    prefactor in Symanzik representation
    """
    V = len(vertices(edges))
    E = len(edges)
    nu = [e[1] for e in edges]  # propagator powers
    nu_tot = sum(nu)
    L = E - V + 1  # loops
    w = nu_tot - L * (D0 - 2 * eps) / 2  # superficial degree of divergence
    num = gamma(w)
    den = prod([gamma(nu[i]) for i in range(E)])
    return num / den


def eps_expansion(res, edges, D0):
    """
    multiplying result by prefactor and expanding in epsilon
    """
    eps_order = len(res)
    terms = 0
    for i, resi in enumerate(res):
        re = resi[0][0]
        im = resi[1][0]
        terms += complex(re, im) * eps**i
    pf = prefactor(edges, D0, eps_order)
    terms = pf * terms
    terms = series(terms, eps, 0, eps_order)
    return terms.evalf(10)

# ========================
#  Prepare kinematic data
# ========================


def to_dict(phase_space_point):
    """
    convert phase_space_point into a dictionary with it being optional whether 'val' is string or not
    """
    n = len(phase_space_point)
    temp = {}
    for i in range(n):
        var, val = phase_space_point[i][0], phase_space_point[i][1]
        temp[var] = str(val)
    return temp


def replace_all(expr, dic):
    """
    replace all elements from a dictionary for expressions such as 's + t' containing more than one variable
    """
    for i, j in dic.items():
        expr = expr.replace(i, j)
    return expr


def replace_sp(replacement_rules, phase_space_point):
    """
    replace sp[i,j] in replacement_rules with the chosen values in phase_space_point
    """
    n = len(replacement_rules)
    phase_space_point = to_dict(phase_space_point)
    temp = []
    for i in range(n):
        var, val = replacement_rules[i][0], replacement_rules[i][1]
        val = replace_all(str(val), phase_space_point)
        val = eval(val)
        temp.append((var, val))
    return temp


def get_m_sqr(edges, phase_space_point):
    """
    list of squared masses
    """
    E = len(edges)
    phase_space_point = to_dict(phase_space_point)
    temp = []
    for e in range(E):
        val = edges[e][2]
        val = replace_all(str(val), phase_space_point)
        val = eval(val)
        temp.append(val)
    return temp


def get_P_uv(edges, phase_space_point, replacement_rules):
    """
    # V x V matrix of scalar products
    """
    V = len(vertices(edges))  # total number of vertices (external and internal)
    P_uv_mat = zeros(V)
    #
    # 1-pt function: all scalar products are zero
    #
    if not replacement_rules:
        return P_uv_mat.tolist()
    #
    # 2-pt or higher
    #
    V_ext = get_V_ext(replacement_rules)  # number of external particles
    n = V_ext - 1  # to simplify for-loops
    #
    # diagonal and upper triangular part
    #
    for i in range(n):
        for j in range(i, n):
            P_uv_mat[i, j] = sp[i, j]
    #
    # copy to upper to lower triangular part
    #
    for i in range(n):
        for j in range(i + 1, n):
            P_uv_mat[j, i] = P_uv_mat[i, j]
    #
    # last entry of column = (-1) * everything above
    #
    for j in range(n):
        colsum = 0
        for i in range(n):
            colsum += P_uv_mat[i, j]
        P_uv_mat[n, j] = (-1) * colsum
    #
    # last entry of row = (-1) * everything to the left
    #
    for i in range(V_ext):
        rowsum = 0
        for j in range(n):
            rowsum += P_uv_mat[i, j]
        P_uv_mat[i, n] = (-1) * rowsum
    #
    # substitute phase space point
    #
    rep = replace_sp(replacement_rules, phase_space_point)
    return P_uv_mat.subs(rep).tolist()

# ======================
#  Tropical integration
# ======================


def prepare_kinematic_data(edges, replacement_rules, phase_space_point):
    """
    returns matrix of scalar products and a list of squared masses
    """
    P_uv = get_P_uv(edges, phase_space_point, replacement_rules)
    m_sqr = get_m_sqr(edges, phase_space_point)
    return P_uv, m_sqr


def tropical_integration(N, D0, Lambda, eps_order, edges, replacement_rules, phase_space_point):
    """
    trop_int is the value of the Feynman integral (without prefactor!)

    Itr is the normalization factor in the tropical measure
    """
    P_uv, m_sqr = prepare_kinematic_data(edges, replacement_rules, phase_space_point)
    edges = [e[:-1] for e in edges]

    P_uv = [[float(e) for e in row] for row in P_uv]
    ft_input = {
        "graph": edges,
        "dimension": D0,
        "scalarproducts": P_uv,
        "masses_sqr": m_sqr,
        "num_eps_terms": eps_order,
        "lambda": Lambda,
        "N": N,
        "seed": 0
    }

    json_str = json.dumps(ft_input)

    if not os.path.exists("./feyntrop"):
        raise RuntimeError("The 'feyntrop' binary file was not found. Please make sure it is compiled and available in the current working directory.")

    print("Starting integration using feyntrop with input:")
    print("Graph with edge weights:", edges)
    print("Dimension:", D0)
    print("Scalarproducts (matrix element (u,v) is the scalar product of ext. momentum flowing into vertices u and v):\n", P_uv)
    print("Squared masses:", m_sqr)
    print("Epsilon order:", eps_order)
    print("Deformation parameter Lambda:", Lambda)
    print("Sample points:", N)

    p = subprocess.Popen(["./feyntrop"], stdout=subprocess.PIPE, stdin=subprocess.PIPE, encoding='utf8')
    out, err = p.communicate(json_str)

    if 0 != p.returncode:
        raise RuntimeError("An error occurred in the feyntrop core-code. No output was generated.")

    output = json.loads(out)

    trop_res = output["integral"]
    Itr = output["IGtr"]

    print("Prefactor: " + str(prefactor(edges, D0, eps_order)) + ".")
    print_res(trop_res)
    return trop_res, Itr
