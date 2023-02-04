import os
import sys
path = os.path.abspath('..')
sys.path.append(path)
from py_feyntrop import *

def save(file, expr):
     file = file + ".txt"
     text_file = open(file, "w")
     text_file.write(str(expr))
     text_file.close()

def print_ex(example):
    string  = "\t# " + example + " #"
    length  = len(string) - 1
    padding = "\t" + "=" * length
    print("\n" + padding)
    print(string)
    print(padding)
    print("")

global_N = int(1e8)

    #========================
    # 2L_4pt_non_planar_muone
    #========================

def ex1():
    example = "2L_4pt_non_planar_muone"
    print_ex(example)
    #
    graph = [((0,1), 1), ((0,4), 1), ((1,5), 1), ((5,2), 1), ((5,3), 1), ((4,3), 1), ((4,2), 1)]
    #
    m, M = 1/200, 1
    S, T = -1/7, -1/3
    U = 2*M + 2*m - S - T 
    momentum_vars = [(p_sqr[0], M), (p_sqr[1], m), (p_sqr[2], m),(s[0,1] , S), (s[1,2] , T), (s[0,2] , U)]
    masses_sqr   = [0, M, m, m, 0, M, 0]
    #
    D0, eps_order, Lambda, N = 6, 5, 1.29, global_N 
    #
    res = tropical_integration(graph, masses_sqr, momentum_vars, D0, eps_order, Lambda, N)[0]
    save(example, res)

    #=================================
    # 2L_5pt_topo_4 (qcd box-pentagon)
    #=================================

def ex2():
    example = "2L_5pt_topo_4"
    print_ex(example)
    #
    graph = [((0,1), 1), ((1,2), 1), ((2,6), 1), ((6,3), 1), ((3,4), 1), ((4,5), 1), ((5,0), 1), ((5,6), 1)]
    #
    m = 1/2
    momentum_vars = [(p_sqr[0], 0), (p_sqr[1], m), (p_sqr[2], m), (p_sqr[3], m), (s[0,1], 2.2), (s[0,2], 2.3), (s[0,3], 2.4), (s[1,2], 2.5), (s[1,3], 2.6), (s[2,3], 2.7)]
    masses_sqr = [0, m, 0, m, 0, m, 0, m]
    #
    D0, eps_order, Lambda, N = 6, 5, 0.28, global_N 
    #
    res = tropical_integration(graph, masses_sqr, momentum_vars, D0, eps_order, Lambda, N)[0]
    save(example, res)

    #================
    # 3L_4pt_envelope
    #================

def ex3():
    example = "3L_4pt_envelope"
    print_ex(example)
    #
    graph = [((0,1), 1), ((1,2), 1), ((2,3), 1), ((3,0), 1), ((0,2), 2), ((1,3), 2)]
    #
    momentum_vars = [(p_sqr[0], 1.1), (p_sqr[1], 1.2), (p_sqr[2], 1.3), (s[0,1] , 2.1), (s[0,2] , 2.2), (s[1,2] , 2.3)]
    masses_sqr = [0.05, 0.06, 0.07, 0.08, 0.09, 0.10]
    #
    D0, eps_order, Lambda, N = 4, 7, 1.24, global_N
    #
    res = tropical_integration(graph, masses_sqr, momentum_vars, D0, eps_order, Lambda, N)[0]
    save(example, res)

    #=========================
    # 4L_0pt_qed_middle_bubble
    #=========================

def ex4():
    example = "4L_0pt_qed_middle_bubble"
    print_ex(example)
    #
    graph = [((0,1), 1), ((1,2), 1), ((2,0), 1), ((0,5), 1), ((1,4), 1), ((2,3), 1), ((3,4), 1), ((4,5), 1), ((5,3), 1)]
    #
    me = 1
    masses_sqr = [me, me, me, 0, 0, 0, me, me, me]
    momentum_vars = [(p_sqr[0], 0)]
    D0, eps_order, Lambda, N = 4, 9, 0, global_N
    #
    res = tropical_integration(graph, masses_sqr, momentum_vars, D0, eps_order, Lambda, N)[0]
    save(example, res)

    #=========
    # 5L_2pt_M
    #=========

def ex5():
    example = "5L_2pt_M"
    print_ex(example)
    #
    graph = [((0,6), 1), ((0,5), 1), ((5,6), 1), ((6,4), 1), ((5,3), 1), ((5,4), 1), ((4,3), 1), ((4,2), 1), ((3,2), 1), ((3,1), 1), ((2,1), 1)]
    #
    masses_sqr = [i for i in range(1,12)]
    momentum_vars = [(p_sqr[0], 100)]
    #
    D0, eps_order, Lambda, N = 3, 11, 0.02, global_N
    #
    res = tropical_integration(graph, masses_sqr, momentum_vars, D0, eps_order, Lambda, N)[0]
    save(example, res)

    #=============
    # RUN EXAMPLES
    #=============

print("\n\n" + "="*35 + " Running examples " + "="*35)
ex1()
ex2()
ex3()
ex4()
ex5()
print("\n\n" + "="*35 + " All examples done "  + "="*35 + "\n")
