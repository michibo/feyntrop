import os
import sys
path = os.path.abspath('..')
sys.path.append(path)
from py_feyntrop import *

    #########
    # Utils #
    #########

# print name of topology 
def print_ex(example):
    string = "\t# Testing " + example + " #"
    length = len(string) - 1
    padding = "\t" + "-" * length
    print("\n" + padding)
    print(string)
    print(padding + "\n")

# read pre-computed results
def from_file(example):
    with open(example + ".txt") as file:
        res = file.read()
        res = eval(res)
    return res

# test if two numbers are close
def is_close(res_1, res_2):
    return round((abs(res_1)+1) / (abs(res_2)+1), 2)

# compare tests with files
def compare_w_file(res_1, example):
    res_2 = from_file(example)
    print("\nratios w.r.t. file (should be close to 1):\n")
    for i in range(len(res_1)):
        re_1 , re_2 = res_1[i][0][0] , res_2[i][0][0]
        im_1 , im_2 = res_1[i][1][0] , res_2[i][1][0]
        is_close_re = is_close(re_1, re_2)
        is_close_im = is_close(im_1, im_2)
        print("\teps^" + str(i) + ": ["  + str(is_close_re) + "] + [" + str(is_close_im) + "] * i")
    print("\n" + "_"*90 + "\n")

    #########
    # Tests #
    #########

def test_0():
    example = "1L_2pt"
    print_ex(example)
    graph = [((0,1), 1), ((0,1), 2)]
    momentum_vars = [(p_sqr[0], 3)]
    masses_sqr = [0.2, 1.2]
    D0, eps_order, Lambda, N = 4, 5, 1, int(1e7)
    res = tropical_integration(graph, masses_sqr, momentum_vars, D0, eps_order, Lambda, N)[0]
    compare_w_file(res, example)

def test_1():
    example = "1L_3pt_eucl"
    print_ex(example)
    graph = [((0, 1), 1), ((1, 2), 1), ((2, 0), 1)]
    momentum_vars = [(p_sqr[0], -4), (p_sqr[1], -5), (s[0,1], -13)]
    masses_sqr = [0.2, 1.2, 2.2]
    D0, eps_order, Lambda, N = 4, 5, 1, int(1e7)
    res = tropical_integration(graph, masses_sqr, momentum_vars, D0, eps_order, Lambda, N)[0]
    compare_w_file(res, example)

def test_2():
    example = "1L_3pt_mink"
    print_ex(example)
    graph = [((0, 1), 1), ((1, 2), 1), ((2, 0), 1)]
    momentum_vars = [(p_sqr[0], 1), (p_sqr[1], 4), (s[0,1], 3)]
    masses_sqr = [0.01, 0.5, 1]
    D0, eps_order, Lambda, N = 4, 5, 1, int(1e7)
    res = tropical_integration(graph, masses_sqr, momentum_vars, D0, eps_order, Lambda, N)[0]
    compare_w_file(res, example)

def test_3():
    example = "2L_2pt_ladder"
    print_ex(example)
    graph  = [((0,2), 1), ((0,3), 1), ((1,2), 1), ((1,3), 1), ((2,3), 1)]
    momentum_vars  = [(p_sqr[0], 3)]
    masses_sqr = [0.1, 0.2, 0.3, 0.4, 0.5]
    D0, eps_order, Lambda, N = 4, 4, 0.5, int(1e7)
    res = tropical_integration(graph, masses_sqr, momentum_vars, D0, eps_order, Lambda, N)[0]
    compare_w_file(res, example)

def test_4():
    example = "2L_4pt_ladder"
    print_ex(example)
    graph = [((0,1), 1), ((0,4), 1), ((1,5), 1), ((5,4), 1), ((5,2), 1), ((4,3), 1), ((3,2), 1)]
    momentum_vars = [(p_sqr[0], 1.2), (p_sqr[1], 1.3), (p_sqr[2], 1.4), (s[0,1], 2.1), (s[0,2], 2.2), (s[1,2], 2.3)]
    masses_sqr = [0] * len(graph)
    D0, eps_order, Lambda, N = 4, 5, 1, int(1e7)
    res = tropical_integration(graph, masses_sqr, momentum_vars, D0, eps_order, Lambda, N)[0]
    compare_w_file(res, example)

    #############
    # Run tests #
    #############

print("\n" + "="*35 + " Running test suite " + "="*35)
test_0()
test_1()
test_2()
test_3()
test_4()
