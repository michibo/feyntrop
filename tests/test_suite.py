import os
import sys
path = os.path.abspath('..')
sys.path.append(path)
from py_feyntrop import *

    #========
    # Utils #
    #========

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
    print("\n" + "_"*80 + "\n")

    #========
    # Tests #
    #========

def test_0():
    example = "1L_2pt"
    print_ex(example)
    edges = [((0,1), 1, 'mm0'), ((0,1), 2, 'mm1')]
    replacement_rule = [(sp[0,0], 'pp')]
    phase_space_point = [('mm0', 0.2), ('mm1', 1.2), ('pp', 3)]
    D0, eps_order, Lambda, N = 4, 5, 1, int(1e7)
    trop_res = tropical_integration(N, D0, Lambda, eps_order, edges, replacement_rule, phase_space_point)[0]
    compare_w_file(trop_res, example)

def test_1():
    example = "1L_3pt_eucl"
    print_ex(example)
    edges = [((0, 1), 1, 'mm0'), ((1, 2), 1, 'mm1'), ((2, 0), 1, 'mm2')]
    replacement_rules = [(sp[0,0], 'pp0'), (sp[1,1], 'pp1'), (sp[0,1], '(s01-pp0-pp1)/2')]
    phase_space_point = [('mm0', 0.2), ('mm1', 1.2), ('mm2', 2.2), ('pp0', -4), ('pp1', -5), ('s01', -13)]
    D0, eps_order, Lambda, N = 4, 5, 1, int(1e7)
    trop_res = tropical_integration(N, D0, Lambda, eps_order, edges, replacement_rules, phase_space_point)[0]
    compare_w_file(trop_res, example)

def test_2():
    example = "1L_3pt_mink"
    print_ex(example)
    edges = [((0, 1), 1, 'mm0'), ((1, 2), 1, 'mm1'), ((2, 0), 1, 'mm2')]
    replacement_rules = [(sp[0,0], 'pp0'), (sp[1,1], 'pp1'), (sp[0,1], '(s01-pp0-pp1)/2')]
    phase_space_point = [('mm0', 0.01), ('mm1', 0.5), ('mm2', 1), ('pp0', 1), ('pp1', 4), ('s01', 3)]
    D0, eps_order, Lambda, N = 4, 5, 1, int(1e7)
    trop_res = tropical_integration(N, D0, Lambda, eps_order, edges, replacement_rules, phase_space_point)[0]
    compare_w_file(trop_res, example)

def test_3():
    example = "2L_2pt_ladder"
    print_ex(example)
    edges  = [((0,2), 1, 'mm0'), ((0,3), 1, 'mm1'), ((1,2), 1, 'mm2'), ((1,3), 1, 'mm3'), ((2,3), 1, 'mm4')]
    replacement_rule = [(sp[0,0], 'pp')]
    phase_space_point = [('mm0', 0.1), ('mm1', 0.2), ('mm2', 0.3), ('mm3', 0.4), ('mm4', 0.5), ('pp', 3)]
    D0, eps_order, Lambda, N = 4, 4, 0.5, int(1e7)
    trop_res = tropical_integration(N, D0, Lambda, eps_order, edges, replacement_rule, phase_space_point)[0]
    compare_w_file(trop_res, example)

def test_4():
    example = "2L_4pt_ladder"
    print_ex(example)
    edges = [((0,1), 1, '0'), ((0,4), 1, '0'), ((1,5), 1, '0'), ((5,4), 1, '0'), ((5,2), 1, '0'), ((4,3), 1, '0'), ((3,2), 1, '0')]
    replacement_rules = [(sp[0,0], 'pp0'), (sp[1,1], 'pp1'), (sp[2,2], 'pp2'), (sp[0,1], '(s01-pp0-pp1)/2'), (sp[0,2], '(s02-pp0-pp2)/2'), (sp[1,2], '(s12-pp1-pp2)/2')]
    phase_space_point = [('pp0', 1.2), ('pp1', 1.3), ('pp2', 1.4), ('s01', 2.1), ('s02', 2.2), ('s12', 2.3)]
    D0, eps_order, Lambda, N = 4, 5, 1, int(1e7)
    trop_res = tropical_integration(N, D0, Lambda, eps_order, edges, replacement_rules, phase_space_point)[0]
    compare_w_file(trop_res, example)

def test_5():
    example = "2L_3pt_tutorial"
    print_ex(example)
    edges = [((0,1), 1, 'mm'), ((1,3), 1, 'mm'), ((3,2), 1, 'mm'), ((2,0), 1, 'mm'), ((0,3), 1, 'mm')]
    replacement_rules = [(sp[0,0], '0'), (sp[1,1], '0'), (sp[0,1], 's01/2')]
    phase_space_point = [('mm', 0.2), ('s01', 1)]
    D0, eps_order, Lambda, N = 2, 5, 7.6, int(1e7)
    trop_res = tropical_integration(N, D0, Lambda, eps_order, edges, replacement_rules, phase_space_point)[0]
    compare_w_file(trop_res, example)

def test_6():
    example = "2L_5pt_box_pent"
    print_ex(example)
    edges = [((0,1), 1, '0'), ((1,2), 1, 'mm'), ((2,6), 1, '0'), ((6,3), 1, 'mm'), ((3,4), 1, '0'), ((4,5), 1, 'mm'), ((5,0), 1, '0'), ((5,6), 1, 'mm')]
    replacement_rules = [
        (sp[0,0], '0'), (sp[1,1], 'mm'), (sp[2,2], 'mm'), (sp[3,3], 'mm'),(sp[0,1], '(s01-mm)/2'), (sp[0,2], '(s02-mm)/2'), (sp[0,3], '(s03-mm)/2'),(sp[1,2], '(s12-2*mm)/2'), (sp[1,3], '(s13-2*mm)/2'), (sp[2,3], '(s23-2*mm)/2')
    ]
    phase_space_point = [('mm',  1/2), ('s01', 2.2), ('s02', 2.3), ('s02', 2.3), ('s03', 2.4), ('s12', 2.5), ('s13', 2.6), ('s23', 2.7)]
    D0, eps_order, Lambda, N = 6, 5, 0.28, int(1e7)
    trop_res = tropical_integration(N, D0, Lambda, eps_order, edges, replacement_rules, phase_space_point)[0]
    compare_w_file(trop_res, example)

    #============
    # Run tests #
    #============

print("\n" + "="*30 + " Running test suite " + "="*30)
test_0()
test_1()
test_2()
test_3()
test_4()
test_5()
test_6()
