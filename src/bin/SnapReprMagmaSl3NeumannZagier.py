#!/usr/bin/python

import os
import sys

this_path, this_file = os.path.split(sys.argv[0])
abs_path = os.path.abspath(this_path)
base_path, this_dir = os.path.split(abs_path)
sys.path.append(base_path)

try:
    from manifold import sl3NeumannZagierType
    from manifold.triangulation import read_triangulation_from_file
    import algebra.magma
except ImportError as e:
    print e
    print
    print "This program was called as       :", sys.argv[0]
    print "Absolute path to this program is :", abs_path
    print "Base path is                     :", base_path
    sys.exit(1)

def get_term_order(polys, pre_vars = [], post_vars = []):
    all_vars = sum([p.variables() for p in polys], [])
    sort_vars = set(all_vars) - set(pre_vars) - set(post_vars)
    sort_vars = list(sort_vars)
    sort_vars.sort()

    return pre_vars + sort_vars + post_vars


def produce_magma_out(trig):

    eqns = sl3NeumannZagierType.produce_all_equations_non_degenerate(trig)
    term_order = get_term_order(eqns, pre_vars = ['t'])
    return algebra.magma.primary_decomposition(eqns, term_order = term_order)

def main():
    
    trig_filename = sys.argv[1]

    if trig_filename[-5:] == '.trig':
        base_filename = trig_filename[:-5]
    else:
        base_filename = trig_filename

    trig = read_triangulation_from_file(trig_filename)

    open(base_filename+'_sl3NeumannZagier.magma','w').write(
        produce_magma_out(trig))

main()
