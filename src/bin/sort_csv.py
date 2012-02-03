#!/usr/bin/python

import os
import sys
import re
import optparse
import traceback
import csv

this_path, this_file = os.path.split(sys.argv[0])
abs_path = os.path.abspath(this_path)
base_path, this_dir = os.path.split(abs_path)
sys.path.append(base_path)

try:
    from manifold.triangulation import triangulation, read_triangulation_from_file
    from algebra.polynomial import polynomial
    from algebra.pari import number, set_pari_precision, set_pari_allowed_error, get_pari_allowed_error, NumericalError
    import manifold.slN
    import algebra.magma
    import algebra.pari
    import hashlib
    import linear_combinations.find_as_multiple_sum_diff as find_lin_comb
except ImportError as e:
    print e
    print
    print "This program was called as       :", sys.argv[0]
    print "Absolute path to this program is :", abs_path
    print "Base path is                     :", base_path
    sys.exit(1)

set_pari_precision(120)

RepresentationHeader = csv.reader(open(sys.argv[1],'rb')).next()
RepresentationReader = csv.DictReader(open(sys.argv[1],'rb'))
OutWriter = csv.DictWriter(open(sys.argv[2],'wb'),RepresentationHeader)


L = [d for d in RepresentationReader if not d['Invariant Trace Field']=='-']

L.sort(key=lambda x:(x['Invariant Trace Field'],number(x['Volume'])))

OutWriter.writerow(dict(zip(RepresentationHeader,RepresentationHeader)))

for d in L:
    OutWriter.writerow(d)
