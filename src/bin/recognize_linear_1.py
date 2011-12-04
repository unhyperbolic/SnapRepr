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
    from algebra.solve_polynomial_equations import solve_polynomial_equations
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

if not len(sys.argv)==4:
    print "Usage: recognize_linear_1 census_volumes.csv representation_volumes.csv linear_combinations.csv"
    sys.exit(1)

def read_census_volumes(csv_file):
    CensusVolumesReader = csv.DictReader(csv_file)

    CensusVolumesDict = []

    for d in CensusVolumesReader:
        d['Volume'] = number(d['Volume'])
        d['Tetrahedra'] = int(d['Tetrahedra'])
        CensusVolumesDict.append(d)

    CensusVolumesDict.sort(key = lambda d: d['Volume'])
    return CensusVolumesDict

def add_submultiples_to_list(D):
    N = []
    for d in D:
        for i in range(1, 1 + int(d['Volume'] / number('0.9'))):
            e = d.copy()
            e['Volume'] = d['Volume'] / number(i)
            e['VolumeFactor'] = i
            N.append(e)
    N.sort(key = lambda d:d['Volume'])
    return N

def find_representative_census_volumes(L):
    def admit(k):
        for l in k:
            if l['VolumeFactor'] == 1:
                return True
        return False

    def reps(k):        
        def all_factor_one(k):
            return [x for x in k if x['VolumeFactor'] == 1]
        def link_complement(k):
            k = [(
                    ('^' in d['Name'],int(d['Name'].split('_')[0].split('^')[0]),d['Name']),
                    d
                 )
                 for d in k if '_' in d['Name']]
            k.sort()
            return [x[1] for x in k[:1]]
        def cusped_census_mfd(k):
            k = [(
                    (d['Tetrahedra'], int(d['Name'][1:])),
                    d
                 )
                 for d in k if not ('_' in d['Name'] or '(' in d['Name'])]
            k.sort()
            return [x[1] for x in k[:1]]
        def closed_census_mfd(k):
            k = [(
                    (d['Tetrahedra'], int(d['Name'][1:].split('(')[0])),
                    d
                 )
                 for d in k if '(' in d['Name']]
            k.sort()
            return [x[1] for x in k[:1]]

        as_cusped_mfd = link_complement(all_factor_one(k))
        if not as_cusped_mfd:
            as_cusped_mfd = cusped_census_mfd(all_factor_one(k))
        if not as_cusped_mfd:
            as_cusped_mfd = link_complement(k)
        if not as_cusped_mfd:
            as_cusped_mfd = cusped_census_mfd(k)
        
        as_closed_mfd = closed_census_mfd(all_factor_one(k))
        if not as_closed_mfd:
            as_closed_mfd = closed_census_mfd(k)

        return as_cusped_mfd + as_closed_mfd

    return [reps(x) for x in L if admit(x)]

def find_equivalent_census_volumes(D):
    vols = []
    prev_vol = number("-1.0")
    mfds = []

    new_list = []

    for d in D:
        if (d['Volume']-prev_vol).abs() > number("1E-30"):
            if mfds:
                new_list.append(mfds)
            mfds = []
            prev_vol = d['Volume']
        mfds.append(d)
    new_list.append(mfds)
    return new_list

CensusVolumesDict = read_census_volumes(open(sys.argv[1],'rb'))
CensusVolumesDict = add_submultiples_to_list(CensusVolumesDict)
CensusVolumesDict = find_equivalent_census_volumes(CensusVolumesDict)
CensusVolumesDict = find_representative_census_volumes(CensusVolumesDict)

print float(CensusVolumesDict[0][0]['Volume'])

find_lin_comb.initialize(CensusVolumesDict)

#for i in CensusVolumesDict:
#    print
#    for j in i:
#        print j

RepresentationHeader = csv.reader(open(sys.argv[2],'rb')).next()
RepresentationReader = csv.DictReader(open(sys.argv[2],'rb'))

OutWriter = csv.DictWriter(open(sys.argv[3],'wb'),RepresentationHeader + ['Linear Combinations'])

for d in RepresentationReader:
    d['Linear Combinations'] = '-'
    if not d['Volume'] == '-':
        v = number(d['Volume'])
        vfloat = float(d['Volume'])
        if vfloat > 1e-10:
            print vfloat, v
            print find_lin_comb.find_as_multiple(v)
            print find_lin_comb.find_as_linear_combination(v)

            

            

    OutWriter.writerow(d)

