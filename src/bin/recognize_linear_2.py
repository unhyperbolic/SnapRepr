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
    from algebra.pari import number, set_pari_precision, set_pari_allowed_error, get_pari_allowed_error, NumericalError, pari_eval
    from algebra.solve_polynomial_equations import solve_polynomial_equations
    import manifold.slN
    import algebra.magma
    import algebra.pari
    import hashlib
except ImportError as e:
    print e
    print
    print "This program was called as       :", sys.argv[0]
    print "Absolute path to this program is :", abs_path
    print "Base path is                     :", base_path
    sys.exit(1)

if not len(sys.argv)==4:
    print "Usage: recognize_linear_2 representation_volumes.csv linear_combinations.csv PRECISION"
    sys.exit(1)

volDict = {
    "4_1" :
        number("2.02988321281930725004240510854904057188337861506059958403497821355319495251648804427294070845651338989172365506271977525"),
    "m003(-3,1)" :
        number("0.942707362776927720921299603092211647590327105766883159014506775752934182774157210312315672643333035804180429759598797803"),
    "m003(-2,3)" :
        number("0.981368828892232088091452189794427068238164321906312438642604199777420481646215962007743465608622970924477159450932847661"),
    "m006(1,3)" :
        number("1.83193118835443803010920702986476822154829874856334426853299623924352603955250953895871302585223021249714884523839239916"),
    "m009" :
        number("2.66674478344905979079671246261065004409838388855263953139317180331572348841561131881281221620442212850291034276448127621"),
    "m004(6,1)":
        number("1.28448530046835444246033708486871125891241877318392809214943701516964507690587972156461082746727624238071994448327461536"),
    "m007(3,2)":
        number("1.58316666062481283616602885187919062839805259468511564120510037138514921358016328417586218367685774975389015282528735199"),
    "9^3_15" :
        number("8.35550214637956599430377376881477535438602570433613182225567019065904209089981277038106183143111474247020120487272572689"),
    "8_16" :
        number("10.5790219168992703258086845500886728985048852795970808175018610806349105463038041532838496386087925737269710231234852770")
    }

precision = int(sys.argv[3])

RepresentationHeader = csv.reader(open(sys.argv[1],'rb')).next() + ['Linear Combinations']

RepresentationReader = csv.DictReader(open(sys.argv[1],'rb'))

OutWriter = csv.DictWriter(open(sys.argv[2],'wb'), RepresentationHeader)
OutWriter.writerow( dict(zip(RepresentationHeader, 
                             RepresentationHeader)))

for d in RepresentationReader:
    d['Linear Combinations'] = '-'
    if d['Volume'] and not d['Volume'] == '-':
        v = number(d['Volume'])
        vfloat = float(d['Volume'])
        if vfloat > 1e-10:

            pari_cmd = ("lindep([%s],%d)" 
                          % (','.join([str(x) 
                                       for x in 
                                       volDict.values() + [v]]),
                             precision))
            s = pari_eval(pari_cmd) 
            if d['Manifold'] == '4_1' and False:
                print s
                print v
                print len(volDict.values() + [v])
            assert s[0] == '['
            assert s[-2:] == ']~'
            s = [int(x) for x in s[1:-2].split(',')]
            k = max([max(s),-min(s)])
            if k < 1000:
                q = '('
                q+= '+'.join(['(%d) * %s' % x for x in zip(s[:-1], volDict.keys()) if x[0]])
                q+= ') / (%d)' % s[-1]

                d['Linear Combinations'] = q

    OutWriter.writerow(d)

