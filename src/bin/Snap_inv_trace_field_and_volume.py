#!/usr/bin/python

# cd triangulations
# Snap_inv_trace_field_and_volume.py -H > census_volumes.csv
# for i in cusped_*/*.trig spun_triangulations/*.trig; do Snap_inv_trace_field_and_volume.py "$i" >> census_volumes.csv; done


import os
import sys
import re
import optparse

this_path, this_file = os.path.split(sys.argv[0])
abs_path = os.path.abspath(this_path)
base_path, this_dir = os.path.split(abs_path)
sys.path.append(base_path)

try:
    from manifold.triangulation import triangulation, read_triangulation_from_file
    from algebra.polynomial import polynomial
    from algebra.pari import number, set_pari_precision, get_pari_error
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

import tempfile
import sys
import subprocess
import re

set_pari_precision(120)

def process_file(trig_file):

    d,in_filename = tempfile.mkstemp()

    tmp_in_file = open(in_filename, 'wb')
    tmp_in_file.write("r f %s\n" % trig_file)
    tmp_in_file.write("set precision 15\n")
    tmp_in_file.write("set digits_printed 120 f\n")
    tmp_in_file.write("print solution_type\n")
    tmp_in_file.write("print volume\n")
    tmp_in_file.write("print complex_volume\n")
    tmp_in_file.write("quit\n")
    tmp_in_file.close()

    l = subprocess.Popen("ulimit -t 1; snap <" + in_filename, shell=True, stdout = subprocess.PIPE)
    s = l.stdout.read()

    if not ("solution type: geometric" in s or "solution type: nongeometric" in s):
        return None


    d,in_filename = tempfile.mkstemp()
    tmp_in_file = open(in_filename, 'wb')
    tmp_in_file.write("r f %s\n" % trig_file)
    tmp_in_file.write("set degree 7\n")
    tmp_in_file.write("set precision 15\n")
    tmp_in_file.write("set digits_printed 120 f\n")
    tmp_in_file.write("compute invariant_trace_field\n")
    tmp_in_file.write("quit\n")
    tmp_in_file.close()

    l = subprocess.Popen("ulimit -t 20; snap <" + in_filename, shell=True, stdout = subprocess.PIPE)
    s2 = l.stdout.read()


    r = re.search("Complex volume: (.*)",s)

    rcvol = '"-"'
    icvol = '"-"'

    if r:
        cvol = number(r.group(1))
        rcvol = str(cvol.real())
        icvol = str(cvol.imag())

    r = re.search("Volume is: (.*)", s)

    if r:
        rcvol = r.group(1)

    f = re.search("Invariant trace field: (.*) \[",s2)
    fp = '"-"'
    fdeg = '"-"'

    if f:
        fp = f.group(1)
        fdeg = str(polynomial(fp).degree())

    
    if re.search(r"\(-?\d+,-?\d+\)", trig_file):
        oc = "closed"
    else:
        oc = "cusped"

    t = read_triangulation_from_file(trig_file)

    return ('"%s",%s,%d,%s,%s,"%s",%s'
            % (trig_file.split('/')[-1].replace('.trig',''),
               oc,
               t.num_tets,
               rcvol, icvol, fp, fdeg))
            
if not len(sys.argv) == 2:
    print "Usage: Snap_inv_trace_field_and_volume.py -H         prints CSV header"
    print "       Snap_inv_trace_field_and_volume.py TRIG_FILE  prints row for CSV file"
    sys.exit(1)

if sys.argv[1] == '-H':
    print 'Name, Type, Tetrahedra, Volume, CS, "Invariant Trace Field", Degree'
else:
    k = process_file(sys.argv[1])
    if k:
        print k

