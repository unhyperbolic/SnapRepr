#!/usr/bin/python
import optparse
import os
import sys
import re

this_path, this_file = os.path.split(sys.argv[0])
abs_path = os.path.abspath(this_path)
base_path, this_dir = os.path.split(abs_path)
sys.path.append(base_path)

try:
    from manifold.triangulation import triangulation
except ImportError as e:
    print e
    print
    print "This program was called as       :", sys.argv[0]
    print "Absolute path to this program is :", abs_path
    print "Base path is                     :", base_path
    sys.exit(1)

if not len(sys.argv) in [2, 3]:
    print "Usage: order_triangulation.py IN_TRIANGULATION_FILE OUT_TRIANGULATION_FILE"
    print "       order_triangulation.py IN_TRIANGULATION"
    sys.exit(1)

in_file = sys.argv[1]

if len(sys.argv) == 2:
    if in_file[-5:].lower() == '.trig':
        out_file_base = in_file[:-5]
    else:
        out_file_base = in_file
    out_file = out_file_base + '_ordered.trig'

    out_file_un = out_file_base + '_unorderable.trig'

else:
    out_file = sys.argv[2]

    out_file_un = out_file


t=triangulation(open(in_file,'r').read())
t.orient()
if t.is_ordered():
    print "Triangulation already ordered"
else:
    orderings=t.find_orderings()
    print "Number of orderings of the triangulation:", len(orderings)
    if not orderings:
        print "No ordering could be found"
        print "Output File:", out_file_un
        open(out_file_un,'w').write(t.to_SnapPea())

        sys.exit(1)
    good_ordering=orderings[0]
    for ordering in orderings:
        if t.is_edge_order_orientation_preserving(ordering):
            print "Found orientation preserving ordering"
            break

    t.check_consistency()
    t.reorder_tets(good_ordering)
    t.check_consistency()
    assert t.is_ordered()

print "Output File:", out_file
open(out_file,'w').write(t.to_SnapPea())
    
