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

inputFileName = sys.argv[1]

if len(sys.argv) == 2:
    # using the program without a OUT_TRIANGULATION_FILE argument

    # chop of the extenstion
    if inputFileName[-5:].lower() == '.trig':
        outFileNameBase = inputFileName[:-5]
    else:
        outFileNameBase = inputFileName
        
    # add "_ordered" or "_unordered" to filename
    outFileNameOrdered = outFileNameBase + '_ordered.trig'
    outFileNameUnordered = outFileNameBase + '_unorderable.trig'

else:
    # using the program with a OUT_TRIANGULATION_FILE argument
    outFileNameOrdered   = sys.argv[2]
    outFileNameUnordered = sys.argv[2]

t = triangulation(open(inputFileName,'r').read())
t.orient()

if t.is_ordered():
    # triangulation is already ordered

    print "Triangulation already ordered"
    outFileName = outFileNameOrdered
    orderings = [ True ]
else:
    # triangulation is not ordered yet
    
    # try to find a canonical ordering
 
    orderings = t.find_orderings()
    print "Number of orderings of the triangulation:", len(orderings)

    if not orderings:
        # No ordering found

        print "No ordering could be found"
        outFileName = outFileNameUnordered

    else:

        # Canonical ordering is the first ordering such that
        # all tetrahedra are oriented the same, if there is no such
        # ordering, we take the first ordering

        canonicalOrdering = orderings[0]

        for ordering in orderings:
            if t.is_edge_order_orientation_preserving(ordering):
                print "Found orientation preserving ordering"
                canonicalOrdering = ordering
                break

        # Apply canonical ordering

        t.check_consistency()
        t.reorder_tets(canonicalOrdering)
        t.check_consistency()
        assert t.is_ordered()

print "Output File:", out_file
open(out_file,'w').write(t.to_SnapPea())
    
if not orderings:
    sys.exit(1)
