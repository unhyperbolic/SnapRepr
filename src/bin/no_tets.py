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


if not len(sys.argv)==2:
    print "Usage: no_tets.py IN_TRIANGULATION_FILE"
    sys.exit(1)

t=triangulation(open(sys.argv[1],'r').read())
print t.num_tets
