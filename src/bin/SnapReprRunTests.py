#!/usr/bin/python
import os
import sys
import re
import optparse
import traceback
import StringIO
import csv
import doctest

this_path, this_file = os.path.split(sys.argv[0])
abs_path = os.path.abspath(this_path)
base_path, this_dir = os.path.split(abs_path)
sys.path.append(base_path)

try:
    import manifold.triangulation
    import manifold.slN
    import algebra.polynomial
    import algebra.magma
    import globalsettings
except ImportError as e:
    print e
    print
    print "This program was called as       :", sys.argv[0]
    print "Absolute path to this program is :", abs_path
    print "Base path is                     :", base_path
    sys.exit(1)

globalsettings.setSetting("testpath", base_path + "tests/")

doctest.testmod(manifold.triangulation)
#doctest.testmod(manifold.slN)
doctest.testmod(algebra.polynomial)
#doctest.testmod(algebra.magma)
