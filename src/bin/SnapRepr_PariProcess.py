#!/usr/bin/python
import os
import sys
import re
import optparse
import traceback
import StringIO
import csv

from fractions import Fraction

this_path, this_file = os.path.split(sys.argv[0])
abs_path = os.path.abspath(this_path)
base_path, this_dir = os.path.split(abs_path)
sys.path.append(base_path)

try:
    from algebra import pari
except ImportError as e:
    print e
    print
    print "This program was called as       :", sys.argv[0]
    print "Absolute path to this program is :", abs_path
    print "Base path is                     :", base_path
    sys.exit(1)

def main():
    
    while True:
        inputStr = sys.stdin.readline()

        sys.stderr.write("got a line")

        if not inputStr:
            break
        inputStr = inputStr.strip()
        inputRaw = eval(inputStr)
        outputRaw = pari.pari_eval(inputRaw)
        print repr(outputRaw)

        sys.stdout.flush()

        sys.stderr.write("printing a line")

main()
