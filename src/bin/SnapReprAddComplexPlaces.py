#!/usr/bin/python
import os
import sys
import re
import optparse
import traceback
import StringIO
import csv

from fractions import Fraction

import mpmath

this_path, this_file = os.path.split(sys.argv[0])
abs_path = os.path.abspath(this_path)
base_path, this_dir = os.path.split(abs_path)
sys.path.append(base_path)

try:
    from csvUtilities.readCensusTable import readCensusTable
    import globalsettings
    import csvUtilities
    import csvUtilities.keyTable
    from linearCombinations import filterRepresentativeMfds, twoTerms
    from linearCombinations.formatLinearCombination import formatLinearCombination
    from utilities.basicAlgorithms import safeDictLookup
    from algebra import pari
except ImportError as e:
    print e
    print
    print "This program was called as       :", sys.argv[0]
    print "Absolute path to this program is :", abs_path
    print "Base path is                     :", base_path
    sys.exit(1)

def main():
    parser = create_parser()
    # Parse command line options

    options, args = parser.parse_args()

    if not len(args) == 1:
        parser.print_help()
        sys.exit(1)

    filename = args[0]

    if filename[-4:] == '.csv':
        outFilenameBase = filename[:-4]
    else:
        outFilenameBase = filename
        
    outFilename = outFilenameBase + '_complexPlaces.csv'

    censusTable = readCensusTable(filename, convertData = False)

    fieldnames = censusTable.header
    if not "NumberOfComplexPlaces" in fieldnames:
        fieldnames = fieldnames + ["NumberOfComplexPlaces"]

    csv_writer = csv.DictWriter(open(outFilename,'w'),
                                fieldnames = fieldnames)

    csv_writer.writerow(dict(zip(fieldnames,fieldnames)))

    for row in censusTable.listOfDicts:
        if row.has_key('InvariantTraceField'):
            if (row['InvariantTraceField'] and 
                not row['InvariantTraceField'] == '-'):

                pariStr = "nfinit(%s).r2" % row['InvariantTraceField']

                row["NumberOfComplexPlaces"] = pari.pari_eval(pariStr)

        csv_writer.writerow(row)

def create_parser():
    
    parser = optparse.OptionParser(
        usage = ("Usage: %prog [options] CSV_FILE"))

    return parser

main()
