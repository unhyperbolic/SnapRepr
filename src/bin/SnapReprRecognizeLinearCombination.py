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
    from linearCombinations import filterRepresentativeMfds, twoTerms
    from linearCombinations.formatLinearCombination import formatLinearCombination
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

    if not args:
        parser.print_help()
        sys.exit(1)

    print "Setting maximal errors..."

    globalsettings.setSetting("maximalError",
                              mpmath.mpf("0.1") ** options.maximalError)
    globalsettings.setSetting("maximalErrorDigits",
                              options.maximalError)

    mpmath.mp.dps = options.maximalError + options.extraDigits

    if options.tableFile:
        censusTableFile = options.tableFile
    else:
        censusTableFile = base_path + "/tests/CensusTable/exampleCensusTable.csv"

    print "Reading census table %s..." % censusTableFile
    
    censusTable = readCensusTable(censusTableFile, sortKey = "Volume")

    print "Finding representatives..."
    representatives = filterRepresentativeMfds.filterRepresentativeMfds(censusTable.listOfDicts)

    print "Initializing C structures..."
    table = twoTerms.censusTable(representatives)
    
    print "Processing files..."
    processCsvFiles(args, table)

def processCsvFiles(files, table):

    for csvFileName in files:
        
        if not '_linComb' in csvFileName:

            if csvFileName[-4:] == '.csv':
                csvOutFilename = csvFileName[:-4] + '_linComb.csv'
            else:
                csvOutFilename = csvFileName + '_linComb.csv'

            print csvFileName, "->", csvOutFilename

            csvTable = csvUtilities.readCensusTable.readCensusTable(csvFileName, 
                                                                    readHeaderFromFile = False)
  
            csv_writer = csv.DictWriter(open(csvOutFilename,'w'),
                                        fieldnames = csvTable.header)

            for row in csvTable.listOfDicts:
                if (row.has_key('Volume') and 
                    not row['Volume'] is None and
                    row['Volume'] > globalsettings.getSetting("maximalError")):
                    row['LinearCombinations'] = ' | '.join(
                        [
                            formatLinearCombination(l) 
                            for l in 
                            table.findAsUpToTwoTerms(row['Volume'])])
                
                csv_writer.writerow(row)
                    
    
def create_parser():
    
    parser = optparse.OptionParser(
        usage = ("Usage: %prog [options] CSV_FILES"))

    parser.add_option("-p", "--precision",
                      type = "int",
                      dest = "maximalError", default = 50,
                      help = "Maximal error in decimal digits")
    parser.add_option("-P", "--extra-digits",
                      type = "int",
                      dest = "extraDigits", default = 10,
                      help = "Digits of precision in intermediate calculation")

    parser.add_option("-t", "--table",
                      type = "str",
                      dest = "tableFile", default = None,
                      help = "Table of census manifold volumes used")

    return parser

main()
