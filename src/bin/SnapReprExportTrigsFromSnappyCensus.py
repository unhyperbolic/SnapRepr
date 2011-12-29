#!/usr/bin/python

import optparse
import sys
import os

import snappy

maxNumTetsSeparateDirs = 30

def cusped_directory_name(trig):
    num_tets = trig.num_tetrahedra()
    if num_tets >= maxNumTetsSeparateDirs:
        return "cusped_%d_higher" % maxNumTetsSeparateDirs
    else:
        return "cusped_%d" % num_tets

def createCuspedDirectory(trig):
    name = cusped_directory_name(trig)

    if os.path.exists(name):
        assert os.path.isdir(name)
    else:
        os.mkdir(name)

    return name

def filename_generic(trig):
    return trig.name()

def filename_morwen_link_factory(no_components, no_crossings):
    basename = "%d" % no_crossings
    if no_components > 1:
        basename += "^%d" % no_components
    def filename_morwen_link(trig, basename = basename):
        return basename + "_" + trig.name()
    return filename_morwen_link

def find_trig_with_least_tets(trig, no_tries):
    m=[trig] + [snappy.Manifold(trig.name()) for x in range(no_tries-2)]
    m.sort(key = lambda x:x.num_tetrahedra())
    return m[0]

def write_files_for_sequence_of_trigs(trigs, 
                                      filename_function = filename_generic,
                                      no_tries = 1, max_tets = 1000):
    for trig in trigs:
        n_trig = find_trig_with_least_tets(trig, no_tries)
        print n_trig.name()
        if n_trig.num_tetrahedra() <= max_tets:
            filename = (
                createCuspedDirectory(trig) + 
                '/' +
                filename_function(trig) + 
                ".trig")
            print filename

            print "saved"
            n_trig.save(filename)

def print_usage():
    print "Usage: export_trigs_from_SnapPy_census.py OrientableCuspedCensus MIN_NO_TETS MAX_NO_TETS"
    print "       export_trigs_from_SnapPy_census.py LinkExteriors RETRIANGULATION_TRIES MIN_CROSSINGS MAX_CROSSINGS MIN_COMPONENT MAX_COMPONENT MAX_NO_TETS"
    print "       export_trigs_from_SnapPy_census.py MorwenLinks   RETRIANGULATION_TRIES MIN_CROSSINGS MAX_CROSSINGS MIN_COMPONENT MAX_COMPONENT MAX_NO_TETS"
    print "       export_trigs_from_SnapPy_census.py KnotExteriors RETRIANGULATION_TRIES MIN_CROSSINGS MAX_CROSSINGS MAX_NO_TETS"

if not len(sys.argv) >= 3: 
    print_usage()
    sys.exit(1)

try:
    args = [int(x) for x in sys.argv[2:]]
except:
    print_usage()
    sys.exit(1)    

if sys.argv[1] == "OrientableCuspedCensus":
    if not len(args) == 2:
        print_usage()
        sys.exit(1)

    min_no_tets, max_no_tets = args
    if min_no_tets > max_no_tets:
        print "MIN number needs to be larger than MAX"
        sys.exit(1)
    if max_no_tets > 8:
        print "CuspedOrientableCensus only goes up to 8 tetrahedra"
        sys.exit(1)
    for i in range(min_no_tets, min(5, max_no_tets)+1):
        write_files_for_sequence_of_trigs(
            [M for M in snappy.OrientableCuspedCensus()[0:301] if M.num_tetrahedra() == i])
    if min_no_tets <= 6 and 6 <= max_no_tets:
        write_files_for_sequence_of_trigs(
            snappy.OrientableCuspedCensus()[301:1263])
    if min_no_tets <= 7 and 7 <= max_no_tets:
        write_files_for_sequence_of_trigs(
            snappy.OrientableCuspedCensus()[1263:4815])
    if min_no_tets <= 8 and 8 <= max_no_tets:
        write_files_for_sequence_of_trigs(
            snappy.OrientableCuspedCensus()[4815:])
        
elif sys.argv[1] == "LinkExteriors":
    if not len(args) == 6:
        print_usage()
        sys.exit(1)
    retriangulation_tries, min_crossings, max_crossings, min_components, max_components, max_tets = args

    if max_components > 5:
        print "LinkExteriors takes at most 5 components"
        sys.exit(1)
    
    for components in range(min_components, max_components + 1):
        def no_crossings(trig):
            return int(trig.name().split('_')[0].split('^')[0])

        write_files_for_sequence_of_trigs(
            [trig for trig in snappy.LinkExteriors(components) if
             no_crossings(trig) >= min_crossings and 
             no_crossings(trig) <= max_crossings],
            no_tries = retriangulation_tries,
            max_tets = max_tets)

elif sys.argv[1] == "MorwenLinks":
    if not len(args) == 6:
        print_usage()
        sys.exit(1)

    retriangulation_tries, min_crossings, max_crossings, min_components, max_components, max_tets = args

    if max_crossings > 14:
        print "MorwenLinks only go up to 14 crossings"
        sys.exit(1)
        
    
    for no_components in range(min_components, max_components + 1):
        for no_crossings in range(min_crossings, max_crossings + 1):
            write_files_for_sequence_of_trigs(
                snappy.MorwenLinks(no_components, no_crossings),
                filename_function = filename_morwen_link_factory(
                    no_components, no_crossings),
                no_tries = retriangulation_tries,
                max_tets = max_tets)

elif sys.argv[1] == "KnotExteriors":
    if not len(args) == 4:
        print_usage()
        sys.exit(1)

    retriangulation_tries, min_crossings, max_crossings, max_tets = args

    if min_crossings < 3:
        print "At least 3 crossings"
        sys.exit(1)
    if max_crossings > 16:
        print "At most 16 crossings"
        sys.exit(1)

    table_crossings_alternating = [
             0,  #  0 crossings
             0,  #  1 crossing
             0,  #  2 crossings
             0,  #  3 crossings
             1,  #  4 crossings
             2,  #  5 crossings
             4,  #  6 crossings
             7,  #  7 crossings
            14,  #  8 crossings
            32,  #  9 crossings
            73,  # 10 crossings
           196,  # 11 crossings
           563,  # 12 crossings
          1851,  # 13 crossings
          6729,  # 14 crossings
         26265,  # 15 crossings
        111528,  # 16 crossings
        491327   # 17 crossings
        ]
    
    table_crossings_nonalternating = [
              0,  #  0 crossings
              0,  #  1 crossing
              0,  #  2 crossings
              0,  #  3 crossings
              0,  #  4 crossings
              0,  #  5 crossings
              0,  #  6 crossings
              0,  #  7 crossings
              0,  #  8 crossings
              3,  #  9 crossings
             11,  # 10 crossings
             53,  # 11 crossings
            238,  # 12 crossings
           1126,  # 13 crossings
           6236,  # 14 crossings
          33672,  # 15 crossings
         201702,  # 16 crossings
        1210608   # 17 crossings
        ]
    


    write_files_for_sequence_of_trigs(
            snappy.AlternatingKnotExteriors()[
                table_crossings_alternating[min_crossings]:table_crossings_alternating[max_crossings+1]],
                no_tries = retriangulation_tries,
                max_tets = max_tets)

    write_files_for_sequence_of_trigs(
            snappy.NonalternatingKnotExteriors()[
                table_crossings_nonalternating[min_crossings]:table_crossings_nonalternating[max_crossings+1]],
                no_tries = retriangulation_tries,
                max_tets = max_tets)

else:
    print_usage()
    sys.exit(1)

