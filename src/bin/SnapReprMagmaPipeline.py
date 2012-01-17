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
    from manifold.triangulation import triangulation, read_triangulation_from_file
    from algebra.polynomial import Polynomial
    from algebra.solvePolynomialEquations import solvePolynomialEquations
    import algebra.mpmathFunctions
    import manifold.slN
    import algebra.magma
    import hashlib
    from utilities import basicAlgorithms
    import globalsettings
    from algebra.exceptions import NumericalError
    import utilities.printNumbers
    from csvUtilities import readCensusTable
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

    globalsettings.setSetting("maximalError",
                              mpmath.mpf("0.1") ** options.maximalError)
    globalsettings.setSetting("maximalErrorDigits",
                              options.maximalError)

    mpmath.mp.dps = options.maximalError + options.extraDigits

    # Exit writing the header for the CSV file
    if options.print_csv_header:
        output = StringIO.StringIO()
        csv_writer = csv.DictWriter(output, fieldnames = readCensusTable.header)
        csv_writer.writerow(dict(zip(readCensusTable.header,
                                     readCensusTable.header)))
        print output.getvalue()[:-1]
        sys.exit(0)

    # Exit with usage information
    if not len(args)==1:
        parser.print_help()
        sys.exit(1)

    # Read input file
    input_file = open(args[0], 'r').read()
    first_line = ''.join(input_file.split('\n')[0:1])

    # Branch depending on whether input file is a triangulation

    if "% triangulation" in first_line.lower():
        # Triangulation
        write_magma_files(options, triangulation_filename = args[0])

    elif "MAGMA=INPUT=FROM=SNAPREPR" in input_file:

        print "This file is supposed to be processed with MAGMA"
        print "Run:"
        print "      magma_out %s" % args[0]
        print "Then:"
        print "      slN_representations.py %s_out" % args[0]

    elif "MAGMA=OUTPUT=FOR=SNAPREPR" in input_file:
        # magma breaks long lines putting a \ at the break
        # we replace those line breaks
        magma_out = ''.join([x.strip() for x in input_file.split('\\\n')])

        if "MAGMA=OUTPUT=STAGE2=FOR=SNAPREPR" in input_file:

            process_magma_output_headers_stage2(options, magma_out, 
                                                magma_filename = args[0])
            
        else:
            process_magma_output_headers(options, magma_out,
                                         magma_filename = args[0])
    else:

        sys.stderr.write("Could not recognize a SnapPea triangulation file or a MAGMA output file ????\n\n")
        
        # XXX:needs to be fixed!

        print '"-","%s","-","-","-","-","-","-","-","-","-","-","-","FILE NOT RECOGNIZED","-","-","-"' % args[0]

def write_magma_files(options, triangulation_filename):
    # The input file is a triangulation
    # produce the Magma file

    # get N for SL(N,C) representations
    N = options.N

    # parse and orient the triangulation
    t = read_triangulation_from_file(triangulation_filename)
    t.orient()

    # determine the base for the output filenames
    # If this is not overwritten by the command line option,
    # this strips the ".trig" extension of the triangulation file
    # and also tries to write the file into ../magma_slN/
    outfile_base = options.magma_base
    if not outfile_base:
        outfile_base = get_outfile_base(triangulation_filename)
        print "Writing files to %s..." % outfile_base

    # if N is even, compute all cohomology obstruction classes
    # the cohomology obstruction is an element in H^2(M, partial M; Z/2)
    if N % 2 == 0:
        H = manifold.slN.get_all_obstruction_classes(t)
    else:
        H = [None] 
    
    # For each cohomology obstruction class h with index c
    for c, h in basicAlgorithms.indexedIterable(H):
        # List all Ptolemy relations with respect to the cohomology class h
        pt_eqns = manifold.slN.get_Ptolemy_relations(t, N, h)

        # List additional relations for fixing decoration
        #fix_decoration_method = "MANUAL_FIX_DECORATION"
        #addl_eqns = manifold.slN.get_additional_eqns_manual(t, N)

        # List additional relations for fixing decoration
        fix_decoration_method = "AUTOMATIC_FIX_DECORATION"
        addl_eqns = manifold.slN.get_additional_eqns_independent(t, N)

        # Add all equations together
        pre_eqns = pt_eqns + addl_eqns

        # Find which Ptolemy coordinates are getting identified
        # id_c_parms is of type equivalence_with_sign which encapsulates
        # an equivalence relationship with signs
        id_c_parms = manifold.slN.get_identified_c_parameters(t, N)

        # identify the Ptolemy coordinates
        eqns       = manifold.slN.identify_c_parameters(pre_eqns, id_c_parms)

        # Append an equation of the form "x * y * t -1" with dummy variable t
        # This is to prevent Ptolemy coordinates from being zero
        eqns.append( manifold.slN.polynomialNonZeroCondition(eqns,'t') )

        # Make the input for magma
        # pre_vars tells the procedure to list t in the term order used
        # for the computation of the Groebner basis
        term_order = get_term_order(eqns, pre_vars = ['t'])
        out = algebra.magma.primary_decomposition(eqns, term_order = term_order)

        # Compute the hash of the ideal
        hash_eqns = hash_ideal(eqns, term_order) 

        # Human readable comment at the beginning of the magma file
        comment = generate_magma_comment(t, N, triangulation_filename, h, id_c_parms, pre_eqns)

        # Header for the magma file consists of a bunch of PRINT statements
        # The results are parsed from the magma output file to recover
        # the triangulation, N, the cohomology obstruction...
        # when processing the output magma produced.
        header = generate_magma_header(t, N, triangulation_filename, c, 
                                       outfile_base, fix_decoration_method,
                                       term_order, hash_eqns)

        # name of magma input file contains 
        # N and the index of the cohomology obstruction
        # outfile = outfile_base + "_sl%d_c%d_hash%s.magma" % (N,c,hash_eqns)
        outfile = outfile_base + "_sl%d_c%d.magma" % (N,c)
        print "          ", outfile
        open(outfile,'w').write(comment + header + out)

def process_magma_output_headers(options, magma_out, magma_filename):

    outputFileNameBase = magma_filename
    if magma_filename[-6:] == '.magma':
        outputFileNameBase = magma_filename[:-6]
    if magma_filename[-10:] == '.magma_out':
        outputFileNameBase = magma_filename[:-10]
    outputFileName = outputFileNameBase + '.magmi'
    
    print "Writing next file to be processed with magma:", outputFileName

    error_condition = None

    try:
        # parse "TERM ORDER       : [ 'a', 'b', 'c' ]
        term_order  = re.search("TERM ORDER       : \[(.*?)\]",
                                magma_out,
                                re.DOTALL).group(1)
        term_order  = [ x.strip()[1:-1] for x in term_order.split(',')]
    except:
        error_condition = "Could not parse TERM ORDER in magma output file"

    prim_decomp = []

    try:
        prim_decomp = algebra.magma.parse_primary_decomposition(magma_out)
    except:
        error_condition = "Could not parse PRIMARY DECOMPOSITION"

    header  = '/* MAGMA=INPUT=FROM=SNAPREPR */\n'
    header += '/* STAGE2=INPUT=FOR=ABSOLUTIZE */\n'
    header += '/* This will print a copy of the old file magma file */\n'
    header += 'print "%s\\n";\n' % '\\n'.join(magma_out.split('\n'))
    header += '/* This will print the new stuff */\n'
    header += 'print "MAGMA=OUTPUT=STAGE2" cat "=FOR=SNAPREPR";\n'
    
    if not error_condition:
        pass
    else:
        print error_condition
        header += 'print "%s";\n' % error_condition

    newMagmaCommands = ""
        
    for index, component in basicAlgorithms.indexedIterable(prim_decomp):

        newMagmaCommands += 'print "VARIETY=NUMBER=%d=BEGINS=HERE";\n' % index

        newMagmaCommands += (
            'print "DIMENSION        : %d";\n' % component.dimension)

        if component.dimension == 0:
            newMagmaCommands += (
                'print "NUMBEROFPOINTS   : %d";\n' % component.numberOfPoints)

            newMagmaCommands += algebra.magma.absolutizeVariety(component,
                                                                term_order)


        newMagmaCommands += 'print "VARIETY=NUMBER=%d=ENDS=HERE";\n' % index

    outputFile = open(outputFileName, 'w')
    outputFile.write(header + newMagmaCommands)

    
    # print "FIRST IDEAL:"
    # print "NUMBER FIELD:", result of computiation
    # print "

def process_magma_output_headers_stage2(options, magma_out, magma_filename):
    error_condition = ""
    output = StringIO.StringIO()

    csvOutputFilenameBase = magma_filename

    if magma_filename[-6:] == '.magmi':
        csvOutputFilenameBase = magma_filename[:-6]
    if magma_filename[-10:] == '.magmi_out':
        csvOutputFilenameBase = magma_filename[:-10]
    csvOutputFileName = csvOutputFilenameBase + '_magmi.csv'

    sys.stderr.write("Output File: %s\n" % csvOutputFileName)

    csv_writer = csv.DictWriter(output, 
                                fieldnames = readCensusTable.header,
                                restval = "")
    csv_dict = { "File": magma_filename}

    # parse the header we produced to recover N, cohomology obstruction class,
    # name of the triangulation, triangulation data, cputime
    
    try:
        c      =   int(re.search("IND OF COH CLASS : (.*)",magma_out).group(1))
        csv_dict['Obstruction Class'] = c
    except:
        error_condition = "Could not parse Cohomology class given in magma output file"

    try:
        cputime = float(re.search("CPUTIME          : (.*)",magma_out).group(1))
        csv_dict['CPUTIME'] = cputime
    except:
        error_condition = "Could not parse CPUTIME in magma output file"
        
    try:
        name   =       re.search("NAME             : (.*)",magma_out).group(1)
        csv_dict['Manifold'] = name
    except:
        error_condition = "Could not parse NAME in magma output file"

    try:
        N      =   int(re.search("N                : (.*)",magma_out).group(1))
        csv_dict['SL(N,C)'] = N
    except:
        error_condition = "Could not parse N in magma output file"

    try:
        t_data =       re.search(('==TRIANGULATION=BEGINS=='
                                  '(.*?)'
                                  '==TRIANGULATION=ENDS=='),
                                 magma_out,
                                 re.DOTALL).group(1)
        t = triangulation(t_data.strip())
        csv_dict["Tetrahedra"] = t.num_tets
        csv_dict["Cusps"] =  t.num_or_cusps
        t.orient()
        if t.is_ordered():
            csv_dict["Ordered"] = "Ordered"
        else:
            csv_dict["Ordered"] = "Unordered"            
    except:
        error_condition = "Could not parse TRIANGULATION in magma output file"

    # Parse primary decomposition from the magma file.
    # This returns a list of objects of type prime_ideal.
    # A PrimeIdeal object is derived from a list, but also
    # has a dimension and numberOfPoints field.
    try:
        prim_decomp = algebra.magma.parse_primary_decomposition(magma_out)
        csv_dict["Number Components"] = len(prim_decomp)
    except Exception as e:
        print "Here..."
        print e
        error_condition = "Could not parse primary decomposition"
        
    if not "PRIMARY=DECOMPOSITION=ENDS=HERE" in magma_out:
        error_condition = 'MAGMA crashed/aborted'
        
    if "All virtual memory has been exhausted" in magma_out:
        error_condition = 'MAGMA out of memory'            

    # Handle various Error Conditions
    if error_condition:
        csv_dict["Warning"] = error_condition
        sys.stderr.write(error_condition+"\n\n")
        csv_writer.writerow(csv_dict)

    # Empty ideal found
    elif not len(prim_decomp):
        csv_dict["Warning"] = "Empty Ideal"
        sys.stderr.write(error_condition+"\n\n")
        csv_writer.writerow(csv_dict)

    else:
        # For every prime ideal prime_ideal in the primary decomposition
        # index_prime_ideal is the index of the prime_ideal in the list
        for index_prime_ideal, prime_ideal in (
            basicAlgorithms.indexedIterable(prim_decomp)):

            csv_dict['Index Component'] = index_prime_ideal
            csv_dict["Dimension"] = prime_ideal.dimension
            if not prime_ideal.numberOfPoints is None:
                csv_dict['Number Solutions'] = prime_ideal.numberOfPoints

            try:
                varietyData = re.search(
                    ('VARIETY=NUMBER=%d=BEGINS=HERE'
                     '(.*?)'
                     'VARIETY=NUMBER=%d=ENDS=HERE') % (index_prime_ideal, 
                                                       index_prime_ideal),
                    magma_out,
                    re.DOTALL).group(1).strip()
            except:
                error_condition = "Could not parse Variety Output %d" % index_prime_ideal

            if prime_ideal.dimension == 0:
                try:
                    csv_dict["PtolemyField"] = (
                        algebra.magma.parseDefiningRelationshipAbsolutizedAlgebraicClosure(
                            varietyData))
                except:
                    error_condition = "Could not parse Field in Variety %d" % index_prime_ideal

                try:
                    cvols = get_complex_volumes(t, N, c, prime_ideal, not_paranoid = not options.paranoid)

                    for i, cvol in basicAlgorithms.indexedIterable(cvols):
                        csv_dict["Index Solution"] = i
                        csv_dict["Volume"] = utilities.printNumbers.printRealNumberAsFixed(cvol.real)
                        csv_dict["CS"] = utilities.printNumbers.printRealNumberAsFixed(cvol.imag)
                        csv_writer.writerow(csv_dict)

                except NumericalError as n:
                    csv_dict["Warning"] = "Numerical Error"
                    csv_dict["Volume"] = "-"
                    csv_dict["CS"] = "-"
                    csv_writer.writerow(csv_dict)
                    sys.stderr.write(traceback.format_exc())
                    sys.stderr.write("\n\n\n")
                except:
                    csv_dict["Warning"] = "UNKNOWN ERROR"
                    csv_dict["Volume"] = "-"
                    csv_dict["CS"] = "-"
                    csv_writer.writerow(csv_dict)
                    sys.stderr.write(traceback.format_exc())
                    sys.stderr.write("\n\n\n")

    resultString = output.getvalue()
    print resultString[:-1]
    open(csvOutputFileName, 'w').write(resultString)
    

def get_complex_volumes(t, N, c, primeIdeal, not_paranoid = False):

    all_cvols = []

    # Find the points in the variety
    # solvePolynomialEquations assumes that prime_ideal is given
    # in Groebner basis
    sys.stderr.write("Solving...\n")

    def conversionFunction(c):
        if isinstance(c, Fraction):
            return mpmath.mpc(c.numerator) / mpmath.mpc(c.denominator)
        else:
            return mpmath.mpc(c)

    primeIdealFloatCoeff = primeIdeal.convertCoefficients(conversionFunction)
    solutions = solvePolynomialEquations(
        primeIdealFloatCoeff,
        polynomialSolver = algebra.mpmathFunctions.PolynomialSolver)
    sys.stderr.write("Solved...\n")

    # Solutions is a list with one element per Galois conjugate
    if not len(solutions) == primeIdeal.numberOfPoints:
        raise NumericalError("Number of solutions doesn't match")

    id_c_parms = manifold.slN.get_identified_c_parameters(t, N)
    if N % 2 == 0:
        H = manifold.slN.get_all_obstruction_classes(t)
        h = H[c]
    else:
        h = None

    # Each element in solutions is a dictionary
    # mapping the variable name to the numerical value
    for solution in solutions:
        # the Ptolemy coordinates are projective
        # We multiply all Ptolemy coordinates by the same random number
        # to avoid taking the logarithm of a negative real number
        solution = manifold.slN.multiply_solution_by_constant(
            solution, 
            mpmath.mpc("1.0","0.32433234644"))

        # We assign a value to each Ptolemy coordinate based on the 
        # solution
        c_parms = manifold.slN.map_solution_to_c_parameters(solution,
                                                            id_c_parms)
        # Consistency check
        if not not_paranoid:
            manifold.slN.check_solution_identification(t, N, c_parms)

        # Consistency check
        manifold.slN.check_solution_on_gluing_equations(
            t, N, c_parms, h)

        # Get Ptolemy cochain
        # P consists of Ptolemy_cochain objects which encapsulate
        #   a sign for orientation and 6 edge parameters c
        P = manifold.slN.get_Ptolemy_cochain(t, N, c_parms, no_check = not_paranoid)

        # convert the 6 edge parameters into an element of the
        # extended Bloch group represented by (w0,w1,w2)
        ws = [manifold.bloch_group.w_triple(x, no_check = not_paranoid) for x in P]

        # convert the triples (w0,w1,w2) into [z;p,q]
        # a different representation of an extended Bloch group element
        zpqs = [manifold.bloch_group.zpq_triple(x, no_check = not_paranoid) for x in ws]

        # Compute the volume function and sum for all [z;p,q]
        vols = [x.volume() for x in zpqs]
        vol = sum(vols)

        # Compute the L function and sum for all [z;p,q]
        cvols = [x.L_function() for x in zpqs]
        cvol = sum(cvols) / 1j

        sys.stderr.write("   %s\n" % (
            utilities.printNumbers.printComplexNumberAsFixed(cvol)))

        # Consistency check
        maxErr = globalsettings.getSetting("maximalError")

        if not abs(vol - cvol.real) < maxErr:
            raise NumericalError(val = [vol, cvol], 
                                 msg = "Vol and Complex Vol not matching")
        all_cvols.append(cvol)
    return all_cvols

#            sys.stderr.write(
#                "Irreducible component of variety No %d is of dimension %d\n" 
#                % (index_prime_ideal, prime_ideal.dimension))
#            print csv_output + '"-","","-","-",%.1f' % cputime
#
#            if options.write_nonzero_dim_ideals:
#                this_outfile = (
#                    outfile_base +
#                    "_sl%d_%d_component%d.magma" % (N,c,no_irreducible_comp))
#                out = magma.ideal_to_magma.genus(prime_ideal)

#                try:
#                    open(this_outfile,'w').write(out)
#                except:
#                    open(path + '/' + this_outfile,'w').write(out)

def get_term_order(polys, pre_vars = [], post_vars = []):
    all_vars = sum([p.variables() for p in polys], [])
    sort_vars = set(all_vars) - set(pre_vars) - set(post_vars)
    sort_vars = list(sort_vars)
    sort_vars.sort()

    return pre_vars + sort_vars + post_vars


def create_parser():
    parser = optparse.OptionParser(
        usage = ("Usage: %prog [options] IN_TRIANGULATION\n" + 
                 "       %prog [options] MAGMA_OUT_FILE"))
    
    parser.add_option("-C", "--paranoid",
                      dest = "paranoid", default = False,
                      action = "store_true",
                      help = "do extensive checking")
    parser.add_option("-b", "--bloch",
                      dest = "display_bloch", default = False,
                      action = "store_true",
                      help = "display elements in the extended Blochgroup")
    parser.add_option("-d", "--d",
                      dest = "display_c", default = False,
                      action = "store_true",
                      help = "display c parameters")
    parser.add_option("-n", "--n",
                      dest = "N", default = 2,
                      type = "int",
                      help = "N of SL(N,C) representation, default 2")
    parser.add_option("-m", "--magma-base",
                      dest = "magma_base", default = None,
                      help = "base filename for MAGMA files")
    parser.add_option("-c", "--csv",
                      dest = "as_csv", default = False,
                      action = "store_true",
                      help = "format results as csv file")
    parser.add_option("-H", "--header-csv",
                      dest = "print_csv_header", default = False,
                      action = "store_true",
                      help = "print csv header")
    parser.add_option("-i", "--ideals",
                      dest = "write_nonzero_dim_ideals", default = False,
                      action = "store_true",
                      help = "Write magma files for non zero dimension ideals")
    parser.add_option("-p", "--precision",
                      type = "int",
                      dest = "maximalError", default = 60,
                      help = "Maximal error in decimal digits")
    parser.add_option("-P", "--extra-digits",
                      type = "int",
                      dest = "extraDigits", default = 20,
                      help = "Digits of precision in intermediate calculation")
    parser.add_option("-e", "--exact",
                      dest = "exact", default = False,
                      action = "store_true",
                      help = "Write exact solutions")
    parser.add_option("-s", "--solution",
                      dest = "print_solution",
                      default = False, action = "store_true",
                      help = "print solutions")
    return parser

def generate_magma_comment(t,N,triangulation_filename,h,id_c_parms, pre_eqns):
    comment  = "/******************************** \n"
    comment += "MAGMA gluing equations for PSL(%d,C) representation\n" % N
    comment += "\n"
    comment += "File: %s \n\n" % triangulation_filename
    if not t.is_ordered():
        comment += "***** Triangulation was not ordered *****\n\n"
    for tet in t.tet_list:
        if tet.positive_orientation:
            comment += "Tetrahedron %d was oriented positively.\n" % tet.index
        else:
            comment += "Tetrahedron %d was oriented negatively.\n" % tet.index    
    comment += "\n\n"
    comment += "H^2(M,boundary M; Z/2) is an obstruction to lifting the\n"
    comment += "representation to SL(2,C).s For this set of equations the\n"
    comment += "following cohomology class was chosen:\n"
    
    comment += "Cohomology class (representative in C^2): %s\n\n" % h
    comment += "Basis for C^2: %s\n\n" % t.get_face_classes()
    comment += "Identified c parameters:\n%s\n\n" % id_c_parms.pretty_print()
    comment += "Equations:\n    %s\n\n" % "\n    ".join(map(str,pre_eqns))
    comment += "********************************/ \n\n\n"
    
    return comment

def generate_magma_header(t, N, triangulation_filename, c, outfile_base,
                          fix_decoration_method, term_order, hash_eqns):
    header  = 'print "MAGMA=OUTPUT" cat "=FOR=SNAPREPR";\n'
    header += '/* MAGMA=INPUT=FROM=SNAPREPR */\n'
    header += 'print "MAGMA BASE       : %s";\n' % outfile_base
    header += 'print "IND OF COH CLASS : %s";\n' % c
    header += 'print "TRIANGULATION    : %s";\n' % (os.path.abspath(triangulation_filename))
    header += 'print "NAME             : %s";\n' % t.name
    header += 'print "N                : %d";\n' % N
    header += 'print "DECORATIONS      : %s";\n' % fix_decoration_method
    header += 'print "TERM ORDER       : %s";\n' % term_order
    header += 'print "HASH OF IDEAL    : %s";\n' % hash_eqns
    header += 'print "==TRIANGULATION=BEGINS==";\n'
    header += ('print "%s";\n' % 
               algebra.magma.quote_string(open(triangulation_filename,'r').read()))
    header += 'print "==TRIANGULATION=ENDS==";\n'
    return header

def get_outfile_base(triangulation_filename):
    base = os.path.abspath(triangulation_filename)
    directory, filename = os.path.split(base)
    parentdir = os.path.dirname(directory)
    if filename[-5:].lower() == '.trig':
        filename = filename[:-5]
    target_dir = os.path.join(parentdir, "magma_slN")
    if os.path.isdir(target_dir):
        outfile_base = os.path.join(target_dir, filename)
    else:
        outfile_base = os.path.join(directory, filename)
        #print "You might want to create a seperate directory for the magma files with"
        #print "   mkdir", target_dir
    return outfile_base

def hash_ideal(eqns, term_order):
    e = [str(p) for p in eqns]
    s = 'SEP'.join(e+[' VAR ']+term_order)
    return hashlib.md5(s).hexdigest()[-6:]

main()
