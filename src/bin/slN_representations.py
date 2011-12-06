#!/usr/bin/python
import os
import sys
import re
import optparse
import traceback
import StringIO
import csv

this_path, this_file = os.path.split(sys.argv[0])
abs_path = os.path.abspath(this_path)
base_path, this_dir = os.path.split(abs_path)
sys.path.append(base_path)

try:
    from manifold.triangulation import triangulation, read_triangulation_from_file
    from algebra.polynomial import Polynomial
    import manifold.slN
    import algebra.magma
    import hashlib
except ImportError as e:
    print e
    print
    print "This program was called as       :", sys.argv[0]
    print "Absolute path to this program is :", abs_path
    print "Base path is                     :", base_path
    sys.exit(1)


def load_pari(options):
    try:
        import algebra.pari as pari
        from algebra.pari import number, set_pari_precision, set_pari_allowed_error, get_pari_allowed_error, NumericalError
        from algebra.solve_polynomial_equations import solve_polynomial_equations

        global pari
        global number, set_pari_precision, set_pari_allowed_error, get_pari_allowed_error, NumericalError
        global solve_polynomial_equations
    except ImportError as e:
        print "Error while loading pari"
        print "It is necessary to compile the cython code, see README"
        print
        print e
        print
        print "This program was called as       :", sys.argv[0]
        print "Absolute path to this program is :", abs_path
        print "Base path is                     :", base_path
        sys.exit(1)

    # Set precision for PARI
    if options.precision and options.max_err:
        set_pari_precision(options.precision)
        set_pari_allowed_error(options.max_err)
    elif options.precision:
        set_pari_precision(options.precision)
        set_pari_allowed_error(options.precision - 6)
    elif options.max_err:
        set_pari_precision(options.max_err + 6)
        set_pari_allowed_error(options.max_err)
    else:
        set_pari_precision(30)
        set_pari_allowed_error(20)

        

CSVHeader = ["Manifold",
             "File",
             "Ordered",
             "Tetrahedra",
             "Cusps",
             "SL(N,C)",
             "Obstruction Class",
             "Index Component",
             "Number Components",
             "Dimension",
             "Additional Eqns",
             "Number Solutions",
             "Index Solution",
             "Warning",
             "Volume",
             "CS",
             "CPUTIME"]

def main():
    parser = create_parser()
    # Parse command line options
    options, args = parser.parse_args()

    # Exit writing the header for the CSV file
    if options.print_csv_header:
        output = StringIO.StringIO()
        csv_writer = csv.DictWriter(output, fieldnames = CSVHeader)
        csv_writer.writerow(dict(zip(CSVHeader, CSVHeader)))
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

    elif "MAGMA=OUTPUT=FOR=SNAPREPR" in input_file:
        # magma breaks long lines putting a \ at the break
        # we replace those line breaks
        load_pari(options)
        magma_out = ''.join([x.strip() for x in input_file.split('\\\n')])

        process_magma_output_headers(options, magma_out, 
                                     magma_filename = args[0])

    elif "MAGMA=INPUT=FROM=SNAPREPR" in input_file:

        print "This file is supposed to be processed with MAGMA"
        print "Run:"
        print "      magma_out %s" % args[0]
        print "Then:"
        print "      slN_representations.py %s_out" % args[0]

    else:

        sys.stderr.write("Could not recognize a SnapPea triangulation file or a MAGMA output file ????\n\n")
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
    for c, h in zip(range(len(H)), H):
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
    error_condition = ""
    output = StringIO.StringIO()
    csv_writer = csv.DictWriter(output, 
                                fieldnames = CSVHeader,
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
                                  '(.*)'
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
    # A prime_ideal object is derived from a list, but also
    # has a dimension and number_of_points field.
    try:
        prim_decomp = algebra.magma.parse_primary_decomposition(magma_out)
        csv_dict["Number Components"] = len(prim_decomp)
    except:
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
        for index_prime_ideal, prime_ideal in zip(
            range(len(prim_decomp)),prim_decomp):

            csv_dict['Index Component'] = index_prime_ideal
            csv_dict["Dimension"] = prime_ideal.dimension
            if not prime_ideal.number_of_points is None:
                csv_dict['Number Solutions'] = prime_ideal.number_of_points

            if prime_ideal.dimension == 0:
                try:
                    cvols = get_complex_volumes(t, N, c, prime_ideal, not_paranoid = not options.paranoid)

                    for i, cvol in zip(range(len(cvols)),cvols):
                        csv_dict["Index Solution"] = i
                        rvol = cvol.real()
                        if rvol.abs() < get_pari_allowed_error():
                            rvol = 0
                        csv_dict["Volume"] = str(rvol)
                        csv_dict["CS"] = str(cvol.imag())
#                        csv_writer.writerow(csv_dict)

                except NumericalError as n:
                    csv_dict["Warning"] = "Numerical Error"
#                    csv_writer.writerow(csv_dict)
                    sys.stderr.write(traceback.format_exc())
                    sys.stderr.write("\n\n\n")
                except:
                    csv_dict["Warning"] = "UNKNOWN ERROR"
#                    csv_writer.writerow(csv_dict)
                    sys.stderr.write(traceback.format_exc())
                    sys.stderr.write("\n\n\n")
            csv_writer.writerow(csv_dict)

    print output.getvalue()[:-1]

def get_complex_volumes(t, N, c, prime_ideal, not_paranoid = False):

    all_cvols = []

    # Find the points in the variety
    # solve_polynomial_equations assumes that prime_ideal is given
    # in Groebner basis
    sys.stderr.write("Solving...\n")
    solutions = solve_polynomial_equations(prime_ideal)
    sys.stderr.write("Solved...\n")

    # Solutions is a list with one element per Galois conjugate
    if not len(solutions) == prime_ideal.number_of_points:
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
            number("1.0+0.32433234644*I"))

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
        cvol = sum(cvols) / pari.I

        sys.stderr.write("   %s\n" % cvol)

        # Consistency check
        if not (vol - cvol.real()).abs() < get_pari_allowed_error():
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
                      dest = "max_err", default = None,
                      help = "Digits of precision")
    parser.add_option("-P", "--internal-precision",
                      type = "int",
                      dest = "precision", default = None,
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
        print "You might want to create a seperate directory for the magma files with"
        print "   mkdir", target_dir
    return outfile_base

def hash_ideal(eqns, term_order):
    e = [str(p) for p in eqns]
    s = 'SEP'.join(e+[' VAR ']+term_order)
    return hashlib.md5(s).hexdigest()[-6:]

main()
