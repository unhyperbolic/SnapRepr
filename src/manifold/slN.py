# This file is structured as follows:
# PART I 
# Procedures needed to produce the Ptolemy variety
# Part II
# Procedures needed to process the solutions of the Ptolemy variety
# Part III
# Procedures to check consistency

from manifold.obstruction_class import cohomology_2_rel_boundary_classes
from manifold.obstruction_class import cohomology_2_rel_boundary_class_to_coeffs
from manifold.bloch_group import Ptolemy_cochain
from algebra.field_p import field_p
from algebra.polynomial import Polynomial
from algebra.pari import pari_eval, pari_eval_bool, number, get_pari_allowed_error, NumericalError

from fractions import Fraction

# PART I
# Procedures needed to produce the Ptolemy variety

class equivalence_with_sign:

    """
    Encapsulates an equivalence relationship between variables with signs,
    e.g., if x is identified with -y, and y with z,
    then x will be identified  with -z automatically.
    This class will fail if x is identified with -x.
    
    The internal representation is as a dictionary of dictionaries:
    d : variable |--->  (sign |-----> list of equivalent variables)

    >>> e=equivalence_with_sign()
    >>> e.identify('a','c')
    >>> e.identify('x','y',-1)
    >>> e.identify('a','d',-1)
    >>> e.identify('a','c')
    >>> e.dict_canonical_representatives()
    {'a': (1, 'a'), 'x': (1, 'x'), 'c': (1, 'a'), 'y': (-1, 'x'), 'd': (-1, 'a')}
    """

    def __init__(self):
        self.d = {}
    def identify(self, var_a, var_b, sign = +1):
        assert sign in [-1, +1]
            
        if not self.d.has_key(var_a):
            self.d[var_a] = { +1 : set([var_a]), -1 : set() }
        if not self.d.has_key(var_b):
            self.d[var_b] = { +1 : set([var_b]), -1 : set() }

        new_val_a = { + 1: self.d[var_a][+1] | self.d[var_b][ sign],
                      - 1: self.d[var_a][-1] | self.d[var_b][-sign]}

        new_minus_a = { +1: new_val_a[-1],
                        -1: new_val_a[+1]}

        for n in [new_val_a, new_minus_a]:
            for x in n[+1]:
                self.d[x] = n
    
    def dict_canonical_representatives(self):
        def canonical_representative(d):
            all_keys = list(d[+1] | d[-1])
            all_keys.sort()
            key = all_keys[0]
            if key in d[+1]:
                return (+1,key)
            else:
                return (-1,key)
        return dict([(k,canonical_representative(v)) for k,v in self.d.items()])

    def dict_canonical_representatives_poly(self):
        def turn_pair_to_Polynomial(p):
            if p[0] == +1:
                return Polynomial.fromVariableName(p[1])
            else:
                return -Polynomial.fromVariableName(p[1])
        return dict([(k,turn_pair_to_Polynomial(v)) for k,v in self.dict_canonical_representatives().items()])

    def pretty_print(self):
        d = self.dict_canonical_representatives()
        keys = d.keys()
        keys.sort()
        def print_pair(key, d = d):
            if d[key][0] == +1:
                return "%s :  %s" % (key, d[key][1])
            else:
                return "%s : -%s" % (key, d[key][1])

        return "\n".join([print_pair(key) for key in keys])

class simplex_coords(tuple):
    """
    Encapsulates integer coordinates (t0,t1,...,tn)
    It is derived from tuple, hence supports len, ...
    It provides two additional methods: 

    __add__ for adding another tuple
    >>> simplex_coords((0,4,1,1)) + (0,1,1,0)
    (0, 5, 2, 1)

    inclusion_onto_boundary
    Consider the i-th face of n+1 simplex. If T is a point on the face,
    then T.inclusion_onto_boundary(i) is the image of T under the inclusion
    of the i-th face of the n+1 simplex.
    >>> simplex_coords((6, 5, 2)).inclusion_onto_boundary(1)
    (6, 0, 5, 2)
    """

    def __new__(cls, z):
        return super(simplex_coords, cls).__new__(cls,z)
    def __add__(self, other):
        return simplex_coords(map(lambda x,y:x+y,self,other))
    def inclusion_onto_boundary(self,i):
        return simplex_coords(self[:i]+(0,)+self[i:])


def tuples_with_fixed_sum(n, s):
    """
    Returns a list of all tuples of length n with sum s
    
    The result is a python native tuple.

    >>> tuples_with_fixed_sum(3, 2)
    [(0, 0, 2), (0, 1, 1), (0, 2, 0), (1, 0, 1), (1, 1, 0), (2, 0, 0)]
    """

    if n == 0:
        if s == 0:
            return [()]
        else:
            return []
    if n == 1:
        return [(s,)]
    r = []
    for i in range(0, s + 1):
        for j in tuples_with_fixed_sum(n - 1, s - i):
            r.append( (i,) + j)
    return r

def simplex_coords_with_fixed_sum(n, s):
    """
    Same as tuples_with_fixed_sum, but returns it as simplex_coords.
    """

    return map(simplex_coords, tuples_with_fixed_sum(n, s))

def c_parameter_var(coords, tet_index):
    """
    Takes simplex_coords and the index of a tetrahedron and return the
    respective variable name for the Ptolemy coordinate.

    >>> c_parameter_var(simplex_coords((1,2,3,1)),5)
    'c_1231_5'
    """

    c = "c_"
    for i in coords:
        if i > 9:
            c += "_"
        c += "%d" % i
    c += "_%d" % tet_index
    return c

def c_parameter(coords, tet_index):
    """
    Same as c_parameter_var but returns it as Polynomial type.
    """
    return Polynomial.fromVariableName(
        c_parameter_var(coords, tet_index))

def get_all_obstruction_classes(t):
    """
    Computes all classes in H^2(M, partial M; Z/2)
    """
    return cohomology_2_rel_boundary_classes(t, field_p(2))

def Z2_to_sign(f):
    """
    Helper function turning the result of get_all_obstruction_classes
    which uses the type field_p into +1, -1 signs.

    >>> Z2_to_sign(field_p(2)(0))
    1
    >>> Z2_to_sign(field_p(2)(1))
    -1
    """
    if f == field_p(2)(0):
        return + 1
    elif f == field_p(2)(1):
        return - 1
    raise Exception, "Z2_to_sign expects field_p(2)(0) or field_p(2)(1)"

def get_Ptolemy_relations(t, N, cohomology_class = None):
    if cohomology_class:
        cohomology_coefficients = (
            cohomology_2_rel_boundary_class_to_coeffs(t, cohomology_class))

    eqns = []
    
    for tet in t.tet_list:
        i = tet.index
        if cohomology_class:
            signs = cohomology_coefficients[tet.index]
            sign_01 = Z2_to_sign(signs[2] + signs[3])
            sign_12 = Z2_to_sign(signs[0] + signs[3])
        else:
            sign_01 = +1
            sign_12 = +1

        sign_01 = Polynomial.constantPolynomial(sign_01)
        sign_12 = Polynomial.constantPolynomial(sign_12)

        for coord in simplex_coords_with_fixed_sum(4, N-2):
            eqns.append(
                - sign_01 * c_parameter(coord+(1,0,0,1),i)
                          * c_parameter(coord+(0,1,1,0),i)
                - sign_12 * c_parameter(coord+(1,1,0,0),i)
                          * c_parameter(coord+(0,0,1,1),i)
                +           c_parameter(coord+(1,0,1,0),i)
                          * c_parameter(coord+(0,1,0,1),i))
    return eqns

def find_independent_integer_vectors(d):
    
    d = dict([ (k, [ Fraction(x,1) for x in v ]) for k, v in d.items() ] )

    keys = d.keys()
    keys.sort()
    
    s = len(d[keys[0]])

    res_keys = []

    m = []

    while len(m) < s:
        found_new_row = False
        for k in keys:
            new_row = d[k]
            for index, old_row in zip(range(len(m)), m):
                new_row = [ x - y * new_row[index]  
                            for x, y in zip(new_row, old_row)]
            
            if new_row[len(m)]:
                m.append([x / new_row[len(m)] for x in new_row])
                res_keys.append(k)
                found_new_row = True
                break
        assert found_new_row
            
    return res_keys


def get_additional_eqns_independent(t, N):
    
    # make sure cusps have indices

    # dictionary associating the decoration change vectors to Ptolemy coordinates
    d = {}

    for i, tet in zip(range(len(t.tet_list)), t.tet_list):
        for coord in simplex_coords_with_fixed_sum(4, N):
            if not N in coord:
                v = [0 for x in range((N-1) * t.num_or_cusps)]
                for j, k in zip(range(4), coord):
                    for p in range(k):
                        v[ tet.cusp_index[j] + p * t.num_or_cusps ] += 1
                d[c_parameter_var(coord, i)] = v

    variables = find_independent_integer_vectors(d)
    
    return [
        Polynomial.fromVariableName(var) - Polynomial.constantPolynomial(1) 
        for var in variables]

def get_additional_eqns_manual(t, N):
    # sets any additional terms to one

    eqns = []

    if N == 2:
        eqns.append(Polynomial("c_1100_0") - 1)
    elif N == 3:
        eqns.append(Polynomial("c_2100_0") - 1)
        eqns.append(Polynomial("c_1110_0") - 1)
    elif N == 4:
        eqns.append(Polynomial("c_3100_0") - 1)
        eqns.append(Polynomial("c_2200_0") - 1)
        eqns.append(Polynomial("c_2110_0") - 1)
    elif N == 5:
        eqns.append(Polynomial("c_4100_0") - 1)
        eqns.append(Polynomial("c_3110_0") - 1)
        eqns.append(Polynomial("c_2111_0") - 1)
        eqns.append(Polynomial("c_2210_0") - 1)
    else:
        raise "Not yet implemented"
           
    return eqns

def get_identified_c_parameters(t, N):
    e = equivalence_with_sign()

    K = simplex_coords_with_fixed_sum(3, N)
    
    M = simplex_coords_with_fixed_sum(4, N)

    for tet in t.tet_list:
        for coord in M:
            if 0 not in coord:
                self_var = c_parameter_var(coord, tet.index)
                e.identify(self_var, self_var, sign = +1)
        for face in range(4):
            for coord in [x.inclusion_onto_boundary(face) for x in K]:
                if not N in coord:
                    self_var = c_parameter_var(coord, tet.index)
                    
                    adj_tet = tet.neighbor_index[face]
                    adj_coord = simplex_coords([coord[x]
                                                for x in tet.gluing[face].inverse()])
                    adj_var = c_parameter_var(adj_coord, adj_tet)

                    new_perm = []

                    for i in range(4):
                        if coord[i] % 2:
                            new_perm.append(tet.gluing[face][i])
                    n=[x for x in new_perm]
                    n.sort()
                    d = dict(zip(n,range(len(n))))
                    new_perm = [d[x] for x in new_perm]
                    l = 0
                    for i in [(0,1),(0,2),(1,2)]:
                        if i[1] < len(new_perm):
                            if new_perm[i[0]] > new_perm[i[1]]:
                                l = l + 1
                            
                    e.identify(self_var, adj_var, sign = (-1) ** l)

    return e

# gets eqns and equivalence relationship
def identify_c_parameters(eqns, e):
    d = e.dict_canonical_representatives_poly()
    return [eqn.substitute(d) for eqn in eqns]

def get_all_variables(poly_list):
    variables = []
    for p in poly_list:
        assert isinstance(p, Polynomial)
        variables = variables + p.variables()
    variables = list(set(variables))
    variables.sort()
    return variables

def polynomialNonZeroCondition(eqns, var):
    variables = get_all_variables(eqns)

    eqns = [e for e in eqns]

    prod = Polynomial.constantPolynomial(1)
    for i in variables + [var]:
        prod = prod * Polynomial.fromVariableName(i)

    return prod - Polynomial.constantPolynomial(1)

# Part II
# Procedures needed to process the solutions of the Ptolemy variety

def map_solution_to_c_parameters(solution, e):
    assert isinstance(solution, dict)
    assert isinstance(e, equivalence_with_sign)

    c_params = {}

    for var, signed_representative in e.dict_canonical_representatives().items():
        sign, representative = signed_representative
        if sign == +1:
            c_params[var] =  solution[representative]
        else:
            c_params[var] = -solution[representative]

    return c_params

def multiply_solution_by_constant(d, n):
    return dict([(k, v*n) for k, v in d.items()])

def get_Ptolemy_cochain(t, N, c_params, no_check = False):
    
    t.check_consistency()

    K = simplex_coords_with_fixed_sum(4, N - 2)

    elements = []
    for i, tet in zip(range(len(t.tet_list)), t.tet_list):
        for coord in K:
            elements.append(
                Ptolemy_cochain(
                    sign = (+1 if tet.positive_orientation else -1),
                    c01 = c_params[c_parameter_var(coord+(1,1,0,0),i)],
                    c02 = c_params[c_parameter_var(coord+(1,0,1,0),i)],
                    c03 = c_params[c_parameter_var(coord+(1,0,0,1),i)],
                    c12 = c_params[c_parameter_var(coord+(0,1,1,0),i)],
                    c13 = c_params[c_parameter_var(coord+(0,1,0,1),i)],
                    c23 = c_params[c_parameter_var(coord+(0,0,1,1),i)],
                    no_check = no_check))
                      
    return elements

# Part III
# Procedures to check consistency



def check_solution_on_gluing_equations(t, N, solution, cohomology_class = None):

    if cohomology_class:
        cohomology_coefficients = (
            cohomology_2_rel_boundary_class_to_coeffs(t,cohomology_class))

    for tet in t.tet_list:
        i = tet.index
        if cohomology_class:
            signs = cohomology_coefficients[tet.index]
            
            sign_01 = Z2_to_sign(signs[2] + signs[3])
            sign_12 = Z2_to_sign(signs[0] + signs[3])

        else:
            sign_01 = +1
            sign_12 = +1

        for coord in simplex_coords_with_fixed_sum(4, N-2):
            if not pari_eval_bool(
                "abs(- (%s) * (%s) * (%s) - (%s) * (%s) * (%s) + (%s) * (%s)) < (%s)" % (
                    sign_01,
                    solution[c_parameter_var(coord+(1,0,0,1),i)],
                    solution[c_parameter_var(coord+(0,1,1,0),i)],
                    sign_12,
                    solution[c_parameter_var(coord+(1,1,0,0),i)],
                    solution[c_parameter_var(coord+(0,0,1,1),i)],
                    solution[c_parameter_var(coord+(1,0,1,0),i)],
                    solution[c_parameter_var(coord+(0,1,0,1),i)],
                    get_pari_allowed_error())):

                raise NumericalError(
                    val =
                    number(eval_this = "- (%s) * (%s) * (%s) - (%s) * (%s) * (%s) + (%s) * (%s)" % (
                            sign_01,solution[c_parameter_var(coord+(1,0,0,1),i)],
                            solution[c_parameter_var(coord+(0,1,1,0),i)],
                            sign_12,solution[c_parameter_var(coord+(1,1,0,0),i)],
                            solution[c_parameter_var(coord+(0,0,1,1),i)],
                            solution[c_parameter_var(coord+(1,0,1,0),i)],
                            solution[c_parameter_var(coord+(0,1,0,1),i)])),
                    msg = "in verify gluing equations")

def check_solution_identification(t, N, solution):    
    e = equivalence_with_sign()

    K = simplex_coords_with_fixed_sum(3, N)

    for tet in t.tet_list:
        for face in range(4):
            for coord in [x.inclusion_onto_boundary(face) for x in K]:
                if not N in coord:
                    self_var = c_parameter_var(coord, tet.index)
                    
                    adj_tet = tet.neighbor_index[face]
                    adj_coord = simplex_coords([coord[x]
                                                for x
                                                in tet.gluing[face].inverse()])
                    adj_var = c_parameter_var(adj_coord, adj_tet)

                    new_perm = []

                    for i in range(4):
                        if coord[i] % 2:
                            new_perm.append(tet.gluing[face][i])
                    n=[x for x in new_perm]
                    n.sort()
                    d = dict(zip(n,range(len(n))))
                    new_perm = [d[x] for x in new_perm]
                    l = 0
                    for i in [(0,1),(0,2),(1,2)]:
                        if i[1] < len(new_perm):
                            if new_perm[i[0]] > new_perm[i[1]]:
                                l = l + 1
                            
                    sign = (-1) ** l

                    if not pari_eval_bool("abs( (%s) * (%s) - (%s) ) < (%s)" % 
                                          (sign,
                                           solution[self_var],
                                           solution[adj_var],
                                           get_pari_allowed_error())):
                        raise NumericalError(val = [self_var, adj_var, sign, solution],
                                             msg = "in identifying equations")


def c_param_elements(t, N, solution):
    
    t.check_consistency()

    c_params = map_c_parameters_back(t, N, solution)

    elements = []
    for i in range(len(t.tet_list)):
        for coord in simplex_coords_with_fixed_sum(4, N - 2):
            elements.append(
                ( (+1 if t.tet_list[i].positive_orientation else -1) ,
                  [c_params[c_parameter_var(coord+(1,1,0,0),i)],
                   c_params[c_parameter_var(coord+(1,0,1,0),i)],
                   c_params[c_parameter_var(coord+(1,0,0,1),i)],
                   c_params[c_parameter_var(coord+(0,1,1,0),i)],
                   c_params[c_parameter_var(coord+(0,1,0,1),i)],
                   c_params[c_parameter_var(coord+(0,0,1,1),i)]]))
                      
    return elements


def complex_volume_old_and_obsolete(t, N, solution):

    complex_volume = 0.0 + 0.0j

    elements = bloch_group_elements(t, N, solution)

    for sign, element in elements:
        complex_volume += sign * bloch_group.L_function(element)

    return complex_volume / 1.0j

def complex_volume2_old_and_obsolete(t, N, solution):

    complex_volume = number("0.0 + 0.0 * I")

    elements = c_param_elements(t, N, solution)

    for sign, element in elements:
        w = c_params_to_w(element[0],element[1],element[2],element[3],element[4],element[5])
        zpq = w_to_zpq(w)
        v = zpq_L_function(zpq)
        
        complex_volume += number(sign) * v

    return complex_volume / number(eval_this="I")

def volume2_old_and_Obsolete(t, N, solution):

    volume = number("0.0")

    elements = c_param_elements(t, N, solution)

    for sign, element in elements:
        w = c_params_to_w(element[0],element[1],element[2],element[3],element[4],element[5])
        zpq = w_to_zpq(w)
        v = zpq_volume(zpq)
        
        volume += number(sign) * v

    return volume
    

def volume_old_and_obsolete(t, N, solution):
    
    volume = 0.0

    elements = bloch_group_elements(t, N, solution)

    for sign, element in elements:
        volume += sign * bloch_group.bloch_volume(element)

    return volume
    
    
class Ptolemy_simplex:
    def __init__(self,c01,c02,c03,c12,c13,c23,sign=+1):
        self.c01 = c01
        self.c02 = c02
        self.c03 = c03
        self.c12 = c12
        self.c13 = c13
        self.c23 = c23
        self.sign = sign
    def check_consistency(self):
        pass
#    def __str__(self):
#        return "Ptolemy_simplex(c01=%s
