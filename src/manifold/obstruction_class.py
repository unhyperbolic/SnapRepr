
#from manifold3D.triangulation import left_out_number,permutation,triangulation

from manifold.triangulation import triangulation
from algebra.matrix import matrix
from algebra.homology import homology_generators


"""
Given a triangulation t, the basis for the chain complex are:
C_3 the tetrahedra itself
C_2 the face classes from t.get_face_classes()
C_1 the edge classes from t.get_edge_classes()

"""

def d3(t):

    """
    the differential C_3 -> C_2 of chains with Z coefficients as a matrix
    """
    
    tets = t.tet_list
    faces = t.get_face_classes()

    def boundary(tet,a_face_class):
        res = 0
        if tet.index==a_face_class.tet1():
            res += +1
            #res += (-1)**a_face_class.face1()
        if tet.index==a_face_class.tet2():
            res += a_face_class.orient()
            #res += ((-1)**a_face_class.face2()) * a_face_class.orient()
        return res

    return matrix([[boundary(i,j) for i in tets] for j in faces],int)

    return m
    
def d2(t):
    """
    the differential C_2 -> C_1 of chains with Z coefficients as a matrix
    """

    faces = t.get_face_classes()
    edges = t.get_edge_classes()

    def boundary(face,an_edge_class):
        res = 0

        for an_edge in an_edge_class:
            if face.tet1()==an_edge.tet():
                if not (face.face1()==an_edge.vert_0() or
                        face.face1()==an_edge.vert_1()):
                    
                    if an_edge.vert_0()>face.face1():
                        vert0 = an_edge.vert_0()-1
                    else:
                        vert0 = an_edge.vert_0()
                    if an_edge.vert_1()>face.face1():
                        vert1 = an_edge.vert_1()-1
                    else:
                        vert1 = an_edge.vert_1()
                    if (vert0,vert1) in [(0,1),(1,2),(2,0)]:
                        res += + (-1) ** face.face1()
                    else:
                        assert (vert0,vert1) in [(1,0),(2,1),(0,2)]
                        res += - (-1) ** face.face1()
        return res

    return matrix([[boundary(i,j) for i in faces] for j in edges],int)

def cohomology_2_rel_boundary(t, field, as_matrix_with_column_vectors = False):
    """
    computes H^2(t,boundary of t; field)
    i.e. t is a manifold with boundary (cusps are cut), coefficients are in field
    
    >>> from algebra.field_p import field_p
    >>> from fractions import Fraction
    >>> Fraction.one = staticmethod(lambda : Fraction(1,1))
    >>> Fraction.zero = staticmethod(lambda : Fraction(0,1))
    >>> def is_zero(m):
    ...     for r in m.values:
    ...         for c in r:
    ...             if not c == 0:
    ...                 return False
    ...     return True
    >>> from manifold3D.triangulation import *
    >>> t = read_triangulation_from_file("lib/triangulations/m003.trig")
    >>> p = d2(t) * d3(t)
    >>> is_zero(p)
    True

    H_1(m003)=H^2(m003,boundary) = Z + Z/5,
    dim(H^2(m003,boundary;Q)) = 1
    
    >>> len(cohomology_2_rel_boundary(t, Fraction))
    1
    >>> t = read_triangulation_from_file("lib/triangulations/m206.trig")
    >>> is_zero(d2(t) * d3(t))
    True
    >>> field = field_p(5)

    H_1(m206)=H^2(m206,boundary) = Z + Z/5,
    dim(H^2(m206,boundary;Z/5)) = 2
    
    >>> len(cohomology_2_rel_boundary(t, field))
    2
    >>> t = read_triangulation_from_file("lib/triangulations/m053.trig")
    >>> is_zero(d2(t) * d3(t))
    True
    >>> len(cohomology_2_rel_boundary(t, Fraction))
    1
    >>> cohomology_2_rel_boundary(t, Fraction, True).no_columns()
    1
    """

    assert isinstance(t,triangulation)
    
    m_d3 = matrix(d3(t), field).transpose()
    m_d2 = matrix(d2(t), field).transpose()
    
    H = homology_generators(m_d2, m_d3, as_matrix_with_column_vectors)

    return H


def cohomology_2_rel_boundary_classes(t, field):
    H = cohomology_2_rel_boundary(t, field,
                                  as_matrix_with_column_vectors = True)

    return [H * v for v in all_vectors(field, H.no_columns())]


def all_vectors(field, length):
    """
    Lists all vectors of field ^ length.

    >>> from algebra.field_p import field_p
    >>> all_vectors(field_p(2),2)
    [[field_p(2)(0), field_p(2)(0)], [field_p(2)(0), field_p(2)(1)], [field_p(2)(1), field_p(2)(0)], [field_p(2)(1), field_p(2)(1)]]

    There are 625 vectors in (Z/5)^4
    
    >>> len(all_vectors(field_p(5),4))
    625
    """

    if length == 0:
        return [ [] ]
    try:
        elems = field.elements()
    except:
        raise Exception("field has no method elements, field is %s of type: %s" % (field, type(field)))

    v = all_vectors(field, length - 1)

    r = []
    for e in elems:
        r = r + [[e]+x for x in v]
    return r

def cohomology_2_rel_boundary_class_to_coeffs(trig, cohomology_class):
    """
    Given a cohomology class as outputed by cohomology_2_rel_boundary
    this will return cohomology_coefficients such that
    cohomology_coefficient[tet_index][face_index] will give the coefficient
    of the face face_index of tetrahedron tet_index

    """
    assert isinstance(trig, triangulation)
    assert isinstance(cohomology_class, list)
    
    face_classes = trig.get_face_classes()
    assert len(cohomology_class) == len(face_classes)

    cohomology_coefficients = [ [None, None, None, None] for tet in trig.tet_list]

    for f_class, coeff in zip(face_classes,cohomology_class):
        cohomology_coefficients[f_class.tet1()][f_class.face1()] = coeff
        if f_class.orient == +1:
            cohomology_coefficients[f_class.tet2()][f_class.face2()] = coeff
        else:
            cohomology_coefficients[f_class.tet2()][f_class.face2()] = -coeff
            
    return cohomology_coefficients


def assign_cohomology_2_rel_boundary_to_tets(trig, cohomology_class):

    """
    Given a cohomology class (cohomology_2_rel_boundary gives a list of
    generators for these), this will store the coefficients in
    cohomology_coefficients of each tetrahedron in the triangulation.
    
    
    >>> from algebra.field_p import field_p
    >>> from fractions import Fraction
    >>> Fraction.one = staticmethod(lambda : Fraction(1,1))
    >>> Fraction.zero = staticmethod(lambda : Fraction(0,1))
    >>> from manifold3D.triangulation import *
    >>> t1 = read_triangulation_from_file("lib/triangulations/m003.trig")
    >>> t2 = read_triangulation_from_file("lib/triangulations/m053.trig")
    >>> t3 = read_triangulation_from_file("lib/triangulations/m032.trig")
    >>> t4 = read_triangulation_from_file("lib/ordered_triangulations/m032_ordered.trig")    
    >>> for field in [Fraction, field_p(2), field_p(5)]:
    ...     for t in [t1, t2,t3,t4]:
    ...         H = cohomology_2_rel_boundary(t, field)
    ...         for h in H:
    ...             assign_cohomology_2_rel_boundary_to_tets(t, h)
    ...             for tet in t.tet_list:
    ...                 c = tet.cohomology_coefficients
    ...                 assert c[0] + c[1] + c[2] + c[3] == field.zero()
    
    """

    assert isinstance(trig, triangulation)
    assert isinstance(cohomology_class, list)
    
    face_classes = trig.get_face_classes()
    assert len(cohomology_class) == len(face_classes)

    cohomology_coefficients = cohomology_2_rel_boundary_coeffs_per_tet(trig, cohomology_class)
    
    for tet in trig.tet_list:
        tet.cohomology_coefficients = cohomology_coefficients[tet.index]
