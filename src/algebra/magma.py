import re
from algebra.polynomial import Polynomial

def quote_string(s):
    r = ""

    for i in s.split('\n'):
        for j in range(1+(len(i)/60)):
            if j == len(i)/60:
                r += i[60*j:60*j+60] + '\\n'
            else:
                r += i[60*j:60*j+60] + '\\\\\\n'
    return r


def ideal_to_magma(polys, term_order):
    """
    Turns a list of polynomials which are considered to be generators of an
    ideal into a format suitable for MAGMA.

    >>> polys = [polynomial("x*y"), polynomial("x*z")]
    >>> print ideal_to_magma(polys)
    P<x, y, z> :=   PolynomialRing(RationalField(), 3);
    I := ideal<P |
          x * y,
          x * z>;
    <BLANKLINE>
    """
    return (  "P<%s> := " % ", ".join(term_order)
            + "  PolynomialRing(RationalField(), %d);\n" % len(term_order)
            + "I := ideal<P |\n"
            +    "      "
            + ",\n      ".join(map(str, polys))
            + ">;\n")


def primary_decomposition(polys, term_order):
    """
    Similar to ideal_to_magma but includes the MAGMA command for Primary
    Decomposition


    >>> polys = [polynomial("x*y"), polynomial("x*z")]
    >>> print primary_decomposition(polys)
    P<x, y, z> :=   PolynomialRing(RationalField(), 3);
    I := ideal<P |
          x * y,
          x * z>;
    <BLANKLINE>
    <BLANKLINE>
    print "PRIMARY=DECOMPOSITION=BEGINS=HERE";
    PrimaryDecomposition(I);
    print "PRIMARY=DECOMPOSITION=ENDS=HERE";
    <BLANKLINE>
    <BLANKLINE>
    <BLANKLINE>
    """
    
    return (  ideal_to_magma(polys, term_order) + "\n\n"
            + 'cputime := Cputime();\n'
            + 'print "PRIMARY=DECOMPOSITION=BEGINS=HERE";\n'
            + "PrimaryDecomposition(I);\n"
            + 'print "PRIMARY=DECOMPOSITION=ENDS=HERE";\n\n\n'
            + 'print "CPUTIME          :", Cputime(cputime);')

def ideal_to_magma_curve(polys, term_order):
    return (  "A<%s> := " % ", ".join(term_order)
            + "  AffineSpace(Rationals(), %d);\n" % len(term_order)
            + "C := Curve(A,[\n"
            +    "      "
            + ",\n      ".join(map(str, polys))
            + "]);\n")
    

def genus(polys, pre_vars = [], post_vars = []):
    """
    Similar to ideal_to_magma but includes the MAGMA command for Genus
    """
    
    return (  ideal_to_magma_curve(polys, pre_vars, post_vars) + "\n\n"
            + 'print "GENUS=BEGINS=HERE";\n'
            + "GeometricGenus(C);\n"
            + 'print "GENUS=DECOMPOSITION=ENDS=HERE";\n\n\n')



class prime_ideal(list):
    def __init__(self,l,dim, number_of_points = None):
        super(prime_ideal,self).__init__(l)
        self.dimension = dim
        if not number_of_points == '':
            self.number_of_points = int(number_of_points)
        else:
            self.number_of_points = None

def parse_primary_decomposition(s_with_backslash):
    """
    >>> p = parse_primary_decomposition(_magma_sample)
    >>> len(p)
    1
    >>> p = p[0]
    >>> len(p)
    3
    >>> str(p[0])
    '-  1  + c_01_0'
    >>> str(p[1])
    '1  + c_02_0 + t'
    >>> str(p[2])
    '1  + t + t^2'
    """

    s = re.search(
        "PRIMARY=DECOMPOSITION=BEGINS=HERE\s*"
        "("
        "\[\s*"
        "(([^\]]*\[[^\]]*\],?\s*)*)"
        "\]\s*"
        ")+"
        "PRIMARY=DECOMPOSITION=ENDS=HERE",
        s_with_backslash,re.DOTALL)

    if not s:
        raise Exception, "Invalid Primary Decomposition"

    s_with_backslash = s.group(2)
    
    s=""
    for line in s_with_backslash.split('\n'):
        line=line.strip()
        if line and line[-1]=='\\':
            line=line[0:-1]
        else:
            line=line+' '
        s = s + line

#    print s

    

    components = re.findall(
        "Ideal of Polynomial ring.*?"
        "Dimension (\d+).*?"
        "(Size of variety over algebraically closed field: (\d+).*?)?"
        "Groebner basis:\s*"
        "\[([^\]]*)\]",
        s)

    print "There"

    print components

    for dimension, tmp, nu, polys in components:
        for p in polys.split(','):
            print "Str"
            print p
            print "Poly"
            print Polynomial.parseFromString(p)
            print "What?"

    components= [
        prime_ideal(
            l = [Polynomial.parseFromString(p) for p in polys.split(',')],
            dim = int(dimension), number_of_points = number_of_points)
        for dimension, tmp, number_of_points, polys in components]

    print "asd"

    print components

    return components

if __name__=='__main__':
    import sys
    import solve_polynomial_equations
    l=open(sys.argv[1],'r').read()
    g_basis=parse_magma_primary_decomposition(l)
    for i in range(len(g_basis)):
        for j in range(len(g_basis[i])):
            print "Magma",g_basis[i][j]
            print "variables",g_basis[i][j].variables()
            print "univariant",g_basis[i][j].is_univariant()
            print "degree",g_basis[i][j].degree()
            if g_basis[i][j].is_univariant():
                print "coefficients",g_basis[i][j].coefficients()
    for a_g_basis in g_basis:
        print solve_polynomial_equations.solve_polynomial_equations(a_g_basis)
    
_magma_sample = """
Magma V2.13-6     Tue May 25 2010 18:05:46 on matthias-laptop [Seed = 531113876]
Type ? for help.  Type <Ctrl>-D to quit.
PRIMARY=DECOMPOSITION=BEGINS=HERE
[
    Ideal of Polynomial ring of rank 3 over Rational Field
    Lexicographical Order
    Variables: c_01_0, c_02_0, t
    Dimension 0, Radical, Prime
    Size of variety over algebraically closed field: 2
    Groebner basis:
    [
        c_01_0 - 1,
        c_02_0 + t + 1,
        t^2 + t + 1
    ]
]
[
    Ideal of Polynomial ring of rank 3 over Rational Field
    Lexicographical Order
    Variables: c_01_0, c_02_0, t
    Dimension 0, Radical, Prime
    Size of variety over algebraically closed field: 2
    Groebner basis:
    [
        c_01_0 - 1,
        c_02_0 + t + 1,
        t^2 + t + 1
    ]
]
PRIMARY=DECOMPOSITION=ENDS=HERE

Total time: 0.839 seconds, Total memory usage: 5.62MB
"""
