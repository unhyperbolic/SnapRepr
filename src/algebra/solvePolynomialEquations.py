import numpy
import math
from algebra.polynomial import Polynomial, uncomparablePrintCoefficientMethod
#from algebra.pari import roots_of_polynomial
#from algebra.pari import random_complex_modulos
#from algebra.pari import number

from fractions import Fraction

import globalsettings

globalsettings.registerSetting("solvePolynomialEquationsLog", True)

class SolverException(Exception):
    def __init__(self, message, poly_hist):
        self.poly_hist = poly_hist
        self.msg = message
    def __str__(self):
        return self.msg + "\nHistory of polynomials:\n" + self.poly_hist

def _printPoly(p):
    return p.printMagma(
        printCoefficientMethod = uncomparablePrintCoefficientMethod)

def filterPoly(polys, skip):
    return [poly for poly in polys if not poly == skip]

def solvePolynomialEquationsExactly(polys,
                                    nf = None, 
                                    variableDict = None):

    # nf is a polynomial in x encoding a number field

    # nf = x^2 + 1

    # variable contains the variables already bound as polynomials in the variable of the number field

    # a : x + 2

    if not variableDict:
        variableDict = {}

    if not polys:
        return variableDict, nf

    univariatePolys = [poly for poly in polys if poly.isUnivariate()]
    linearPolys = [poly for poly in polys if poly.isLinear()]

    if linearPolys:
        linearPoly = linearPolys[0]

        factor, constant = linearPoly.getCoefficients(Fraction)
        variable = linearPoly.variables()[0]

        if not isinstance(constant, Polynomial):
            constant = Polynomial.constantPolynomial(constant)

        factor = Polynomial.constantPolynomial(-Fraction(1)/factor)

        variableDict[variable] = factor * constant

        return solvePolynomialEquationsExactly(
            filterPoly(polys, linearPoly),
            nf,
            variableDict)

    if univariatePolys:
        univariatePoly = univariatePolys[0]
        variable = univariatePoly.variables()[0]

        if not nf:
            nf = univariatePoly.substitute({variable: Polynomial.fromVariableName('x')})
        else:

            univariatePoly = univariatePoly.subsitute({variable: Polynomial.fromVariableName('y')})

            pariNf = nfinit(nf)
            nf, transform = rnfequation(pariNf, univariatePoly)

            for k in variableDict.keys():
                variableDict[k] = variableDict[k].substitute({'x' : transform})

        variableDict[variable] = Polynomial.fromVariableName('x')
            
        return solvePolynomialEquationsExactly(
            filterPoly(polys, univariatePoly),
            nf,
            variableDict)

    raise Exception, "Should never get here"


# fills free variables with random values

def solvePolynomialEquations(polys,
                             polynomialSolver,
                             free_dim = 0,
                             with_poly_history = False,
                             poly_history="",
                             variable_dict = { },
                             non_linear_equation_encountered=False):
    
#    polys = [polynomial.convertCoefficients(number) for polynomial in polys]

    if globalsettings.getSetting("solvePolynomialEquationsLog"):
        poly_history += '\n\n\n\n'+'\n'.join(map(_printPoly,polys))+'\n\n============\n'

    if not polys:
        assert free_dim == 0
        if with_poly_history:
            return [(variable_dict,poly_history)]
        else:
            return [variable_dict]
    solutions=[]
    for i in polys:
        assert isinstance(i,Polynomial)
        
    univariate_polys = [ poly for poly in polys if poly.isUnivariate() ]

    if globalsettings.getSetting("solvePolynomialEquationsLog"):
        poly_history=poly_history+'\n\n'+str(map(_printPoly,univariate_polys))+'\n'
    
    if univariate_polys:
        univariate_poly = univariate_polys[0]
        #    print univariate_poly
        variable_name = univariate_poly.variables()[0]
        if globalsettings.getSetting("solvePolynomialEquationsLog"):
            poly_history = poly_history + '\n\nSolving for %s\n' % variable_name

        try:
            sol = polynomialSolver(univariate_poly)
            if globalsettings.getSetting("solvePolynomialEquationsLog"):
                poly_history = poly_history+'\n'+str(sol)+'\n'
        except Exception as e:
            raise SolverException("Error in find_complex_roots when solving: " +
                                  str(univariate_poly) + " " + repr(e),
                                  poly_history)

        assert len(sol)==univariate_poly.degree()
        #if not len(sol)==1:
        #    if non_linear_equation_encountered:
        #        raise SolverException(
        #            "Encountered second non linear equation: " +
        #            str(univariate_poly),
        #            poly_history)
        #    
        #    non_linear_equation_encountered = True
    else:
        if free_dim == 0:
            raise SolverException("No univariate polynomial left",
                                  poly_history)
        else:
            univariate_poly = None
            variable_name = polys[-1].variables()[0]
            sol = [random_complex_modulos()]
            if globalsettings.getSetting("solvePolynomialEquationsLog"):
                poly_history += "In pick random solution for %s:\n %s\n\n" % (variable_name, sol)

            free_dim = free_dim - 1
        
    for value in sol:
        new_variable_dict = dict(variable_dict)
        new_variable_dict[variable_name] = value
        new_polys = [
            poly.substitute(
                { variable_name : Polynomial.constantPolynomial(value) })
            for poly in polys if not poly is univariate_poly]
        new_solutions = solvePolynomialEquations(
            new_polys,
            polynomialSolver = polynomialSolver,
            free_dim = free_dim,
            with_poly_history = with_poly_history,
            poly_history = poly_history,
            variable_dict = new_variable_dict)
        solutions = solutions + new_solutions
        
    return solutions
