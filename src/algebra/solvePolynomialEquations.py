import numpy
import math
from algebra.polynomial import Polynomial, uncomparablePrintCoefficientMethod
from algebra import pari

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

def _filterPoly(polys, skip):
    return [poly for poly in polys if not poly == skip]

def exactSolutionsToNumerical(
        variableDict, nf, coeffConversion, polynomialSolver):

    # convert coefficients of nf and solve

    if nf: # if number field given
        nfSolutions = polynomialSolver(
            nf.convertCoefficients(coeffConversion))
    else:  # if no number field given, only solution is zero
        nfSolutions = [ coeffConversion(0) ]

    # convert all the polynomials in the variable dict
    variableDict = dict(
        [ (k, p.convertCoefficients(coeffConversion))
          for k, p in variableDict.items() ])

    def computeVariableDict(nfSolution,
                            variableDict = variableDict):
        def computeNumeric(
                p,
                substituteDict = {'x' : 
                                   Polynomial.constantPolynomial(
                                       nfSolution)}):
            
            c = p.substitute(substituteDict)
            
            assert c.isConstant()
            return c.getConstant()

        return dict(
            [(var, computeNumeric(val)) 
             for var, val in variableDict.items()])

    return [computeVariableDict(sol) for sol in nfSolutions]

def solvePolynomialEquationsExactly(polys, timeout = None):

    def conversionFunction(c):
        return Polynomial.constantPolynomial(c)

    polys = [ poly.convertCoefficients(conversionFunction) for poly in polys ]
    
    return _solvePolynomialEquationsExactly(polys,
                                            nf = None, variableDict = {},
                                            timeout = timeout)

def _transformVariableDict(variableDict, newExpressionForX, nf):
    return dict( [(k, v.substitute( {'x': newExpressionForX}) % nf)
                   for k, v in variableDict.items()] )

def _transformCoefficientsOfPolynomials(polys, newExpressionForX, nf):
    def substitute(p, newExpressionForX = newExpressionForX):
        return p.substitute( {'x': newExpressionForX} ) % nf
    return [ poly.convertCoefficients(substitute)
             for poly in polys]

def _setValueInPolynomials(polys, variable, value, nf = None):
    res = [
        poly.substitute(
            {variable:Polynomial.constantPolynomial(value)})
        for poly in polys]

    if nf:
        res = [poly.convertCoefficients(lambda x: x % nf) for poly in res]

    return res

def _solveExactlyOverNumberField(univariatePoly, nf, timeout):
    
    variable = univariatePoly.variables()[0]

    def convertXtoY(p):
        return p.substitute({'x' : Polynomial.fromVariableName('y')})

    univariatePoly = univariatePoly.convertCoefficients(convertXtoY)
    univariatePoly = univariatePoly.substitute(
        { variable : Polynomial.constantPolynomial(
            Polynomial.fromVariableName('x'))})

    if not nf:
        assert univariatePoly.isConstant()
        newSolution       = Polynomial.fromVariableName('x')
        newNf             = univariatePoly.getConstant()
        newExpressionForX = Polynomial.constantPolynomial(0)
    else:
        nf = convertXtoY(nf)

        pariStr = "PRIAVTEsEONF = rnfequation(nfinit(%s), %s, 1)" % (
            nf, univariatePoly)

        print pariStr
        print timeout
        r = pari.pari_eval(pariStr, timeout = timeout)
        # print r

        newNf              = Polynomial.parseFromMagma(pari.pari_eval(
                "PRIAVTEsEONF[1]", timeout = timeout))
        newExpressionForX  = Polynomial.parseFromMagma(pari.pari_eval(
                "PRIAVTEsEONF[2].pol", timeout = timeout))
        factor             = int(pari.pari_eval(
                "PRIAVTEsEONF[3]", timeout = timeout))
        newSolution = (
            Polynomial.fromVariableName('x')
            - Polynomial.constantPolynomial(factor) * newExpressionForX)

    return newSolution, newNf, newExpressionForX

def _convertToMonicNf(nf, timeout):

    pariStr = "PRIVATEconvertToMonicNf = nfinit(%s, 3)" % nf.printMagma()
    print pariStr
    print timeout
    r       = pari.pari_eval(pariStr, timeout = timeout)
    nf      = Polynomial.parseFromMagma(
        pari.pari_eval("PRIVATEconvertToMonicNf[1].pol", timeout = timeout))
    newExpressionForX = Polynomial.parseFromMagma(
        pari.pari_eval("PRIVATEconvertToMonicNf[2].pol", timeout = timeout))

    return nf, newExpressionForX

def _inverseOfConstantPolynomial(p):
    assert p.isConstant()
    constant = p.getConstant()
    invConstant = Fraction(1,1) / constant
    return Polynomial.constantPolynomial(invConstant)

def _solvePolynomialEquationsExactlyHandleLinear(
        polys,
        linearPoly, variable,
        nf, variableDict,
        timeout):
    
    factor, constant = linearPoly.getCoefficients()

    assert isinstance(factor, Polynomial) 
    assert isinstance(constant, Polynomial)

    newSolution = -constant * _inverseOfConstantPolynomial(factor)

    variableDict[variable] = newSolution
    
    return _solvePolynomialEquationsExactly(
        polys = _setValueInPolynomials(polys, variable, newSolution),
        nf = nf,
        variableDict = variableDict,
        timeout = timeout)

def _solvePolynomialEquationsExactlyHandleNonMonicNf(
        polys, nf, variableDict,
        timeout):

    nf, newExpressionForX = _convertToMonicNf(nf, timeout = timeout)

    return _solvePolynomialEquationsExactly(
        polys = _transformCoefficientsOfPolynomials(
            polys,
            newExpressionForX = newExpressionForX,
            nf = nf),
        nf = nf,    
        variableDict = _transformVariableDict(
            variableDict,
            newExpressionForX = newExpressionForX,
            nf = nf),
        timeout = timeout)
            
def _solvePolynomialEquationsExactlyHandleUnivariate(
        polys,
        univariatePoly, variable,
        nf, variableDict,
        timeout):

    newSolution, newNf, newExpressionForX = _solveExactlyOverNumberField(
            univariatePoly, nf, timeout = timeout)
    
    variableDict = _transformVariableDict(
        variableDict,
        newExpressionForX = newExpressionForX,
        nf = newNf)
    
    polys = _transformCoefficientsOfPolynomials(
        polys,
        newExpressionForX = newExpressionForX,
        nf = newNf)
    
    variableDict[variable] = newSolution
    polys = _setValueInPolynomials(polys, variable, newSolution, 
                                   nf = newNf)
    
    return _solvePolynomialEquationsExactly(
        polys,
        newNf,
        variableDict,
        timeout = timeout)
    
def _hasIntegralCoefficients(poly):
    coeffs = poly.getCoefficients()
    for coeff in coeffs:
        if not isinstance(coeff, int):
            assert isinstance(coeff, Fraction)
            if not coeff.denominator == 1:
                return False
    return True

def _solvePolynomialEquationsExactly(polys,
                                     nf = None, 
                                     variableDict = None,
                                     timeout = None):


    # nf is a polynomial in x encoding a number field

    # nf = x^2 + 1

    # variable contains the variables already bound as polynomials in the variable of the number field

    # a : x + 2
    if not polys:
        return variableDict, nf

    #print "=================== Enter _solvePolynomialEquationsExactly"

    #if nf:
    #    print "Number field: ", nf.printMagma()
    
    #print "Polynomials:"

    #for i in polys:
    #    print "         ", i.printMagma()

    linearPolys = [poly for poly in polys if poly.isLinear()]

    if linearPolys:
        linearPoly = linearPolys[0]

        return _solvePolynomialEquationsExactlyHandleLinear(
            polys = _filterPoly(polys, linearPoly),
            linearPoly = linearPoly,
            variable = linearPoly.variables()[0],
            nf = nf,
            variableDict = variableDict,
            timeout = timeout)

    univariatePolys = [poly for poly in polys if poly.isUnivariate()]
    if univariatePolys:

        if nf and not (nf.isMonic() and _hasIntegralCoefficients(nf)):
            return _solvePolynomialEquationsExactlyHandleNonMonicNf(
                polys, nf, variableDict, timeout = timeout)
        
        #if not nf:
        #    nf = (  Polynomial.fromVariableName('x') 
        #          - Polynomial.constantPolynomial(1))
            
        univariatePoly = univariatePolys[0]

        return _solvePolynomialEquationsExactlyHandleUnivariate(
            polys = _filterPoly(polys, univariatePoly),
            univariatePoly = univariatePoly,
            variable = univariatePoly.variables()[0],
            nf = nf,
            variableDict = variableDict,
            timeout = timeout)

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
