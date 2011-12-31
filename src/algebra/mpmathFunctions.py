import mpmath
from algebra.polynomial import Polynomial

def PolynomialSolver(polynomial):
    assert isinstance(polynomial, Polynomial)
    assert polynomial.coefficientType(mpmath.mpc) == mpmath.mpc
    return [
        mpmath.mpc(x) 
        for x in mpmath.polyroots(polynomial.getCoefficients(mpmath.mpc))]
