import mpmath
from algebra.polynomial import Polynomial

def PolynomialSolver(polynomial):
    assert isinstance(polynomial, Polynomial)
    assert polynomial.coefficientType(mpmath.mpc) == mpmath.mpc
    return [mpmath.mpc(x) for x in mpmath.polyroots(polynomial.getCoefficients(mpmath.mpc))]

def complexModulusSquared(z):
    return z.real * z.real + z.imag * z.imag

def _dilogPowerSeries(z, rsquare):
    
    noTerms = mpmath.mp.prec * 2 / 1
    if rsquare < 0.25:
        noTerms = mpmath.mp.prec * 2 / 2
        if rsquare < 0.125:
            noTerms = mpmath.mp.prec * 2 / 3
            if rsquare < 0.0625:
                noTerms = mpmath.mp.prec * 2 / 4

    noTerms += 5
    terms = noTerms * [ 0 ]

    PowerOfZ = 1

    for i in range(1, noTerms + 1):
        PowerOfZ *= z
        terms[noTerms - i] = PowerOfZ / i ** 2

    return sum(terms)


def dilog(z):
    """
    

    >>> mpmath.mp.dps = 100

    Test power series

    >>> mpmath.nstr(dilog(mpmath.mpc("0.3", "0.1")),90)

    Test outside of disk of radius sqrt(2)

    >>> mpmath.nstr(dilog(mpmath.mpc("10.0", "0.5")),90)

    Test inside the anulus around 0 with radii sqrt(1/2) and sqrt(2)
    and outside of anulus around 1

    >>> mpmath.nstr(dilog(mpmath.mpc("-0.5", "0.8")),90)

    Test intersection of two anulii

    >>> mpmath.nstr(dilog(mpmath.mpc("0.5", "0.7")),90)
    >>> mpmath.nstr(dilog(mpmath.mpc("0.5", "-0.8")),90)
    >>> mpmath.nstr(dilog(mpmath.mpc("0.5", "0.86")),90)

    Test other random values

    >>> mpmath.nstr(dilog(mpmath.mpc("2.3", "4.5")),90)
    >>> mpmath.nstr(dilog(mpmath.mpc("-1.2", "3.4")),90)
    >>> mpmath.nstr(dilog(mpmath.mpc("-1.2", "-1.2")),90)
    >>> mpmath.nstr(dilog(mpmath.mpc("1.0", "0.3")),90)



    """

    # http://arxiv.org/pdf/hep-th/9408113v2 (page 31, L. Euler 1768)
    # http://www.math.ualberta.ca/~mlalin/dilogarithm.pdf

    #######
    # First cover cases where z is outside the anulus around 0 
    # and radii sqrt(1/2) and sqrt(2)

    # if |z| > sqrt(2) use
    #
    #    dilog(z) + dilog(1/z) = - pi**2 / 6 - log(-z)**2 /2
    #

    rsquare = complexModulusSquared(z)

    if rsquare > 2:
        return - mpmath.pi ** 2 / 6 - mpmath.log(-z) ** 2 / 2 - dilog(1 / z)

    # if |z| < sqrt(1/2) use the power series for the dilogarithm

    if rsquare < 0.5:
        return _dilogPowerSeries(z, rsquare)

    #######
    # Cover cases where z is outside the anulus around 1
    # and radii sqrt(1/2) and sqrt(2)

    #     dilog(z) + dilog(1-z) = pi**2 / 6 - ln(z) * ln(1-z)

    OneMinusZAbsSquared = complexModulusSquared(1 - z)
    if OneMinusZAbsSquared < 0.5 or OneMinusZAbsSquared > 2:
        return (
              mpmath.pi ** 2 / 6 
            - mpmath.log(z) * mpmath.log(1 - z) 
            - dilog(1 - z))

    ### Remaining case: z is in the intersection of the two anulii
    
    # dilog(z)  = 2 (dilog(sqrt(z)) + dilog(-sqrt(z)))

    sqrtZ = mpmath.sqrt(z)

    return 2 * (dilog(sqrtZ) + dilog(-sqrtZ))
