import mpmath
import globalsettings

def lessThanMaxErr(r):
    return abs(r) < globalsettings.getSetting("maximalError")

def _printComplexNumberAsFixed(z):
    assert isinstance(z, mpmath.mpc)

    r = printRealNumberAsFixed(z.real)

    if lessThanMaxErr(z.imag):
        return r

    i = printRealNumberAsFixed(z.imag)
    if not i[0] in '+-':
        i = '+' + i

    return r + ' ' + i + 'j'

def printComplexNumberAsFixed(z):
    return '(%s)' % _printComplexNumberAsFixed(z)

def printRealNumberAsFixed(r):

    assert isinstance(r, mpmath.mpf)
    
    if lessThanMaxErr(r):
        return '0'
    else:
        maxErrDigits = globalsettings.getSetting("maximalErrorDigits")
        s = mpmath.nstr(r,
                        maxErrDigits,
                        min_fixed = -mpmath.inf, max_fixed = mpmath.inf)
        a, b = s.split('.')

        # mpmath chops 1.2000000 to 1.2, so we are adding '0' to b

        alen = len(a)
        if a[0] in '+-':
            alen -= 1
        
        b += (maxErrDigits - alen + len(b)) * '0'

        # b should never be more than maximalErrorDigits

        b = b[:maxErrDigits - alen]

        return a + '.' + b
