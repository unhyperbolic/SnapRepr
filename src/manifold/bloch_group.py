#from cmath import exp,pi,phase
#import cmath
import re
import math
import algebra.polynomial
from algebra import mpmathFunctions
from algebra.pari import *
import globalsettings
import mpmath

def mpmathDilog(z):
    return mpmath.polylog(2,z)

def pariDilog(z):
    assert isinstance(z, mpmath.mpc)

    pariStr = 'myDilog = dilog((%s) + (%s) * I)' % (
        mpmath.nstr(z.real, mpmath.mp.dps),
        mpmath.nstr(z.imag, mpmath.mp.dps))

    pari_eval(pariStr)
    
    return mpmath.mpc(
        pari_eval('real(myDilog)').replace(' E','E'),
        pari_eval('imag(myDilog)').replace(' E','E'))

myDilog = pariDilog

class DegenerateSimplexException(Exception):
    pass
    
class PtolemyCochain:
    def __init__(self, sign, c01, c02, c03, c12, c13, c23, no_check = False,
                 ptolemy_index = None, tet_index = None):
        self.sign = sign
        self.c01 = c01
        self.c02 = c02
        self.c03 = c03
        self.c12 = c12
        self.c13 = c13
        self.c23 = c23

        self.ptolemy_index = ptolemy_index
        self.tet_index = tet_index

        if not no_check:
            self.checkConsistency()

    def checkConsistency(self):
        
        c01c23 = self.c01 * self.c23
        c02c13 = self.c02 * self.c13
        c03c12 = self.c03 * self.c12

        maxErr = globalsettings.getSetting("maximalError")

        if not ( (abs(- c03c12 - c01c23 + c02c13) < maxErr) or
                 (abs(- c03c12 + c01c23 + c02c13) < maxErr) or
                 (abs(  c03c12 - c01c23 + c02c13) < maxErr) or
                 (abs(  c03c12 + c01c23 + c02c13) < maxErr)):

            raise NumericalError(
                val=[ (- c03c12 - c01c23 + c02c13).abs(),
                      (- c03c12 + c01c23 + c02c13).abs(),
                      (c03c12 - c01c23 + c02c13).abs(),
                      (  c03c12 + c01c23 + c02c13).abs()],
                msg = "PtolemyCochain(%s,%s,%s,%s,%s,%s)" %
                ( self.c01, self.c02, self.c03,
                  self.c12, self.c13, self.c23))

def logOfSquare(c):
    return mpmath.log(c)

    csquare = c * c

    maxErr = globalsettings.getSetting("maximalError")

    if not ((csquare.real > 0 or
             abs(csquare.imag) > maxErr)):
        raise NumericalError(c, msg = "logOfSqaure near branch cut")

    return mpmath.log(csquare) / 2
 
class w_triple:
    def __init__(self, P, no_check = False):
        assert isinstance(P, PtolemyCochain)
        self.sign = P.sign

        self.ptolemy_index = P.ptolemy_index
        self.tet_index = P.tet_index

        self.w0 = logOfSquare(P.c03) + logOfSquare(P.c12) - logOfSquare(P.c02) - logOfSquare(P.c13)
        self.w1 = logOfSquare(P.c02) + logOfSquare(P.c13) - logOfSquare(P.c01) - logOfSquare(P.c23)

        if no_check:
            self.w2 = - (self.w0 + self.w1)
            return

        self.w2 = logOfSquare(P.c01) + logOfSquare(P.c23) - logOfSquare(P.c03) - logOfSquare(P.c12)

        w0_e = mpmath.exp(self.w0)
        w1_e = mpmath.exp(-self.w1)

        errs = [abs(1 - w0_e - w1_e),
                abs(1 + w0_e - w1_e),
                abs(1 - w0_e + w1_e),
                abs(1 + w0_e + w1_e)]

        if not (errs[0] < globalsettings.getSetting("maximalError") or
                errs[1] < globalsettings.getSetting("maximalError") or
                errs[2] < globalsettings.getSetting("maximalError") or
                errs[3] < globalsettings.getSetting("maximalError")):
            raise NumericalError(val = errs,
                                 msg = "w triple (%s,%s,%s) from %s %s %s %s %s %s" % 
                                 (self.w0, self.w1, self.w2, P.c01, P.c02, P.c03, P.c12, P.c13, P.c23))
            
        if not abs(self.w0+self.w1+self.w2) < globalsettings.getSetting("maximalError"):
            raise NumericalError(val = (self.w0+self.w1+self.w2).abs(),
                                 msg = "w triple (%s,%s,%s) not sum to 0 from %s %s %s %s %s %s" % 
                                 (self.w0, self.w1, self.w2, P.c01, P.c02, P.c03, P.c12, P.c13, P.c23))

def mpmathRoundToInt(z):
    return int(mpmath.nint(z.real))

class zpq_triple:
    def __init__(self, w, no_check = False):
        assert isinstance(w, w_triple)
        self.sign = w.sign

        self.ptolemy_index = w.ptolemy_index
        self.tet_index = w.tet_index
        
        def pqCandidates(z, w = w):

            PiI = mpmath.pi * 1j

            if abs(1 - z) < globalsettings.getSetting("maximalError"):
                return (z, 0, 0), 1

            p = mpmathRoundToInt( (w.w0 - mpmath.log(  z)) / PiI)
            q = mpmathRoundToInt( (w.w1 + mpmath.log(1-z)) / PiI)
            err = abs(w.w2 + mpmath.log(z) - mpmath.log(1-z) + p * PiI + q * PiI)

            return (z, p, q), err

        candidate1, err1 = pqCandidates(z =   mpmath.exp(w.w0))
        candidate2, err2 = pqCandidates(z = - mpmath.exp(w.w0))

        if err1 < err2:
            candidate, err = candidate1, err1
        else:
            candidate, err = candidate2, err2

        if err > globalsettings.getSetting("maximalError"):
            raise NumericalError(val = err,
                                 msg = "err1 %s %s" % (candidate, err))

        self.z, self.p, self.q = candidate

        self.zp = 1 / (1 - self.z)
        self.zpp = 1 - 1 / self.z

        return

        if abs(w.w0) < globalsettings.getSetting("maximalError"):
            err1 = 1
        else:
            z1 = mpmath.exp(w.w0)
            p1 = int( ((w.w0 - mpmath.log(z1)) / PiI).real )
            q1 = int( ((w.w1 + mpmath.log(1 - z1)) / PiI).real )
            err1 = abs((w.w2 + mpmath.log(z1) - mpmath.log(1 - z1) + p1 * PiI + q1 * PiI))

        z2 = - mpmath.exp(w.w0)
        p2 = int( ((w.w0 - mpmath.log(z2)) / PiI).real )
        q2 = int( ((w.w1 + mpmath.log(1 - z2)) / PiI).real )
        err2 = abs((w.w2 + mpmath.log(z2) - mpmath.log(1 - z2) + p2 * PiI + q2 * PiI))

        if err1 < err2:
            if not err1 < globalsettings.getSetting("maximalError"):
                raise NumericalError(val = err1,
                                     msg = "err1 %s %s %s %s %s %s" % (
                        z1, p1, q1, w.w0, w.w1, w.w2))
            self.z = z1
            self.p = p1
            self.q = q1
        else:
            if not err2 < globalsettings.getSetting("maximalError"):
                raise NumericalError(val = err2,
                                     msg = "err1 %s %s %s %s %s %s" % (
                        z2, p2, q2, w.w0, w.w1, w.w2))
            self.z = z2
            self.p = p2
            self.q = q2

    def __repr__(self):
        return "[%s;%d,%d]" % (str(self.z), self.p, self.q)

    def cross_ratios(self):
        ext_str = ""
        for i in self.ptolemy_index:
            if i > 9:
                ext_str += "_"
            ext_str += "%d" %i
        ext_str += "_%d" % self.tet_index

        return {
            "z_" + ext_str : self.z,
            "zp_" + ext_str : self.zp,
            "zpp_" + ext_str : self.zpp
            }

    def L_function(self):
        p = self.p
        q = self.q
        z = self.z
        PiI = mpmath.pi * 1j

        val= (
              myDilog(z)
            + (mpmath.log(z) + p * PiI) * ( mpmath.log(1 - z) + q * PiI) / 2
            - mpmath.pi ** 2 / 6)

        if self.sign == -1:
            return -val
        else:
            return val
    
    def volume(self):
        z = self.z
        val = (  mpmath.arg(1 - z) * mpmath.log(abs(z))
               + myDilog(z).imag)
        if self.sign == -1:
            return -val
        else:
            return val
        
