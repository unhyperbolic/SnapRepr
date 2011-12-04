from cmath import exp,pi,phase
import cmath
import re
import math
import algebra.polynomial
from algebra.pari import *

class DegenerateSimplexException(Exception):
    pass
    
class Ptolemy_cochain:
    def __init__(self,sign, c01, c02, c03, c12, c13, c23, no_check = False):
        self.sign = sign
        self.c01 = c01
        self.c02 = c02
        self.c03 = c03
        self.c12 = c12
        self.c13 = c13
        self.c23 = c23

        if no_check:
            return
        
        c01c23 = c01 * c23
        c02c13 = c02 * c13
        c03c12 = c03 * c12

        if not ( ((- c03c12 - c01c23 + c02c13).abs() < get_pari_allowed_error()) or
                 ((- c03c12 + c01c23 + c02c13).abs() < get_pari_allowed_error()) or
                 ((  c03c12 - c01c23 + c02c13).abs() < get_pari_allowed_error()) or
                 ((  c03c12 + c01c23 + c02c13).abs() < get_pari_allowed_error())):

            raise NumericalError(
                val=[ (- c03c12 - c01c23 + c02c13).abs(),
                      (- c03c12 + c01c23 + c02c13).abs(),
                      (c03c12 - c01c23 + c02c13).abs(),
                      (  c03c12 + c01c23 + c02c13).abs()],
                msg = "Ptolemy_cochain(%s,%s,%s,%s,%s,%s)" %
                (c01, c02, c03, c12, c13, c23))

def log_square(c):
    csquare = c * c

    if not ((csquare.real() > 0 or
             csquare.imag().abs() > get_pari_allowed_error())):
        raise NumericalError(c, msg = "log_square near branch cut")

    return csquare.log() / 2
 
class w_triple:
    def __init__(self, P, no_check = False):
        assert isinstance(P, Ptolemy_cochain)
        self.sign = P.sign

        self.w0 = log_square(P.c03) + log_square(P.c12) - log_square(P.c02) - log_square(P.c13)
        self.w1 = log_square(P.c02) + log_square(P.c13) - log_square(P.c01) - log_square(P.c23)

        if no_check:
            self.w2 = - (self.w0 + self.w1)
            return

        self.w2 = log_square(P.c01) + log_square(P.c23) - log_square(P.c03) - log_square(P.c12)

        w0_e = self.w0.exp()
        w1_e = (-self.w1).exp()

        errs = [(1 - w0_e - w1_e).abs(),
                (1 + w0_e - w1_e).abs(),
                (1 - w0_e + w1_e).abs(),
                (1 + w0_e + w1_e).abs()]

        if not (errs[0] < get_pari_allowed_error() or
                errs[1] < get_pari_allowed_error() or
                errs[2] < get_pari_allowed_error() or
                errs[3] < get_pari_allowed_error()):
            raise NumericalError(val = errs,
                                 msg = "w triple (%s,%s,%s) from %s %s %s %s %s %s" % 
                                 (self.w0, self.w1, self.w2, P.c01, P.c02, P.c03, P.c12, P.c13, P.c23))
            
        if not (self.w0+self.w1+self.w2).abs() < get_pari_allowed_error():
            raise NumericalError(val = (self.w0+self.w1+self.w2).abs(),
                                 msg = "w triple (%s,%s,%s) not sum to 0 from %s %s %s %s %s %s" % 
                                 (self.w0, self.w1, self.w2, P.c01, P.c02, P.c03, P.c12, P.c13, P.c23))
            

class zpq_triple:
    def __init__(self, w, no_check = False):
        assert isinstance(w, w_triple)
        self.sign = w.sign

        self.sign = w.sign
        
        if w.w0.abs() < get_pari_allowed_error():
            err1 = 1
        else:
            z1 = w.w0.exp()
            p1 = int( ((w.w0 - z1.log()) / PiI).real() )
            q1 = int( ((w.w1 + (1 - z1).log()) / PiI).real() )
            err1 = ((w.w2 + z1.log() - (1 - z1).log() + p1 * PiI + q1 * PiI)).abs()

        z2 = - w.w0.exp()
        p2 = int( ((w.w0 - z2.log()) / PiI).real() )
        q2 = int( ((w.w1 + (1 - z2).log()) / PiI).real() )
        err2 = ((w.w2 + z2.log() - (1 - z2).log() + p2 * PiI + q2 * PiI)).abs()

        if err1 < err2:
            if not err1 < get_pari_allowed_error():
                raise NumericalError(val = err1,
                                     msg = "err1 %s %s %s %s %s %s" % (
                        z1, p1, q1, w.w0, w.w1, w.w2))
            self.z = z1
            self.p = p1
            self.q = q1
        else:
            if not err2 < get_pari_allowed_error():
                raise NumericalError(val = err1,
                                     msg = "err1 %s %s %s %s %s %s" % (
                        z2, p2, q2, w.w0, w.w1, w.w2))
            self.z = z2
            self.p = p2
            self.q = q2

    def L_function(self):
        val= (
              self.z.dilog()
            + (self.z.log() + self.p * PiI) * ( (1 - self.z).log() + self.q * PiI) / 2
            - Pi ** 2 / 6)

        if self.sign == -1:
            return -val
        else:
            return val
    
    def volume(self):
        val = (1 - self.z).arg() * self.z.abs().log() + self.z.dilog().imag()
        if self.sign == -1:
            return -val
        else:
            return val
        
