import fractions

lib_err = False

try:
    import _pari
except:
    lib_err = True

class NumericalError(Exception):
    def __init__(self, val, msg):
        self.val = val
        self.msg = msg
    def __str__(self):
        return "NumericalError(val = %s) : %s" % (self.val, self.msg)

def pari_eval(s):
    assert isinstance(s, str), "pari_eval requires string, but got %s" % s
    if lib_err:
        raise Exception, "pari not available"
    return _pari._pari_eval(s)

def pari_eval_bool(s):
    assert isinstance(s, str),\
        "pari_eval_bool requires string, but got %s" % s
    result = pari_eval(s)
    assert result in ("0","1"),\
        "expression in pari_eval_bool not boolean, argument was %s" % s
    return result == "1"

_precision = 0
_error = 0

def get_pari_precision():
    return _precision

def set_pari_allowed_error(precision):
    global _error
    _error = number(eval_this = "1 / 10^%d" % precision)

def get_pari_allowed_error():
    return _error

def get_pari_error():
    return _error

def set_pari_precision(precision):
    global _precision

    assert isinstance(precision, int) and precision > 9,\
        "set_pari_precision needs integer at least 10"

    pari_eval("default(realprecision,%d)" % precision)

    _precision = precision

class number(object):
    def __init__(self, val = None, eval_this = None):
        if val == None:
            if eval_this:
                self.val = pari_eval(eval_this)
            else:
                self.val = "0.0"
        elif isinstance(val,str):
            self.val = val
        elif isinstance(val,number):
            self.val = val.val
        elif (isinstance(val,int) or 
              isinstance(val,long) or 
              isinstance(val,fractions.Fraction)):
            self.val=str(val)
        else:
            raise Exception, "constructor to pari.number got %s" % val
    
    @staticmethod
    def do_op(a,op,b):
        if not isinstance(a,number):
            a=number(a)
        if not isinstance(b,number):
            b=number(b)
        return number(pari_eval("(%s)%s(%s)" % (a.val,op,b.val)))

    def __str__(self):
        return self.val.replace(" E","E")

    def dilog(self):
        return number(eval_this = "dilog(%s)" % self.val)
    def log(self):
        return number(eval_this = "log(%s)" % self.val)
    def exp(self):
        return number(eval_this = "exp(%s)" % self.val)
    def abs(self):
        return number(eval_this = "abs(%s)" % self.val)
    def real(self):
        return number(eval_this = "real(%s)" % self.val)
    def imag(self):
        return number(eval_this = "imag(%s)" % self.val)
    def arg(self):
        return number(eval_this = "arg(%s)" % self.val)

    def __int__(self):
        return int(pari_eval("round(%s)" % self.val))
        
    def __add__(self,other):
        return number.do_op(self,'+',other)
    def __sub__(self,other):
        return number.do_op(self,'-',other)
    def __neg__(self):
        return number.do_op(number("0"),'-',self)
    def __mul__(self,other):
        return number.do_op(self,'*',other)
    def __div__(self,other):
        return number.do_op(self,'/',other)
    def __pow__(self,other):
        return number.do_op(self,"^",other)
    def __radd__(self,other):
        return number.do_op(other,'+',self)
    def __rsub__(self,other):
        return number.do_op(other,'-',self)
    def __rmul__(self,other):
        return number.do_op(other,'*',self)
    def __rdiv__(self,other):
        return number.do_op(other,'/',self)
    def __mod__(self,other):
        return number.do_op(self,'%',other)
    def __abs__(self):
        return number(pari_eval("abs(%s)" % self.val))
        
    def __repr__(self):
        return "number('%s')" % self.val
    def __eq__(self,other):
        if not isinstance(other,number):
            other = number(other)
        return pari_eval_bool("(%s) == (%s)" % (self.val, other.val))
    def __cmp__(self,other):
        if not isinstance(other,number):
            other = number(other)
        if pari_eval_bool("(%s) < (%s)" % (self.val,other.val)):
            return -1
        if pari_eval_bool("(%s) > (%s)" % (self.val,other.val)):
            return +1
        return 0

    def pretty_print(self):
        real_part = number(eval_this = "real(%s)" % self.val)
        imag_part = number(eval_this = "imag(%s)" % self.val)
        
        sign_real = "+" if pari_eval_bool("%s > 0" % real_part.val) else "-"
        sign_imag = "+" if pari_eval_bool("%s > 0" % imag_part.val) else "-"

        real_part = real_part.abs()
        imag_part = imag_part.abs()

        non_zero_real = pari_eval_bool("%s > 10^(5-%s)" % (real_part.val, _precision))
        non_zero_imag = pari_eval_bool("%s > 10^(5-%s)" % (imag_part.val, _precision))

        if non_zero_real:
            real_part = sign_real + " " + real_part.val + " " * (5 + _precision - len(real_part.val))
        else:
            real_part = " " * (7 + _precision)

        if non_zero_imag:
            imag_part = sign_imag + " " + imag_part.val + " " * (5 + _precision - len(imag_part.val)) + " * I"
        else:
            imag_part = " "
        

        if not (non_zero_real or non_zero_imag):
            return " 0"

        else:
            return real_part + imag_part
        
def roots_of_polynomial(pol):
    assert pol.isUnivariate()
    s=pari_eval("polroots(%s)" % str(pol))
    assert s[0]=='['
    assert s[-2:]==']~'
    sol=[number(x) for x in s[1:-2].split(',')]
    assert len(sol)==pol.degree()
    return sol

def random():
    return number(eval_this="random()")

def random_complex_modulos():
    real_part = (100000 - (random() % 200000)) / 300000
    imag_part = 2 * Pi * I * ((random() % 100000) / 100000)
    r = (real_part+imag_part).exp()
    return r

Pi = number("Pi")
I = number("I")
PiI = number("Pi * I")

if not lib_err:
    set_pari_precision(15)
