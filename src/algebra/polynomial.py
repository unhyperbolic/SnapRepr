import re
import operator
#from fractions import Fraction

#from algebra.pari import *

#mathTypes

#typesafeOperatorPolicy(a, b, op):

#getImaginaryUnit

def parseIntCoefficient(s):
    coeff, rest = re.match('([0-9]*)(.*)',s).groups()

    if coeff:
        coeff = int(coeff)
    else:
        coeff = None
    return coeff, rest

_operators = {
    '+' : operator.add,
    '-' : operator.sub,
    '*' : operator.mul,
    '^' : operator.pow
    }

def _applyOperator(op, l, r):
    return _operators[op](l,r)

_operatorPrecedence = {
    None : 0,
    '+' : 1,
    '-' : 1,
    '*' : 2,
    '^' : 3
    }
                
def _parsePolynomial(s, parseCoefficient = parseIntCoefficient):

    operandStack = [Polynomial()]
    operatorStack = []

    def debugPrint(s):
        print "=" * 75
        print "Remaining string : ", s
        print "Operator Stack   : ", operatorStack
        print "Operand Stack    : ", operandStack

    def evalPrecedingOperatorsOnStack(operator = None):
        while operatorStack:
            topOperator = operatorStack[-1]
            
            if topOperator == '(':
                return
            
            if (_operatorPrecedence[topOperator] <
                _operatorPrecedence[operator]):
                return
            
            topOperator = operatorStack.pop()
            r = operandStack.pop()
            l = operandStack.pop()

            operandStack.append(
                _applyOperator(topOperator, l, r))

    def processNextToken(s):
        s = s.lstrip()

        operand, rest = parseCoefficient(s)
        if operand:
            operandStack.append(operand)
            return rest

        nextChar, rest = s[0], s[1:]
        
        if nextChar in _operators.keys():
            operator = nextChar
            evalPrecedingOperatorsOnStack(operator)
            operatorStack.append(operator)
            return rest

        if nextChar in '()':
            parenthesis = nextChar
            if parenthesis == '(':
                operatorStack.append('(')
            else:
                evalPrecedingOperatorsOnStack()
                assert operatorStack.pop() == '('
            return rest

        raise "WTF?"

    s = s.strip()
    while s:
        # debugPrint(s)
        s = processNextToken(s)

    # debugPrint(s)
    evalPrecedingOperatorsOnStack(None)

    # debugPrint(s)

    assert not operatorStack
    assert (len(operandStack) == 1
            or (
                len(operandStack) == 2 and
                operandStack[0] == Polynomial())) 

    return operandStack[-1]

# This implements our type policy for multiplication
# Allowed:
#      multiplying two objects of same type
#      multiplying by an integer (from either left or right)
#
# If multiplying with an integer, the integer is cast to the type
# of the other factor.

def _typesafeMultiply(a, b):
    if type(a) == type(b):
        return a * b
    if type(a) == int:
        return type(b)(a) * b
    if type(b) == int:
        return a * type(a)(b)
    raise Exception, "Cannot multiply two different types"

def _combineDicts(listOfDicts, combineFunction):
    result = {}
    for aDict in listOfDicts:
        for k, v in a.Dict.items():
            if result.has_key(k):
                result[k] = combineFunction(result[k], v)
            else:
                result[k] = v
    return result

def _dictToOrderedTupleOfPairs(d):
    l = d.items()
    l.sort(key = lambda x:x[0])
    return tuple(l)

class Monomial(object):

    # Given a string as variable name, constructs the corresponding monomial
    @classmethod
    def fromVariableName(var):
        return Monomial(1, ((var, 1),))

    # Constructor
    def __init__(self, coefficient, vars):
        """
        >>> Monomial(2, (('a', 2), ('b', 3)))
        """

        self._coefficient = _coefficient

        if isinstance(vars, dict):
            self._vars = _dictToOrderedTupleOfPairs(vars)
        else:
            assert isinstance(vars, tuple)
            self._vars = vars

    # Returns the coefficient
    def getCoefficient(self):
        return self._coefficient

    # Returns the type of the coefficient
    def coefficientType(self):
        return type(self._coefficient)

    # Returns a tuple of pairs (variableName, Exponent)
    def getVars(self):
        return self._vars

    # Returns the list of variables
    def variables(self):
        return [var[0] for var in self._vars if var[1] > 0]

    # Returns the degree of the monomial
    def degree(self):
        return sum([var[1] for var in self._vars])
        
    # Multiply two monomials
    def __mul__(self, other):
        
        assert isinstance(other, Monomial)

        # Determine the coefficient
        coefficient = _typesafeMultiply(self._coefficient, other._coefficient)

        # Compute the variables
        varDict = _combineDicts([self._vars,other._vars],
                                lambda x,y: x + y)

        return Monomial(coefficient, varDict)

    # Check whether two monomials are equal
    def __eq__(self,other):

        assert isinstance(other, monomial)

        return (
            self._coefficient == other._coefficient and
            self._vars == other._vars)

class Polynomial(object):
    
    def __init__(self, monomials = ()):

        # combine monomials with the same variables and exponents
        # and bring them into canonical order

        assert isinstance(monomials, tuple)

        # create for each monomial a dictionary
        # with key being the variables and exponents
        # and value being the coefficient

        listOfVarsCoeffDicts = [
            { monomial.getVars() : monomial.getCoeffcient() }
            for monomial in monomials]

        # combine the dictionaries using sum
        combinedVarsCoeffDict = _combineDicts(listOfVarsCoeffDicts,
                                              lambda x, y: x + y)

        # turn dictionary into a list of pairs (vars, coefficient)
        # in canonical order
        orderedTupleOfVarsCoeffPairs = _dictToOrderedTupleOfPairs(
            combinedVarsCoeffDict)

        # turn pairs into monomials, skip trivial monomials
        combinedMonomials = [
            Monomial(coefficient, vars)
            for vars, coefficient in orderedTupleOfVarsCoeffPairs
            if not coefficient == 0]

        # convert to tuple
        self._monomials = tuple(combinedMonomials)

    def __eq__(self, other):
        return self._monomials == other._monomials

    def __add__(self, other):
        assert isinstance(other, Polynomial)
        return Polynomial(self._monomials + other._monomials)

    def __mul__(self, other):
        monomials = []
        
        for m in self._monomials:
            for n in self._monomials:
                monomials.append(m * n)
                
        return Polynomial(monomials)


class OldPolynomial(object): ### allows multiplication with integers
    def __init__(self,terms=[]):
        """
        >>> p1 = Polynomial('3 * t * t + t ^ 6 + x * t * y')
        >>> p2 = Polynomial('t * x * y + 3 * t^2 + t^6')
        >>> p1 == p2
        True
        >>> str(p1 + p2)
        '2 * t * x * y + 6 * t^2 + 2 * t^6'
        >>> str(p1 - p2)
        ''
        >>> str(p1  * p2)
        't^2 * x^2 * y^2 + 6 * t^3 * x * y + 9 * t^4 + 2 * t^7 * x * y + 6 * t^8 + t^12'
        >>> p1 == p1 ** 2
        False
        >>> p3 = Polynomial('x+1')
        >>> p4 = p3 ** 3
        >>> str(p4)
        '1  + 3 * x + 3 * x^2 + x^3'
        >>> p5 = p4.substitute_numerical({'x':4.5})
        >>> str(p5)
        '166.375'
        >>> p5.is_constant()
        True
        >>> p5.get_constant()
        166.375
        >>> str(p4.substitute({'x':p4}))
        '8 + 36 * x + 90 * x^2 + 147 * x^3 + 171 * x^4 + 144 * x^5 + 87 * x^6 + 36 * x^7 + 9 * x^8 + x^9'
        >>> p1.variables()
        ['t', 'x', 'y']
        >>> p1.is_univariate()
        False
        >>> p4.is_univariate()
        True
        >>> p1.leading_coefficient()
        Traceback (most recent call last):
        ...
        AssertionError
        >>> p4.leading_coefficient()
        1
        >>> p6 = Polynomial('1+x^2')
        >>> str(p4 % p6)
        '- 2 + 2 * x'
        >>> str(Polynomial('4+3*x').make_monic())
        '(4/3) + x'
        >>> str(p5 + 2 * p3 + p2 * 3)
        '168.375 + 3 * t * x * y + 9 * t^2 + 3 * t^6 + 2 * x'
        """

        if isinstance(terms, Fraction):
            self.terms=[(terms,)]
        if isinstance(terms,int):
            self.terms=[(terms,)]        
        if isinstance(terms,Polynomial):
            self.terms=terms.terms
        if isinstance(terms,str):
            # does the expression have any +, * or ^ sign
            def eval_str(s):
                found_plus_minus=None
                found_times=None
                found_power=None
                s=s.strip()
                if not s:
                    return Polynomial()
                opening_parenthesis=0
                for i in range(len(s)):
                    if s[i]=='(':
                        opening_parenthesis+=1
                    if s[i]==')':
                        opening_parenthesis-=1
                    if not opening_parenthesis:
                        if s[i] in ['+','-']:

                            found_plus_minus=i
                        if s[i]=='*':
                            found_times=i
                        if s[i]=='^':
                            found_power=i
                if not found_plus_minus==None:
                    if s[found_plus_minus]=='+':
                        return eval_str(s[:found_plus_minus])+eval_str(s[found_plus_minus+1:])
                    if s[found_plus_minus]=='-':
                        return eval_str(s[:found_plus_minus])-eval_str(s[found_plus_minus+1:])
                if not found_times==None:
                    l1=eval_str(s[:found_times])
                    l2=eval_str(s[found_times+1:])
                    return l1*l2
                if not found_power==None:
                    return eval_str(s[:found_power])**int(s[found_power+1:])
                if s[0]=='(':
                    assert s[-1]==')'
                    return eval_str(s[1:-1])
                if s in ['i','I']:
                    return Polynomial([(1j,)])
                if re.match(r'[_A-Za-z][_A-Za-z0-9]*$',s):
                    return Polynomial([(1,(s,1))])
                if re.match(r'[0-9]+$',s):
                    return Polynomial([(int(s),)])
                if re.match(r'[0-9]+/[0-9]+$',s):
                    return Polynomial([(Fraction(s),)])
                r=re.match(r'([0-9]+\.[0-9]+)(j?)$',s)
                if r:
                    if r.group(2):
                        return Polynomial([(complex(0,float(r.group(1))),)])
                    else:
                        return Polynomial([(float(r.group(1)),)])
                    
                raise Exception, "while parsing Polynomial %s" % s

            self.terms=eval_str(terms).terms
        if isinstance(terms,list):
            self.terms=terms
            
        assert isinstance(self.terms,list)
        for term in self.terms:
            assert isinstance(term,tuple)
            for j in term[1:]:
                assert isinstance(j,tuple)
                assert isinstance(j[0],str)
                assert isinstance(j[1],int)
                assert j[1]>=0

    @classmethod
    def type_name(cls):
        return "Polynomial"

    def constructor_argument(self):
        return str(self)

    def __eq__(self,other):
        if isinstance(other,type(None)):
            return False
        if isinstance(other,Polynomial):
            return self.terms==other.terms
        if not self.is_constant():
            return False
        return self.get_constant()==other
        

    def __add__(self,other):
        if not isinstance(other, Polynomial):
            other = Polynomial(other)
        return Polynomial(self.terms+other.terms).simplify()

    def __radd__(self,other):
        return self + other

    def __pow__(self,other):
        if other==0:
            return Polynomial([(1,)])
        if other % 2:
            return self * self ** (other - 1)
        return (self*self) ** (other/2)
        
    def __mul__(self, other):
        if not isinstance(other, Polynomial):
            return self * Polynomial(other)
        
        assert isinstance(other, Polynomial)
        new_terms=[]
        for i in self.terms:
            for j in other.terms:
                new_terms.append((i[0]*j[0],)+i[1:]+j[1:])
        return Polynomial(new_terms).simplify()

    def __rmul__(self, other):
        return self * Polynomial(other)

    def __sub__(self,other):
        return self+other*Polynomial([(-1,)])

    def __repr__(self):
        return "Polynomial("+repr(self.terms)+")"

    def __str__(self):
        res=""
        for term in self.terms:
            def to_power(pair):
                if pair[1]==1:
                    return pair[0]
                else:
                    return pair[0]+"^"+str(pair[1])
            coeff=[]
            if term[0] in [+1,-1]:
                if term[0]==+1:
                    res = res + ' + '
                else:
                    res = res + ' - '
            else:
                if (isinstance(term[0],int) or
                    isinstance(term[0],long) or
                    isinstance(term[0],float) or
                    isinstance(term[0],Fraction)) and term[0]<0:
                    res = res + ' - '
                    coeff = [str(abs(term[0]))]
                elif isinstance(term[0],Polynomial):
                    res = res + ' + '
                    coeff = ['( ' + str(term[0]) + ' )']
                else:
                    res = res + ' + '
                    coeff = ['('+str(term[0])+')']
                    
            res = res + ' * '.join(coeff+map(to_power,term[1:]))
            if not (coeff or term[1:]):
                res = res + ' 1 '
        res=res.strip()
        if res and res[0]=='+':
            return res[1:].strip()
        return res

    def simplify(self):
        new_terms={}
        for term in self.terms:
            coeff=term[0]
            vars=term[1:]
            new_vars={}
            for var in vars:
                if new_vars.has_key(var[0]):
                    new_vars[var[0]]=new_vars[var[0]]+var[1]
                else:
                    new_vars[var[0]]=var[1]
            all_vars=new_vars.keys()
            all_vars.sort()
            new_vars=tuple(
                [(var,new_vars[var]) for var in all_vars if new_vars[var]])
            if new_terms.has_key(new_vars):
                new_terms[new_vars]=new_terms[new_vars]+coeff
            else:
                new_terms[new_vars]=coeff
        all_terms=new_terms.keys()
        all_terms.sort()
        return Polynomial([(new_terms[vars],)+vars for vars in all_terms if not new_terms[vars]==0])
    
    def substitute_numerical(self,var,value=None):
        if isinstance(var,dict):
            p=Polynomial(self)
            for var,value in var.items():
                p=p.substitute_numerical(var,value)
            return p
            
        def handle_vars(v,var=var,value=value):
            coeff=v[0]
            for the_var,the_exp in v[1:]:
                if the_var==var:
                    coeff = coeff*value**the_exp
            return (coeff,)+tuple([x for x in v[1:] if not x[0]==var])
        return Polynomial(map(handle_vars,self.terms)).simplify()

    def substitute(self,var1,var2=None):
        """
        Substitute can take two arguments or a dictionary
        It can substitute a variable by another variable or a Polynomial.

        >>> p = Polynomial("x^3+y^3")

        Subsitute x by y
        >>> str(p.substitute("x","y"))
        '2 * y^3'
        >>> str(p.substitute({"x":"y"}))
        '2 * y^3'
        >>> str(p.substitute("x",Polynomial("y")))
        '2 * y^3'
        >>> str(p.substitute({"x":Polynomial("y")}))
        '2 * y^3'

        Subsitute x by a+1 and y by b+1
        >>> str(p.substitute({"x":Polynomial("a+1"),
        ...                  "y":Polynomial("b+1")}))
        '2 + 3 * a + 3 * a^2 + a^3 + 3 * b + 3 * b^2 + b^3'
        """
        
        if isinstance(var1,dict):
            assert var2==None
            d=var1
        else:
            d={var1:var2}

        def new_poly(term,d=d):
            new_term=[term[0]]
            poly=Polynomial([(1,)])
            for var,exp in term[1:]:
                if d.has_key(var):
                    if isinstance(d[var],str):
                        new_term.append((d[var],exp))
                    elif isinstance(d[var],Polynomial):
                        poly=poly*(d[var]**exp)
                    else:
                        raise Exception, 'substitute: must be variable or Polynomial'
                else:
                    new_term.append((var,exp))
            if poly==Polynomial([(1,)]):
                return [tuple(new_term)]
            else:
                return (Polynomial([tuple(new_term)])*poly).terms
            
        res=[]
        for term in self.terms:
            res=res+new_poly(term)
        return Polynomial(res).simplify()
        raise TypeError, " in Polynomial substitute"
            
    def variables(self):
        vars=[]
        for term in self.terms:
            for var,exp in term[1:]:
                vars.append(var)
        vars=list(set(vars))
        vars.sort()
        return vars

    def is_linear(self):
        return self.degree()==1

    def get_linear_and_constant_term(self):
        assert self.is_linear()
        assert self.is_univariate()
        return (self.terms[1][0],self.terms[0][0])
        
    def is_univariate(self):
        return len(self.variables())==1

    def is_constant(self):
        return len(self.variables())==0

    def get_constant(self):
        assert self.is_constant()
        assert len(self.terms)<=1
        if self.terms:
            return self.terms[0][0]
        else:
            return 0

    def degree(self,variable=None):
        max=0
        for term in self.terms:
            deg=0
            for var,exp in term[1:]:
                if (not variable) or var==variable:
                    deg=deg+exp            
            if deg>max:
                max=deg
        return max

    def leading_coefficient(self):
        assert self.is_univariate()
        lead_coeff=None
        highest_deg=-1
        for term in self.terms:
            this_deg=0
            if len(term)==2:
                this_deg=term[1][1]
            if this_deg>highest_deg:
                highest_deg=this_deg
                lead_coeff=term[0]
        return lead_coeff
    
    def coefficients(self):
        assert self.is_univariate()
        coeffs=[0 for i in range(self.degree()+1)]
        for term in self.terms:
            if len(term)==1:
                coeffs[0]=term[0]
            else:
                coeffs[term[1][1]]=term[0]
        return coeffs

    def __mod__(self,other):
        assert isinstance(other,Polynomial)
        assert other.is_univariate()
        other=other.make_monic()
        deg=other.degree()
        var=other.variables()[0]
        highest_deg=self.degree(var)
        if highest_deg < deg:
            return self
        for term in self.terms:
            # found the term with highest power of var
            if (var,highest_deg) in term: 
                p=self+other*Polynomial([
                    (-term[0],(var,highest_deg-deg))])
                # force this term to be zero!
                return Polynomial(
                    [t for t in p.terms if not t[1:]==term[1:]]) % other
            
        return self
        
    def __hash__(self):
        return hash(tuple(self.terms))
