def is_field_p_element(f, p = None):
    """
    >>> is_field_p_element(field_p(5)(3),5)
    True
    >>> is_field_p_element(field_p(5)(3),3)
    False
    >>> is_field_p_element(field_p(3)(2))
    True
    >>> is_field_p_element(field_p(3))
    False
    >>> is_field_p_element(5)
    False
    """
    
    if p == None:
        return isinstance(f, _field_p)
    
    return isinstance(f, _field_p) and f.get_p() == p

def is_field_p_class(f, p = None):
    """
    >>> is_field_p_class(field_p(3))
    True
    >>> is_field_p_class(field_p(3),3)
    True
    >>> is_field_p_class(field_p(3),5)
    False
    >>> is_field_p_class(field_p(3)(2))
    False
    >>> is_field_p_class(4)
    False
    """

    if p == None:
        try:
            return issubclass(f, _field_p)
        except:
            return False

    return f == field_p(p)

_field_cache = {}

def field_p(p):
    """
    creates a class representing Z/p, i.e. instances will be representing
    elements in Z/p.

    Z/5
    >>> field = field_p(5)

    3 * 4 in Z/5 is 2
    >>> field(3) * field(4)
    field_p(5)(2)

    As integer such that 0 <= x < p
    >>> int(field(3) * field(4))
    2
    """
    
    assert isinstance(p,int)
    assert p >= 2
    for i in range(2, p-1):
        assert p % i, "field_p(%d): %d is not prime" % (p,p)
    
    if not _field_cache.has_key(p):
        _field_cache[p] = type('field_%d' % p, (_field_p,), dict(p=p))

    return _field_cache[p]

_division_cache = {}

class _field_p(object):
    @classmethod
    def one(cls):
        return cls(1)

    @classmethod
    def zero(cls):
        return cls(0)

    @classmethod
    def elements(cls):
        return [cls(x) for x in range(cls.p)]

    @classmethod
    def type_name(cls):
        return "field_p(%d)" % cls.p
    
    def __init__(self,x):
        if isinstance(x,_field_p):
            self.p=x.p
            self.x=x.x
        else:
            assert isinstance(x,int)
            assert self.p>=2
            self.x= x % self.p

    def constructor_argument(self):
        return self.x

    def __add__(self,other):
        assert isinstance(other,_field_p)
        assert self.p==other.p
        return type(self)((self.x + other.x) % self.p)
    def __sub__(self,other):
        assert isinstance(other,_field_p)
        assert self.p==other.p
        return type(self)((self.x - other.x) % self.p)
    def __mul__(self,other):
        assert isinstance(other,_field_p)
        assert self.p==other.p
        return type(self)((self.x * other.x) % self.p)
    def __neg__(self):
        return type(self)((-self.x) % self.p)
    def __eq__(self,other):
        if isinstance(other,_field_p):
            assert self.p==other.p
            return (self.x-other.x) % self.p == 0
        else:
            assert isinstance(other,int)
            return (self.x-other) % self.p == 0
    def __repr__(self):
        return "field_p(%d)(%d)" % (self.p,self.x)
    def __str__(self):
        return "%d" % self.x
    def __hash__(self):
        return hash((self.x,self.p))

    def __div__(self,other):
        assert isinstance(other,_field_p)
        assert self.p==other.p
        assert not other.x % self.p == 0
        if not _division_cache.has_key(self.p):
            _division_cache[self.p]=range(self.p)
            for i in range(1,self.p):
                for j in range(1,self.p):
                    if i * j % self.p == 1:
                        _division_cache[self.p][i]=type(self)(j)
        return _division_cache[self.p][other.x % self.p]*self
        
    def get_p(self):
        return self.p

    def __int__(self):
        return self.x % self.p
