def identity_matrix(size,the_type):
    return matrix([[the_type.one() if i==j else the_type.zero()
                    for j in range(size)]
                   for i in range(size)],the_type)

class matrix(object):
    """
    Represents a matrix
    
    >>> matrix([[1,2,3],[2,3,5]],int)
    matrix([
            [ 1, 2, 3],
            [ 2, 3, 5]
           ], the_type = int)

    The coefficients are converted to field (through field constructor)
    >>> from algebra.field_p import field_p
    >>> field = field_p(3)
    >>> matrix([[2,1,2],[0,2,1]], field)
    matrix([
            [ 2, 1, 2],
            [ 0, 2, 1]
           ], the_type = field_p(3))

    >>> from algebra.polynomial import polynomial
    >>> matrix([['x*x','x'],['a','b']],polynomial)
    matrix([
            [ 'x^2', 'x'],
            [ 'a', 'b']
           ], the_type = polynomial)
    >>> matrix([[polynomial('x'),polynomial('a*a')]],polynomial)
    matrix([
            [ 'x', 'a^2']
           ], the_type = polynomial)
    
    """
    def _add_rows(x):
        return reduce(lambda x,y:map(lambda x,y:x+y,x,y),x)
    def _mult(m,n):
        def mult_vec_mat(x,m=n): 
            return matrix._add_rows(
                   map(
                         lambda x,y:map(lambda x,y=y:x*y,x),
                         m,x
                      ) 
                                   )
        return map(mult_vec_mat,m)
    _mult = staticmethod(_mult)
    _add_rows = staticmethod(_add_rows)

    def __init__(self,values,the_type=None):
        
        if isinstance(values,matrix):
            if the_type==None:
                self.values=[[x for x in y] for y in values.values]
                self.the_type=values.the_type
            else:
                self.values=[[the_type(x) for x in y] for y in values.values]
                self.the_type=the_type
        else:
            assert isinstance(values,list)
            # assert len(values) >= 1
            if len(values) > 0:
                assert isinstance(values[0],list)
            # assert len(values[0]) >= 1
            if the_type==None:
                the_type=type(values[0][0])
            self.the_type=the_type
            for row in values:
                assert len(row)==len(values[0])
            self.values=[[the_type(x) for x in y] for y in values]
            
    def __mul__(self,other):
        if isinstance(other,matrix):
            assert self.no_columns() == other.no_rows()            
            return matrix(matrix._mult(self.values,other.values),self.the_type)
        if isinstance(other,list):
            assert self.no_columns() == len(other)
            if len(other)== 0:
                return []
            def inner_product(v1,v2=other):
                return reduce(lambda x, y: x+y,
                              map(lambda x, y: x*y, v1,v2))
            return map(inner_product,self.values)
        return other.__matrix_action__(self)
    def __add__(self,other):
        new_values=[[None for x in y] for y in self.values]
        for i in range(len(self.values)):
            for j in range(len(self.values[0])):
                new_values[i][j]=self.values[i][j]+other.values[i][j]
        return matrix(new_values)
    def __sub__(self,other):
        new_values=[[None for x in y] for y in self.values]
        for i in range(len(self.values)):
            for j in range(len(self.values[0])):
                new_values[i][j]=self.values[i][j]-other.values[i][j]
        return matrix(new_values)
        
    def __eq__(self,other):
        if other==None:
            return False
        for i1,i2 in zip(self.values,other.values):
            for j1,j2 in zip(i1,i2):
                if not j1==j2:
                    return False
        return True

    def __repr__(self):
        try:
            type_str = self.the_type.type_name()
        except:
            if self.the_type == int:
                type_str = 'int'
            else:
                type_str = '%s' % self.the_type

        def convert_row(r):
            try:
                l = [repr(v.constructor_argument()) for v in r]
            except:
                l = [repr(v) for v in r]
                
            return '[ ' + ', '.join(l) + ']'

        return (
                 "matrix(["
            +  "\n        "
            + ",\n        ".join(map(convert_row,self.values))
            +  "\n       ], the_type = %s)" % type_str)

        return convert_row(self.values[0])

        return repr(val.constructor_argument())
                
        return ("matrix([\n        %s\n       ],the_type=%s)" %
                    (',\n        '.join([
                        ['[ %s ]' % ', '.join([repr(val.constructor_argument())
                                               for val in row])]
                         for row in self.values
                        ]),
                     self.type_str))


        try:
            return ("matrix([\n        %s\n       ],the_type=%s)" %
                    (',\n        '.join([
                        ['[ %s ]' % ', '.join([val.constructor_argument()
                                               for val in row])]
                         for row in self.values
                        ]),
                     self.type_str))
        except:
            return "matrix([\n        %s\n      ],the_type=%s)" % (',\n        '.join([repr(row) for row in self.values]),type_str)

    def determinant(self):
        v=self.values
        if len(v)==2:
            return v[0][0]*v[1][1]-v[1][0]*v[0][1]
        if len(v)==3:
            return ( v[0][0]*v[1][1]*v[2][2]
                    +v[0][1]*v[1][2]*v[2][0]
                    +v[0][2]*v[1][0]*v[2][1]
                    -v[0][2]*v[1][1]*v[2][0]
                    -v[0][0]*v[1][2]*v[2][1]
                    -v[0][1]*v[1][0]*v[2][2])
        raise Exception, "Not implemented"

    def transpose(self):
        if len(self.values) == 0:
            return matrix([[]],
                          self.the_type)
        return matrix([[self.values[j][i] for j in range(len(self.values))]
                       for i in range(len(self.values[0]))],
                      self.the_type)

    def swap_columns(self,i,j):
        for r in self.values:
            t=r[i]
            r[i]=r[j]
            r[j]=t

    def swap_rows(self,i,j):
        t=self.values[i]
        self.values[i]=self.values[j]
        self.values[j]=t

    def add_to_column_p_times_column(self,j,p,i):
        for r in self.values:
            r[j]=r[j]+p*r[i]

    def add_to_row_p_times_row(self,j,p,i):
        self.values[j]=[self.values[j][t]+p*self.values[i][t]
                        for t in range(len(self.values[i]))]

    def multiply_column(self,i,p):
        for t in range(len(self.values)):
            self.values[t][i]=p*self.values[t][i]

    def multiply_row(self,i,p):
        self.values[i]=[p*t for t in self.values[i]]

    def no_rows(self):
        return len(self.values)

    def no_columns(self):
        return len(self.values[0])

class complex2x2_matrix(matrix):
    def __init__(self,a,b=None,c=None,d=None):
       if b==None and c==None and d==None:
            self.values=a
       else:
            self.values=[[a,b],[c,d]]
    def __mul__(self,other):
        if isinstance(other,complex2x2_matrix):
            return complex2x2_matrix(
                     super(complex2x2_matrix,self).__mul__(other).values
                  )
        if isinstance(other,complex):
            return ((self.values[0][0]*other+self.values[0][1])/
                    (self.values[1][0]*other+self.values[1][1]))
        return other.__matrix_action__(self)
    def identity():
        return complex2x2_matrix(complex(1,0),complex(0,0),
                                  complex(0,0),complex(1,0))
    identity=staticmethod(identity)
      
