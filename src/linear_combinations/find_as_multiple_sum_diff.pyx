import algebra.pari
from algebra.pari import number

# gets a list of list of dicts

CensusVols = []
cdef double fl_volumes[100000]
cdef int no_volumes

algebra.pari.set_pari_precision(120)

cdef extern from "math.h":
    double fabs(double x)

def initialize(CensVols):
    global CensusVols
    global no_volumes
    global fl_volumes

    assert len(CensVols) < 100000

    CensusVols = CensVols

    no_volumes = 0

    for d in CensVols:
        fl_volumes[no_volumes] = float(d[0]['Volume'])
        no_volumes += 1

def print_db():
    cdef int i
    for i in range(0,no_volumes):
        print fl_volumes[i]

def find_as_multiple(target, factor_range = 20, precision = number("1E-50")):
    cdef int i
    cdef double c
    cdef double c_target
    cdef int j
    cdef int c_factor_range

    list = []

    c_factor_range = factor_range + 1
    c_target = float(target.val)

    for i in range(1,c_factor_range):
        c = c_target / i
        for j in range(0, no_volumes):
            if fabs(fl_volumes[j] - c) < 1E-10:
                if (CensusVols[j][0]['Volume'] - target / number(i)).abs() < precision:
                    list.append([(i,CensusVols[j])])
    return list

def find_as_sum(target, factors = [1, 1], precision = number("1E-50")):
    cdef int start, end
    cdef double c_target
    cdef int f1, f2
    
    result = []
    f1 = factors[0]
    f2 = factors[1]
    c_target = float(target.val)

    start = 0
    end = no_volumes - 1
    
    while (start < no_volumes and end >= 0
           and f1 * fl_volumes[start] + 1E-9 < f2 * fl_volumes[end]):
        if fabs(f1 * fl_volumes[start] + f2 * fl_volumes[end]
                - c_target) < 1E-9:
            if (number(factors[0]) * CensusVols[start][0]['Volume'] +
                number(factors[1]) * CensusVols[end][0]['Volume'] -
                target).abs() < precision:
                result.append(
                    [(factors[0], CensusVols[start]),
                     (factors[1], CensusVols[end])])
        if f1 * fl_volumes[start] + f2 * fl_volumes[end] > c_target:
            end -= 1
        else:
            start += 1
    return result

def find_as_difference(target, factors = [1, 1], precision = number("1E-50")):
    cdef int start, end
    cdef double c_target
    cdef int f1, f2
    
    result = []
    f1 = factors[0]
    f2 = factors[1]
    c_target = float(target.val)

    start = 0
    end = 0
    
    while (start < no_volumes and end < no_volumes):
        if not start == end:
            if fabs(f1 * fl_volumes[start] - f2 * fl_volumes[end]
                    - c_target) < 1E-9:
                if (number(factors[0]) * CensusVols[start][0]['Volume'] -
                    number(factors[1]) * CensusVols[end][0]['Volume'] -
                    target).abs() < precision:
                    result.append(
                        [(factors[0], CensusVols[start]),
                         (-factors[1], CensusVols[end])])
        if f1 * fl_volumes[start] - f2 * fl_volumes[end] > c_target:
            end += 1
        else:
            start += 1
    return result


def find_as_linear_combination(target, factor_range = 20, precision = number("1E-50")):
    result = []
    for i in range(1,factor_range+1):
        for j in range(1,factor_range+1):
            result += find_as_sum(target, factors = [i,j], precision = precision)

    for i in range(1,factor_range+1):
        for j in range(1,factor_range+1):
            result += find_as_difference(target, factors = [i,j], precision = precision)
    return result
