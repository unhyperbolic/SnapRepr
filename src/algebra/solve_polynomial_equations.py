import numpy
import math
from algebra.polynomial import Polynomial
from algebra.pari import roots_of_polynomial
from algebra.pari import random_complex_modulos
from algebra.pari import number

class SolverException(Exception):
    def __init__(self, message, poly_hist):
        self.poly_hist = poly_hist
        self.msg = message
    def __str__(self):
        return self.msg + "\nHistory of polynomials:\n" + self.poly_hist

# fills free variables with random values

def solve_polynomial_equations(polys,
                               free_dim = 0,
                               with_poly_history = False,
                               poly_history="",
                               variable_dict={},
                               non_linear_equation_encountered=False):

    polys = [polynomial.convertCoefficients(lambda x:number(x)) for polynomial in polys]
    
    poly_history += '\n\n\n\n'+'\n'.join(map(str,polys))+'\n\n============\n'
    if not polys:
        assert free_dim == 0
        if with_poly_history:
            return [(variable_dict,poly_history)]
        else:
            return [variable_dict]
    solutions=[]
    for i in polys:
        assert isinstance(i,Polynomial)
        
    univariate_polys = [poly for poly in polys if poly.isUnivariate()]
    
    poly_history=poly_history+'\n\n'+str(univariate_polys)+'\n'
    
    if univariate_polys:
        univariate_poly = univariate_polys[0]
        #    print univariate_poly
        variable_name = univariate_poly.variables()[0]
        poly_history = poly_history + '\n\nSolving for %s\n' % variable_name

        try:
            sol = roots_of_polynomial(univariate_poly)
            poly_history = poly_history+'\n'+str(sol)+'\n'
        except Exception as e:
            raise SolverException("Error in find_complex_roots when solving: " +
                                  str(univariate_poly) + " " + repr(e),
                                  poly_history)

        assert len(sol)==univariate_poly.degree()
        #if not len(sol)==1:
        #    if non_linear_equation_encountered:
        #        raise SolverException(
        #            "Encountered second non linear equation: " +
        #            str(univariate_poly),
        #            poly_history)
        #    
        #    non_linear_equation_encountered = True
    else:
        if free_dim == 0:
            raise SolverException("No univariate polynomial left",
                                  poly_history)
        else:
            univariate_poly = None
            variable_name = polys[-1].variables()[0]
            sol = [random_complex_modulos()]
            poly_history += "In pick random solution for %s:\n %s\n\n" % (variable_name, sol)

            free_dim = free_dim - 1
        
    for value in sol:
        new_variable_dict = dict(variable_dict)
        new_variable_dict[variable_name] = value
        new_polys = [
            poly.substitute(
                {variable_name:Polynomial.constantPolynomial(value)})
            for poly in polys if not poly is univariate_poly]
        new_solutions = solve_polynomial_equations(new_polys, free_dim,
                                                   with_poly_history,
                                                   poly_history,
                                                   new_variable_dict)
        solutions = solutions + new_solutions
        
    return solutions
