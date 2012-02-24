import operator

from triangulation import left_out_number
from algebra.polynomial import Monomial, Polynomial
from manifold.slN import polynomialNonZeroCondition

class NeumannZagierTypeEquation(object):

    def __init__(self, left, right = None):
        self.left = left

        if right is None:
            right = Monomial.constantMonomial(1)

        self.right = right

    def __repr__(self):
        return ("NeumannZagierEquation(%s,%s)" % (
                str(self.left), str(self.right)))

    def toPolynomial(self):

        def processOneSide(origLeftMonomial):
            
            newLeft = Monomial.constantMonomial(origLeftMonomial.getCoefficient())
            newRight = Monomial.constantMonomial(1)

            for var, expo in origLeftMonomial.getVars():
                switch_sides = 'Inv' in var
                var = var.replace('Inv','')

                if 'zpprime' in var:
                    newLeftMultiply = (
                        Monomial.constantMonomial(-1) *
                        Monomial.fromVariableName(
                            var.replace('zpprime','OneMinusz')))
                    newRightMultiply = (
                        Monomial.fromVariableName(
                            var.replace('zpprime','z')))
                elif 'zprime' in var:
                    newLeftMultiply = (
                        Monomial.constantMonomial(1))
                    newRightMultiply = (
                        Monomial.fromVariableName(
                            var.replace('zprime','OneMinusz')))
                else:
                    newLeftMultiply = (
                        Monomial.fromVariableName(var))
                    newRightMultiply = (
                        Monomial.constantMonomial(1))
                    
                if switch_sides:
                    newLeft = newLeft * newRightMultiply**expo
                    newRight = newRight * newLeftMultiply**expo
                else:
                    newLeft = newLeft * newLeftMultiply**expo
                    newRight = newRight * newRightMultiply**expo
            
            return newLeft, newRight
        
        l1, r1 = processOneSide(self.left)
        l2, r2 = processOneSide(self.right)

        newLeft = Polynomial((l1*r2,))
        newRight = Polynomial((l2*r1,))


        substDict = dict(
            [
                (var, 
                 Polynomial(
                        (
                            Monomial.constantMonomial(1),
                            Monomial.constantMonomial(-1)*
                            Monomial.fromVariableName(
                                var.replace('OneMinusz','z')))))
                for var in newLeft.variables() + newRight.variables()
                if 'OneMinusz' in var])

        return (
              newLeft.substitute(substDict) 
            - newRight.substitute(substDict))

def get_edge_parameter(vert0, vert1, small_tet_vert, tet):
    
    if vert0 > vert1:
        nvert0, nvert1 = vert1, vert0
    else:
        nvert0, nvert1 = vert0, vert1

    if (nvert0, nvert1) in [(0,1),(2,3)]:
        b = "z"
    if (nvert0, nvert1) in [(0,3),(1,2)]:
        b = "zprime"
    if (nvert0, nvert1) in [(0,2),(1,3)]:
        b = "zpprime"

    if not tet.positive_orientation:
        b += "Inv"

    return Monomial.fromVariableName(b+"_%d_%d" % (small_tet_vert,tet.index))


def produce_edge_equation_for_edge_class(edgeClass,trig):

    return NeumannZagierTypeEquation(
        reduce(operator.__mul__,
               [get_edge_parameter(e.vert_0(),e.vert_1(),e.vert_0(),trig.tet_list[e.tet()]) for e in edgeClass]))

def produce_edge_equations(trig):

    return [produce_edge_equation_for_edge_class(e,trig) for e in trig.get_edge_classes(both_orientations=True)]

def produce_face_equation_for_face(tet, face):
    
    def one_mini_edge(small_tet_index, tet = tet, face = face):
        vert0 = left_out_number([small_tet_index,face])
        vert1 = left_out_number([small_tet_index,face,vert0])
        return get_edge_parameter(vert0, vert1, small_tet_index, tet)

    return reduce(operator.__mul__,
                  [ one_mini_edge(small_tet_index)
                    for small_tet_index in range(4) 
                    if not small_tet_index == face])

def produce_face_equation_for_face_class(faceClass, trig):
    
    return NeumannZagierTypeEquation(
        produce_face_equation_for_face(trig.tet_list[faceClass.tet1()],
                                       faceClass.face1()),
        produce_face_equation_for_face(trig.tet_list[faceClass.tet2()],
                                       faceClass.face2()))

def produce_face_equations(triangulation):

    return [produce_face_equation_for_face_class(f, triangulation) for f in triangulation.get_face_classes()]

def produce_all_equations(trig):
    return produce_face_equations(trig) + produce_edge_equations(trig)


def produce_all_equations_non_degenerate(trig):
    eqnsNZ = produce_all_equations(trig)
    eqnsPolys = [e.toPolynomial() for e in eqnsNZ]
    return eqnsPolys + [
        polynomialNonZeroCondition(eqnsPolys,
                                   var = 't',
                                   andNonOne = True)]
