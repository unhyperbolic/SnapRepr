import operator

from manifold.triangulation import left_out_number,permutation
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


def produce_edge_equation_for_edge_class(edgeClass, trig):

    return NeumannZagierTypeEquation(
        reduce(operator.__mul__,
               [get_edge_parameter(e.vert_0(),e.vert_1(),e.vert_0(),trig.tet_list[e.tet()])
                for e in edgeClass]))

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
    return produce_face_equations(trig) + produce_edge_equations(trig) + produce_peripheral_equations(trig)

def produce_peripheral_equation_for_cusp(trig,
                                         curve_coeffs,
                                         cusp_index,
                                         dist):

    monomialLeft = Monomial.constantMonomial(1)
    monomialRight = Monomial.constantMonomial(1)

    for tet in trig.tet_list:
        for vNearCusp in range(4):
            if tet.cusp_index[vNearCusp] == cusp_index:
                for vAwayCusp in range(4):
                    if not vNearCusp == vAwayCusp:
                        edge_param = get_edge_parameter(
                            vNearCusp,vAwayCusp,
                            vNearCusp if dist == 0
                            else vAwayCusp,
                            tet)
                        for c in range(2):

                            def FLOW(a,b):
                                if (a<0 and b<0) or (a>0 and b>0):
                                    return 0
                                max_abs = min(abs(a), abs(b))
                                if a<0:
                                    return -max_abs
                                else:
                                    return max_abs

                            face = left_out_number([vNearCusp, vAwayCusp])
                            other_face = left_out_number([vNearCusp, vAwayCusp, face])
                            if permutation([vNearCusp, vAwayCusp, face, other_face]).is_odd:
                                face, other_face = other_face, face

                            #face, other_face such that with vNearCusp and vAwayCusp forms positively oriented tet

                            flow = FLOW(tet.peripheral_curves[c][0][vNearCusp][face],
                                        tet.peripheral_curves[c][0][vNearCusp][other_face])
                            
                            if not flow == 0:
                                print flow, edge_param

                            expo = curve_coeffs[c] * flow

                            if expo > 0:
                                monomialLeft *= edge_param ** expo
                            else:
                                monomialRight *= edge_param ** -expo


    return NeumannZagierTypeEquation(monomialLeft, monomialRight)


def produce_peripheral_equations_for_cusp(trig, curve_coeffs, cusp_index):
    
    return [
        produce_peripheral_equation_for_cusp(trig, curve_coeffs, cusp_index, dist) for dist in range(3-1)]

def produce_peripheral_equations_for_curve(trig, curve_coeffs):

    return reduce(operator.__add__,
                  [produce_peripheral_equations_for_cusp(trig, curve_coeffs, cusp_index)
                   for cusp_index in range(trig.num_or_cusps)])

def produce_peripheral_equations(trig):

    return produce_peripheral_equations_for_curve(trig,(1,0)) # just for meridian

def produce_all_equations_non_degenerate(trig):
    eqnsNZ = produce_all_equations(trig)
    eqnsPolys = [e.toPolynomial() for e in eqnsNZ]
    return eqnsPolys + [
        polynomialNonZeroCondition(eqnsPolys,
                                   var = 't',
                                   andNonOne = True)]

