import copy
import operator
import globalsettings

def left_out_number(l):
    """
    >>> left_out_number([0,1,3])
    2
    """
    
    i = 0
    while i in l:
        i = i + 1
    return i

class edge(list):
    def __init__(self,tetrahedron,vert_0,vert_1):
        super(edge,self).__init__((tetrahedron,vert_0,vert_1))
    def __repr__(self):
        return "edge(tet=%d,%d,%d)" % (self[0],self[1],self[2])
    def __hash__(self):
        return hash(tuple(self))
    def flip(self):
        return edge(self[0],self[2],self[1])
    def tet(self):
        return self[0]
    def vert_0(self):
        return self[1]
    def vert_1(self):
        return self[2]

class face_class(list):
    """
    represents one face in the triangulation,
    the orientation is the one induced from the tetrahedron
    not the orientation induced from the ordering
    in particular, the boundary of a tetrahedron is represented by 4 face classes, all with positive coefficients
    (if the tetrahedron is not identified with itself and is tet1)
    a face has in the triangulation is obtained by identifying 
       face1 of tet1 and face2 of tet2
    orient = +1 means that the induced orientations from tet1 and tet2 agree
    orient = -1 means they do not
    if tet1 and tet2 have the same orientation, then orient = -1
    in this case, the boundary of tet1 and tet2 would not contain this face
    """
    
    def __init__(self,tet1,face1,tet2,face2,orient):
        super(face_class,self).__init__((tet1,face1,tet2,face2,orient))
    def __repr__(self):
        return "face_class(tet1=%d,face1=%d,tet2=%d,face2=%d,orient=%d)" % tuple(self)
                
    def __hash__(self):
        return hash(tuple(self))
    def tet1(self):
        return self[0]
    def face1(self):
        return self[1]
    def tet2(self):
        return self[2]
    def face2(self):
        return self[3]
    def orient(self):
        return self[4]

class permutation(list):
    def __init__(self,p=[0,1,2,3]):
        """
        permutation("0132")
        permutation([0,1,3,2])
        """
        if isinstance(p,str):
            super(permutation,self).__init__(map(int,p))
        else:
            super(permutation,self).__init__(p)
    def inverse(self):
        inv=[0 for i in range(len(self))]
        for i in range(len(self)):
            inv[self[i]]=i
        return permutation(inv)
    def is_odd(self):
        l = 0
        for i in [(0,1),(0,2),(0,3),(1,2),(1,3),(2,3)]:
            if self[i[0]] > self[i[1]]:
                l = l + 1
        return bool(l % 2)
    def is_even(self):
        return not self.is_odd()
    def sign(self):
        if self.is_odd():
            return -1
        else:
            return +1
    def __mul__(self,other):
        assert isinstance(other,permutation)
        assert len(self)==len(other)
        return permutation([self[other[x]] for x in range(len(self))])
    def __repr__(self):
        return "permutation(%s)" % super(permutation,self).__repr__()

class tetrahedron:
    def __init__(self,index,neighbor_index,gluing,
                 cusp_index=[-1,-1,-1,-1],positive_orientation=True,
                 peripheral_curves = None):
        self.index=index
        self.neighbor_index=neighbor_index
        self.cusp_index=cusp_index
        self.gluing=map(permutation,gluing)
        self.positive_orientation = bool(positive_orientation)
        self.edge_orientation=[[0 for y in range(4)] for x in range(4)]
        if peripheral_curves:
            self.peripheral_curves = peripheral_curves
        else:
            self.peripheral_curves = [[[[0 
                                         for x in range(4)]
                                        for y in range(4)]
                                       for z in range(2)]
                                      for a in range(2)]
        # marks whether edge from x and y is oriented positive or not

    def __repr__(self):
        return ("tetrahedron(index=%d,neighbor_index=%s,gluing=%s,"
                "cusp_index=%s,positive_orientation=%s)"
                % (self.index,self.neighbor_index,self.gluing,
                   self.cusp_index,self.positive_orientation))

    def is_edge_orientation_consistent_on_face(self,face):
        f=face+(face[0],) # make cyclic
        signs = set()
        for j in range(3):
            signs = signs | set([self.edge_orientation[f[j]][f[j+1]]])
        return (0 in signs) or (signs == set([+1,-1]))
            
    def is_edge_orientation_consistent(self):
        for i in [(1,2,3),(0,2,3),(0,1,3),(0,1,2)]: # for each face
            if not self.is_edge_orientation_consistent_on_face(i):
                return False
        return True

def read_triangulation_from_file(filename):
    """
    >>> dir = globalsettings.getSetting("testTriangulationsPath")
    >>> t=read_triangulation_from_file(dir + "/m003.trig")
    >>> len(t)
    2
    >>> t.num_or_cusps
    1
    >>> t.num_nonor_cusps
    0
    >>> t[0].gluing[1]
    permutation([2, 1, 0, 3])
    >>> t.cusp_structure()
    [0]
    >>> t.num_or_cusps
    1
    >>> t.num_nonor_cusps
    0
    >>> [x.positive_orientation for x in t]
    [True, True]
    >>> t.orient()
    >>> [x.positive_orientation for x in t]
    [True, True]
    >>> a=triangulation(t.to_SnapPea())
    >>> a[1].gluing[2]
    permutation([0, 3, 2, 1])
    >>> a[1].gluing[2].sign()
    -1
    >>> t.is_ordered()
    False
    >>> len(t.get_edge_classes())
    2
    >>> len(t.get_face_classes())
    4
    >>> t.find_orderings()
    []
    >>> t.check_consistency()
    >>> t[0].neighbor_index[0]=0
    >>> t.check_consistency()
    Traceback (most recent call last):
    ...
    AssertionError
        
    """
    
    
    return triangulation(open(filename,'r').read())
            
class triangulation:
    def reorder_orient(self):
        self.orient()
        perms = [
            permutation() if t.positive_orientation else permutation([0,1,3,2])
            for t in self.tet_list]
        self.reorder_tets(perms)
    
    def orient(self):

        oriented_tets = [set([])]

        def orient_tet(tet,
                       positive_orientation=True,
                       tets = self.tet_list,
                       oriented_tets = oriented_tets):

            if not tet.index in oriented_tets[0]:
                tet.positive_orientation = positive_orientation
                oriented_tets[0] = oriented_tets[0] | set([tet.index])
                for face_index in range(4):
                    orient_tet(tets[tet.neighbor_index[face_index]],
                               bool(positive_orientation ^ tet.gluing[face_index].is_even()))

        orient_tet(self.tet_list[0])
    
    def allTetsPositiveOrientation(self):
        self.orient()
        return reduce(
            operator.and_,
            [tet.positive_orientation for tet in self.tet_list],
            True)
        
    def cusp_structure(self):
        self.num_or_cusps=0
        self.num_nonor_cusps=0
        
        def one_vertex(tet_index,vert_index,cusp_index,tet_list=self.tet_list):
            if tet_list[tet_index].cusp_index[vert_index]==-1:
                tet_list[tet_index].cusp_index[vert_index]=cusp_index
                for face_index in range(4):
                    if not face_index==vert_index:
                        one_vertex(tet_list[tet_index].neighbor_index[face_index],
                                   tet_list[tet_index].gluing[face_index][vert_index],
                                   cusp_index)
        for tet in self.tet_list:
            for vert_index in range(4):
                tet.cusp_index[vert_index]=-1
                
        for tet in self.tet_list:
            for vert_index in range(4):
                if tet.cusp_index[vert_index]==-1:
                    one_vertex(tet.index,vert_index,self.num_or_cusps)
                    self.num_or_cusps += 1

        def get_euler(cusp_index,
                      tet_list=self.tet_list,
                      edge_classes=self.get_edge_classes(both_orientations=True)):
            verts=[]
            no_verts=0
            triangs=0
            for tet in tet_list:
                for vert_index in range(4):
                    if tet.cusp_index[vert_index]==cusp_index:
                        # for every vertex of every tetrahedron belonging to the cusp
                        triangs += 1
                        for vert2_index in range(4):
                            if not vert2_index==vert_index:
                                for edge_class in edge_classes:
                                    if edge(tet.index,
                                            vert_index,
                                            vert2_index) in edge_class:
                                        if not edge_class in verts:
                                            verts.append(edge_class)
                                            no_verts += 1
            assert triangs % 2 == 0, "odd number of triangles"
            return no_verts - (triangs/2)

        self.cusp_shape=[('torus',0.0,0.0) for x in range(self.num_or_cusps)]
        return [get_euler(x) for x in range(self.num_or_cusps)]

    def __len__(self):
        return len(self.tet_list)

    def __getitem__(self,k):
        return self.tet_list[k]

    def __init__(self,s=None):
        self.tet_list=[]

        if isinstance(s, triangulation):
            for k, v in s.__dict__.items():
                self.__dict__[k] = copy.deepcopy(v)
        elif isinstance(s,str):
            l=s.split('\n')
            self.comment_line=l[0].split()
            self.name=l[1]
            self.solution_type=l[2]
            self.orientation=l[3]
            self.cs=l[4]
            l=' '.join(l[5:])
            l=l.split()
            self.num_or_cusps=int(l[0])
            self.num_nonor_cusps=int(l[1])

            self.num_or_cusps = self.num_or_cusps + self.num_nonor_cusps
            self.num_nonor_cusps = 0

            self.cusp_shape=[]
            l=l[2:]
            for i in range(self.num_or_cusps+self.num_nonor_cusps):
                self.cusp_shape.append((l[0],float(l[1]),float(l[2])))
                l=l[3:]
            self.num_tets=int(l[0])
            l=l[1:]
            for i in range(self.num_tets):
                positive_orientation=True
                if (len(self.comment_line)>2+i and 
                    self.comment_line[1]=='orientations:' and 
                    self.comment_line[2+i]=='negative'):
                    positive_orientation=False
                if (len(self.comment_line)>3+i and 
                    self.comment_line[2]=='orientations:' and 
                    self.comment_line[3+i]=='negative'):
                    positive_orientation=False
                self.tet_list.append(
                    tetrahedron(index = i,
                                neighbor_index = map(int,l[0:4]),
                                gluing = l[4:8],
                                cusp_index = map(int,l[8:12]),
                                peripheral_curves = 
                                 [[[[int(l[12+16*merdLong+8*handedness+4*vert+face])
                                     for face in range(4)]
                                    for vert in range(4)]
                                   for handedness in range(2)]
                                  for merdLong in range(2)],
                                positive_orientation = positive_orientation))
                l=l[4*3+4*16+2:]
        else:
            raise ValueError, "Argument to constructor needs to be a string or a triangulation"

    def to_SnapPea(self):
        out = "% Triangulation orientations:"
        for i in self.tet_list:
            if i.positive_orientation:
                out = out + " positive"
            else:
                out = out + " negative"
        out = out + "\n"
        out = out + "%s\nnot_attempted 0.0000\n" % self.name
        out = out + "unknown_orientability\nCS_unknown\n\n"
        out = out + "%d %d\n" % (self.num_or_cusps,self.num_nonor_cusps)
        for c in self.cusp_shape:
            out = out + "      %s %f %f\n" % (c[0],c[1],c[2])
        out = out + "\n%d\n" % self.num_tets
        for tet in self.tet_list:
            out = out + "%5d%5d%5d%5d\n " % tuple(tet.neighbor_index)
            for perm in tet.gluing:
                out = out + " %d%d%d%d" % tuple(perm)
            out = out + "\n"
            out = out + "%5d%5d%5d%5d\n" % tuple(tet.cusp_index)
            out = out + 4 * (" "+16 * "  0"+"\n")
            out = out + "   0.1 0.1\n\n"
        return out

    def check_consistency(self):
        assert len(self.tet_list) == self.num_tets
        for tet_index in range(self.num_tets):
            for face_index in range(4):
                tet=self.tet_list[tet_index]
                new_tet_index=tet.neighbor_index[face_index]
                new_face_index=tet.gluing[face_index][face_index]
                new_tet=self.tet_list[new_tet_index]
                assert (
                    self.tet_list[new_tet_index].positive_orientation ==
                    tet.positive_orientation ^
                    tet.gluing[face_index].is_even())
                
                for vert_index in range(4):
                    new_vert_index=tet.gluing[face_index][vert_index]
                    assert (tet_index==
                            new_tet.neighbor_index[new_face_index])
                    assert (face_index==
                            new_tet.gluing[new_face_index][new_face_index])
                    assert (vert_index==
                            new_tet.gluing[new_face_index][new_vert_index])

    # vertex i of tet j will be named vertex perm[j][i]
    def reorder_tets(self,perms_or_edge_ordering):
        assert isinstance(perms_or_edge_ordering,list)
        if isinstance(perms_or_edge_ordering[0],permutation):
            assert len(perms_or_edge_ordering)==self.num_tets
            perms=perms_or_edge_ordering
        else:
            perms=self.turn_edge_ordering_into_relabel_tets(perms_or_edge_ordering)
        assert len(perms)==self.num_tets
        
        inv_perms=[perm.inverse() for perm in perms]

        def new_tet(tet,s=self,inv_perms=inv_perms,perms=perms):
            return tetrahedron(
                tet.index,
                neighbor_index=[
                    tet.neighbor_index[inv_perms[tet.index][i]]
                    for i in range(4)],
                gluing=[
                    permutation(
                        [
                        perms[tet.neighbor_index[inv_perms[tet.index][i]]]
                             [tet.gluing[inv_perms[tet.index][i]]
                                        [inv_perms[tet.index][j]]]
                        for j in range(4)])
                    for i in range(4)],
                cusp_index=[tet.cusp_index[inv_perms[tet.index][i]]
                            for i in range(4)],
                positive_orientation=(tet.positive_orientation
                                      ^ perms[tet.index].is_odd()))
        
        new_tet_list = [new_tet(tet) for tet in self.tet_list]
        self.tet_list=new_tet_list
                    
    def __repr__(self):
        return str(self.tet_list)
    
    def is_ordered(self):
        for i in self.tet_list:
            for j in range(4):
                g=i.gluing[j]
                g=g[:j]+g[j+1:]
                if g[0]>g[1]:
                    return False
                if g[1]>g[2]:
                    return False
        return True
    
    def get_next_edge(self,tet,v0,v1,face):
        next_tet=self.tet_list[tet].neighbor_index[face]
        next_vert0=self.tet_list[tet].gluing[face][v0]
        next_vert1=self.tet_list[tet].gluing[face][v1]
        next_face=(self.tet_list[tet].gluing[face]
                                            [left_out_number([v0,v1,face])])
        return next_tet,next_vert0,next_vert1,next_face

    def get_face_classes(self):
        face_classes = []
        processed_faces = set()
        for tet in self.tet_list:
            for f in range(4):
                if not (tet.index,f) in processed_faces:
                    #gluing = [(x-1 if x > tet.gluing[f][f] else x) for x in tet.gluing[f] if not x == tet.gluing[f][f]]
                    #if gluing in [ [0,1,2], [1,2,0], [2,0,1] ]:
                    #    orient = +1
                    #else:
                    #    orient = -1
                    orient = tet.gluing[f].sign() # does not work
                    face_classes.append(face_class(tet.index,
                                                   f,
                                                   tet.neighbor_index[f],
                                                   tet.gluing[f][f],
                                                   orient))
                    processed_faces.add((tet.index,f))
                    processed_faces.add((tet.neighbor_index[f],
                                         tet.gluing[f][f]))
        return face_classes
    
    def get_edge_class(self,tet,vert0,vert1):
        face = 0
        while face in [vert0,vert1]: face = face+1
        edge_class=[]
        while not edge(tet,vert0,vert1) in edge_class:
            edge_class.append(edge(tet,vert0,vert1))
            tet,vert0,vert1,face=self.get_next_edge(tet,vert0,vert1,face)
        return edge_class
    
    def get_edge_classes(self,both_orientations=False):
        """
        If the triangulation is ordered and both_orientations is False,
        then the edge classes will have the orientation induced by the ordering
        on each edge.
        I.e. the edge in each edge class will be 01, 02, 03, 12, 13, 23.
        """
        
        if both_orientations:
            the_edges=[(0,1),(0,2),(0,3),(1,2),(1,3),(2,3),
                       (1,0),(2,0),(3,0),(2,1),(3,1),(3,2)]
        else:
            the_edges=[(0,1),(0,2),(0,3),(1,2),(1,3),(2,3)]
       
        processed_edges=set()
        edge_classes=[]
        for i in range(self.num_tets):
            for j, k in the_edges:
                if not edge(i,j,k) in processed_edges:
                    e=self.get_edge_class(i,j,k)
                    e.sort()
                    edge_classes.append(e)
                    processed_edges=processed_edges | set(e)
                    if not both_orientations:
                        processed_edges=processed_edges | set([ed.flip() for ed in e])
        self.edge_classes=edge_classes
        return edge_classes

    def turn_edge_ordering_into_relabel_tets(self,edge_ordering):
        edge_ordering_tet=[[] for i in range(self.num_tets)]
        for i in edge_ordering:
            for j in i:
                assert isinstance(j,edge)
                edge_ordering_tet[j[0]].append(j[1:])

        def new_order_for_tet(l):
            k=map(lambda x:x[1],l)
            perm=[0,0,0,0]
            for i in k:
                perm[i]=perm[i]+1
            return permutation(perm)
                    
        res=[new_order_for_tet(edge_ordering_tet[i])
             for i in range(self.num_tets)]
        return res

    def is_relabel_tets_orientation_preserving(self,perms):
        return reduce(lambda x,y:x and y,[p.is_odd() for p in perms])

    def is_edge_order_orientation_preserving(self,edge_order):
        return self.is_relabel_tets_orientation_preserving(
            self.turn_edge_ordering_into_relabel_tets(edge_order))

    def find_orderings(self):
        """
        >>> dir = globalsettings.getSetting("testTriangulationsPath")
        >>> t=read_triangulation_from_file(dir + "/m053.trig")
        >>> len(t.find_orderings())
        4
        >>> t.is_edge_order_orientation_preserving(t.find_orderings()[0])
        False
        >>> t.is_edge_order_orientation_preserving(t.find_orderings()[1])
        False
        >>> t.is_edge_order_orientation_preserving(t.find_orderings()[2])
        False
        >>> t.reorder_tets(t.find_orderings()[0])
        >>> t.check_consistency()
        >>> [x.positive_orientation for x in t]
        [True, True, True, False]
        >>> t.is_ordered()
        True
        >>> t.orient()
        >>> [x.positive_orientation for x in t]
        [True, True, True, False]
        >>> t.reorder_orient()
        >>> [x.positive_orientation for x in t]
        [True, True, True, True]
        >>> t.check_consistency()
        """
        e_one=self.get_edge_classes()
        return self.orient_one_edge([],e_one)

    def orient_one_edge(self,edge_classes_oriented,edge_classes_left):
        res=[]
        if not edge_classes_left:
            return [edge_classes_oriented]
            # found a valid ordering of the edges
        # try orient edge_classes_left[0] one way first then the other
        for func in [lambda x:x,lambda x:x.flip()]:
            new_edge_class = map(func,edge_classes_left[0])
            # orient tetrahedron edges
            for i in new_edge_class:
                self.tet_list[i[0]].edge_orientation[i[1]][i[2]] = +1
                self.tet_list[i[0]].edge_orientation[i[2]][i[1]] = -1
            # check in if whether tetrahedra are consistently oriented
            consistent=True
            for i in set(map(lambda x:x[0],new_edge_class)):
                # get all tetrahedra involved
                consistent=(consistent and
                            self.tet_list[i].is_edge_orientation_consistent())
            if consistent:
                res=(res+
                     self.orient_one_edge(
                         edge_classes_oriented + [new_edge_class],
                         edge_classes_left[1:]))
            # disorient tetrahedron edges
            for i in new_edge_class:
                self.tet_list[i[0]].edge_orientation[i[1]][i[2]] = 0
                self.tet_list[i[0]].edge_orientation[i[2]][i[1]] = 0
        return res


