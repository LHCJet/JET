from libcpp.string cimport string
from libcpp.vector cimport vector
cimport cjetjet as jj
import math


cdef class lorentz_vector:
    cdef public double px, py, pz, E

    def __init__(self, double px, double py, double pz, double E):
        self.px = px
        self.py = py
        self.pz = pz
        self.E = E

    def p_abs(self):
        return math.sqrt(self.px*self.px+self.py*self.py+self.pz*self.pz)

cdef class jet_definition:
    EConeDefinition = 1
    EtConeDefinition = 2
    EtAlphaConeDefinition = 3

cdef class jet:
    cdef public lorentz_vector vector
    cdef public list particles
    cdef public double jet_function

    def __init__(self, v, p, jf):
        self.vector = v
        self.particles = p
        self.jet_function = jf

cdef class jet_finder:
    cdef jj.JetFinder * jet_finder_ptr
    cdef jj.JetDefinition * jet_definition_ptr
    cdef list vl

    def __cinit__(self, list vl, jd, double beta, double alpha=1.0):
        cdef vector[jj.Vector] particles
        for v in vl:
            particles.push_back(jj.Vector(v.px, v.py, v.pz, v.p_abs()))

        if jd == jet_definition.EConeDefinition:
            self.jet_definition_ptr = <jj.JetDefinition * >new jj.EConeDefinition(beta)
        elif jd == jet_definition.EtConeDefinition:
            self.jet_definition_ptr = <jj.JetDefinition * >new jj.EtConeDefinition(beta)
        elif jd == jet_definition.EtAlphaConeDefinition:
            self.jet_definition_ptr = <jj.JetDefinition * >new jj.EtAlphaConeDefinition(alpha, beta)

        print(self.jet_definition_ptr.description())
        self.vl = vl
        self.jet_finder_ptr = new jj.JetFinder(particles, self.jet_definition_ptr)

    def jets(self):
        cjets = self.jet_finder_ptr.jets()
        jet_list = []
        for i in range(cjets.size()):
            v = cjets[i].jet()
            jet_v = lorentz_vector(v.px(), v.py(), v.pz(), v.E())
            indices = cjets[i].content()
            jet_p = []
            for j in indices:
                jet_p.append(self.vl[j])
            jet_list.append(jet(jet_v, jet_p, cjets[i].jet().jetFunction()))
        return jet_list
