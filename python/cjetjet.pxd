from libcpp.string cimport string
from libcpp.vector cimport vector

cdef extern from "Vector.h" namespace "JETJet":
    cdef cppclass Vector:
        Vector()
        Vector(double, double, double, double)
        double px()
        double py()
        double pz()
        double E()
        double jetFunction()


cdef extern from "JetDefinition.h" namespace "JETJet":
    cdef cppclass JetDefinition:
        JetDefinition(double)
        string description()
    cdef cppclass EConeDefinition:
        EConeDefinition(double)
        string description()
    cdef cppclass EtConeDefinition:
        EtConeDefinition(double)
        string description()
    cdef cppclass EtAlphaConeDefinition:
        EtAlphaConeDefinition(double,double)
        string description()


cdef extern from "JetFinder.h" namespace "JETJet":
    cdef cppclass Jet:
        const Vector & jet()
        const vector[uint] & content()

    cdef cppclass JetFinder:
        JetFinder(vector[Vector] & particles, JetDefinition * jetDefinition)
        vector[Jet] jets()

