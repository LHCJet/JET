import jetjet
import sys

particles = []

for l in open(sys.argv[1]):
    p = [float(x) for x in l.split()]
    particles.append(jetjet.lorentz_vector(p[0], p[1], p[2], p[3]))

jf = jetjet.jet_finder(particles, jetjet.jet_definition.EtConeDefinition, 6.0)

jet_list = jf.jets()

for j in jet_list:
    v = j.vector
    print("%f, %f, %f, %f, %f" % (v.px, v.py, v.pz, v.E, j.jet_function))
