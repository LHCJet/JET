#ifndef JETFINDER_H
#define JETFINDER_H

#include "Vector.h"
#include "JetCone.h"
#include "JetDefinition.h"

#include <vector>

namespace SlowJet {
class Jet
{
public:
    Jet(IndexList & c, const VectorList & p)
        : m_jet{}, m_content{c}, m_jetFunction(0) {
        PArray jetP = JetDefinition::instance()->sumP(c, p);
        m_jet = Vector(jetP[0], jetP[1], jetP[2], jetP[3]);
        m_jetFunction = JetDefinition::instance()->jetFunction(jetP);
    }
    Jet() : m_jet{}, m_content{}, m_jetFunction(0) {}
    ~Jet() {}
    const Vector & jet() const {return m_jet;}
    const IndexList & content() const {return m_content;}
    double jetFunction() const {return m_jetFunction;}

private:
    Vector m_jet;
    IndexList m_content;
    double m_jetFunction;

};

typedef std::vector<Jet> JetList;

class JetFinder
{
public:
    JetFinder(const VectorList & particles) : m_particles{particles}, m_cones{} {
        m_distances = NULL;
        m_particlesRemaining = new int[particles.size()];
        for (unsigned int i; i < particles.size(); i++) {
            m_particlesRemaining[i] = 1;
        }
    }
    ~JetFinder() {
    }
    JetList jets();
private:
    Jet findOneJet();
    Jet protoJet(VectorList subset);
    VectorList m_particles;
    JetConeList m_cones;
    double ** m_distances;
    int * m_particlesRemaining;
    int * m_conesRemaining;
};

} // namespace SlowJet;

#endif // JETFINDER_H
