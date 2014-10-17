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
    Jet(VectorList & c)
        : m_jet{}, m_content{c}, m_jetFunction(0) {
        PArray jetP = JetDefinition::instance()->sumP(m_content);
        m_jet = Vector(jetP[0], jetP[1], jetP[2], jetP[3]);
        m_jetFunction = JetDefinition::instance()->jetFunction(jetP);
    }

    ~Jet() {}
    const Vector & jet() const {return m_jet;}
    const VectorList & content() const {return m_content;}
    double jetFunction() const {return m_jetFunction;}

private:
    Vector m_jet;
    VectorList m_content;
    double m_jetFunction;

};

typedef std::vector<Jet> JetList;

class JetFinder
{
public:
    JetFinder(const VectorList & particles) : m_particles{particles}, m_cones{} {}
    ~JetFinder() {}
    JetList jets();
private:
    Jet findOneJet();
    Jet protoJet(VectorList subset);
    VectorList m_particles;
    JetConeList m_cones;
};

} // namespace SlowJet;

#endif // JETFINDER_H
