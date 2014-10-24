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
        : m_jet{}, m_content{c} {
        PArray jetP = JetDefinition::instance()->sumP(c, p);
        m_jet = Vector(jetP[0], jetP[1], jetP[2], jetP[3]);
    }
    Jet() : m_jet{}, m_content{} {}
    ~Jet() {}
    const Vector & jet() const {return m_jet;}
    const IndexList & content() const {return m_content;}

private:
    Vector m_jet;
    IndexList m_content;
};

typedef std::vector<Jet> JetList;

class JetFinder
{
public:
    JetFinder(const VectorList & particles) : m_particles{particles}, m_cones{} {}
    ~JetFinder() {
    }
    JetList jets();
private:
    Jet findOneJet();
    Jet protoJet(VectorList subset);
    VectorList m_particles;
    JetConeList m_cones;
};

} // namespace SlowJet;

#endif // JETFINDER_H
