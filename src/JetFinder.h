/* Any copyright is dedicated to the Public Domain.
 *  * http://creativecommons.org/publicdomain/zero/1.0/ */

#ifndef JETFINDER_H
#define JETFINDER_H

#include "Vector.h"
#include "JetCone.h"
#include "JetDefinition.h"

#include <vector>

namespace JETJet {
class Jet
{
public:
    Jet(IndexList & c)
        : m_jet{}, m_content{c} {}
    Jet() : m_jet{}, m_content{} {}
    ~Jet() {}
    void setVector (const Vector & v) { m_jet = v; }
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
    JetFinder(VectorList & particles, JetDefinition * jetDefinition);
    ~JetFinder() {}
    JetList jets();
private:
    Jet findOneJet();
    Jet protoJet(VectorList subset);
    VectorList m_particles;
    JetConeList m_cones;
    JetDefinition * m_jetDefinition;
};

} // namespace JETJet;

#endif // JETFINDER_H
