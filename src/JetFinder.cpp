#include "JetFinder.h"
#include "JetDefinition.h"
#include "Debug.h"
#include <algorithm>

namespace SlowJet
{

JetList JetFinder::jets()
{
    DEBUG_MSG("Total number of particles: " << m_particles.size());
    m_cones = JetDefinition::instance()->generateCones(m_particles);
    DEBUG_MSG("Total number of cones: " << m_cones.size());
    DEBUG_MSG("====================Start Clustering...====================");
    JetList results{};
    while (not m_particles.empty()) {
        Jet j = findOneJet();
        results.push_back(j);
        DEBUG_MSG("Jet:" << j.jet().prettyOutput() << ' ' << j.content().size());
        for (auto pt : j.content()) {
            m_particles.erase(std::remove(m_particles.begin(), m_particles.end(), pt), m_particles.end());
            DEBUG_MSG("Number of particles left: " << m_particles.size());
            if (not m_cones.empty()) {
                m_cones.erase(std::remove_if(m_cones.begin(), m_cones.end(), [&](const JetCone & c) {
                    return std::find(c.boundary().begin(), c.boundary().end(), pt) != c.boundary().end();
                }), m_cones.end());
                DEBUG_MSG("Number of cones left: " << m_cones.size());
            }
        }
    }
    return results;
}

Jet JetFinder::findOneJet()
{
    Jet jet;
    if (m_cones.empty()) {
        VectorList pt{m_particles.back()};
        m_particles.pop_back();
        Jet single(pt);
        return single;
    }

    double maxJetFunction = -100000;

    for (auto cone : m_cones) {
        Vector c = cone.center();
        double r = cone.radius();
        VectorList b = cone.boundary();
        VectorList inner(m_particles.size());
        auto it = std::copy_if(m_particles.begin(), m_particles.end(), inner.begin(), [&](const Vector & pt) {
            return JetDefinition::instance()->zt(c, pt) > r and std::find(b.begin(),b.end(), pt) == b.end();
        });
        inner.resize(std::distance(inner.begin(), it));
        DEBUG_MSG(inner.size() << " particles inside this cone");

        VectorList subset{};
        subset.reserve(inner.size() + b.size());
        subset.insert(subset.end(),inner.begin(),inner.end());
        subset.insert(subset.end(),b.begin(),b.end());
        Jet protoJet(subset);
        if (protoJet.jetFunction() > maxJetFunction) {
            jet = protoJet;
        }

        if (b.size() == 3) {
            for (unsigned int i =0; i < 2; i++) {
                for (unsigned int j = i+1; j < 3; j++) {
                    subset.clear();
                    subset.insert(subset.end(),inner.begin(),inner.end());
                    subset.push_back(b[i]);
                    subset.push_back(b[j]);
                    Jet protoJet(subset);
                    if (protoJet.jetFunction() > maxJetFunction) {
                        jet = protoJet;
                    }
                }
            }
        } else if (b.size() == 2) {
            for (unsigned int i = 0; i < 2; i++) {
                subset.clear();
                subset.insert(subset.end(),inner.begin(),inner.end());
                subset.push_back(b[i]);
                Jet protoJet(subset);
                if (protoJet.jetFunction() > maxJetFunction) {
                    jet = protoJet;
                }
            }
        } else {
            DEBUG_MSG("The boundary contains: " << b.size() << " particles. Something is wrong!");
        }
    }
    return jet;
}

}
