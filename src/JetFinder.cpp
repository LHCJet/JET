#include "JetFinder.h"
#include "JetDefinition.h"
#include "Debug.h"
#include <algorithm>

namespace SlowJet
{

JetList JetFinder::jets()
{
    DEBUG_MSG("Total number of particles: " << m_particles.size());
    m_distances = JetDefinition::instance()->generateDistanceTable(m_particles);
    m_cones = JetDefinition::instance()->generateCones(m_particles, m_distances);
    DEBUG_MSG("Total number of cones: " << m_cones.size());
    DEBUG_MSG("====================Start Clustering...====================");
    JetList results{};
    int n = 0;
    do {
        Jet j = findOneJet();
        results.push_back(j);
        DEBUG_MSG("Jet:" << j.jet().prettyOutput() << ' ' << j.content().size());
//        for (auto pt : j.content()) {
//            m_particles.erase(std::remove(m_particles.begin(), m_particles.end(), pt), m_particles.end());
//            DEBUG_MSG("Number of particles left: " << m_particles.size());
//            if (not m_cones.empty()) {
//                m_cones.erase(std::remove_if(m_cones.begin(), m_cones.end(), [&](const JetCone & c) {
//                    return std::find(c.boundary().begin(), c.boundary().end(), pt) != c.boundary().end();
//                }), m_cones.end());
//                DEBUG_MSG("Number of cones left: " << m_cones.size());
//            }
//        }
        m_cones.erase(std::remove_if(m_cones.begin(), m_cones.end(), [&](JetCone & c) {
            return c.discarded();
        }), m_cones.end());
        //DEBUG_MSG("Number of cones left: " << m_cones.size());
        std::cout << m_cones.size() << " cones left" << std::endl;
        n = 0;
        for (unsigned int i = 0; i < m_particles.size(); i ++) {
            n += m_discarded[i];
        }
        DEBUG_MSG(n << " particles left");
    } while (n > 0);
    return results;
}

Jet JetFinder::findOneJet()
{
    Jet jet;
    if (m_cones.empty()) {
        //FIXME: no cone
        for (unsigned int i = 0; i < m_particles.size(); i++ ) {
            if (m_discarded[i] > 0) {
                VectorList pt{m_particles[i]};
                m_discarded[i] = 0;
                Jet single(pt);
                return single;
            }
        }
        DEBUG_MSG("Should not reach here");
    }

    double maxJetFunction = -100000;

    JetCone jetCone;
    int keep = -1;
    //for (unsigned int i_cone = 0; i_cone < m_cones.size(); i_cone++) {
    for (auto & cone : m_cones) {
        keep = -1;
        if (cone.dirty()) {
            Vector c = cone.center();
            double r = cone.radius();
            std::vector<unsigned int > b = cone.boundary();
            //VectorList inner(cone.indices().size());
            VectorList inner{};

            for (unsigned int i = 0; i < cone.indices().size(); i++) {
                if (m_discarded[cone.indices()[i]] > 0) {
                    inner.push_back(m_particles[cone.indices()[i]]);
                }
            }

            DEBUG_MSG(inner.size() << " particles inside this cone");

            VectorList subset{};
            subset.reserve(inner.size() + b.size());
            subset.insert(subset.end(),inner.begin(),inner.end());
            for (auto i: b) {
                subset.push_back(m_particles[i]);
            }
            Jet protoJet(subset);

            //if (protoJet.jetFunction() > cone.jetFunction()) {
                //DEBUG_MSG("====Cache result: " << protoJet.jetFunction());
                cone.setJetFunction(protoJet.jetFunction());
            //}
            if (protoJet.jetFunction() > maxJetFunction) {
                DEBUG_MSG("Jet Function Value:" << protoJet.jetFunction() << " " << protoJet.content().size());
                jet = protoJet;
                jetCone = cone;
                maxJetFunction = jet.jetFunction();
            }

            if (b.size() == 3) {
                for (unsigned int i =0; i < 2; i++) {
                    for (unsigned int j = i+1; j < 3; j++) {
                        subset.clear();
                        subset.insert(subset.end(),inner.begin(),inner.end());
                        subset.push_back(m_particles[b[i]]);
                        subset.push_back(m_particles[b[j]]);
                        Jet protoJet(subset);
                        keep = b[3 - i - j]; // the vector not included;
                        if (protoJet.jetFunction() > cone.jetFunction()) {
                            cone.setJetFunction(protoJet.jetFunction());
                            cone.setKeep(keep);
                        }
                        if (protoJet.jetFunction() > maxJetFunction) {
                            DEBUG_MSG("Jet Function Value:" << protoJet.jetFunction() << " " << protoJet.content().size());
                            jetCone = cone;
                            jet = protoJet;
                            maxJetFunction = jet.jetFunction();
                        }
                    }
                }
            } else if (b.size() == 2) {
                for (unsigned int i = 0; i < 2; i++) {
                    subset.clear();
                    subset.insert(subset.end(),inner.begin(),inner.end());
                    subset.push_back(m_particles[b[i]]);
                    Jet protoJet(subset);
                    keep = b[2 - i]; // the vector not included;
                    if (protoJet.jetFunction() > cone.jetFunction()) {
                        cone.setJetFunction(protoJet.jetFunction());
                        cone.setKeep(keep);
                    }
                    if (protoJet.jetFunction() > maxJetFunction) {
                        DEBUG_MSG("Jet Function Value:" << protoJet.jetFunction() << " " << protoJet.content().size());
                        jetCone = cone;
                        jet = protoJet;
                        maxJetFunction = jet.jetFunction();
                    }
                }
            } else {
                DEBUG_MSG("The boundary contains: " << b.size() << " particles. Something is wrong!");
            }
        }
        else {
            DEBUG_MSG("Cached results" << cone.jetFunction());
            if (maxJetFunction < cone.jetFunction()) {
                jetCone = cone;
                maxJetFunction = cone.jetFunction();
            }
        }
    }
    //VectorList inner(cone.indices().size());
    VectorList inner{};
    DEBUG_MSG("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
    DEBUG_MSG(jetCone.content().size());
    DEBUG_MSG("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");

    for (unsigned int i = 0; i < jetCone.content().size(); i++) {
        DEBUG_MSG("!!particle in cone: " << jetCone.content()[i]);
        if (m_discarded[jetCone.content()[i]] > 0) {
            inner.push_back(m_particles[jetCone.content()[i]]);
        }
        DEBUG_MSG("Passed");
    }

    jet = Jet(inner);
    DEBUG_MSG("Find cone " << jet.content().size() << " " << maxJetFunction);
    //std::vector<unsigned int > content(jetCone.indices());
    std::vector<unsigned int > content{};
    for (auto i : jetCone.content()) {
        if (m_discarded[i] > 0) {
            content.push_back(i);
        }
    }
    for (auto i: content) {
        DEBUG_MSG("particle: " << i);
        m_discarded[i] = 0;
        if (not m_cones.empty()) {
            int count = 0;
            for (auto & c: m_cones) {
                if (not c.discarded()) {
                    DEBUG_MSG("size: " << c.indices().size());
                    //for (auto j: c.boundary()) {
                    //    MSG("boundary: " << j);
                    //}
                    //for (auto j: c.indices()) {
                    //    MSG("inner: " << j);
                    //}
                    if (std::find(c.boundary().begin(), c.boundary().end(), i) != c.boundary().end()) {
                        DEBUG_MSG("On the boundary");
                        c.discard();
                        if (c.discarded()) {
                            DEBUG_MSG("REALLY DISCARDED");
                        }
                    }
                    if (std::find(c.indices().begin(), c.indices().end(), i) != c.indices().end()) {
                        DEBUG_MSG("Dirty Cone");
                        //MSG(std::find(c.indices().begin(), c.indices().end(), i))
                        c.markDirty();
                        count += 1;
                    } /*else {
                        MSG("Clean!");
                    }*/
                    DEBUG_MSG("========================================================================");
                } else {
                    DEBUG_MSG("Already discarded!");
                }
            }
            //m_cones.erase(std::remove_if(m_cones.begin(), m_cones.end(), [&](const JetCone & c) {
            //    return std::find(c.boundary().begin(), c.boundary().end(), i) != c.boundary().end();
            //}), m_cones.end());
            std::cout << count << " dirty cones" << std::endl;
        }
    }
    //std::cout << "continue" << std::endl;
    return jet;
}

}
