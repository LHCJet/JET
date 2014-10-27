#include "JetFinder.h"
#include "JetDefinition.h"
#include "Debug.h"
#include <algorithm>

namespace SlowJet
{
JetFinder::JetFinder(VectorList & particles, JetDefinition * jetDefinition)
    : m_particles{particles}, m_cones{}
{
    m_jetDefinition = jetDefinition;
    for (auto & p : m_particles) {
        p.setJetFunction(m_jetDefinition->jetFunction(p.fourVector()));
    }
}

JetList JetFinder::jets()
{
    DEBUG_MSG("Total number of particles: " << m_particles.size());
    m_cones = m_jetDefinition->generateCones(m_particles);
    DEBUG_MSG("Total number of cones: " << m_cones.size());
    DEBUG_MSG("====================Start Clustering...====================");
    JetList results{};
    int n = 0;
    do {
        Jet j = findOneJet();
        results.push_back(j);
        DEBUG_MSG("Jet:" << j.jet().prettyOutput() << ' ' << j.content().size());
        n = 0;
        DEBUG_MSG("===================Remaining Cones==========================");
        for (unsigned int i = 0; i < m_cones.size(); i++) {
            if (not m_cones[i].discarded()) {
                n += 1;
            }
        }
        DEBUG_MSG(n << " cones left");
    } while (n > 0);

    for (unsigned int i = 0; i < m_particles.size(); i++) {
        if (not m_particles[i].discarded()) {
            IndexList ids{i};
            Jet single(ids, m_particles);
            single.setVector(m_particles[i]);
            results.push_back(single);
        }
    }

    return results;
}

Jet JetFinder::findOneJet()
{

    double maxJetFunction = -100000;

    JetCone jetCone;
    int keep = -1;
    //for (unsigned int i_cone = 0; i_cone < m_cones.size(); i_cone++) {
    unsigned int n = 0;
    for (auto & cone : m_cones) {
        keep = -1;
        if (cone.discarded()) {
            continue;
        }
        if (cone.dirty()) {
            Vector c = cone.center();
            double r = cone.radius();
            double jf = 0;
            IndexList b = cone.boundary();
            //VectorList inner(cone.indices().size());
            IndexList inner{};
            for (auto i : cone.indices()) {
                if (not m_particles[i].discarded()) {
                    inner.push_back(i);
                }
            }

            DEBUG_MSG(inner.size() << " particles inside this cone");

            IndexList subset{};
            subset.reserve(inner.size() + b.size());
            subset.insert(subset.end(),inner.begin(),inner.end());
            subset.insert(subset.end(),b.begin(),b.end());

            PArray pj = m_jetDefinition->sumP(subset, m_particles);
            jf = m_jetDefinition->jetFunction(pj);

            cone.setJetFunction(jf);
            cone.setKeep(-1);

            if (jf > maxJetFunction) {
                DEBUG_MSG("Jet Function Value:" << jf << " " << subset.size());
                jetCone = cone;
                maxJetFunction = jf;
            }

            if (b.size() == 3) {
                for (unsigned int i =0; i < 2; i++) {
                    for (unsigned int j = i+1; j < 3; j++) {
                        subset.clear();
                        subset.insert(subset.end(),inner.begin(),inner.end());
                        subset.push_back(b[i]);
                        subset.push_back(b[j]);
                        pj = m_jetDefinition->sumP(subset, m_particles);
                        jf = m_jetDefinition->jetFunction(pj);
                        keep = b[3 - i - j]; // the vector not included;
                        if (jf > cone.jetFunction()) {
                            cone.setJetFunction(jf);
                            cone.setKeep(keep);
                        }
                        if (jf > maxJetFunction) {
                            DEBUG_MSG("Jet Function Value:" << jf << " " << subset.size());
                            jetCone = cone;
                            maxJetFunction = jf;
                        }
                    }
                }
            }
        } else {
            DEBUG_MSG("Cached results" << cone.jetFunction());
            if (maxJetFunction < cone.jetFunction()) {
                jetCone = cone;
                maxJetFunction = cone.jetFunction();
            }
        }
    }

    int single = -1;

    for(unsigned int i = 0; i < m_particles.size(); i++) {
        if (not m_particles[i].discarded() and m_particles[i].jetFunction() > maxJetFunction) {
            maxJetFunction = m_particles[i].jetFunction();
            single = i;
        }
    }

    IndexList content{};

    Vector vj{};

    if (single >= 0) {
        content.push_back((unsigned int)single);
        vj = m_particles[single];
    } else {
        for (auto i : jetCone.content()) {
            if (not m_particles[i].discarded()) {
                content.push_back(i);
            }
        }
        PArray pj = m_jetDefinition->sumP(content,m_particles);
        vj = Vector(pj[0],pj[1],pj[2],pj[3]);
        vj.setJetFunction(maxJetFunction);
    }


    Jet jet = Jet(content, m_particles);
    jet.setVector(vj);
    DEBUG_MSG("Find cone " << jet.content().size() << " " << maxJetFunction);
    //std::vector<unsigned int > content(jetCone.indices());
    for (auto i: content) {
        DEBUG_MSG("particle: " << i);
        m_particles[i].discard();
        DEBUG_MSG("Discard " << m_particles[i].associatedBoundaries().size() << " cones!");
        DEBUG_MSG(m_particles[i].associatedInners().size() << " dirty cones!");
        for (auto j: m_particles[i].associatedBoundaries()) {
            m_cones[j].discard();
        }
        for (auto j: m_particles[i].associatedInners()) {
            m_cones[j].markDirty();
        }
    }
    //std::cout << "continue" << std::endl;
    return jet;
}

}
