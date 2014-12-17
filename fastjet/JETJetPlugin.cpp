#include "JETJetPlugin.h"
#include "Vector.h"
#include "JetDefinition.h"
#include "JetFinder.h"

#include <fastjet/ClusterSequence.hh>
#include <cmath>

FASTJET_BEGIN_NAMESPACE

std::string JETJetPlugin::description() const
{
    return "Test JETJet";
JETJetPlugin::JETJetPlugin(Algorithm alg, double beta, double alpha)
    :JetDefinition::Plugin(), m_beta(beta), m_alpha(alpha)
{
    m_definition = NULL;
    switch(alg) {
        case ECone:
        {
            m_definition = new JETJet::EConeDefinition(m_beta);
            break;
        }
        case EtCone:
        {
            m_definition = new JETJet::EtConeDefinition(m_beta);
            break;
        }
        case EtAlphaCone:
        {
            m_definition = new JETJet::EtAlphaConeDefinition(m_alpha, m_beta);
            break;
        }
        default:
            break;
    }
}

void JETJetPlugin::run_clustering(ClusterSequence & clust_seq) const
{
    const std::vector<PseudoJet> & initial_particles = clust_seq.jets();
    JETJet::VectorList input_particles;
    int jet_i, jet_j, jet_k;
    double dij = 0;
    double px,py,pz,E;
    for (const auto & p : initial_particles) {
        px = p.px();
        py = p.py();
        pz = p.pz();
        E = sqrt(px*px+py*py+pz*pz);
        input_particles.push_back(JETJet::Vector(px,py,pz,E));
    }
    JETJet::JetFinder jf(input_particles, m_definition);
    JETJet::JetList jets = jf.jets();
    for (const auto & jet : jets) {
        JETJet::IndexList cst = jet.content();
        jet_k = cst[0];
        if (cst.size() > 1) {
            for (unsigned int j = 1; j < cst.size(); j++) {
                jet_i = jet_k;
                jet_j = cst[j];
                clust_seq.plugin_record_ij_recombination(jet_i, jet_j, dij, jet_k);
            }
        }
        // we have merged all the jet's particles into a single object, so now
        // "declare" it to be a beam (inclusive) jet.
        // [NB: put a sensible looking d_iB just to be nice...]
        double d_iB = clust_seq.jets()[jet_k].perp2();
        clust_seq.plugin_record_iB_recombination(jet_k, d_iB);
    }
}

FASTJET_END_NAMESPACE
