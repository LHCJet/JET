#ifndef JETJETPLUGIN_H
#define JETJETPLUGIN_H

#include <fastjet/JetDefinition.hh>
#include "JetDefinition.h"

FASTJET_BEGIN_NAMESPACE

class JETJetPlugin : public JetDefinition::Plugin {
public:
    enum Algorithm {
        ECone = 0,
        EtCone,
        EtAlphaCone
    };
    JETJetPlugin(Algorithm alg, double beta, double alpha=1.0);
    ~JETJetPlugin();
    virtual std::string description () const;
    virtual void run_clustering(ClusterSequence &) const;
    virtual double R() const {return m_beta;}
private:
    JETJet::JetDefinition * m_definition;
    double m_beta;
    double m_alpha;
};

FASTJET_END_NAMESPACE

#endif // JETJETPLUGIN_H
