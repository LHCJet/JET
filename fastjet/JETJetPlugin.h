#ifndef JETJETPLUGIN_H
#define JETJETPLUGIN_H

#include <fastjet/JetDefinition.hh>

FASTJET_BEGIN_NAMESPACE

class JETJetPlugin : public JetDefinition::Plugin {
public:
    JETJetPlugin(double beta) : m_beta(beta) {}
    virtual std::string description () const;
    virtual void run_clustering(ClusterSequence &) const;
    virtual double R() const {return m_beta;}
private:
    double m_beta;
};

FASTJET_END_NAMESPACE

#endif // JETJETPLUGIN_H
