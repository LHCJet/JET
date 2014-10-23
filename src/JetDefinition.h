#ifndef JETDEFINITION_H
#define JETDEFINITION_H

#include "Vector.h"
#include "JetCone.h"

namespace SlowJet {

class JetDefinition
{
public:
    static JetDefinition * instance();
    double jetFunction(const PArray &);
    double zt(const Vector &, const Vector &);
    void setBeta(double);
    double beta() const {return m_beta;}
    PArray sumP(const VectorList &);
    JetConeList generateCones(const VectorList &, double **);
    double ** generateDistanceTable(const VectorList &);

private:
    JetDefinition() : m_betaMin(2.0), m_beta(2.0), m_b(0.5) {}
    ~JetDefinition() {}
    JetCone findCone(const Vector &, const Vector &, const Vector &);
    JetCone findCone(const Vector &, const Vector &);
    static JetDefinition * m_instance;
    double m_beta;
    double m_betaMin;
    double m_b;
};

}

#endif // JETDEFINITION_H
