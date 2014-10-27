#ifndef JETDEFINITION_H
#define JETDEFINITION_H

#include "Vector.h"
#include "JetCone.h"

namespace SlowJet {

class JetDefinition
{
public:
    JetDefinition(double beta) : m_beta(beta) {}
    ~JetDefinition() {}
    double zt(const Vector &, const Vector &);
    double beta() const {return m_beta;}
    PArray sumP(const IndexList &, const VectorList &);
    JetConeList generateCones(VectorList &);
    virtual double jetFunction(const PArray &) const = 0;

protected:
    virtual double fiducialBoundary(const Vector & pt) const = 0;
    virtual Vector fiducialCenter(const Vector & pt) const = 0;
    virtual double coneBoundary(const PArray & center) const = 0;
    static JetDefinition * m_instance;
    double m_beta;

private:
    JetCone findCone(const Vector &, const Vector &, const Vector &);
    JetCone findCone(const Vector &, const Vector &);
    double ** generateDistanceTable(const VectorList &);
};

class EtConeDefinition : public JetDefinition
{
public:
    EtConeDefinition(double beta);
    virtual double jetFunction(const PArray &) const;

protected:
    virtual double fiducialBoundary(const Vector & pt) const;
    virtual Vector fiducialCenter(const Vector & pt) const;
    virtual double coneBoundary(const PArray & center) const;

private:
    double m_b;
    double m_cos2th;
};

}

#endif // JETDEFINITION_H
