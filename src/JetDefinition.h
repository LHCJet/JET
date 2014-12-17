/* Any copyright is dedicated to the Public Domain.
 *  * http://creativecommons.org/publicdomain/zero/1.0/ */

#ifndef JETDEFINITION_H
#define JETDEFINITION_H

#include "Vector.h"
#include "JetCone.h"

namespace JETJet {

class JetDefinition
{
public:
    JetDefinition(double beta);
    ~JetDefinition() {}
    double zt(const Vector &, const Vector &);
    double beta() const {return m_beta;}
    PArray sumP(const IndexList &, const VectorList &);
    JetConeList generateCones(VectorList &);
    virtual double jetFunction(const PArray &) const = 0;
    virtual std::string description() const = 0;

protected:
    virtual double fiducialBoundary(const Vector & pt) const = 0;
    virtual Vector fiducialCenter(const Vector & pt) const = 0;
    virtual double coneBoundary(const PArray & center) const = 0;
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
    virtual std::string description() const;

protected:
    virtual double fiducialBoundary(const Vector & pt) const;
    virtual Vector fiducialCenter(const Vector & pt) const;
    virtual double coneBoundary(const PArray & center) const;

private:
    double m_b;
    double m_cos2th;
};

class EtAlphaConeDefinition : public JetDefinition
{
public:
    EtAlphaConeDefinition(double alpha, double beta);
    virtual double jetFunction(const PArray &) const;
    virtual std::string description() const;

protected:
    virtual double fiducialBoundary(const Vector & pt) const;
    virtual Vector fiducialCenter(const Vector & pt) const;
    virtual double coneBoundary(const PArray & center) const;

private:
    double m_b;
    double m_alpha;
    double m_cos2th;
};

class EConeDefinition : public JetDefinition
{
public:
    EConeDefinition(double beta);
    virtual double jetFunction(const PArray &) const;
    virtual std::string description() const;

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
