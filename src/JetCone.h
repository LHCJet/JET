#ifndef JETCONE_H
#define JETCONE_H

#include "Vector.h"

namespace SlowJet {

class JetCone {
public:
    JetCone() : m_boundary{}, m_center{}, m_radius(0) {}
    ~JetCone() {}
    void setBoundary(const VectorList & b) { m_boundary = b; }
    void setCenter(const Vector & c) { m_center = c; }
    void setRadius(double r) {m_radius = r; }
    const VectorList & boundary() const { return m_boundary; }
    const Vector & center() const { return m_center; }
    double radius() const {return m_radius; }
    bool operator == (const JetCone & d) const {
        return m_boundary == d.boundary();
    }
private:
    VectorList m_boundary;
    Vector m_center;
    double m_radius;
};

typedef std::vector<JetCone> JetConeList;

}

#endif // JETCONE_H
