/* Any copyright is dedicated to the Public Domain.
 *  * http://creativecommons.org/publicdomain/zero/1.0/ */

#ifndef JETCONE_H
#define JETCONE_H

#include "Vector.h"
#include "Debug.h"

namespace JETJet {

class JetCone {
public:
    JetCone() : m_boundary{}, m_indices{}, m_center{}, m_radius(0), m_jetFunction(-100000), m_keep(-1), m_discarded(false), m_dirty(true) {}
    ~JetCone() {}
    void setBoundary(const IndexList & b) { m_boundary = b; }
    void setCenter(const Vector & c) { m_center = c; }
    void setRadius(double r) {m_radius = r; }
    void setJetFunction(double j) {m_jetFunction = j; m_dirty = false;}
    const IndexList & boundary() const { return m_boundary; }
    const Vector & center() const { return m_center; }
    double radius() const {return m_radius; }
    double jetFunction() const {return m_jetFunction;}
    void addIndex(unsigned int i) { m_indices.push_back(i); }
    void setKeep(unsigned int k) { m_keep = k; }
    const IndexList & indices() const { return m_indices; }
    IndexList content() const {
        IndexList res(m_indices);
        for (unsigned int i = 0; i < m_boundary.size(); i++) {
            if ((m_keep < 0) or (m_boundary[i] != static_cast<unsigned int>(m_keep))) {
                res.push_back(m_boundary[i]);
            }
        }
        return res;
    }
    bool operator == (const JetCone & d) const {
        return m_boundary == d.boundary();
    }
    bool discarded() const {return m_discarded;}
    void discard() { m_discarded = true; }
    bool dirty() const { return m_dirty; }
    void markDirty() { m_dirty = true; }
private:
    IndexList m_boundary;
    IndexList m_indices;
    Vector m_center;
    double m_radius;
    double m_jetFunction;
    int m_keep;
    bool m_discarded;
    bool m_dirty;
};

typedef std::vector<JetCone > JetConeList;

}

#endif // JETCONE_H
