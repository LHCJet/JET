#ifndef JETCONE_H
#define JETCONE_H

#include "Vector.h"
#include "Debug.h"

namespace SlowJet {

class JetCone {
public:
    JetCone() : m_boundary{}, m_center{}, m_radius(0), m_indices{}, m_jetFunction(-100000), m_keep(-1), m_discarded(false) {}
    ~JetCone() {}
    void setBoundary(const std::vector<unsigned int >& b) { m_boundary = b; }
    void setCenter(const Vector & c) { m_center = c; }
    void setRadius(double r) {m_radius = r; }
    void setJetFunction(double j) {m_jetFunction = j; m_dirty = false;}
    const std::vector<unsigned int > & boundary() const { return m_boundary; }
    const Vector & center() const { return m_center; }
    double radius() const {return m_radius; }
    double jetFunction() const {return m_jetFunction;}
    void addIndex(unsigned int i) { m_indices.push_back(i); }
    void setKeep(unsigned int k) {
        m_keep = k;
        DEBUG_MSG("set exclude: " << m_keep);
    }
    const std::vector<unsigned int > & indices() const { return m_indices; }
    std::vector<unsigned int > content() const {
        std::vector<unsigned int> res(m_indices);
        DEBUG_MSG("exclude: " << m_keep);
        for (unsigned int i = 0; i < m_boundary.size(); i++) {
            DEBUG_MSG(m_boundary[i]);
            if (m_boundary[i] != m_keep) {
                res.push_back(m_boundary[i]);
            }
        }
        return res;
    }
    bool operator == (const JetCone & d) const {
        return m_boundary == d.boundary();
    }
    bool discarded() const {return m_discarded;}
    void discard() {
        DEBUG_MSG("DISCARD CONE");
        m_discarded = true; }
    bool dirty() const { return m_dirty; }
    void markDirty() { m_dirty = true; }
private:
    std::vector<unsigned int > m_boundary;
    std::vector<unsigned int > m_indices;
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
