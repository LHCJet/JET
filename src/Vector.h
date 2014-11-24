/* Any copyright is dedicated to the Public Domain.
 *  * http://creativecommons.org/publicdomain/zero/1.0/ */

#ifndef VECTOR_H
#define VECTOR_H

#include <array>
#include <vector>
#include <string>
#include <cmath>

namespace JETJet {

typedef std::array<double, 4 > PArray;
typedef std::vector<unsigned int > IndexList;

class Vector {
public:
    Vector() :m_p{{0,0,0,0}}, m_p_normalized{{0,0,0,0}}, m_associatedInners{}, m_associatedBoundaries{}, m_discarded(false), m_jetFunction(0) {}
    Vector(const double px, const double py, const double pz, const double E);
    const PArray & fourVector() const {return m_p;}
    double pt() const { return sqrt(m_p[0]*m_p[0]+m_p[1]*m_p[1]); }
    double eta() const {
        double pAbs = sqrt(m_p[0]*m_p[0]+m_p[1]*m_p[1]+m_p[2]*m_p[2]);
        return 0.5*log((pAbs+m_p[2])/(pAbs-m_p[2]));
    }
    double y() const { return 0.5*log((m_p[3]+m_p[2])/(m_p[3]-m_p[2])); }
    double phi() const {
        double arg = atan2(m_p[1],m_p[0]);
        if (arg < 0) {
            arg += 2*M_PI;
        }
        return arg;
    }
    const PArray & normalizedFourVector() const {return m_p_normalized;}
    void addAssociatedInners(unsigned int c) { m_associatedInners.push_back(c); }
    void addAssociatedBoundaries(unsigned int c) { m_associatedBoundaries.push_back(c); }
    const IndexList & associatedInners() const { return m_associatedInners; }
    const IndexList & associatedBoundaries() const { return m_associatedBoundaries; }
    bool operator == (const Vector & d) const {
        return m_p == d.fourVector();
    }
    std::string prettyOutput() const;

    bool discarded() { return m_discarded; }
    void discard() { m_discarded = true; }

    double jetFunction() const { return m_jetFunction; }
    void setJetFunction(double j) { m_jetFunction = j; }

private:
    void normalizeVector();
    PArray m_p;
    PArray m_p_normalized;
    IndexList m_associatedInners;
    IndexList m_associatedBoundaries;
    bool m_discarded;
    double m_jetFunction;
};

typedef std::vector<Vector> VectorList;
} // namespace JETJet

#endif // VECTOR_H
