#ifndef VECTOR_H
#define VECTOR_H

#include <array>
#include <vector>
#include <string>
#include <cmath>

namespace SlowJet {

typedef std::array<double, 4> PArray;

class Vector {
public:
    Vector() :m_p{0,0,0,0}, m_p_normalized{0,0,0,0} {}
    Vector(const double px, const double py, const double pz, const double E);
    const PArray & fourVector() const {return m_p;}
    const double pt() const { return sqrt(m_p[0]*m_p[0]+m_p[1]*m_p[1]); }
    const double eta() const {
        double pAbs = sqrt(m_p[0]*m_p[0]+m_p[1]*m_p[1]+m_p[2]*m_p[2]);
        return 0.5*log((pAbs+m_p[2])/(pAbs-m_p[2]));
    }
    const double y() const { return 0.5*log((m_p[3]+m_p[2])/(m_p[3]-m_p[2])); }
    const double phi() const {
        double arg = atan2(m_p[1],m_p[0]);
        if (arg < 0) {
            arg += 2*M_PI;
        }
        return arg;
    }
    const PArray & normalizedFourVector() const {return m_p_normalized;}
    bool operator == (const Vector & d) const {
        return m_p == d.fourVector();
    }
    std::string prettyOutput() const;
private:
    void normalizeVector();
    PArray m_p;
    PArray m_p_normalized;
};

typedef std::vector<Vector> VectorList;
} // namespace SlowJet

#endif // VECTOR_H
