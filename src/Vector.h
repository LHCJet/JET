#ifndef VECTOR_H
#define VECTOR_H

#include <array>
#include <vector>
#include <string>

namespace SlowJet {

typedef std::array<double, 4> PArray;

class Vector {
public:
    Vector() :m_p{0,0,0,0}, m_p_normalized{0,0,0,0} {}
    Vector(const double px, const double py, const double pz, const double E);
    const PArray & fourVector() const {return m_p;}
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
