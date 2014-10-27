#include "Vector.h"
#include "JetDefinition.h"
#include <cmath>
#include <sstream>

namespace SlowJet {

Vector::Vector(const double px, const double py, const double pz, const double E)
    :m_p{px, py, pz, E}, m_p_normalized{0, 0, 0, 0}, m_associatedInners{}, m_associatedBoundaries{}, m_discarded(false), m_jetFunction(0)
{
    normalizeVector();
}

void Vector::normalizeVector()
{
    double p = sqrt(m_p[0]*m_p[0]+m_p[1]*m_p[1]+m_p[2]*m_p[2]);
    for(int i = 0; i < 4; i++)
        m_p_normalized[i] = m_p[i]/p;
}

std::string Vector::prettyOutput() const
{
    std::ostringstream stringStream;
    stringStream.precision(10);
    stringStream  << eta() << ' ' << phi() << ' ' << pt();// << ' ' << m_p[3];
    return stringStream.str();
}

}
