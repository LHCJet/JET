#include "JetDefinition.h"
#include "Debug.h"
#include <numeric>
#include <algorithm>
#include <functional>
#include <cmath>

namespace SlowJet {

JetDefinition * JetDefinition::m_instance = 0;

JetDefinition * JetDefinition::instance()
{
    if (!m_instance) {
        m_instance = new JetDefinition;
    }
    return m_instance;
}

void JetDefinition::setBeta(double beta)
{
    if (beta < m_betaMin){
        DEBUG_MSG("beta is too small:" << beta);
        return;
    }
    DEBUG_MSG("New beta: " << beta);
    m_beta = beta;
    m_b = 1.0 - 1.0 / beta;
}

double JetDefinition::zt(const Vector & pt1, const Vector & pt2)
{
    PArray A = pt1.normalizedFourVector();
    PArray B = pt2.normalizedFourVector();
    return A[0]*B[0]+A[1]*B[1]+A[2]*B[2];
}

PArray JetDefinition::sumP(const VectorList & particles)
{
    PArray jetP{0,0,0,0};
    for(auto pt: particles) {
        PArray p = pt.fourVector();
        std::transform(p.begin(), p.end(), jetP.begin(), jetP.begin(),
                std::plus<double>());
    }
    return jetP;
}

double JetDefinition::jetFunction(const PArray & jetP)
{
    double Et2 = jetP[3]*jetP[3] - jetP[2]*jetP[2];
    return (1-m_beta)*Et2 + m_beta*(jetP[0]*jetP[0]+jetP[1]*jetP[1]);
}

JetCone JetDefinition::findCone(const Vector & pt1, const Vector & pt2, const Vector & pt3)
{
    JetCone cone;
    PArray A = pt1.normalizedFourVector();
    PArray B = pt2.normalizedFourVector();
    PArray C = pt3.normalizedFourVector();
    double px = (A[2]*(B[1]-C[1])+B[2]*C[1]-B[1]*C[2]+A[1]*(C[2]-B[2]));
    double py = (A[2]*(C[0]-B[0])+B[0]*C[2]-B[2]*C[0]+A[0]*(B[2]-C[2]));
    double pz = (A[1]*(B[0]-C[0])+B[1]*C[0]-B[0]*C[1]+A[0]*(C[1]-B[1]));
    Vector center(px, py, pz, 1.0);
    double radius = zt(center,pt1);
    cone.setBoundary(VectorList{pt1, pt2, pt3});
    if (radius > 0) {
        cone.setCenter(center);
        cone.setRadius(radius);
    } else {
        cone.setCenter(Vector(-px, -py, -pz, 1.0));
        cone.setRadius(-radius);
    }
    return cone;
}

JetCone JetDefinition::findCone(const Vector & pt1, const Vector & pt2)
{
    JetCone cone;
    PArray A = pt1.normalizedFourVector();
    PArray B = pt2.normalizedFourVector();
    double px = (A[0]+B[0])/2.0;
    double py = (A[1]+B[1])/2.0;
    double pz = (A[2]+B[2])/2.0;
    Vector center(px, py, pz, 1.0);
    double radius = zt(center,pt1);
    cone.setBoundary(VectorList{pt1, pt2});
    cone.setCenter(center);
    cone.setRadius(radius);
    return cone;
}

JetConeList JetDefinition::generateCones(VectorList & particles)
{
    JetConeList cones{};
    double cos2th = 2*m_b*m_b-1;
    for (unsigned int i = 0; i < particles.size() - 2; i++) {
        PArray pp = particles[i].normalizedFourVector();
        Vector xp(pp[0],pp[1],cos2th*pp[2],pp[3]);
        for (unsigned int j = i+1; j < particles.size() - 1; j++) {
            if (zt(xp, particles[j]) < cos2th/sqrt(1-(1-cos2th*cos2th)*pp[2]*pp[2])) {
                continue;
            }
            for (unsigned int k = j+1; k < particles.size(); k++) {
                if (zt(xp, particles[k]) < cos2th/sqrt(1-(1-cos2th*cos2th)*pp[2]*pp[2])) {
                    continue;
                }
                JetCone cone = findCone(particles[i], particles[j], particles[k]);
                PArray p = cone.center().normalizedFourVector();
                if (cone.radius() > sqrt(1+(1/m_b/m_b-1)*p[2]*p[2])*m_b) {
                    cones.push_back(cone);
                }
            }
        }
    }
    for (unsigned int i = 0; i < particles.size() - 1; i++) {
        for (unsigned int j = i+1; j < particles.size(); j++) {
            JetCone cone = findCone(particles[i], particles[j]);
            PArray p = cone.center().normalizedFourVector();
            if (cone.radius() > sqrt(1+(1/m_b/m_b-1)*p[2]*p[2])*m_b) {
                cones.push_back(cone);
            }
        }
    }
    DEBUG_MSG(cones.size() << " cones generated!");
    return cones;
}

}
