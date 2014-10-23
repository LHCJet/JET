#include "JetDefinition.h"
#include "Debug.h"
#include <numeric>
#include <algorithm>
#include <functional>
#include <cmath>
#include <array>

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

PArray JetDefinition::sumP(const IndexList & indices, const VectorList & particles)
{
    PArray jetP{0,0,0,0};
    for(auto i: indices) {
        PArray p = particles[i].fourVector();
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
    //cone.setBoundary(VectorList{pt1, pt2, pt3});
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
    //cone.setBoundary(VectorList{pt1, pt2});
    cone.setCenter(center);
    cone.setRadius(radius);
    return cone;
}

double ** JetDefinition::generateDistanceTable(const VectorList & particles)
{
    unsigned int n = particles.size();
    double ** distances  = new double *[n];
    double cos2th = 2*m_b*m_b-1;
    for (unsigned int i = 0; i < n; i++) {
        distances[i] = new double[n];
        PArray p = particles[i].normalizedFourVector();
        Vector v(p[0],p[1],p[2]*cos2th, p[3]);
        for (unsigned int j = 0; j < n; j++) {
            if (zt(v, particles[j]) > cos2th/sqrt(1-(1-cos2th*cos2th)*p[2]*p[2])) {
                distances[i][j] = zt(particles[i],particles[j]);
            } else {
                distances[i][j] = -1.0;
            }
            //printf ("%f ", distances[i][j]);
        }
        //printf ("\n");
    }
    return distances;
}

JetConeList JetDefinition::generateCones(VectorList & particles, double ** distances)
{
    JetConeList cones{};
    double cos2th = 2*m_b*m_b-1;
    generateDistanceTable(particles);
    //std::cout << cos2th << std::endl;
    unsigned int cone_index = 0;
    for (unsigned int i = 0; i < particles.size() - 2; i++) {
        unsigned int count_n =0;
        unsigned int count_n2 =0;
        PArray pp = particles[i].normalizedFourVector();
        Vector xp(pp[0],pp[1],cos2th*pp[2],pp[3]);
        for (unsigned int j = i+1; j < particles.size() - 1; j++) {
            if (distances[i][j] < 0) {
                continue;
            }
            count_n++;
            for (unsigned int k = j+1; k < particles.size(); k++) {
                if (distances[i][k] < 0 or distances[j][k] < 0) {
                    continue;
                }
                JetCone cone = findCone(particles[i], particles[j], particles[k]);
                Vector c = cone.center();
                PArray p = c.normalizedFourVector();
                if (cone.radius() > sqrt(1+(1/m_b/m_b-1)*p[2]*p[2])*m_b) {
                    for (unsigned int l = 0; l < particles.size(); l++) {
                        if (l == i or l == j or l == k) {
                            particles[l].addAssociatedBoundaries(cone_index);
                            continue;
                        }
                        if (distances[i][l] < 0) {
                            continue;
                        }
                        if (zt(c, particles[l]) >= cone.radius()) {
                            particles[l].addAssociatedInners(cone_index);
                            cone.addIndex(l);
                        }
                    }
                    cone.setBoundary (std::vector<unsigned int > {i,j,k});
                    cones.push_back(cone);
                    cone_index++;
                    DEBUG_MSG("Cone contains: " << cone.indices().size() << " particles");
                    //count_n2++;
                }
            }
        }
//        std::cout << count_n << " particles in the fiducial region." << std::endl;
    }
    for (unsigned int i = 0; i < particles.size() - 1; i++) {
        for (unsigned int j = i+1; j < particles.size(); j++) {
            if (distances[i][j] < 0) {
                continue;
            }
            JetCone cone = findCone(particles[i], particles[j]);
            Vector c = cone.center();
            PArray p = c.normalizedFourVector();
            if (cone.radius() > sqrt(1+(1/m_b/m_b-1)*p[2]*p[2])*m_b) {
                for (unsigned int l = 0; l < particles.size(); l++) {
                    if (l == i or l == j) {
                        particles[l].addAssociatedBoundaries(cone_index);
                        continue;
                    }
                    if (distances[i][l] < 0) {
                        continue;
                    }
                    if (zt(c, particles[l]) >= cone.radius()) {
                        particles[l].addAssociatedInners(cone_index);
                        cone.addIndex(l);
                    }
                }
                cone.setBoundary (std::vector<unsigned int > {i,j});
                cones.push_back(cone);
                cone_index++;
                DEBUG_MSG("Cone contains: " << cone.indices().size() << " particles");
            }
        }
    }
    //DEBUG_MSG(cones.size() << " cones generated!");
    std::cout << cones.size() << " cones generated!" << std::endl;
    return cones;
}

}
