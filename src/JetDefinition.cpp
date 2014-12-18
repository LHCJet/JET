/* Any copyright is dedicated to the Public Domain.
 *  * http://creativecommons.org/publicdomain/zero/1.0/ */

#include "JetDefinition.h"
#include "Debug.h"
#include <numeric>
#include <algorithm>
#include <functional>
#include <cmath>
#include <array>
#include <sstream>

namespace JETJet {

JetDefinition::JetDefinition(double beta)
    : m_beta(beta)
{
    if (beta < 1.0){
        m_beta = 1.0;
        DEBUG_MSG("beta is too small:" << beta << ", set it to 1.");
    }
    std::cout << "================================================================================" << std::endl;
    std::cout << "==   JETJet: Implementation of JET algorithms described in arXiv: 1411.3705   ==" << std::endl;
    std::cout << "================================================================================" << std::endl;
    std::cout << std::endl;
}

double JetDefinition::zt(const Vector & pt1, const Vector & pt2)
{
    PArray A = pt1.normalizedFourVector();
    PArray B = pt2.normalizedFourVector();
    return A[0]*B[0]+A[1]*B[1]+A[2]*B[2];
}

PArray JetDefinition::sumP(const IndexList & indices, const VectorList & particles)
{
    PArray jetP{{0,0,0,0}};
    for(auto i: indices) {
        PArray p = particles[i].fourVector();
        jetP[0] += p[0];
        jetP[1] += p[1];
        jetP[2] += p[2];
        jetP[3] += p[3];
    }
    return jetP;
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
    for (unsigned int i = 0; i < n; i++) {
        distances[i] = new double[n];
        Vector v = fiducialCenter(particles[i]);
        double boundary = fiducialBoundary(particles[i]);
        for (unsigned int j = 0; j < n; j++) {
            if (zt(v, particles[j]) > boundary) {
                distances[i][j] = zt(particles[i],particles[j]);
            } else {
                distances[i][j] = -1.1;
            }
        }
    }
    return distances;
}

JetConeList JetDefinition::generateCones(VectorList & particles)
{
    JetConeList cones{};
    double ** distances = generateDistanceTable(particles);
    unsigned int cone_index = 0;
    for (unsigned int i = 0; i < particles.size() - 2; i++) {
        for (unsigned int j = i+1; j < particles.size() - 1; j++) {
            if (distances[i][j] < -1.0) {
                continue;
            }
            for (unsigned int k = j+1; k < particles.size(); k++) {
                if (distances[i][k] < -1.0 or distances[j][k] < -1.0) {
                    continue;
                }
                JetCone cone = findCone(particles[i], particles[j], particles[k]);
                Vector c = cone.center();
                PArray p = c.normalizedFourVector();
                if (cone.radius() > coneBoundary(p)) {
                    for (unsigned int l = 0; l < particles.size(); l++) {
                        if (distances[i][l] < -1.0) {
                            continue;
                        }
                        if (l == i or l == j or l == k) {
                            particles[l].addAssociatedBoundaries(cone_index);
                            continue;
                        }
                        if (zt(c, particles[l]) >= cone.radius()) {
                            particles[l].addAssociatedInners(cone_index);
                            cone.addIndex(l);
                        }
                    }
                    cone.setBoundary(IndexList{i,j,k});
                    cones.push_back(cone);
                    cone_index++;
                    DEBUG_MSG("Cone contains: " << cone.indices().size() << " particles");
                }
            }
        }
    }
    for (unsigned int i = 0; i < particles.size() - 1; i++) {
        for (unsigned int j = i+1; j < particles.size(); j++) {
            if (distances[i][j] < -1.0) {
                continue;
            }
            JetCone cone = findCone(particles[i], particles[j]);
            Vector c = cone.center();
            PArray p = c.normalizedFourVector();
            if (cone.radius() > coneBoundary(p)) {
                for (unsigned int l = 0; l < particles.size(); l++) {
                    if (distances[i][l] < -1.0) {
                        continue;
                    }
                    if (l == i or l == j) {
                        particles[l].addAssociatedBoundaries(cone_index);
                        continue;
                    }
                    if (zt(c, particles[l]) >= cone.radius()) {
                        particles[l].addAssociatedInners(cone_index);
                        cone.addIndex(l);
                    }
                }
                cone.setBoundary(IndexList{i,j});
                cones.push_back(cone);
                cone_index++;
                DEBUG_MSG("Cone contains: " << cone.indices().size() << " particles");
            }
        }
    }

    for (unsigned int i = 0; i < particles.size(); i++) {
        delete [] distances[i];
    }
    delete [] distances;

    DEBUG_MSG(cones.size() << " cones generated!");
    //std::cout << cones.size() << " cones generated!" << std::endl;
    return cones;
}

EtConeDefinition::EtConeDefinition(double beta)
    :JetDefinition(beta), m_b(0), m_cos2th(0)
{
    DEBUG_MSG("New beta: " << m_beta);
    m_b = sqrt(1.0 - 1.0/beta);
    m_cos2th = 2*m_b*m_b-1;
}

std::string EtConeDefinition::description() const
{
    std::ostringstream desc;
    desc << "JET algorithm (EtCone definition) with beta = " << m_beta;
    return desc.str();
}

double EtConeDefinition::jetFunction(const PArray & jetP) const
{
    double Et = sqrt(jetP[3]*jetP[3] - jetP[2]*jetP[2]);
    return ((1-m_beta)*Et*Et + m_beta*(jetP[0]*jetP[0]+jetP[1]*jetP[1]))/Et;
}

Vector EtConeDefinition::fiducialCenter(const Vector & pt) const
{
    PArray p = pt.fourVector();
    return Vector(p[0],p[1],p[2]*m_cos2th,p[3]);
}

double EtConeDefinition::fiducialBoundary(const Vector & pt) const
{
    PArray p = pt.normalizedFourVector();
    return m_cos2th/sqrt(1-(1 - m_cos2th*m_cos2th)*p[2]*p[2]);
}

double EtConeDefinition::coneBoundary(const PArray & center) const
{
    return m_b*sqrt(1+(1/m_b/m_b -1)*center[2]*center[2]);
}

EtAlphaConeDefinition::EtAlphaConeDefinition(double alpha, double beta)
    :JetDefinition(beta), m_b(0), m_alpha(alpha), m_cos2th(0)
{
    if (alpha < 0 or alpha > 2) {
        DEBUG_MSG("Wrong alpha, reset it to 1");
        m_alpha = 1;
    }
    DEBUG_MSG("New alpha: " << m_alpha << ", beta: " << m_beta);
    //ToCheck: whether this is still true in the forward region
    if (alpha < 1) {
        m_b = sqrt(1.0 - 1.0 / beta);
    } else {
        double v_maximize = sqrt(alpha/(alpha-2)*(beta -1)/beta);
        if ( v_maximize < 1.0) {
            m_b = sqrt(alpha * (2.0 - alpha) * (1.0 - 1.0 / beta));
        } else {
            m_b = 1.0 - 0.5 * alpha / beta;
        }
    }
    m_cos2th = 2*m_b*m_b-1;
}

std::string EtAlphaConeDefinition::description() const
{
    std::ostringstream desc;
    desc << "JET algorithm (EtAlphaCone definition) with alpha = " << m_alpha << " and beta = " << m_beta;
    return desc.str();
}

double EtAlphaConeDefinition::jetFunction(const PArray & jetP) const
{
    double Et = sqrt(jetP[3]*jetP[3] - jetP[2]*jetP[2]);
    return pow(Et,m_alpha)*(1 - m_beta*(jetP[3]*jetP[3]-jetP[0]*jetP[0]-jetP[1]*jetP[1]-jetP[2]*jetP[2])/(Et*Et));
}

Vector EtAlphaConeDefinition::fiducialCenter(const Vector & pt) const
{
    PArray p = pt.fourVector();
    return Vector(p[0],p[1],p[2]*m_cos2th,p[3]);
}

double EtAlphaConeDefinition::fiducialBoundary(const Vector & pt) const
{
    PArray p = pt.normalizedFourVector();
    return m_cos2th/sqrt(1-(1 - m_cos2th*m_cos2th)*p[2]*p[2]);
}

double EtAlphaConeDefinition::coneBoundary(const PArray & center) const
{
    return m_b*sqrt(1+(1/m_b/m_b -1)*center[2]*center[2]);
}

EConeDefinition::EConeDefinition(double beta)
    :JetDefinition(beta), m_b(0), m_cos2th(0)
{
    DEBUG_MSG("New beta: " << m_beta);
    m_b = sqrt(1.0 - 1.0 / beta);
    m_cos2th = 2*m_b*m_b-1;
}

std::string EConeDefinition::description() const
{
    std::ostringstream desc;
    desc << "JET algorithm (ECone definition) with beta = " << m_beta;
    return desc.str();
}

double EConeDefinition::jetFunction(const PArray & jetP) const
{
    return (1-m_beta)*jetP[3]+ m_beta*(jetP[0]*jetP[0]+jetP[1]*jetP[1]+jetP[2]*jetP[2])/jetP[3];
}

Vector EConeDefinition::fiducialCenter(const Vector & pt) const
{
    return pt;
}

double EConeDefinition::fiducialBoundary(const Vector &) const
{
    return m_cos2th;
}

double EConeDefinition::coneBoundary(const PArray &) const
{
    return m_b;
}

}
