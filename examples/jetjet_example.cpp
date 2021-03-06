/* Any copyright is dedicated to the Public Domain.
 *  * http://creativecommons.org/publicdomain/zero/1.0/ */

#include "Vector.h"
#include "JetDefinition.h"
#include "JetFinder.h"
#include "Debug.h"
#include <iostream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <cmath>

int main()
{
    JETJet::VectorList input_particles;
    std::string line;
    while (std::getline(std::cin, line)) {
        std::istringstream linestream(line);
        // take substrings to avoid problems when there are extra "pollution"
        // characters (e.g. line-feed).
        if (line.substr(0,1) == "#") {continue;}
        double px,py,pz,E;
        linestream >> px >> py >> pz >> E;
        E = sqrt(px*px+py*py+pz*pz);
        input_particles.push_back(JETJet::Vector(px,py,pz,E));
    }
    DEBUG_MSG("================================================================================");
    JETJet::JetDefinition * jetDefinition = new JETJet::EtConeDefinition(6.0);
    JETJet::JetFinder jf(input_particles, jetDefinition);
    JETJet::JetList jets = jf.jets();
    std::sort(jets.begin(),jets.end(),[](const JETJet::Jet & a, const JETJet::Jet & b){
            return a.jet().pt() > b.jet().pt(); // Descendent sorting
    });
    for (unsigned i = 0; i < jets.size(); i++) {
        if (jets[i].jet().pt() < 5) {
            continue;
        }
        JETJet::PArray p = jets[i].jet().fourVector();
        std::cout << std::setprecision(12) << i  << " "
             << p[0] << " "
             << p[1] << " "
             << p[2] << " "
             << p[3] << std::endl;
        JETJet::IndexList cst = jets[i].content();
        for (unsigned j = 0; j < cst.size() ; j++) {
            std::cout << " " << j << " "
               << input_particles[cst[j]].y() << " "
               << input_particles[cst[j]].phi() << " "
               << input_particles[cst[j]].pt() << std::endl;
        }
        std::cout << "#END" << std::endl;
    }
    //JETJet::JetConeList cones = JETJet::JetDefinition::instance()->generateCones(input_particles);
    //for (auto cone : cones) {
    //    JETJet::VectorList b = cone.boundary();
    //    for (auto v : b) {
    //        DEBUG_MSG(v.prettyOutput());
    //    }
    //    DEBUG_MSG("================================================================================");
    //    DEBUG_MSG(cone.center().prettyOutput());
    //    DEBUG_MSG(cone.radius());
    //    DEBUG_MSG("================================================================================");
    //}
    return 0;
}

