#include "Vector.h"
#include "JetCone.h"
#include "JetDefinition.h"
#include "JetFinder.h"
#include "Debug.h"
#include <iostream>
#include <sstream>
#include <vector>

int main()
{
    SlowJet::VectorList input_particles;
    double px, py , pz, E;
    std::string line;
    while (std::getline(std::cin, line)) {
        std::istringstream linestream(line);
        // take substrings to avoid problems when there are extra "pollution"
        // characters (e.g. line-feed).
        if (line.substr(0,1) == "#") {continue;}
        double px,py,pz,E;
        linestream >> px >> py >> pz >> E;
        input_particles.push_back(SlowJet::Vector(px,py,pz,E));
    }
#ifdef DEBUG
    for (auto pt : input_particles) {
        DEBUG_MSG("Input particles: " << pt.prettyOutput());
    }
#endif
    DEBUG_MSG("================================================================================");
    SlowJet::JetDefinition::instance()->setBeta(12.0);
    SlowJet::JetFinder jf(input_particles);
    SlowJet::JetList jets = jf.jets();
    int count = 0;
    for (auto jet : jets) {
        std::cout << count << ' ' <<  jet.jet().prettyOutput() << ' ' << jet.content().size() << std::endl;
        count++;
    }
    //SlowJet::JetConeList cones = SlowJet::JetDefinition::instance()->generateCones(input_particles);
    //for (auto cone : cones) {
    //    SlowJet::VectorList b = cone.boundary();
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

