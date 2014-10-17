#include "Vector.h"
#include "JetCone.h"
#include "JetDefinition.h"
#include "JetFinder.h"
#include "Debug.h"
#include <iostream>
#include <vector>

int main()
{
    SlowJet::VectorList input_particles;
    double px, py , pz, E;
    while (std::cin >> px >> py >> pz >> E) {
        input_particles.push_back(SlowJet::Vector(px,py,pz,E));
    }
    for (auto pt : input_particles) {
        DEBUG_MSG("Input particles: " << pt.prettyOutput());
    }
    DEBUG_MSG("================================================================================");
    SlowJet::JetDefinition::instance()->setBeta(15.0);
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

