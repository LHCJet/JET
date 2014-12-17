/* Any copyright is dedicated to the Public Domain.
 *  * http://creativecommons.org/publicdomain/zero/1.0/ */

#include "JETJetPlugin.h"
#include <fastjet/JetDefinition.hh>
#include <fastjet/PseudoJet.hh>
#include <fastjet/ClusterSequence.hh>
#include <iostream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <stdio.h>

int main()
{
    std::vector<fastjet::PseudoJet> input_particles;
    std::string line;
    while (std::getline(std::cin, line)) {
        std::istringstream linestream(line);
        // take substrings to avoid problems when there are extra "pollution"
        // characters (e.g. line-feed).
        if (line.substr(0,1) == "#") {continue;}
        double px,py,pz,E;
        linestream >> px >> py >> pz >> E;
        input_particles.push_back(fastjet::PseudoJet(px,py,pz,E));
    }
    fastjet::JETJetPlugin * plugin = new fastjet::JETJetPlugin(fastjet::JETJetPlugin::EtCone, 6.0);
    fastjet::JetDefinition jet_def(plugin);

    // run the jet clustering with the above jet definition
    //----------------------------------------------------------
    fastjet::ClusterSequence clust_seq(input_particles, jet_def);


    // get the resulting jets ordered in pt
    //----------------------------------------------------------
    double ptmin = 5.0;
    std::vector<fastjet::PseudoJet> inclusive_jets = fastjet::sorted_by_pt(clust_seq.inclusive_jets(ptmin));


    // tell the user what was done
    //  - the description of the algorithm used
    //  - extract the inclusive jets with pt > 5 GeV
    //    show the output as•
    //      {index, rap, phi, pt}
    //----------------------------------------------------------
    std::cout << "Ran " << jet_def.description() << std::endl;

    // label the columns
    printf("%5s %15s %15s %15s\n","jet #", "rapidity", "phi", "pt");
    // print out the details for each jet
    for (unsigned int i = 0; i < inclusive_jets.size(); i++) {
        printf("%5u %15.8f %15.8f %15.8f\n",
            i, inclusive_jets[i].rap(), inclusive_jets[i].phi(),
            inclusive_jets[i].perp());
    }
    return 0;
}

