/*
  Copyright (c) 2015-2016 Lester Hedges <lester.hedges+vmmc@gmail.com>

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#ifdef ISOTROPIC
#error patchy_disc.cpp cannot be linked to isotropic VMMC library!
#endif

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <csignal>


#include "System.h"
#include "VMMC.h"
//include "parse_input.h"


#ifndef M_PI
    #define M_PI 3.1415926535897932384626433832795
#endif

#include <nlohmann/json.hpp>
#include <fstream>
using json = nlohmann::json;
json ReadJsonFromFile(std::string file_name) {
    try {
        return json::parse(std::ifstream{file_name, std::ios::in});
    } catch (json::parse_error& e) {
        std::cerr << "JSON parse exception : " << e.what() << std::endl;
    } catch (std::ifstream::failure& e) {
        std::cerr << "Stream exception : " << e.what() << std::endl;
    } catch (std::exception& e) {
        std::cerr << "Exception : " << e.what() << std::endl;
    } catch (...) {
        std::cerr << "Unknown error" << std::endl;
    }
    return {};
}


int stop = 0;
void gbl_terminate(int arg) {
	/**
	 * The next time the signal is intercepted, the default handler will be invoked
	 * in order to avoid making the program unkillable.
	 */
	signal(arg, SIG_DFL);
	fprintf(stderr, "# Caught SIGNAL %d; setting stop = 1\n", arg);
	stop = 1;
}

int main(int argc, char** argv){
	if(argc == 1){
        fprintf(stderr, "Usage: %s input.json\n", argv[0]);
        exit(EXIT_FAILURE);
    }

	/**
	 * Here we handle a few SIG* signals: whenever we intercept one of these signals
	 * the program will quit the main loop and die as gracefully as possible.
	 */
	signal(SIGTERM, gbl_terminate);
	signal(SIGABRT, gbl_terminate);
	signal(SIGINT, gbl_terminate);
	signal(SIGUSR2, gbl_terminate);


	//input_file input;
	//loadInputFile(&input, argv[1]);
	//if(input.state == ERROR) exit(1);
    //getInputDouble(&input, "density", &density, 1);


    // Simulation parameters.
    unsigned int dimension = 2;                     // dimension of simulation box
    double interactionRange = 1.0;                  // interaction range. #diameter of patch (in units of particle diameter)
    double interactionEnergy = 1.0;                 // pair interaction energy scale (in units of kBT)
    unsigned int maxInteractions = 36;              // maximum number of interactions per particle

	/**
	 * Load the input file and Initialise
	 */
    json inp;
    inp = ReadJsonFromFile(argv[1]);

    //std::cout << std::setw(4) << inp << std::endl;

    unsigned int nTypes=inp["types"];
    //init topology struct
    Top top(nTypes);

    unsigned int nParticles = 0;
    unsigned int t1, t2, counter;

    //init Ni
    for (unsigned int i=0;i<nTypes;i++){
        t1 = inp["particle_numbers"][i]["type"];
        top.Ni[t1] = inp["particle_numbers"][i]["N"];
        nParticles += top.Ni[inp["particle_numbers"][i]["type"]];
    }

    //init patches
    double dummy;
    for (unsigned int i=0;i<nTypes;i++){
        t1 = inp["patches"][i]["type"];
        top.nPatches[t1] = inp["patches"][i]["nPatches"];
        //std::cout << inp["patches"][i]["angles"] << std::endl;
        for (unsigned int j=0;j<top.nPatches[t1];j++){
            dummy = inp["patches"][i]["angles"][j];
            top.patchAngles[t1].push_back(M_PI * (dummy / 180.0));
        }
    }

    //init pair_coeff
    counter = 0;
    for (unsigned int i=0;i<nTypes;i++){
        for (unsigned int j=0;j<nTypes;j++){
            t1 = inp["pair_coeff"][counter]["type1"];
            t2 = inp["pair_coeff"][counter]["type2"];
            //std::cout << t1 << t2 << std::endl;
            top.epsilon[t1][t2] = inp["pair_coeff"][counter]["epsilon"];
            top.delta[t1][t2] = inp["pair_coeff"][counter]["delta"];
            top.sigma[t1][t2] = inp["pair_coeff"][counter]["sigma"];
            top.sigma_p[t1][t2] = inp["pair_coeff"][counter]["sigma_p"];
            top.rcut[t1][t2] = inp["pair_coeff"][counter]["rcut"];
            counter += 1;

            if (interactionEnergy<top.epsilon[t1][t2]) interactionEnergy=top.epsilon[t1][t2];
            if (interactionRange<top.rcut[t1][t2]) interactionRange=top.rcut[t1][t2];
        }
    }
    //interactionRange = 0.1;
    double box_size = inp["box_size"];
    double density = nParticles * (M_PI/4) / (box_size*box_size);                           // particle number density density
    printf("density = %5.4f\n", density);

    unsigned int output_every = inp["output_every"];
    unsigned int sweeps = inp["sweeps"];
    std::string last_conf = "last_conf.xyz";
    //if inp.contains("last_conf") last_conf = inp["last_conf"];
    std::string trajectory = "trajectory.xyz";
    //if inp.contains("trajectory") trajectory = inp["trajectory"];

    //exit(1);

    //double baseLength;                              // base length of simulation box

    // Data structures.
    std::vector<Particle> particles(nParticles);    // particle container
    CellList cells;                                 // cell list
    bool isIsotropic[nParticles];                   // whether the potential of each particle is isotropic
    // Resize particle container.
    particles.resize(nParticles);

    // // Work out base length of simulation box (particle diameter is one).
    // if (dimension == 2) baseLength = std::pow((nParticles*M_PI)/(4.0*density), 1.0/2.0);
    // else baseLength = std::pow((nParticles*M_PI)/(6.0*density), 1.0/3.0);
    std::vector<double> boxSizeVec;
    for (unsigned int i=0;i<dimension;i++)
        boxSizeVec.push_back(box_size);
    // Initialise simulation box object.
    Box box(boxSizeVec);

    // Initialise input/output class,
    InputOutput io;

    // Create VMD script.
    io.vmdScript(boxSizeVec);

    // Initialise cell list.
    cells.setDimension(dimension);
    cells.initialise(box.boxSize, interactionRange);

    // Initialise the patchy disc model.
    GaussianPatchyDisc patchyDisc(box, particles, cells, top, maxInteractions, interactionEnergy, interactionRange);

    // Initialise random number generator.
    MersenneTwister rng;

    // Initialise particle initialisation object.
    Initialise initialise;

    // Generate a random particle configuration.
    //std::cout << top.Ni[0] <<", "<< top.Ni[1] << std::endl;
    initialise.random(top, particles, cells, box, rng, false);

    // double p1[2] = {1., 1.};
    // double p2[2] = {2.4, 1.};
    // double o1[2] = {1., 0.};
    // double o2[2] = {1., 0.};
    // std::cout << patchyDisc.computePairEnergy(0, p1, o1, 1, p2, o2) << std::endl;


    // Initialise data structures needed by the VMMC class.
    double coordinates[dimension*nParticles];
    double orientations[dimension*nParticles];

    // Copy particle coordinates and orientations into C-style arrays.
    for (unsigned int i=0;i<nParticles;i++)
    {
        for (unsigned int j=0;j<dimension;j++)
        {
            coordinates[dimension*i + j] = particles[i].position[j];
            orientations[dimension*i + j] = particles[i].orientation[j];
        }

        // Set all particles as anisotropic.
        isIsotropic[i] = false;
    }

    // Initialise the VMMC callback functions.
    using namespace std::placeholders;
    vmmc::CallbackFunctions callbacks;
#ifndef ISOTROPIC
    callbacks.energyCallback =
        std::bind(&GaussianPatchyDisc::computeEnergy, patchyDisc, _1, _2, _3);
    callbacks.pairEnergyCallback =
        std::bind(&GaussianPatchyDisc::computePairEnergy, patchyDisc, _1, _2, _3, _4, _5, _6);
    callbacks.interactionsCallback =
        std::bind(&GaussianPatchyDisc::computeInteractions, patchyDisc, _1, _2, _3, _4);
    callbacks.postMoveCallback =
        std::bind(&GaussianPatchyDisc::applyPostMoveUpdates, patchyDisc, _1, _2, _3);
#else
    callbacks.energyCallback =
        std::bind(&GaussianPatchyDisc::computeEnergy, patchyDisc, _1, _2);
    callbacks.pairEnergyCallback =
        std::bind(&GaussianPatchyDisc::computePairEnergy, patchyDisc, _1, _2, _3, _4);
    callbacks.interactionsCallback =
        std::bind(&GaussianPatchyDisc::computeInteractions, patchyDisc, _1, _2, _3);
    callbacks.postMoveCallback =
        std::bind(&GaussianPatchyDisc::applyPostMoveUpdates, patchyDisc, _1, _2);
#endif

    // Initialise VMMC object.
    vmmc::VMMC vmmc(nParticles, dimension, coordinates, orientations,
        0.15, 0.2, 0.5, 0.5, maxInteractions, &boxSizeVec[0], isIsotropic, false, callbacks);

    // save init conf and print init energy (is not needed if restarting)
    long long int starting_step = 0;
    long long int curr_step;
    printf("sweeps = %10lld, energy = %5.4f\n", starting_step, patchyDisc.getEnergy());
    io.appendXyzTrajectory(trajectory, starting_step, box, particles, true, true);
    // Execute the simulation.
    for(curr_step = starting_step; curr_step < sweeps && !stop; curr_step++)
    {
        vmmc += nParticles;
        if(curr_step > 0 && (curr_step % (output_every)) == 0)
        {
            io.appendXyzTrajectory(trajectory, curr_step, box, particles, false, true);
            printf("sweeps = %10lld, energy = %5.4f\n", curr_step, patchyDisc.getEnergy());
        }
    }
    printf("Writing the last conf to `%s`.", last_conf.c_str());
    io.appendXyzTrajectory(last_conf, curr_step, box, particles, true, true);
    std::cout << "\nComplete!\n";

    // We're done!
    return (EXIT_SUCCESS);
}
