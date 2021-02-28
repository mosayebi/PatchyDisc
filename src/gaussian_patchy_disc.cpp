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
#include <random>



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
    std::cerr << "[ERROR] Failed to parse json!" << std::endl;
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


    // Initialisation
    unsigned int nParticles = 0;
    unsigned int nTypes = 1;
    Top top;
    double box_size;

    std::string name="initialisation";
    unsigned int seed = std::random_device{}();
    if (inp[name].contains("seed")) seed = inp[name]["seed"];
    // figure out initialisation mode
    std::string init_mode = inp[name]["mode"];
    if (inp[name].contains(init_mode))
    {
        if (init_mode == "random_conf")
        {
            nTypes =inp[name][init_mode]["types"];
            top.setSize(nTypes);
            //init Ni
            for (unsigned int i=0;i<nTypes;i++){
                unsigned int t1 = inp[name][init_mode]["particle_numbers"][i]["type"];
                top.Ni[t1] = inp[name][init_mode]["particle_numbers"][i]["N"];
                nParticles += top.Ni[t1];
            }

            box_size = inp[name][init_mode]["box_size"];
            // double density = nParticles * (M_PI/4) / (box_size*box_size);   // area fraction assuming all sigma=1.
            printf("# Number density is %5.4f\n", nParticles/box_size/box_size );
        }
        else if (init_mode == "init_conf")
        {
            //TODO
        }
    }
    else
    {
        std::cerr << "[ERROR] Initialisation mode is not defined!\n";
        exit(EXIT_FAILURE);
    }

    // Initialise topology
    // patches
    unsigned int t1, t2;
    name="topology";
    for (unsigned int i=0;i<nTypes;i++){
        t1 = inp[name]["patches"][i]["type"];
        top.nPatches[t1] = inp[name]["patches"][i]["nPatches"];
        //std::cout << inp["patches"][i]["angles"] << std::endl;
        for (unsigned int j=0;j<top.nPatches[t1];j++){
            double dummy = inp[name]["patches"][i]["angles"][j];
            top.patchAngles[t1].push_back(M_PI * (dummy / 180.0));
        }
    }
    //init pair_coeff
    unsigned int counter = 0;
    for (unsigned int i=0;i<nTypes;i++){
        for (unsigned int j=0;j<nTypes;j++){
            t1 = inp[name]["pair_coeff"][counter]["type1"];
            t2 = inp[name]["pair_coeff"][counter]["type2"];
            top.epsilon[t1][t2] = inp[name]["pair_coeff"][counter]["epsilon"];
            top.delta[t1][t2] = inp[name]["pair_coeff"][counter]["delta"];
            top.sigma[t1][t2] = inp[name]["pair_coeff"][counter]["sigma"];
            top.sigma_p[t1][t2] = inp[name]["pair_coeff"][counter]["sigma_p"];
            top.rcut[t1][t2] = inp[name]["pair_coeff"][counter]["rcut"];
            counter += 1;

            if (interactionEnergy<top.epsilon[t1][t2]) interactionEnergy=top.epsilon[t1][t2];
            if (interactionRange<top.rcut[t1][t2]) interactionRange=top.rcut[t1][t2];
        }
    }


    unsigned int output_every = 1000;
    name = "output";
    if (inp[name].contains("output_every")) output_every = inp[name]["output_every"];
    unsigned int sweeps = inp[name]["sweeps"];
    std::string last_conf = "last_conf.xyz";
    if (inp[name].contains("last_conf")) last_conf = inp[name]["last_conf"];
    std::string trajectory = "trajectory.xyz";
    if (inp[name].contains("trajectory")) trajectory = inp[name]["trajectory"];

    //exit(1);

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
    rng.setSeed(seed);
    std::cout << "# Random seed set to "<< rng.getSeed()<< std::endl;

    // Initialise particle initialisation object.
    Initialise initialise;

    // Generate a random particle configuration.
    //std::cout << top.Ni[0] <<", "<< top.Ni[1] << std::endl;
    initialise.random(top, particles, cells, box, rng, false);

    // for (double t=0; t<=360; t+=1)
    // {
    //     double p1[2] = {10., 10.};
    //     double p2[2] = {10.84, 10.84};
    //     double o1[2] = {cos(30./180*M_PI), sin(30./180*M_PI)};
    //     double o2[2] = {cos(t/180*M_PI), sin(t/180*M_PI)};
    //     double E = patchyDisc.computePairEnergy(0, p1, o1, 1, p2, o2);
    //     std::cout << t << " " << E <<"     "<< o2[0]<< " "<<o2[1] << std::endl << std::endl;
    // }
    // exit(1);


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
    vmmc::VMMC vmmc(rng, nParticles, dimension, coordinates, orientations,
        0.15, 0.2, 0.5, 0.5, maxInteractions, &boxSizeVec[0], isIsotropic, false, callbacks);

    // Initialise single particle move object.
    // Todo rng should be an arguman.
    //SingleParticleMove MC(&patchyDisc, 0.2, 0.1, 0.5, false);

    // save init conf and print init energy (is not needed if restarting)
    long long int starting_step = 0;
    long long int curr_step;
    printf("sweeps = %10lld, energy = %5.4f\n", starting_step, patchyDisc.getEnergy());
    io.appendXyzTrajectory(trajectory, starting_step, box, particles, true, true);
    // Execute the simulation.
    for(curr_step = starting_step; curr_step < sweeps && !stop; curr_step++)
    {
        for (unsigned int i=0; i<nParticles/2; i++)
        {
            vmmc ++;
            //MC ++;
        }
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
