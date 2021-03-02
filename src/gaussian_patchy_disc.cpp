#ifdef ISOTROPIC
#error gaussian_patchy_disc.cpp cannot be linked to isotropic VMMC library!
#endif

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include "System.h"


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
    std::string inputFileName = argv[1];

	//Here we handle a few SIG* signals: whenever we intercept one of these signals
	//the program will quit the main loop and die as gracefully as possible.
	signal(SIGTERM, gbl_terminate);
	signal(SIGABRT, gbl_terminate);
	signal(SIGINT, gbl_terminate);
	signal(SIGUSR2, gbl_terminate);

	//Load the input file and Initialise
    inp = ReadJsonFromFile(inputFileName);
    parseInitialisationBlock();
    parseTopologyBlock();
    parseSimulationParamBlock();


    std::vector<Particle> particles(nParticles);    // particle container
    bool isIsotropic[nParticles];                   // whether the potential of each particle is isotropic
    particles.resize(nParticles);                   // Resize particle container.
    Box box(boxSizeVec);
    printf("# Number density is %5.4f\n", nParticles/box.Volume);

    // Initialise cell list.
    CellList cells;
    cells.setDimension(dimension);
    cells.initialise(box.boxSize, interactionRange);

    // Initialise random number generator.
    MersenneTwister rng;
    rng.setSeed(seed);
    std::cout << "# Random seed set to "<< rng.getSeed()<< std::endl;

    // Initialise input/output class,
    InputOutput io;
    // Create VMD script.
    io.vmdScript(boxSizeVec);

    // Initialise the patchy disc model.
    GaussianPatchyDisc patchyDisc(box, particles, cells, top, maxInteractions, interactionEnergy, interactionRange);
    // if (interaction == "GaussianPatchyc")
    // {
    //     patchyDisc = GaussianPatchyDisc(box, particles, cells, top, maxInteractions, interactionEnergy, interactionRange);
    //     hasFiniteRepulsion = true;
    // }
    // else if (interaction == "GaussianPatchyDiscHR")
    // {
    //     patchyDisc = GaussianPatchyDiscHR(box, particles, cells, top, maxInteractions, interactionEnergy, interactionRange);
    //     hasFiniteRepulsion = false;
    // }
    // else if (interaction == "GaussianPatchyDiscHRSW")
    // {
    //     patchyDisc = GaussianPatchyDiscHRSW(box, particles, cells, top, maxInteractions, interactionEnergy, interactionRange);
    //     hasFiniteRepulsion = false;
    // }
    // else
    // {
    //     std::cerr << "[ERROR] Invalid interaction ("<<interaction<<"!\n";
    //     exit(EXIT_FAILURE);
    // }


    // Initialise particle initialisation object.
    Initialise initialise;

    if (init_mode == "from_random_conf")
    {
        //Generate a random particle configuration.
        initialise.random(top, particles, cells, box, rng, false);
    }
    else if (init_mode == "from_init_conf")
    {
        loadInitialConfiguration(init_conf, box, particles, cells);
    }

    // // For debug
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
        0.15, 0.2, 0.5, 0.5, maxInteractions, &boxSizeVec[0], isIsotropic, hasFiniteRepulsion, callbacks);
    // Initialise single particle move object.
    SingleParticleMove MC(rng, &patchyDisc, 0.2, 0.1, 0.5, false);

    // append init_conf details depending on restart_step_counter
    if (restart_step_counter)
    {
        starting_step = 0;
        bool clearFile = true;
        io.appendXyzTrajectory(trajectory, starting_step, box, particles, clearFile, true);
        io.appendLog(log_file, starting_step, patchyDisc.getEnergy(), clearFile);
    }
    else
    {
        bool clearFile = false;
        if (starting_step % output_every == 0)
        {
            io.appendXyzTrajectory(trajectory, starting_step, box, particles, clearFile, true);
            io.appendLog(log_file, starting_step, patchyDisc.getEnergy(), clearFile);
        }
    }

    // Execute the simulation.
    for(curr_step = starting_step; curr_step < sweeps && !stop; curr_step++)
    {
        for (unsigned int i=0; i<nParticles/2; i++)
        {
            vmmc ++;
            MC ++;
        }
        if(curr_step > 0 && (curr_step % (output_every)) == 0)
        {
            io.appendXyzTrajectory(trajectory, curr_step, box, particles, false, true);
            io.appendLog(log_file, curr_step, patchyDisc.getEnergy(), false);
        }
    }
    printf("Writing the last conf to `%s`.", last_conf.c_str());
    io.appendXyzTrajectory(last_conf, curr_step, box, particles, true, true);
    std::cout << "\nComplete!\n";

    // We're done!
    return (EXIT_SUCCESS);
}
