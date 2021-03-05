#ifdef ISOTROPIC
#error gaussian_patchy_disc.cpp cannot be linked to isotropic VMMC library!
#endif

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include "System.h"


int main(int argc, char** argv)
{
    std::cout << "# PatchyDisc compiled on " << COMPILED_ON << std::endl;
    std::cout << "# GIT information: branch '"<< GIT_BRANCH <<"', commit hash '" << GIT_COMMIT_HASH << "'" <<std::endl << std::endl;

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

	//Load the input file and Initialise (inp is a json object)
    inp = ReadJsonFromFile(inputFileName);
    parseInputFile();

    std::vector<Particle> particles(nParticles);    // particle container
    bool isIsotropic[nParticles];                   // whether the potential of each particle is isotropic
    particles.resize(nParticles);                   // Resize particle container.
    Box box(boxSizeVec);
    printf("# Number of species %7d\n", nTypes);
    printf("# Number of particles %7d\n", nParticles);
    printf("# Number density is %7.4f\n", nParticles/box.Volume);

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
    hasFiniteRepulsion = true;
    if (interaction != "GaussianPatchyDisc") std::cout<< "# [NotImplementedError] Currently, changing the interaction is not allowed from the input_file. You can do so by editing the source code.\n# Using the default GaussianPatchyDisc interaction" << std::endl;
    // if (interaction == "GaussianPatchyDisc")
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



    if (init_mode == "from_random_conf")
    {
        // Initialise particle initialisation object.
        Initialise initialise;
        //Generate a random particle configuration.
        initialise.random(top, particles, cells, box, rng, false);
    }
    else if (init_mode == "from_init_conf")
    {
        loadInitialConfiguration(init_conf, box, particles, cells);
        std::cout<< "# Starting configuration is loaded from '"<< last_conf<<"'. Setting current step to "<< starting_step <<"."<< std::endl;
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
    // SingleParticleMove MC(rng, &patchyDisc, 0.2, 0.1, 0.5, false);

    // Initialisation depending on the restart_step_counter
    if (restart_step_counter)
    {
        starting_step = 0;
        bool clearFile = true;
        std::cout<< "# Resetting the step counter to zero. Owerwiting '"<< trajectory << "' and '"<< log_file<<"'.\n" << std::endl;
        io.appendXyzTrajectory(trajectory, starting_step, box, particles, clearFile, true);
        io.appendLog(log_file, starting_step, patchyDisc.getEnergy(), clearFile);
    }
    else
    {
        bool clearFile = false;
        std::cout<< "# Resuming the simulation; Appending to '"<< trajectory << "' and '"<< log_file<<"'.\n"<< std::endl;
        if (starting_step % output_every == 0)
        {
            io.appendXyzTrajectory(trajectory, starting_step, box, particles, clearFile, true);
            io.appendLog(log_file, starting_step, patchyDisc.getEnergy(), clearFile);
        }
    }

    // Execute the simulation.
    for(curr_step = starting_step; curr_step < sweeps && !stop; curr_step++)
    {
        for (unsigned int i=0; i<nParticles; i++)
        {
            vmmc ++;
            // MC ++;   // both moves cannot be used at the same time. BUG?
        }
        if(curr_step > 0 && (curr_step % (output_every)) == 0)
        {
            io.appendXyzTrajectory(trajectory, curr_step, box, particles, false, true);
            io.appendLog(log_file, curr_step, patchyDisc.getEnergy(), false);
        }
    }
    printf("\nWriting the last conf to `%s`.", last_conf.c_str());
    io.appendXyzTrajectory(last_conf, curr_step, box, particles, true, true);
    std::cout << "\nComplete!\n";

    // We're done!
    return (EXIT_SUCCESS);
}
