#ifndef _SYSTEM_H
#define _SYSTEM_H

#include "Box.h"
#include "CellList.h"
#include "Initialise.h"
#include "InputOutput.h"
#include "Model.h"
#include "Particle.h"
#include "Top.h"

#include "Utils.h"

#include "SingleParticleMove.h"
#include "VMMC.h"

#include "GaussianPatchyDisc.h"
#include "GaussianPatchyDiscHR.h"
#include "GaussianPatchyDiscHRSW.h"

#include <random>
#include <vector>

#ifndef M_PI
    #define M_PI 3.1415926535897932384626433832795
#endif

int stop = 0;

// Simulation parameters.
unsigned int dimension = 2;                     // dimension of simulation box
double interactionRange = 1.0;                  // interaction range. #diameter of patch (in units of particle diameter)
double interactionEnergy = 1.0;                 // pair interaction energy scale (in units of kBT)
unsigned int maxInteractions = 36;              // maximum number of interactions per particle


nlohmann::json inp;                             // inputfile json
unsigned int nParticles;
unsigned int nTypes;
Top top;
std::vector<double> boxSizeVec (dimension);
std::string init_conf="init_conf.xyz";
std::string trajectory = "trajectory.xyz";
std::string last_conf = "last_conf.xyz";
std::string init_mode = "from_init_conf";

unsigned int output_every = 1000;
unsigned int sweeps = 10000;
long long int starting_step = 0;
long long int curr_step = 0;
bool restart_step_counter=true;

std::string interaction = "GaussainPatchyDisk";
bool hasFiniteRepulsion = true;

unsigned int seed = std::random_device{}();

void getParamsFromInitConf(std::string fileName)
{
    std::ifstream dataFile;
    std::vector<double> dummy (dimension);
    unsigned int type, max_type=0;

    dataFile.open(fileName.c_str());

    // Check that the file is valid.
    if (dataFile.good())
    {
        dataFile >> nParticles;
        dataFile >> starting_step;
        for (unsigned int j=0;j<dimension;j++) dataFile >> boxSizeVec[j];
        // find out nTypes
        for (unsigned int i=0;i<nParticles;i++)
        {
            dataFile >> type;
            if (type>max_type) max_type=type;
            for (unsigned int j=0;j<dimension;j++) dataFile >> dummy[j];
            for (unsigned int j=0;j<dimension;j++) dataFile >> dummy[j];
        }
        nTypes = max_type + 1;
    }
    else
    {
        std::cerr << "[ERROR] getParamsFromInitConf(): Invalid init_conf file ('"<< fileName <<"')!\n";
        exit(EXIT_FAILURE);
    }
    // Close file stream.
    dataFile.close();
}


void parseInitialisationBlock()
{
    std::string name="initialisation";
    //Zzstd::cout << name;
    // figure out initialisation mode

    if (inp[name].contains("mode")) init_mode.assign(inp[name]["mode"]);
    //init_mode::assign()
    if (inp[name].contains(init_mode))
    {
        if (init_mode == "from_random_conf")
        {
            nTypes =inp[name][init_mode]["types"];
            top.setSize(nTypes);
            for (unsigned int i=0;i<nTypes;i++){
                unsigned int t1 = inp[name][init_mode]["particle_numbers"][i]["type"];
                top.Ni[t1] = inp[name][init_mode]["particle_numbers"][i]["N"];
                nParticles += top.Ni[t1];
                //std::cout<< nParticles << std::endl;
            }
            for (unsigned int k=0; k<dimension; k++)
            {
                boxSizeVec[k] = inp[name][init_mode]["box"][k];
            }
        }
        else if (init_mode == "from_init_conf")
        {
            init_conf.assign(inp[name][init_mode]["init_conf"]);
            restart_step_counter = inp[name][init_mode]["restart_step_counter"];
            getParamsFromInitConf(init_conf);
            top.setSize(nTypes);
            // do not forget to initialise top.Ni when reading the init_conf. 
        }
    }
    else
    {
        std::cerr << "[ERROR] Initialisation mode is not defined!\n";
        exit(EXIT_FAILURE);
    }
}

void parseTopologyBlock()
{
    // Initialise topology
    // patches
    unsigned int t1, t2;
    std::string name="topology";
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
}

void parseSimulationParamBlock()
{
    std::string name = "simulation parameters";
    if (inp[name].contains("output_every")) output_every = inp[name]["output_every"];
    if (inp[name].contains("last_conf")) last_conf.assign(inp[name]["last_conf"]);
    if (inp[name].contains("trajectory")) trajectory.assign(inp[name]["trajectory"]);
    if (inp[name].contains("seed")) seed = inp[name]["seed"];
    if (inp.contains("interaction")) interaction.assign(inp["interaction"]);
    sweeps = inp[name]["sweeps"];
}



#endif
