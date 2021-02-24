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

  Adapted from PatchyDisc.cpp
*/

#include <cmath>
#include <cstdlib>
#include <iostream>

#include "Box.h"
#include "CellList.h"
#include "Particle.h"
#include "Top.h"
#include "GaussianPatchyDisc.h"

#ifndef M_PI
    #define M_PI 3.1415926535897932384626433832795
#endif

GaussianPatchyDisc::GaussianPatchyDisc(
    Box& box_,
    std::vector<Particle>& particles_,
    CellList& cells_,
    Top& top_,
    unsigned int maxInteractions_,
    double interactionEnergy_,
    double interactionRange_) :
    Model(box_, particles_, cells_, maxInteractions_, interactionEnergy_, interactionRange_),
    top(top_)
{
#ifdef ISOTROPIC
    std::cerr << "[ERROR] GaussianPatchyDisc: Cannot be used with isotropic VMMC library!\n";
    exit(EXIT_FAILURE);
#endif

    // std::cout << "nPatches = "<<nPatches<<"\n";
    // std::cout << "(GaussianPatchyDisk) nPatchesTypes = "<<top.nTypes<<"\n";
    // std::cout << "(GaussianPatchyDisk) epsilon[0][0] = "<<top.epsilon[0][0]<<"\n";

    // Check dimensionality.
    if (box.dimension != 2)
    {
        std::cerr << "[ERROR] GaussianPatchyDisc: Model only valid in two dimensions!\n";
        exit(EXIT_FAILURE);
    }

    for (unsigned int i=0;i<particles.size();i++)
    {
        idx2type.push_back(particles[i].type);
    }

    // // Populate lookup tables.
    // unsigned int np;
    // // Resize rotation matrix lookup tables.
    // cosTheta.resize(particles.size());
    // sinTheta.resize(particles.size());
    // for (unsigned int i=0;i<particles.size();i++)
    // {
    //     idx2type.push_back(particles[i].type);
    //     // Work out the angle between patches.
    //     np = top.nPatches[particles[i].type];
    //     patchSeparation = 2.0*M_PI/np;
    //     cosTheta[i].resize(np);
    //     sinTheta[i].resize(np);
    //     for (unsigned int j=0;j<np;j++)
    //     {
    //         cosTheta[i][j] = cos(i*patchSeparation);
    //         sinTheta[i][j] = sin(i*patchSeparation);
    //     }
    // }
}

double GaussianPatchyDisc::computePairEnergy(unsigned int particle1, const double* position1,
    const double* orientation1, unsigned int particle2, const double* position2, const double* orientation2)
{
    unsigned int t1 = idx2type[particle1];
    unsigned int t2 = idx2type[particle2];

    // Separation vector.
    std::vector<double> sep(2);
    std::vector<double> patch_sep(2);

    // Calculate disc separation.
    sep[0] = position1[0] - position2[0];
    sep[1] = position1[1] - position2[1];

    // Enforce minimum image.
    box.minimumImage(sep);

    // Calculate squared norm of vector.
    double normSqd = sep[0]*sep[0] + sep[1]*sep[1];

    //Discs overlap.
    if (normSqd < 1) //top.rcut_sq[t1][t2]) 
    {
        return INF;
    }

    // Total interaction energy sum.
    double energy = 0;

    // Test interactions between all patch top.
    for (unsigned int i=0;i<top.nPatches[t1];i++)
    {
        // Compute position of patch i on first disc.
        std::vector<double> coord1(2);
        coord1[0] = position1[0] + 0.5*(orientation1[0]*top.cosPatchAngles[t1][i] - orientation1[1]*top.sinPatchAngles[t1][i]);
        coord1[1] = position1[1] + 0.5*(orientation1[0]*top.sinPatchAngles[t1][i] + orientation1[1]*top.cosPatchAngles[t1][i]);

        // Enforce periodic boundaries.
        box.periodicBoundaries(coord1);

        for (unsigned int j=0;j<top.nPatches[t2];j++)
        {
            // Compute position of patch j on second disc.
            std::vector<double> coord2(2);
            coord2[0] = position2[0] + 0.5*(orientation2[0]*top.cosPatchAngles[t2][i] - orientation2[1]*top.sinPatchAngles[t2][i]);
            coord2[1] = position2[1] + 0.5*(orientation2[0]*top.sinPatchAngles[t2][i] + orientation2[1]*top.cosPatchAngles[t2][i]);

            // Enforce periodic boundaries.
            box.periodicBoundaries(coord2);

            // Calculate patch separation.
            patch_sep[0] = coord1[0] - coord2[0];
            patch_sep[1] = coord1[1] - coord2[1];

            // Enforce minimum image.
            box.minimumImage(patch_sep);

            // Calculate squared norm of vector.
            normSqd = patch_sep[0]*patch_sep[0] + patch_sep[1]*patch_sep[1];

            // Patches interact.
            //std::cout << normSqd << std::endl;
            if (normSqd < squaredCutOffDistance)
                energy -= top.epsilon[t1][t2];

        }
    }

    return energy;
}

unsigned int GaussianPatchyDisc::computeInteractions(unsigned int particle,
    const double* position, const double* orientation, unsigned int* interactions)
{
    // Interaction counter.
    unsigned int nInteractions = 0;

    // Check all neighbouring cells including same cell.
    for (unsigned int i=0;i<cells.getNeighbours();i++)
    {
        // Cell index.
        unsigned int cell = cells[particles[particle].cell].neighbours[i];

        // Check all particles within cell.
        for (unsigned int j=0;j<cells[cell].tally;j++)
        {
            // Index of neighbouring particle.
            unsigned int neighbour = cells[cell].particles[j];

            // Make sure the particles are different.
            if (neighbour != particle)
            {
                // Calculate pair energy.
                double energy = computePairEnergy(particle, position, orientation,
                                neighbour, &particles[neighbour].position[0],
                                &particles[neighbour].orientation[0]);

                // Particles interact.
                if (energy < 0 )
                {
                    if (nInteractions == maxInteractions)
                    {
                        std::cerr << "[ERROR] GaussianPatchyDisc: Maximum number of interactions exceeded!\n";
                        exit(EXIT_FAILURE);
                    }

                    interactions[nInteractions] = neighbour;
                    nInteractions++;
                }
            }
        }
    }

    return nInteractions;
}
