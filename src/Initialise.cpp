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

#include <cstdlib>
#include <iostream>
//#include <numeric>


#include "Box.h"
#include "CellList.h"
#include "Particle.h"
#include "Initialise.h"
#include "MersenneTwister.h"
#include "Top.h"

Initialise::Initialise()
{
}

void Initialise::random(unsigned int nParticles, unsigned int maxParticles, Top& top, std::vector<Particle>& particles, CellList& cells, Box& box, MersenneTwister& rng, bool isSpherocylinder)
{
    std::cout<< "# Generating random starting configuration..."<< std::endl;

    if (isSpherocylinder && (box.dimension != 3)){
        std::cerr << "[ERROR] Initialise: Spherocylindrical boundary only valid for three dimensional simulation box!\n";
        exit(EXIT_FAILURE);
    }

    // Copy box dimensions.
    boxSize = box.boxSize;
    unsigned int sum = 0;
    sum = std::accumulate(top.Ni.begin(), top.Ni.end(), 0);
    if (sum != nParticles){
        std::cerr << "[ERROR] Initialise: Invalid number of particles! ("<< nParticles <<" != "<< sum << ")\n";
        exit(EXIT_FAILURE);
    }
    if (maxParticles != particles.size()){
        std::cerr << "[ERROR] Initialise: Invalid max number of particles! ("<< particles.size() <<" != "<< maxParticles << ")\n";
        exit(EXIT_FAILURE);
    }
    unsigned int idx=0;
    // Current number of attempted particle insertions.
    unsigned int nTrials = 0;
    for (unsigned int k=0;k<top.Ni.size();k++)
    {
        //std::cout << k <<": "<< top.Ni[k]<<std::endl;
        for (unsigned i=0;i<top.Ni[k];i++)
        {

            // Whether particle overlaps.
            bool isOverlap = true;

            // Temporary vector.
            std::vector<double> vec(box.dimension);

            particles[idx].ghost = false;

            // Set particle index.
            particles[idx].index = idx;

            // Set particle type.
            particles[idx].type = k;

            // Keep trying to insert particle until there is no overlap.
            while (isOverlap)
            {
                nTrials++;

                // Generate a random position.
                for (unsigned int j=0;j<box.dimension;j++)
                    vec[j] = rng()*box.boxSize[j];

                particles[idx].position = vec;

                // Generate a random orientation.
                for (unsigned int j=0;j<box.dimension;j++)
                    vec[j] = rng.normal();

                // Calculate vector norm.
                double norm = 0;
                for (unsigned int j=0;j<box.dimension;j++)
                    norm += vec[j]*vec[j];
                norm = sqrt(norm);

                // Convert orientation to a unit vector.
                for (unsigned int j=0;j<box.dimension;j++)
                    vec[j] /= norm;

                particles[idx].orientation = vec;

                // Calculate the particle's cell index.
                particles[idx].cell = cells.getCell(particles[idx]);

                // Enforce spherocylindrical boundary.
                if (isSpherocylinder)
                {
                    // Make sure particle lies within the spherocylinder.
    #ifndef ISOTROPIC
                    if (!outsideSpherocylinder(idx, &particles[idx].position[0], &particles[idx].orientation[0]))
    #else
                    if (!outsideSpherocylinder(idx, &particles[idx].position[0]))
    #endif
                    {
                        // See if there is any overlap between particles.
                        isOverlap = checkOverlap(particles[idx], particles, cells, box, top);
                    }
                    else isOverlap = true;
                }
                else
                {
                    // See if there is any overlap between particles.
                    isOverlap = checkOverlap(particles[idx], particles, cells, box, top);
                }

                // Check trial limit isn't exceeded.
                if (nTrials == MAX_TRIALS)
                {
                    std::cerr << "[ERROR] Initialise: Maximum number of trial insertions reached. Only "<< idx <<" particles could be added.\n";
                    exit(EXIT_FAILURE);
                }
            }

            // Update cell list.
            cells.initCell(particles[idx].cell, particles[idx]);
            idx += 1;
        }
    }
    // create all ghost particles
    for (unsigned int i=nParticles; i<maxParticles; i++)
    {
        std::vector<double> vec(box.dimension);
        // Generate a random position.
        for (unsigned int j=0;j<box.dimension;j++)
            vec[j] = rng()*box.boxSize[j];
        particles[idx].position = vec;

        // Generate a random orientation.
        for (unsigned int j=0;j<box.dimension;j++)
            vec[j] = rng.normal();
        double norm = 0;
        for (unsigned int j=0;j<box.dimension;j++)
            norm += vec[j]*vec[j];
        norm = sqrt(norm);
        for (unsigned int j=0;j<box.dimension;j++)
            vec[j] /= norm;
        particles[idx].orientation = vec;
        particles[idx].cell = cells.getCell(particles[idx]);
        particles[idx].type = 0;
        particles[idx].ghost = true;
        particles[idx].index = true;
        // Update cell list.
        cells.initCell(particles[idx].cell, particles[idx]);
        idx += 1;
    }

std::cout<< "# Random starting configuration is generated after "<< nTrials <<" insertion trials."<< std::endl;
std::cout<< "# "<< maxParticles - nParticles << " ghost particles are added to the simulation box"<< std::endl;
}

#ifndef ISOTROPIC
bool Initialise::outsideSpherocylinder(unsigned int particle, const double* position, const double* orientation)
#else
bool Initialise::outsideSpherocylinder(unsigned int particle, const double* position)
#endif
{
    // Centre of sphere or circle.
    std::vector<double> centre(3);

    // Separation vector.
    std::vector<double> sep(3);

    // Squared radius of spherical cap (minus squared radius of particle).
    double radiusSqd = 0.25*(boxSize[0] - 1)*(boxSize[0] - 1);

    // Squared norm of separation vector.
    double normSqd = 0;

    // Initialise x and y coordinates of sphere centres.
    centre[0] = 0.5*boxSize[0];
    centre[1] = 0.5*boxSize[0];

    // Check whether particle lies in lower cap.
    if (position[2] < 0.5*boxSize[0])
    {
        centre[2] = 0.5*boxSize[0];

        // Calculate separation.
        for (unsigned int i=0;i<3;i++)
        {
            sep[i] = position[i] - centre[i];
            normSqd += sep[i]*sep[i];
        }

        // Particle lies outside of cap.
        if (normSqd > radiusSqd) return true;
    }
    else
    {
        // Check whether particle lies in upper cap.
        if (position[2] > (boxSize[2] - 0.5*boxSize[0]))
        {
            centre[2] = boxSize[2] - 0.5*boxSize[0];

            // Calculate separation.
            for (unsigned int i=0;i<3;i++)
            {
                sep[i] = position[i] - centre[i];
                normSqd += sep[i]*sep[i];
            }

            // Particle lies outside of cap.
            if (normSqd > radiusSqd) return true;
        }
        else
        {
            // Calculate separation.
            for (unsigned int i=0;i<2;i++)
            {
                sep[i] = position[i] - centre[i];
                normSqd += sep[i]*sep[i];
            }

            // Particle lies outside of cylinder.
            if (normSqd > radiusSqd) return true;
        }
    }

    // Inside spherocylinder.
    return false;
}

bool Initialise::checkOverlap(Particle& particle, std::vector<Particle>& particles, CellList& cells, Box& box, Top& top)
{
    unsigned int cell, neighbour;

    // Check all neighbouring cells including same cell.
    for (unsigned int i=0;i<cells.getNeighbours();i++)
    {
        cell = cells[particle.cell].neighbours[i];

        // Check all particles within cell.
        for (unsigned int j=0;j<cells[cell].tally;j++)
        {
            neighbour = cells[cell].particles[j];

            // Make sure particles are different.
            if (neighbour != particle.index)
            {
                // Particle separtion vector.
                std::vector<double> sep(box.dimension);

                // Compute separation.
                for (unsigned int k=0;k<box.dimension;k++)
                    sep[k] = particle.position[k] - particles[neighbour].position[k];

                // Compute minimum image.
                box.minimumImage(sep);

                double normSqd = 0;

                // Calculate squared norm of vector.
                for (unsigned int k=0;k<box.dimension;k++)
                    normSqd += sep[k]*sep[k];

                // Overlap if normSqd is less than particle diameter (box is scaled in diameter units).
                if (normSqd < top.sigma[particle.type][particles[neighbour].type]*top.sigma[particle.type][particles[neighbour].type]) return true; //TODO: use sigma_ij*sigma_ij instead
            }
        }
    }

    // If we get this far, no overlaps.
    return false;
}
