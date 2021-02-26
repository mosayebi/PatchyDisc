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

    //init all variables and look-up tables
    for (unsigned int i=0;i<particles.size();i++)
    {
        idx2type.push_back(particles[i].type);
    }
    rcut_sq.resize(top.nTypes, std::vector<double>(top.nTypes));
    lj_shift.resize(top.nTypes, std::vector<double>(top.nTypes));
    twoSigmapSq.resize(top.nTypes, std::vector<double>(top.nTypes));
    twoSigmapSq.resize(top.nTypes, std::vector<double>(top.nTypes));
    sigmaSq.resize(top.nTypes, std::vector<double>(top.nTypes));
    for (unsigned int i=0;i<top.nTypes;i++)
    {
        for (unsigned int j=0;j<top.nTypes;j++)
        {
            rcut_sq[i][j] = top.rcut[i][j] * top.rcut[i][j];
            lj_shift[i][j] = pow(top.sigma[i][j]/top.rcut[i][j], 12) - pow(top.sigma[i][j]/top.rcut[i][j], 6);
            lj_shift[i][j] *= -4 * top.epsilon[i][j];
            twoSigmapSq[i][j] = 2 * top.sigma_p[i][j] * top.sigma_p[i][j];
            sigmaSq[i][j] = top.sigma[i][j] * top.sigma[i][j];
         }
    }
// std::cout << "init done!" << std::endl;
// std::cout << std::endl;
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
    sep[0] = position2[0] - position1[0];
    sep[1] = position2[1] - position1[1];

    // Enforce minimum image.
    box.minimumImage(sep);

    // Calculate squared norm of vector.
    double normSqd = sep[0]*sep[0] + sep[1]*sep[1];

    // Particles interact.
    double energyLJ;
    if (normSqd < sigmaSq[t1][t2])
    {
        double r2Inv = sigmaSq[t1][t2] / normSqd;
        double r6Inv = r2Inv*r2Inv*r2Inv;
        energyLJ = 4.0*top.epsilon[t1][t2]*((r6Inv*r6Inv) - r6Inv) + lj_shift[t1][t2];
        return energyLJ;
    }
    else if (normSqd < rcut_sq[t1][t2])
    {
        double r2Inv = sigmaSq[t1][t2] / normSqd;
        double r6Inv = r2Inv*r2Inv*r2Inv;
        energyLJ = 4.0*top.epsilon[t1][t2]*((r6Inv*r6Inv) - r6Inv) + lj_shift[t1][t2];
    }
    else return 0;

    // Calculate the angular modulation of the LJ interaction.
    // examine all patch pairs and take the one that has the maximum value.
    double max_modulation = 0.0;
    double r12Angle, r21Angle, modulation;
    // TODO: fast atan2 approx:https://gist.github.com/volkansalma/2972237
    // TODO: atan2 table
    r12Angle = atan2(sep[1], sep[0]);
    r21Angle = r12Angle + M_PI;

    for (unsigned int i=0;i<top.nPatches[t1];i++)
    {
        double p1Angle = atan2(orientation1[1], orientation1[0]);
        //std::cout<< "orientation1 angle " << p1Angle/M_PI*180 << std::endl;
        p1Angle += top.patchAngles[t1][i];
        //std::cout << "1 (patch "<<i<<"): o1=("<<orientation1[0]<<","<<orientation1[1]<<"),   angle1="<<p1Angle/M_PI*180<<std::endl;

        // std::vector<double> coord1(2);
        // coord1[0] = position1[0] + 0.5*(orientation1[0]*cosPatchAngles[t1][i] - orientation1[1]*sinPatchAngles[t1][i]);
        // coord1[1] = position1[1] + 0.5*(orientation1[0]*sinPatchAngles[t1][i] + orientation1[1]*cosPatchAngles[t1][i]);
        // // Enforce periodic boundaries.
        // box.periodicBoundaries(coord1);

        for (unsigned int j=0;j<top.nPatches[t2];j++)
        {
            double p2Angle = atan2(orientation2[1], orientation2[0]);
            //std::cout<< "orientation2 angle " << p1Angle/M_PI*180 << std::endl;
            p2Angle += top.patchAngles[t2][j];
            //std::cout << "2 (patch "<<j<<"): o2=(" <<orientation2[0]<<","<<orientation2[1]<<"),   angle2="<<p2Angle/M_PI*180<<std::endl;
            p1Angle = p1Angle - r12Angle;
            p2Angle = p2Angle - r21Angle;
            while (p1Angle > M_PI) p1Angle -= 2*M_PI;
            while (p2Angle < -M_PI) p2Angle += 2*M_PI;
            //std::cout << p1Angle << std::endl;

            modulation = (p1Angle*p1Angle)/twoSigmapSq[t1][t2] + (p2Angle*p2Angle)/twoSigmapSq[t1][t2];
            modulation = exp(-modulation);
            if (max_modulation < modulation) max_modulation = modulation;

            // // Compute position of patch j on second disc.
            // std::vector<double> coord2(2);
            // coord2[0] = position2[0] + 0.5*(orientation2[0]*top.cosPatchAngles[t2][i] - orientation2[1]*top.sinPatchAngles[t2][i]);
            // coord2[1] = position2[1] + 0.5*(orientation2[0]*top.sinPatchAngles[t2][i] + orientation2[1]*top.cosPatchAngles[t2][i]);
            // // Enforce periodic boundaries.
            // box.periodicBoundaries(coord2);
            // // Calculate patch separation.
            // patch_sep[0] = coord1[0] - coord2[0];
            // patch_sep[1] = coord1[1] - coord2[1];
            // // Enforce minimum image.
            // box.minimumImage(patch_sep);
            // // Calculate squared norm of vector.
            // normSqd = patch_sep[0]*patch_sep[0] + patch_sep[1]*patch_sep[1];
            // // Patches interact.
            // //std::cout << normSqd << std::endl;
            // if (normSqd < squaredCutOffDistance)
            //     energy -= top.epsilon[t1][t2];

        }
    }
    //std::cout<<particle1<<"<>"<<particle2<< energyLJ << "," << max_modulation << std::endl;
    return energyLJ*max_modulation;
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
                std::vector<double> sep(box.dimension);

                // Compute separation.
                for (unsigned int k=0;k<box.dimension;k++)
                    sep[k] = position[k] - particles[neighbour].position[k];

                // Enforce minimum image.
                box.minimumImage(sep);

                double normSqd = 0;

                // Calculate squared norm of vector.
                for (unsigned int k=0;k<box.dimension;k++)
                    normSqd += sep[k]*sep[k];

                // Particles interact.

                if (normSqd < rcut_sq[idx2type[particle]][idx2type[neighbour]])
                {
                    if (nInteractions == maxInteractions)
                    {
                        std::cerr << "[ERROR] Model: Maximum number of interactions exceeded!\n";
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