#include <cmath>
#include <cstdlib>
#include <iostream>

#include "Box.h"
#include "CellList.h"
#include "Particle.h"
#include "Top.h"
#include "Utils.h"
#include "GaussianPatchyDisc.h"
#include "GaussianPatchyDiscHRSW.h"


#ifndef M_PI
    #define M_PI 3.1415926535897932384626433832795
#endif

GaussianPatchyDiscHRSW::GaussianPatchyDiscHRSW(
    Box& box_,
    std::vector<Particle>& particles_,
    CellList& cells_,
    Top& top_,
    unsigned int maxInteractions_,
    double interactionEnergy_,
    double interactionRange_) :
    GaussianPatchyDisc(box_, particles_, cells_, top_, maxInteractions_, interactionEnergy_, interactionRange_)
{
    hasFiniteRepulsion = false;
    std::cout << "# Initialised GaussianPatchyDiscHRSW interaction; Derived from GaussianPatchyDisc" << std::endl;
}

double GaussianPatchyDiscHRSW::computePairEnergy(unsigned int particle1, const double* position1,
    const double* orientation1, unsigned int particle2, const double* position2, const double* orientation2)
{
    // Early exis if either of particles is a ghost
    if ( idx2ifGhost[particle1] || idx2ifGhost[particle2] ) return INF;

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
        return INF;
    }
    else if (normSqd < rcut_sq[t1][t2])
    {
        energyLJ = - top.epsilon[t1][t2];
    }
    else return 0;

    // Calculate the angular modulation of the LJ interaction.
    // examine all patch pairs and take the one that has the maximum value.
    double max_modulation, p1Angle, p2Angle;
    double r12Angle, r21Angle, modulation;
    max_modulation = 0.0;
    // to use the 10-20 % faster you can use atan2_approximation function. 
    // TODO: atan2 table
    r12Angle = atan2(sep[1], sep[0]);
    r21Angle = r12Angle + M_PI;

    for (unsigned int i=0;i<top.nPatches[t1];i++)
    {
        p1Angle = atan2(orientation1[1], orientation1[0]);
        p1Angle += top.patchAngles[t1][i];
        p1Angle = p1Angle - r12Angle;

        // std::vector<double> coord1(2);
        // coord1[0] = position1[0] + 0.5*(orientation1[0]*cosPatchAngles[t1][i] - orientation1[1]*sinPatchAngles[t1][i]);
        // coord1[1] = position1[1] + 0.5*(orientation1[0]*sinPatchAngles[t1][i] + orientation1[1]*cosPatchAngles[t1][i]);
        // // Enforce periodic boundaries.
        // box.periodicBoundaries(coord1);

        for (unsigned int j=0;j<top.nPatches[t2];j++)
        {
            //double p2Angle;
            p2Angle = atan2(orientation2[1], orientation2[0]);
            p2Angle += top.patchAngles[t2][j];
            p2Angle = p2Angle - r21Angle;
            while (p1Angle > M_PI) p1Angle -= 2*M_PI;
            while (p2Angle < -M_PI) p2Angle += 2*M_PI;

            modulation = (p1Angle*p1Angle)/twoSigmapSq[t1][t2] + (p2Angle*p2Angle)/twoSigmapSq[t1][t2];
            modulation = exp(-modulation);
            // if (modulation>=max_x){modulation = 0;} else {
            //    modulation = expTable[ (unsigned int) (modulation/dx) ];}
            if (max_modulation < modulation) max_modulation = modulation;

        }
    }
    return energyLJ*max_modulation;
}