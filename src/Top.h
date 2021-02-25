#ifndef _TOP_H
#define _TOP_H

#include <vector>

/*! \file Top.h
    \brief A simple pair data type.
*/

//! Structure containing attributes for an individual particle.
struct Top
{
    //! Default constructor.
    Top(unsigned int nTypes_);

    unsigned int nTypes;                 //!< no Particle type.

    std::vector<unsigned int> Ni;

    std::vector<unsigned int> nPatches;
    std::vector<std::vector<double>> patchAngles;

    std::vector<std::vector<double>> epsilon;
    std::vector<std::vector<double>> delta;
    std::vector<std::vector<double>> sigma;
    std::vector<std::vector<double>> sigma_p;
    std::vector<std::vector<double>> rcut;
};

#endif  /* _TOP_H */