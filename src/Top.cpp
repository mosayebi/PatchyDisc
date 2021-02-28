#include "Top.h"

Top::Top()
{
}

void Top::setSize(unsigned int nTypes_)
{

    nTypes = nTypes_;

    epsilon.resize(nTypes, std::vector<double>(nTypes));
    delta.resize(nTypes, std::vector<double>(nTypes));
    sigma.resize(nTypes, std::vector<double>(nTypes));
    sigma_p.resize(nTypes, std::vector<double>(nTypes));
    rcut.resize(nTypes, std::vector<double>(nTypes));

    Ni.resize(nTypes);
    nPatches.resize(nTypes);
    patchAngles.resize(nTypes);
}
