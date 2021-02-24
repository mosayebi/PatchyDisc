#include "Top.h"

Top::Top(unsigned int nTypes_) :
    nTypes(nTypes_),
    Ni(std::vector<unsigned int>(nTypes)),
    nPatches(std::vector<unsigned int>(nTypes)),
    patchAngles(std::vector<std::vector<double>>(nTypes)),
    cosPatchAngles(std::vector<std::vector<double>>(nTypes)),
    sinPatchAngles(std::vector<std::vector<double>>(nTypes)),
    epsilon(std::vector<std::vector<double>>(nTypes, std::vector<double>(nTypes))),
    delta(std::vector<std::vector<double>>(nTypes, std::vector<double>(nTypes))),
    sigma(std::vector<std::vector<double>>(nTypes, std::vector<double>(nTypes))),
    sigma_p(std::vector<std::vector<double>>(nTypes, std::vector<double>(nTypes))),
    rcut(std::vector<std::vector<double>>(nTypes, std::vector<double>(nTypes))),
    rcut_sq(std::vector<std::vector<double>>(nTypes, std::vector<double>(nTypes)))
{
}