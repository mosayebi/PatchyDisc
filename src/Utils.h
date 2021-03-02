
#ifndef _UTILS_H
#define _UTILS_H

#include <nlohmann/json.hpp>
#include <fstream>
#include <string>


double atan2_approximation(double, double);
nlohmann::json ReadJsonFromFile(std::string);


#endif