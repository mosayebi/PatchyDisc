
/*
#include <cmath>
#include <iostream>
#include <string>
#include "parse_input.h"
#include <signal.h>
#include <csignal>
#include <nlohmann/json.hpp>
#include <fstream>

#ifndef M_PI
    #define M_PI 3.1415926535897932384626433832795
#endif

//using json = nlohmann::json;

nlohmann::json ReadJsonFromFile(std::string file_name) {
    try {
        return nlohmann::json::parse(std::ifstream{file_name, std::ios::in});
    }
		catch (json::parse_error& e) {
			std::cerr << "JSON parse exception : " << e.what() << std::endl;
		} catch (std::ifstream::failure& e) {
			std::cerr << "Stream exception : " << e.what() << std::endl;
		} catch (std::exception& e) {
			std::cerr << "Exception : " << e.what() << std::endl;
		} catch (...) {
			std::cerr << "Unknown error" << std::endl;
		}
    	return {};
}
**/
