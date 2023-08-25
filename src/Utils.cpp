#include <cmath>
#include <iostream>
#include <csignal>

#include "Utils.h"

#ifndef M_PI
    #define M_PI 3.1415926535897932384626433832795
#endif


double atan2_approximation(double y, double x)
{
    //http://pubs.opengroup.org/onlinepubs/009695399/functions/atan2.html
    //Volkan SALMA

    const double ONEQTR_PI = M_PI / 4.0;
	const double THRQTR_PI = 3.0 * M_PI / 4.0;
	double r, angle;
	double abs_y = fabs(y) + 1e-10f;      // kludge to prevent 0/0 condition
	if ( x < 0.0f )
	{
		r = (x + abs_y) / (abs_y - x);
		angle = THRQTR_PI;
	}
	else
	{
		r = (x - abs_y) / (x + abs_y);
		angle = ONEQTR_PI;
	}
	angle += (0.1963f * r * r - 0.9817f) * r;
	if ( y < 0.0f )
		return( -angle );     // negate if in quad III or IV
	else
		return( angle );
}


//using json = nlohmann::json;
nlohmann::json ReadJsonFromFile(std::string fileName)
{
    try
		{
			return nlohmann::json::parse(std::ifstream{fileName, std::ios::in});
		}
		catch (nlohmann::json::parse_error& e) {
			std::cerr << "JSON parse exception : " << e.what() << std::endl;
		}
		catch (std::ifstream::failure& e) {
			std::cerr << "Stream exception : " << e.what() << std::endl;
		}
		catch (std::exception& e) {
			std::cerr << "Exception : " << e.what() << std::endl;
		}
		catch (...) {
			std::cerr << "Unknown error" << std::endl;
		}
		std::cerr << "[ERROR] Failed to parse json!" << std::endl;
		return {};
}

