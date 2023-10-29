#ifndef PHOTON_MATH_COMMON_MATH_HPP
#define PHOTON_MATH_COMMON_MATH_HPP

#define MAX(A, B) ((A) > (B) ? (A) : (B))
#define MIN(A, B) ((A) < (B) ? (A) : (B))

#define FLOAT_INFINITY std::numeric_limits<float>::infinity();
#define PI32 3.141592653589793238462643f 

// TODO: Replace std library.
#include <random>
#include <cmath>

inline float degrees_to_radians(float degrees)
{
    return degrees * PI32 / 180.0f;
}

inline float random_float()
{
	static std::uniform_real_distribution<float> distribution(0.0f, 1.0f);
	static std::mt19937 generator;
	return distribution(generator);
}

inline float random_float(float min, float max)
{
	return min + (max-min) * random_float();
}

inline double random_double()
{
	static std::uniform_real_distribution<double> distribution(0.0f, 1.0f);
	static std::mt19937 generator;
	return distribution(generator);
}

inline double random_double(double min, double max)
{
	return min + (max - min) * random_double();
}

#endif
