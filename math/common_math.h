#ifndef MATH_COMMON_MATH_H
#define MATH_COMMON_MATH_H

#include <stdint.h>
#include <math.h>
#include <limits>
#include <random>

#define MATH_ASSERT(expression) if(!(expression)) {*(volatile int*)0 = 0;}

#define MAX(A, B) ((A) > (B) ? (A) : (B))
#define MIN(A, B) ((A) < (B) ? (A) : (B))

#define PI32_SP 3.141592653589793238462643f // Single-precision float representation of π
#define PI32 3.141592653589793238462643 // Double-precision float representation of π
#define FLOAT_INFINITY std::numeric_limits<float>::infinity()
#define DOUBLE_INFINITY std::numeric_limits<double>::infinity()

inline float degrees_to_radians(float degrees)
{
    return degrees * PI32_SP / 180.0f;
}

inline double degrees_to_radians(double degrees)
{
    return degrees * PI32 / 180.0;
}

inline float random_float()
{
	static std::uniform_real_distribution<float> distribution(0.0f, 1.0f);
	static std::mt19937 generator;
	return distribution(generator);
}

inline float random_float(float min, float max)
{
	return min + (max - min) * random_float();
}

inline double random_double()
{
	static std::uniform_real_distribution<double> distribution(0.0, 1.0);
	static std::mt19937 generator;
	return distribution(generator);
}

inline double random_double(double min, double max)
{
	return min + (max - min) * random_double();
}

inline float clamp(float x, float min, float max)
{
    if (x < min) return min;
    if (x > max) return max;
    return x;
}

// Source: https://entropymine.com/imageworsener/srgbformula/
inline float linear_to_srgb(float value)
{
    if(value < 0.0f)
    {
        value = 0.0f;
    }
    
    if(value > 1.0f)
    {
        value = 1.0f;
    }
    
    float result = value * 12.92f;
    if(value > 0.0031308f)
    {
        result = 1.055f * std::pow(value, 1.0f/2.4f) - 0.055f;
    }
    
    return result;
}

#endif
