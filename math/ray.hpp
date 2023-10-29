#ifndef PHOTON_MATH_RAY_HPP
#define PHOTON_MATH_RAY_HPP

#include "common_math.hpp"
#include "vec3.hpp"
#include "point3.hpp"

// NOTE: direction must be normalized.
struct ray_t
{
	ray_t() = default;

	ray_t(const point3_t& p, const vec3_t& v)
		:origin(p), direction(v)
	{
	}

	point3_t at(float t) const
	{
		return origin + t * direction;
	}

	// float t_min, t_max;
	point3_t origin;
	vec3_t direction;
};

#endif
