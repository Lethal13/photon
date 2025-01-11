#ifndef MATH_RAY3F_HPP
#define MATH_RAY3F_HPP

#include "vec3f.hpp"
#include "point3f.hpp"

struct ray3f_t
{
	ray3f_t() = default;

	ray3f_t(const point3f_t& o, const vec3f_t& d)
		: origin(o), direction(d)
	{
	}

	point3f_t at(float t) const
	{
		return origin + (t * direction);
	}

	ray3f_t& operator=(const ray3f_t& r)
	{
		this->origin = r.origin;
		this->direction = r.direction;

		return *this;
	}

	point3f_t origin;
	vec3f_t direction;
};

#endif
