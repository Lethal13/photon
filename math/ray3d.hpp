#ifndef MATH_RAY3D_HPP
#define MATH_RAY3D_HPP

#include "vec3d.hpp"
#include "point3d.hpp"

struct ray3d_t
{
	ray3d_t() = default;

	ray3d_t(const point3d_t& o, const vec3d_t& d)
		: origin(o), direction(d)
	{
	}

	point3d_t at(double t) const
	{
		return (origin + t * direction);
	}

	ray3d_t& operator=(const ray3d_t& r)
	{
		this->origin = r.origin;
		this->direction = r.direction;

		return *this;
	}

	point3d_t origin;
	vec3d_t direction;
};

#endif
