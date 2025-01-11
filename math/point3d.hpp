#ifndef MATH_POINT3D_HPP
#define MATH_POINT3D_HPP

#include "common_math.h"
#include "vec3d.hpp"

struct point3d_t
{
	point3d_t()
		: x(0.0), y(0.0), z(0.0)
	{
	}

	point3d_t(double x, double y, double z)
		: x(x), y(y), z(z)
	{
		MATH_ASSERT(!has_NaNs());
	}

	point3d_t(double scalar)
		: x(scalar), y(scalar), z(scalar)
	{
		MATH_ASSERT(!has_NaNs());
	}

	point3d_t(const point3d_t& p)
		: x(p.x), y(p.y), z(p.z)
	{
		MATH_ASSERT(!has_NaNs());
	}

	point3d_t operator-() {
		return point3d_t(-x, -y, -z);
	}

	point3d_t& operator=(const point3d_t& p)
	{
		this->x = p.x;
		this->y = p.y;
		this->z = p.z;

		return *this;
	}

	point3d_t& operator+=(const vec3d_t& v)
	{
		this->x += v.x;
		this->y += v.y;
		this->z += v.z;

		return *this;
	}

	point3d_t& operator-=(const vec3d_t& v)
	{
		this->x -= v.x;
		this->y -= v.y;
		this->z -= v.z;

		return *this;
	}

    const double operator[](uint32_t i) const
    {
        MATH_ASSERT(i >= 0 && i < 3);
        return (&x)[i];
    }

    double& operator[](uint32_t i)
    {
        MATH_ASSERT(i >= 0 && i < 3);
        return (&x)[i];
    }

    int32_t has_NaNs() const
    {
        return !(x == x) || !(y == y) || !(z == z);
    }

	double x, y, z;
};

inline point3d_t operator+(const point3d_t& p, const vec3d_t& v)
{
    return point3d_t(p.x + v.x, p.y + v.y, p.z + v.z);
}

inline point3d_t operator+(const vec3d_t& p, const point3d_t& v)
{
    return point3d_t(v.x + p.x, v.y + p.y, v.z + p.z);
}

inline point3d_t operator-(const point3d_t& p, const vec3d_t& v)
{
    return point3d_t(p.x - v.x, p.y - v.y, p.z - v.z);
}

inline point3d_t operator-(const vec3d_t& p, const point3d_t& v)
{
    return point3d_t(v.x - p.x, v.y - p.y, v.z - p.z);
}

inline vec3d_t operator-(const point3d_t& a, const point3d_t& b)
{
	return vec3d_t(a.x - b.x, a.y - b.y, a.z - b.z);
}

inline double min_component(const point3d_t& p)
{
    return MIN(p.x, MIN(p.y, p.z));
}

inline double max_component(const point3d_t& p)
{
    return MAX(p.x, MAX(p.y, p.z));
}

#endif
