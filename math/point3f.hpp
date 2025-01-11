#ifndef MATH_POINT3F_HPP
#define MATH_POINT3F_HPP

#include "common_math.h"
#include "vec3f.hpp"

struct point3f_t
{
	point3f_t()
		: x(0.0), y(0.0), z(0.0)
	{
	}

	point3f_t(float x, float y, float z)
		: x(x), y(y), z(z)
	{
		MATH_ASSERT(!has_NaNs());
	}

	point3f_t(float scalar)
		: x(scalar), y(scalar), z(scalar)
	{
		MATH_ASSERT(!has_NaNs());
	}

	point3f_t(const point3f_t& p)
		: x(p.x), y(p.y), z(p.z)
	{
		MATH_ASSERT(!has_NaNs());
	}

	point3f_t operator-() {
		return point3f_t(-x, -y, -z);
	}

	point3f_t& operator=(const point3f_t& p)
	{
		this->x = p.x;
		this->y = p.y;
		this->z = p.z;

		return *this;
	}

	point3f_t& operator+=(const vec3f_t& v)
	{
		this->x += v.x;
		this->y += v.y;
		this->z += v.z;

		return *this;
	}

	point3f_t& operator-=(const vec3f_t& v)
	{
		this->x -= v.x;
		this->y -= v.y;
		this->z -= v.z;

		return *this;
	}

    const float operator[](uint32_t i) const
    {
        MATH_ASSERT(i >= 0 && i < 3);
        return (&x)[i];
    }

    float& operator[](uint32_t i)
    {
        MATH_ASSERT(i >= 0 && i < 3);
        return (&x)[i];
    }

    int32_t has_NaNs() const
    {
        return !(x == x) || !(y == y) || !(z == z);
    }

	float x, y, z;
};

inline point3f_t operator+(const point3f_t& p, const vec3f_t& v)
{
    return point3f_t(p.x + v.x, p.y + v.y, p.z + v.z);
}

inline point3f_t operator+(const vec3f_t& p, const point3f_t& v)
{
    return point3f_t(v.x + p.x, v.y + p.y, v.z + p.z);
}

inline point3f_t operator-(const point3f_t& p, const vec3f_t& v)
{
    return point3f_t(p.x - v.x, p.y - v.y, p.z - v.z);
}

inline point3f_t operator-(const vec3f_t& p, const point3f_t& v)
{
    return point3f_t(v.x - p.x, v.y - p.y, v.z - p.z);
}

inline vec3f_t operator-(const point3f_t& a, const point3f_t& b)
{
	return vec3f_t(a.x - b.x, a.y - b.y, a.z - b.z);
}

inline float min_component(const point3f_t& p)
{
    return MIN(p.x, MIN(p.y, p.z));
}

inline float max_component(const point3f_t& p)
{
    return MAX(p.x, MAX(p.y, p.z));
}

#endif
