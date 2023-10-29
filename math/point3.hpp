#ifndef PHOTON_MATH_POINT3_HPP
#define PHOTON_MATH_POINT3_HPP

#include "common_math.hpp"
#include "vec3.hpp"

struct point3_t
{
	point3_t() = default;

	point3_t(float a, float b, float c)
		:x(a), y(b), z(c)
	{
		PHOTON_ASSERT(!has_NaNs());
	}

	point3_t(float scalar)
		:x(scalar), y(scalar), z(scalar)
	{
		PHOTON_ASSERT(!has_NaNs());
	}

	point3_t(const point3_t& p)
		:x(p.x), y(p.y), z(p.z)
	{
	}

	point3_t operator-()
	{
		return point3_t(-x, -y, -z);
	}

	point3_t& operator=(const point3_t& p)
	{
		this->x = p.x;
		this->y = p.y;
		this->z = p.z;

		return *this;
	}

	point3_t& operator+=(const vec3_t& v)
	{
		this->x += v.x;
		this->y += v.y;
		this->z += v.z;
		PHOTON_ASSERT(!has_NaNs());
		
		return *this;
	}

	point3_t& operator-=(const vec3_t& v)
	{
		this->x -= v.x;
		this->y -= v.y;
		this->z -= v.z;
		PHOTON_ASSERT(!has_NaNs());

		return *this;
	}

    const float operator[](uint32_t i) const
    {
        PHOTON_ASSERT(i >= 0 && i < 3);
        return (&x)[i];
    }

    float& operator[](int32_t i)
    {
        PHOTON_ASSERT(i >= 0 && i < 3);
        return (&x)[i];
    }

    const float operator[](int32_t i) const
    {
        PHOTON_ASSERT(i >= 0 && i < 3);
        return (&x)[i];
    }

    int32_t has_NaNs() const
    {
        return !(x == x) || !(y == y) || !(z == z);
    }

    uint32_t near_zero() const
    {
        // Return true if the vector is close to zero in all dimensions.
        constexpr auto s = 1e-8;
        return (fabs(x) < s) && (fabs(y) < s) && (fabs(z) < s);
    }

	float x, y, z;
};

inline point3_t operator+(const point3_t& p, const vec3_t& v)
{
    return point3_t(p.x + v.x, p.y + v.y, p.z + v.z);
}

inline point3_t operator+(const vec3_t& p, const point3_t& v)
{
    return point3_t(v.x + p.x, v.y + p.y, v.z + p.z);
}

inline point3_t operator-(const point3_t& p, const vec3_t& v)
{
    return point3_t(p.x - v.x, p.y - v.y, p.z - v.z);
}

inline point3_t operator-(const vec3_t& p, const point3_t& v)
{
    return point3_t(v.x - p.x, v.y - p.y, v.z - p.z);
}

inline vec3_t operator-(const point3_t& a, const point3_t& b)
{
	return vec3_t(a.x - b.x, a.y - b.y, a.z - b.z);
}

inline float min_component(const point3_t& p)
{
    return MIN(p.x, MIN(p.y, p.z));
}

inline float max_component(const point3_t& p)
{
    return MAX(p.x, MAX(p.y, p.z));
}

#endif
