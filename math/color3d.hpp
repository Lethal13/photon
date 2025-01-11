#ifndef MATH_COLOR3D_HPP
#define MATH_COLOR3D_HPP

#include "common_math.h"

struct color3d_t
{
	color3d_t()
		: r(0.0), g(0.0), b(0.0)
	{
	}

	color3d_t(double r, double g, double b)
		: r(r), g(g), b(b)
	{
		MATH_ASSERT(!has_NaNs());
	}

	color3d_t(double scalar)
		: r(scalar), g(scalar), b(scalar)
	{
		MATH_ASSERT(!has_NaNs());
	}

	color3d_t(const color3d_t& c)
		: r(c.r), g(c.g), b(c.b)
	{
	}

	color3d_t& operator=(const color3d_t& c)
	{
		this->r = c.r;
		this->g = c.g;
		this->b = c.b;

		return *this;
	}

	color3d_t& operator+=(const color3d_t& c)
	{
		this->r += c.r;
		this->g += c.g;
		this->b += c.b;
		MATH_ASSERT(!has_NaNs());
		
		return *this;
	}

	color3d_t& operator*=(double s)
	{
		this->r *= s;
		this->g *= s;
		this->b *= s;
		MATH_ASSERT(!has_NaNs());

		return *this;
	}

    int32_t has_NaNs() const
    {
        return !(r == r) || !(g == g) || !(b == b);
    }

	double r, g, b;
};


inline color3d_t operator+(const color3d_t& a, const color3d_t& b)
{
	return color3d_t(a.r + b.r, a.g + b.g, a.b + b.b);
}

inline color3d_t operator*(const color3d_t& a, const color3d_t& b)
{
	return color3d_t(a.r * b.r, a.g * b.g, a.b * b.b);
}

inline color3d_t operator*(double a, const color3d_t& c)
{
	return color3d_t(a * c.r, a * c.g, a * c.b);
}

inline color3d_t operator*(const color3d_t& c, double a)
{
	return color3d_t(a * c.r, a * c.g, a * c.b);
}

#endif
