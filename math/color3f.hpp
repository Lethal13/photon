#ifndef MATH_COLOR3F_HPP
#define MATH_COLOR3F_HPP

#include "common_math.h"

struct color3f_t
{
	color3f_t()
		: r(0.0), g(0.0), b(0.0)
	{
	}

	color3f_t(float r, float g, float b)
		: r(r), g(g), b(b)
	{
		MATH_ASSERT(!has_NaNs());
	}

	color3f_t(float scalar)
		: r(scalar), g(scalar), b(scalar)
	{
		MATH_ASSERT(!has_NaNs());
	}

	color3f_t(const color3f_t& c)
		: r(c.r), g(c.g), b(c.b)
	{
	}

	color3f_t& operator=(const color3f_t& c)
	{
		this->r = c.r;
		this->g = c.g;
		this->b = c.b;

		return *this;
	}

	color3f_t& operator+=(const color3f_t& c)
	{
		this->r += c.r;
		this->g += c.g;
		this->b += c.b;
		MATH_ASSERT(!has_NaNs());
		
		return *this;
	}

	color3f_t& operator*=(float s)
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

	float r, g, b;
};

inline color3f_t operator+(const color3f_t& a, const color3f_t& b)
{
	return color3f_t(a.r + b.r, a.g + b.g, a.b + b.b);
}

inline color3f_t operator*(const color3f_t& a, const color3f_t& b)
{
	return color3f_t(a.r * b.r, a.g * b.g, a.b * b.b);
}

inline color3f_t operator*(float a, const color3f_t& c)
{
	return color3f_t(a * c.r, a * c.g, a * c.b);
}

inline color3f_t operator*(const color3f_t& c, float a)
{
	return color3f_t(a * c.r, a * c.g, a * c.b);
}

#endif
