#ifndef PHOTON_CORE_COLOR_HPP
#define PHOTON_CORE_COLOR_HPP

struct color_t
{
	color_t() = default;
	color_t(double red, double green, double blue)
		:r(red), g(green), b(blue)
	{
	}

	color_t(double scalar)
		:r(scalar), g(scalar), b(scalar)
	{
	}

	color_t(const color_t& c)
		:r(c.r), g(c.g), b(c.b)
	{
	}

	color_t& operator=(const color_t& c)
	{
		this->r = c.r;
		this->g = c.g;
		this->b = c.b;

		return *this;
	}

	color_t& operator+=(const color_t& c)
	{
		this->r += c.r;
		this->g += c.g;
		this->b += c.b;

		return *this;
	}

	color_t& operator*=(const color_t& c)
	{
		this->r *= c.r;
		this->g *= c.g;
		this->b *= c.b;

		return *this;
	}

	color_t& operator*=(double s)
	{
		this->r *= s;
		this->g *= s;
		this->b *= s;

		return *this;
	}

	color_t& operator/=(double s)
	{
		PHOTON_ASSERT(s != 0.0);
		s = 1.0 / s;
		
		this->r *= s;
		this->g *= s;
		this->b *= s;

		return *this;
	}

	double r, g, b;
};

inline color_t operator+(const color_t& a, const color_t& b)
{
	return color_t(a.r + b.r, a.g + b.g, a.b + b.b);
}

inline color_t operator*(const color_t& a, const color_t& b)
{
	return color_t(a.r * b.r, a.g * b.g, a.b * b.b);
}

inline color_t operator*(double scalar, const color_t& c)
{
	return color_t(scalar * c.r, scalar * c.g, scalar * c.b);
}

#endif
