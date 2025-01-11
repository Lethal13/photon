#ifndef MATH_VEC3D_HPP
#define MATH_VEC3D_HPP

#include "common_math.h"

struct vec3d_t
{
	vec3d_t()
		: x(0.0), y(0.0), z(0.0)
	{
	}

	vec3d_t(double x, double y, double z)
		: x(x), y(y), z(z)
	{
		MATH_ASSERT(!has_NaNs());
	}

	vec3d_t(double scalar)
		: x(scalar), y(scalar), z(scalar)
	{
		MATH_ASSERT(!has_NaNs());
	}

	vec3d_t(const vec3d_t& v)
		: x(v.x), y(v.y), z(v.z)
	{
		MATH_ASSERT(!has_NaNs());
	}

	vec3d_t operator-() const
	{
		return vec3d_t(-x, -y , -z);
	}

	vec3d_t& operator=(const vec3d_t& v)
	{
		this->x = v.x;
		this->y = v.y;
		this->z = v.z;
		MATH_ASSERT(!has_NaNs());

		return *this;
	}

	vec3d_t& operator+=(const vec3d_t& v)
	{
		this->x += v.x;
		this->y += v.y;
		this->z += v.z;
		MATH_ASSERT(!has_NaNs());

		return *this;
	}

	vec3d_t& operator-=(const vec3d_t& v)
	{
		this->x -= v.x;
		this->y -= v.y;
		this->z -= v.z;
		MATH_ASSERT(!has_NaNs());

		return *this;
	}

	vec3d_t& operator*=(double scalar)
	{
		this->x *= scalar;
		this->y *= scalar;
		this->z *= scalar;
		MATH_ASSERT(!has_NaNs());

		return *this;
	}

	vec3d_t& operator/=(double scalar)
	{
		// TODO: Add asset scalar != 0.0
		scalar = 1.0 / scalar;
		this->x *= scalar;
		this->y *= scalar;
		this->z *= scalar;

		return *this;
	}

	int32_t has_NaNs() const
	{
		return !(x == x) || !(y == y) || !(z == z);
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

    uint32_t near_zero() const
    {
        // Return true if the vector is close to zero in all dimensions.
        constexpr double s = 1e-8;
        return (fabs(x) < s) && (fabs(y) < s) && (fabs(z) < s);
    }
	
	double x, y, z;
};

inline vec3d_t operator+(const vec3d_t& a, const vec3d_t& b)
{
	return vec3d_t(a.x + b.x, a.y + b.y, a.z + b.z);
}

inline vec3d_t operator-(const vec3d_t& a, const vec3d_t& b)
{
	return vec3d_t(a.x - b.x, a.y - b.y, a.z - b.z);
}

inline vec3d_t operator*(const vec3d_t& v, double scalar)
{
	return vec3d_t(v.x * scalar, v.y * scalar, v.z * scalar);
}

inline vec3d_t operator*(double scalar, const vec3d_t& v)
{
	return vec3d_t(v.x * scalar, v.y * scalar, v.z * scalar);
}

inline vec3d_t operator/(const vec3d_t& v, double scalar)
{
	// TODO: Add asset scalar != 0.0
    scalar = 1.0 / scalar;
    return vec3d_t(v.x * scalar, v.y * scalar, v.z * scalar);
}

inline double magnitude(const vec3d_t& v)
{
	return sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
}

inline double magnitude_squared(const vec3d_t& v)
{
	return v.x * v.x + v.y * v.y + v.z * v.z;
}

inline vec3d_t normalize(const vec3d_t& v)
{
    return v / magnitude(v);
}

inline double dot(const vec3d_t& a, const vec3d_t& b)
{
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

inline vec3d_t cross(const vec3d_t& a, const vec3d_t& b)
{
    return vec3d_t(a.y * b.z - a.z * b.y,
				   a.z * b.x - a.x * b.z,
				   a.x * b.y - a.y * b.x);
}

inline double min_component(const vec3d_t& v)
{
    return MIN(v.x, MIN(v.y, v.z));
}

inline double max_component(const vec3d_t& v)
{
    return MAX(v.x, MAX(v.y, v.z));
}

inline vec3d_t permute(const vec3d_t& v, uint32_t x, uint32_t y, uint32_t z)
{
    return vec3d_t(v[x], v[y], v[z]);
}

inline vec3d_t project(const vec3d_t& a, const vec3d_t& b)
{
    return (b * dot(a, b) / dot(b, b));
}

inline vec3d_t reject(const vec3d_t& a, const vec3d_t& b)
{
    return (a - b * (dot(a, b) / dot(b, b)));
}

inline vec3d_t random_vec3d()
{
	return vec3d_t(random_double(), random_double(), random_double());
}

inline vec3d_t random(double min, double max)
{
	return vec3d_t(random_double(min, max), random_double(min, max), random_double(min, max));
}

inline vec3d_t random_unit_vec3d()
{
    while (true)
	{
        vec3d_t v = random(-1.0, 1.0);
        double lensq = magnitude_squared(v);
        if (1e-160 < lensq && lensq <= 1.0)
            return v / sqrt(lensq);
    }
}

inline vec3d_t random_on_hemisphere(const vec3d_t& normal)
{
    vec3d_t on_unit_sphere = random_unit_vec3d();
    if (dot(on_unit_sphere, normal) > 0.0) // In the same hemisphere as the normal
        return on_unit_sphere;
    else
        return -on_unit_sphere;
}

inline vec3d_t reflect(const vec3d_t& v, const vec3d_t& n)
{
    return v - 2 * dot(v, n) * n;
}

inline vec3d_t refract(const vec3d_t& uv, const vec3d_t& n, double etai_over_etat)
{
    double cos_theta = fmin(dot(-uv, n), 1.0);
    vec3d_t r_out_perp =  etai_over_etat * (uv + cos_theta*n);
    vec3d_t r_out_parallel = -sqrt(fabs(1.0 - magnitude_squared(r_out_perp))) * n;
    return r_out_perp + r_out_parallel;
}

inline vec3d_t random_vec3d_in_unit_disk()
{
    while (true)
	{
        vec3d_t v = vec3d_t(random_double(-1.0, 1.0), random_double(-1.0, 1.0), 0.0);
        if (magnitude_squared(v) < 1)
            return v;
    }
}

#endif
