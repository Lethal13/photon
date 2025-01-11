#ifndef MATH_VEC3F_HPP
#define MATH_VEC3F_HPP

#include "common_math.h"

struct vec3f_t
{
	vec3f_t()
		: x(0.0f), y(0.0f), z(0.0f)
	{
	}

	vec3f_t(float x, float y, float z)
		: x(x), y(y), z(z)
	{
		MATH_ASSERT(!has_NaNs());
	}

	vec3f_t(float scalar)
		: x(scalar), y(scalar), z(scalar)
	{
		MATH_ASSERT(!has_NaNs());
	}

	vec3f_t(const vec3f_t& v)
		: x(v.x), y(v.y), z(v.z)
	{
		MATH_ASSERT(!has_NaNs());
	}

	vec3f_t operator-() const
	{
		return vec3f_t(-x, -y , -z);
	}

	vec3f_t& operator=(const vec3f_t& v)
	{
		this->x = v.x;
		this->y = v.y;
		this->z = v.z;
		MATH_ASSERT(!has_NaNs());

		return *this;
	}

	vec3f_t& operator+=(const vec3f_t& v)
	{
		this->x += v.x;
		this->y += v.y;
		this->z += v.z;
		MATH_ASSERT(!has_NaNs());

		return *this;
	}

	vec3f_t& operator-=(const vec3f_t& v)
	{
		this->x -= v.x;
		this->y -= v.y;
		this->z -= v.z;
		MATH_ASSERT(!has_NaNs());

		return *this;
	}

	vec3f_t& operator*=(float scalar)
	{
		this->x *= scalar;
		this->y *= scalar;
		this->z *= scalar;
		MATH_ASSERT(!has_NaNs());

		return *this;
	}

	vec3f_t& operator/=(float scalar)
	{
		// TODO: Add asset scalar != 0.0f
		scalar = 1.0f / scalar;
		this->x *= scalar;
		this->y *= scalar;
		this->z *= scalar;

		return *this;
	}


	int32_t has_NaNs() const
	{
		return !(x == x) || !(y == y) || !(z == z);
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

    uint32_t near_zero() const
    {
        // Return true if the vector is close to zero in all dimensions.
        constexpr double s = 1e-8;
        return (fabs(x) < s) && (fabs(y) < s) && (fabs(z) < s);
    }

	float x, y, z;
};

inline vec3f_t operator+(const vec3f_t& a, const vec3f_t& b)
{
	return vec3f_t(a.x + b.x, a.y + b.y, a.z + b.z);
}

inline vec3f_t operator-(const vec3f_t& a, const vec3f_t& b)
{
	return vec3f_t(a.x - b.x, a.y - b.y, a.z - b.z);
}

inline vec3f_t operator*(const vec3f_t& v, float scalar)
{
	return vec3f_t(v.x * scalar, v.y * scalar, v.z * scalar);
}

inline vec3f_t operator*(float scalar, const vec3f_t& v)
{
    return vec3f_t(v.x * scalar, v.y * scalar, v.z * scalar);
}

inline vec3f_t operator/(const vec3f_t& v, float scalar)
{
	// TODO: Add asset scalar != 0.0f
    scalar = 1.0f / scalar;
    return vec3f_t(v.x * scalar, v.y * scalar, v.z * scalar);
}

inline float magnitude(const vec3f_t& v)
{
	return sqrtf(v.x * v.x + v.y * v.y + v.z * v.z);
}

inline float magnitude_squared(const vec3f_t& v)
{
	return v.x * v.x + v.y * v.y + v.z * v.z;
}

inline vec3f_t normalize(const vec3f_t& v)
{
    return v / magnitude(v);
}

inline float dot(const vec3f_t& a, const vec3f_t& b)
{
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

inline vec3f_t cross(const vec3f_t& a, const vec3f_t& b)
{
    return vec3f_t(a.y * b.z - a.z * b.y,
				   a.z * b.x - a.x * b.z,
				   a.x * b.y - a.y * b.x);
}

inline float min_component(const vec3f_t& v)
{
    return MIN(v.x, MIN(v.y, v.z));
}

inline float max_component(const vec3f_t& v)
{
    return MAX(v.x, MAX(v.y, v.z));
}

inline vec3f_t permute(const vec3f_t& v, uint32_t x, uint32_t y, uint32_t z)
{
    return vec3f_t(v[x], v[y], v[z]);
}

inline vec3f_t project(const vec3f_t& a, const vec3f_t& b)
{
    return (b * dot(a, b) / dot(b, b));
}

inline vec3f_t reject(const vec3f_t& a, const vec3f_t& b)
{
    return (a - b * (dot(a, b) / dot(b, b)));
}

inline vec3f_t random_vec3f()
{
	return vec3f_t(random_float(), random_float(), random_float());
}

inline vec3f_t random(float min, float max)
{
	return vec3f_t(random_float(min, max), random_float(min, max), random_float(min, max));
}

inline vec3f_t random_unit_vec3f()
{
    while (true)
	{
        vec3f_t v = random(-1.0f, 1.0f);
        float lensq = magnitude_squared(v);
        if (1e-12 < lensq && lensq <= 1.0f)
            return v / sqrtf(lensq);
    }
}

inline vec3f_t random_on_hemisphere(const vec3f_t& normal)
{
    vec3f_t on_unit_sphere = random_unit_vec3f();
    if (dot(on_unit_sphere, normal) > 0.0) // In the same hemisphere as the normal
        return on_unit_sphere;
    else
        return -on_unit_sphere;
}

inline vec3f_t reflect(const vec3f_t& v, const vec3f_t& n)
{
    return v - 2 * dot(v, n) * n;
}

inline vec3f_t refract(const vec3f_t& uv, const vec3f_t& n, float etai_over_etat)
{
    float cos_theta = fmin(dot(-uv, n), 1.0f);
    vec3f_t r_out_perp =  etai_over_etat * (uv + cos_theta*n);
    vec3f_t r_out_parallel = -sqrtf(fabs(1.0f - magnitude_squared(r_out_perp))) * n;
    return r_out_perp + r_out_parallel;
}

inline vec3f_t random_vec3f_in_unit_disk()
{
    while (true)
	{
        vec3f_t v = vec3f_t(random_float(-1.0f, 1.0f), random_float(-1.0f, 1.0f), 0.0f);
        if (magnitude_squared(v) < 1)
            return v;
    }
}

#endif
