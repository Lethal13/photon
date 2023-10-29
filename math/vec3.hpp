#ifndef PHOTON_MATH_VEC3_HPP
#define PHOTON_MATH_VEC3_HPP

#include "common_math.hpp"

struct vec3_t
{
	vec3_t() = default;

	vec3_t(float a, float b, float c)
		:x(a), y(b), z(c)
	{
		PHOTON_ASSERT(!has_NaNs());
	}

	vec3_t(float scalar)
		:x(scalar), y(scalar), z(scalar)
	{
		PHOTON_ASSERT(!has_NaNs());
	}

	vec3_t(const vec3_t& v)
		:x(v.x), y(v.y), z(v.z)
	{
	}

	vec3_t operator-()
	{
		return vec3_t(-x, -y, -z);
	}

	vec3_t& operator=(const vec3_t& v)
	{
		this->x = v.x;
		this->y = v.y;
		this->z = v.z;

		return *this;
	}

	vec3_t& operator+=(const vec3_t& v)
	{
		this->x += v.x;
		this->y += v.y;
		this->z += v.z;
		PHOTON_ASSERT(!has_NaNs());

		return *this;
	}

	vec3_t& operator-=(const vec3_t& v)
	{
		this->x -= v.x;
		this->y -= v.y;
		this->z -= v.z;
		PHOTON_ASSERT(!has_NaNs());

		return *this;
	}

	vec3_t& operator*=(float scalar)
	{
		this->x *= scalar;
		this->y *= scalar;
		this->z *= scalar;
		PHOTON_ASSERT(!has_NaNs());

		return *this;
	}

	vec3_t& operator/=(float scalar)
	{
        PHOTON_ASSERT(scalar != 0.0f);
		scalar = 1.0f / scalar;

		this->x *= scalar;
		this->y *= scalar;
		this->z *= scalar;

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

inline vec3_t operator+(const vec3_t& a, const vec3_t& b)
{
    return vec3_t(a.x + b.x, a.y + b.y, a.z + b.z);
}

inline vec3_t operator-(const vec3_t& a, const vec3_t& b)
{
    return vec3_t(a.x - b.x, a.y - b.y, a.z - b.z);
}

inline vec3_t operator*(const vec3_t& v, float scalar)
{
    return vec3_t(v.x * scalar, v.y * scalar,  v.z * scalar);
}

inline vec3_t operator*(float scalar, const vec3_t& v)
{
    return vec3_t(v.x * scalar, v.y * scalar, v.z * scalar);
}

inline vec3_t operator/(const vec3_t& v, float scalar)
{
    PHOTON_ASSERT(scalar != 0.0f);
    scalar = 1.0f / scalar;

    return vec3_t(v.x * scalar, v.y * scalar, v.z * scalar);
}

inline float magnitude(const vec3_t& v)
{
    return std::sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
}

inline float magnitude_squared(const vec3_t& v)
{
    return v.x * v.x + v.y * v.y + v.z * v.z;
}

inline vec3_t normalize(const vec3_t& v)
{
    return v / magnitude(v);
}

inline float dot(const vec3_t& a, const vec3_t& b)
{
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

inline vec3_t cross(const vec3_t& a, const vec3_t& b)
{
    return vec3_t(a.y * b.z - a.z * b.y,
                a.z * b.x - a.x * b.z,
                a.x * b.y - a.y * b.x);
}

inline float min_component(const vec3_t& v)
{
    return MIN(v.x, MIN(v.y, v.z));
}

inline float max_component(const vec3_t& v)
{
    return MAX(v.x, MAX(v.y, v.z));
}

inline vec3_t permute(const vec3_t& v, uint32_t x, uint32_t y, uint32_t z)
{
    return vec3_t(v[x], v[y], v[z]);
}

inline vec3_t project(const vec3_t& a, const vec3_t& b)
{
    return (b * dot(a, b) / magnitude_squared(b));
}

inline vec3_t reject(const vec3_t& a, const vec3_t& b)
{
    return (a - b * (dot(a, b) / magnitude_squared(b)));
}

// NOTE: Assumes the normal is already normalized.
inline vec3_t reflect(const vec3_t& v, const vec3_t& n)
{
    // vec3_t n_temp = normalize(n);
    return v - 2 * dot(v, n) * n;
}

inline vec3_t refract(const vec3_t& v, const vec3_t& n, float refractive_index_ratio) 
{
#if 0
    float cos_theta = fmin(ods::dot(-v, n), 1.0f);
    ods::vec3_t r_out_perp =  refractive_index_ratio * (v + cos_theta * n);
    ods::vec3_t r_out_parallel = -sqrt(fabs(1.0f - ods::dot(r_out_perp, r_out_perp))) * n;
    return r_out_perp + r_out_parallel;
#else
    vec3_t out;
    float dot_product = dot(n, v);
    float k = 1.f - refractive_index_ratio * refractive_index_ratio * (1.f - dot_product * dot_product);
    if(k < 0.0f)
        out = vec3_t(0.0f, 0.0f, 0.0f);
    else
        out = refractive_index_ratio * v - (refractive_index_ratio * dot_product + sqrtf(k)) * n;
    return out;
#endif
}

inline vec3_t random()
{
	return vec3_t(random_float(), random_float(), random_float());
}

inline vec3_t random(float min, float max)
{
	return vec3_t(random_float(min, max), random_float(min, max), random_float(min, max));
}

inline vec3_t random_in_unit_sphere()
{
	while(1)
	{
		vec3_t v = random(-1.0f, 1.0f);
		if(dot(v, v) >= 1) continue;
		return v;
	}
}

inline vec3_t random_in_unit_disk()
{
	while(1)
	{
		vec3_t v = vec3_t(random_float(-1.0f, 1.0f), random_float(-1.0f, 1.0f), 0.0f);
		if(dot(v, v) >= 1) continue;
		return v;
	}
}

inline vec3_t random_unit_vector()
{
	return normalize(random_in_unit_sphere());
}

// TODO: Implement vec3_t ==
#if 0
inline int32_t operator==(const vec3_t& a, const vec3_t& b)
{
    return ((a.x == b.x) && (a.y == b.y) && (a.z == b.z));
}
#endif

#endif
