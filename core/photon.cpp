#include "photon.h"

typedef struct
{
	point3_t intersection_point;
	vec3_t normal;
	float t;
	uint32_t hit_material_index;
	uint32_t object_found;
} intersection_record_t;

double linear_to_srgb(double value)
{
    if(value < 0.0)
    {
        value = 0.0f;
    }
    
    if(value > 1.0f)
    {
        value = 1.0f;
    }
    
    double result = value * 12.92f;
    if(value > 0.0031308f)
    {
        result = 1.055 * std::pow(value, 1.0/2.4) - 0.055;
    }
    
    return result;
}

inline double clamp(double x, double min, double max)
{
    if (x < min) return min;
    if (x > max) return max;
    return x;
}

// ray-sphere intersection.
float hit_sphere(const point3_t& center, float radius, const ray_t& r, float t_min, float t_max)
{
    vec3_t oc = r.origin - center;
    float a = dot(r.direction, r.direction);
    float b = 2.0f * dot(oc, r.direction);
    float c = dot(oc, oc) - radius*radius;
    float discriminant = b * b - 4.0f * a * c;

    if(discriminant < 0.0f) return -1.0f;
    float sqrtd = std::sqrt(discriminant);

    // Find the nearest root that lies in the acceptable range.
    float root_a = (-b + sqrtd) / (2.0f * a);
    float root_b = (-b - sqrtd) / (2.0f * a);

    float t = root_a;

    if((root_b > t_min) && (root_b < root_a))
    {
        t = root_b;
        // return t;
    }

    if((t > t_min) && (t < t_max))
    {
        return t;
    }

    return -1.0f;
}

// https://www.scratchapixel.com/lessons/3d-basic-rendering/ray-tracing-rendering-a-triangle/barycentric-coordinates.html
// Möller–Trumbore intersection algorithm
// ray-triangle intersection.
float hit_triangle(const point3_t vertices[3], const ray_t& r, float t_min, float t_max)
{
	vec3_t edge1, edge2, h, s, q;
	float a, f, u, v;

	edge1.x = vertices[1].x - vertices[0].x;
	edge1.y = vertices[1].y - vertices[0].y;
	edge1.z = vertices[1].z - vertices[0].z;

	edge2.x = vertices[2].x - vertices[0].x;
	edge2.y = vertices[2].y - vertices[0].y;
	edge2.z = vertices[2].z - vertices[0].z;

	h = cross(r.direction, edge2);
	a = dot(edge1, h);

	static constexpr float epsilon = 1e-5f;
	// NOTE: Check if Ray and Triangle are parallel.
	if(a > -epsilon && a < epsilon) return -1.0f;

	f = 1.0f / a;
	s.x = r.origin.x - vertices[0].x;
	s.y = r.origin.y - vertices[0].y;
	s.z = r.origin.z - vertices[0].z;
	u = f * dot(s, h);

	// NOTE: Check if the intersection is outside of the triangle.
	if(u < 0.0f || u > 1.0f) return -1.0f;

	q = cross(s, edge1);
	v = f * dot(r.direction, q);

	// NOTE: Check if the intersection is outside of the triangle.
	if(v < 0.0f || (u + v) > 1.0f) return -1.0f;

	float t = f * dot(edge2, q);
	if(t > epsilon) return t;
	return -1.0f;
}

color_t cast_ray(world_t *world, ray_t *ray, uint32_t max_bounces)
{
	color_t pixel_color = {};
	float t_min = 0.001f;
	float t_max = FLOAT_INFINITY;
	color_t attenuation = {1.0f};

	float closest_so_far = FLOAT_INFINITY;

	for(uint32_t bounce = 0; bounce < max_bounces; ++bounce)
	{
		intersection_record_t intersection = {0};

		for(uint32_t sphere_index = 0; sphere_index < world->total_spheres; ++sphere_index)
		{
			sphere_t *sphere = &world->spheres[sphere_index];
			float t = hit_sphere(sphere->center, sphere->radius, *ray, t_min, t_max);

			if(t <= 0.0f) continue;
		
			if(t < closest_so_far) 
			{
				closest_so_far = t;
				intersection.intersection_point = ray->at(t);
				intersection.t = t;
				intersection.hit_material_index = sphere->material_index;
				intersection.normal = (ray->at(closest_so_far) - sphere->center) / sphere->radius;
				intersection.object_found = 1;
			}
		}

		for(uint32_t triangle_index = 0; triangle_index < world->total_triangles; ++triangle_index)
		{
			triangle_t *triangle = &world->triangles[triangle_index];
			float t = hit_triangle(triangle->vertices, *ray, t_min, t_max);

			if(t <= 0.0f) continue;

			if(t < closest_so_far)
			{
				vec3_t normal;
				vec3_t edge1, edge2;
				edge1.x = triangle->vertices[1].x - triangle->vertices[0].x;
				edge1.y = triangle->vertices[1].y - triangle->vertices[0].y;
				edge1.z = triangle->vertices[1].z - triangle->vertices[0].z;

				edge2.x = triangle->vertices[2].x - triangle->vertices[0].x;
				edge2.y = triangle->vertices[2].y - triangle->vertices[0].y;
				edge2.z = triangle->vertices[2].z - triangle->vertices[0].z;

				normal.x = edge1.y * edge2.z - edge1.z * edge2.y;
				normal.y = edge1.z * edge2.x - edge1.x * edge2.z;
				normal.z = edge1.x * edge2.y - edge1.y * edge2.x;

				normal = normalize(normal);

				closest_so_far = t;
				intersection.intersection_point = ray->at(t);
				intersection.t = t;
				intersection.hit_material_index = triangle->material_index;
				intersection.normal = normal;
				intersection.object_found = 1;
			}
		}

		if(intersection.object_found)
		{
			material_t *material = &world->materials[intersection.hit_material_index];

			switch(material->material)
			{
				case diffuse:
				{
					vec3_t normal = intersection.normal;
					if(dot(ray->direction, normal) > 0.0f) 
					{
						normal = -intersection.normal;
					}

					attenuation *= material->base_color;
				
					point3_t target = intersection.intersection_point + normal + random_unit_vector();
					ray->origin = intersection.intersection_point;
					ray->direction = target - intersection.intersection_point;

					// Catch the case where random_unit_vector() is opposite to normal, so the direction is 0.
					if(ray->direction.near_zero()) ray->direction = normal;
				} break;
			};

			// pixel_color.r = 255;
		}
		else
		{
			// pixel_color.b = 255;
			vec3_t unit_direction = normalize(ray->direction);
			float t = 0.5f * (unit_direction.y + 1.0f);
			color_t c = (1.0f - t) * color_t(1.0f, 1.0f, 1.0f) + t * color_t(0.5f, 0.7f, 1.0f);
			pixel_color = attenuation * c;
// end_point:
			break;
		}
	}
	return pixel_color;
}

void render_tile(framebuffer_t *buffer, world_t *world, camera_t *camera, settings_t *settings, work_order *order)
{

	for(uint32_t j = order->y_min; j < order->y_max; ++j)
	{
		uint32_t offset = (order->x_min + j * buffer->width) * 4;
		uint8_t *pixels = (uint8_t*)buffer->pixels + offset;
		for(uint32_t i = order->x_min; i < order->x_max; ++i)
		{
			color_t final_color = {};
			double contribution = 1.0 / settings->samples_per_pixel;
			for(uint32_t sample = 0; sample < settings->samples_per_pixel; ++sample)
			{
				ray_t ray = camera->get_ray((float)i, (float)j);
				final_color += cast_ray(world, &ray, settings->max_bounces);
			}

			final_color *= contribution;

			final_color.r = linear_to_srgb(clamp(final_color.r, 0.0, 0.999));
			final_color.g = linear_to_srgb(clamp(final_color.g, 0.0, 0.999));
			final_color.b = linear_to_srgb(clamp(final_color.b, 0.0, 0.999));

			float alpha = 255.0f;
			*pixels++ = (uint8_t)(255 * final_color.b);
			*pixels++ = (uint8_t)(255 * final_color.g);
			*pixels++ = (uint8_t)(255 * final_color.r);
			*pixels++ = (uint8_t)alpha;
		}
	}
}

// raytrace(framebuffer_t *buffer, world_t *world, camera_t *camera, settings_t *settings, work_queue *queue)
RAYTRACER(raytrace)
{
	PHOTON_ASSERT(buffer != 0 && world != 0 && settings != 0 && queue != 0);

	uint64_t work_order_index = settings->lock_add(&queue->work_order_index, 1);
	if(work_order_index >= queue->total_work_orders) return 0;

	work_order *order = queue->work_orders + work_order_index;

	render_tile(buffer, world, camera, settings, order);

	settings->lock_add(&queue->tiles_retired_counter, 1);

	return 1;
}
