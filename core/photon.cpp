#include "photon.h"

typedef struct
{
	point3f_t intersection_point;
	vec3f_t normal;
	float t;
	uint32_t hit_material_index;
	bool object_found;
} intersection_record_t;

static float reflectance(float cosine, float refraction_index)
{
	// Use Schlick's approximation for reflectance.
	float r0 = (1.0f - refraction_index) / (1.0f + refraction_index);
	r0 = r0*r0;
	return r0 + (1.0f - r0) * (float)pow((1.0f - cosine), 5);
}

float hit_sphere(const point3f_t& center, float radius, const ray3f_t& ray, float t_min, float t_max)
{
    vec3f_t oc = center - ray.origin;  // This is correct as is
    float a = magnitude_squared(ray.direction);
    float half_b = dot(ray.direction, oc);  // Changed to half_b
    float c = dot(oc, oc) - radius * radius;
    float discriminant = half_b * half_b - a * c;
    
    if(discriminant < 0.0f)
        return -1.0f;
    
    float sqrtd = sqrtf(discriminant);
    float root = (half_b - sqrtd) / a;
    if(root < t_min || root > t_max)
    {
        root = (half_b + sqrtd) / a;
        if(root < t_min || root > t_max)
            return -1.0f;
    }
    return root;
}

// Möller–Trumbore intersection algorithm
// ray-triangle intersection.
float hit_triangle(const point3f_t vertices[3], const ray3f_t& r, float t_min, float t_max)
{
	vec3f_t edge1, edge2, h, s, q;
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

color3f_t compute_color_it32(world_t *world, ray3f_t *r, uint32_t max_bounces, info_t *info)
{
	color3f_t accumulated_color(1.0f, 1.0f, 1.0f);
	float t_min = 0.001f;
	float t_max = FLOAT_INFINITY;

	for(uint32_t bounce = 0; bounce < max_bounces; ++bounce)
	{
        ++info->bounces_computed;
		intersection_record_t intersection = {0};
		float closest_so_far = FLOAT_INFINITY;

		for(uint32_t sphere_index = 0; sphere_index < world->total_spheres; ++sphere_index)
		{
			sphere_t *sphere = &world->spheres[sphere_index];
			float t = hit_sphere(sphere->center, sphere->radius, *r, t_min, t_max);

			if(t <= 0.0f) continue;

			if(t < closest_so_far)
			{
				closest_so_far = t;
				intersection.intersection_point = r->at(t);
				intersection.t = t;
				intersection.hit_material_index = sphere->material_index;
				intersection.normal = (r->at(t) - sphere->center) / sphere->radius;
				intersection.object_found = true;
			}
		}

		for(uint32_t triangle_index = 0; triangle_index < world->total_triangles; ++triangle_index)
		{
			triangle_t *triangle = &world->triangles[triangle_index];

			float t = hit_triangle(triangle->vertices, *r, t_min, t_max);

			if(t <= 0.0f) continue;

			if(t < closest_so_far)
			{
                vec3f_t normal;
				vec3f_t edge1;
				vec3f_t edge2;
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
				intersection.intersection_point = r->at(t);
				intersection.t = t;
				intersection.hit_material_index = triangle->material_index;
				intersection.normal = normal;
				intersection.object_found = true;
			}
		}

		if(intersection.object_found)
		{
			material_t *material = &world->materials[intersection.hit_material_index];
			switch(material->material)
			{
				case MAT_DIFFUSE:
				{
					accumulated_color = accumulated_color * material->base_color;
					bool front_face = dot(r->direction, intersection.normal) < 0;
					intersection.normal = front_face ? intersection.normal : -intersection.normal;

					vec3f_t direction = intersection.normal + random_unit_vec3f();

					if(direction.near_zero())
						direction = intersection.normal;

					*r = ray3f_t(intersection.intersection_point, direction);
				} break;
				case MAT_METAL:
				{
					material->roughness = material->roughness < 1.0f ? material->roughness : 1.0f;
					accumulated_color = accumulated_color * material->base_color;
					bool front_face = dot(r->direction, intersection.normal) < 0;
					intersection.normal = front_face ? intersection.normal : -intersection.normal;
					vec3f_t direction = normalize(reflect(r->direction, intersection.normal)) + (material->roughness * random_unit_vec3f());
					if(dot(direction, intersection.normal) > 0.0f)
						*r = ray3f_t(intersection.intersection_point, direction);
					else 
						goto end;
				} break;
				case MAT_DIELECTRIC:
				{
					accumulated_color = accumulated_color * material->base_color;
					bool front_face = dot(r->direction, intersection.normal) < 0;
					intersection.normal = front_face ? intersection.normal : -intersection.normal;
					float ri = front_face ? (1.0f / material->refraction_index) : material->refraction_index;

					vec3f_t normalized_direction = normalize(r->direction);
					float cos_theta = fmin(dot(-normalized_direction, intersection.normal), 1.0f);
					float sin_theta = sqrtf(1.0f - cos_theta * cos_theta);
					
					bool cannot_refract = ri * sin_theta > 1.0f;
					vec3f_t direction;

					if(cannot_refract || reflectance(cos_theta, ri) > random_float())
						direction = reflect(normalized_direction, intersection.normal);
					else
						direction = refract(normalized_direction, intersection.normal, ri);

					*r = ray3f_t(intersection.intersection_point, direction);
				} break;
				default:
				break;
			}
		}
		else
		{
			vec3f_t unit_direction = normalize(r->direction);
            float a = 0.5f * (unit_direction.y + 1.0f);
			color3f_t background_color = (1.0f - a) * color3f_t(1.0f, 1.0f, 1.0f) + a * color3f_t(0.5f, 0.7f, 1.0f);
            // color3f_t background_color = color3f_t(0.145f, 0.169f, 0.216f);
			return accumulated_color * background_color;
		}
	}
end:
	return color3f_t(0.0f, 0.0f, 0.0f);
}

color3f_t compute_color_rec32(world_t *world, const ray3f_t& r, uint32_t max_bounces, info_t *info)
{
	if(max_bounces <= 0)
		return color3f_t(0.0f, 0.0f, 0.0f);

    ++info->bounces_computed;

	float t_min = 0.001f;
	float t_max = FLOAT_INFINITY;
	intersection_record_t intersection = {0};
	float closest_so_far = FLOAT_INFINITY;
	for(uint32_t sphere_index = 0; sphere_index < world->total_spheres; ++sphere_index)
	{
		sphere_t *sphere = &world->spheres[sphere_index];
		float t = hit_sphere(sphere->center, sphere->radius, r, t_min, t_max);

		if(t <= 0.0f) continue;

		if(t < closest_so_far)
		{
			closest_so_far = t;
			intersection.intersection_point = r.at(t);
			intersection.t = t;
			intersection.hit_material_index = sphere->material_index;
			intersection.normal = (r.at(t) - sphere->center) / sphere->radius;
			intersection.object_found = true;
		}
	}

    for(uint32_t triangle_index = 0; triangle_index < world->total_triangles; ++triangle_index)
    {
        triangle_t *triangle = &world->triangles[triangle_index];

        float t = hit_triangle(triangle->vertices, r, t_min, t_max);

        if(t <= 0.0f) continue;

        if(t < closest_so_far)
        {
            vec3f_t normal;
            vec3f_t edge1;
            vec3f_t edge2;
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
            intersection.intersection_point = r.at(t);
            intersection.t = t;
            intersection.hit_material_index = triangle->material_index;
            intersection.normal = normal;
            intersection.object_found = true;
        }
    }

	if(intersection.object_found)
	{
		material_t *material = &world->materials[intersection.hit_material_index];
		switch(material->material)
		{
			case MAT_DIFFUSE:
			{
				bool front_face = dot(r.direction, intersection.normal) < 0;
				intersection.normal = front_face ? intersection.normal : -intersection.normal;
				vec3f_t direction = intersection.normal + random_unit_vec3f();
				if(direction.near_zero())
					direction = intersection.normal;

				return material->base_color * compute_color_rec32(world, ray3f_t(intersection.intersection_point, direction), max_bounces - 1, info);
			} break;
			case MAT_METAL:
			{
				material->roughness = material->roughness < 1.0f ? material->roughness : 1.0f;
				bool front_face = dot(r.direction, intersection.normal) < 0;
				intersection.normal = front_face ? intersection.normal : -intersection.normal;
				vec3f_t direction = normalize(reflect(r.direction, intersection.normal)) + (material->roughness * random_unit_vec3f());
				if(dot(direction, intersection.normal) > 0.0f)
					return material->base_color * compute_color_rec32(world, ray3f_t(intersection.intersection_point, direction), max_bounces - 1, info);
				else
					return color3f_t(0.0f, 0.0f, 0.0f);
			} break;
			case MAT_DIELECTRIC:
			{
				bool front_face = dot(r.direction, intersection.normal) < 0;
				intersection.normal = front_face ? intersection.normal : -intersection.normal;

				float ri = front_face ? (1.0f / material->refraction_index) : material->refraction_index;
				vec3f_t normalized_direction = normalize(r.direction);

				float cos_theta = fmin(dot(-normalized_direction, intersection.normal), 1.0f);
				float sin_theta = sqrtf(1.0f - cos_theta * cos_theta);
				
				bool cannot_refract = ri * sin_theta > 1.0f;
				vec3f_t direction;

				if(cannot_refract || reflectance(cos_theta, ri) > random_float())
					direction = reflect(normalized_direction, intersection.normal);
				else
					direction = refract(normalized_direction, intersection.normal, ri);

				return material->base_color * compute_color_rec32(world, ray3f_t(intersection.intersection_point, direction), max_bounces - 1, info);
			} break;
			default:
			break;
		}
	}

    vec3f_t unit_direction = normalize(r.direction);
    float a = 0.5f * (unit_direction.y + 1.0f);
    return (1.0f - a) * color3f_t(1.0f, 1.0f, 1.0f) + a * color3f_t(0.5f, 0.7f, 1.0f);
}

RAYTRACE_FUNCTION_SP(raytracer_iterative_sp)
{
	uint32_t *pixels = (uint32_t*)framebuffer->data;
	for(uint32_t j = 0; j < framebuffer->height; ++j)
	{
		for(uint32_t i = 0; i < framebuffer->width; ++i)
		{
			float contribution = 1.0f / camera->samples_per_pixel;
			color3f_t pixel_color(0.0f, 0.0f, 0.0f);
			for(uint32_t sample = 0; sample < camera->samples_per_pixel; ++sample)
			{
				ray3f_t camera_ray = get_camera_ray(camera, i, j);
				pixel_color += compute_color_it32(world, &camera_ray, camera->max_bounces, info);
			}

			pixel_color *= contribution;
			pixel_color.r = linear_to_srgb(pixel_color.r);
			pixel_color.g = linear_to_srgb(pixel_color.g);
			pixel_color.b = linear_to_srgb(pixel_color.b);


			float a = 1.0f;
			// Translate the [0,1] component values to the byte range [0,255].
			int32_t ia = int(256 * clamp(a, 0.0f, 0.999f));
			int32_t ir = int(256 * clamp(pixel_color.r, 0.0f, 0.999f));
			int32_t ig = int(256 * clamp(pixel_color.g, 0.0f, 0.999f));
			int32_t ib = int(256 * clamp(pixel_color.b, 0.0f, 0.999f));

			// Pack components into a single 32-bit ARGB value.
			uint32_t final_value = 0;
			final_value = ia << 24 | ir << 16 | ig << 8 | ib << 0;
			*pixels++ = final_value;
		}
	}
	return 0;
}

RAYTRACE_FUNCTION_SP(raytracer_recursive_sp)
{
	uint32_t *pixels = (uint32_t*)framebuffer->data;
	for(uint32_t j = 0; j < framebuffer->height; ++j)
	{
		for(uint32_t i = 0; i < framebuffer->width; ++i)
		{
			float contribution = 1.0f / camera->samples_per_pixel;
			color3f_t pixel_color(0.0f, 0.0f, 0.0f);
			for(uint32_t sample = 0; sample < camera->samples_per_pixel; ++sample)
			{
				ray3f_t camera_ray = get_camera_ray(camera, i, j);
				pixel_color += compute_color_rec32(world, camera_ray, camera->max_bounces, info);
			}

			pixel_color *= contribution;
			pixel_color.r = linear_to_srgb(pixel_color.r);
			pixel_color.g = linear_to_srgb(pixel_color.g);
			pixel_color.b = linear_to_srgb(pixel_color.b);

			float a = 1.0f;
			// Translate the [0,1] component values to the byte range [0,255].
			int32_t ia = int(256 * clamp(a, 0.0f, 0.999f));
			int32_t ir = int(256 * clamp(pixel_color.r, 0.0f, 0.999f));
			int32_t ig = int(256 * clamp(pixel_color.g, 0.0f, 0.999f));
			int32_t ib = int(256 * clamp(pixel_color.b, 0.0f, 0.999f));


			// Pack components into a single 32-bit ARGB value.
			uint32_t final_value = 0;
			final_value = ia << 24 | ir << 16 | ig << 8 | ib << 0;
			*pixels++ = final_value;
		}
	}
	return 0;
}

RAYTRACE_FUNCTION_ORDER_SP(raytracer_iterative_queue_sp)
{
    world_t *world = work_order->world;
    framebuffer_t *framebuffer = &work_order->framebuffer;

    uint32_t row_width = framebuffer->width;
    uint32_t starting_offset = work_order->x_min + (work_order->y_min * row_width);
    uint32_t *pixels = (uint32_t*)framebuffer->data + starting_offset;

	for(uint32_t j = work_order->y_min; j < work_order->y_max; ++j)
	{
		for(uint32_t i = work_order->x_min; i < work_order->x_max; ++i)
		{
			float contribution = 1.0f / camera->samples_per_pixel;
			color3f_t pixel_color(0.0f, 0.0f, 0.0f);
			for(uint32_t sample = 0; sample < camera->samples_per_pixel; ++sample)
			{
				ray3f_t camera_ray = get_camera_ray(camera, i, j);
				pixel_color += compute_color_it32(world, &camera_ray, camera->max_bounces, info);
			}

			pixel_color *= contribution;
			pixel_color.r = linear_to_srgb(pixel_color.r);
			pixel_color.g = linear_to_srgb(pixel_color.g);
			pixel_color.b = linear_to_srgb(pixel_color.b);


			float a = 1.0f;
			// Translate the [0,1] component values to the byte range [0,255].
			int32_t ia = int(256 * clamp(a, 0.0f, 0.999f));
			int32_t ir = int(256 * clamp(pixel_color.r, 0.0f, 0.999f));
			int32_t ig = int(256 * clamp(pixel_color.g, 0.0f, 0.999f));
			int32_t ib = int(256 * clamp(pixel_color.b, 0.0f, 0.999f));

			// Pack components into a single 32-bit ARGB value.
			uint32_t final_value = 0;
			final_value = ia << 24 | ir << 16 | ig << 8 | ib << 0;
			*pixels++ = final_value;
		}
        pixels += (row_width - (work_order->x_max - work_order->x_min));
	}
	return 0;
}

RAYTRACE_FUNCTION_QUEUE_SP(raytracer_queue_sp)
{
	uint64_t work_order_index = lock_add(&work_queue->work_order_index, 1);

	if(work_order_index >= work_queue->work_orders_count)
        return 0;

	work_order_t *work_order = work_queue->work_orders + work_order_index;

    info_t info = {0};
    raytracer_iterative_queue_sp(work_order, &work_queue->camera, &info);

	lock_add(&work_queue->tile_retired_count, 1);
	lock_add(&work_queue->bounces_computed, info.bounces_computed);

	return 1;
}
