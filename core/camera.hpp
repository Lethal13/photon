#ifndef PHOTON_CORE_CAMERA_HPP
#define PHOTON_CORE_CAMERA_HPP

#include "../math/vec3.hpp"
#include "../math/point3.hpp"
#include "../math/ray.hpp"

struct camera_t
{

	camera_t() = default;

	camera_t(point3_t look_from, point3_t look_at, vec3_t vup, float vfov, uint32_t image_width, uint32_t image_height)
	{
		position = look_from;
		direction = normalize(look_at - look_from);

		camera_z = normalize(direction);
		camera_x = normalize(cross(camera_z, vup));
		camera_y = cross(camera_x, camera_z);

		// Focus distance = the distance between the projection point and the plane where everything is in perfect focus.
		// Focal length = is the distance between the projection point and the image plane.
		// focal_length = magnitude(look_at - look_from);
		focal_length = 1.0f;

		// TODO: Add fov support.
		// float theta = degrees_to_radians(vfov);
		// float h = tan(theta/2);

		// film_distance is the distance between the film_plane and the pinhole point or lense.
		film_distance = 1.0f;

		film_width = 2.0f;
		film_height = 2.0f;

		if(image_width >= image_height)
		{
			film_aspect_ratio = (float)image_width / (float)image_height;
			film_width *= film_aspect_ratio;
		}
		else
		{
			film_aspect_ratio = (float)image_height / (float)image_width;
			film_height *= film_aspect_ratio;
		}

		// Calculate the vectors across the horizontal and down the vertical viewport edges.
		film_u = vec3_t(film_width, 0.0f, 0.0f);
		film_v = vec3_t(0.0f, -film_height, 0.0f);

		// Calculate the horizontal and vertical delta vectors from pixel to pixel.
		pixel_delta_u = film_u / (float)image_width;
		pixel_delta_v = film_v / (float)image_height;

		// Calculate the location of the upper left pixel.
		// film_upper_left = position + vec3_t(0.0f, 0.0f, focal_length) - film_u/2.0f - film_v/2.0f;
		film_upper_left = position + vec3_t(0.0f, 0.0f, film_distance) - film_u/2.0f - film_v/2.0f;
		// NOTE: We add 0.5f, in order to be at the center of the pixel, and not in the upper left.
		pixel00_loc = film_upper_left + 0.5f * (pixel_delta_u + pixel_delta_v);
	}

	ray_t get_ray(float x, float y) const
	{
		point3_t pixel_center = pixel00_loc + x * pixel_delta_u + y * pixel_delta_v;
		vec3_t ray_direction = position - pixel_center;
		ray_t ray{ position, ray_direction};

		return ray;
	}

	point3_t position;
	vec3_t direction;
	vec3_t camera_x, camera_y, camera_z;
	float film_aspect_ratio;
	float film_width, film_height;
	vec3_t film_u, film_v;
	vec3_t pixel_delta_u, pixel_delta_v;
	float focal_length, film_distance;

	point3_t film_upper_left, pixel00_loc;
};

#endif
