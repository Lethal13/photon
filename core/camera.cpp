#include "camera.h"
#include "./../math/ray3f.hpp"
#include "./../math/ray3d.hpp"

void initialize_camera32(float aspect_ratio, uint32_t image_width, uint32_t image_height, const point3f_t& look_from, const point3f_t& look_at, const vec3f_t& v_up, float defocus_angle, float focus_dist, uint32_t samples_per_pixel, uint32_t max_bounces, camera32_t *camera)
{
	camera->camera_center = look_from;
	camera->defocus_angle = defocus_angle;
	camera->focus_dist = focus_dist; // magnitude(look_from - look_at)

	camera->vfov = 20.0f;
	float theta = degrees_to_radians(camera->vfov);
	float h = tan(theta / 2);
	float viewport_height = 2.0f * h * camera->focus_dist;
	float viewport_width = viewport_height * (float(image_width)/image_height);

	// Calculate the u, v, w unit basis vectors for the camera coordinate frame
	vec3f_t w = normalize(look_from - look_at);
	vec3f_t u = normalize(cross(v_up, w));
	vec3f_t v = cross(w, u);

	// Calculate the vectors across the horizontal and down the vertical viewport edges
	vec3f_t viewport_u = viewport_width * u; // Vector across viewport horizontal edge
	vec3f_t viewport_v = viewport_height * -v; // Vector down viewport vertical edge

	// Calculate the horizontal and vertical delta vectors from pixel to pixel
	camera->pixel_delta_u = viewport_u / (float)image_width;
	camera->pixel_delta_v = viewport_v / (float)image_height;

	// Calculate the location of the upper left pixel
	camera->viewport_upper_left = camera->camera_center - (camera->focus_dist * w) - viewport_u / 2 - viewport_v / 2;
	camera->pixel00_loc = camera->viewport_upper_left + 0.5f * (camera->pixel_delta_u + camera->pixel_delta_v);

	// Calculate the camera defocus disk basis vectors
	float defocus_radius = camera->focus_dist * tan(degrees_to_radians(camera->defocus_angle / 2));
	camera->defocus_disk_u = u * defocus_radius;
	camera->defocus_disk_v = v * defocus_radius;

	camera->samples_per_pixel = samples_per_pixel;
	camera->max_bounces = max_bounces;
}

ray3f_t get_camera_ray(camera32_t *camera, int s, int t)
{
	// Construct a camera ray originating from the defocus disk and directed at a randomly
	// sampled point around the pixel location i, j.

	// Creates a vector to a random point in the [-.5,-.5]-[+.5,+.5] unit square.
	vec3f_t offset = vec3f_t(random_float() - 0.5f, random_float() - 0.5f, 0.0f);

	point3f_t pixel_center = camera->pixel00_loc + ((s + offset.x) * camera->pixel_delta_u) + ((t + offset.y) * camera->pixel_delta_v);
	vec3f_t tmp = random_vec3f_in_unit_disk();
	point3f_t ray_origin = (camera->defocus_angle <= 0) ? camera->camera_center : camera->camera_center + (tmp.x * camera->defocus_disk_u) + (tmp.y * camera->defocus_disk_v);
	vec3f_t ray_direction = pixel_center - ray_origin;
	ray3f_t ray(ray_origin, ray_direction);
	return ray;
}

void initialize_camera64(double aspect_ratio, uint32_t image_width, uint32_t image_height, const point3d_t& look_from, const point3d_t& look_at, const vec3d_t& v_up, double defocus_angle, double focus_dist, uint32_t samples_per_pixel, uint32_t max_bounces, camera64_t *camera)
{
	camera->defocus_angle = defocus_angle;
	camera->focus_dist = focus_dist;

	camera->camera_center = look_from;

	camera->vfov = 20.0;
	double theta = degrees_to_radians(camera->vfov);
	double h = tan(theta / 2);
	double viewport_height = 2.0 * h * camera->focus_dist;
	double viewport_width = viewport_height * (float(image_width)/image_height);

	// Calculate the u, v, w unit basis vectors for the camera coordinate frame
	vec3d_t w = normalize(look_from - look_at);
	vec3d_t u = normalize(cross(v_up, w));
	vec3d_t v = cross(w, u);

	// Calculate the vectors across the horizontal and down the vertical viewport edges
	vec3d_t viewport_u = viewport_width * u; // Vector across viewport horizontal edge
	vec3d_t viewport_v = viewport_height * -v; // Vector down viewport vertical edge

	// Calculate the horizontal and vertical delta vectors from pixel to pixel
	camera->pixel_delta_u = viewport_u / (double)image_width;
	camera->pixel_delta_v = viewport_v / (double)image_height;

	// Calculate the location of the upper left pixel
	camera->viewport_upper_left = camera->camera_center - (camera->focus_dist * w) - viewport_u / 2 - viewport_v / 2;
	camera->pixel00_loc = camera->viewport_upper_left + 0.5 * (camera->pixel_delta_u + camera->pixel_delta_v);

	// Calculate the camera defocus disk basis vectors
	double defocus_radius = camera->focus_dist * tan(degrees_to_radians(camera->defocus_angle / 2));
	camera->defocus_disk_u = u * defocus_radius;
	camera->defocus_disk_v = v * defocus_radius;

	camera->samples_per_pixel = samples_per_pixel;
	camera->max_bounces = max_bounces;
}

ray3d_t get_camera_ray(camera64_t *camera, int s, int t)
{
	// Construct a camera ray originating from the defocus disk and directed at a randomly
	// sampled point around the pixel location i, j.

	// Creates a vector to a random point in the [-.5,-.5]-[+.5,+.5] unit square.
	vec3d_t offset = vec3d_t(random_double() - 0.5, random_double() - 0.5, 0.0);

	point3d_t pixel_center = camera->pixel00_loc + ((s + offset.x) * camera->pixel_delta_u) + ((t + offset.y) * camera->pixel_delta_v);
	vec3d_t tmp = random_vec3d_in_unit_disk();
	point3d_t ray_origin = (camera->defocus_angle <= 0) ? camera->camera_center : camera->camera_center + (tmp.x * camera->defocus_disk_u) + (tmp.y * camera->defocus_disk_v);
	vec3d_t ray_direction = pixel_center - ray_origin;
	ray3d_t ray(ray_origin, ray_direction);
	return ray;
}

