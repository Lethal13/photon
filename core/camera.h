#ifndef CORE_CAMERA_H
#define CORE_CAMERA_H

#include "./../math/vec3f.hpp"
#include "./../math/vec3d.hpp"
#include "./../math/point3f.hpp"
#include "./../math/point3d.hpp"

typedef struct
{
    float aspect_ratio;
    uint32_t width;
    uint32_t height;
    float defocus_angle;
    float focus_dist;
    uint32_t samples_per_pixel;
    uint32_t max_bounces;
} camera_setup_t;

typedef struct
{
	point3f_t camera_center;
	vec3f_t pixel_delta_u; // Offset to pixel to the right
	vec3f_t pixel_delta_v; // Offset to pixel below
	point3f_t viewport_upper_left;
	point3f_t pixel00_loc; // Location of pixel 0, 0
	uint32_t samples_per_pixel; // Count of random samples for each pixel
	uint32_t max_bounces; // Maximum number of ray bounces into scene
	float vfov; // Vertical view angle (field of view)
	float defocus_angle; // Variation angle of rays through each pixel
	float focus_dist; // Distance from camera look_from point to plane of perfect focus.
	vec3f_t defocus_disk_u; // Defocus disk horizontal radius
	vec3f_t defocus_disk_v; // Defocus disk vertical radius
} camera32_t;

typedef struct
{
	point3d_t camera_center;
	vec3d_t pixel_delta_u;
	vec3d_t pixel_delta_v;
	point3d_t viewport_upper_left;
	point3d_t pixel00_loc;
	uint32_t samples_per_pixel; // Count of random samples for each pixel
	uint32_t max_bounces; // Maximum number of ray bounces into scene
	double vfov; // Vertical view angle (field of view)
	double defocus_angle; // Variation angle of rays through each pixel
	double focus_dist; // Distance from camera look_from point to plane of perfect focus.
	vec3d_t defocus_disk_u; // Defocus disk horizontal radius
	vec3d_t defocus_disk_v; // Defocus disk vertical radius
} camera64_t;

#endif
