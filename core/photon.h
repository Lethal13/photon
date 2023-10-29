#ifndef CORE_PHOTON_H
#define CORE_PHOTON_H

#include <stdint.h>

#if PHOTON_DEBUG
#define PHOTON_ASSERT(expression) if(!(expression)) { *(volatile int*)0 = 0; }
#else
#define PHOTON_ASSERT(expression)
#endif

#define ARRAY_SIZE(array) (sizeof(array)/sizeof((array)[0]))

#define LOCK_ADD(name) uint64_t name(uint64_t volatile *value, uint64_t add)
typedef LOCK_ADD(lock_add_function);

#include "../math/vec3.hpp"
#include "../math/point3.hpp"
#include "../math/ray.hpp"

#include "camera.hpp"
#include "color.hpp"

typedef struct
{
	void *data;
	uint32_t total_size;
	uint32_t width;
	uint32_t height;
	uint32_t header_size;
} ppm_image;

typedef struct
{
	void *pixels;
	uint32_t width;
	uint32_t height;
	uint32_t bytes_per_pixel;
	uint32_t pitch;
} framebuffer_t;

typedef struct
{
	uint32_t max_bounces;
	uint32_t samples_per_pixel;
	uint32_t total_cpu_cores;
	uint32_t cache_line_size;
	lock_add_function *lock_add;

	uint32_t tile_width;
	uint32_t tile_height;
	uint32_t total_tiles;

} settings_t;

enum material_type
{
	diffuse = 0,
};

typedef struct
{
	material_type material;
	color_t base_color;
	float roughness;
	float ior;
} material_t;

typedef struct
{
	point3_t center;
	float radius;
	uint32_t material_index;
} sphere_t;

typedef struct
{
	sphere_t *spheres;
	uint32_t total_spheres;
	material_t *materials;
	uint32_t total_materials;
} world_t;

typedef struct
{
	uint32_t x_min;
	uint32_t x_max;
	uint32_t y_min;
	uint32_t y_max;
} work_order;

typedef struct
{
	work_order *work_orders;
	volatile uint32_t total_work_orders;
	volatile uint64_t work_order_index;
	volatile uint64_t tiles_retired_counter;
} work_queue;

#define RAYTRACER(name) uint64_t name(framebuffer_t *buffer, world_t *world, camera_t *camera, settings_t *settings, work_queue *queue)
typedef RAYTRACER(raytrace_function);

#endif
