#ifndef CORE_PHOTON_H
#define CORE_PHOTON_H

#include <stdint.h>
#include "./../math/vec3d.hpp"
#include "./../math/vec3f.hpp"
#include "./../math/point3d.hpp"
#include "./../math/point3f.hpp"
#include "./../math/ray3d.hpp"
#include "./../math/ray3f.hpp"
#include "./../math/color3d.hpp"
#include "./../math/color3f.hpp"
#include "./camera.cpp"

typedef enum
{
    MAT_UNDEFINED = 0,
    MAT_DIFFUSE = 1,
    MAT_METAL = 2,
    MAT_DIELECTRIC = 3
} material_type;

typedef enum
{
    IMG_UNDEFINED = 0,
    IMG_PPM = 1
} image_format;

typedef struct
{
	image_format format;
	uint32_t total_size;
	uint32_t width;
	uint32_t height;
	uint32_t header_size;
	void *data;
} image_t;

typedef struct
{
	uint32_t width;
	uint32_t height;
	uint32_t bytes_per_pixel;
	void *data;
} framebuffer_t;

typedef struct
{
	point3f_t center;
	float radius;
	uint32_t material_index;
} sphere_t;

typedef struct
{
    point3f_t vertices[3];
    uint32_t material_index;
} triangle_t;

typedef struct
{
	material_type material;
	color3f_t base_color;
	float roughness;
	float refraction_index;
} material_t;

typedef struct
{
	sphere_t *spheres;
	uint32_t total_spheres;
    triangle_t *triangles;
    uint32_t total_triangles;
	material_t *materials;
	uint32_t total_materials;
} world_t;

typedef struct
{
    volatile uint64_t bounces_computed;
} info_t;

typedef struct
{
    world_t *world;
    framebuffer_t framebuffer;
    uint32_t x_min;
    uint32_t x_max;
    uint32_t y_min;
    uint32_t y_max;
} work_order_t;

typedef struct
{
    uint32_t work_orders_count;
    work_order_t *work_orders;

    uint32_t samples_per_pixel;
    uint64_t bounces_computed;

    volatile uint64_t work_order_index;
    volatile uint64_t tile_retired_count;

    camera_t camera;
} work_queue_t;

typedef uint64_t (*lock_add_fn)(uint64_t volatile*, uint64_t);

#define RAYTRACE_FUNCTION_SP(name) uint32_t name(framebuffer_t *framebuffer, world_t *world, camera_t *camera, info_t *info)
typedef RAYTRACE_FUNCTION_SP(raytrace_function_sp);

#define RAYTRACE_FUNCTION_QUEUE_SP(name) uint32_t name(work_queue_t *work_queue, lock_add_fn lock_add)
typedef RAYTRACE_FUNCTION_QUEUE_SP(raytrace_function_queue_sp);

#define RAYTRACE_FUNCTION_ORDER_SP(name) uint32_t name(work_order_t *work_order, camera_t *camera, info_t *info)
typedef RAYTRACE_FUNCTION_ORDER_SP(raytrace_function_order_sp);

typedef struct
{
    work_queue_t *queue;
    lock_add_fn lock_add;
    raytrace_function_queue_sp *raytracer;
} thread_params_t;

#endif
