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
} sphere32_t;

typedef struct
{
	point3d_t center;
	double radius;
	uint32_t material_index;
} sphere64_t;

typedef struct
{
    point3f_t vertices[3];
    uint32_t material_index;
} triangle32_t;

typedef struct
{
    point3d_t vertices[3];
    uint32_t material_index;
} triangle64_t;

typedef struct
{
	material_type material;
	color3f_t base_color;
	float roughness;
	float refraction_index;
} material32_t;

typedef struct
{
	material_type material;
	color3d_t base_color;
	double roughness;
	double refraction_index;
} material64_t;

typedef struct
{
	sphere32_t *spheres;
	uint32_t total_spheres;
    triangle32_t *triangles;
    uint32_t total_triangles;
	material32_t *materials;
	uint32_t total_materials;
} world32_t;

typedef struct
{
	sphere64_t *spheres;
	uint32_t total_spheres;
    triangle64_t *triangles;
    uint32_t total_triangles;
	material64_t *materials;
	uint32_t total_materials;
} world64_t;

typedef struct
{
    volatile uint64_t bounces_computed;
} info_t;

typedef struct
{
    world32_t *world;
    framebuffer_t framebuffer;
    uint32_t x_min;
    uint32_t x_max;
    uint32_t y_min;
    uint32_t y_max;
} work_order32_t;

typedef struct
{
    uint32_t work_orders_count;
    work_order32_t *work_orders;

    uint32_t samples_per_pixel;
    uint64_t bounces_computed;

    volatile uint64_t work_order_index;
    volatile uint64_t tile_retired_count;

    camera32_t camera;
} work_queue32_t;

typedef struct
{
    world64_t *world;
    framebuffer_t framebuffer;
    uint32_t x_min;
    uint32_t x_max;
    uint32_t y_min;
    uint32_t y_max;
} work_order64_t;

typedef struct
{
    uint32_t work_orders_count;
    work_order64_t *work_orders;

    uint32_t samples_per_pixel;
    uint64_t bounces_computed;

    volatile uint64_t work_order_index;
    volatile uint64_t tile_retired_count;

    camera64_t camera;
} work_queue64_t;

typedef uint64_t (*lock_add_fn)(uint64_t volatile*, uint64_t);

#define RAYTRACE_FUNCTION_SP(name) uint32_t name(framebuffer_t *framebuffer, world32_t *world, camera32_t *camera, info_t *info)
typedef RAYTRACE_FUNCTION_SP(raytrace_function_sp);

#define RAYTRACE_FUNCTION_DP(name) uint32_t name(framebuffer_t *framebuffer, world64_t *world, camera64_t *camera, info_t *info)
typedef RAYTRACE_FUNCTION_DP(raytrace_function_dp);

#define RAYTRACE_FUNCTION_QUEUE_SP(name) uint32_t name(work_queue32_t *work_queue, lock_add_fn lock_add)
typedef RAYTRACE_FUNCTION_QUEUE_SP(raytrace_function_queue_sp);

#define RAYTRACE_FUNCTION_QUEUE_DP(name) uint32_t name(work_queue64_t *work_queue, lock_add_fn lock_add)
typedef RAYTRACE_FUNCTION_QUEUE_DP(raytrace_function_queue_dp);

#define RAYTRACE_FUNCTION_ORDER_SP(name) uint32_t name(work_order32_t *work_order, camera32_t *camera, info_t *info)
typedef RAYTRACE_FUNCTION_ORDER_SP(raytrace_function_order_sp);

#define RAYTRACE_FUNCTION_ORDER_DP(name) uint32_t name(work_order64_t *work_order, camera64_t *camera, info_t *info)
typedef RAYTRACE_FUNCTION_ORDER_DP(raytrace_function_order_dp);

typedef struct
{
    work_queue32_t *queue;
    lock_add_fn lock_add;
    raytrace_function_queue_sp *raytracer;
} thread_params32_t;

typedef struct
{
    work_queue64_t *queue;
    lock_add_fn lock_add;
    raytrace_function_queue_dp *raytracer;
} thread_params64_t;

#endif
