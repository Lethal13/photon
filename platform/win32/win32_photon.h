#ifndef PLATFORM_WIN32_WIN32_PHOTON_H
#define PLATFORM_WIN32_WIN32_PHOTON_H

#include "./../../core/photon.h"
#include <stdio.h>
#include <Windows.h>
#include <intrin.h>

typedef struct
{
	HMODULE library;
	raytrace_function_sp *raytracer_iterative_sp;
	raytrace_function_sp *raytracer_recursive_sp;
    raytrace_function_queue_sp *raytracer_queue_sp;
} win32_photon;

static win32_photon win32_load_photon_dll(const char *dll);
static void win32_unload_photon_dll(win32_photon *photon);
static void PLATFORM_WRITE_FILE(const char *filename, uint32_t memory_size, void *memory);
static void convert_framebuffer_to_image(void *data, uint32_t width, uint32_t height, image_t *image);

inline uint64_t get_os_timer_freq(void);
inline uint64_t read_os_timer(void);
inline uint64_t read_cpu_timer(void);
inline uint64_t estimate_cpu_freq(void);

uint64_t lock_add(uint64_t volatile *value, uint64_t add);
DWORD WINAPI worker_thread(void *lpParameter);
void create_worker_thread(thread_params_t *params);

void setup_camera_sp(const camera_setup_t& setup, camera_t* camera);
void run_raytracer_test_sp(const char *test_name, raytrace_function_sp raytracer_func, framebuffer_t *framebuffer,
        world_t *world, camera_t *camera, const uint64_t cpu_freq);
void run_raytracer_queue_test_sp(const char *test_name, raytrace_function_queue_sp raytracer_func, framebuffer_t *framebuffer,
        world_t *world, camera_t *camera, const uint64_t cpu_freq, uint32_t cpu_cores);

#endif
