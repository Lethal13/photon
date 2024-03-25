#ifndef PLATFORMS_WIN32_PHOTON
#define PLATFORMS_WIN32_PHOTON

// TODO list:

#include "../../core/photon.h"
#include "../../vendors/obj.h"
#include "../../core/test_worlds.hpp"
#include <Windows.h>
#include <intrin.h>
#include <sstream>
#include <cstring>

struct win32_raytrace_code
{
	HMODULE library;
	raytrace_function *raytracer;
	add_filter_function *add_filter;
};

struct win32_wrapper
{
	framebuffer_t *framebuffer;
	world_t *world;
	camera_t *camera;
	settings_t *settings;
	work_queue *queue;
	raytrace_function *raytracer;
	add_filter_function *add_filter;
};

static uint64_t get_os_timer_freq(void);
static uint64_t read_os_timer(void);
inline uint64_t read_cpu_timer(void);
static uint64_t estimate_cpu_time_freq(void); // returns an estimation of cpu clocks per second.

win32_raytrace_code load_raytracer_library(char *dll);
uint32_t get_total_cpu_cores(void);
uint32_t get_cache_line_size();
LOCK_ADD(win32_lock_add);
DWORD WINAPI worker_thread(void *lpParameter);
void create_worker_thread(void *parameter);
void PLATFORM_WRITE_FILE(const char* file_name, uint32_t memory_size, void* memory);
void convert_framebuffer_to_ppm_image(void *data, uint32_t width, uint32_t height, ppm_image *image);

#endif
