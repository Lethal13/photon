#ifndef PLATFORMS_WIN32_PHOTON
#define PLATFORMS_WIN32_PHOTON

#include "../../core/photon.h"
#include <Windows.h>
#include <sstream>
#include <cstring>

struct win32_raytrace_code
{
	HMODULE library;
	raytrace_function *raytracer;
};

struct win32_wrapper
{
	framebuffer_t *buffer;
	world_t *world;
	camera_t *camera;
	settings_t *settings;
	work_queue *queue;
	raytrace_function *raytracer;
};


win32_raytrace_code load_raytracer_library(char *dll);
uint32_t get_total_cpu_cores(void);
uint32_t get_cache_line_size();
LOCK_ADD(win32_lock_add);
DWORD WINAPI worker_thread(void *lpParameter);
void create_worker_thread(void *parameter);
void PLATFORM_WRITE_FILE(const char* file_name, uint32_t memory_size, void* memory);
void convert_framebuffer_to_ppm_image(void *data, uint32_t width, uint32_t height, ppm_image *image);

#endif