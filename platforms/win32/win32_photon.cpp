#include "win32_photon.h"
#include <stdio.h>

static uint64_t get_os_timer_freq(void)
{
	LARGE_INTEGER value;
	QueryPerformanceFrequency(&value); // ticks per second.
	return value.QuadPart;
}

static uint64_t read_os_timer(void)
{
	LARGE_INTEGER value;
	QueryPerformanceCounter(&value); // get the current value of ticks.
	return value.QuadPart;
}

inline uint64_t read_cpu_timer(void)
{
	return __rdtsc();
}

static uint64_t estimate_cpu_time_freq(void)
{
	uint64_t milliseconds_to_wait = 100;
	uint64_t os_freq = get_os_timer_freq();

	uint64_t cpu_start = read_cpu_timer();
	uint64_t os_start = read_os_timer();
	uint64_t os_end = 0;
	uint64_t os_elapsed = 0;
	uint64_t os_wait_time = os_freq * milliseconds_to_wait / 1000;
	while(os_elapsed < os_wait_time)
	{
		os_end = read_os_timer();
		os_elapsed = os_end - os_start;
	}
	
	uint64_t cpu_end = read_cpu_timer();
	uint64_t cpu_elapsed = cpu_end - cpu_start;
	
	uint64_t cpu_freq = 0;
	if(os_elapsed)
	{
		cpu_freq = os_freq * cpu_elapsed / os_elapsed;
	}
	
	return cpu_freq;
}

// TODO: Collapse all memory allocations, into one big memory allocation at startup.
// TODO: Add the case that library was not loaded properly.
win32_raytrace_code load_raytracer_library(char *dll)
{
	win32_raytrace_code result = {0};
	result.library =  LoadLibraryA(dll);
	result.raytracer = (raytrace_function*)GetProcAddress(result.library, "raytrace");

	return result;
}

uint32_t get_total_cpu_cores(void)
{
	SYSTEM_INFO system_info;
	GetSystemInfo(&system_info);
	uint32_t total_cpu_cores = system_info.dwNumberOfProcessors;

	return total_cpu_cores;
}

// https://learn.microsoft.com/en-us/windows/win32/api/sysinfoapi/nf-sysinfoapi-getlogicalprocessorinformation
// https://stackoverflow.com/questions/150294/how-to-programmatically-get-the-cpu-cache-line-size-in-c
uint32_t get_cache_line_size()
{
	uint32_t result = 64;
	using CpuInfo = SYSTEM_LOGICAL_PROCESSOR_INFORMATION;
	DWORD length;
	GetLogicalProcessorInformation(0, &length);

	CpuInfo *processor_information = (CpuInfo*)malloc(length);
	GetLogicalProcessorInformation(processor_information, &length);

	DWORD count = length / sizeof(CpuInfo);
	for(DWORD i = 0; i < count; ++i)
	{
		// This will be true for multiple returned caches, we need just one
		if(processor_information[i].Relationship == RelationCache)
		{
			result = processor_information[i].Cache.LineSize;
			break;
		}
	}

	free(processor_information);
	return result;
}

LOCK_ADD(win32_lock_add)
{
	uint64_t result = InterlockedExchangeAdd64((volatile int64_t*)value, add);
	return result;
}

DWORD WINAPI worker_thread(void *lpParameter)
{
	win32_wrapper *wrapper = (win32_wrapper*)lpParameter;
	while(wrapper->raytracer(wrapper->framebuffer, wrapper->world, wrapper->camera, wrapper->settings, wrapper->queue))
	{
	}

	return 0;
}

void create_worker_thread(void *parameter)
{
	DWORD thread_id;
	HANDLE thread_handle = CreateThread(0, 0, worker_thread, parameter, 0, &thread_id);
	CloseHandle(thread_handle);
}

void PLATFORM_WRITE_FILE(const char* file_name, uint32_t memory_size, void* memory)
{
    HANDLE file_handle = CreateFile(file_name, GENERIC_WRITE, 0, 0, OPEN_ALWAYS, 0, 0);

    if (file_handle != INVALID_HANDLE_VALUE)
    {
        DWORD bytes_written;

        //NOTE: File read successfully.
        if (WriteFile(file_handle, memory, memory_size, &bytes_written, 0))
        {
            // result = (bytes_written == memory_size);
        }
        else
        {
            // TODO: Logging.
        }

        CloseHandle(file_handle);
    }
    else
    {
        //TODO: Logging.
    }
}

void convert_framebuffer_to_ppm_image(void *data, uint32_t width, uint32_t height, ppm_image *image)
{
	// TODO: Find a way to convert framebuffer to image without using std library.
	std::ostringstream oss;
    oss << "P6\n" << width << ' ' << height << '\n' << "255\n";
	std::string ppm_header = oss.str();

	uint32_t ppm_header_size = (uint32_t)strlen(ppm_header.c_str());
	uint32_t total_size = (uint32_t)(ppm_header_size + width * height * sizeof(uint8_t) * 3);

	void *image_data = malloc(total_size);

	memset(image_data, 0, total_size);
	memcpy(image_data, ppm_header.c_str(), ppm_header_size);

	uint8_t *output_data = ((uint8_t*)image_data) + ppm_header_size;
	uint32_t current_pixel = 0;
	for(uint32_t i = 0; i < width * height; ++i)
	{
		uint32_t pixel = ((uint32_t*)data)[i];
		uint8_t red = (uint8_t)((0x00FF0000 & pixel) >> 16);
		uint8_t green = (uint8_t)((0x0000FF00 & pixel) >> 8);
		uint8_t blue = (uint8_t)((0x000000FF & pixel) >> 0);

		output_data[current_pixel] = red;
		output_data[current_pixel + 1] = green;
		output_data[current_pixel + 2] = blue;
		current_pixel += 3;
	}

	image->data = image_data;
	image->width = width;
	image->height = height;
	image->header_size = ppm_header_size;
	image->total_size = total_size;
}

void create_work_queue(work_queue *queue, settings_t *settings, world_t *world, framebuffer_t *framebuffer)
{
	uint32_t tile_width = (framebuffer->width + settings->total_cpu_cores - 1) / settings->total_cpu_cores;
	uint32_t tile_height = tile_width;

	uint32_t tile_count_x = (framebuffer->width + tile_width - 1) / tile_width;
	uint32_t tile_count_y = (framebuffer->height + tile_height -1) / tile_height;
	uint32_t total_tiles = tile_count_x * tile_count_y;

	settings->tile_width = tile_width;
	settings->tile_height = tile_height;
	settings->total_tiles = total_tiles;

	queue->total_work_orders = 0;
	queue->work_orders = (work_order*)malloc(sizeof(work_order) * total_tiles);

	for(uint32_t tile_y = 0; tile_y < tile_count_y; ++tile_y)
	{
		uint32_t y_min = tile_y * tile_height;
		uint32_t y_max = y_min + tile_height;

		if(y_max > framebuffer->height)
		{
			y_max = framebuffer->height;
		}

		for(uint32_t tile_x = 0; tile_x < tile_count_x; ++tile_x)
		{
			uint32_t x_min = tile_x * tile_width;
			uint32_t x_max = x_min + tile_width;

			if(x_max > framebuffer->width)
			{
				x_max = framebuffer->width;
			}

			work_order *order = queue->work_orders + queue->total_work_orders++;
			PHOTON_ASSERT(queue->total_work_orders <= total_tiles);

			order->x_min = x_min;
			order->x_max = x_max;
			order->y_min = y_min;
			order->y_max = y_max;
		}
	}

	PHOTON_ASSERT(queue->total_work_orders == total_tiles);
}

// TODO(Alexandris): Add asserts.
static void create_world_test1(framebuffer_t *framebuffer, world_t *world, camera_t *camera, settings_t *settings)
{
	framebuffer->width = 800;
	framebuffer->height = 400;
	framebuffer->bytes_per_pixel = 4;
	framebuffer->pitch = 500 * framebuffer->bytes_per_pixel;
	framebuffer->pixels = malloc(framebuffer->width * framebuffer->height * framebuffer->bytes_per_pixel);

	uint32_t total_spheres = 1;
	uint32_t total_materials = 1;

	world->total_spheres = total_spheres;
	world->spheres = (sphere_t*)malloc(sizeof(sphere_t) * total_spheres);
	world->spheres[0].center = point3_t(0.0f, 0.0f, -2.0f);
	world->spheres[0].radius = 1.0f;
	world->spheres[0].material_index = 0;

	world->total_materials = total_materials;
	world->materials = (material_t*)malloc(sizeof(material_t) * total_materials);
	world->materials[0].material = diffuse;
	world->materials[0].base_color = color_t{1.0f, 0.0f, 0.0f};
	world->materials[0].roughness = 0.0f;
	world->materials[0].ior = 0.0f;

	settings->samples_per_pixel = 100;
	settings->max_bounces = 10;
	settings->total_cpu_cores = get_total_cpu_cores();
	settings->total_cpu_cores = 4;
	settings->cache_line_size = get_cache_line_size();
	settings->lock_add = win32_lock_add;

	*camera = camera_t(point3_t{0.0f, 0.0f, 0.0f} , point3_t{0.0f, 0.0f, -1.0f} , vec3_t{0.0f, 1.0f, 0.0f}, 60.0f, framebuffer->width, framebuffer->height);
}

int main(int argc, char **argv)
{
	uint64_t cpu_freq = estimate_cpu_time_freq();
	uint64_t cpu_start = read_cpu_timer();
	win32_raytrace_code raytrace_code = load_raytracer_library("photon.dll");

	framebuffer_t framebuffer = {0};
	world_t world = {0};
	camera_t camera = {};
	settings_t settings = {0};

	create_world_test1(&framebuffer, &world, &camera, &settings);

	work_queue queue = {0};
	create_work_queue(&queue, &settings, &world, &framebuffer);

	win32_wrapper wrapper = {0};
	wrapper.framebuffer = &framebuffer;
	wrapper.world = &world;
	wrapper.camera = &camera;
	wrapper.settings = &settings;
	wrapper.queue = &queue;
	wrapper.raytracer = raytrace_code.raytracer;

	printf("Configuration: %d cores with %dx%d tiles (%dKByte/tile)\n", settings.total_cpu_cores,
			settings.total_tiles, settings.tile_width, settings.tile_height);
	printf("Quality: %d samples/pixel, %d bounces (max) per ray.\n", settings.samples_per_pixel,
			settings.max_bounces);

	win32_lock_add(&queue.work_order_index, 0);

	for(uint32_t i = 0; i < settings.total_cpu_cores; ++i)
	{
		create_worker_thread(&wrapper);
	}

#if 1
	while(queue.tiles_retired_counter < settings.total_tiles)
#else
	while(queue.work_order_index < queue.total_work_orders) // Experiment.
#endif
	{
		if(raytrace_code.raytracer(&framebuffer, &world, &camera, &settings, &queue))
		{
			fprintf(stderr, "\rRaycasting %d%%...", 100 * (uint32_t)queue.tiles_retired_counter / settings.total_tiles);
			fflush(stdout);
		}
		
	}


	ppm_image image = {0};
	convert_framebuffer_to_ppm_image(framebuffer.pixels, framebuffer.width, framebuffer.height, &image);
	PLATFORM_WRITE_FILE("test.ppm", image.total_size, image.data);

	uint64_t cpu_end = read_cpu_timer();

	printf("\nTotal execution time: %.2f(secs).\n", ((double)cpu_end - cpu_start)/(double)cpu_freq);

	free(queue.work_orders);
	free(world.materials);
	free(world.spheres);
	free(image.data);
	free(framebuffer.pixels);

	return 0;
}
