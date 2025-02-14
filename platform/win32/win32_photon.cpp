#include "win32_photon.h"
#include "./../../core/test_world.cpp"

static win32_photon win32_load_photon_dll(const char *dll)
{
	win32_photon result = {0};
	result.library = LoadLibraryA(dll);
	result.raytracer_iterative_sp = (raytrace_function_sp*)GetProcAddress(result.library, "raytracer_iterative_sp");
	result.raytracer_recursive_sp = (raytrace_function_sp*)GetProcAddress(result.library, "raytracer_recursive_sp");
	result.raytracer_queue_sp = (raytrace_function_queue_sp*)GetProcAddress(result.library, "raytracer_queue_sp");
	return result;
}

static void win32_unload_photon_dll(win32_photon *photon)
{
	if(photon)
	{
		if(photon->library)
		{
			FreeLibrary(photon->library);
			photon->library = 0;
			photon->raytracer_iterative_sp = 0;
			photon->raytracer_recursive_sp = 0;
            photon->raytracer_queue_sp = 0;
		}
	}
}

static void PLATFORM_WRITE_FILE(const char *file_name, uint32_t memory_size, void *memory)
{
	HANDLE file_handle = CreateFile(file_name, GENERIC_WRITE, 0, 0, OPEN_ALWAYS, 0, 0); 

	if(file_handle != INVALID_HANDLE_VALUE)
	{
		DWORD bytes_written;

		if(WriteFile(file_handle, memory, memory_size, &bytes_written, 0))
		{
			// result = (bytes_written == memory_size)
			CloseHandle(file_handle);
		}
		else
		{
			// TODO(Alexandris): Logging.
		}
	}
	else
	{
		// TODO(Alexandris): Logging.
	}

}

static void convert_framebuffer_to_image(void *data, uint32_t width, uint32_t height, image_t *image)
{
	switch(image->format)
	{
		case IMG_PPM:
		{
			char ppm_header[64]; // TODO(Alexandris): Assert that ppm_header < 64.
			snprintf(ppm_header, sizeof(ppm_header), "P6\n%d %d\n255\n", width, height);

			uint32_t ppm_header_size = (uint32_t)strlen(ppm_header);
			uint32_t total_size = (uint32_t)(ppm_header_size + width * height * sizeof(uint8_t) * 3);

			// TODO(Alexandris): Probably we should not allocate memory inside the function :-)
			void *image_data = malloc(total_size);
			memset(image_data, 0, total_size);
			memcpy(image_data, ppm_header, ppm_header_size);

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
		}  break;
		default:
		break;
	}
}

inline uint64_t get_os_timer_freq(void)
{
	LARGE_INTEGER freq;
	QueryPerformanceFrequency(&freq);
	return freq.QuadPart;
}

inline uint64_t read_os_timer(void)
{
	LARGE_INTEGER value;
	QueryPerformanceCounter(&value);
	return value.QuadPart;
}

inline uint64_t read_cpu_timer(void)
{
	return __rdtsc();
}

inline uint64_t estimate_cpu_freq(void)
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

uint64_t lock_add(uint64_t volatile *value, uint64_t add)
{
	uint64_t result = InterlockedExchangeAdd64((volatile int64_t*)value, add);
	return result;
}

DWORD WINAPI worker_thread(void *lpParameter)
{
    thread_params_t *params = (thread_params_t*)lpParameter;
	while(params->raytracer(params->queue, params->lock_add))
	{
	}

	return 0;
}

void create_worker_thread(thread_params_t *params)
{
	DWORD thread_id;
	HANDLE thread_handle = CreateThread(0, 0, worker_thread, params, 0, &thread_id);
	CloseHandle(thread_handle);
}

// Setup functions for SP and DP cameras
void setup_camera_sp(const camera_setup_t& setup, camera_t* camera)
{
    point3f_t look_from(13.0f, 2.0f, 3.0f);
    point3f_t look_at(0.0f, 0.0f, 0.0f);
    vec3f_t v_up(0.0f, 1.0f, 0.0f);
    
    initialize_camera32(setup.aspect_ratio, 
                       setup.width, 
                       setup.height, 
                       look_from, 
                       look_at, 
                       v_up, 
                       setup.defocus_angle, 
                       setup.focus_dist, 
                       setup.samples_per_pixel, 
                       setup.max_bounces, 
                       camera);
}

void run_raytracer_test_sp(const char *test_name, raytrace_function_sp raytracer_func, framebuffer_t *framebuffer,
        world_t *world, camera_t *camera, const uint64_t cpu_freq)
{
    info_t info = {0};
    uint64_t start_cpu_timer = __rdtsc();

    raytracer_func(framebuffer, world, camera, &info);

    // Convert and save image
    image_t image = {};
    image.format = IMG_PPM;
    convert_framebuffer_to_image(framebuffer->data, framebuffer->width, framebuffer->height, &image);

    char filename[256];
    snprintf(filename, sizeof(filename), "test_%s.ppm", test_name);
    PLATFORM_WRITE_FILE(filename, image.total_size, image.data);
    free(image.data);

    // Performance measurements
    uint64_t total_tsc_elapsed = __rdtsc() - start_cpu_timer;

    if(cpu_freq)
    {
        double time_ms = 1000.0 * (double)total_tsc_elapsed / (double)cpu_freq;
        printf("\n%s Raytracer total time: %0.4fms (estimated CPU freq %llu)\n", 
               test_name, time_ms, cpu_freq);
        printf("Bounces Computed: %lld bounces\n", info.bounces_computed);
        printf("Performance: %fms/bounce\n", 
               time_ms / (double)info.bounces_computed);
    }
}

void run_raytracer_queue_test_sp(const char *test_name, raytrace_function_queue_sp raytracer_func, framebuffer_t *framebuffer,
        world_t *world, camera_t *camera, const uint64_t cpu_freq, uint32_t cpu_cores)
{
    uint64_t start_cpu_timer = __rdtsc();

    work_queue_t work_queue = {0};

    uint32_t tile_width = (framebuffer->width  +  (cpu_cores - 1)) / cpu_cores;
    uint32_t tile_height = tile_width;
    uint32_t tile_count_x = (framebuffer->width  + tile_width - 1) / tile_width;
    uint32_t tile_count_y = (framebuffer->height + tile_height - 1)/ tile_height;
    uint32_t total_tiles = tile_count_x * tile_count_y;

    work_queue.work_orders = (work_order_t*)malloc(sizeof(work_order_t) * total_tiles);

    for(uint32_t tile_y = 0; tile_y < tile_count_y; ++tile_y)
    {
        uint32_t y_min = tile_y * tile_height;
        uint32_t y_max = y_min + tile_height;
        if(y_max > framebuffer->height)
            y_max = framebuffer->height;

        for(uint32_t tile_x = 0; tile_x < tile_count_x; ++tile_x)
        {
            uint32_t x_min = tile_x * tile_width;
            uint32_t x_max = x_min + tile_height;
            if(x_max > framebuffer->width)
                x_max = framebuffer->width;

            work_order_t *work_order = work_queue.work_orders + work_queue.work_orders_count++;
            // TODO(alexandris): add assert
            // PHOTON_ASSERT(work_queue.work_orders_count <= total_tiles);

            work_order->world = world;
            work_order->framebuffer = *framebuffer;
            work_order->x_min = x_min;
            work_order->x_max = x_max;
            work_order->y_min = y_min;
            work_order->y_max = y_max;
        }

    }

    work_queue.camera = *camera;

    // TODO(alexandris): add assert
    // PHOTON_ASSERT(work_queue.work_orders_count == total_tiles);

	printf("\nConfiguration: %d cores with %d %dx%d tiles (%dKByte/tile)\n",
			cpu_cores, total_tiles, tile_width, tile_height,
			tile_width * tile_height * 3 / 1024);
	printf("Quality: %d samples/pixel, %d bounces (max) per ray\n",
			camera->samples_per_pixel, camera->max_bounces);

	lock_add(&work_queue.work_order_index, 0);

    thread_params_t thread_params = {0};
    thread_params.queue = &work_queue;
    thread_params.lock_add = lock_add;
    thread_params.raytracer = raytracer_func;

	for(uint32_t core = 1; core < cpu_cores; ++core)
	{
		create_worker_thread(&thread_params);
	}

	while(work_queue.tile_retired_count < total_tiles)
	{
		if(raytracer_func(&work_queue, lock_add))
		{
			fprintf(stderr, "\rRaycasting %d%%...", 100 * (uint32_t)work_queue.tile_retired_count / total_tiles);
			fflush(stdout);
		}
	}

    // Convert and save image
    image_t image = {};
    image.format = IMG_PPM;
    convert_framebuffer_to_image(framebuffer->data, framebuffer->width, framebuffer->height, &image);

    char filename[256];
    snprintf(filename, sizeof(filename), "test_%s.ppm", test_name);
    PLATFORM_WRITE_FILE(filename, image.total_size, image.data);

    free(image.data);
    free(work_queue.work_orders);

    // Performance measurements
    uint64_t total_tsc_elapsed = __rdtsc() - start_cpu_timer;

    if(cpu_freq)
    {
        double time_ms = 1000.0 * (double)total_tsc_elapsed / (double)cpu_freq;
        printf("\n%s Raytracer total time: %0.4fms (estimated CPU freq %llu)\n", 
               test_name, time_ms, cpu_freq);
        printf("Bounces Computed: %lld bounces\n", work_queue.bounces_computed);
        printf("Performance: %fms/bounce\n", 
               time_ms / (double)work_queue.bounces_computed);
    }

}

int main(int argc, char **argv)
{
	win32_photon photon = win32_load_photon_dll("photon.dll");
	printf("Photon Library was loaded successfully.\n");

    // TODO(Alexandris): Support optional user input for choosing which raytracer or/and which test case

	float aspect_ratio = 16.0f / 9.0f;
	uint32_t image_width = 1200;
	uint32_t image_height = uint32_t(image_width / aspect_ratio);
	image_height = (image_height < 1) ? 1 : image_height;

    uint32_t samples_per_pixel = 10;
    uint32_t max_bounces = 50;

	framebuffer_t framebuffer = {0};
	framebuffer.width = image_width;
	framebuffer.height = image_height;
	framebuffer.bytes_per_pixel = sizeof(uint32_t);
	framebuffer.data = malloc(framebuffer.bytes_per_pixel * framebuffer.width * framebuffer.height);

    // Initialize common setup
    camera_setup_t camera_setup = {0};
    camera_setup.aspect_ratio = aspect_ratio;
    camera_setup.width = framebuffer.width;
    camera_setup.height = framebuffer.height;
    camera_setup.defocus_angle = 0.6f;
    camera_setup.focus_dist = 10.0f;
    camera_setup.samples_per_pixel = samples_per_pixel;
    camera_setup.max_bounces = max_bounces;

	camera_t camera;
    setup_camera_sp(camera_setup, &camera);

	world_t world = {0};
	init_world(&world);

    uint64_t cpu_freq = estimate_cpu_freq();

    // Run raytracing tests
    run_raytracer_test_sp("iterative_sp", photon.raytracer_iterative_sp, 
                         &framebuffer, &world, &camera, cpu_freq);
    
    
    run_raytracer_test_sp("recursive_sp", photon.raytracer_recursive_sp, 
                         &framebuffer, &world, &camera, cpu_freq);
    

	free(world.spheres);
	free(world.materials);

    create_world_case_obj(&world,
        aspect_ratio,
        framebuffer.width,
        framebuffer.height,
        50,
        10,
        &camera);


    run_raytracer_test_sp("iterative_sp_obj", photon.raytracer_iterative_sp, 
                         &framebuffer, &world, &camera, cpu_freq);

    run_raytracer_queue_test_sp("iterative_queue_sp_obj", photon.raytracer_queue_sp, &framebuffer,
            &world, &camera, cpu_freq, 4);

    free(world.triangles);
	free(world.materials);

	free(framebuffer.data);
	win32_unload_photon_dll(&photon);
	return 0;
}
