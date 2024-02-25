#ifndef CORE_TEST_WORLDS_H
#define CORE_TEST_WORLDS_H

/*
 * create_world_case1: Creates an image with a sphere that uses ideal diffuse material.
 * create_world_case2: Renders an object from an .obj file.
 */

// TODO(Alexandris): Add asserts.
static void create_world_case1(framebuffer_t *framebuffer, world_t *world, camera_t *camera, settings_t *settings)
{
	framebuffer->width = 800;
	framebuffer->height = 400;
	framebuffer->bytes_per_pixel = 4;
	framebuffer->pitch = 500 * framebuffer->bytes_per_pixel;
	framebuffer->pixels = malloc(framebuffer->width * framebuffer->height * framebuffer->bytes_per_pixel);

	constexpr uint32_t total_spheres = 1;
	constexpr uint32_t total_materials = 1;

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

	/*
	settings->samples_per_pixel = 100;
	settings->max_bounces = 10;
	settings->total_cpu_cores = get_total_cpu_cores();
	settings->total_cpu_cores = 4;
	settings->cache_line_size = get_cache_line_size();
	settings->lock_add = win32_lock_add;
	*/

	*camera = camera_t(point3_t{0.0f, 0.0f, 0.0f} , point3_t{0.0f, 0.0f, -1.0f} , vec3_t{0.0f, 1.0f, 0.0f}, 60.0f, framebuffer->width, framebuffer->height);
}

static void create_world_case2(framebuffer_t *framebuffer, world_t *world, camera_t *camera, settings_t *settings)
{
	framebuffer->width = 800;
	framebuffer->height = 400;
	framebuffer->bytes_per_pixel = 4;
	framebuffer->pitch = 500 * framebuffer->bytes_per_pixel;
	framebuffer->pixels = malloc(framebuffer->width * framebuffer->height * framebuffer->bytes_per_pixel);

	char const *filename = "monkey_smooth.obj";
	ObjResult obj_result;
	memset(&obj_result, 0, sizeof(ObjResult));
	ObjInfo obj_info;
	memset(&obj_info, 0, sizeof(ObjInfo));
// 	int32_t result =
// 	TODO(Alexandris): read value returned from read_obj_file(...) function and add paths accordinly to that.
	read_obj_file(filename, &obj_result, &obj_info);

	constexpr uint32_t total_materials = 1;

	world->total_materials = total_materials;
	world->materials = (material_t*)malloc(sizeof(material_t) * total_materials);
	world->materials[0].material = diffuse;
	world->materials[0].base_color = color_t{1.0f, 0.0f, 0.0f};
	world->materials[0].roughness = 0.0f;
	world->materials[0].ior = 0.0f;

	world->total_triangles = (uint32_t)(obj_result.index_buffer_size / 3);
	world->triangles = (triangle_t*)malloc(sizeof(triangle_t) * world->total_triangles);

	size_t index = 0;
	constexpr float z_offset = -2.2f;

	uint16_t *index_buffer = (uint16_t*)obj_result.index_buffer;

	for(size_t i = 0; i < world->total_triangles; i++)
	{
		world->triangles[i].vertices[0].x = obj_result.vertex_buffer[index_buffer[index]].point.x;
		world->triangles[i].vertices[0].y = obj_result.vertex_buffer[index_buffer[index]].point.y;
		world->triangles[i].vertices[0].z = obj_result.vertex_buffer[index_buffer[index]].point.z + z_offset;

		world->triangles[i].vertices[1].x = obj_result.vertex_buffer[index_buffer[index + 1]].point.x;
		world->triangles[i].vertices[1].y = obj_result.vertex_buffer[index_buffer[index + 1]].point.y;
		world->triangles[i].vertices[1].z = obj_result.vertex_buffer[index_buffer[index + 1]].point.z + z_offset;

		world->triangles[i].vertices[2].x = obj_result.vertex_buffer[index_buffer[index + 2]].point.x;
		world->triangles[i].vertices[2].y = obj_result.vertex_buffer[index_buffer[index + 2]].point.y;
		world->triangles[i].vertices[2].z = obj_result.vertex_buffer[index_buffer[index + 2]].point.z + z_offset;

		world->triangles[i].material_index = 0;

		index += 3;
	}

	PHOTON_ASSERT(index == obj_result.index_buffer_size);

	free(obj_result.data);

	/*
	settings->samples_per_pixel = 100;
	settings->max_bounces = 10;
	settings->total_cpu_cores = get_total_cpu_cores();
	settings->total_cpu_cores = 4;
	settings->cache_line_size = get_cache_line_size();
	settings->lock_add = win32_lock_add;
	*/

	*camera = camera_t(point3_t{0.0f, 0.0f, 0.0f} , point3_t{0.0f, 0.0f, -1.0f} , vec3_t{0.0f, 1.0f, 0.0f}, 60.0f, framebuffer->width, framebuffer->height);
}


#endif
