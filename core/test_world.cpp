#include "test_world.h"

#include "../vendor/obj.h"

void init_world(world_t *world)
{
    uint32_t total_spheres = 485;
    uint32_t total_materials = 485;

    world->total_spheres = total_spheres;
    world->spheres = (sphere_t *)malloc(sizeof(sphere_t) * total_spheres);
    world->total_triangles = 0;
    world->triangles = 0;
    world->total_materials = total_materials;
    world->materials = (material_t *)malloc(sizeof(material_t) * total_materials);

    memset(world->spheres, 0, sizeof(sphere_t) * total_spheres);
    memset(world->materials, 0, sizeof(material_t) * total_materials);

    uint32_t current_material_index = 0;
    uint32_t current_sphere_index = 0;

    // Ground material and sphere
    // Ground sphere
    world->materials[current_material_index].material = MAT_DIFFUSE;
    world->materials[current_material_index].base_color = color3f_t(0.5f, 0.5f, 0.5f);
    world->materials[current_material_index].roughness = 0.0f;
    world->materials[current_material_index].refraction_index = 0.0f;

    world->spheres[current_sphere_index].center = point3f_t(0.0f, -1000.0f, 0.0f);
    world->spheres[current_sphere_index].radius = 1000.0f;
    world->spheres[current_sphere_index].material_index = current_material_index;
    
	current_material_index++;
    current_sphere_index++;

    // Small random spheres
    for (int a = -11; a < 11; ++a)
	{
        for (int b = -11; b < 11; ++b)
		{
            float choose_mat = random_float();
            point3f_t center = {
                a + 0.9f * random_float(),
                0.2f,
                b + 0.9f * random_float()
			};
            
			vec3f_t tmp = center - point3f_t(4.0f, 0.2f, 0.0f);

            if (dot(tmp, tmp) > 0.9f)
			{
                if (choose_mat < 0.8f)
				{
                    // Diffuse
					vec3f_t v1 = random_vec3f();
					vec3f_t v2 = random_vec3f();
                    world->materials[current_material_index].material = MAT_DIFFUSE;
                    world->materials[current_material_index].base_color = color3f_t(v1.x * v2.x, v1.y * v2.y, v1.z * v2.z);
                    world->materials[current_material_index].roughness = 0.0f;
                    world->materials[current_material_index].refraction_index = 0.0f;
				} else if (choose_mat < 0.95f) {
                    // Metal
					vec3f_t v1 = random(0.5f, 1.0f);
                    world->materials[current_material_index].material = MAT_METAL;
                    world->materials[current_material_index].base_color = color3f_t(v1.x, v1.y, v1.z);
                    world->materials[current_material_index].roughness = random_float(0.0f, 0.5f);
                    world->materials[current_material_index].refraction_index = 0.0f;
                } else {
                    // Glass
					world->materials[current_material_index].material = MAT_DIELECTRIC;
                    world->materials[current_material_index].base_color = color3f_t(1.0f, 1.0f, 1.0f);
                    world->materials[current_material_index].roughness = 0.0f;
                    world->materials[current_material_index].refraction_index = 1.5f;
                }
                world->spheres[current_sphere_index].center = center;
                world->spheres[current_sphere_index].radius = 0.2f;
                world->spheres[current_sphere_index].material_index = current_material_index;

                current_material_index++;
                current_sphere_index++;
            }
        }
    }

    // Add three large spheres
    // Three larger spheres
    world->materials[current_material_index].material = MAT_DIELECTRIC;
    world->materials[current_material_index].base_color = color3f_t(0.8f, 0.9f, 1.0f);
    world->materials[current_material_index].roughness = 0.0f;
    world->materials[current_material_index].refraction_index = 1.5f;

    world->spheres[current_sphere_index].center = point3f_t(0.0f, 1.0f, 0.0f);
    world->spheres[current_sphere_index].radius = 1.0f;
    world->spheres[current_sphere_index].material_index = current_material_index;
    
	current_material_index++;
    current_sphere_index++;

    world->materials[current_material_index].material = MAT_DIFFUSE;
    world->materials[current_material_index].base_color = color3f_t(0.4f, 0.2f, 0.1f);
    world->materials[current_material_index].roughness = 0.2f;
    world->materials[current_material_index].refraction_index = 0.0f;
    
	world->spheres[current_sphere_index].center = point3f_t(-4.0f, 1.0f, 0.0f);
    world->spheres[current_sphere_index].radius = 1.0f;
    world->spheres[current_sphere_index].material_index = current_material_index;
    
	current_material_index++;
    current_sphere_index++;

    world->materials[current_material_index].material = MAT_METAL;
    world->materials[current_material_index].base_color = color3f_t(0.7f, 0.6f, 0.5f);
    world->materials[current_material_index].roughness = 0.0f;
    world->materials[current_material_index].refraction_index = 0.0f;
    
	world->spheres[current_sphere_index].center = point3f_t(4.0f, 1.0f, 0.0f);
    world->spheres[current_sphere_index].radius = 1.0f;
    world->spheres[current_sphere_index].material_index = current_material_index;
    
	world->total_materials = ++current_material_index;
    world->total_spheres = ++current_sphere_index;
}

void create_world_case_obj(world_t *world,
        float aspect_ratio,
        uint32_t framebuffer_width,
        uint32_t framebuffer_height,
        uint32_t samples_per_pixel,
        uint32_t max_bounces,
        camera_t *camera)
{
	// PHOTON_ASSERT(world != 0 && settings != 0);

    point3f_t look_from(0.0f, 0.0f, 8.0f);
    point3f_t look_at(0.0f, 0.0f, -1.0f);
    vec3f_t v_up(0.0f, 1.0f, 0.0f);

    camera_setup_t camera_setup = {0};
    camera_setup.aspect_ratio = aspect_ratio;
    camera_setup.width = framebuffer_width;
    camera_setup.height = framebuffer_height;
    camera_setup.defocus_angle = 0.0f;
    camera_setup.focus_dist = 10.0f;
    camera_setup.samples_per_pixel = samples_per_pixel;
    camera_setup.max_bounces = max_bounces;
    
    initialize_camera32(camera_setup.aspect_ratio, 
                       camera_setup.width, 
                       camera_setup.height, 
                       look_from, 
                       look_at, 
                       v_up, 
                       camera_setup.defocus_angle, 
                       camera_setup.focus_dist, 
                       camera_setup.samples_per_pixel, 
                       camera_setup.max_bounces, 
                       camera);

	char const *filename = "monkey_smooth.obj";
	obj_result_t obj_result;
	memset(&obj_result, 0, sizeof(obj_result_t));
	obj_info_t obj_info;
	memset(&obj_info, 0, sizeof(obj_info_t));
	read_obj_file(filename, &obj_result, &obj_info);

	constexpr uint32_t materials_count = 1;

    world->total_spheres = 0;
	world->total_materials = materials_count;
	world->materials = (material_t*)malloc(sizeof(material_t) * materials_count);

	// PHOTON_ASSERT((obj_result.index_buffer_size % 3) == 0);

	world->total_triangles = (uint32_t)(obj_result.index_buffer_size / 3);
	world->triangles = (triangle_t*)malloc(sizeof(triangle_t) * world->total_triangles);

	material_t material1 = {};
	material1.material = MAT_DIFFUSE;
	material1.base_color = color3f_t(0.753f, 0.647f, 0.435f);

	world->materials[0] = material1;

	size_t index = 0;
	float z_offset = 0.0f;

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

	// PHOTON_ASSERT(index == obj_result.index_buffer_size);

	free(obj_result.data);
}
