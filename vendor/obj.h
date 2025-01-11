#ifndef OBJ_H
#define OBJ_H

#include <stdio.h> // fread, printf
#include <stdint.h> // int32_t, ...
#include <stdlib.h> // malloc, free
#include <string.h> // memset
#include <assert.h> // assert

// http://paulbourke.net/dataformats/obj/

// TODO: Add asserts where is needed.
//
//
//
/*
 * obj_info.index_type: 0 = uint16_t buffer, 1 = uint32_t buffer, 2 = no index buffer
	char const *filename = "assets/viking_room.obj";
	ObjResult obj_result = {0};
	// memset(&obj_result, 0, sizeof(ObjResult));
	ObjInfo obj_info = {0};
	// memset(&obj_info, 0, sizeof(ObjInfo));
	// obj_info.index_type = 2;
	int32_t result = read_obj_file(filename, &obj_result, &obj_info);

	for(uint32_t i = 0; i < obj_result.vertex_buffer_size; ++i)
	{
		printf("v %f %f %f \n", obj_result.vertex_buffer[i].point.x, obj_result.vertex_buffer[i].point.y, obj_result.vertex_buffer[i].point.z);
	}

	for(uint32_t i = 0; i < obj_result.index_buffer_size; ++i)
	{
	}

	for(uint32_t i = 0; i < obj_result.material_buffer_size; ++i)
	{
		obj_result.material_buffer[i];
	}
	free(obj_result.data);
 */

typedef struct
{
    float x, y, z;
} obj_point3_t;

typedef struct
{
    float u, v;
} obj_texture_coordinate_t;

typedef struct
{
    float x, y, z;
} obj_normal_t;

typedef struct
{
	uint32_t v;
	uint32_t uv;
	uint32_t n;
	uint32_t m;
} obj_face_t;

typedef struct
{
    float red, green, blue;
} obj_color_t;

typedef struct
{
	obj_color_t ambient; // Ka
	obj_color_t diffuse; // Kd
	obj_color_t specular; // Ks
	obj_color_t emissive; // Ke
	obj_color_t transmission_filter_color; // Tf
	float specular_exponent; // Ns
	float ior; // Ni
	float transparency; // d
	uint32_t illumination_mode; // illum
	// char texture_map[50]; // map_Ka, map_Kd, map_Ks, etc...
} obj_material_t;

typedef struct
{
    obj_point3_t point;
    obj_texture_coordinate_t texture_coordinate;
    obj_normal_t normal;
    uint32_t material_index;
} obj_vertex_t;

typedef struct
{
	uint32_t total_vertices;
	uint32_t total_points;
	uint32_t total_uvs;
	uint32_t total_normals;
	uint32_t total_materials;
	// char material[128];
    uint8_t is_initialized;
    uint8_t index_type; // 0 = uint16_t, 1 = uint32_t, 2 = no index buffer
} obj_info_t;

typedef struct
{
	void *data;
	obj_vertex_t *vertex_buffer;
	void *index_buffer;
	obj_material_t *material_buffer;
	uint32_t vertex_buffer_size;
	uint32_t index_buffer_size;
	uint32_t material_buffer_size;
} obj_result_t;

// data should point at the first digit
static float parse_float(const char *data, uint32_t *consumed_size) 
{
    if (!data || !*data)
        return 0.0f;

    const char *start = data;
    float result = 0.0f;
    bool is_negative = *data == '-';
    if (is_negative)
        ++data;

    // Parse integer part
    while (*data >= '0' && *data <= '9')
    {
        result = result * 10.0f + (*data - '0');
        data++;
    }

    // Parse decimal part
    if (*data == '.')
    {
        data++;
        float factor = 0.1f;
        while (*data >= '0' && *data <= '9') {
            result += (*data - '0') * factor;
            factor *= 0.1f;
            data++;
        }
    }

    if (consumed_size)
        *consumed_size = (uint32_t)(data - start);

    return is_negative ? -result : result;
}

static int32_t parse_integer(const char *data, uint32_t *consumed_size)
{
    if (!data || !*data)
        return 0;

    const char *start = data;
    int32_t result = 0;
    bool is_negative = *data == '-';
    if (is_negative)
        ++data;

    while (*data >= '0' && *data <= '9')
    {
        result = result * 10 + (*data - '0');
        data++;
    }

    if (consumed_size)
        *consumed_size = (uint32_t)(data - start);

    return is_negative ? -result : result;
}

static int32_t read_obj_file(const char *filename, obj_result_t *obj_result, obj_info_t *obj_info)
{
	assert(obj_result != 0 && obj_info != 0);

	obj_info->is_initialized = 0;

	FILE *file;
	fopen_s(&file, filename, "rb");
	if(!file)
	{
		printf("ERROR: Failed to open file: %s\n", filename);
		return 1;
	}

	fseek(file, 0L, SEEK_END);
	uint32_t file_size = ftell(file);
	fseek(file, 0L, SEEK_SET);

	char *data = (char*)malloc(sizeof(char) * file_size);
	fread(data, file_size, 1, file);
	char *current_ptr = data;
	uint32_t current_size = file_size;

	uint32_t total_faces = 0;
	char material_file[64];

	while(current_size)
	{
		if(current_ptr[0] == 'v' && current_ptr[1] == ' ')
		{
			++obj_info->total_points;
			obj_info->is_initialized |= 0x1;
		}
		else if(current_ptr[0] == 'v' && current_ptr[1] == 't')
		{
			++obj_info->total_uvs;
			obj_info->is_initialized |= 0x2;
		}
		else if(current_ptr[0] == 'v' && current_ptr[1] == 'n')
		{
			++obj_info->total_normals;
			obj_info->is_initialized |= 0x4;
		}
		else if(current_ptr[0] == 'f' && current_ptr[1] == ' ')
		{
			++total_faces;
		}
		else if(current_ptr[0] == 'u' && current_ptr[1] == 's' && current_ptr[2] == 'e' && current_ptr[3] == 'm' && current_ptr[4] == 't' && current_ptr[5] == 'l')
		{
			++obj_info->total_materials;
		}
		else if(current_ptr[0] == 'm' && current_ptr[1] == 't' && current_ptr[2] == 'l' && current_ptr[3] == 'l' && current_ptr[4] == 'i' && current_ptr[5] == 'b')
		{
			current_ptr += 7;
			current_size -= 7;
			uint32_t length = 0;
			char *temp = current_ptr;
			while(*temp++ != '\n')
			{
				++length;
			}
			memcpy(material_file, current_ptr, length);
			material_file[length] = '\0';
		}

		while(current_ptr[0] != '\n')
		{
			++current_ptr;
			--current_size;

			if(!current_size) break;
		}

		++current_ptr;
		--current_size;
	}

	if(obj_info->total_materials) ++obj_info->total_materials;

	uint32_t face_multiplier = 10;
	void *memory = malloc(obj_info->total_points * sizeof(obj_point3_t) + obj_info->total_uvs * sizeof(obj_texture_coordinate_t) + obj_info->total_normals * sizeof(obj_normal_t) + total_faces * face_multiplier * sizeof(obj_face_t));

	obj_point3_t *points = (obj_point3_t*)memory;
	obj_texture_coordinate_t *texture_coordinates = (obj_texture_coordinate_t*)((uint8_t*)memory + obj_info->total_points * sizeof(obj_point3_t));
	obj_normal_t *normals = (obj_normal_t*)((uint8_t*)memory + obj_info->total_points * sizeof(obj_point3_t) + obj_info->total_uvs * sizeof(obj_texture_coordinate_t));
	obj_face_t *faces = (obj_face_t*)((uint8_t*)memory + obj_info->total_points * sizeof(obj_point3_t) + obj_info->total_uvs * sizeof(obj_texture_coordinate_t) + obj_info->total_normals * sizeof(obj_normal_t));

	fseek(file, 0, SEEK_SET);
	uint32_t v_index = 0;
	uint32_t uv_index = 0;
	uint32_t n_index = 0;
	uint32_t m_index = 0;
	uint32_t f_index = 0;

	current_ptr = data;
	current_size = file_size;
	uint32_t consumed_size;

	// TODO: Add aseert to check if faces_counter exceeds temp_faces array size.
	obj_face_t temp_faces[16];
	int32_t faces_counter;
	uint32_t start;

	while(current_size)
	{
		if(current_ptr[0] == 'v' && current_ptr[1] == ' ')
		{
			current_ptr += 2;
			current_size -= 2;

			points[v_index].x = parse_float(current_ptr, &consumed_size);
			current_ptr += consumed_size + 1; // the +1 removes the '/' between values.
			current_size -= consumed_size + 1;
			points[v_index].y = parse_float(current_ptr, &consumed_size);
			current_ptr += consumed_size + 1; // the +1 removes the '/' between values.
			current_size -= consumed_size + 1;
			points[v_index].z = parse_float(current_ptr, &consumed_size);
			current_ptr += consumed_size;
			current_size -= consumed_size;
			++v_index;
		}
		else if(current_ptr[0] == 'v' && current_ptr[1] == 't')
		{
			current_ptr += 3;
			current_size -= 3;

			texture_coordinates[uv_index].u = parse_float(current_ptr, &consumed_size);
			current_ptr += consumed_size + 1;
			current_size -= consumed_size + 1;
			texture_coordinates[uv_index].v = parse_float(current_ptr, &consumed_size);
			current_ptr += consumed_size;
			current_size -= consumed_size;
			++uv_index;
		}
		else if(current_ptr[0] == 'v' && current_ptr[1] == 'n')
		{
			current_ptr += 3;
			current_size -= 3;

			normals[n_index].x = parse_float(current_ptr, &consumed_size);
			current_ptr += consumed_size + 1;
			current_size -= consumed_size + 1;
			normals[n_index].y = parse_float(current_ptr, &consumed_size);
			current_ptr += consumed_size + 1;
			current_size -= consumed_size + 1;
			normals[n_index].z = parse_float(current_ptr, &consumed_size);
			current_ptr += consumed_size;
			current_size -= consumed_size;
			++n_index;
		}
		else if(current_ptr[0] == 'u' && current_ptr[1] == 's' && current_ptr[2] == 'e' && current_ptr[3] == 'm' && current_ptr[4] == 't' && current_ptr[5] == 'l')
		{
			++m_index;
		}
		else if(current_ptr[0] == 'f' && current_ptr[1] == ' ')
		{
			current_ptr += 2;
			current_size -= 2;

			switch(obj_info->is_initialized)
			{
				case 0x7: // vertex + uv + normal
				{
					faces_counter = 0;
					start = 0;

					while(*current_ptr != '\n')
					{
						temp_faces[faces_counter].v = parse_integer(current_ptr, &consumed_size);
						current_ptr += consumed_size + 1;
						current_size -= consumed_size + 1;
						temp_faces[faces_counter].uv = parse_integer(current_ptr, &consumed_size);
						current_ptr += consumed_size + 1;
						current_size -= consumed_size + 1;
						temp_faces[faces_counter].n = parse_integer(current_ptr, &consumed_size);
						current_ptr += consumed_size;
						current_size -= consumed_size;

						if(*current_ptr == ' ' || *current_ptr == '\r')
						{
							++current_ptr;
							--current_size;
						}

						temp_faces[faces_counter].m = m_index;
						++faces_counter;
					}

					int32_t c_index = 2;

					for(;;)
					{
						faces[f_index].v = temp_faces[0].v;
						faces[f_index].uv = temp_faces[0].uv;
						faces[f_index].n = temp_faces[0].n;
						faces[f_index].m = temp_faces[0].m;
						++f_index;

						faces[f_index].v = temp_faces[start + 1].v;
						faces[f_index].uv = temp_faces[start + 1].uv;
						faces[f_index].n = temp_faces[start + 1].n;
						faces[f_index].m = temp_faces[start + 1].m;
						++f_index;

						faces[f_index].v = temp_faces[start + 2].v;
						faces[f_index].uv = temp_faces[start + 2].uv;
						faces[f_index].n = temp_faces[start + 2].n;
						faces[f_index].m = temp_faces[start + 2].m;
						++f_index;

						++c_index;
						if(c_index != faces_counter)
                            ++start;
						else
                            break;
					}
				} break;
				case 0x5: // vertex + normal
				{
					faces_counter = 0;
					start = 0;

					while(*current_ptr != '\n')
					{
						temp_faces[faces_counter].v = parse_integer(current_ptr, &consumed_size);
						current_ptr += consumed_size + 2;
						current_size -= consumed_size + 2;
						temp_faces[faces_counter].n = parse_integer(current_ptr, &consumed_size);
						current_ptr += consumed_size;
						current_size -= consumed_size;

						if(*current_ptr == ' ' || *current_ptr == '\r')
						{
							++current_ptr;
							--current_size;
						}

						temp_faces[faces_counter].m = m_index;
						++faces_counter;
					}

					int32_t c_index = 2;

					for(;;)
					{
						faces[f_index].v = temp_faces[0].v;
						faces[f_index].n = temp_faces[0].n;
						faces[f_index].m = temp_faces[0].m;
						++f_index;

						faces[f_index].v = temp_faces[start + 1].v;
						faces[f_index].n = temp_faces[start + 1].n;
						faces[f_index].m = temp_faces[start + 1].m;
						++f_index;

						faces[f_index].v = temp_faces[start + 2].v;
						faces[f_index].n = temp_faces[start + 2].n;
						faces[f_index].m = temp_faces[start + 2].m;
						++f_index;

						++c_index;
						if(c_index != faces_counter) ++start;
						else break;
					}
				} break;
				case 0x3: // vertex + uv
				{
					// TODO: Implement this scenario.
				} break;
				default: // vertex
					faces_counter = 0;
					start = 0;

					while(*current_ptr != '\n')
					{
						temp_faces[faces_counter].v = parse_integer(current_ptr, &consumed_size);
						current_ptr += consumed_size;
						current_size -= consumed_size;

						if(*current_ptr == ' ' || *current_ptr == '\r')
						{
							++current_ptr;
							--current_size;
						}

						temp_faces[faces_counter].m = m_index;
						++faces_counter;
					}

					int32_t c_index = 2;

					for(;;)
					{
						faces[f_index].v = temp_faces[0].v;
						faces[f_index].n = temp_faces[0].n;
						faces[f_index].m = temp_faces[0].m;
						++f_index;

						faces[f_index].v = temp_faces[start + 1].v;
						faces[f_index].n = temp_faces[start + 1].n;
						faces[f_index].m = temp_faces[start + 1].m;
						++f_index;

						faces[f_index].v = temp_faces[start + 2].v;
						faces[f_index].n = temp_faces[start + 2].n;
						faces[f_index].m = temp_faces[start + 2].m;
						++f_index;

						++c_index;
						if(c_index != faces_counter) ++start;
						else break;
					}
			}
		}

		while(*current_ptr != '\n')
		{
			++current_ptr;
			--current_size;

			if(!current_size) break;
		}

		++current_ptr;
		--current_size;
	}

	if(obj_info->index_type == 2) // No index buffer.
	{
		obj_result->data = malloc(f_index * sizeof(obj_vertex_t) + obj_info->total_materials * sizeof(obj_material_t));
		memset(obj_result->data, 0, f_index * sizeof(obj_vertex_t) + obj_info->total_materials * sizeof(obj_material_t)); // NOTE: Maybe its not needed.
		obj_result->vertex_buffer = (obj_vertex_t*)obj_result->data;
		obj_result->vertex_buffer_size = f_index;
		obj_result->material_buffer = (obj_material_t*)((char*)obj_result->data + f_index * sizeof(obj_vertex_t));
		obj_result->material_buffer_size = obj_info->total_materials;

		switch(obj_info->is_initialized)
		{
			case 0x7: // vertex + uv + normal
			{
				uint32_t counter = 0;
				for(uint32_t i = 0; i < f_index; ++i)
				{
					v_index = faces[i].v - 1;
					uv_index = faces[i].uv - 1;
					n_index = faces[i].n - 1;
					m_index = faces[i].m - 1;
					
					obj_result->vertex_buffer[i].point = points[v_index];
					obj_result->vertex_buffer[i].texture_coordinate = texture_coordinates[uv_index];
					obj_result->vertex_buffer[i].normal = normals[n_index];
					obj_result->vertex_buffer[i].material_index = m_index;
					++counter;
				}
			} break;
			case 0x5: // vertex + normal
			{
				uint32_t counter = 0;
				for(uint32_t i = 0; i < f_index; ++i)
				{
					v_index = faces[i].v - 1;
					n_index = faces[i].n - 1;
					m_index = faces[i].m - 1;
					
					obj_result->vertex_buffer[i].point = points[v_index];
					obj_result->vertex_buffer[i].normal = normals[n_index];
					obj_result->vertex_buffer[i].material_index = m_index;
					++counter;
				}
			} break;
			case 0x3: // vertex + uv
			{
				uint32_t counter = 0;
				for(uint32_t i = 0; i < f_index; ++i)
				{
					v_index = faces[i].v - 1;
					uv_index = faces[i].uv - 1;
					m_index = faces[i].m - 1;
					
					obj_result->vertex_buffer[i].point = points[v_index];
					obj_result->vertex_buffer[i].material_index = m_index;
					++counter;
				}
			} break;
			default: // vertex
				uint32_t counter = 0;
				for(uint32_t i = 0; i < f_index; ++i)
				{
					v_index = faces[i].v - 1;
					m_index = faces[i].m - 1;
					
					obj_result->vertex_buffer[i].point = points[v_index];
					obj_result->vertex_buffer[i].texture_coordinate = texture_coordinates[uv_index];
					obj_result->vertex_buffer[i].material_index = m_index;
					++counter;
				}
		}
	}
	else
	{
		if(obj_info->index_type == 0) // uint16_t index_buffer
		{
			obj_result->data = malloc(obj_info->total_points * sizeof(obj_vertex_t) + f_index * sizeof(uint16_t) + obj_info->total_materials * sizeof(obj_material_t));
			memset(obj_result->data, 0, obj_info->total_points * sizeof(obj_vertex_t) + f_index * sizeof(uint16_t) + obj_info->total_materials * sizeof(obj_material_t));

			obj_result->vertex_buffer = (obj_vertex_t*)obj_result->data;
			obj_result->vertex_buffer_size = obj_info->total_points;

			obj_result->index_buffer = (char*)obj_result->data + obj_info->total_points * sizeof(obj_vertex_t);
			obj_result->index_buffer_size = f_index;

			obj_result->material_buffer = (obj_material_t*)((char*)obj_result->data + obj_info->total_points * sizeof(obj_vertex_t) + f_index * sizeof(uint16_t));
			obj_result->material_buffer_size = obj_info->total_materials;

			switch(obj_info->is_initialized)
			{
				case 0x7: // vertex + uv + normal
				{
					uint32_t counter = 0;
					for(uint32_t i = 0; i < f_index; ++i)
					{
						v_index = faces[i].v - 1;
						uv_index = faces[i].uv - 1;
						n_index = faces[i].n - 1;
						m_index = faces[i].m;

						obj_result->vertex_buffer[v_index].point = points[v_index];
						obj_result->vertex_buffer[v_index].texture_coordinate = texture_coordinates[uv_index];
						obj_result->vertex_buffer[v_index].normal = normals[n_index];
						obj_result->vertex_buffer[v_index].material_index = m_index;

						((uint16_t*)obj_result->index_buffer)[counter++] = (uint16_t)v_index;
					}
				} break;
				case 0x5: // vertex + normal
				{
					uint32_t counter = 0;
					for(uint32_t i = 0; i < f_index; ++i)
					{
						v_index = faces[i].v - 1;
						n_index = faces[i].n - 1;
						m_index = faces[i].m;

						obj_result->vertex_buffer[v_index].point = points[v_index];
						obj_result->vertex_buffer[v_index].normal = normals[n_index];
						obj_result->vertex_buffer[v_index].material_index = m_index;

						((uint16_t*)obj_result->index_buffer)[counter++] = (uint16_t)v_index;
					}
				} break;
				case 0x3: // vertex + uv
				{
					uint32_t counter = 0;
					for(uint32_t i = 0; i < f_index; ++i)
					{
						v_index = faces[i].v - 1;
						uv_index = faces[i].uv - 1;
						m_index = faces[i].m;

						obj_result->vertex_buffer[v_index].point = points[v_index];
						obj_result->vertex_buffer[v_index].texture_coordinate = texture_coordinates[uv_index];
						obj_result->vertex_buffer[v_index].material_index = m_index;

						((uint16_t*)obj_result->index_buffer)[counter++] = (uint16_t)v_index;
					}
				} break;
				default: // vertex
					uint32_t counter = 0;
					for(uint32_t i = 0; i < f_index; ++i)
					{
						v_index = faces[i].v - 1;
						m_index = faces[i].m;

						obj_result->vertex_buffer[v_index].point = points[v_index];
						obj_result->vertex_buffer[v_index].material_index = m_index;

						((uint16_t*)obj_result->index_buffer)[counter++] = (uint16_t)v_index;
					}
			}
		}
		else // uint32_t index_buffer
		{
			obj_result->data = malloc(obj_info->total_points * sizeof(obj_vertex_t) + f_index * sizeof(uint32_t) + obj_info->total_materials * sizeof(obj_material_t));
			memset(obj_result->data, 0, obj_info->total_points * sizeof(obj_vertex_t) + f_index * sizeof(uint32_t) + obj_info->total_materials * sizeof(obj_material_t));

			obj_result->vertex_buffer = (obj_vertex_t*)obj_result->data;
			obj_result->vertex_buffer_size = obj_info->total_points;

			obj_result->index_buffer = (char*)obj_result->data + obj_info->total_points * sizeof(obj_vertex_t);
			obj_result->index_buffer_size = f_index;

			obj_result->material_buffer = (obj_material_t*)((char*)obj_result->data + obj_info->total_points * sizeof(obj_vertex_t) + f_index * sizeof(uint32_t));
			obj_result->material_buffer_size = obj_info->total_materials;

			switch(obj_info->is_initialized)
			{
				case 0x7: // vertex + uv + normal
				{
					uint32_t counter = 0;
					for(uint32_t i = 0; i < f_index; ++i)
					{
						v_index = faces[i].v - 1;
						uv_index = faces[i].uv - 1;
						n_index = faces[i].n - 1;
						m_index = faces[i].m;

						obj_result->vertex_buffer[v_index].point = points[v_index];
						obj_result->vertex_buffer[v_index].texture_coordinate = texture_coordinates[uv_index];
						obj_result->vertex_buffer[v_index].normal = normals[n_index];
						obj_result->vertex_buffer[v_index].material_index = m_index;

						((uint32_t*)obj_result->index_buffer)[counter++] = (uint32_t)v_index;
					}
				} break;
				case 0x5: // vertex + normal
				{
					uint32_t counter = 0;
					for(uint32_t i = 0; i < f_index; ++i)
					{
						v_index = faces[i].v - 1;
						n_index = faces[i].n - 1;
						m_index = faces[i].m;

						obj_result->vertex_buffer[v_index].point = points[v_index];
						obj_result->vertex_buffer[v_index].normal = normals[n_index];
						obj_result->vertex_buffer[v_index].material_index = m_index;

						((uint32_t*)obj_result->index_buffer)[counter++] = (uint32_t)v_index;
					}
				} break;
				case 0x3: // vertex + uv
				{
					uint32_t counter = 0;
					for(uint32_t i = 0; i < f_index; ++i)
					{
						v_index = faces[i].v - 1;
						uv_index = faces[i].uv - 1;
						m_index = faces[i].m;

						obj_result->vertex_buffer[v_index].point = points[v_index];
						obj_result->vertex_buffer[v_index].texture_coordinate = texture_coordinates[uv_index];
						obj_result->vertex_buffer[v_index].material_index = m_index;

						((uint32_t*)obj_result->index_buffer)[counter++] = (uint32_t)v_index;
					}
				} break;
				default: // vertex
					uint32_t counter = 0;
					for(uint32_t i = 0; i < f_index; ++i)
					{
						v_index = faces[i].v - 1;
						m_index = faces[i].m;

						obj_result->vertex_buffer[v_index].point = points[v_index];
						obj_result->vertex_buffer[v_index].material_index = m_index;

						((uint32_t*)obj_result->index_buffer)[counter++] = (uint32_t)v_index;
					}
			}
		}
	}

	free(data);
	fclose(file);

	if(obj_info->total_materials)
	{
		fopen_s(&file, material_file, "rb");
		if (!file)
		{
			printf("Failed to open mtl file: %s\n", material_file);
			return 1;
		}

		fseek(file, 0L, SEEK_END);
		file_size = ftell(file);
		fseek(file, 0L, SEEK_SET);

		data = (char*)malloc(sizeof(char) * file_size);
		// TODO: Maybe read in chunks not all file at once.
		fread(data, file_size, 1, file);
		current_ptr = data;
		current_size = file_size;

		m_index = 0;

		while(current_size)
		{
			if(current_ptr[0] == 'n' && current_ptr[1] == 'e' && current_ptr[2] == 'w' && current_ptr[3] == 'm' && current_ptr[4] == 't' && current_ptr[5] == 'l')
			{
				++m_index;
			}
			else if(current_ptr[0] == 'K' && current_ptr[1] == 'a')
			{
				current_ptr += 3;
				current_size -= 3;
				
				obj_result->material_buffer[m_index].ambient.red = parse_float(current_ptr, &consumed_size);
				current_ptr += consumed_size + 1; // removes ' ' between values.
				current_size -= consumed_size + 1;
				obj_result->material_buffer[m_index].ambient.green = parse_float(current_ptr, &consumed_size);
				current_ptr += consumed_size + 1; // removes ' ' between values.
				current_size -= consumed_size + 1;
				obj_result->material_buffer[m_index].ambient.blue = parse_float(current_ptr, &consumed_size);
				current_ptr += consumed_size;
				current_size -= consumed_size;
			}
			else if(current_ptr[0] == 'K' && current_ptr[1] == 'd')
			{
				current_ptr += 3;
				current_size -= 3;
				
				obj_result->material_buffer[m_index].diffuse.red = parse_float(current_ptr, &consumed_size);
				current_ptr += consumed_size + 1; // removes ' ' between values.
				current_size -= consumed_size + 1;
				obj_result->material_buffer[m_index].diffuse.green = parse_float(current_ptr, &consumed_size);
				current_ptr += consumed_size + 1; // removes ' ' between values.
				current_size -= consumed_size + 1;
				obj_result->material_buffer[m_index].diffuse.blue = parse_float(current_ptr, &consumed_size);
				current_ptr += consumed_size;
				current_size -= consumed_size;
			}
			else if(current_ptr[0] == 'K' && current_ptr[1] == 's')
			{
				current_ptr += 3;
				current_size -= 3;
				
				obj_result->material_buffer[m_index].specular.red = parse_float(current_ptr, &consumed_size);
				current_ptr += consumed_size + 1; // removes ' ' between values.
				current_size -= consumed_size + 1;
				obj_result->material_buffer[m_index].specular.green = parse_float(current_ptr, &consumed_size);
				current_ptr += consumed_size + 1; // removes ' ' between values.
				current_size -= consumed_size + 1;
				obj_result->material_buffer[m_index].specular.blue = parse_float(current_ptr, &consumed_size);
				current_ptr += consumed_size;
				current_size -= consumed_size;
			}
			else if(current_ptr[0] == 'K' && current_ptr[1] == 'e')
			{
				current_ptr += 3;
				current_size -= 3;
				
				obj_result->material_buffer[m_index].emissive.red = parse_float(current_ptr, &consumed_size);
				current_ptr += consumed_size + 1; // removes ' ' between values.
				current_size -= consumed_size + 1;
				obj_result->material_buffer[m_index].emissive.green = parse_float(current_ptr, &consumed_size);
				current_ptr += consumed_size + 1; // removes ' ' between values.
				current_size -= consumed_size + 1;
				obj_result->material_buffer[m_index].emissive.blue = parse_float(current_ptr, &consumed_size);
				current_ptr += consumed_size;
				current_size -= consumed_size;
			}
			else if(current_ptr[0] == 'N' && current_ptr[1] == 's')
			{
				current_ptr += 3;
				current_size -= 3;
				
				obj_result->material_buffer[m_index].specular_exponent = parse_float(current_ptr, &consumed_size);
				current_ptr += consumed_size;
				current_size -= consumed_size;
			}
			else if(current_ptr[0] == 'd' && current_ptr[1] == ' ')
			{
				current_ptr += 2;
				current_size -= 2;
				
				obj_result->material_buffer[m_index].transparency = parse_float(current_ptr, &consumed_size);
				current_ptr += consumed_size;
				current_size -= consumed_size;
			}
			else if(current_ptr[0] == 'T' && current_ptr[1] == 'r')
			{
				current_ptr += 3;
				current_size -= 3;
				
				obj_result->material_buffer[m_index].transparency = parse_float(current_ptr, &consumed_size);
				obj_result->material_buffer[m_index].transparency = 1 - obj_result->material_buffer[m_index].transparency;
				current_ptr += consumed_size;
				current_size -= consumed_size;
			}
			else if(current_ptr[0] == 'T' && current_ptr[1] == 'f')
			{
				// NOTE: We do not support Transmission Filter Color for now.
			}
			else if(current_ptr[0] == 'N' && current_ptr[1] == 'i')
			{
				current_ptr += 3;
				current_size -= 3;
				
				obj_result->material_buffer[m_index].ior= parse_float(current_ptr, &consumed_size);
				current_ptr += consumed_size;
				current_size -= consumed_size;
			}
			else if(current_ptr[0] == 'i' && current_ptr[1] == 'l' && current_ptr[2] == 'l' && current_ptr[3] == 'u' && current_ptr[4] == 'm')
			{
				current_ptr += 6;
				current_size -= 6;

				obj_result->material_buffer[m_index].illumination_mode = parse_integer(current_ptr, &consumed_size);
				current_ptr += consumed_size;
				current_size -= consumed_size;
			}
			// NOTE: For now we do not support texture maps.
			else if(current_ptr[0] == 'm' && current_ptr[1] == 'a' && current_ptr[2] == 'p' && current_ptr[3] == '_')
			{
				// TODO: For now we dont support external texture data.
			}


			while(current_ptr[0] != '\n')
			{
				++current_ptr;
				--current_size;

				if(!current_size) break;
			}

			++current_ptr;
			--current_size;
		}

		fclose(file);
	}

	free(memory);

	return 0;
}

#endif
