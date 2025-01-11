#ifndef CORE_TEST_WORLD_H
#define CORE_TEST_WORLD_H

#include "photon.h"

void init_world(world32_t *world);

void init_world(world64_t *world);

void create_world_case_obj(world32_t *world,
        float aspect_ratio,
        uint32_t framebuffer_width,
        uint32_t framebuffer_height,
        uint32_t samples_per_pixel,
        uint32_t max_bounces,
        camera32_t *camera);

void create_world_case_obj(world64_t *world,
        float aspect_ratio,
        uint32_t framebuffer_width,
        uint32_t framebuffer_height,
        uint32_t samples_per_pixel,
        uint32_t max_bounces,
        camera64_t *camera);

#endif
