/* =============================================================================
// PROJECT CHRONO - http://projectchrono.org
//
// Copyright (c) 2019 projectchrono.org
// All rights reserved.
//
// Use of this source code is governed by a BSD-style license that can be found
// in the LICENSE file at the top level of the distribution and at
// http://projectchrono.org/license-chrono.txt.
//
// =============================================================================
// Authors: Asher Elmquist
// =============================================================================
//
// RT kernels for coloring upon ray not intersecting anything
//
// ============================================================================= */

#include <math_constants.h>
#include <optixu/optixu_aabb.h>
#include "ray_utils.h"

using namespace optix;

rtDeclareVariable(PerRayData_camera, prd_camera, rtPayload, );
rtDeclareVariable(PerRayData_lidar, prd_lidar, rtPayload, );
rtDeclareVariable(PerRayData_radar, prd_radar, rtPayload, );
rtDeclareVariable(optix::Ray, ray, rtCurrentRay, );
rtDeclareVariable(float3, default_color, , );
rtDeclareVariable(float, default_depth, , );
rtDeclareVariable(float, default_range, , );

// environment map
rtTextureSampler<float4, 2> environment_map;
rtDeclareVariable(int, has_environment_map, , );

RT_PROGRAM void camera_miss() {
    if (has_environment_map) {
        float theta = atan2f(ray.direction.x, ray.direction.y);
        float phi = asinf(ray.direction.z);
        float tex_x = theta / (2 * M_PIf);
        float tex_y = phi / (M_PIf) + 0.5;

        prd_camera.color = make_float3(tex2D(environment_map, tex_x, tex_y));

    } else {
        float theta = atan2f(ray.direction.x, ray.direction.y);
        float phi = asinf(ray.direction.z);
        float tex_x = theta / (2 * M_PIf);
        float tex_y = phi / (M_PIf) + 0.5;

        prd_camera.color = make_float3(tex2D(environment_map, tex_x, tex_y));
        // prd_camera.color.x = powf(prd_camera.color.x / 255.0f, 2.2f) * 255.0f;
        // prd_camera.color.y = powf(prd_camera.color.y / 255.0f, 2.2f) * 255.0f;
        // prd_camera.color.z = powf(prd_camera.color.z / 255.0f, 2.2f) * 255.0f;
    }
    if (prd_camera.mode == GLOBAL_ILLUMINATION) {
        if (prd_camera.depth == 1) {
            prd_camera.normal = normalize(-ray.direction);
            prd_camera.albedo = prd_camera.color;
        }
        prd_camera.color *= prd_camera.contribution_to_firsthit;
        // terminate the ray when it hit the environment map
        prd_camera.contribution_to_firsthit = make_float3(0.0f);
    }
}

RT_PROGRAM void lidar_miss() {
    prd_lidar.range = -1;
    prd_lidar.intensity = 0.f;
}

RT_PROGRAM void radar_miss(){
    prd_radar.range = -1;
    prd_radar.rcs = 0;
}
