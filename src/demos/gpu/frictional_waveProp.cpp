// =============================================================================
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
// Authors: Luning Fang
// =============================================================================
// settling phase of goldenburg test
// =============================================================================

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <ctime>
#include <iomanip>
#include <algorithm>
#include <stdlib.h>

#include "chrono_thirdparty/filesystem/path.h"
#include "chrono_gpu/physics/ChSystemGpu.h"
#include "chrono/core/ChGlobal.h"
#include "chrono/core/ChVector.h"
#include "chrono/utils/ChUtilsSamplers.h"
#include "chrono/utils/ChUtilsInputOutput.h"

#include "demos/gpu/GoldenburgHelpers.cpp"

using namespace chrono;
using namespace chrono::gpu;

bool verbose = true;

// Show command line usage
void ShowUsage(std::string name) {
    std::cout << "usage: " + name + " <radius (cm)> " + " <fric coeff mu> " + " <F_ext_ratio * mg> " << std::endl;
}

int main(int argc, char* argv[]) {
    if (argc != 4) {
        ShowUsage(argv[0]);
        return 1;
    }

    // sphere input parameters
    double sphere_radius = std::stof(argv[1]);
    double friction_coeff = std::stof(argv[2]);
    double F_ext_ratio = std::stof(argv[3]);

    double sphere_density = 7.8;

    // 31 spheres per layer, total of 15 layers
    int x_dim_num = 61;
    int z_dim_num = 15;

    // box dim
    float box_X = x_dim_num * sphere_radius * 2;
    float box_Y = box_X;
    float box_Z = z_dim_num * sphere_radius * 2;

    // material based parameter
    double sphere_volume = 4.0f / 3.0f * CH_C_PI * sphere_radius * sphere_radius * sphere_radius;
    double sphere_mass = sphere_density * sphere_volume;

    double gravity = 1000.0f;

    // stiffness
    double kt_over_kn = 1.0f;
    double kn = 30000 * sphere_mass * gravity / sphere_radius;
    double kt = kt_over_kn * kn;
    // damping parameters
    double gamma_n = 5000;
    double gamma_t = 5000;

    printf("radius, %e, density, %e, mu, %e, {kn, kt, gn, gt} = {%e, %e, %e, %e}\n", sphere_radius, sphere_density,
           friction_coeff, kn, kt, gamma_n, gamma_t);

    // friction coefficient
    double mu_s_s2s = friction_coeff;
    double mu_s_s2w = friction_coeff;

    // set gravity
    double grav_X = 0.0f;
    double grav_Y = 0.0f;
    double grav_Z = -1000.0f;
    double grav_mag = std::abs(grav_Z);


    // time integrator for Milad
    // double time_settling = 0.5f;
    // double time_Fduration = 0.5f;
    // double time_extF = 0.5f;
    // double step_size = 5e-7;

    double time_settling = 0.5f;
    double time_Fduration = 0.5f;
    double time_extF = 0.5f;
    double step_size = 1e-7;

    double time_end = time_settling + time_Fduration + time_extF;

    // setup simulation gran_sys
    ChSystemGpu gran_sys(sphere_radius, sphere_density, ChVector<float>(box_X, box_Y, box_Z));

    // intialize particle positions
    std::vector<ChVector<float>> initialPos = initializePositions(x_dim_num, z_dim_num, sphere_radius);
    int numSpheres = initialPos.size();
    std::cout << "number of spheres: " << initialPos.size() << "\n";
    gran_sys.SetParticles(initialPos);

    // normal force model
    gran_sys.SetKn_SPH2SPH(kn);
    gran_sys.SetKn_SPH2WALL(kn);
    gran_sys.SetGn_SPH2SPH(gamma_n);
    gran_sys.SetGn_SPH2WALL(gamma_n);

    // tangential force model
    gran_sys.SetFrictionMode(CHGPU_FRICTION_MODE::MULTI_STEP);
    gran_sys.SetKt_SPH2SPH(kt);
    gran_sys.SetKt_SPH2WALL(kt);
    gran_sys.SetGt_SPH2SPH(gamma_t);
    gran_sys.SetGt_SPH2WALL(gamma_t);
    gran_sys.SetStaticFrictionCoeff_SPH2SPH(mu_s_s2s);
    gran_sys.SetStaticFrictionCoeff_SPH2WALL(mu_s_s2w);

    // set gravity
    gran_sys.SetGravitationalAcceleration(ChVector<float>(grav_X, grav_Y, grav_Z));

    // Set the position of the BD
    gran_sys.SetBDFixed(true);

    // set time integrator
    gran_sys.SetTimeIntegrator(CHGPU_TIME_INTEGRATOR::CHUNG);
    gran_sys.SetFixedStepSize(step_size);

    // record contact info
    gran_sys.Initialize();
    int top_center_sphereID = findTopCenterSphereID(gran_sys, numSpheres);

    // set output directory
    char out_dir[100];
    sprintf(out_dir, "61_radius_%.0ecm_mu_%.1e", sphere_radius, friction_coeff);

    // create folder for outputs
    if (!filesystem::create_directory(filesystem::path(out_dir))) {
        std::cout << "Error creating directory " << out_dir << std::endl;
        return -1;
    }

    double curr_time = 0;
    int currframe = 0;
    double frame_size = 0.05;

    int print_force_per_steps = (int)(frame_size / step_size);

    // initialize values that I want to keep track of
    double sysKE, avgKE;
    std::vector<force_over_x> pos_force_array;
    double force_left, force_right, diff;  // tracker to make sure force profile is symmetric

    clock_t start = std::clock();


    char output_pos_filename[500];

    // PHASE ONE: SETTLING
    while (curr_time < time_settling) {
        gran_sys.AdvanceSimulation(frame_size);
        curr_time += frame_size;

        pos_force_array =
            getSortedBoundaryForces(gran_sys, numSpheres, sphere_radius, kn, gamma_n, sphere_mass, -box_Z / 2.0f);

        force_left = pos_force_array.at(0).force;
        force_right = pos_force_array.at(x_dim_num - 1).force;
        diff = std::abs((force_left - force_right) / force_left);

        sysKE = getSystemKE(sphere_radius, sphere_density, gran_sys, numSpheres);
        avgKE = sysKE / numSpheres / (sphere_mass * gravity * sphere_radius);

        std::cout << curr_time << ", " << diff << ", " << avgKE << ", " << gran_sys.GetParticlePosition(452).x() <<  std::endl;
        currframe++;
        sprintf(output_pos_filename, "%s/output_pos_step%05d.csv", out_dir, currframe);
        gran_sys.WriteParticleFile(std::string(output_pos_filename));

    }

    // output reaction force at settling stage
    char settling_reactionF_filename[500];
    sprintf(settling_reactionF_filename, "%s/settling_pos_force.csv", out_dir);
    std::ofstream forcestream(std::string(settling_reactionF_filename), std::ios::out);

    for (int i = 0; i < pos_force_array.size(); i++) {
        forcestream << std::setprecision(7) << pos_force_array.at(i).pos << ", " << pos_force_array.at(i).force << "\n";
    }
    forcestream.close();
    // write position at settling stage
    char settling_pos_filename[500];
    sprintf(settling_pos_filename, "%s/settling_position.csv", out_dir);
    gran_sys.WriteParticleFile(std::string(settling_pos_filename));

    if (verbose){
        std::cout << "finish PHASE ONE, settling phase" << std::endl;
    }


    // ************************************
    // PHASE TWO: APPLYING FORCE GRADUALLY
    // ************************************
    int F_ext_ratio_array_size = std::round(time_Fduration / step_size) + 1;

    double *F_ext_ratio_array = (double*) malloc(F_ext_ratio_array_size * sizeof(double));

    double slope = F_ext_ratio / time_Fduration;

    for (int i = 0; i < F_ext_ratio_array_size; i++) {
        F_ext_ratio_array[i] = slope * step_size * i;

    }
    int force_counter = 0;
    while (curr_time < time_settling + time_Fduration && force_counter < F_ext_ratio_array_size) {
        gran_sys.AdvanceSimulation(step_size);
        curr_time += step_size;

        gran_sys.setWavePropagationParameters(top_center_sphereID, F_ext_ratio_array[force_counter], grav_mag);
        force_counter++;

        if (force_counter % print_force_per_steps == 0) {
            pos_force_array = getSortedBoundaryForces(gran_sys, numSpheres, sphere_radius, kn, gamma_n, sphere_mass,
                                                      -box_Z / 2.0f);

            force_left = pos_force_array.at(0).force;
            force_right = pos_force_array.at(x_dim_num - 1).force;
            diff = std::abs((force_left - force_right) / force_left);

            sysKE = getSystemKE(sphere_radius, sphere_density, gran_sys, numSpheres);
            avgKE = sysKE / numSpheres / (sphere_mass * gravity * sphere_radius);
            std::cout << curr_time << ", " << diff << ", " << avgKE << ", " << F_ext_ratio_array[force_counter - 1] << ", " << gran_sys.GetParticlePosition(452).x()
                      << std::endl;
        }
    }
    free(F_ext_ratio_array);
    // PHASE THREE: KEEPING EXTERNAL FORCE CONSTANT
    gran_sys.setWavePropagationParameters(top_center_sphereID, F_ext_ratio, grav_mag);
    while (curr_time < time_end) {
        gran_sys.AdvanceSimulation(frame_size);
        curr_time += frame_size;

        sysKE = getSystemKE(sphere_radius, sphere_density, gran_sys, numSpheres);
        avgKE = sysKE / numSpheres / (sphere_mass * gravity * sphere_radius);

        pos_force_array =
            getSortedBoundaryForces(gran_sys, numSpheres, sphere_radius, kn, gamma_n, sphere_mass, -box_Z / 2.0f);

        force_left = pos_force_array.at(0).force;
        force_right = pos_force_array.at(x_dim_num - 1).force;
        diff = std::abs((force_left - force_right) / force_left);

        std::cout << curr_time << ", " << diff << ", " << avgKE << ", " << gran_sys.GetParticlePosition(452).x() <<  std::endl;

        sprintf(output_pos_filename, "%s/output_pos_step%05d.csv", out_dir, currframe);
        gran_sys.WriteParticleFile(std::string(output_pos_filename));

        currframe++;
    }

    char outForceFile_string[500];
    sprintf(outForceFile_string, "%s/ending_pos_force_F_%.1fmg.csv", out_dir, F_ext_ratio);
    std::ofstream outstream(std::string(outForceFile_string), std::ios::out);

    for (int i = 0; i < pos_force_array.size(); i++) {
        outstream << std::setprecision(7) << pos_force_array.at(i).pos << ", " << pos_force_array.at(i).force << "\n";
    }
    outstream.close();

    char end_pos_filename[500];
    sprintf(end_pos_filename, "%s/ending_position_F_%.1fmg.csv", out_dir, F_ext_ratio);
    gran_sys.WriteParticleFile(std::string(end_pos_filename));

    clock_t end_time = std::clock();
    double computation_time = ((double)(end_time - start)) / CLOCKS_PER_SEC;
    std::cout << "Time: " << computation_time << " seconds" << std::endl;
    return 0;
}
