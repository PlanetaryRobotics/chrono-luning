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
#include "chrono_thirdparty/filesystem/path.h"
#include "chrono_granular/api/ChApiGranularChrono.h"
#include "chrono_granular/physics/ChGranular.h"
#include "chrono/utils/ChUtilsSamplers.h"
#include "chrono/utils/ChUtilsInputOutput.h"
#include "chrono_granular/utils/ChCudaMathUtils.cuh"
#include "demos/granular/ChGranularDemoUtils.hpp"
#include "demos/granular/GoldenburgHelpers.hpp"
#include "chrono/core/ChStream.h"
#include "chrono/core/ChVector.h"
#include "chrono_thirdparty/filesystem/path.h"

using namespace chrono;
using namespace chrono::granular;

// Show command line usage
void ShowUsage(std::string name) {
    std::cout << "usage: " + name + " <radius (cm)> " + " <fric coeff mu> " + " <F_ext_ratio * mg> "  << std::endl;
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
    int x_dim_num = 31;
    int z_dim_num = 15;


    // box dim
    float box_X = x_dim_num * sphere_radius * 2;
    float box_Y = box_X;
    float box_Z = z_dim_num * sphere_radius * 2;

    // material based parameter
    double sphere_volume = 4.0f/3.0f * CH_C_PI * sphere_radius * sphere_radius * sphere_radius;
    double sphere_mass = sphere_density * sphere_volume;

    double gravity = 980.0f;

    // stiffness
    double kt_over_kn = 0.8f;
    double kn = 30000 * sphere_mass * gravity/sphere_radius;
    double kt = kt_over_kn * kn;
    // damping parameters
    double gamma_n = 50000;
    double gamma_t = kt_over_kn * gamma_n;

    printf("radius, %e, density, %e, mu, %e, {kn, kt, gn, gt} = {%e, %e, %e, %e}\n", sphere_radius, sphere_density, friction_coeff, kn, kt, gamma_n, gamma_t);

    // friction coefficient
    double mu_s_s2s = friction_coeff;
    double mu_s_s2w = friction_coeff;

    // set gravity
    double grav_X = 0.0f;
    double grav_Y = 0.0f;
    double grav_Z = -980.0f;
    double grav_mag = std::abs(grav_Z);

    // time integrator
    // double time_settling = 5.0f;
    // double time_Fduration = 2.0f;
    // double time_extF = 4.0f;
    // double time_end = time_settling + time_Fduration + time_extF;
    // double step_size = 5e-6;

    // time integrator for testing
    double time_settling = 0.2f;
    double time_Fduration = 1.0f;
    double time_extF = 1.0f;
    double time_end = time_settling + time_Fduration + time_extF;
    double step_size = 5e-6;



    // setup simulation gran_sys
    ChSystemGranularSMC gran_sys(sphere_radius, sphere_density,
                                 make_float3(box_X, box_Y, box_Z));

    ChGranularSMC_API apiSMC;
    apiSMC.setGranSystem(&gran_sys);

    // intialize particle positions
    std::vector<ChVector<float>> initialPos = initializePositions(x_dim_num, z_dim_num, sphere_radius);
    int numSpheres = initialPos.size();
    std::cout << "number of spheres: " << initialPos.size() << "\n";
    apiSMC.setElemsPositions(initialPos);

    float psi_T = 32.0f;
    float psi_L = 256.0f;
    gran_sys.setPsiFactors(psi_T, psi_L);

    // normal force model
    gran_sys.set_K_n_SPH2SPH(kn);
    gran_sys.set_K_n_SPH2WALL(2*kn);
    gran_sys.set_Gamma_n_SPH2SPH(gamma_n);
    gran_sys.set_Gamma_n_SPH2WALL(gamma_n);
    
    // tangential force model
    gran_sys.set_friction_mode(GRAN_FRICTION_MODE::MULTI_STEP);
    gran_sys.set_K_t_SPH2SPH(kt);
    gran_sys.set_K_t_SPH2WALL(kt);
    gran_sys.set_Gamma_t_SPH2SPH(gamma_t);
    gran_sys.set_Gamma_t_SPH2WALL(gamma_t);
    gran_sys.set_static_friction_coeff_SPH2SPH(mu_s_s2s);
    gran_sys.set_static_friction_coeff_SPH2WALL(mu_s_s2w);

    // set gravity
    gran_sys.set_gravitational_acceleration(grav_X, grav_Y, grav_Z);

    // Set the position of the BD
    gran_sys.set_BD_Fixed(true);

    // set time integrator
    gran_sys.set_timeIntegrator(GRAN_TIME_INTEGRATOR::FORWARD_EULER);
    gran_sys.set_fixed_stepSize(step_size);

    // record contact info
    gran_sys.setRecordingContactInfo(true);
    gran_sys.initialize();
    int top_center_sphereID = findTopCenterSphereID(gran_sys, numSpheres);

    // set output directory
    char out_dir[100];
    sprintf(out_dir, "31_radius_%.0ecm_mu_%.1e", sphere_radius, friction_coeff);

    // create folder for outputs
    if (!filesystem::create_directory(filesystem::path(out_dir))){
        std::cout << "Error creating directory " << out_dir << std::endl;
        return -1;
    }

    double curr_time = 0;
    int currframe = 0;
    double frame_size = 0.05;
    int print_force_per_steps = (int)(frame_size/step_size);

    // initialize values that I want to keep track of
    double sysKE, avgKE;
    std::vector<force_over_x> pos_force_array;
    double force_left, force_right, diff; // tracker to make sure force profile is symmetric

    clock_t start = std::clock();

    // PHASE ONE: SETTLING
	while (curr_time < time_settling) {
        gran_sys.advance_simulation(frame_size);
        curr_time += frame_size;

        pos_force_array = getSortedBoundaryForces(gran_sys, numSpheres, sphere_radius, 2*kn, gamma_n, sphere_mass, -box_Z/2.0f);
                
        force_left = pos_force_array.at(0).force;
        force_right = pos_force_array.at(x_dim_num-1).force;
        diff = std::abs((force_left - force_right)/force_left);

        sysKE = getSystemKE(sphere_radius, sphere_density, apiSMC, numSpheres);
        avgKE = sysKE/numSpheres/(sphere_mass*gravity*sphere_radius);

        std::cout << curr_time << ", "  << diff << ", " << avgKE << std::endl;
        currframe++;
    }

    // output reaction force at settling stage
    char settling_reactionF_filename[100];
    sprintf(settling_reactionF_filename, "%s/settling_pos_force.csv", out_dir);
    std::ofstream forcestream(std::string(settling_reactionF_filename), std::ios::out);

    for (int i = 0; i < pos_force_array.size(); i ++){
        forcestream << std::setprecision(7) << pos_force_array.at(i).pos << ", " << pos_force_array.at(i).force << "\n";
    }
    forcestream.close();
    // write position at settling stage
    char settling_pos_filename[100];
    sprintf(settling_pos_filename, "%s/settling_position", out_dir);
    gran_sys.writeFile(std::string(settling_pos_filename));
    // write contact network at settling stage
    char settling_contact_filename[100];
    sprintf(settling_contact_filename, "%s/settling_contact_forces", out_dir, currframe);
    gran_sys.writeContactInfoFile(std::string(settling_contact_filename));

    // PHASE TWO: APPLYING FORCE GRADUALLY
    int F_ext_ratio_array_size = std::round(time_Fduration/step_size)+1;
    double F_ext_ratio_array[F_ext_ratio_array_size];
    double slope = F_ext_ratio/time_Fduration;
    for (int i = 0; i < F_ext_ratio_array_size; i++){
        F_ext_ratio_array[i] = slope * step_size * i;
    }
    int force_counter = 0;
    while (curr_time < time_settling + time_Fduration && force_counter < F_ext_ratio_array_size){
        gran_sys.advance_simulation(step_size);
        curr_time += step_size;

        gran_sys.setWavePropagationParameters(top_center_sphereID, F_ext_ratio_array[force_counter], grav_mag);
        force_counter++;

        if (force_counter%print_force_per_steps == 0){

            pos_force_array = getSortedBoundaryForces(gran_sys, numSpheres, sphere_radius, 2*kn, gamma_n, sphere_mass, -box_Z/2.0f);

            force_left = pos_force_array.at(0).force;
            force_right = pos_force_array.at(x_dim_num-1).force;
            diff = std::abs((force_left - force_right)/force_left);

            sysKE = getSystemKE(sphere_radius, sphere_density, apiSMC, numSpheres);
            avgKE = sysKE/numSpheres/(sphere_mass*gravity*sphere_radius);
            std::cout << curr_time << ", "  << diff << ", " << avgKE << ", " <<  F_ext_ratio_array[force_counter-1] << std::endl;
        }

    }

    // PHASE THREE: KEEPING EXTERNAL FORCE CONSTANT
    gran_sys.setWavePropagationParameters(top_center_sphereID, F_ext_ratio, grav_mag);
	while (curr_time < time_end) {

        gran_sys.advance_simulation(frame_size);
        curr_time += frame_size;

        sysKE = getSystemKE(sphere_radius, sphere_density, apiSMC, numSpheres);
        avgKE = sysKE/numSpheres/(sphere_mass*gravity*sphere_radius);

        pos_force_array = getSortedBoundaryForces(gran_sys, numSpheres, sphere_radius, 2*kn, gamma_n, sphere_mass, -box_Z/2.0f);

        force_left = pos_force_array.at(0).force;
        force_right = pos_force_array.at(x_dim_num-1).force;
        diff = std::abs((force_left - force_right)/force_left);

        std::cout << curr_time << ", "  << diff << ", " << avgKE << std::endl;
        currframe++;

    }

    char outForceFile_string[100];
    sprintf(outForceFile_string, "%s/ending_pos_force_F_%.1fmg.csv", out_dir, F_ext_ratio);
    std::ofstream outstream(std::string(outForceFile_string), std::ios::out);

    for (int i = 0; i < pos_force_array.size(); i ++){
        outstream << std::setprecision(7) << pos_force_array.at(i).pos << ", " << pos_force_array.at(i).force << "\n";
    }
    outstream.close();


    char end_pos_filename[100];
    sprintf(end_pos_filename, "%s/ending_position_F_%.1fmg", out_dir, F_ext_ratio);
    gran_sys.writeFile(std::string(end_pos_filename));

    char end_contact_filename[100];
    sprintf(end_contact_filename, "%s/ending_contact_forces_network_F_%.1fmg", out_dir, F_ext_ratio);
    gran_sys.writeContactInfoFile(std::string(end_contact_filename));

	clock_t end_time = std::clock();
	double computation_time = ((double)(end_time - start)) / CLOCKS_PER_SEC;
	std::cout << "Time: " << computation_time << " seconds" << std::endl;
    return 0;
}

