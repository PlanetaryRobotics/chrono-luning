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

struct force_over_x{
    double pos;
    double force;
};

// Show command line usage
void ShowUsage(std::string name) {
    std::cout << "usage: " + name + " <radius (cm)>  "  + " <friction coefficient mu> " << std::endl;
}

int main(int argc, char* argv[]) {


    if (argc != 3) {
        ShowUsage(argv[0]);
        return 1;
    }

    // sphere input parameters
    double sphere_radius = std::stof(argv[1]);
    double friction_coeff = std::stof(argv[2]);
    double sphere_density = 7.8;
    double step_size = 5e-6;

    // 61 spheres per layer, total of 15 layers
    int x_dim_num = 122;
    int z_dim_num = 30;

    // box dim
    float box_X = x_dim_num * sphere_radius;
    float box_Y = box_X;
    float box_Z = z_dim_num * sphere_radius;

    // material based parameter
    double sphere_volume = 4.0f/3.0f * CH_C_PI * sphere_radius * sphere_radius * sphere_radius;
    double sphere_mass = sphere_density * sphere_volume;

    double gravity = 980.0f;
//    double kt_over_kn = 0.8; used this before, with mu = 0.2, sort of worked
    double kt_over_kn = 2.0f/7.0f;
    double kn = 30000 * sphere_mass * gravity/sphere_radius;
    double kt = kt_over_kn * kn;
    // damping parameters
    double gamma_n = 50000;
    float gamma_t = kt_over_kn * gamma_n;

    printf("radius, %e, density, %e, mu, %e, {kn, kt, gn, gt} = {%e, %e, %e, %e}\n", sphere_radius, sphere_density, friction_coeff, kn, kt, gamma_n, gamma_t);

    // friction coefficient
    float mu_s_s2s = friction_coeff;
    float mu_s_s2w = friction_coeff;

    // set gravity
    float grav_X = 0.0f;
    float grav_Y = 0.0f;
    float grav_Z = -980.0f;

    // time integrator
    float time_end = 5.0f;

    // setup simulation gran_sys
    ChSystemGranularSMC gran_sys(sphere_radius, sphere_density,
                                 make_float3(box_X, box_Y, box_Z));

    ChGranularSMC_API apiSMC;
    apiSMC.setGranSystem(&gran_sys);

    std::vector<ChVector<float>> initialPos = initializePositions(x_dim_num, z_dim_num, sphere_radius);

    int numSpheres = initialPos.size();

    std::cout << "number of spheres: " << initialPos.size();

    apiSMC.setElemsPositions(initialPos);

    float psi_T = 32.0f;
    float psi_L = 256.0f;
    gran_sys.setPsiFactors(psi_T, psi_L);


    // normal force model
    gran_sys.set_K_n_SPH2SPH(kn);
    gran_sys.set_K_n_SPH2WALL(2*kn);
    gran_sys.set_Gamma_n_SPH2SPH(gamma_n);
    gran_sys.set_Gamma_n_SPH2WALL(gamma_n);
    
    //tangential force model
    gran_sys.set_friction_mode(GRAN_FRICTION_MODE::MULTI_STEP);
    gran_sys.set_K_t_SPH2SPH(kt);
    gran_sys.set_K_t_SPH2WALL(kt);
    gran_sys.set_Gamma_t_SPH2SPH(gamma_t);
    gran_sys.set_Gamma_t_SPH2WALL(gamma_t);
    gran_sys.set_static_friction_coeff_SPH2SPH(mu_s_s2s);
    gran_sys.set_static_friction_coeff_SPH2WALL(mu_s_s2w);

    gran_sys.set_gravitational_acceleration(grav_X, grav_Y, grav_Z);

    // Set the position of the BD
    gran_sys.set_BD_Fixed(true);
    gran_sys.set_timeIntegrator(GRAN_TIME_INTEGRATOR::FORWARD_EULER);
    gran_sys.set_fixed_stepSize(step_size);

    gran_sys.initialize();


    char out_dir[100];
    sprintf(out_dir, "settling_radius_%.0ecm_mu_%.1e", sphere_radius, friction_coeff);

    //create folder for outputs
    if (!filesystem::create_directory(filesystem::path(out_dir))){
        std::cout << "Error creating directory " << out_dir << std::endl;
        return -1;
    }


    float curr_time = 0;
    int currframe = 0;
    double frame_size = 0.01;
    int print_force_per_steps = (int)(frame_size/step_size);

    // initialize values that I want to keep track of
    float sysKE;
    clock_t start = std::clock();

    //gran_sys.writeFile("initialPosition");

	while (curr_time < time_end) {
        gran_sys.advance_simulation(frame_size);
        curr_time += frame_size;

        std::vector<ChVector<double>> normalForces;
        std::vector<int> particlesInContact;
        calculateBoundaryForces(gran_sys, numSpheres, sphere_radius, 2*kn, gamma_n, sphere_mass, -box_Z/2.0f, normalForces, particlesInContact);
        
        if (normalForces.size() != 0) {
            char outForceFile_string[100];
            sprintf(outForceFile_string, "%s/pos_force%06d.csv", out_dir, currframe);
            std::ofstream outstream(std::string(outForceFile_string), std::ios::out);
        
            std::vector<force_over_x> pos_force_array;
            force_over_x myData;
            for (int i = 0; i < normalForces.size(); i ++){
                myData.pos = gran_sys.getPosition(particlesInContact.at(i)).x;
                myData.force = normalForces.at(i).z();
                pos_force_array.push_back(myData);
            }

            // sort data
            std::sort(pos_force_array.begin(), pos_force_array.end(), [](auto const &a, auto const &b) {return a.pos < b.pos;});
            for (int i = 0; i < pos_force_array.size(); i ++){
                outstream << std::setprecision(7) << pos_force_array.at(i).pos << ", " << pos_force_array.at(i).force << "\n";
            }
            outstream.close();
                
            double force_left = pos_force_array.at(0).force;
            double force_right = pos_force_array.at(60).force;
            double diff = std::abs((force_left - force_right)/force_left);

            sysKE = getSystemKE(sphere_radius, sphere_density, apiSMC, numSpheres);
            double avgKE = sysKE/numSpheres/(sphere_mass*gravity*sphere_radius);

            std::cout << curr_time << ", write file " << outForceFile_string << ", diff in force: " << diff << ", avgKE " << avgKE << std::endl;


        }
        currframe++;
    }

    char settling_pos_filename[100];
    sprintf(settling_pos_filename, "%s/settling_position.csv", out_dir);
    gran_sys.writeFile(std::string(settling_pos_filename));

	clock_t end_time = std::clock();
	double computation_time = ((double)(end_time - start)) / CLOCKS_PER_SEC;
	std::cout << "Time: " << computation_time << " seconds" << std::endl;
    return 0;
}

