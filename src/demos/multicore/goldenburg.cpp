// =============================================================================
// PROJECT CHRONO - http://projectchrono.org
//
// Copyright (c) 2014 projectchrono.org
// All rights reserved.
//
// Use of this source code is governed by a BSD-style license that can be found
// in the LICENSE file at the top level of the distribution and at
// http://projectchrono.org/license-chrono.txt.
//
// =============================================================================
// Authors: Luning Fang
// =============================================================================
//
// Chrono::Multicore test program using SMC method for frictional contact.
//
// The model simulated here consists of a number of spherical objects falling
// in a fixed container.
//
// The global reference frame has Z up.
//
// If available, OpenGL is used for run-time rendering. Otherwise, the
// simulation is carried out for a pre-defined duration and output files are
// generated for post-processing with POV-Ray.
// =============================================================================


#include <cstdio>
#include <vector>
#include <cmath>
#include <stdio.h>
#include <ctime>
#include <sys/time.h>


#include "chrono_multicore/physics/ChSystemMulticore.h"
#include "chrono_thirdparty/filesystem/path.h"

#include "chrono/core/ChTimer.h"
#include "chrono/ChConfig.h"
#include "chrono/utils/ChUtilsCreators.h"
#include "chrono/utils/ChUtilsInputOutput.h"

#include "chrono/utils/ChUtilsInputOutput.h"


#include "chrono_opengl/ChVisualSystemOpenGL.h"
#include "demos/multicore/goldenburg_helpers.cpp"
using namespace chrono;
using namespace chrono::collision;



// Show command line usage
void ShowUsage(std::string name) {
    std::cout << "usage: " + name + " <fric coeff mu> "  + " <stiffness ratio > 3e4>"  + " <step size> " +  " <test_folder> " << std::endl;
}



// -----------------------------------------------------------------------------
// Create the system, specify simulation parameters, and run simulation loop.
// -----------------------------------------------------------------------------
int main(int argc, char* argv[]) {
    GetLog() << "Copyright (c) 2017 projectchrono.org\nChrono version: " << CHRONO_VERSION << "\n\n";


    if (argc != 5) {
        ShowUsage(argv[0]);
        return 1;
    }    


    // input parameters
    double mu = std::stof(argv[1]);
    kn_ratio = std::stof(argv[2]);
    time_step = std::stof(argv[3]);

    char test_dir[300];
    sprintf(test_dir, "%s_kn_%.0e_dt_%.0e", argv[4], kn_ratio, time_step);


    if (!filesystem::create_directory(filesystem::path(test_dir))) {
        std::cout << "Error creating directory " << test_dir << std::endl;
        return -1;
    }
    std::cout << "Create test folder: " << test_dir << std::endl;
    std::cout << "Running tests of kn = " << kn_ratio << "mg/R, and step size of " << time_step << std::endl;

    // set subfolder directory, mu=0, 0.1 and 0.5
    char subtest_dir[400];
    sprintf(subtest_dir, "%s/mu_%.1e", test_dir, mu);

    // create folder for outputs
    if (!filesystem::create_subdirectory(filesystem::path(subtest_dir))) {
        std::cout << "Error creating sub directory " << subtest_dir << std::endl;
        return -1;
    }
    std::cout << "Create test subfolder: " << subtest_dir << std::endl;


    // Simulation parameters for debugging 
    // ---------------------
    double time_settle = 0.5;
    double time_F_dur = 0.03;
    double time_end = 1;
    // Simulation parameters for actual tests
    // double time_settle = 0.5;
    // double time_F_dur = 0.03;
    // double time_end = 1;


    double output_rate = 0.01; // write output contact force every 0.02 seconds
    int output_per_step = (int)(output_rate/time_step);

    // variables for timing
    struct timeval start;
    struct timeval end;
    
    // Create system
    // -------------
    ChSystemMulticoreSMC msystem;

    // Set number of threads
    msystem.SetNumThreads(8);

    // Set gravitational acceleration
    msystem.Set_G_acc(ChVector<>(0, 0, -gravity));

    // Set solver parameters
    msystem.GetSettings()->collision.narrowphase_algorithm = ChNarrowphase::Algorithm::MPR;
    msystem.GetSettings()->collision.bins_per_axis = vec3(4, 4, 4);
    msystem.GetSettings()->solver.contact_force_model = ChSystemSMC::ContactForceModel::Hooke;
    msystem.GetSettings()->solver.tangential_displ_mode = ChSystemSMC::TangentialDisplacementModel::MultiStep;
    msystem.GetSettings()->solver.use_material_properties = false;

    // Create the fixed and moving bodies
    // ----------------------------------
    AddContainer(&msystem, sphere_radius, mu);
    AddFallingBalls(&msystem, sphere_radius, mu);

    // Perform the simulation

    // Run simulation for specified time
    int num_steps = (int)std::ceil(time_end / time_step);
    double curr_time = 0;
    auto creporter = chrono_types::make_shared<ContactReporter>();

    bool useVis = false;

    #ifdef CHRONO_OPENGL
        opengl::ChVisualSystemOpenGL vis;
    #endif

    #ifdef CHRONO_OPENGL
        if (useVis == true){
            vis.AttachSystem(&msystem);
            vis.SetWindowTitle("Balls SMC");
            vis.SetWindowSize(1280, 720);
            vis.SetRenderMode(opengl::WIREFRAME);
            vis.Initialize();
            vis.SetCameraPosition(ChVector<>(0, -10, 2), ChVector<>(0, 0, 0));
            vis.SetCameraVertical(CameraVerticalDir::Z);
    
        }
    #endif


    int curr_step = 0;
    gettimeofday(&start, NULL);

    for (curr_step = 0; curr_step < num_steps; curr_step++) {
        msystem.DoStepDynamics(time_step);


        #ifdef CHRONO_OPENGL
            if (useVis == true){
                vis.Render();
            }
        #endif

        curr_time += time_step;

        if (curr_step%output_per_step == 0){
            gettimeofday(&end, NULL);
            double KE_ratio = calcKE(&msystem)/(gravity * sphere_mass * sphere_radius );
            std::cout << "t = " << curr_time << ", KE/mgR: " <<  KE_ratio << ", simulation of " << output_rate << "sec took " <<  time_diff(&start, &end) << " cpu seconds. " << std::endl;
            gettimeofday(&start, NULL);

            // print out contact force and write particle positions
            WrtieOutputInfo(&msystem, subtest_dir, int(curr_step/output_per_step));

            if (KE_ratio < 1e-7){
                break;
            }            
        }
        if (curr_time > time_settle){
            break;
        }
    }
    

    // ************************************
    // PHASE TWO: APPLYING FORCE GRADUALLY
    // ************************************
    int F_ext_ratio_array_size = std::round(time_F_dur / time_step) + 1;
    double *F_ext_ratio_array = (double*) malloc(F_ext_ratio_array_size * sizeof(double));
    double slope = F_ext_ratio / time_F_dur;

    for (int i = 0; i < F_ext_ratio_array_size; i++) {
        F_ext_ratio_array[i] = slope * time_step * i;
    }

    int force_counter = 0;
    while (curr_time < time_settle + time_F_dur && force_counter < F_ext_ratio_array_size) {
        msystem.DoStepDynamics(time_step);
        curr_step++;

        #ifdef CHRONO_OPENGL
            if (useVis == true){
                vis.Render();
            }
        #endif


        curr_time += time_step;

        double external_force = -F_ext_ratio_array[force_counter] * top_middle->GetMass() * gravity;
        // apply gradual force
        top_middle->Empty_forces_accumulators();
        top_middle->Accumulate_force(ChVector<double>(0, 0, external_force), top_middle->GetPos(), false);

        force_counter++;

        if (curr_step % output_per_step == 0){
            std::cout << "t = " << curr_time << ", applied force = " << external_force <<  ", KE/mgR: " << calcKE(&msystem)/(gravity * sphere_mass * sphere_radius )<< std::endl;
            // print out contact force and write particle positions
            WrtieOutputInfo(&msystem, subtest_dir, int(curr_step/output_per_step));
        }

    }

    // free the memory because i'm a good samaritan 
    free(F_ext_ratio_array);

    while(curr_time < time_end){

        // apply constant force
        top_middle->Empty_forces_accumulators();
        top_middle->Accumulate_force(ChVector<double>(0, 0, -F_ext_ratio * top_middle->GetMass() * gravity), top_middle->GetPos(), false);
        msystem.DoStepDynamics(time_step);

        #ifdef CHRONO_OPENGL
            if (useVis == true){
                vis.Render();
            }           
        #endif

        curr_time += time_step;


        if (curr_step % output_per_step == 0){
            std::cout << "t = " << curr_time << ", KE/mgR: " << calcKE(&msystem)/(gravity * sphere_mass * sphere_radius )<< std::endl;
            // print out contact force and write particle positions
            WrtieOutputInfo(&msystem, subtest_dir, int(curr_step/output_per_step));

        }
        curr_step++;
    }



    // #endif

    return 0;
}