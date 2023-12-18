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
// Single wheel test for moon ranger rover
// Units in cgs 
// =============================================================================

#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <cmath>
#include <fstream>

#include "chrono/core/ChGlobal.h"
#include "chrono/physics/ChSystemSMC.h"
#include "chrono/physics/ChBody.h"
#include "chrono/physics/ChForce.h"
#include "chrono/timestepper/ChTimestepper.h"
#include "chrono/utils/ChUtilsSamplers.h"
#include "chrono/utils/ChUtilsCreators.h"
#include "chrono/assets/ChSphereShape.h"

#include "chrono/ChConfig.h"
#include "chrono/physics/ChSystemSMC.h"
#include "chrono/physics/ChBody.h"
#include "chrono/physics/ChInertiaUtils.h"
#include "chrono/physics/ChLinkMotorRotationAngle.h"
#include "chrono/utils/ChUtilsGeometry.h"
#include "chrono/utils/ChUtilsInputOutput.h"
#include "chrono/assets/ChTriangleMeshShape.h"
#include "chrono/core/ChTimer.h"

#include "chrono_gpu/physics/ChSystemGpu.h"
#include "chrono_gpu/utils/ChGpuJsonParser.h"
#include "chrono_gpu/utils/ChGpuVisualization.h"

#include "chrono_thirdparty/filesystem/path.h"

using namespace chrono;
using namespace chrono::gpu;
using namespace chrono::geometry;

// Output frequency in frames per second
float out_fps = 50;

double wheel_AngVel = 0.523599;  // prescribed angular velocity

// radius of the wheel (in cm)
double wheel_radius = 9;

ChVector<> wheel_IniVel(0.0, 0.0, 0.0);


int main(int argc, char* argv[]) {


    //double total_cart_load = 4235.5;  // one fourth the mass of moon ranger? unit in gram

     double total_cart_load = 2300.f/4.f;  // one fourth the mass of IRIS, unit in gram

    double slope = 10;
    double slope_rad = slope * CH_C_DEG_TO_RAD;
    double grav = 980;
    ChVector<float> gravity(-sin(slope_rad) * grav, 0, -cos(slope_rad) * grav);

    // directory for wheel obj
    //std::string wheel_obj = "C:/Users/fang/Documents/NSF_Collaboration/CMU_MoonRanger/data/wheel.obj";
    std::string wheel_obj = "C:/Users/fang/Documents/NSF_Collaboration/CMU_MoonRanger/data/SIMPLIFIED_IRIS_WHEEL_V2_blender.obj";
    //std::string wheel_obj = "C:/Users/fang/Downloads/luning_wheel.obj";

    // input json file
    std::string inputJson = "C:/Users/fang/Documents/NSF_Collaboration/CMU_MoonRanger/data/singleWheel.json";

    // set up parameters
    ChGpuSimulationParameters params;
    if (!ParseJSON(inputJson, params)) {
        std::cout << "ERROR: reading input file " << inputJson << std::endl;
        return 1;
    }

    // Output directory and subfolders
    std::string out_dir = "C:/Users/fang/Documents/NSF_Collaboration/CMU_MoonRanger/outputs/";
    filesystem::create_directory(filesystem::path(out_dir));
    out_dir = out_dir + params.output_dir;
    filesystem::create_directory(filesystem::path(out_dir));

    // Create particle system
    ChSystemGpuMesh gpu_sys(params.sphere_radius, params.sphere_density,
                            ChVector<float>(params.box_X, params.box_Y, params.box_Z));

    float integration_step = params.step_size;


    // Sample particle positions
    chrono::utils::PDSampler<float> sampler(2.01f * params.sphere_radius);
    double fill_bottom = -params.box_Z / 2.0f + 2 * params.sphere_radius;
    double fill_top = fill_bottom + 15;

    // leave margin at edges of sampling
    ChVector<> hdims(params.box_X / 2 - params.sphere_radius, params.box_Y / 2 - params.sphere_radius, 0);
    ChVector<> center(0, 0, fill_bottom + 2.0 * params.sphere_radius);
    std::vector<ChVector<float>> body_points;

    // create random particles for each layer 
    while (center.z() < fill_top) {
        // You can uncomment this line to see a report on particle creation process.
        std::cout << "Create layer at " << center.z() << std::endl;
        auto points = sampler.SampleBox(center, hdims);
        body_points.insert(body_points.end(), points.begin(), points.end());
        center.z() += 2.02 * params.sphere_radius;
    }
    std::cout << body_points.size() << " particles sampled." << std::endl;

    // Create MBS single wheel system first
    ChSystemSMC sysMBS;
    // Set gravity
    sysMBS.Set_G_acc(gravity);

    /* ground body */
    auto ground = chrono_types::make_shared<ChBody>();
    ground->SetIdentifier(-1);
    ground->SetPos(VNULL);
    ground->SetBodyFixed(true);
    ground->SetCollide(false);
    sysMBS.AddBody(ground);

    /* wheel body */
    // load wheel mesh (this is for computing mass and inertia only)
    auto trimesh = chrono_types::make_shared<ChTriangleMeshConnected>();
    //double scale_ratio = 100.0;  // convert m to cm
    //double scale_ratio = 0.1f;  // convert mm to cm
     double scale_ratio = 0.1f;  // convert mm to cm

    trimesh->LoadWavefrontMesh(wheel_obj, true, true);
    trimesh->Transform(ChVector<>(0, 0, 0), ChMatrix33<>(scale_ratio));  // scale to a different size
    trimesh->RepairDuplicateVertexes(1e-9);                              // if meshes are not watertight

    // Compute mass inertia from mesh
    double mmass;
    double mdensity = 2.72; // wheel density?
    ChVector<> mcog;
    ChMatrix33<> minertia;
    trimesh->ComputeMassProperties(true, mmass, mcog, minertia);
    ChMatrix33<> principal_inertia_rot;
    ChVector<> principal_I;
    ChInertiaUtils::PrincipalInertia(minertia, principal_I, principal_inertia_rot);

    ChVector<> wheel_IniPos(-params.box_X / 2 + wheel_radius * 1.2, 0.0, fill_top + wheel_radius);

    // Set the abs orientation, position and velocity
    auto wheel = chrono_types::make_shared<ChBodyAuxRef>();
    ChQuaternion<> wheel_Rot = Q_from_Euler123(ChVector<double>(0, 0, 0));

    // Set the COG coordinates to barycenter, without displacing the REF reference.
    wheel->SetFrame_COG_to_REF(ChFrame<>(mcog, wheel_Rot));

    // Set inertia
    //double wheel_mass = mmass * mdensity;
    double wheel_mass = 17.f;
    wheel->SetMass(wheel_mass);
    std::cout << "wheel mass: " << wheel_mass << " g" << std::endl;
    wheel->SetInertiaXX(mdensity * principal_I);
    wheel->SetPos_dt(wheel_IniVel); // initial wheel velocity
    wheel->SetWvel_loc(ChVector<>(0.0, 0.0, 0.0));  // initial angular velocity

    // Set the absolute position of the body:
    wheel->SetFrame_REF_to_abs(ChFrame<>(wheel_IniPos, QUNIT));
    sysMBS.AddBody(wheel);
    wheel->SetBodyFixed(true); // initially wheel fixed until particles are settled
    wheel->SetCollide(false);

    // rotational motor for the wheel
    auto motor = chrono_types::make_shared<ChLinkMotorRotationAngle>();

    motor->SetName("wheel_motor");
    // free to move in x and z direction, rotation is locked
    // z direction is the motorized direction
    motor->SetSpindleConstraint(false, false, true, true, true);
    motor->SetAngleFunction(chrono_types::make_shared<ChFunction_Ramp>(0, wheel_AngVel));

    motor->Initialize(wheel, ground, ChFrame<>(wheel->GetPos(), Q_from_AngX(-CH_C_PI_2)));
    sysMBS.AddLink(motor);
    motor->SetDisabled(true);


    /////////////////////////////////////////
    // Properties of the granular material //
    /////////////////////////////////////////
    gpu_sys.SetParticles(body_points);

    gpu_sys.SetKn_SPH2SPH(params.normalStiffS2S);
    gpu_sys.SetKn_SPH2WALL(params.normalStiffS2W);
    gpu_sys.SetKn_SPH2MESH(params.normalStiffS2M);

    gpu_sys.SetGn_SPH2SPH(params.normalDampS2S);
    gpu_sys.SetGn_SPH2WALL(params.normalDampS2W);
    gpu_sys.SetGn_SPH2MESH(params.normalDampS2M);

    gpu_sys.SetKt_SPH2SPH(params.tangentStiffS2S);
    gpu_sys.SetKt_SPH2WALL(params.tangentStiffS2W);
    gpu_sys.SetKt_SPH2MESH(params.tangentStiffS2M);

    gpu_sys.SetGt_SPH2SPH(params.tangentDampS2S);
    gpu_sys.SetGt_SPH2WALL(params.tangentDampS2W);
    gpu_sys.SetGt_SPH2MESH(params.tangentDampS2M);

    gpu_sys.SetCohesionRatio(params.cohesion_ratio);

    gpu_sys.SetRollingCoeff_SPH2SPH(params.rolling_friction_coeffS2S);
    gpu_sys.SetRollingCoeff_SPH2WALL(params.rolling_friction_coeffS2W);
    gpu_sys.SetRollingCoeff_SPH2MESH(params.rolling_friction_coeffS2M);

    gpu_sys.SetGravitationalAcceleration(gravity);

    gpu_sys.SetFixedStepSize(params.step_size);
    gpu_sys.SetFrictionMode(CHGPU_FRICTION_MODE::MULTI_STEP); // friction model
    gpu_sys.SetRollingMode(CHGPU_ROLLING_MODE::SCHWARTZ); // rolling friction model
    gpu_sys.SetTimeIntegrator(CHGPU_TIME_INTEGRATOR::CENTERED_DIFFERENCE);

    gpu_sys.SetStaticFrictionCoeff_SPH2SPH(params.static_friction_coeffS2S);
    gpu_sys.SetStaticFrictionCoeff_SPH2WALL(params.static_friction_coeffS2W);
    gpu_sys.SetStaticFrictionCoeff_SPH2MESH(params.static_friction_coeffS2M);

    gpu_sys.SetParticleOutputMode(params.write_mode);
    gpu_sys.SetVerbosity(params.verbose);
    gpu_sys.SetBDFixed(true);

    // In the settling run we disable the mesh.
    gpu_sys.EnableMeshCollision(false);

    gpu_sys.Initialize();

    unsigned int out_steps = (unsigned int)(1 / (out_fps * integration_step));
    unsigned int total_frames = (unsigned int)(params.time_end * out_fps);
    int currframe = 0;
    unsigned int curr_step = 0;

    // Add a ball mesh to the GPU system
    gpu_sys.AddMesh(wheel_obj, VNULL, ChMatrix33<float>(scale_ratio), wheel_mass);
    //gpu_sys.UseMeshNormals(true);
    gpu_sys.Initialize();

    bool enableCollision = false;

    // values to be reported
    ChVector<double> motor_force;
    ChVector<double> motor_torque;
    ChVector<double> wheel_angvel;
    ChVector<double> wheel_velocity;



    for (double t = 0; t < (double)params.time_end; t += integration_step, curr_step++) {
        gpu_sys.ApplyMeshMotion(0, wheel->GetPos(), wheel->GetRot(), wheel->GetPos_dt(), wheel->GetWvel_par());

        ChVector<> gran_force;
        ChVector<> gran_torque;
        gpu_sys.CollectMeshContactForces(0, gran_force, gran_torque);

        wheel->Empty_forces_accumulators();
        wheel->Accumulate_force(gran_force, wheel->GetPos(), false);

        double external_load = (total_cart_load - wheel_mass) * grav;
        
        // ramp up the external load, external load stays constant after 1 second.
        if (t < 1.0)
            external_load = t * external_load;


        ChVector<double> external_force(-sin(slope_rad) * external_load, 0, -cos(slope_rad) * external_load);

        wheel->Accumulate_force(external_force, wheel->GetPos(), false);
        wheel->Accumulate_torque(gran_torque, false);

        // enable mesh and terrain collision after 0.2 seconds
        if (t > 0.05 && enableCollision == false) {
            gpu_sys.EnableMeshCollision(true);
            enableCollision = true;
            wheel->SetBodyFixed(false);
            motor->SetDisabled(false);
            wheel_IniPos.z() = gpu_sys.GetMaxParticleZ() + params.sphere_radius + 0.1 + wheel_radius;
            wheel->SetPos(wheel_IniPos);
            wheel->SetPos_dt(VNULL);
            wheel->SetPos_dtdt(VNULL);
        }


        if (curr_step % out_steps == 0) {
            std::cout << "Output frame " << currframe + 1 << " of " << total_frames << ", time = " << t << " sec, ";
            char filename[300];
            char mesh_filename[300];

            sprintf(filename, "%s/step%06d.csv", out_dir.c_str(), currframe);
            sprintf(mesh_filename, "%s/step%06d_mesh", out_dir.c_str(), currframe++);
            gpu_sys.WriteParticleFile(std::string(filename));
            gpu_sys.WriteMeshes(std::string(mesh_filename));

            motor_force = motor->Get_react_force();
            motor_torque = motor->Get_react_torque();
            wheel_angvel = wheel->GetWvel_par();
            wheel_velocity = wheel->GetPos_dt();

            double slip_ratio = 1 - wheel_velocity.x() / (wheel_AngVel * wheel_radius);


            std::cout << "wheel_velo, " << wheel_velocity.x() / 100 << " m/s, slip_ratio, " << slip_ratio << ", DBP, "
                      << motor_force.z() / 1e5 << " N, torque: " << motor->GetMotorTorque() / 1e7 << " Nm" << std::endl;
        }

        gpu_sys.AdvanceSimulation(integration_step);
        sysMBS.DoStepDynamics(integration_step);
    }

    return 0;
}