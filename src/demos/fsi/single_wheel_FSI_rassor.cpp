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
// Author: Wei Hu, Radu Serban
// =============================================================================

#include <cassert>
#include <cstdlib>
#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <cmath>
#include <fstream>
#include <stdio.h>
#include <cstring>

#include "chrono/ChConfig.h"
#include "chrono/physics/ChSystemSMC.h"
#include "chrono/physics/ChBody.h"
#include "chrono/physics/ChInertiaUtils.h"
#include "chrono/physics/ChLinkMotorRotationAngle.h"
#include "chrono/utils/ChUtilsGeometry.h"
#include "chrono/utils/ChUtilsInputOutput.h"
#include "chrono/assets/ChTriangleMeshShape.h"
#include "chrono/core/ChTimer.h"
#include "chrono_fsi/ChSystemFsi.h"
#include "chrono_thirdparty/filesystem/path.h"
#ifdef CHRONO_VSG
    #include "chrono_fsi/visualization/ChFsiVisualizationVSG.h"
#endif



// Chrono namespaces
using namespace chrono;
using namespace chrono::fsi;
using namespace chrono::geometry;


// Physical properties of terrain particles
double iniSpacing;
double kernelLength;
double density;
double slope_angle; // Terrain slope 

// Dimension of the terrain container
double smalldis = 1.0e-9;
//double bxDim = 5.0 + smalldis;
double bxDim = 2.0 + smalldis;
double byDim = 0.5 + smalldis;   // byDim depending on the wheel width
double bzDim = 0.15 + smalldis;

// Size of the wheel
double wheel_radius = 0.25;
double wheel_wide = 0.205;

double wheel_slip = 0.0;
//double wheel_vel = -0.05;
//double wheel_AngVel = -0.7; // for rTot = 250mm, 0.4 rad/s ~ 0.1 m/s linear velocity

double wheel_vel = -2.4;
double wheel_AngVel = -33.5;  // for rTot = 250mm, 0.4 rad/s ~ 0.1 m/s linear velocity

double total_mass = 17.5;

// Initial Position of wheel
ChVector<> wheel_IniPos(bxDim / 2 - wheel_radius * 1.2, 0.0, wheel_radius + bzDim/2.0 + 0.01);
//ChVector<> wheel_IniVel(0.0, 0.0, -5.0f);
 ChVector<> wheel_IniVel(0.0, 0.0, 0.0f);

// Simulation time and stepsize
//double total_time = 5.0;
//double total_time = 60.0;
double total_time = 10;
double dT;

// linear actuator and angular actuator
auto actuator = chrono_types::make_shared<ChLinkLinActuator>();
auto motor = chrono_types::make_shared<ChLinkMotorRotationAngle>();

// Save data as csv files to see the results off-line using Paraview
bool output = true;
int out_fps = 100;

// Output directories and settings
std::string out_dir = "/root/sbel/outputs/FSI_Single_Wheel_Test_RealSlope_mode_slope";

// Enable/disable run-time visualization (if Chrono::OpenGL is available)
bool render = false;
float render_fps = 100;

// Verbose terminal output
bool verbose_fsi = true;
bool verbose = true;

std::string wheel_obj = "C:/Users/fang/Documents/RASSOR/test_32k_triangles.obj";
//std::string drum_BCE_csvfile = "C:/Users/fang/Documents/RASSOR/preprocessing/drum_BCE.csv";
std::string drum_BCE_csvfile = "C:/Users/fang/Documents/RASSOR/preprocessing/single_drum_0.01.csv";

//std::string drum_BCE_csvfile = "C:/Users/fang/Documents/RASSOR/DEM_solver/ViperWheelSimple.csv";

//------------------------------------------------------------------
// Create the objects of the MBD system. Rigid bodies, and if FSI,
// their BCE representation are created and added to the systems
//------------------------------------------------------------------
void CreateSolidPhase(ChSystemSMC& sysMBS, ChSystemFsi& sysFSI) {
    // Common contact material
    auto cmaterial = chrono_types::make_shared<ChMaterialSurfaceSMC>();
    cmaterial->SetYoungModulus(1e8);
    cmaterial->SetFriction(0.9f);
    cmaterial->SetRestitution(0.4f);
    cmaterial->SetAdhesion(0);

    // Create a container -- always FIRST body in the system
    auto ground = chrono_types::make_shared<ChBody>();
    ground->SetIdentifier(-1);
    ground->SetBodyFixed(true);
    sysMBS.AddBody(ground);

    ground->GetCollisionModel()->ClearModel();
    chrono::utils::AddBoxContainer(ground, cmaterial, ChFrame<>(), ChVector<>(bxDim, byDim, bzDim), 0.1,
        ChVector<int>(0, 0, -1), false);
    ground->GetCollisionModel()->BuildModel();
    ground->SetCollide(false);

    // Add BCE particles attached on the walls into FSI system
    sysFSI.AddBoxContainerBCE(ground, ChFrame<>(), ChVector<>(bxDim, byDim, bzDim), ChVector<int>(2, 0, -1));

    // Create the wheel -- always SECOND body in the system
    auto trimesh = chrono_types::make_shared<ChTriangleMeshConnected>();
    double scale_ratio = 1e-3;
    trimesh->LoadWavefrontMesh(GetChronoDataFile(wheel_obj), false, true);
    trimesh->Transform(ChVector<>(0, 0, 0), ChMatrix33<>(scale_ratio));  // scale to a different size
    trimesh->RepairDuplicateVertexes(1e-9);                              // if meshes are not watertight

    // Compute mass inertia from mesh
    double mmass;
    double mdensity = 1500.0;
    ChVector<> mcog;
    ChMatrix33<> minertia;
    trimesh->ComputeMassProperties(true, mmass, mcog, minertia);
    ChMatrix33<> principal_inertia_rot;
    ChVector<> principal_I;
    ChInertiaUtils::PrincipalInertia(minertia, principal_I, principal_inertia_rot);

    // look into principal inertia computation!!! why it's so small!!!!
    principal_I.x() = principal_I.z();
    mcog = ChVector<>(0.0, 0.0, 0.0);

    // Set the abs orientation, position and velocity
    auto wheel = chrono_types::make_shared<ChBodyAuxRef>();
    ChQuaternion<> wheel_Rot = Q_from_Euler123(ChVector<double>(0, 0, 0));

    // Set the COG coordinates to barycenter, without displacing the REF reference.
    // Make the COG frame a principal frame.
    wheel->SetFrame_COG_to_REF(ChFrame<>(mcog, principal_inertia_rot));

    // Set inertia
    wheel->SetMass(total_mass * 1.0 / 2.0);
    wheel->SetInertiaXX(mdensity * principal_I);
    std::cout << "principal inertia: " << std::endl;
    std::cout << mdensity * principal_I << std::endl;
    wheel->SetPos_dt(wheel_IniVel);
    wheel->SetWvel_loc(ChVector<>(0.0, 0.0, 0.0));  // set an initial anular velocity (rad/s)

    // Set the absolute position of the body:
    wheel->SetFrame_REF_to_abs(ChFrame<>(ChVector<>(wheel_IniPos), ChQuaternion<>(wheel_Rot)));
    sysMBS.AddBody(wheel);

    wheel->SetBodyFixed(false);
    wheel->GetCollisionModel()->ClearModel();
    wheel->GetCollisionModel()->AddTriangleMesh(cmaterial, trimesh, false, false, VNULL, ChMatrix33<>(1), 0.005);
    wheel->GetCollisionModel()->BuildModel();
    wheel->SetCollide(false);

    // Add this body to the FSI system
     std::vector<ChVector<>> BCE_wheel;
     // read the drum BCE from csv file, add points to the BCE_wheel vector
     std::ifstream file(drum_BCE_csvfile);
     std::string line;
     std::getline(file, line); // skip the first line
     while (std::getline(file, line))
     {
		 std::stringstream iss(line);
		 std::string val;
		 std::vector<double> values;
         while (std::getline(iss, val, ','))
         {
			 values.push_back(std::stod(val));
		 }
		 BCE_wheel.push_back(ChVector<>(values[0], values[1], values[2]));
         //std::cout << values[0] << " " << values[1] << " " << values[2] << std::endl;
	 }
     sysFSI.AddPointsBCE(wheel, BCE_wheel, ChFrame<>(), true);
     sysFSI.AddFsiBody(wheel);

    std::cout << "Added BCE particles for rassor wheel" << std::endl;

    // Create the chassis -- always THIRD body in the system
    auto chassis = chrono_types::make_shared<ChBody>();
    chassis->SetMass(total_mass/2.0f);
    chassis->SetPos(wheel->GetPos());
    chassis->SetCollide(false);
    chassis->SetBodyFixed(false);

    // Add geometry of the chassis.
    chassis->GetCollisionModel()->ClearModel();
    chrono::utils::AddBoxGeometry(chassis.get(), cmaterial, ChVector<>(0.1, 0.1, 0.1), ChVector<>(0, 0, 0));
    chassis->GetCollisionModel()->BuildModel();
    sysMBS.AddBody(chassis);

    // Create the axle -- always FOURTH body in the system
    auto axle = chrono_types::make_shared<ChBody>();
    axle->SetMass(total_mass * 1.0 / 2.0);
    axle->SetPos(wheel->GetPos());
    axle->SetCollide(false);
    axle->SetBodyFixed(false);

    // Add geometry of the axle.
    axle->GetCollisionModel()->ClearModel();
    chrono::utils::AddSphereGeometry(axle.get(), cmaterial, 0.5, ChVector<>(0, 0, 0));
    axle->GetCollisionModel()->BuildModel();
    sysMBS.AddBody(axle);

    // Connect the chassis to the containing bin (ground) through a translational joint and create a linear actuator.
    auto prismatic1 = chrono_types::make_shared<ChLinkLockPrismatic>();
    prismatic1->Initialize(ground, chassis, ChCoordsys<>(chassis->GetPos(), Q_from_AngY(CH_C_PI_2)));
    prismatic1->SetName("prismatic_chassis_ground");
    sysMBS.AddLink(prismatic1);

    double velocity = wheel_vel;  // wheel_AngVel * wheel_radius * (1.0 - wheel_slip);
    auto actuator_fun = chrono_types::make_shared<ChFunction_Ramp>(0.0, velocity);

    actuator->Initialize(ground, chassis, false, ChCoordsys<>(chassis->GetPos(), QUNIT),
                         ChCoordsys<>(chassis->GetPos() + ChVector<>(1, 0, 0), QUNIT));
    actuator->SetName("actuator");
    actuator->SetDistanceOffset(1);
    actuator->SetActuatorFunction(actuator_fun);
    sysMBS.AddLink(actuator);



    // Connect the axle to the chassis through a vertical translational joint.
    auto prismatic2 = chrono_types::make_shared<ChLinkLockPrismatic>();
    prismatic2->Initialize(chassis, axle, ChCoordsys<>(chassis->GetPos(), QUNIT));
    prismatic2->SetName("prismatic_axle_chassis");
    sysMBS.AddLink(prismatic2);

    // Connect the wheel to the axle through a engine joint.
    motor->SetName("engine_wheel_axle");
    motor->Initialize(wheel, axle, ChFrame<>(wheel->GetPos(),
        chrono::Q_from_AngAxis(-CH_C_PI / 2.0, ChVector<>(1, 0, 0))));
    motor->SetAngleFunction(chrono_types::make_shared<ChFunction_Ramp>(0, wheel_AngVel));
    //sysMBS.AddLink(motor);
}

// =============================================================================

int main(int argc, char* argv[]) {
    // The path to the Chrono data directory
    //SetChronoDataPath(CHRONO_DATA_DIR);
    SetChronoDataPath("");

    // Create the MBS and FSI systems
    ChSystemSMC sysMBS;
    ChSystemFsi sysFSI(&sysMBS);

    sysFSI.SetVerbose(verbose_fsi);

    // Use JSON file to set the FSI parameters
    std::string inputJson = "C:/Users/fang/Documents/RASSOR/input_json/demo_FSI_SingleDrum.json";

    total_mass = 17.5;
    slope_angle = 0;
    out_dir = "C:/Users/fang/Documents/RASSOR/fsi_outputs/debug";
    //wheelfilename = "C:/Users/fang/source/nasa_singlewheel/wheel_obj/withGrousers.obj";

    std::cout << "total_mass: " << total_mass << "\n";
    std::cout << "slope_angle: " << slope_angle << "\n";
    std::cout << "wheel_angvel: " << wheel_AngVel << "\n";
    std::cout << "out_dir: " << out_dir << "\n";
    std::cout << "out_fps: " << out_fps << "\n";

    // Create oputput directories
    if (!filesystem::create_subdirectory(filesystem::path(out_dir))) {
        std::cerr << "Error creating directory " << out_dir << std::endl;
        return 1;
    }
    if (!filesystem::create_subdirectory(filesystem::path(out_dir + "/particles"))) {
        std::cerr << "Error creating directory " << out_dir + "/particles" << std::endl;
        return 1;
    }
    if (!filesystem::create_subdirectory(filesystem::path(out_dir + "/fsi"))) {
        std::cerr << "Error creating directory " << out_dir + "/fsi" << std::endl;
        return 1;
    }
    std::cout << "Reading json file" << std::endl;

    sysFSI.ReadParametersFromFile(inputJson);
    std::cout << "Done reading json file" << std::endl;
    
    double gravity_G = sysFSI.Get_G_acc().z();
    ChVector<> gravity = ChVector<>(gravity_G * sin(slope_angle), 0, gravity_G * cos(slope_angle));
    sysMBS.Set_G_acc(gravity);
    sysFSI.Set_G_acc(gravity);

    // Get the simulation stepsize
    dT = sysFSI.GetStepSize();

    // Get the initial particle spacing
    iniSpacing = sysFSI.GetInitialSpacing();

    // Get the SPH kernel length
    kernelLength = sysFSI.GetKernelLength();

    // Set the terrain container size
    sysFSI.SetContainerDim(ChVector<>(bxDim, byDim, bzDim));

    // Set SPH discretization type, consistent or inconsistent
    sysFSI.SetDiscreType(false, false);

    // Set wall boundary condition
    sysFSI.SetWallBC(BceVersion::ADAMI);

    // Set rigid body boundary condition
    sysFSI.SetRigidBodyBC(BceVersion::ADAMI);

    // Set cohsion of the granular material
    sysFSI.SetCohesionForce(0.0);

    // Setup the SPH method
    sysFSI.SetSPHMethod(FluidDynamics::WCSPH);

    // Set up the periodic boundary condition (if not, set relative larger values)
    ChVector<> cMin(-bxDim / 2 * 10, -byDim / 2 - 0.5 * iniSpacing, -bzDim * 10);
    ChVector<> cMax( bxDim / 2 * 10,  byDim / 2 + 0.5 * iniSpacing, bzDim * 10);
    sysFSI.SetBoundaries(cMin, cMax);

    // Initialize the SPH particles
    ChVector<> boxCenter(0.0, 0.0, 0.0);
    ChVector<> boxHalfDim(bxDim / 2, byDim / 2, bzDim / 2);
    sysFSI.AddBoxSPH(boxCenter, boxHalfDim);

    // Create Solid region and attach BCE SPH particles
    CreateSolidPhase(sysMBS, sysFSI);

    // Set simulation data output length
    sysFSI.SetOutputLength(0);

    // Construction of the FSI system must be finalized before running
    sysFSI.Initialize();

    auto wheel = sysMBS.Get_bodylist()[1];
    ChVector<> force = actuator->Get_react_force();
    ChVector<> torque = motor->Get_react_torque();
    ChVector<> w_pos = wheel->GetPos();
    ChVector<> w_vel = wheel->GetPos_dt();
    ChVector<> angvel = wheel->GetWvel_loc();

    ChVector<> w_pos_init = wheel->GetPos();

    // Save wheel mesh
    //ChTriangleMeshConnected wheel_mesh;
    //wheel_mesh.LoadWavefrontMesh(GetChronoDataFile(wheel_obj), false, true);
    //wheel_mesh.RepairDuplicateVertexes(1e-9);

    // Write the information into a txt file
    std::ofstream myFile;
    std::ofstream myDBP_Torque;
    if (output) {
        myFile.open(out_dir + "/results.txt", std::ios::trunc);
        myDBP_Torque.open(out_dir + "/DBP_Torque.txt", std::ios::trunc);
    }

    // Create a run-tme visualizer
    std::shared_ptr<ChFsiVisualization> fsi_vis;
    // attach it to the FSI system
    fsi_vis = chrono_types::make_shared<ChFsiVisualizationVSG>(&sysFSI);

    if (render) {
        fsi_vis->SetTitle("Chrono::FSI single wheel demo");
        fsi_vis->SetSize(1280, 720);

        fsi_vis->AddCamera(ChVector<>(0, -5 * byDim, 5 * bzDim), ChVector<>(0, 0, 0));
        fsi_vis->SetCameraMoveScale(0.05f);
        fsi_vis->EnableBoundaryMarkers(true);
        fsi_vis->EnableBoundaryMarkers(false);
        fsi_vis->EnableRigidBodyMarkers(true);
        fsi_vis->SetRenderMode(ChFsiVisualization::RenderMode::SOLID);
        fsi_vis->SetParticleRenderMode(ChFsiVisualization::RenderMode::SOLID);
        fsi_vis->Initialize();
    }

    // Start the simulation
    unsigned int output_steps = (unsigned int)round(1 / (out_fps * dT));
    unsigned int render_steps = (unsigned int)round(1 / (render_fps * dT));

    double time = 0.0;
    int current_step = 0;


    while (time < total_time) {
        // Get the infomation of the wheel
        force = actuator->Get_react_force();
        torque = motor->Get_react_torque();
        w_pos = wheel->GetPos();
        w_vel = wheel->GetPos_dt();
        angvel = wheel->GetWvel_loc();

        if (time < 0.1)
        {
            w_pos_init = wheel->GetPos();
        }

        //}

        if (output) {
            myFile << time << "\t" << w_pos.x() << "\t" << w_pos.y() << "\t" << w_pos.z() << "\t" << w_vel.x() << "\t"
                << w_vel.y() << "\t" << w_vel.z() << "\t" << angvel.x() << "\t" << angvel.y() << "\t" << angvel.z()
                << "\t" << force.x() << "\t" << force.y() << "\t" << force.z() << "\t" << torque.x() << "\t"
                << torque.y() << "\t" << torque.z() << "\n";
            myDBP_Torque << time << "\t" << force.x() << "\t" << torque.z() << "\n";
        }

        if (output && current_step % output_steps == 0) {
         //if (output) {

            std::cout << "-------- Output" << std::endl;
            sysFSI.PrintParticleToFile(out_dir + "/particles");
            sysFSI.PrintFsiInfoToFile(out_dir + "/fsi", time);
            static int counter = 0;        
            std::cout << "time: " << time << std::endl;
            std::cout << "  wheel position:         " << w_pos << std::endl;
            std::cout << "  wheel linear velocity:  " << w_vel << std::endl;
            std::cout << "  wheel angular velocity: " << angvel << std::endl;
            std::cout << "  drawbar pull:           " << force << std::endl;
            std::cout << "  wheel torque:           " << torque << std::endl;
        
        
        
        }

        // Call the FSI solver
        sysFSI.DoStepDynamics_FSI();
        time += dT;
        current_step++;

        if (w_pos.x() + wheel_radius > bxDim / 2.0f) {
            break;
        }

        if (time > 0.4168 && time < 0.42) {
            std::cout << "time: " << time << std::endl;
            std::cout << "  wheel position:         " << w_pos << std::endl;
            std::cout << "  wheel linear velocity:  " << w_vel << std::endl;
            std::cout << "  wheel angular velocity: " << angvel << std::endl;
            std::cout << "  drawbar pull:           " << force << std::endl;
            std::cout << "  wheel torque:           " << torque << std::endl;

            sysMBS.EnableSolverMatrixWrite(true, "C:/Users/fang/Documents/RASSOR/build/bin/Release/matrix_solve");

        }

        if (render && current_step % render_steps == 0) {
            if (!fsi_vis->Render())
                break;
        }

        if (time > 0.43) {
            break;
        }

    }

    if (output) {
        myFile.close();
        myDBP_Torque.close();
    }

    return 0;
}