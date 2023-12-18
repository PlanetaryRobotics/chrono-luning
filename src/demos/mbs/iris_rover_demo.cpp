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
// Authors: Alessandro Tasora, Radu Serban
// =============================================================================
//
// Demo code about Iris rover on deformable terrain
//
// =============================================================================

#include "chrono/physics/ChSystemSMC.h"
#include "chrono/physics/ChSystemNSC.h"



#include "chrono/physics/ChBodyEasy.h"
#include "chrono/physics/ChInertiaUtils.h"
#include "chrono/assets/ChTexture.h"
#include "chrono/assets/ChTriangleMeshShape.h"
#include "chrono/geometry/ChTriangleMeshConnected.h"
#include "chrono/physics/ChLinkMotorRotationSpeed.h"
#include "chrono_irrlicht/ChVisualSystemIrrlicht.h"
#include "chrono_thirdparty/filesystem/path.h"

#include "chrono_models/robot/iris/Iris.h"


// Use the namespaces of Chrono
using namespace chrono;
using namespace chrono::geometry;
using namespace chrono::irrlicht;

int main(int argc, char* argv[]) {
    GetLog() << "Copyright (c) 2017 projectchrono.org\nChrono version: " << CHRONO_VERSION << "\n\n";

    // wheel position relative to COM of the chassis
    ChVector<double> LF_wheel_pos( 11.5,  11, -3.5);
    ChVector<double> RF_wheel_pos( 11.5, -11, -3.5);
    ChVector<double> LR_wheel_pos(-11.5,  11, -3.5);
    ChVector<double> RR_wheel_pos(-11.5, -11, -3.5);

    std::vector<ChVector<double>> wheel_pos_list;
    wheel_pos_list.push_back(LF_wheel_pos);
    wheel_pos_list.push_back(LR_wheel_pos);
    wheel_pos_list.push_back(RF_wheel_pos);
    wheel_pos_list.push_back(RR_wheel_pos);



    std::vector<std::shared_ptr<ChLinkMotorRotationSpeed>> motor_list;

    double clearance = 4.5;

    ChVector<double> chassis_pos(-5, 0, 0);
    ChVector<double> chassis_dim(25, 17.5, 14);

    double total_mass = 2300;
    double wheel_mass = 28.44;
    double chassis_mass = total_mass - 4 * wheel_mass;

    double motor_ang_velo = CH_C_PI / 5.0;

    // Create a Chrono physical system
    ChSystemNSC sys;

    sys.Set_G_acc(ChVector<>(0, 0, -980));

    // Create all the rigid bodies.

    //collision::ChCollisionModel::SetDefaultSuggestedEnvelope(0.025);
    //collision::ChCollisionModel::SetDefaultSuggestedMargin(0.025);

    // - Create a floor

    auto floor_mat = chrono_types::make_shared<ChMaterialSurfaceNSC>();
    //floor_mat->SetYoungModulus(2.0e5f);
    floor_mat->SetFriction(0.4f);
    auto floor = chrono_types::make_shared<ChBodyEasyBox>(500, 500, 5, 1000, true, true, floor_mat);
    floor->SetPos(ChVector<>(0, 0, chassis_pos.z() - chassis_dim.z() * 0.5 - clearance * 2.5));
    floor->SetBodyFixed(true);
    floor->GetVisualShape(0)->SetTexture(GetChronoDataFile("textures/concrete.jpg"));
    floor->GetVisualShape(0)->SetOpacity(0.8f);
    floor->SetCollide(true);
    sys.Add(floor);


    // - Create a falling item with triangle mesh shape

    auto chassis =
        chrono_types::make_shared<ChBodyEasyBox>(chassis_dim.x(), chassis_dim.y(), chassis_dim.z(), 1, true, false);

    chassis->SetPos(chassis_pos);
    chassis->GetVisualShape(0)->SetColor(ChColor(0.3f, 0.3f, 0.6f));
    sys.Add(chassis);
    chassis->SetBodyFixed(false);
    chassis->SetCollide(false);
    chassis->SetMass(chassis_mass);
    for (int j = 0; j < 4; ++j) {


            //// flip mesh and mesh shape
            // Shared contact material for all meshes
            auto mesh_mat = chrono_types::make_shared<ChMaterialSurfaceNSC>();
            mesh_mat->SetFriction(0.5f);
            //mesh_mat->SetYoungModulus(1.0e7f);
            mesh_mat->SetRestitution(0.2f);

            auto mesh = ChTriangleMeshConnected::CreateFromWavefrontFile(
                "C:/Users/fang/Documents/NSF_Collaboration/CMU_MoonRanger/data/iris_wheel.obj", false, true);
            mesh->Transform(ChVector<>(0, 0, 0), ChMatrix33<>(0.1));  // scale to a different size

             //auto mesh = ChTriangleMeshConnected::CreateFromWavefrontFile(
             //    "C:/Users/fang/Documents/NSF_Collaboration/CMU_MoonRanger/data/wheel.obj", false, true);
             //mesh->Transform(ChVector<>(0, 0, 0), ChMatrix33<>(100));  // scale to a different size



            if (j == 2 || j == 3) {
                ChQuaternion wheel_rot = Q_from_AngAxis(CH_C_PI, VECT_Z);
                mesh->Transform(ChVector<>(0, 0, 0), ChMatrix33(wheel_rot));
            }


            mesh->RepairDuplicateVertexes(1e-9);  // if meshes are not watertight

            // compute mass inertia from mesh
            double mass;
            ChVector<> cog;
            ChMatrix33<> inertia;
            double density = 1.78;
            mesh->ComputeMassProperties(true, mass, cog, inertia);
            ChMatrix33<> principal_inertia_rot;
            ChVector<> principal_I;
            ChInertiaUtils::PrincipalInertia(inertia, principal_I, principal_inertia_rot);

            // Create a shared visual model containing a visualizatoin mesh
            auto mesh_shape = chrono_types::make_shared<ChTriangleMeshShape>();
            mesh_shape->SetMesh(mesh);
            mesh_shape->SetMutable(false);
            mesh_shape->SetColor(ChColor(1.0f, 0.5f, 0.5f));
            mesh_shape->SetBackfaceCull(true);

            auto vis_model = chrono_types::make_shared<ChVisualModel>();
            vis_model->AddShape(mesh_shape);

            auto wheel = chrono_types::make_shared<ChBodyAuxRef>();

            // Set the COG coordinates to barycenter, without displacing the REF reference.
            // Make the COG frame a principal frame.
            wheel->SetFrame_COG_to_REF(ChFrame<>(cog, principal_inertia_rot));

            // Set inertia
            wheel->SetMass(wheel_mass);
            wheel->SetInertiaXX(density * principal_I);

            // Set the absolute position of the body:
            wheel->SetFrame_REF_to_abs(
                ChFrame<>(wheel_pos_list.at(j) + chassis_pos  ));
            sys.Add(wheel);

            wheel->GetCollisionModel()->ClearModel();
            wheel->GetCollisionModel()->AddTriangleMesh(mesh_mat, mesh, false, false, VNULL, ChMatrix33<>(1), 0.005);
            wheel->GetCollisionModel()->BuildModel();
            wheel->SetCollide(true);

            wheel->AddVisualModel(vis_model);

            //wheel->SetBodyFixed(true);


            // now add the rotational motor between wheel and chassis
            // TODO: add the pointer to a list of motors later so i can control them 
            std::shared_ptr<ChLinkMotorRotationSpeed> motor;
            motor = chrono_types::make_shared<ChLinkMotorRotationSpeed>();

            motor->Initialize(wheel, chassis, ChFrame<>(wheel_pos_list.at(j) + chassis_pos, Q_ROTATE_Z_TO_Y));


            // TURN RIGHT

   //         if (j == 0 || j == 1) {
			//	motor->SetMotorFunction(
			//		chrono_types::make_shared<ChFunction_Const>(motor_ang_velo));  // actually, default function type
			//}
			//else {
			//	motor->SetMotorFunction(
			//		chrono_types::make_shared<ChFunction_Const>(-motor_ang_velo));  // actually, default function type
			//}


			motor->SetMotorFunction(
            	chrono_types::make_shared<ChFunction_Const>(motor_ang_velo));  // actually, default function type

            sys.AddLink(motor);

            motor_list.push_back(motor);



    }

    // Shared contact material for falling objects
    //auto obj_mat = chrono_types::make_shared<ChMaterialSurfaceSMC>();
    //obj_mat->SetFriction(0.2f);

    
    std::vector<std::shared_ptr<ChBodyAuxRef>> rock;
    std::vector<std::string> rock_meshfile = {"robot/curiosity/rocks/rock1.obj",
                                              "robot/curiosity/rocks/rock1.obj",  //
                                              "robot/curiosity/rocks/rock1.obj"};

    std::vector<ChVector<>> rock_pos = {
        ChVector<>(chassis_pos.x() + 1.f * chassis_dim.x(), chassis_pos.y() + chassis_dim.y()/2.5, 0), 
        ChVector<>(chassis_pos.x() + 3.f * chassis_dim.x(), chassis_pos.y() - chassis_dim.y()/1.8, 0)};
    std::vector<double> rock_scale = {
        7,  8,   //
        7};
    double rock_density = 3;

    auto rock_mat = chrono_types::make_shared<ChMaterialSurfaceNSC>();

    rock_mat->SetFriction(0.4f);
    //rock_mat->SetRestitution(0.5f);
    //rock_mat->SetYoungModulus(1.0e7f);

    for (int i = 0; i < 2; i++) {
            auto mesh =
                ChTriangleMeshConnected::CreateFromWavefrontFile(GetChronoDataFile(rock_meshfile[i]), false, true);
            mesh->Transform(ChVector<>(0, 0, 0), ChMatrix33<>(rock_scale[i]));

            double mass;
            ChVector<> cog;
            ChMatrix33<> inertia;
            mesh->ComputeMassProperties(true, mass, cog, inertia);
            ChMatrix33<> principal_inertia_rot;
            ChVector<> principal_I;
            ChInertiaUtils::PrincipalInertia(inertia, principal_I, principal_inertia_rot);

            auto body = chrono_types::make_shared<ChBodyAuxRef>();
            sys.Add(body);
            body->SetFrame_REF_to_abs(ChFrame<>(ChVector<>(rock_pos[i]), QUNIT));
            body->SetFrame_COG_to_REF(ChFrame<>(cog, principal_inertia_rot));
            body->SetMass(mass * rock_density);
            body->SetInertiaXX(rock_density * principal_I);

            body->GetCollisionModel()->ClearModel();
            body->GetCollisionModel()->AddTriangleMesh(rock_mat, mesh, false, false, VNULL, ChMatrix33<>(1), 0.005);
            body->GetCollisionModel()->BuildModel();
            body->SetCollide(true);

            auto mesh_shape = chrono_types::make_shared<ChTriangleMeshShape>();
            mesh_shape->SetMesh(mesh);
            mesh_shape->SetBackfaceCull(true);
            body->AddVisualShape(mesh_shape);

            rock.push_back(body);
            
    }







    // Create the Irrlicht visualization system
    auto vis = chrono_types::make_shared<ChVisualSystemIrrlicht>();
    vis->AttachSystem(&sys);
    vis->SetWindowSize(1280, 720);
    vis->SetWindowTitle("Iris rover model");
    vis->Initialize();
    vis->AddLogo();
    vis->AddSkyBox();
    vis->SetCameraVertical(CameraVerticalDir::Z);
    vis->AddCamera(ChVector<>(70, 30, 20));
    vis->AddTypicalLights();
    vis->AddLight(ChVector<>(0, 20, 20), 90, ChColor(0.5f, 0.5f, 0.5f));
    vis->ShowInfoPanel(true);

    //vis->AddTypicalLights(ChVector<>(-10, -10, 20), ChVector<>(-10, 10, 20), 90, 90, 40, 40);
    //vis->AddLightWithShadow(ChVector<>(1.5, 5.5, -2.5), ChVector<>(0, 0, 0), 3, 2.2, 7.2, 40, 512,
    //    ChColor(0.8f, 0.8f, 1.0f));
    //vis->EnableShadows();

        // Default camera uses Z up

    ////application.SetContactsDrawMode(ContactsDrawMode::CONTACT_DISTANCES);

    // Simulation loop
    double end_time = 20;

    double frame_step = 0.01;
    double step_size = 5e-5;
    int image_out_steps = (int)std::ceil(frame_step / step_size);
    int num_steps = 0;

    // create a folder to save images
    std::string out_dir = "iris_rover_demo_rock";
    filesystem::create_directory(filesystem::path(out_dir));

    bool render = false;

    while (vis->Run()) {
        vis->BeginScene(true, true, ChColor(0.55f, 0.63f, 0.75f));
        vis->Render();
        vis->EndScene();
        sys.DoStepDynamics(step_size);

        num_steps++;
        

    //    if (sys.GetChTime() > 2) {
    //        motor_list.at(0)->SetMotorFunction(
				//chrono_types::make_shared<ChFunction_Const>(-motor_ang_velo));  // actually, default function type

    //        motor_list.at(1)->SetMotorFunction(
    //            chrono_types::make_shared<ChFunction_Const>(-motor_ang_velo));  // actually, default function type


    //    }

    //    if (sys.GetChTime() > 7) {
    //        motor_list.at(0)->SetMotorFunction(
    //            chrono_types::make_shared<ChFunction_Const>(motor_ang_velo));  // actually, default function type

    //        motor_list.at(1)->SetMotorFunction(
    //            chrono_types::make_shared<ChFunction_Const>(motor_ang_velo));  // actually, default function type
    //    }

        std::cout << "time: " << sys.GetChTime() << ", chassis pos: " << chassis->GetPos().x()
				  << ", " << chassis->GetPos().y() << ", " << chassis->GetPos().z() << std::endl;

        //std::cout << "time: " << sys.GetChTime() << ", motor spped: [" << motor_list.at(0)->GetMotorRot_dt() << ", "
        //          << motor_list.at(1)->GetMotorRot_dt() << ", " << motor_list.at(2)->GetMotorRot_dt() << ", "
        //          << motor_list.at(3)->GetMotorRot_dt() << " ]" << std::endl;



        // write images every image_out_steps steps
        if (num_steps % image_out_steps == 0 && render) {
			char filename[100];
			sprintf(filename, "%s/img_%04d.jpg", out_dir.c_str(), num_steps / image_out_steps);
			vis->WriteImageToFile(filename);
		}


        if (sys.GetChTime() > end_time) {
            return 0;
        }
    }

    return 0;
}
