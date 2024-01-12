// =============================================================================
// PROJECT CHRONO - http://projectchrono.org
//
// Copyright (c) 2021 projectchrono.org
// All right reserved.
//
// Use of this source code is governed by a BSD-style license that can be found
// in the LICENSE file at the top level of the distribution and at
// http://projectchrono.org/license-chrono.txt.
//
// =============================================================================
// Authors: Luning Bakke
// =============================================================================
// Demo to show Iris Rover operated on Rigid Terrain
// This Demo includes operation to spawn an Iris rover, control wheel speed
// =============================================================================
// TODO: different wheel types that uses different meshes

#include "chrono_models/robot/iris/Iris.h"

#include "chrono/physics/ChSystemNSC.h"
#include "chrono/physics/ChBodyEasy.h"

#ifdef CHRONO_POSTPROCESS
    #include "chrono_postprocess/ChGnuPlot.h"
#endif

#include "chrono/assets/ChVisualSystem.h"
#ifdef CHRONO_IRRLICHT
    #include "chrono_irrlicht/ChVisualSystemIrrlicht.h"
using namespace chrono::irrlicht;
#endif

using namespace chrono;
using namespace chrono::iris;

// -----------------------------------------------------------------------------

// Define Iris rover wheel type
IrisWheelType wheel_type = IrisWheelType::RealWheel;

// Simulation time step
double time_step = 1e-3;

// -----------------------------------------------------------------------------
// contact material, define friction coefficienn, cor, etc.
std::shared_ptr<ChMaterialSurface> CustomWheelMaterial(ChContactMethod contact_method) {
    float mu = 0.4f;   // coefficient of friction
    float cr = 0.1f;   // coefficient of restitution
    float Y = 2e7f;    // Young's modulus
    float nu = 0.3f;   // Poisson ratio
    float kn = 2e5f;   // normal stiffness
    float gn = 40.0f;  // normal viscous damping
    float kt = 2e5f;   // tangential stiffness
    float gt = 20.0f;  // tangential viscous damping

    switch (contact_method) {
        case ChContactMethod::NSC: {
            auto matNSC = chrono_types::make_shared<ChMaterialSurfaceNSC>();
            matNSC->SetFriction(mu);
            matNSC->SetRestitution(cr);
            return matNSC;
        }
        case ChContactMethod::SMC: {
            auto matSMC = chrono_types::make_shared<ChMaterialSurfaceSMC>();
            matSMC->SetFriction(mu);
            matSMC->SetRestitution(cr);
            matSMC->SetYoungModulus(Y);
            matSMC->SetPoissonRatio(nu);
            matSMC->SetKn(kn);
            matSMC->SetGn(gn);
            matSMC->SetKt(kt);
            matSMC->SetGt(gt);
            return matSMC;
        }
        default:
            return std::shared_ptr<ChMaterialSurface>();
    }
}

int main(int argc, char* argv[]) {
    GetLog() << "Copyright (c) 2017 projectchrono.org\nChrono version: " << CHRONO_VERSION << "\n\n";

    // Create the Chrono system with gravity in the negative Z direction
    ChSystemNSC sys;
    sys.Set_G_acc(ChVector<>(0, 0, -9.81));

    sys.SetCollisionSystemType(ChCollisionSystem::Type::BULLET);
    ChCollisionModel::SetDefaultSuggestedEnvelope(0.0025);
    ChCollisionModel::SetDefaultSuggestedMargin(0.0025);

    // Create the ground.
    auto ground_mat = chrono_types::make_shared<ChMaterialSurfaceNSC>();
    auto ground = chrono_types::make_shared<ChBodyEasyBox>(30, 30, 1, 1000, true, true, ground_mat);
    ground->SetPos(ChVector<>(0, 0, -0.5));
    ground->SetBodyFixed(true);
    ground->GetVisualShape(0)->SetTexture(GetChronoDataFile("textures/concrete.jpg"), 60, 45);
    sys.Add(ground);

    // Construct a Iris rover and the asociated driver
    // anguler velocity reaches 0.4 rad/s after 1 second of simulation
    auto driver = chrono_types::make_shared<IrisSpeedDriver>(1.0, 0.4f);

    Iris iris(&sys, wheel_type);
    iris.SetSpeedDriver(driver);
    iris.SetChassisFixed(false);
    iris.SetWheelContactMaterial(CustomWheelMaterial(ChContactMethod::NSC));

    iris.SetChassisVisualization(false);

    // Initialize the rover at position (0, 0, 0.14)
    iris.Initialize(ChFrame<>(ChVector<>(0, 0, 0.14), QUNIT));

    std::cout << "iris total mass: " << iris.GetRoverMass() << std::endl;
    std::cout << "  chassis:        " << iris.GetChassis()->GetBody()->GetMass() << std::endl;
    std::cout << "  wheel:          " << iris.GetWheel(IrisWheelID::LF)->GetBody()->GetMass() << std::endl;
    std::cout << std::endl;

    // Create the run-time visualization interface
    auto vis = chrono_types::make_shared<ChVisualSystemIrrlicht>();
    vis->AttachSystem(&sys);
    vis->SetCameraVertical(CameraVerticalDir::Z);
    vis->SetWindowSize(800, 600);
    vis->SetWindowTitle("Iris Rover on Rigid Terrain");
    vis->Initialize();
    vis->AddLogo();
    vis->AddSkyBox();
    vis->AddCamera(ChVector<>(3, 3, 1));
    vis->AddTypicalLights();
    vis->EnableContactDrawing(ContactsDrawMode::CONTACT_DISTANCES);
    vis->EnableShadows();

    // Simulation loop
    while (vis->Run()) {
        vis->BeginScene();
        vis->Render();
        vis->EndScene();

        // Set current steering angle
        double time = iris.GetSystem()->GetChTime();
        // Update the rover controls
        iris.Update();

        std::cout << sys.GetChTime() << ", " << iris.GetChassisVel() << ", " << std::endl;
        sys.DoStepDynamics(time_step);
    }

    return 0;
}
