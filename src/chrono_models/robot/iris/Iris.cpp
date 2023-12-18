// =============================================================================
// PROJECT CHRONO - http://projectchrono.org
//
// Copyright (c) 2021 projectchrono.org
// All rights reserved.
//
// Use of this source code is governed by a BSD-style license that can be found
// in the LICENSE file at the top level of the distribution and at
// http://projectchrono.org/license-chrono.txt.
//
// =============================================================================
// Authors: Luning Bakke
// =============================================================================
//
// CMU IRIS Lunar Rover Model Class.
// This class contains model for CMU's IRIS lunar rover for NASA's CLPS mission
// exploration mission.
//
// Luning TODO:
// 
// User can specify wheel obj files and mass properties
// Move wheel frame info to private member of iris class
// Chassis mesh
// 
// =============================================================================
//
// =============================================================================

#include <cmath>

#include "chrono/assets/ChBoxShape.h"
#include "chrono/physics/ChBodyEasy.h"
#include "chrono/assets/ChCylinderShape.h"
#include "chrono/assets/ChSphereShape.h"
#include "chrono/assets/ChTexture.h"
#include "chrono/assets/ChTriangleMeshShape.h"

#include "chrono/motion_functions/ChFunction_Setpoint.h"

#include "chrono/physics/ChLinkMotorRotationAngle.h"
#include "chrono/physics/ChLinkMotorRotationSpeed.h"
#include "chrono/physics/ChLinkMotorRotationTorque.h"
#include "chrono/physics/ChShaftsBody.h"

#include "chrono/physics/ChInertiaUtils.h"

#include "chrono_models/robot/iris/Iris.h"

namespace chrono {
namespace iris {


// =============================================================================

// Default contact material for rover parts
std::shared_ptr<ChMaterialSurface> DefaultContactMaterial(ChContactMethod contact_method) {
    float mu = 0.4f;   // coefficient of friction
    float cr = 0.0f;   // coefficient of restitution
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

// Add a rotational speed motor between two bodies at the given position and orientation
// (expressed in and relative to the chassis frame).
std::shared_ptr<ChLinkMotorRotationSpeed> AddMotorSpeed(std::shared_ptr<ChBody> body1,
                                                        std::shared_ptr<ChBody> body2,
                                                        std::shared_ptr<IrisChassis> chassis,
                                                        const ChVector<>& joint_pos,
                                                        const ChQuaternion<>& joint_rot) {
    // Express relative frame in global
    ChFrame<> X_GC = chassis->GetBody()->GetFrame_REF_to_abs() * ChFrame<>(joint_pos, joint_rot);

    // Create motor (actuated DOF about Z axis of X_GC frame)
    auto motor = chrono_types::make_shared<ChLinkMotorRotationSpeed>();
    motor->Initialize(body1, body2, X_GC);
    chassis->GetBody()->GetSystem()->AddLink(motor);

    return motor;
}

// Add a rotational angle motor between two bodies at the given position and orientation
// (expressed in and relative to the chassis frame).
std::shared_ptr<ChLinkMotorRotationAngle> AddMotorAngle(std::shared_ptr<ChBody> body1,
                                                        std::shared_ptr<ChBody> body2,
                                                        std::shared_ptr<IrisChassis> chassis,
                                                        const ChVector<>& rel_pos,
                                                        const ChQuaternion<>& rel_rot) {
    // Express relative frame in global
    ChFrame<> X_GC = chassis->GetBody()->GetFrame_REF_to_abs() * ChFrame<>(rel_pos, rel_rot);

    // Create motor (actuated DOF about Z axis of X_GC frame)
    auto motor = chrono_types::make_shared<ChLinkMotorRotationAngle>();
    motor->Initialize(body1, body2, X_GC);
    chassis->GetBody()->GetSystem()->AddLink(motor);

    return motor;
}

// Add a rotational torque motor between two bodies at the given position and orientation
// (expressed in and relative to the chassis frame).
std::shared_ptr<ChLinkMotorRotationTorque> AddMotorTorque(std::shared_ptr<ChBody> body1,
                                                          std::shared_ptr<ChBody> body2,
                                                          std::shared_ptr<IrisChassis> chassis,
                                                          const ChVector<>& rel_pos,
                                                          const ChQuaternion<>& rel_rot) {
    // Express relative frame in global
    ChFrame<> X_GC = chassis->GetBody()->GetFrame_REF_to_abs() * ChFrame<>(rel_pos, rel_rot);

    // Create motor (actuated DOF about Z axis of X_GC frame)
    auto motor = chrono_types::make_shared<ChLinkMotorRotationTorque>();
    motor->Initialize(body1, body2, X_GC);
    chassis->GetBody()->GetSystem()->AddLink(motor);

    return motor;
}
// =============================================================================

// Base class for all Iris Part, position relative to the chassis frame
IrisPart::IrisPart(const std::string& name,
                     const ChFrame<>& rel_pos,
                     std::shared_ptr<ChMaterialSurface> mat,
                     bool collide)
    : m_name(name), m_pos(rel_pos), m_mat(mat), m_collide(collide), m_visualize(true) {}

void IrisPart::Construct(ChSystem* system) {
    m_body = std::shared_ptr<ChBodyAuxRef>(system->NewBodyAuxRef());
    m_body->SetNameString(m_name + "_body");
    m_body->SetMass(m_mass);
    m_body->SetInertiaXX(m_inertia);
    m_body->SetFrame_COG_to_REF(m_cog);

    // Add visualization shape
    if (m_visualize) {
        auto vis_mesh_file = GetChronoDataFile("robot/iris/obj/" + m_mesh_name + ".obj");
        auto trimesh_vis = geometry::ChTriangleMeshConnected::CreateFromWavefrontFile(vis_mesh_file, true, true);
        trimesh_vis->Transform(m_mesh_xform.GetPos(), m_mesh_xform.GetA());  // translate/rotate/scale mesh
        trimesh_vis->RepairDuplicateVertexes(1e-9);                          // if meshes are not watertight

        auto trimesh_shape = chrono_types::make_shared<ChTriangleMeshShape>();
        trimesh_shape->SetMesh(trimesh_vis);
        trimesh_shape->SetName(m_mesh_name);
        trimesh_shape->SetMutable(false);
        m_body->AddVisualShape(trimesh_shape);
    }

    // Add collision shape
    if (m_collide) {
        auto col_mesh_file = GetChronoDataFile("robot/iris/" + m_mesh_name + ".obj");
        auto trimesh_col = geometry::ChTriangleMeshConnected::CreateFromWavefrontFile(col_mesh_file, false, false);
        trimesh_col->Transform(m_mesh_xform.GetPos(), m_mesh_xform.GetA());  // translate/rotate/scale mesh
        trimesh_col->RepairDuplicateVertexes(1e-9);                          // if meshes are not watertight

        m_body->GetCollisionModel()->ClearModel();
        m_body->GetCollisionModel()->AddTriangleMesh(m_mat, trimesh_col, false, false, VNULL, ChMatrix33<>(1), 0.005);
        m_body->GetCollisionModel()->BuildModel();
        m_body->SetCollide(m_collide);
    }

    system->AddBody(m_body);
}

void IrisPart::CalcMassProperties(double density) {
    auto mesh_filename = GetChronoDataFile("robot/iris/col/" + m_mesh_name + ".obj");
    auto trimesh_col = geometry::ChTriangleMeshConnected::CreateFromWavefrontFile(mesh_filename, false, false);
    trimesh_col->Transform(m_mesh_xform.GetPos(), m_mesh_xform.GetA());  // translate/rotate/scale mesh
    trimesh_col->RepairDuplicateVertexes(1e-9);                          // if meshes are not watertight

    double vol;
    ChVector<> cog_pos;
    ChMatrix33<> cog_rot;
    ChMatrix33<> inertia;
    trimesh_col->ComputeMassProperties(true, vol, cog_pos, inertia);
    ChInertiaUtils::PrincipalInertia(inertia, m_inertia, cog_rot);
    m_mass = density * vol;
    m_inertia *= density;
    m_cog = ChFrame<>(cog_pos, cog_rot);
}

void IrisPart::Initialize(std::shared_ptr<ChBodyAuxRef> chassis) {
    Construct(chassis->GetSystem());

    // Set absolute position
    ChFrame<> X_GC = chassis->GetFrame_REF_to_abs() * m_pos;
    m_body->SetFrame_REF_to_abs(X_GC);
}

// =============================================================================

// Rover Chassis
IrisChassis::IrisChassis(const std::string& name, std::shared_ptr<ChMaterialSurface> mat)
    : IrisPart(name, ChFrame<>(VNULL, QUNIT), mat, false) {
    m_mesh_name = "iris_chassis";
    m_color = ChColor(1.0f, 1.0f, 1.0f);
    CalcMassProperties(165);
}

void IrisChassis::Initialize(ChSystem* system, const ChFrame<>& pos) {
    Construct(system);

    m_body->SetFrame_REF_to_abs(pos);
}

// =============================================================================

// Iris Wheel
IrisWheel::IrisWheel(const std::string& name,
                       const ChFrame<>& rel_pos,
                       std::shared_ptr<ChMaterialSurface> mat
                       )
    : IrisPart(name, rel_pos, mat, true) {

    m_color = ChColor(0.4f, 0.7f, 0.4f);
    CalcMassProperties(800);
}


// =============================================================================

// =============================================================================

// Rover model
Iris::Iris(ChSystem* system) : m_system(system), m_chassis_fixed(false) {
    // Set default collision model envelope commensurate with model dimensions.
    // Note that an SMC system automatically sets envelope to 0.
    auto contact_method = m_system->GetContactMethod();
    if (contact_method == ChContactMethod::NSC) {
        collision::ChCollisionModel::SetDefaultSuggestedEnvelope(0.01);
        collision::ChCollisionModel::SetDefaultSuggestedMargin(0.005);
    }

    // Create the contact materials
    m_default_material = DefaultContactMaterial(contact_method);
    m_wheel_material = DefaultContactMaterial(contact_method);

    Create();
}

void Iris::Create() {
    // create rover chassis
    m_chassis = chrono_types::make_shared<IrisChassis>("chassis", m_default_material);

    // initilize rover wheels
    double wx = 0.5618 + 0.08;
    double wy = 0.2067 + 0.32 + 0.0831;
    double wz = 0.0;

    m_wheels[LF] = chrono_types::make_shared<IrisWheel>("wheel_LF", ChFrame<>(ChVector<>(+wx, +wy, wz), QUNIT),
                                                           m_wheel_material);
    m_wheels[RF] = chrono_types::make_shared<IrisWheel>("wheel_RF", ChFrame<>(ChVector<>(+wx, -wy, wz), QUNIT),
                                                           m_wheel_material);
    m_wheels[LB] = chrono_types::make_shared<IrisWheel>("wheel_LB", ChFrame<>(ChVector<>(-wx, +wy, wz), QUNIT),
                                                           m_wheel_material);
    m_wheels[RB] = chrono_types::make_shared<IrisWheel>("wheel_RB", ChFrame<>(ChVector<>(-wx, -wy, wz), QUNIT),
                                                           m_wheel_material);

    m_wheels[LF]->m_mesh_xform = ChFrame<>(VNULL, Q_from_AngZ(CH_C_PI));
    m_wheels[LB]->m_mesh_xform = ChFrame<>(VNULL, Q_from_AngZ(CH_C_PI));


}

void Iris::Initialize(const ChFrame<>& pos) {
    assert(m_driver);

    m_chassis->Initialize(m_system, pos);
    m_chassis->GetBody()->SetBodyFixed(m_chassis_fixed);

    for (int i = 0; i < 4; i++) {
        m_wheels[i]->Initialize(m_chassis->GetBody());
    }

    // wheel relative position w.r.t chassis COM
    double w_lx =  11.5;
    double w_ly =  11;
    double w_lz = -3.5;

    ChVector<> wheel_rel_pos[] = {
        ChVector<>(+w_lx, +w_ly, w_lz),  // LF
        ChVector<>(+w_lx, -w_ly, w_lz),  // RF
        ChVector<>(-w_lx, +w_ly, w_lz),  // LB
        ChVector<>(-w_lx, -w_ly, w_lz)   // RB
    };

    // Orientation of driving motors.
    // TODO: Flip coordinate system so positive speed is forward.
    ChQuaternion<> driving_rot[] = {
        QUNIT,                 // LF
        QUNIT,                 // RF
        Q_from_AngX(CH_C_PI),  // LB
        Q_from_AngX(CH_C_PI)   // RB
    };


    ChQuaternion<> z2x = Q_from_AngY(CH_C_PI_2);

    for (int i = 0; i < 4; i++) {


        auto steer_rod = chrono_types::make_shared<ChBodyEasyBox>(0.1, 0.1, 0.1, 1000, true, false);
        steer_rod->SetPos(m_wheels[i]->GetPos());
        steer_rod->SetBodyFixed(false);
        m_system->Add(steer_rod);

        ChQuaternion<> z2y;
        z2y.Q_from_AngAxis(CH_C_PI / 2, ChVector<>(1, 0, 0));

        m_drive_motor_funcs[i] = chrono_types::make_shared<ChFunction_Const>();
        m_drive_motors[i] =
        AddMotorSpeed(m_chassis->GetBody(), m_wheels[i]->GetBody(), m_chassis, wheel_rel_pos[i], z2y);
        m_drive_motors[i]->SetMotorFunction(m_drive_motor_funcs[i]);


    }

}

void Iris::SetDriver(std::shared_ptr<IrisDriver> driver) {
    m_driver = driver;
    m_driver->iris = this;
}

void Iris::SetWheelContactMaterial(std::shared_ptr<ChMaterialSurface> mat) {
    for (auto& wheel : m_wheels)
        wheel->m_mat = mat;
}

void Iris::SetChassisFixed(bool fixed) {
    m_chassis_fixed = fixed;
}

void Iris::SetChassisVisualization(bool state) {
    m_chassis->SetVisualize(state);
}

void Iris::SetWheelVisualization(bool state) {
    for (auto& wheel : m_wheels)
        wheel->SetVisualize(state);
}

ChVector<> Iris::GetWheelContactForce(IrisWheelID id) const {
    return m_wheels[id]->GetBody()->GetContactForce();
}

ChVector<> Iris::GetWheelContactTorque(IrisWheelID id) const {
    return m_wheels[id]->GetBody()->GetContactTorque();
}

ChVector<> Iris::GetWheelAppliedForce(IrisWheelID id) const {
    return m_wheels[id]->GetBody()->GetAppliedForce();
}

ChVector<> Iris::GetWheelAppliedTorque(IrisWheelID id) const {
    return m_wheels[id]->GetBody()->GetAppliedTorque();
}

double Iris::GetWheelTracTorque(IrisWheelID id) const {
    if (m_driver->GetDriveMotorType() == IrisDriver::DriveMotorType::TORQUE)
        return 0;

    return m_drive_motors[id]->GetMotorTorque();
}

double Iris::GetRoverMass() const {
    double tot_mass = m_chassis->GetBody()->GetMass();
    for (int i = 0; i < 4; i++) {
        tot_mass += m_wheels[i]->GetBody()->GetMass();
    }
    return tot_mass;
}

double Iris::GetWheelMass() const {
    return m_wheels[0]->GetBody()->GetMass();
}

void Iris::Update() {
    double time = m_system->GetChTime();
    m_driver->Update(time);

    for (int i = 0; i < 4; i++) {
        // Extract driver inputs
        double driving = m_driver->drive_speeds[i];


        // Set motor functions
        m_drive_motor_funcs[i]->SetSetpoint(driving, time);
    }
}

// =============================================================================

IrisDriver::IrisDriver()
    : drive_speeds({0, 0, 0, 0}), iris(nullptr) {}


IrisSpeedDriver::IrisSpeedDriver(double time_ramp, double speed) : m_ramp(time_ramp), m_speed(speed) {}

void IrisSpeedDriver::Update(double time) {
    double speed = m_speed;
    if (time < m_ramp)
        speed = m_speed * (time / m_ramp);
    drive_speeds = {speed, speed, speed, speed};
}

}  // namespace iris
}  // namespace chrono
