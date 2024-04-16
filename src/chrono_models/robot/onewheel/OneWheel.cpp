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
// Authors: Jason Zhou, Radu Serban, Sidney Nimako
// =============================================================================
//
// NASA VIPER Lunar Rover Model Class.
// This class contains model for NASA's VIPER lunar rover for NASA's 2024 Moon
// exploration mission.
//
// =============================================================================
//
// RADU TODO:
// - Recheck kinematics of mechanism (for negative lift angle)
// - Forces and torques are reported relative to the part's centroidal frame.
//   Likely confusing for a user since all bodies are ChBodyAuxRef!
// - Consider using a torque motor instead of driveshafts
//   (for a driver that uses torque control)
//
// =============================================================================

#include <cmath>

#include <iostream>

#include "chrono/assets/ChVisualShapeBox.h"
#include "chrono/physics/ChBodyEasy.h"
#include "chrono/assets/ChVisualShapeCylinder.h"
#include "chrono/assets/ChVisualShapeSphere.h"
#include "chrono/assets/ChTexture.h"
#include "chrono/assets/ChVisualShapeTriangleMesh.h"

#include "chrono/motion_functions/ChFunction_Setpoint.h"

#include "chrono/physics/ChLinkMotorRotationAngle.h"
#include "chrono/physics/ChLinkMotorRotationSpeed.h"
#include "chrono/physics/ChLinkMotorRotationTorque.h"
#include "chrono/physics/ChShaftsBody.h"

#include "chrono/physics/ChInertiaUtils.h"

#include "chrono_models/robot/onewheel/onewheel.h"

namespace chrono {
namespace onewheel {


// Default contact material for rover parts
std::shared_ptr<ChMaterialSurface> DefaultContactMaterial(ChContactMethod contact_method, SkidSteerParameters m_params) {
    float mu = m_params.wheel_mu;  // coefficient of friction
    float cr = m_params.cr;  // coefficient of restitution
    float Y = m_params.Y;    // Young's modulus
    float nu = m_params.nu;  // Poisson ratio
    float kn = m_params.kn;  // normal stiffness
    float gn = m_params.gn;  // normal viscous damping
    float kt = m_params.kt;  // tangential stiffness
    float gt = m_params.gt;  // tangential viscous damping

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
                                                        std::shared_ptr<OneWheelChassis> chassis,
                                                        const ChVector<>& rel_pos,
                                                        const ChQuaternion<>& rel_rot) {
    // Express relative frame in global
    ChFrame<> X_GC = chassis->GetBody()->GetFrame_REF_to_abs() * ChFrame<>(rel_pos, rel_rot);

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
                                                        std::shared_ptr<OneWheelChassis> chassis,
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
                                                          std::shared_ptr<OneWheelChassis> chassis,
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

// Base class for all OneWheel Part
OneWheelPart::OneWheelPart(const std::string& name,
                             const ChFrame<>& rel_pos,
                             std::shared_ptr<ChMaterialSurface> mat,
                             bool collide)
    : m_name(name), m_pos(rel_pos), m_mat(mat), m_collide(collide), m_visualize(true) {}

void OneWheelPart::Construct(ChSystem* system) {
    m_body = chrono_types::make_shared<ChBodyAuxRef>();
    m_body->SetNameString(m_name + "_body");
    m_body->SetMass(m_mass);
    m_body->SetInertiaXX(m_inertia);
    m_body->SetFrame_COG_to_REF(m_cog);

    // Add visualization shape
    if (m_visualize) {
        //TODO: Fully parametrize
        auto vis_mesh_file = m_mesh_name;
        auto trimesh_vis = geometry::ChTriangleMeshConnected::CreateFromWavefrontFile(vis_mesh_file, true, true);

        // scale mesh
        trimesh_vis->Transform(m_mesh_xform.GetPos(), m_mesh_xform.GetA());  // translate/rotate/scale mesh
        trimesh_vis->RepairDuplicateVertexes(1e-9);                          // if meshes are not watertight

        auto trimesh_shape = chrono_types::make_shared<ChVisualShapeTriangleMesh>();
        trimesh_shape->SetMesh(trimesh_vis);
        trimesh_shape->SetName(m_mesh_name);
        trimesh_shape->SetMutable(false);
        trimesh_shape->SetColor(m_color);
        m_body->AddVisualShape(trimesh_shape);

    }
    // Add collision shape
    if (m_collide) {
        auto col_mesh_file = m_mesh_name;

        auto trimesh_col = geometry::ChTriangleMeshConnected::CreateFromWavefrontFile(col_mesh_file, false, false);

        trimesh_col->Transform(m_mesh_xform.GetPos(), m_mesh_xform.GetA());  // translate/rotate/scale mesh
        trimesh_col->RepairDuplicateVertexes(1e-9);                          // if meshes are not watertight

        auto shape = chrono_types::make_shared<ChCollisionShapeTriangleMesh>(m_mat, trimesh_col, false, false, 0.005);
        m_body->AddCollisionShape(shape);
        m_body->SetCollide(m_collide);
    }

    system->AddBody(m_body);
}

void OneWheelPart::CalcMassProperties(double density) {
    auto mesh_filename = m_mesh_name;
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

void OneWheelPart::Initialize(std::shared_ptr<ChBodyAuxRef> chassis) {
    Construct(chassis->GetSystem());

    // Set absolute position
    ChFrame<> X_GC = chassis->GetFrame_REF_to_abs() * m_pos;
    m_body->SetFrame_REF_to_abs(X_GC);
}

// =============================================================================

// Rover Chassis
OneWheelChassis::OneWheelChassis(const std::string& name, 
                                    std::shared_ptr<ChMaterialSurface> mat, 
                                    std::string mesh_name, double mass)
    : OneWheelPart(name, ChFrame<>(VNULL, QUNIT), mat, false) {
    m_mesh_name = mesh_name;
    m_mass = mass;
    m_color = ChColor(1.0f, 1.0f, 1.0f);

    m_inertia = ChVector<>(1e-2, 0.014, 0.015);  // TODO: ask heather what to put for inertia?

    m_visualize = false;
    m_collide = false;
}

void OneWheelChassis::Initialize(ChSystem* system, const ChFrame<>& pos) {
    Construct(system);

    m_body->SetFrame_REF_to_abs(pos);
}

// =============================================================================

// OneWheel Wheel
OneWheelWheel::OneWheelWheel(const std::string& name,
                               const ChFrame<>& rel_pos,
                               std::shared_ptr<ChMaterialSurface> mat,
                               OneWheelWheelType wheel_type,
                               std::string wheel_mesh_name,
                               double mass
                               )
    : OneWheelPart(name, rel_pos, mat, true) {
    switch (wheel_type) {
        case OneWheelWheelType::RealWheel:
            m_mesh_name = wheel_mesh_name;
            break;
        case OneWheelWheelType::SimpleWheel:
            m_mesh_name = wheel_mesh_name;
            break;
        case OneWheelWheelType::CylWheel:
            m_mesh_name = wheel_mesh_name;
            break;
    }

    m_color = ChColor(0.4f, 0.7f, 0.4f);
    m_mass = mass;                                 // weight of the wheel
    m_inertia = ChVector<double>(8.74e-4, 8.77e-4, 16.81e-4);  // principal inertia

    ChMatrix33<> A;
    A.setZero();
    A(0, 0) = -0.243;
    A(2, 0) = 0.97;
    A(0, 1) = 0.97;
    A(2, 1) = 0.243;
    A(1, 2) = 1;

    m_cog = ChFrame<>(ChVector<>(0, 0.00426, 0), A.Get_A_quaternion());
}

// =============================================================================

// Rover model
OneWheel::OneWheel(ChSystem* system, OneWheelWheelType wheel_type, std::string fp)
    : m_system(system), m_chassis_fixed(false) {
    // Set default collision model envelope commensurate with model dimensions.
    // Note that an SMC system automatically sets envelope to 0.
    auto contact_method = m_system->GetContactMethod();
    if (contact_method == ChContactMethod::NSC) {
        ChCollisionModel::SetDefaultSuggestedEnvelope(0.01);
        ChCollisionModel::SetDefaultSuggestedMargin(0.005);
    }

    // Create the contact materials
    m_default_material = DefaultContactMaterial(contact_method, m_params);
    m_wheel_material = DefaultContactMaterial(contact_method, m_params);

    char* f = new char[fp.length() + 1];
    strcpy(f, fp.c_str());

    ChStreamInAsciiFile mfilei(f);
    ChArchiveInJSON marchivei(mfilei);
    SkidSteerParameters skid_steer_params;
    marchivei >> CHNVP(skid_steer_params);
    m_params = skid_steer_params;

    Create(wheel_type);
}


void OneWheel::Create(OneWheelWheelType wheel_type) {
    // initilize rover wheels
    const double wheel_x = m_params.wheel_x;
    const double wheel_y = m_params.wheel_y;
    const double wheel_z = m_params.wheel_z;
    std::cout << "Making chassis" << std::endl;
    std::cout << m_params.chassis_mass << std::endl;
    
    // create rover chassis
    m_chassis = chrono_types::make_shared<OneWheelChassis>("chassis", m_default_material, m_params.chassis_mesh_name, m_params.chassis_mass);

    m_wheel chrono_types::make_shared<OneWheelWheel>(
        "wheel", ChFrame<>(ChVector<>(+wheel_x, +wheel_y, wheel_z), QUNIT), m_wheel_material, wheel_type, m_params.wheel_mesh_name, m_params.wheel_mass);
}

void OneWheel::Initialize(const ChFrame<>& pos) {
    double wheel_x = m_params.wheel_x;
    double wheel_y = m_params.wheel_y;
    double wheel_z = m_params.wheel_z;

    assert(m_driver);
   
    m_chassis->Initialize(m_system, pos);
  
    m_chassis->GetBody()->SetBodyFixed(m_chassis_fixed);
 

    m_wheel->Initialize(m_chassis->GetBody());

    ChVector<> wheel_rel_pos = ChVector<>(+wheel_x, +wheel_y, wheel_z);

    ChQuaternion<> z2x = Q_from_AngY(CH_C_PI_2);

    ChQuaternion<> z2y;
    z2y.Q_from_AngAxis(CH_C_PI / 2, ChVector<>(1, 0, 0));
    m_drive_motor_funcs = chrono_types::make_shared<ChFunction_Setpoint>();
    m_drive_motors = AddMotorSpeed(m_chassis->GetBody(), m_wheel->GetBody(), m_chassis, wheel_rel_pos[i], z2y);
    m_drive_motors->SetMotorFunction(m_drive_motor_funcs[i]);
}

void OneWheel::SetWheelContactMaterial(std::shared_ptr<ChMaterialSurface> mat) {
    m_wheel->m_mat = mat;
}

void OneWheel::SetChassisFixed(bool fixed) {
    m_chassis_fixed = fixed;
}

void OneWheel::SetChassisVisualization(bool state) {
    m_chassis->SetVisualize(state);
}

void OneWheel::SetWheelVisualization(bool state) {
    m_wheel->SetVisualize(state);
}

ChVector<> OneWheel::GetWheelContactForce() const {
    return m_wheel>GetBody()->GetContactForce();
}

ChVector<> OneWheel::GetWheelContactTorque() const {
    return m_wheel->GetBody()->GetContactTorque();
}

ChVector<> OneWheel::GetWheelAppliedForce() const {
    return m_wheel->GetBody()->GetAppliedForce();
}

ChVector<> OneWheel::GetWheelAppliedTorque() const {
    return m_wheel->GetBody()->GetAppliedTorque();
}

double OneWheel::GetWheelTracTorque() const {
    return m_drive_motors->GetMotorTorque();
}

double OneWheel::GetRoverMass() const {
    double tot_mass = m_chassis->GetBody()->GetMass();
    m_wheel->GetBody()->GetMass();
    return tot_mass;
}

double OneWheel::GetWheelMass() const {
    return m_wheel->GetBody()->GetMass();
}

void OneWheel::Update() {
    double time = m_system->GetChTime();
    m_driver->Update(time);

    double driving = m_driver->drive_speeds[i];

    m_drive_motor_funcs->SetSetpoint(driving, time);
}
// =============================================================================

OneWheelSpeedDriver::OneWheelSpeedDriver(double time_ramp, double speed) : m_ramp(time_ramp), m_speed(speed) {}

void OneWheelSpeedDriver::Update(double time) {
    double speed = m_speed;
    if (time < m_ramp)
        speed = m_speed * (time / m_ramp);
    drive_speeds = {speed, speed, speed, speed};
}

}  // namespace OneWheel
}  // namespace chrono
