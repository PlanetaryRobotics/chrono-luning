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
// Authors: Sidney Nimako
// =============================================================================
//

#ifndef IRIS_H
#define IRIS_H

#include <string>
#include <array>

#include "chrono/assets/ChColor.h"
#include "chrono/physics/ChLinkMotorRotation.h"
#include "chrono/physics/ChSystem.h"
#include "chrono/physics/ChShaft.h"

#include "chrono_models/ChApiModels.h"

#include "chrono/serialization/ChArchive.h"
#include "chrono/serialization/ChArchiveJSON.h"

namespace chrono {

/// Namespace with classes for the SkidSteer model.
namespace skidsteer {

/// @addtogroup robot_models_skidsteer
/// @{

/// SkidSteer wheel/suspension identifiers.
enum SkidSteerWheelID {
    LF = 0,  ///< left front
    RF = 1,  ///< right front
    LB = 2,  ///< left back
    RB = 3   ///< right back
};

/// SkidSteer wheel type.
enum class SkidSteerWheelType {
    RealWheel,    ///< actual geometry of the skidsteer wheel
    SimpleWheel,  ///< simplified wheel geometry
    CylWheel      ///< cylindrical wheel geometry
};

// -----------------------------------------------------------------------------
/// SkidSteer rover parameters structure.
class SkidSteerParameters {
  public:
    // Geometry
    double wheel_x;
    double wheel_y;
    double wheel_z;
    double chassis_dim_x;
    double chassis_dim_y;
    double chassis_dim_z;
    // Contact material
    float mu;
    float cr;
    float Y;
    float nu;
    float kn;
    float gn;
    float kt;
    float gt;
    // Mesh
    std::string chassis_mesh;
    std::string wheel_mesh;
    // Inertia
    double m_chassis;
    double i_chassis;
    double m_wheel;
    double i_wheel;

    SkidSteerParameters() {}
    SkidSteerParameters(double wheel_x,
                        double wheel_y,
                        double wheel_z,
                        double chassis_dim_x,
                        double chassis_dim_y,
                        double chassis_dim_z,
                        float mu,
                        float cr,
                        float Y,
                        float nu,
                        float kn,
                        float gn,
                        float kt,
                        float gt,
                        std::string chassis_mesh,
                        std::string wheel_mesh,
                        double m_chassis,
                        double i_chassis,
                        double m_wheel,
                        double i_wheel)
        : wheel_x(wheel_x),
          wheel_y(wheel_y),
          wheel_z(wheel_z),
          chassis_dim_x(chassis_dim_x),
          chassis_dim_y(chassis_dim_y),
          chassis_dim_z(chassis_dim_z),
          mu(mu),
          cr(cr),
          Y(Y),
          nu(nu),
          kn(kn),
          gn(gn),
          kt(kt),
          gt(gt),
          chassis_mesh(chassis_mesh),
          wheel_mesh(wheel_mesh),
          m_chassis(m_chassis),
          i_chassis(i_chassis),
          m_wheel(m_wheel),
          i_wheel(i_wheel) {}

    virtual ~SkidSteerParameters() {}

    virtual void ArchiveOut(ChArchiveOut& marchive) {
        marchive << CHNVP(wheel_x);
        marchive << CHNVP(wheel_y);
        marchive << CHNVP(wheel_z);
        marchive << CHNVP(chassis_dim_x);
        marchive << CHNVP(chassis_dim_y);
        marchive << CHNVP(chassis_dim_z);
        marchive << CHNVP(mu);
        marchive << CHNVP(cr);
        marchive << CHNVP(Y);
        marchive << CHNVP(nu);
        marchive << CHNVP(kn);
        marchive << CHNVP(gn);
        marchive << CHNVP(kt);
        marchive << CHNVP(gt);
        marchive << CHNVP(chassis_mesh);
        marchive << CHNVP(wheel_mesh);
        marchive << CHNVP(m_chassis);
        marchive << CHNVP(i_chassis);
        marchive << CHNVP(m_wheel);
        marchive << CHNVP(i_wheel);
    }

    virtual void ArchiveIn(ChArchiveIn& marchive) {
        marchive >> CHNVP(wheel_x);
        marchive >> CHNVP(wheel_y);
        marchive >> CHNVP(wheel_z);
        marchive >> CHNVP(chassis_dim_x);
        marchive >> CHNVP(chassis_dim_y);
        marchive >> CHNVP(chassis_dim_z);
        marchive >> CHNVP(mu);
        marchive >> CHNVP(cr);
        marchive >> CHNVP(Y);
        marchive >> CHNVP(nu);
        marchive >> CHNVP(kn);
        marchive >> CHNVP(gn);
        marchive >> CHNVP(kt);
        marchive >> CHNVP(gt);
        marchive >> CHNVP(chassis_mesh);
        marchive >> CHNVP(wheel_mesh);
        marchive >> CHNVP(m_chassis);
        marchive >> CHNVP(i_chassis);
        marchive >> CHNVP(m_wheel);
        marchive >> CHNVP(i_wheel);
    }
};

static SkidSteerParameters skidsteer_params;
// -----------------------------------------------------------------------------
/// Base class definition for all skidsteer parts.
/// skidsteer Rover Parts include Chassis, Steering, Upper Suspension Arm, Bottom Suspension Arm and Wheel.
class CH_MODELS_API SkidSteerPart {
  public:
    SkidSteerPart(const std::string& name,                 ///< part name
                  const ChFrame<>& rel_pos,                ///< position relative to chassis frame
                  std::shared_ptr<ChMaterialSurface> mat,  ///< contact material
                  bool collide                             ///< enable collision?
    );
    virtual ~SkidSteerPart() {}

    /// Return the name of the part.
    const std::string& GetName() const { return m_name; }

    /// Set the name of the part.
    void SetName(const std::string& name) { m_name = name; }

    /// Enable/disable visualization.
    void SetVisualize(bool state) { m_visualize = state; }

    /// Enable/disable collision.
    void SetCollide(bool state) { m_collide = state; }

    /// Initialize the rover part by attaching it to the specified chassis body.
    void Initialize(std::shared_ptr<ChBodyAuxRef> chassis);

    /// Return the ChBody of the corresponding SkidSteer part.
    std::shared_ptr<ChBodyAuxRef> GetBody() const { return m_body; }

    /// Return the position of the SkidSteer part.
    /// This is the absolute location of the part reference frame.
    const ChVector<>& GetPos() const { return m_body->GetFrame_REF_to_abs().GetPos(); }

    /// Return the rotation of the SkidSteer part.
    /// This is the orientation wrt the global frame of the part reference frame.
    const ChQuaternion<>& GetRot() const { return m_body->GetFrame_REF_to_abs().GetRot(); }

    /// Return the linear velocity of the SkidSteer part.
    /// This is the absolute linear velocity of the part reference frame.
    const ChVector<>& GetLinVel() const { return m_body->GetFrame_REF_to_abs().GetPos_dt(); }

    /// Return the angular velocity of the SkidSteer part.
    /// This is the absolute angular velocity of the part reference frame.
    const ChVector<> GetAngVel() const { return m_body->GetFrame_REF_to_abs().GetWvel_par(); }

    /// Return the linear acceleration of the SkidSteer part.
    /// This is the absolute linear acceleration of the part reference frame.
    const ChVector<>& GetLinAcc() const { return m_body->GetFrame_REF_to_abs().GetPos_dtdt(); }

    /// Return the angular acceleratino of the SkidSteer part.
    /// This is the absolute angular acceleratin of the part reference frame.
    const ChVector<> GetAngAcc() const { return m_body->GetFrame_REF_to_abs().GetWacc_par(); }

  protected:
    /// Utility function for calculating mass properties using the part's collision mesh.
    void CalcMassProperties(double density);

    /// Construct the part body.
    void Construct(ChSystem* system);

    std::string m_name;                        ///< part name
    std::shared_ptr<ChBodyAuxRef> m_body;      ///< part rigid body
    std::shared_ptr<ChMaterialSurface> m_mat;  ///< contact material

    std::string m_mesh_name;  ///< visualization mesh name
    ChFrame<> m_mesh_xform;   ///< mesh transform (translate, rotate, scale)
    ChColor m_color;          ///< visualization asset color

    ChFrame<> m_pos;       ///< relative position wrt the chassis
    double m_mass;         ///< mass
    ChVector<> m_inertia;  ///< principal moments of inertia
    ChFrame<> m_cog;       ///< COG frame (relative to body frame)

    bool m_visualize;  ///< part visualization flag
    bool m_collide;    ///< part collision flag
};

/// SkidSteer rover Chassis.
class CH_MODELS_API SkidSteerChassis : public SkidSteerPart {
  public:
    SkidSteerChassis(const std::string& name,                ///< part name
                     std::shared_ptr<ChMaterialSurface> mat  ///< contact material
    );
    ~SkidSteerChassis() {}

    /// Initialize the chassis at the specified (absolute) position.
    void Initialize(ChSystem* system, const ChFrame<>& pos);
};

/// SkidSteer rover Wheel.
class CH_MODELS_API SkidSteerWheel : public SkidSteerPart {
  public:
    SkidSteerWheel(const std::string& name,                 ///< part name
                   const ChFrame<>& rel_pos,                ///< position relative to chassis frame
                   std::shared_ptr<ChMaterialSurface> mat,  ///< contact material
                   SkidSteerWheelType wheel_type            ///< wheel type
    );
    ~SkidSteerWheel() {}

    friend class SkidSteer;
};

class SkidSteer;

// -----------------------------------------------------------------------------
class CH_MODELS_API SkidSteerSpeedDriver {
  public:
    SkidSteerSpeedDriver(double time_ramp, double speed);
    ~SkidSteerSpeedDriver() {}

    void Update(double time);

    double m_ramp;
    double m_speed;

    SkidSteer* skidsteer;  ///< associated SkidSteer rover

    std::array<double, 4> drive_speeds;  ///< angular speeds for drive motors

    friend class SkidSteer;
};

/// SkidSteer rover class.
/// This class encapsulates the location and rotation information of all SkidSteer parts wrt the chassis.
/// This class should be the entry point to create a complete rover.
class CH_MODELS_API SkidSteer {
  public:
    SkidSteer(ChSystem* system, SkidSteerWheelType wheel_type = SkidSteerWheelType::RealWheel);

    ~SkidSteer() {}

    /// Get the containing system.
    ChSystem* GetSystem() const { return m_system; }

    /// Set wheel contact material.
    void SetWheelContactMaterial(std::shared_ptr<ChMaterialSurface> mat);

    /// Fix the chassis to ground.
    void SetChassisFixed(bool fixed);

    /// Enable/disable visualization of the rover chassis (default: true).
    void SetChassisVisualization(bool state);

    /// Enable/disable visualization of rover wheels (default: true).
    void SetWheelVisualization(bool state);

    void SetSpeedDriver(std::shared_ptr<SkidSteerSpeedDriver> driver) {
        m_driver = driver;
        m_driver->skidsteer = this;
    }

    /// Initialize the SkidSteer rover at the specified position.
    void Initialize(const ChFrame<>& pos);

    /// Get the rover chassis.
    std::shared_ptr<SkidSteerChassis> GetChassis() const { return m_chassis; }

    /// Get all rover wheels.
    std::array<std::shared_ptr<SkidSteerWheel>, 4> GetWheels() const { return m_wheels; }

    /// Get the specified rover wheel.
    std::shared_ptr<SkidSteerWheel> GetWheel(SkidSteerWheelID id) const { return m_wheels[id]; }

    /// Get chassis position.
    ChVector<> GetChassisPos() const { return m_chassis->GetPos(); }

    /// Get chassis orientation.
    ChQuaternion<> GetChassisRot() const { return m_chassis->GetRot(); }

    /// Get chassis linear velocity.
    ChVector<> GetChassisVel() const { return m_chassis->GetLinVel(); }

    /// Get chassis linear acceleration.
    ChVector<> GetChassisAcc() const { return m_chassis->GetLinAcc(); }

    /// Get wheel speed.
    ChVector<> GetWheelLinVel(SkidSteerWheelID id) const { return m_wheels[id]->GetLinVel(); }

    /// Get wheel angular velocity.
    ChVector<> GetWheelAngVel(SkidSteerWheelID id) const { return m_wheels[id]->GetAngVel(); }

    /// Get wheel contact force.
    ChVector<> GetWheelContactForce(SkidSteerWheelID id) const;

    /// Get wheel contact torque.
    ChVector<> GetWheelContactTorque(SkidSteerWheelID id) const;

    /// Get wheel total applied force.
    ChVector<> GetWheelAppliedForce(SkidSteerWheelID id) const;

    /// Get wheel tractive torque - if DC control set to off
    double GetWheelTracTorque(SkidSteerWheelID id) const;

    /// Get wheel total applied torque.
    ChVector<> GetWheelAppliedTorque(SkidSteerWheelID id) const;

    /// Get total rover mass.
    double GetRoverMass() const;

    /// Get total wheel mass.
    double GetWheelMass() const;

    /// Get drive motor function.
    /// This will return an empty pointer if the associated driver uses torque control.
    std::shared_ptr<ChFunction_Setpoint> GetDriveMotorFunc(SkidSteerWheelID id) const {
        return m_drive_motor_funcs[id];
    }

    /// Get drive motor.
    /// This will return an empty pointer if the associated driver uses torque control.
    std::shared_ptr<ChLinkMotorRotation> GetDriveMotor(SkidSteerWheelID id) const { return m_drive_motors[id]; }

    /// SkidSteer update function.
    /// This function must be called before each integration step.
    void Update();

  private:
    /// Create the rover parts.
    void Create(SkidSteerWheelType wheel_type);

    ChSystem* m_system;  ///< pointer to the Chrono system

    bool m_chassis_fixed;  ///< fix chassis to ground

    std::shared_ptr<SkidSteerChassis> m_chassis;              ///< rover chassis
    std::array<std::shared_ptr<SkidSteerWheel>, 4> m_wheels;  ///< rover wheels (LF, RF, LR, RB)

    std::array<std::shared_ptr<ChLinkMotorRotation>, 4> m_drive_motors;  ///< drive motors

    std::array<std::shared_ptr<ChFunction_Setpoint>, 4> m_drive_motor_funcs;  ///< drive motor functions

    std::shared_ptr<SkidSteerSpeedDriver> m_driver;  ///< rover driver system

    std::shared_ptr<ChMaterialSurface> m_default_material;  ///< common contact material for all non-wheel parts
    std::shared_ptr<ChMaterialSurface> m_wheel_material;    ///< wheel contact material (shared across limbs)

    static const double m_max_steer_angle;  ///< maximum steering angle
};

// -----------------------------------------------------------------------------

/// Concrete SkidSteer speed driver.
/// This driver applies the same angular speed (ramped from 0 to a prescribed value) to all wheels.

/// @} robot_models_skidsteer

}  // namespace skidsteer
}  // namespace chrono

#endif
