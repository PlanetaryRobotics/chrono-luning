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
// Authors: Jason Zhou, Radu Serban
// =============================================================================
//
// NASA VIPER Lunar Rover Model Class.
// This class contains model for NASA's VIPER lunar rover for NASA's 2024 Moon
// exploration mission.
//
// =============================================================================

#ifndef IRIS_H
#define IRIS_H

#include <string>
#include <array>

#include "chrono/assets/ChColor.h"
#include "chrono/physics/ChLinkMotorRotation.h"
#include "chrono/physics/ChSystem.h"
#include "chrono/physics/ChShaft.h"

#include "chrono_models/ChApiModels.h"

namespace chrono {

/// Namespace with classes for the Iris model.
namespace iris {

/// @addtogroup robot_models_iris
/// @{

/// Iris wheel/suspension identifiers.
enum IrisWheelID {
    LF = 0,  ///< left front
    RF = 1,  ///< right front
    LB = 2,  ///< left back
    RB = 3   ///< right back
};

/// Iris wheel type.
enum class IrisWheelType {
    RealWheel,    ///< actual geometry of the iris wheel
    SimpleWheel,  ///< simplified wheel geometry
    CylWheel      ///< cylindrical wheel geometry
};

// -----------------------------------------------------------------------------

/// Base class definition for all iris parts.
/// iris Rover Parts include Chassis, Steering, Upper Suspension Arm, Bottom Suspension Arm and Wheel.
class CH_MODELS_API IrisPart {
  public:
    IrisPart(const std::string& name,                 ///< part name
              const ChFrame<>& rel_pos,                ///< position relative to chassis frame
              std::shared_ptr<ChMaterialSurface> mat,  ///< contact material
              bool collide                             ///< enable collision?
    );
    virtual ~IrisPart() {}

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

    /// Return the ChBody of the corresponding Iris part.
    std::shared_ptr<ChBodyAuxRef> GetBody() const { return m_body; }

    /// Return the position of the Iris part.
    /// This is the absolute location of the part reference frame.
    const ChVector<>& GetPos() const { return m_body->GetFrame_REF_to_abs().GetPos(); }

    /// Return the rotation of the Iris part.
    /// This is the orientation wrt the global frame of the part reference frame.
    const ChQuaternion<>& GetRot() const { return m_body->GetFrame_REF_to_abs().GetRot(); }

    /// Return the linear velocity of the Iris part.
    /// This is the absolute linear velocity of the part reference frame.
    const ChVector<>& GetLinVel() const { return m_body->GetFrame_REF_to_abs().GetPos_dt(); }

    /// Return the angular velocity of the Iris part.
    /// This is the absolute angular velocity of the part reference frame.
    const ChVector<> GetAngVel() const { return m_body->GetFrame_REF_to_abs().GetWvel_par(); }

    /// Return the linear acceleration of the Iris part.
    /// This is the absolute linear acceleration of the part reference frame.
    const ChVector<>& GetLinAcc() const { return m_body->GetFrame_REF_to_abs().GetPos_dtdt(); }

    /// Return the angular acceleratino of the Iris part.
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

/// Iris rover Chassis.
class CH_MODELS_API IrisChassis : public IrisPart {
  public:
    IrisChassis(const std::string& name,                ///< part name
                 std::shared_ptr<ChMaterialSurface> mat  ///< contact material
    );
    ~IrisChassis() {}

    /// Initialize the chassis at the specified (absolute) position.
    void Initialize(ChSystem* system, const ChFrame<>& pos);
};

/// Iris rover Wheel.
class CH_MODELS_API IrisWheel : public IrisPart {
  public:
    IrisWheel(const std::string& name,                 ///< part name
               const ChFrame<>& rel_pos,                ///< position relative to chassis frame
               std::shared_ptr<ChMaterialSurface> mat,  ///< contact material
               IrisWheelType wheel_type                ///< wheel type
    );
    ~IrisWheel() {}

    friend class Iris;
};


// -----------------------------------------------------------------------------
class CH_MODELS_API IrisSpeedDriver {
  public:
    IrisSpeedDriver(double time_ramp, double speed);
    ~IrisSpeedDriver() {}

    void Update(double time);

    double m_ramp;
    double m_speed;

    Iris* iris;  ///< associated Iris rover

    std::array<double, 4> drive_speeds;  ///< angular speeds for drive motors

    friend class Iris;
};

/// Iris rover class.
/// This class encapsulates the location and rotation information of all Iris parts wrt the chassis.
/// This class should be the entry point to create a complete rover.
class CH_MODELS_API Iris {
  public:
    Iris(ChSystem* system, IrisWheelType wheel_type = IrisWheelType::RealWheel);

    ~Iris() {}

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

    void SetSpeedDriver(std::shared_ptr<IrisSpeedDriver> driver) {
        m_driver = driver;
        m_driver->iris = this;
    }

    /// Initialize the Iris rover at the specified position.
    void Initialize(const ChFrame<>& pos);

    /// Get the rover chassis.
    std::shared_ptr<IrisChassis> GetChassis() const { return m_chassis; }

    /// Get all rover wheels.
    std::array<std::shared_ptr<IrisWheel>, 4> GetWheels() const { return m_wheels; }

    /// Get the specified rover wheel.
    std::shared_ptr<IrisWheel> GetWheel(IrisWheelID id) const { return m_wheels[id]; }


    /// Get chassis position.
    ChVector<> GetChassisPos() const { return m_chassis->GetPos(); }

    /// Get chassis orientation.
    ChQuaternion<> GetChassisRot() const { return m_chassis->GetRot(); }

    /// Get chassis linear velocity.
    ChVector<> GetChassisVel() const { return m_chassis->GetLinVel(); }

    /// Get chassis linear acceleration.
    ChVector<> GetChassisAcc() const { return m_chassis->GetLinAcc(); }

    /// Get wheel speed.
    ChVector<> GetWheelLinVel(IrisWheelID id) const { return m_wheels[id]->GetLinVel(); }

    /// Get wheel angular velocity.
    ChVector<> GetWheelAngVel(IrisWheelID id) const { return m_wheels[id]->GetAngVel(); }

    /// Get wheel contact force.
    ChVector<> GetWheelContactForce(IrisWheelID id) const;

    /// Get wheel contact torque.
    ChVector<> GetWheelContactTorque(IrisWheelID id) const;

    /// Get wheel total applied force.
    ChVector<> GetWheelAppliedForce(IrisWheelID id) const;

    /// Get wheel tractive torque - if DC control set to off
    double GetWheelTracTorque(IrisWheelID id) const;

    /// Get wheel total applied torque.
    ChVector<> GetWheelAppliedTorque(IrisWheelID id) const;

    /// Get total rover mass.
    double GetRoverMass() const;

    /// Get total wheel mass.
    double GetWheelMass() const;

    /// Get drive motor function.
    /// This will return an empty pointer if the associated driver uses torque control.
    std::shared_ptr<ChFunction_Setpoint> GetDriveMotorFunc(IrisWheelID id) const { return m_drive_motor_funcs[id]; }


    /// Get drive motor.
    /// This will return an empty pointer if the associated driver uses torque control.
    std::shared_ptr<ChLinkMotorRotation> GetDriveMotor(IrisWheelID id) const { return m_drive_motors[id]; }


    /// Iris update function.
    /// This function must be called before each integration step.
    void Update();

  private:
    /// Create the rover parts.
    void Create(IrisWheelType wheel_type);

    ChSystem* m_system;  ///< pointer to the Chrono system

    bool m_chassis_fixed;  ///< fix chassis to ground

    std::shared_ptr<IrisChassis> m_chassis;                     ///< rover chassis
    std::array<std::shared_ptr<IrisWheel>, 4> m_wheels;         ///< rover wheels (LF, RF, LR, RB)

    std::array<std::shared_ptr<ChLinkMotorRotation>, 4> m_drive_motors;  ///< drive motors

    std::array<std::shared_ptr<ChFunction_Setpoint>, 4> m_drive_motor_funcs;  ///< drive motor functions


    std::shared_ptr<IrisSpeedDriver> m_driver;  ///< rover driver system

    std::shared_ptr<ChMaterialSurface> m_default_material;  ///< common contact material for all non-wheel parts
    std::shared_ptr<ChMaterialSurface> m_wheel_material;    ///< wheel contact material (shared across limbs)

    static const double m_max_steer_angle;  ///< maximum steering angle

};

// -----------------------------------------------------------------------------



/// Concrete Iris speed driver.
/// This driver applies the same angular speed (ramped from 0 to a prescribed value) to all wheels.


/// @} robot_models_iris

}  // namespace iris
}  // namespace chrono

#endif
