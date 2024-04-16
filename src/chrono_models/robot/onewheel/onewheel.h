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

#pragma once

#include <string>
#include <array>

#include "chrono/assets/ChColor.h"
#include "chrono/physics/ChLinkMotorRotation.h"
#include "chrono/physics/ChSystem.h"
#include "chrono/physics/ChShaft.h"

#include "chrono_models/ChApiModels.h"

#include "chrono/serialization/ChArchive.h"
#include "chrono/serialization/ChArchiveJSON.h"

#include "chrono_models/robot/skidsteer/SkidSteerParameters.h"

namespace chrono {

/// Namespace with classes for the OneWheel model.
namespace onewheel {

/// @addtogroup robot_models_OneWheel
/// @{

/// OneWheel wheel/suspension identifiers.
// enum OneWheelWheelID {
//     LF = 0,  ///< left front
//     RF = 1,  ///< right front
//     LB = 2,  ///< left back
//     RB = 3   ///< right back
// };

/// OneWheel wheel type.
enum class OneWheelWheelType {
    RealWheel,    ///< actual geometry of the OneWheel wheel
    SimpleWheel,  ///< simplified wheel geometry
    CylWheel      ///< cylindrical wheel geometry
};
// -----------------------------------------------------------------------------
/// Base class definition for all OneWheel parts.
/// OneWheel Rover Parts include Chassis, Steering, Upper Suspension Arm, Bottom Suspension Arm and Wheel.
class CH_MODELS_API OneWheelPart {
  public:
    OneWheelPart(const std::string& name,                 ///< part name
                  const ChFrame<>& rel_pos,                ///< position relative to chassis frame
                  std::shared_ptr<ChMaterialSurface> mat,  ///< contact material
                  bool collide                             ///< enable collision?
    );
    virtual ~OneWheelPart() {}

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

    /// Return the ChBody of the corresponding OneWheel part.
    std::shared_ptr<ChBodyAuxRef> GetBody() const { return m_body; }

    /// Return the position of the OneWheel part.
    /// This is the absolute location of the part reference frame.
    const ChVector<>& GetPos() const { return m_body->GetFrame_REF_to_abs().GetPos(); }

    /// Return the rotation of the OneWheel part.
    /// This is the orientation wrt the global frame of the part reference frame.
    const ChQuaternion<>& GetRot() const { return m_body->GetFrame_REF_to_abs().GetRot(); }

    /// Return the linear velocity of the OneWheel part.
    /// This is the absolute linear velocity of the part reference frame.
    const ChVector<>& GetLinVel() const { return m_body->GetFrame_REF_to_abs().GetPos_dt(); }

    /// Return the angular velocity of the OneWheel part.
    /// This is the absolute angular velocity of the part reference frame.
    const ChVector<> GetAngVel() const { return m_body->GetFrame_REF_to_abs().GetWvel_par(); }

    /// Return the linear acceleration of the OneWheel part.
    /// This is the absolute linear acceleration of the part reference frame.
    const ChVector<>& GetLinAcc() const { return m_body->GetFrame_REF_to_abs().GetPos_dtdt(); }

    /// Return the angular acceleratino of the OneWheel part.
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

/// OneWheel rover Chassis.
class CH_MODELS_API OneWheelChassis : public OneWheelPart {
  public:
    OneWheelChassis(const std::string& name,                ///< part name
                     std::shared_ptr<ChMaterialSurface> mat,  ///< contact material
                     std::string mesh_name,                 /// < Mesh name
                     double mass                            /// < Mass of the chassis
    );
    ~OneWheelChassis() {}

    /// Initialize the chassis at the specified (absolute) position.
    void Initialize(ChSystem* system, const ChFrame<>& pos);
};

/// OneWheel rover Wheel.
class CH_MODELS_API OneWheelWheel : public OneWheelPart {
  public:
    OneWheelWheel(const std::string& name,                 ///< part name
                   const ChFrame<>& rel_pos,                ///< position relative to chassis frame
                   std::shared_ptr<ChMaterialSurface> mat,  ///< contact material
                   OneWheelWheelType wheel_type,           ///< wheel type
                   std::string wheel_mesh_name,             /// < Path to wheel mesh file
                   double mass                            /// < Mass of the wheel
    );
    ~OneWheelWheel() {}

    friend class OneWheel;
};

class OneWheel;

// -----------------------------------------------------------------------------
class CH_MODELS_API OneWheelSpeedDriver {
  public:
    OneWheelSpeedDriver(double time_ramp, double speed);
    ~OneWheelSpeedDriver() {}

    void Update(double time);

    double m_ramp;
    double m_speed;

    OneWheel* onewheel;  ///< associated OneWheel rover

    std::array<double> drive_speeds;  ///< angular speeds for drive motors

    friend class OneWheel;
};

/// OneWheel rover class.
/// This class encapsulates the location and rotation information of all OneWheel parts wrt the chassis.
/// This class should be the entry point to create a complete rover.
class CH_MODELS_API OneWheel {
  public:
    OneWheel(ChSystem* system,
              OneWheelWheelType wheel_type = OneWheelWheelType::RealWheel, std::string fp = ""); //TODO: Change order for required parameters

    ~OneWheel() {}

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

    void SetSpeedDriver(std::shared_ptr<OneWheelSpeedDriver> driver) {
        m_driver = driver;
        m_driver->onewheel = this;
    }

    /// Initialize the OneWheel rover at the specified position.
    void Initialize(const ChFrame<>& pos);

    /// Get the rover chassis.
    std::shared_ptr<OneWheelChassis> GetChassis() const { return m_chassis; }

    /// Get all rover wheels.
    std::shared_ptr<OneWheelWheel> GetWheels() const { return m_wheel; }

    /// Get wheel speed.
    ChVector<> GetWheelLinVel(OneWheelWheelID id) const { return m_wheel->GetLinVel(); }

    /// Get wheel angular velocity.
    ChVector<> GetWheelAngVel(OneWheelWheelID id) const { return m_wheel->GetAngVel(); }

    /// Get wheel contact force.
    ChVector<> GetWheelContactForce(OneWheelWheelID id) const;

    /// Get wheel contact torque.
    ChVector<> GetWheelContactTorque(OneWheelWheelID id) const;

    /// Get wheel total applied force.
    ChVector<> GetWheelAppliedForce(OneWheelWheelID id) const;

    /// Get wheel tractive torque - if DC control set to off
    double GetWheelTracTorque(OneWheelWheelID id) const;

    /// Get wheel total applied torque.
    ChVector<> GetWheelAppliedTorque(OneWheelWheelID id) const;

    /// Get total rover mass.
    double GetRoverMass() const;

    /// Get total wheel mass.
    double GetWheelMass() const;

    /// Get drive motor function.
    /// This will return an empty pointer if the associated driver uses torque control.
    std::shared_ptr<ChFunction_Setpoint> GetDriveMotorFunc(OneWheelWheelID id) const {
        return m_drive_motor_funcs[id];
    }

    /// Get drive motor.
    /// This will return an empty pointer if the associated driver uses torque control.
    std::shared_ptr<ChLinkMotorRotation> GetDriveMotor(OneWheelWheelID id) const { return m_drive_motors[id]; }

    /// OneWheel update function.
    /// This function must be called before each integration step.
    void Update();

    /// Update the rover parameters. DEPRECATED
    void UpdateParameters(const OneWheelParameters& new_params);

    /// Get the rover parameters. DEPRECATED
    const OneWheelParameters& GetParameters() const { return m_params; }


    SkidSteerParameters m_params;  ///< rover parameters

  private:
    /// Create the rover parts.
    void Create(OneWheelWheelType wheel_type);

    ChSystem* m_system;  ///< pointer to the Chrono system

    bool m_chassis_fixed;  ///< fix chassis to ground


    std::shared_ptr<OneWheelWheel> m_wheel;  
    std::shared_ptr<ChLinkMotorRotation> m_drive_motors; 

    std::shared_ptr<ChFunction_Setpoint> m_drive_motor_funcs;  ///< drive motor functions

    std::shared_ptr<OneWheelSpeedDriver> m_driver;  ///< rover driver system

    std::shared_ptr<ChMaterialSurface> m_default_material;  ///< common contact material for all non-wheel parts
    std::shared_ptr<ChMaterialSurface> m_wheel_material;    ///< wheel contact material (shared across limbs)

    static const double m_max_steer_angle;  ///< maximum steering angle
};

// -----------------------------------------------------------------------------

/// Concrete OneWheel speed driver.
/// This driver applies the same angular speed (ramped from 0 to a prescribed value) to all wheels.

/// @} robot_models_OneWheel

}  // namespace OneWheel
}  // namespace chrono