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
// Authors: Jason Zhou
// =============================================================================
//
// NASA VIPER Lunar Rover Model Class.
// This class contains model for NASA's VIPER lunar rover for NASA's 2024 Moon
// exploration mission.
//
// =============================================================================

#ifndef VIPER_H
#define VIPER_H

#include <string>
#include <array>

#include "chrono/assets/ChColor.h"
#include "chrono/physics/ChLinkMotorRotationSpeed.h"
#include "chrono/physics/ChSystem.h"
#include "chrono/physics/ChShaft.h"
#include "chrono/physics/ChShaftsGear.h"
#include "chrono/physics/ChShaftsBody.h"

#include "chrono_models/ChApiModels.h"

namespace chrono {

/// Namespace with classes for the Viper model.
namespace viper {

/// @addtogroup robot_models_viper
/// @{

/// Viper wheel/suspension identifiers.
enum WheelID {
    LF = 0,  ///< left front
    RF = 1,  ///< right front
    LB = 2,  ///< left back
    RB = 3   ///< right back
};

/// Viper turning signal.
enum class TurnSig {
    LEFT = +1,   ///< left turn signal
    RIGHT = -1,  ///< right turn signal
    HOLD = 0     ///< hold signal
};

/// Viper wheel type.
enum class Wheel_Type {
    RealWheel,    ///< actual geometry of the Viper wheel
    SimpleWheel,  ///< simplified wheel geometry
    CylWheel      ///< cylindrical wheel geometry
};

/// Base class definition of the Viper Rover Part.
/// Viper Rover Parts include Chassis, Steering, Upper Suspension Arm, Bottom Suspension Arm and Wheel.
/// This class encapsulates base fields and functions.
class CH_MODELS_API Viper_Part {
  public:
    Viper_Part(const std::string& name, const ChFrame<>& rel_pos, std::shared_ptr<ChMaterialSurface> mat, bool collide);
    virtual ~Viper_Part() {}

    /// Return the name of the part.
    const std::string& GetName() const { return m_name; }

    /// Set the name of the part.
    void SetName(const std::string& name) { m_name = name; }

    /// Enable/disable visualization.
    void SetVisualize(bool state);

    /// Enable/disable collision.
    void SetCollide(bool state);

    /// Initialize the rover part by attaching it to the specified chassis body.
    void Initialize(std::shared_ptr<ChBodyAuxRef> chassis);

    /// Return the ChBody of the corresponding Viper part.
    std::shared_ptr<ChBodyAuxRef> GetBody() const { return m_body; }

    /// Return the position of the Viper part.
    /// This is the absolute location of the part reference frame.
    const ChVector<>& GetPos() const { return m_body->GetFrame_REF_to_abs().GetPos(); }

    /// Return the rotation of the Viper part.
    /// This is the orientation wrt the global frame of the part reference frame.
    const ChQuaternion<>& GetRot() const { return m_body->GetFrame_REF_to_abs().GetRot(); }

    /// Return the linear velocity of the Viper part.
    /// This is the absolute linear velocity of the part reference frame.
    const ChVector<>& GetLinVel() const { return m_body->GetFrame_REF_to_abs().GetPos_dt(); }

    /// Return the angular velocity of the Viper part.
    /// This is the absolute angular velocity of the part reference frame.
    const ChVector<> GetAngVel() const { return m_body->GetFrame_REF_to_abs().GetWvel_par(); }

    /// Return the linear acceleration of the Viper part.
    /// This is the absolute linear acceleration of the part reference frame.
    const ChVector<>& GetLinAcc() const { return m_body->GetFrame_REF_to_abs().GetPos_dtdt(); }

    /// Return the angular acceleratino of the Viper part.
    /// This is the absolute angular acceleratin of the part reference frame.
    const ChVector<> GetAngAcc() const { return m_body->GetFrame_REF_to_abs().GetWacc_par(); }

  protected:
    /// Complete construction of the part.
    void Construct(ChSystem* system);

    std::string m_name;                        ///< part name
    std::shared_ptr<ChBodyAuxRef> m_body;      ///< part rigid body
    std::shared_ptr<ChMaterialSurface> m_mat;  ///< contact material (shared among all shapes)

    std::string m_mesh_name;  ///< visualization mesh name
    ChColor m_color;          ///< visualization asset color

    ChFrame<> m_pos;   ///< relative position wrt chassis
    double m_density;  ///< part density

    bool m_visualize;  ///< part visualization flag
    bool m_collide;    ///< part collision flag
};

/// Viper rover Chassis.
class CH_MODELS_API Viper_Chassis : public Viper_Part {
  public:
    Viper_Chassis(const std::string& name, std::shared_ptr<ChMaterialSurface> mat);
    ~Viper_Chassis() {}

    /// Initialize the chassis at the specified (absolute) position.
    void Initialize(ChSystem* system, const ChFrame<>& pos);
};

/// Viper rover Wheel.
class CH_MODELS_API Viper_Wheel : public Viper_Part {
  public:
    Viper_Wheel(const std::string& name,
                const ChFrame<>& rel_pos,
                std::shared_ptr<ChMaterialSurface> mat,
                Wheel_Type wheel_type);
    ~Viper_Wheel() {}

    friend class ViperRover;
};

/// The upper arm of the Viper rover suspension.
class CH_MODELS_API Viper_Up_Arm : public Viper_Part {
  public:
    Viper_Up_Arm(const std::string& name,
                 const ChFrame<>& rel_pos,
                 std::shared_ptr<ChMaterialSurface> mat,
                 const int& side);  ///< indicate which side of the suspension 0->L, 1->R
    ~Viper_Up_Arm() {}
};

/// The bottom arm of the Viper rover suspension.
class CH_MODELS_API Viper_Bottom_Arm : public Viper_Part {
  public:
    Viper_Bottom_Arm(const std::string& name,
                     const ChFrame<>& rel_pos,
                     std::shared_ptr<ChMaterialSurface> mat,
                     const int& side);  ///< indicate which side of the suspension 0->L, 1->R
    ~Viper_Bottom_Arm() {}
};

/// Steering rod of the Viper rover.
/// The steering rod is connected to the steering cyl, this link is controlled steering.
/// There are two connecting rods on the steering rod, linking to upper and bottom arms of the suspension.
class CH_MODELS_API Viper_Steer : public Viper_Part {
  public:
    Viper_Steer(const std::string& name,
                const ChFrame<>& rel_pos,
                std::shared_ptr<ChMaterialSurface> mat,
                const int& side);  ///< indicate which side of the rover 0->L, 1->R
    ~Viper_Steer() {}
};

/// Viper rover class.
/// This class encapsulates the location and rotation information of all Viper parts wrt the chassis.
/// This class should be the entry point to create a complete rover.
class CH_MODELS_API ViperRover {
  public:
    ViperRover(ChSystem* system, Wheel_Type wheel_type = Wheel_Type::RealWheel);

    ~ViperRover();

    /// Get the containing system.
    ChSystem* GetSystem() const { return m_system; }

    /// Enable/disable DC motor control.
    void SetDCControl(bool dc_control);

    /// Set motor stall torque for the specified wheel.
    /// This value is used only if DC motor control is enabled.
    void SetMotorStallTorque(double torque, WheelID id);

    /// Set DC motor no load speed.
    /// This value is used only if DC motor control is enabled.
    void SetMotorNoLoadSpeed(double rad_speed, WheelID id);

    /// Set motor speed.
    /// This value is used only if DC motor control is disabled.
    void SetMotorSpeed(double rad_speed, WheelID id);

    /// Set lift motor speed
    void SetLiftMotorSpeed(double rad_speed, WheelID id);

    /// Set wheel contact material.
    void SetWheelContactMaterial(std::shared_ptr<ChMaterialSurface> mat);

    /// Fix the chassis to ground.
    void SetChassisFixed(bool fixed);

    /// Enable/disable visualization of the rover chassis (default: true).
    void SetChassisVisualization(bool state);

    /// Enable/disable visualization of rover wheels (default: true).
    void SetWheelVisualization(bool state);

    /// Enable/disable visualization of rover suspensions (default: true).
    void SetSuspensionVisualization(bool state);

    /// Initialize the Viper rover at the specified position.
    void Initialize(const ChFrame<>& pos);

    /// Get the chassis part.
    std::shared_ptr<Viper_Chassis> GetChassis() const { return m_chassis; }

    /// Get the wheel part.
    std::shared_ptr<Viper_Wheel> GetWheel(WheelID id) const { return m_wheels[id]; }

    /// Get the steering part.
    std::shared_ptr<Viper_Steer> GetSteering(WheelID id) const { return m_steers[id]; }

    /// Get the upper arm part.
    std::shared_ptr<Viper_Up_Arm> GetUpArm(WheelID id) const { return m_up_suss[id]; }

    /// Get the bottom arm body.
    std::shared_ptr<Viper_Bottom_Arm> GetBottomArm(WheelID id) const { return m_bts_suss[id]; }

    /// Get chassis position.
    ChVector<> GetChassisPos() const { return m_chassis->GetPos(); }

    /// Get chassis orientation.
    ChQuaternion<> GetChassisRot() const { return m_chassis->GetRot(); }

    /// Get chassis linear velocity.
    ChVector<> GetChassisVel() const { return m_chassis->GetLinVel(); }

    /// Get chassis linear acceleration.
    ChVector<> GetChassisAcc() const { return m_chassis->GetLinAcc(); }

    /// Get wheel speed.
    ChVector<> GetWheelLinVel(WheelID id) const { return m_wheels[id]->GetLinVel(); }

    /// Get wheel angular velocity.
    ChVector<> GetWheelAngVel(WheelID id) const { return m_wheels[id]->GetAngVel(); }

    /// Get wheel contact force.
    ChVector<> GetWheelContactForce(WheelID id) const;

    /// Get wheel contact torque.
    ChVector<> GetWheelContactTorque(WheelID id) const;

    /// Get wheel total applied force.
    ChVector<> GetWheelAppliedForce(WheelID id) const;

    /// Get wheel tractive torque - if DC control set to off
    double GetWheelTracTorque(WheelID id) const;

    /// Get wheel total applied torque.
    ChVector<> GetWheelAppliedTorque(WheelID id) const;

    /// Get total rover mass.
    double GetRoverMass() const;

    /// Get total wheel mass.
    double GetWheelMass() const;

    /// Get main motor function.
    std::shared_ptr<ChFunction_Const> GetMainMotorFunc(WheelID id);

    /// Get steer motor function.
    std::shared_ptr<ChFunction_Const> GetSteerMotorFunc(WheelID id);

    /// Get main motor link.
    std::shared_ptr<ChLinkMotorRotationSpeed> GetMainMotorLink(WheelID id);

    /// Get steer motor link.
    std::shared_ptr<ChLinkMotorRotationSpeed> GetSteerMotorLink(WheelID id);

    /// Set viper turning signal left/right/hold.
    void SetTurn(TurnSig id, double turn_speed = 0.0);

    /// Get viper turning angle - ranges from -CH_C_PI/3 to CH_C_PI/3.
    double GetTurnAngle() const;

    /// Get viper turning state - HOLD, LEFT, OR RIGHT.
    TurnSig GetTurnState() const;

    /// A viper status check and update function.
    /// Note: this function needs to be included in the main simulation loop.
    void Update();

    /// Sub-update function to update DC motor limit control
    void UpdateDCMotorControl();

    /// Sub-update function to update steering control
    void UpdateSteeringControl();

  private:
    /// Create the rover parts.
    void Create(Wheel_Type wheel_type);

    ChSystem* m_system;  ///< pointer to the Chrono system

    bool m_chassis_fixed;     ///< fix chassis to ground
    bool m_dc_motor_control;  ///< use DC motor controller

    std::shared_ptr<Viper_Chassis> m_chassis;                     ///< rover chassis
    std::array<std::shared_ptr<Viper_Wheel>, 4> m_wheels;         ///< rover wheels (LF, RF, LR, RB)
    std::array<std::shared_ptr<Viper_Up_Arm>, 4> m_up_suss;       ///< rover upper suspensions (LF, RF, LR, RB)
    std::array<std::shared_ptr<Viper_Bottom_Arm>, 4> m_bts_suss;  ///< rover bottom suspensions (LF, RF, LR, RB)
    std::array<std::shared_ptr<Viper_Steer>, 4> m_steers;         ///< rover steering stand (LF, RF, LR, RB)
    std::array<std::shared_ptr<ChBody>, 4> m_steers_rod;          ///< rover steering rod (LF, RF, LR, RB)

    std::array<std::shared_ptr<ChLinkMotorRotationSpeed>, 4> m_motors;        ///< drive motors
    std::array<std::shared_ptr<ChLinkMotorRotationSpeed>, 4> m_steer_motors;  ///< steering motors
    std::array<std::shared_ptr<ChLinkMotorRotationSpeed>, 8> m_lift_motors;   ///< lifting motors

    // DC Motor Test
    std::array<std::shared_ptr<ChShaft>, 4> m_power_shafts;
    std::array<std::shared_ptr<ChShaft>, 4> m_driven_shafts;
    std::array<std::shared_ptr<ChShaftsGear>, 4> m_shaft_gears;

    std::array<double, 4> m_stall_torque;   ///< stall torque of the motors
    std::array<double, 4> m_no_load_speed;  ///< no load speed of the motors

    TurnSig m_turn_state;  ///< turning state of the rover

    std::array<std::shared_ptr<ChFunction_Const>, 4> m_motors_func;        ///< motor angular speed func
    std::array<std::shared_ptr<ChFunction_Const>, 4> m_steer_motors_func;  ///< steering motor angular speed func
    std::array<std::shared_ptr<ChFunction_Const>, 4> m_lift_motors_func;   ///< lifting motor angular speed func

    std::array<std::shared_ptr<ChLinkTSDA>, 4> m_sus_springs;  ///< suspension springs

    std::shared_ptr<ChMaterialSurface> m_default_material;  ///< common contact material for all non-wheel parts
    std::shared_ptr<ChMaterialSurface> m_wheel_material;    ///< wheel contact material (shared across limbs)
};

/// @} robot_models_viper

}  // namespace viper
}  // namespace chrono
#endif
