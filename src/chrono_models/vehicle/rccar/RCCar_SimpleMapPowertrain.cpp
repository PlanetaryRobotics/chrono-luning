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
// Authors: Radu Serban, Jayne Henry, Luning Fang
// =============================================================================
//
// Simple powertrain model for the RCCar vehicle.
// - based on torque-speed engine maps
// - both power and torque limited
// - no torque converter
// - simple gear-shifting model (in automatic mode)
//
// =============================================================================

#include "chrono_models/vehicle/rccar/RCCar_SimpleMapPowertrain.h"

using namespace chrono::vehicle;
using namespace chrono;

namespace chrono {
namespace vehicle {
namespace rccar {

const double rpm2rads = CH_C_PI / 30;

RCCar_SimpleMapPowertrain::RCCar_SimpleMapPowertrain(const std::string& name) : ChSimpleMapPowertrain(name) {}

double RCCar_SimpleMapPowertrain::GetMaxEngineSpeed() {
    return Kv_rating * supply_voltage * max_voltage_ratio * rpm2rads;
}

/// specify ratio of maximum engine voltage
void RCCar_SimpleMapPowertrain::SetMaxMotorVoltageRatio(double maxVoltageRatio){
    max_voltage_ratio = maxVoltageRatio;        
}

/// specify stall torque
void RCCar_SimpleMapPowertrain::SetStallTorque(double stall_torque){
    stallTorque = stall_torque;
}


void RCCar_SimpleMapPowertrain::SetEngineTorqueMaps(ChFunction_Recorder& map0, ChFunction_Recorder& mapF) {
    double max_rpm = Kv_rating * supply_voltage * max_voltage_ratio;
    // since this is a model of motor and ESC combination, we assume a linear relationship.
    // while brushless motors dont follow linear torque-speed relationship, most hobby electronic
    // speed controllers control the motor such that it approximately follows a linear relationship
    mapF.AddPoint(0,stallTorque); //stall torque //TODO
    mapF.AddPoint(max_rpm * rpm2rads, 0);  // no load speed

    // N-m and rad/s
    map0.AddPoint(0, 0);
    // map0.AddPoint(.1 * max_rpm * rpm2rads, 0);
    // map0.AddPoint(max_rpm * rpm2rads, -stallTorque);  // TODO, currently a guess
    // map0.AddPoint(.1 * max_rpm * rpm2rads, 0);
    map0.AddPoint(max_rpm * rpm2rads, 0);  // TODO, currently a guess

}

void RCCar_SimpleMapPowertrain::SetGearRatios(std::vector<double>& fwd, double& rev) {
    rev = -1.0 / 3;
    fwd.push_back(1.0 / 3);
}

void RCCar_SimpleMapPowertrain::SetShiftPoints(std::vector<std::pair<double, double>>& shift_bands) {
    double max_rpm = Kv_rating * supply_voltage;

    shift_bands.push_back(std::pair<double, double>(0, 10 * max_rpm * rpm2rads)); //never shifts
}

}  // end namespace rccar
}  // namespace vehicle
}  // namespace chrono
