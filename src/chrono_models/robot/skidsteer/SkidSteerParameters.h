#pragma once

#include "chrono/serialization/ChArchive.h"
#include "chrono/serialization/ChArchiveJSON.h"

#include <stdio.h>
#include <iostream>

using namespace chrono;
class SkidSteerParameters {
  public:
    float mu_wheel;

    // Wheel geometry
    float wheel_rad;
    float wheel_width;
    float m_wheel;

    float wheel_x;
    float wheel_y;
    float wheel_z;

    // Wheel material parameters
    float wheel_mu;
    float cr;
    float Y;
    float nu;
    float kn;
    float gn;
    float kt;
    float gt;

    // Chassis intertia
    float chassis_mass;

    float chassis_dim_x;
    float chassis_dim_y;
    float chassis_dim_z;

    // Model
    std::string wheel_mesh_name;
    std::string chassis_mesh_name;

    SkidSteerParameters(){};

    virtual ~SkidSteerParameters() {}

    virtual void ArchiveOut(ChArchiveOut& marchive);
    virtual void ArchiveIn(ChArchiveIn& marchive);
};