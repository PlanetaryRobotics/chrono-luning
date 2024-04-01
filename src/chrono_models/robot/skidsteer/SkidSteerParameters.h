#pragma once

#include "chrono/serialization/ChArchive.h"
#include "chrono/serialization/ChArchiveJSON.h"

#include <stdio.h>
#include <iostream>

using namespace chrono;
class SkidSteerParameters {
  public:
    // Wheel geometry
    float wheel_rad;
    float wheel_width;
    float wheel_mass;

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

    virtual void ArchiveOut(ChArchiveOut& marchive) {
      marchive << CHNVP(wheel_rad);
      marchive << CHNVP(wheel_width);
      marchive << CHNVP(wheel_mass);
      marchive << CHNVP(wheel_x);
      marchive << CHNVP(wheel_y);
      marchive << CHNVP(wheel_z);
      marchive << CHNVP(wheel_mu);
      marchive << CHNVP(cr);
      marchive << CHNVP(Y);
      marchive << CHNVP(nu);
      marchive << CHNVP(kn);
      marchive << CHNVP(gn);
      marchive << CHNVP(kt);
      marchive << CHNVP(gt);
      marchive << CHNVP(chassis_mass);
      marchive << CHNVP(chassis_dim_x);
      marchive << CHNVP(chassis_dim_y);
      marchive << CHNVP(chassis_dim_z);
      marchive << CHNVP(wheel_mesh_name);
      marchive << CHNVP(chassis_mesh_name);
    }
    virtual void ArchiveIn(ChArchiveIn& marchive) {
      marchive >> CHNVP(wheel_rad);
      marchive >> CHNVP(wheel_width);
      marchive >> CHNVP(wheel_mass);
      marchive >> CHNVP(wheel_x);
      marchive >> CHNVP(wheel_y);
      marchive >> CHNVP(wheel_z);
      marchive >> CHNVP(wheel_mu);
      marchive >> CHNVP(cr);
      marchive >> CHNVP(Y);
      marchive >> CHNVP(nu);
      marchive >> CHNVP(kn);
      marchive >> CHNVP(gn);
      marchive >> CHNVP(kt);
      marchive >> CHNVP(gt);
      marchive >> CHNVP(chassis_mass);
      marchive >> CHNVP(chassis_dim_x);
      marchive >> CHNVP(chassis_dim_y);
      marchive >> CHNVP(chassis_dim_z);
      marchive >> CHNVP(wheel_mesh_name);
      marchive >> CHNVP(chassis_mesh_name);

  }
};