#pragma once
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <ctime>
#include <iomanip>
#include <algorithm>
#include "chrono_thirdparty/filesystem/path.h"
#include "chrono_granular/api/ChApiGranularChrono.h"
#include "chrono_granular/physics/ChGranular.h"
#include "chrono/utils/ChUtilsInputOutput.h"
#include "chrono_granular/utils/ChCudaMathUtils.cuh"
#include "demos/granular/ChGranularDemoUtils.hpp"
#include "chrono/core/ChStream.h"
#include "chrono/core/ChVector.h"
#include "chrono_thirdparty/filesystem/path.h"

using namespace chrono;
using namespace chrono::granular;

struct force_over_x {
    double pos;
    double force;
};

double getMass(double rad, double density) {
    double volume = 4.0f / 3.0f * CH_C_PI * std::pow(rad, 3);
    double mass = volume * density;
    return mass;
}

// calculate kinetic energy of the system
double getSystemKE(double rad, double density, ChGranularSMC_API& apiSMC, int numSpheres) {
    double sysKE = 0.0f;
    double sphere_KE;
    ChVector<float> angularVelo_f;
    ChVector<float> velo_f;

    double mass = getMass(rad, density);
    double inertia = 0.4f * mass * std::pow(rad, 2);

    for (int i = 0; i < numSpheres; i++) {
        angularVelo_f = apiSMC.getAngularVelo(i);
        velo_f = apiSMC.getVelo(i);
        ChVector<double> angularVelo = ChVector<double>(angularVelo_f.x(), angularVelo_f.y(), angularVelo_f.z());
        ChVector<double> velo = ChVector<double>(velo_f.x(), velo_f.y(), velo_f.z());
        sphere_KE = 0.5f * mass * velo.Length2() + 0.5f * inertia * angularVelo.Length2();
        sysKE = sysKE + sphere_KE;
    }
    return sysKE;
}

void calculateBoundaryForces(ChSystemGranularSMC& gran_sys,
                             int numSpheres,
                             double rad,
                             double kn,
                             double gn,
                             double mass,
                             double bottom_plate_position,
                             std::vector<ChVector<double>>& normalForces,
                             std::vector<int>& particlesInContact) {
    double3 velo;
    double penetration;
    double force_multiplier;
    double3 contact_normal = make_double3(0.0, 0.0, 1.0);
    for (int i = 0; i < numSpheres; i++) {
        double3 pos_f = gran_sys.getPositionDouble(i);
        double3 pos = make_double3(pos_f.x, pos_f.y, pos_f.z);

        // check if it's in contact with the bottom boundary
        if (pos.z - rad < bottom_plate_position) {
            penetration = std::abs(pos.z - rad - bottom_plate_position);
            force_multiplier = sqrt(penetration / rad);
            double3 Fn = kn * penetration * contact_normal;

            float3 velo_f = gran_sys.getVelocity(i);
            double3 velo = make_double3(velo_f.x, velo_f.y, velo_f.z);
            double3 rel_vel = velo;

            double projection = Dot(rel_vel, contact_normal);

            // add damping
            Fn = Fn + -1. * gn * projection * contact_normal * mass;
            Fn = Fn * force_multiplier;
            double3 F_el = kn * penetration * contact_normal * force_multiplier;
            double3 F_d = -1. * gn * projection * contact_normal * mass * force_multiplier;

            normalForces.push_back(ChVector<double>(Fn.x, Fn.y, Fn.z));
            particlesInContact.push_back(i);
        }
    }
}

// return boundary force from floor sorted
std::vector<force_over_x> getSortedBoundaryForces(ChSystemGranularSMC& gran_sys,
                                                  int numSpheres,
                                                  double rad,
                                                  double kn,
                                                  double gn,
                                                  double mass,
                                                  double bottom_plate_position) {
    std::vector<ChVector<double>> normalForces;
    std::vector<int> particlesInContact;
    calculateBoundaryForces(gran_sys, numSpheres, rad, kn, gn, mass, bottom_plate_position, normalForces,
                            particlesInContact);

    std::vector<force_over_x> pos_force_array;
    force_over_x myData;
    for (int i = 0; i < normalForces.size(); i++) {
        myData.pos = gran_sys.getPosition(particlesInContact.at(i)).x;
        myData.force = normalForces.at(i).z();
        pos_force_array.push_back(myData);
    }

    // sort data
    std::sort(pos_force_array.begin(), pos_force_array.end(),
              [](auto const& a, auto const& b) { return a.pos < b.pos; });

    return pos_force_array;
}

// initialize velocity, dimenstion of the slab: x_dim_num * radius by y_dim_num * radius
std::vector<ChVector<float>> initializePositions(int x_dim_num, int z_dim_num, float radius) {
    std::vector<ChVector<float>> pos;
    double z = (-(double)z_dim_num + 1.0f) * radius;
    double y = 0;
    double z_diff = std::sqrt(3.1f) * radius;  // add some difference and see what happens
    double x;
    int layers = 0;
    int total_layers = (int)(z_dim_num);
    while (z <= ((double)z_dim_num - 3.0f) * radius - z_diff) {
        x = (-(double)x_dim_num + 1.0f) * radius;
        while (x <= ((double)x_dim_num - 1.0f) * radius) {
            ChVector<float> position(x, y, z);
            pos.push_back(position);
            x = x + 2 * radius;
        }
        layers = layers + 1;

        if (layers == total_layers) {
            break;
        }

        x = (-(double)x_dim_num + 2.0f) * radius;
        z = z + z_diff;
        while (x <= ((double)x_dim_num - 1.0f) * radius) {
            ChVector<float> position(x, y, z);
            pos.push_back(position);
            x = x + 2 * radius;
        }

        layers = layers + 1;
        if (layers == total_layers) {
            break;
        }

        z = z + z_diff;
    }
    return pos;
}

// initialize velocity, dimenstion of the slab: x_dim_num * radius by y_dim_num * radius
std::vector<ChVector<float>> initializeSpherePosition(int x_dim_num, int z_dim_num, float radius) {
    std::vector<ChVector<float>> pos;
    double z = (-(double)z_dim_num - 1.0f) * radius;
    for (int j = 0; j < z_dim_num; j++) {
        float x = j % 2 == 0 ? (-x_dim_num + 1) * radius : (-x_dim_num + 2) * radius;
        for (int i = 0; i < x_dim_num - (j) % 2; i++) {
            pos.push_back(ChVector<>(x, 0, z));
            x += 2 * radius;
        }
        z += std::sqrt(3) * radius;
    }
    //    std::cout << " this row " << pos.back().x() << " " << pos.back().y() << " " << pos.back().z() << std::endl;

    std::cout << "size pos=" << pos.size();
    return pos;
}

// find top center sphere ID
int findTopCenterSphereID(ChSystemGranularSMC& gran_sys, int numSpheres) {
    double max_z = gran_sys.get_max_z();
    float3 particlePos;
    for (int i = 0; i < numSpheres; i++) {
        particlePos = gran_sys.getPosition(i);
        if (std::abs(particlePos.x) < 0.001 && std::abs(particlePos.z - max_z) < 0.01) {
            std::cout << "top sphere ID: " << i << ", pos: " << particlePos.x << ", " << particlePos.z << std::endl;
            return i;
            break;
        }
    }
    printf("ERROR! can not find the particle, check if the setup is correct.\n");
    return -1;
}

// find indexes of particles that are mirrored,
// output index of particles on left and right
void findSymmetricIndices(ChSystemGranularSMC& gran_sys,
                          int numSpheres,
                          std::vector<int>& left_index,
                          std::vector<int>& right_index) {
    for (int i = 0; i < numSpheres; i++) {
        float3 pos_L = gran_sys.getPosition(i);
        if (pos_L.x < -0.05) {
            left_index.push_back(i);  // add left side particle to left array

            // search and find corresponding right ones
            for (int j = 0; j < numSpheres; j++) {
                float3 pos_R = gran_sys.getPosition(j);
                if (std::abs(pos_L.x + pos_R.x) < 1e-2 && std::abs(pos_L.z - pos_R.z) < 1e-2) {
                    right_index.push_back(j);
                    break;
                }
            }
        }
    }
    if (left_index.size() != right_index.size()) {
        std::cout << "left index not equal to the right ones! " << std::endl;
    }
}

// find particles at x = 0, center column of the setup
void findCenterIndices(ChSystemGranularSMC& gran_sys, int numSpheres, std::vector<int>& center_index) {
    for (int i = 0; i < numSpheres; i++) {
        float3 pos = gran_sys.getPosition(i);
        if (std::abs(pos.x) < 0.001) {
            center_index.push_back(i);  // add left side particle to left array
        }
    }
}

// given particle ID on the left, find mirrored particle index on the right
int findMirroredParticleID(std::vector<int> left_index, std::vector<int> right_index, int ID_left) {
    std::vector<int>::iterator it;
    it = std::find(left_index.begin(), left_index.end(), ID_left);
    if (it == left_index.end()) {
        std::cout << "ERROR: given index not found on the left";
        return -1;
    } else {
        int pos = it - left_index.begin();
        return right_index.at(pos);
    }
}

// print debug info
// # of contacts, pos_x, pos_z, penetration, velo_x, velo_z, F_el and F_damp
void printParticleInfo(ChSystemGranularSMC& gran_sys,
                       int myParticleIndex,
                       int numSpheres,
                       float radius,
                       float kn,
                       float gn,
                       float mass,
                       float bottom_plate_position,
                       ChStreamOutAsciiFile& stream) {
    int numContacts = 0;
    float3 mySpherePos = gran_sys.getPosition(myParticleIndex);
    for (int i = 0; i < numSpheres; i++) {
        float3 theirSpherePos = gran_sys.getPosition(i);
        float dist = Length(mySpherePos - theirSpherePos);

        if (dist > radius && dist < 2 * radius) {
            numContacts++;
        }
    }

    // get penetration etc
    float3 contact_normal = make_float3(0.0f, 0.0f, 1.0f);
    float penetration = std::abs(mySpherePos.z - radius - bottom_plate_position);

    // elastic component of contact force
    float force_multiplier = sqrt(penetration / radius);
    float3 Fn = kn * penetration * contact_normal * force_multiplier;

    float3 velo = gran_sys.getVelocity(myParticleIndex);
    // if (myParticleIndex == 200){
    //     printf("print velocity: %e, %e\n", velo.x, velo.z);
    // }
    float3 rel_vel = velo;
    float projection = Dot(rel_vel, contact_normal);

    // damping component of contact force
    float3 F_d = -1. * gn * projection * contact_normal * mass * force_multiplier;

    stream << numContacts << ", " << mySpherePos.x << ", " << mySpherePos.z << ", " << penetration << ", " << velo.x
           << ", " << velo.z << ", " << Fn.z << ", " << F_d.z << "\n";
}

// return number of contacts per particle
int findNumContacts(ChSystemGranularSMC& gran_sys, int numSpheres, int myParticleIndex, float radius) {
    int numContacts = 0;
    float3 mySpherePos = gran_sys.getPosition(myParticleIndex);
    for (int i = 0; i < numSpheres; i++) {
        float3 theirSpherePos = gran_sys.getPosition(i);
        float dist = Length(mySpherePos - theirSpherePos);

        if (dist > radius && dist < 2 * radius) {
            numContacts++;
        }
    }
    return numContacts;
}
