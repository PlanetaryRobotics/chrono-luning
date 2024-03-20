#include "chrono/serialization/ChArchive.h"
#include "chrono/serialization/ChArchiveJSON.h"

#include <stdio.h>
#include <iostream>

using namespace chrono;

class EnvironmentParameters {
  public:
    // Soil parameters
    float mu;
    float mu_wall;  // Friction coefficient with the wall
    float CoR;
    float E;

    float w_r;       // Rover wheel driving speed (rad/s)
    float z_offset;  // Height to spawn the rover

    float time_end;
    unsigned int fps;
    unsigned int report_frequency;

    virtual ~EnvironmentParameters() {}

    virtual void ArchiveOut(ChArchiveOut& marchive) {
        marchive << CHNVP(mu);
        marchive << CHNVP(mu_wall);
        marchive << CHNVP(CoR);
        marchive << CHNVP(E);
        marchive << CHNVP(w_r);
        marchive << CHNVP(z_offset);
        marchive << CHNVP(time_end);
        marchive << CHNVP(fps);
        marchive << CHNVP(report_frequency);
    }

    virtual void ArchiveIn(ChArchiveIn& marchive) {
        marchive >> CHNVP(mu);
        marchive >> CHNVP(mu_wall);
        marchive >> CHNVP(CoR);
        marchive >> CHNVP(E);
        marchive >> CHNVP(w_r);
        marchive >> CHNVP(z_offset);
        marchive >> CHNVP(time_end);
        marchive >> CHNVP(fps);
        marchive >> CHNVP(report_frequency);
    }
}