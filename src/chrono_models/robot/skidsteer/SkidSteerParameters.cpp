#include "SkidSteerParameters.h"
virtual void SkidSteerParameters::ArchiveOut(ChArchiveOut& marchive) {
    marchive << CHNVP(mu_wheel);
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
    marchive << CHNVP(model_file_path);
}

virtual void SkidSteerParameters::ArchiveIn(ChArchiveIn& marchive) {
    marchive >> CHNVP(mu_wheel);
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
    marchive >> CHNVP(model_file_path);
}
