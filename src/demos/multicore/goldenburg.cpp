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
// Authors: Radu Serban
// =============================================================================
//
// Chrono::Multicore test program using SMC method for frictional contact.
//
// The model simulated here consists of a number of spherical objects falling
// in a fixed container.
//
// The global reference frame has Z up.
//
// If available, OpenGL is used for run-time rendering. Otherwise, the
// simulation is carried out for a pre-defined duration and output files are
// generated for post-processing with POV-Ray.
// =============================================================================

// unit cgs

#include <cstdio>
#include <vector>
#include <cmath>
#include <stdio.h>
#include <ctime>
#include <sys/time.h>


#include "chrono_multicore/physics/ChSystemMulticore.h"

#include "chrono/core/ChTimer.h"
#include "chrono/ChConfig.h"
#include "chrono/utils/ChUtilsCreators.h"
#include "chrono/utils/ChUtilsInputOutput.h"




#ifdef CHRONO_OPENGL
#include "chrono_opengl/ChVisualSystemOpenGL.h"
#endif

using namespace chrono;
using namespace chrono::collision;

// Tilt angle (about global Y axis) of the container.
double tilt_angle = 0;

// Number of balls: (2 * count_X + 1) * (2 * count_Y + 1)
float sphere_radius = 1;
int count_X = 60;
int count_Y = 60;
int count_Z = 15;


// int count_X = 2;
// int count_Y = 2;
// int count_Z = 3;


// int count_X = 1;
// int count_Y = 1;
// int count_Z = 1;


// Material properties (same on bin and balls)
float Y = 2e8f;
float cr = 0.01f;
float damp = 1e4f;

auto top_middle = chrono_types::make_shared <ChBody>();

double rho = 7.4;
double gravity = 1000;
//double kn_ratio = 3e7;  // step size: 1e-7 for kn ratio 3e7

// double kn_ratio = 3e5;  // step size: 1e-6 for kn ratio 3e5
// double time_step = 1e-6; // 


double kn_ratio = 3e5;  // step size: 1e-6 for kn ratio 3e5
double time_step = 1e-6; // 

double F_ext_ratio = 100.f;

// -----------------------
// Output kinetic energy
// -----------------------
double calcKE(ChSystemMulticoreSMC* sys){
    double KE = 0;
    for (auto body: sys->Get_bodylist()){


        ChVector<> eng_trn = 0.5 * body->GetMass() * body->GetPos_dt() * body->GetPos_dt();
        ChVector<> eng_rot = 0.5 * body->GetInertiaXX() * body->GetWvel_par() * body->GetWvel_par();

        double KE_trn = eng_trn.x() + eng_trn.y() + eng_trn.z();
        double KE_rot = eng_rot.x() + eng_rot.y() + eng_rot.z();
        double KE_tot = KE_trn + KE_rot;

        KE = KE + KE_tot;

    }
    return KE;
}

// -----------------------------------------------------------------------------
// Callback class for contact reporting
// -----------------------------------------------------------------------------
class ContactReporter : public ChContactContainer::ReportContactCallback {
  public:
    ContactReporter(std::shared_ptr<ChBody> obj1){
        m_obj1 = obj1;
        m_csv << "pos" << "force" << std::endl;
    }

    void getNormalForce(std::vector<double>& normal_F_array){normal_F_array = contact_force_z;}

    void writeNormalForces(const std::string& filename){        
        m_csv.write_to_file(filename);
        std::cout << "Successfully write csv file! " << std::endl;
    }

  private:
    virtual bool OnReportContact(const ChVector<>& pA,
                                 const ChVector<>& pB,
                                 const ChMatrix33<>& plane_coord,
                                 const double& distance,
                                 const double& eff_radius,
                                 const ChVector<>& cforce,
                                 const ChVector<>& ctorque,
                                 ChContactable* modA,
                                 ChContactable* modB) override {


        if (modA == m_obj1.get()){
            auto bodyB = static_cast<ChBody*>(modB);
            const ChVector<>& nrm = plane_coord.Get_A_Xaxis();
            ChVector<> force_global = plane_coord * cforce;
                if (std::abs(std::abs(nrm.z()) - 1) < 1e-5){
                    // report contact force from the plate in normal direction
                    // contact_force_z.push_back(cforce.z());
                    // contact_force_z.push_back(cforce.z());

                    if (bodyB->GetId() > count_Y + 1){
                        fprintf(stderr, "contact pairs out of range! new sphere in contact %d \n", bodyB->GetId());
                        return -1;
                    }


                    contact_force_z.at( bodyB->GetId() - 1) = force_global.z();
                    m_csv << std::scientific << pA.x() << force_global.z() << std::endl;
                }

        }

        if (modA == top_middle.get()){
            ChVector<> force_global = plane_coord * cforce;
            std::cout << "top midlle contact force: " << force_global.x() << ", " << force_global.y() << ", " << force_global.z() << std::endl;
        }
        return true;
    }

    std::shared_ptr<ChBody> m_obj1;
    std::shared_ptr<ChBody> m_obj2;
    std::vector<double> contact_force_z = std::vector<double>(count_Y+1);
    utils::CSV_writer m_csv;
};

float time_diff(struct timeval *start, struct timeval *end)
{
    return (end->tv_sec - start->tv_sec) + 1e-6*(end->tv_usec - start->tv_usec);
}


// -----------------------------------------------------------------------------
// Create a bin consisting of five boxes attached to the ground.
// -----------------------------------------------------------------------------
void AddContainer(ChSystemMulticoreSMC* sys, double sphere_radius, double mu) {
    // IDs for the container
    int binId = -200;


    // Create a common material
    auto mat = chrono_types::make_shared<ChMaterialSurfaceSMC>();

    double volume = 4.f/3.f * CH_C_PI * sphere_radius * sphere_radius * sphere_radius;
    double mass = rho * volume;

    mat->SetYoungModulus(1e7);
    mat->SetKn(kn_ratio * mass * gravity/sphere_radius);
    mat->SetKt(kn_ratio * mass * gravity/sphere_radius);

    mat->SetFriction(mu);
    mat->SetRestitution(cr);

    // mat->SetGn(damp);
    // mat->SetGt(damp);
    // mat->SetKn(3000*9.81*1/0.1);

    // Create the containing bin (4 x 4 x 1)
    auto bin = std::shared_ptr<ChBody>(sys->NewBody());
    bin->SetIdentifier(binId);
    bin->SetMass(10000);
    bin->SetPos(ChVector<>(0, 0, 0));
    bin->SetRot(Q_from_AngY(tilt_angle));
    bin->SetCollide(true);
    bin->SetBodyFixed(true);

    // 61 spheres per layer, total of 15 layers
    double box_hdim = (count_X + 1) * sphere_radius;
    ChVector<> hdim(box_hdim, box_hdim, count_Z * 2 * sphere_radius);
    double hthick = 2*sphere_radius;


    bin->GetCollisionModel()->ClearModel();
    utils::AddBoxGeometry(bin.get(), mat, ChVector<>(hdim.x(), hdim.y(), hthick),  ChVector<>(0, 0, -hthick));
    utils::AddBoxGeometry(bin.get(), mat, ChVector<>(hthick, hdim.y(), hdim.z()),  ChVector<>(-hdim.x() - hthick, 0, hdim.z()));
    utils::AddBoxGeometry(bin.get(), mat, ChVector<>(hthick, hdim.y(), hdim.z()),  ChVector<>( hdim.x() + hthick, 0, hdim.z()));
    utils::AddBoxGeometry(bin.get(), mat, ChVector<>(hdim.x(), hthick, hdim.z()),  ChVector<>(0, -hdim.y() - hthick, hdim.z()));
    utils::AddBoxGeometry(bin.get(), mat, ChVector<>(hdim.x(), hthick, hdim.z()),  ChVector<>(0,  hdim.y() + hthick, hdim.z()));

    bin->GetCollisionModel()->BuildModel();

    sys->AddBody(bin);

    fprintf(stderr, "Add container\n");

    
}

// -----------------------------------------------------------------------------
// Create the falling spherical objects in a uniform rectangular grid.
// -----------------------------------------------------------------------------
void AddFallingBalls(ChSystemMulticoreSMC* sys, double sphere_radius, double mu) {
    // Common material
    // Create the falling balls
    int ballId = 0;
    double volume = 4.f/3.f * CH_C_PI * sphere_radius * sphere_radius * sphere_radius;
    double mass = rho * volume;
    double spacing = std::sqrt(3) * sphere_radius;

    auto ballMat = chrono_types::make_shared<ChMaterialSurfaceSMC>();

    ballMat->SetYoungModulus(1e7);
    ballMat->SetKn(kn_ratio * mass * gravity/sphere_radius);
    ballMat->SetKt(kn_ratio * mass * gravity/sphere_radius);
    ballMat->SetFriction(mu);
    ballMat->SetRestitution(cr);


    ChVector<> inertia = (2.0 / 5.0) * mass * sphere_radius * sphere_radius * ChVector<>(1, 1, 1);
    for (int iz = 0; iz < count_Z; iz++){
        if ( iz%2 == 0){
            for (int ix = -count_X/2; ix <= count_X/2; ix++) {
                ChVector<> pos((2*sphere_radius) * ix, 0,  sphere_radius + spacing * iz);
                auto ball = std::shared_ptr<ChBody>(sys->NewBody());
                ball->SetIdentifier(ballId++);
                ball->SetMass(mass);
                ball->SetInertiaXX(inertia);
                ball->SetPos(pos);
                ball->SetRot(ChQuaternion<>(1, 0, 0, 0));
                ball->SetBodyFixed(false);
                ball->SetCollide(true);

                ball->GetCollisionModel()->ClearModel();
                utils::AddSphereGeometry(ball.get(), ballMat, sphere_radius);
                ball->GetCollisionModel()->BuildModel();

                sys->AddBody(ball);

                if( ix == 0 && iz == count_Z - 1){
                    top_middle = ball; // get pointer to the top middle sphere
                }
                        
                }
        }else{
                for (int ix = -count_X/2; ix <= count_X/2-1; ix++) {

                    ChVector<> pos( (2*sphere_radius) * ix + sphere_radius, 0, sphere_radius + spacing * iz);

                    auto ball = std::shared_ptr<ChBody>(sys->NewBody());
                    ball->SetIdentifier(ballId++);
                    ball->SetMass(mass);
                    ball->SetInertiaXX(inertia);
                    ball->SetPos(pos);
                    ball->SetRot(ChQuaternion<>(1, 0, 0, 0));
                    ball->SetBodyFixed(false);
                    ball->SetCollide(true);

                    ball->GetCollisionModel()->ClearModel();
                    utils::AddSphereGeometry(ball.get(), ballMat, sphere_radius);
                    ball->GetCollisionModel()->BuildModel();

                    sys->AddBody(ball);


                    
                }
        }

    }
}

// -----------------------------------------------------------------------------
// Create the system, specify simulation parameters, and run simulation loop.
// -----------------------------------------------------------------------------
int main(int argc, char* argv[]) {
    GetLog() << "Copyright (c) 2017 projectchrono.org\nChrono version: " << CHRONO_VERSION << "\n\n";

    // Simulation parameters
    // ---------------------
    double mu = 0.5;
    double time_settle = 0.1;
    double time_F_dur = 0.03;
    double time_end = 0.3;

    // variables for timing
    struct timeval start;
    struct timeval end;
    


    uint max_iteration = 100;
    real tolerance = 1e-3;

    double out_fps = 100;

    // Create system
    // -------------

    ChSystemMulticoreSMC msystem;

    // Set number of threads
    msystem.SetNumThreads(8);

    // Set gravitational acceleration
    msystem.Set_G_acc(ChVector<>(0, 0, -gravity));

    // Set solver parameters
    msystem.GetSettings()->solver.max_iteration_bilateral = max_iteration;
    msystem.GetSettings()->solver.tolerance = tolerance;

    msystem.GetSettings()->collision.narrowphase_algorithm = ChNarrowphase::Algorithm::MPR;
    msystem.GetSettings()->collision.bins_per_axis = vec3(4, 4, 4);

    // The following two lines are optional, since they are the default options. They are added for future reference,
    // i.e. when needed to change those models.
    // msystem.GetSettings()->solver.contact_force_model = ChSystemSMC::ContactForceModel::Hertz;
    msystem.GetSettings()->solver.contact_force_model = ChSystemSMC::ContactForceModel::Hooke;
    msystem.GetSettings()->solver.tangential_displ_mode = ChSystemSMC::TangentialDisplacementModel::MultiStep;

    msystem.GetSettings()->solver.use_material_properties = true;

    // Create the fixed and moving bodies
    // ----------------------------------
    AddContainer(&msystem, sphere_radius, mu);
    AddFallingBalls(&msystem, sphere_radius, mu);

    // 


    // Perform the simulation

    // Run simulation for specified time
    int num_steps = (int)std::ceil(time_end / time_step);
    double curr_time = 0;
    double normal_force = 0;
    auto creporter = chrono_types::make_shared<ContactReporter>(msystem.Get_bodylist().at(0));

    // #ifdef CHRONO_OPENGL
    //     opengl::ChVisualSystemOpenGL vis;
    //     vis.AttachSystem(&msystem);
    //     vis.SetWindowTitle("Balls SMC");
    //     vis.SetWindowSize(1280, 720);
    //     vis.SetRenderMode(opengl::WIREFRAME);
    //     vis.Initialize();
    //     vis.SetCameraPosition(ChVector<>(0, -10, 2), ChVector<>(0, 0, 0));
    //     vis.SetCameraVertical(CameraVerticalDir::Z);

    // while (true) {
    //     if (vis.Run()) {
    //         msystem.DoStepDynamics(time_step);
    //         vis.Render();
    //         // Print cumulative contact force on container bin.
    //         real3 frc = msystem.GetBodyContactForce(0);
    //         std::cout << "t = " << msystem.GetChTime() << ", KE: " << calcKE(&msystem) << std::endl;
    //         msystem.GetContactContainer()->ReportAllContacts(creporter);
    //         std::vector<double> normalF;
    //         creporter->getNormalForce(normalF);

    //         std::cout << ", F_l = " << normalF.at(0) << ", F_r = "  << normalF.at(count_Y) << std::endl;

    //     } else {
    //         break;
    //     }
    // }
    // #else

    bool useVis = false;
    opengl::ChVisualSystemOpenGL vis;

    if (useVis == true){
        vis.AttachSystem(&msystem);
        vis.SetWindowTitle("Balls SMC");
        vis.SetWindowSize(1280, 720);
        vis.SetRenderMode(opengl::WIREFRAME);
        vis.Initialize();
        vis.SetCameraPosition(ChVector<>(0, -10, 2), ChVector<>(0, 0, 0));
        vis.SetCameraVertical(CameraVerticalDir::Z);
   
    }

    int curr_step = 0;
    clock_t cpu_time = clock();
    gettimeofday(&start, NULL);

    for (curr_step = 0; curr_step < num_steps; curr_step++) {
        msystem.DoStepDynamics(time_step);

        if (useVis == true){
            vis.Render();
        }

        curr_time += time_step;

        int output_per_step = 1000; // write output to console per 1000 time step


        
        if (curr_step%output_per_step == 0){
            gettimeofday(&end, NULL);

            auto body = msystem.Get_bodylist().at(2);
            std::cout << "t = " << msystem.GetChTime() << ", KE/mgR: " << calcKE(&msystem)/(gravity * body->GetMass() * sphere_radius ) << ", simulation of " << output_per_step * time_step << "sec took " <<  time_diff(&start, &end) << std::endl;

            gettimeofday(&start, NULL);

            
        }


        if (curr_time > time_settle){
            break;
        }

    }
    msystem.GetContactContainer()->ReportAllContacts(creporter);
    // std::vector<double> normalF;
    // creporter->getNormalForce(normalF);

    // std::cout << ", F_l = " << normalF.at(0) << ", F_r = "  << normalF.at(count_Y) << std::endl;

    const std::string contact_file = "contacts_settled.csv";
    creporter->writeNormalForces(contact_file);
    

    // ************************************
    // PHASE TWO: APPLYING FORCE GRADUALLY
    // ************************************
    int F_ext_ratio_array_size = std::round(time_F_dur / time_step) + 1;

    double *F_ext_ratio_array = (double*) malloc(F_ext_ratio_array_size * sizeof(double));

    double slope = F_ext_ratio / time_F_dur;

    for (int i = 0; i < F_ext_ratio_array_size; i++) {
        F_ext_ratio_array[i] = slope * time_step * i;

    }
    int force_counter = 0;
    while (curr_time < time_settle + time_F_dur && force_counter < F_ext_ratio_array_size) {
        msystem.DoStepDynamics(time_step);


        if (useVis == true){
            vis.Render();
        }


        curr_time += time_step;

        double external_force = -F_ext_ratio_array[force_counter] * top_middle->GetMass() * gravity;


        // apply gradual force
        top_middle->Empty_forces_accumulators();
        top_middle->Accumulate_force(ChVector<double>(0, 0, external_force), top_middle->GetPos(), false);

        force_counter++;

        std::cout << "t = " << msystem.GetChTime() << ", apply force = " << F_ext_ratio_array[force_counter] <<  ", KE/mgR: " << calcKE(&msystem)/(gravity * top_middle->GetMass() * sphere_radius )<< std::endl;


    }
    free(F_ext_ratio_array);

    int frame_counter = 0;
    while(curr_time < time_end){

        // apply constant force
        top_middle->Empty_forces_accumulators();
        top_middle->Accumulate_force(ChVector<double>(0, 0, -F_ext_ratio * top_middle->GetMass() * gravity), top_middle->GetPos(), false);


        msystem.DoStepDynamics(time_step);
        if (useVis == true){
            vis.Render();
        }

        curr_time += time_step;

        frame_counter++;

        if (frame_counter%50000 == 0){

            auto creporter = chrono_types::make_shared<ContactReporter>(msystem.Get_bodylist().at(0));


            msystem.GetContactContainer()->ReportAllContacts(creporter);

            char end_filename[500];
            sprintf(end_filename, "contacts_wF_%04d.csv", int(frame_counter/1000));

            creporter->writeNormalForces(std::string(end_filename));

            std::cout << "t = " << msystem.GetChTime() << ", KE/mgR: " << calcKE(&msystem)/(gravity * top_middle->GetMass() * sphere_radius )<< std::endl;


        }

    }



    // #endif

    return 0;
}