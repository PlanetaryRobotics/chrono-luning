#include <cstdio>
#include <vector>
#include <cmath>
#include <stdio.h>
#include <ctime>
#include <sys/time.h>


#include "chrono_multicore/physics/ChSystemMulticore.h"
#include "chrono_thirdparty/filesystem/path.h"

#include "chrono/core/ChTimer.h"
#include "chrono/ChConfig.h"
#include "chrono/utils/ChUtilsCreators.h"
#include "chrono/utils/ChUtilsInputOutput.h"



using namespace chrono;


// Number of balls: (2 * count_X + 1) * (2 * count_Y + 1)
float sphere_radius = 0.01;
double sphere_mass = 0.01;
int count_X = 60;
int count_Y = 60;
int count_Z = 15;

// Material properties (same on bin and balls)
float Y = 2e8f;
// float cr = 1e-1;
float cr = 0.01;

auto top_middle = chrono_types::make_shared <ChBody>();

double rho = 7.4;
double gravity = 10;


double kn_ratio = 4e6;  // step size: 1e-6 for kn ratio 3e5
double time_step = 1e-7; // 

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


// Write to a CSV file pody position, orientation, and (optionally) linear and
// angular velocity. Optionally, only active bodies are processed.
// -----------------------------------------------------------------------------
void WriteBodies(ChSystemMulticoreSMC* system, const std::string& filename){

    utils::CSV_writer csv(",");
    csv << "x" << "y" << "z" << "vx"  << "vy" << "vz" << "wx" << "wy" << "wz" << std::endl; 

    for (auto body : system->Get_bodylist()) {
        if (body->IsActive())
        csv << body->GetPos() << body->GetPos_dt() << body->GetWvel_loc() << std::endl;
    }

    csv.write_to_file(filename);
}


// -----------------------------------------------------------------------------
// Callback class for contact reporting
// -----------------------------------------------------------------------------
class ContactReporter : public ChContactContainer::ReportContactCallback {
  public:
    ContactReporter(){
        m_csv_all << "bi" << "bj" << "Fn" << "Ft1"  << "Ft2" << "nx" << "ny" << "nz" << "x" << "y" << "z" << "dist" << std::endl; 
    }
    // compute force in global direction
    void writeContactForcesAll(const std::string& filename){        
        m_csv_all.write_to_file(filename);
        std::cout << "Successfully write all contact forces: " <<  filename  << std::endl;
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

        ChVector<double> normal = plane_coord.Get_A_Xaxis();
        m_csv_all << modA->GetPhysicsItem()->GetIdentifier() << modB->GetPhysicsItem()->GetIdentifier() << cforce  << normal << (pA+pB)/2 << distance<< std::endl;
        
        return true;
    }

    // std::shared_ptr<ChBody> m_obj1;
    // std::shared_ptr<ChBody> m_obj2;
    // std::vector<double> contact_force_z = std::vector<double>(count_Y+1);
    // utils::CSV_writer m_csv;
    utils::CSV_writer m_csv_all; // csv stream for all contacts
};


float time_diff(struct timeval *start, struct timeval *end)
{
    return (end->tv_sec - start->tv_sec) + 1e-6*(end->tv_usec - start->tv_usec);
}


// -----------------------------------------------------------------------------
// Create a bin consisting of five boxes attached to the ground.
// -----------------------------------------------------------------------------
void AddContainer(ChSystemMulticoreSMC* sys, double sphere_radius, double mu) {
    // 61 spheres per layer, total of 15 layers
    double box_hdim = (count_X + 1) * sphere_radius;
    ChVector<> hdim(box_hdim, box_hdim, count_Z * 2 * sphere_radius);
    double hthick = 2*sphere_radius;


    // Create a common material
    auto mat = chrono_types::make_shared<ChMaterialSurfaceSMC>();

    mat->SetKn(kn_ratio * sphere_mass * gravity/sphere_radius);
    mat->SetKt(kn_ratio * sphere_mass * gravity/sphere_radius);
    mat->SetFriction(mu);
    mat->SetRestitution(cr);


    // Create the containing bin (4 x 4 x 1)
    auto left_wall = std::shared_ptr<ChBody>(sys->NewBody());
    left_wall->SetMass(100);
    left_wall->SetIdentifier(908);
    // left_wall->SetPos(ChVector<>(-hdim.x() - hthick, 0, hdim.z()));
    left_wall->SetCollide(true);
    left_wall->SetBodyFixed(true);
    left_wall->GetCollisionModel()->ClearModel();
    utils::AddBoxGeometry(left_wall.get(), mat, ChVector<>(hthick, hdim.y(), hdim.z()),  ChVector<>(-hdim.x() - hthick, 0, hdim.z()));
    left_wall->GetCollisionModel()->BuildModel();

    sys->AddBody(left_wall);

    auto right_wall = std::shared_ptr<ChBody>(sys->NewBody());
    right_wall->SetMass(100);
    right_wall->SetCollide(true);
    right_wall->SetBodyFixed(true);
    right_wall->GetCollisionModel()->ClearModel();
    right_wall->SetIdentifier(909);
    utils::AddBoxGeometry(right_wall.get(), mat, ChVector<>(hthick, hdim.y(), hdim.z()),  ChVector<>( hdim.x() + hthick, 0, hdim.z()));
    right_wall->GetCollisionModel()->BuildModel();
    sys->AddBody(right_wall);



    // Create the containing bin (4 x 4 x 1)
    auto floor = std::shared_ptr<ChBody>(sys->NewBody());
    floor->SetMass(100);
    floor->SetCollide(true);
    floor->SetBodyFixed(true);
    floor->GetCollisionModel()->ClearModel();
    floor->SetIdentifier(910);
    utils::AddBoxGeometry(floor.get(), mat, ChVector<>(hdim.x(), hdim.y(), hthick),  ChVector<>(0, 0, -hthick));
    floor->GetCollisionModel()->BuildModel();
    sys->AddBody(floor);

}

// -----------------------------------------------------------------------------
// Create the falling spherical objects in a uniform rectangular grid.
// -----------------------------------------------------------------------------
void AddFallingBalls(ChSystemMulticoreSMC* sys, double sphere_radius, double mu) {
    // Common material
    // Create the falling balls
    int ballId = 0;
    double spacing = std::sqrt(3) * sphere_radius;

    auto ballMat = chrono_types::make_shared<ChMaterialSurfaceSMC>();

    ballMat->SetYoungModulus(1e7);
    ballMat->SetKn(kn_ratio * sphere_mass * gravity/sphere_radius);
    ballMat->SetKt(kn_ratio * sphere_mass * gravity/sphere_radius);
    ballMat->SetFriction(mu);
    ballMat->SetRestitution(cr);


    ChVector<> inertia = (2.0 / 5.0) * sphere_mass * sphere_radius * sphere_radius * ChVector<>(1, 1, 1);
    for (int iz = 0; iz < count_Z; iz++){
        if ( iz%2 == 0){
            for (int ix = -count_X/2; ix <= count_X/2; ix++) {
                ChVector<> pos((2*sphere_radius) * ix, 0,  sphere_radius + spacing * iz);
                auto ball = std::shared_ptr<ChBody>(sys->NewBody());
                ball->SetIdentifier(ballId++);
                ball->SetMass(sphere_mass);
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
                    ball->SetMass(sphere_mass);
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



// given output info time id, write contact force and particle positions
//int(curr_step/output_per_step)
void WrtieOutputInfo(ChSystemMulticoreSMC* sys, const std::string& subtest_dir, int time_id){
    auto creporter = chrono_types::make_shared<ContactReporter>();
    sys->GetContactContainer()->ReportAllContacts(creporter);

    char contact_force_filename[500];
    char particles_filename[500];
    sprintf(contact_force_filename, "%s/stepforce%03d.csv", subtest_dir.c_str(), time_id);
    creporter->writeContactForcesAll(contact_force_filename);
    sprintf(particles_filename, "%s/stepdata_sphere_%03d.csv", subtest_dir.c_str(), time_id);
    WriteBodies(sys, particles_filename);
}

