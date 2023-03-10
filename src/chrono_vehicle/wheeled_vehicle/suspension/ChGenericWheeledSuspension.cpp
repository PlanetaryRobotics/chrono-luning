// =============================================================================
// PROJECT CHRONO - http://projectchrono.org
//
// Copyright (c) 2023 projectchrono.org
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
// Base class for a generic wheeled vehicle suspension. Except for the spindle
// bodies and the associated revolute joints, the topology of such a suspension
// is completely arbitrary and left to derived classes.
//
// The suspension subsystem is modeled with respect to a right-handed frame,
// with X pointing towards the front, Y to the left, and Z up (ISO standard).
// The suspension reference frame is assumed to be always aligned with that of
// the vehicle.  When attached to a chassis, only an offset is provided.
//
// All point locations are assumed to be given for the left half of the
// suspension and will be mirrored (reflecting the y coordinates) to construct
// the right side.
//
// =============================================================================

#include <algorithm>

#include "chrono/assets/ChCylinderShape.h"
#include "chrono/assets/ChPointPointShape.h"

#include "chrono_vehicle/wheeled_vehicle/suspension/ChGenericWheeledSuspension.h"

namespace chrono {
namespace vehicle {

// -----------------------------------------------------------------------------

ChGenericWheeledSuspension::ChGenericWheeledSuspension(const std::string& name) : ChSuspension(name) {}

ChGenericWheeledSuspension::~ChGenericWheeledSuspension() {
    auto sys = m_spindle[0]->GetSystem();
    if (sys) {
        for (auto& body : m_bodies)
            sys->Remove(body.second.body);
        for (auto& joint : m_joints)
            ChChassis::RemoveJoint(joint.second.joint);
        for (auto& dist : m_dists)
            sys->Remove(dist.second.dist);
        for (auto& tsda : m_tsdas)
            sys->Remove(tsda.second.tsda);

        m_bodies.clear();
        m_joints.clear();
        m_dists.clear();
        m_tsdas.clear();
    }
}

// -----------------------------------------------------------------------------

void ChGenericWheeledSuspension::DefineBody(const std::string& name,
                                            bool mirrored,
                                            const ChVector<>& pos,
                                            const ChQuaternion<>& rot,
                                            double mass,
                                            const ChVector<>& inertia_moments,
                                            const ChVector<>& inertia_products,
                                            std::shared_ptr<ChVehicleGeometry> geometry) {
    Body b;
    b.body = nullptr;
    b.pos = pos;
    b.rot = rot.GetNormalized();
    b.mass = mass;
    b.inertia_moments = inertia_moments;
    b.inertia_products = inertia_products;
    b.geometry = geometry;

    if (!mirrored) {
        m_bodies.insert({{name, -1}, b});
    } else {
        m_bodies.insert({{name, 0}, b});
        m_bodies.insert({{name, 1}, b});
    }
}

void ChGenericWheeledSuspension::DefineJoint(const std::string& name,
                                             bool mirrored,
                                             ChVehicleJoint::Type type,
                                             BodyIdentifier body1,
                                             BodyIdentifier body2,
                                             const ChVector<>& pos,
                                             const ChQuaternion<>& rot,
                                             std::shared_ptr<ChVehicleBushingData> bdata) {
    Joint j;
    j.joint = nullptr;
    j.type = type;
    j.body1 = body1;
    j.body2 = body2;
    j.pos = pos;
    j.rot = rot.GetNormalized();
    j.bdata = bdata;

    if (!mirrored) {
        m_joints.insert({{name, -1}, j});
    } else {
        m_joints.insert({{name, 0}, j});
        m_joints.insert({{name, 1}, j});
    }
}

void ChGenericWheeledSuspension::DefineDistanceConstraint(const std::string& name,
                                                          bool mirrored,
                                                          BodyIdentifier body1,
                                                          BodyIdentifier body2,
                                                          const ChVector<>& point1,
                                                          const ChVector<>& point2) {
    DistanceConstraint d;
    d.dist = nullptr;
    d.body1 = body1;
    d.body2 = body2;
    d.point1 = point1;
    d.point2 = point2;

    if (!mirrored) {
        m_dists.insert({{name, -1}, d});
    } else {
        m_dists.insert({{name, 0}, d});
        m_dists.insert({{name, 1}, d});
    }
}

void ChGenericWheeledSuspension::DefineTSDA(const std::string& name,
                                            bool mirrored,
                                            BodyIdentifier body1,
                                            BodyIdentifier body2,
                                            const ChVector<>& point1,
                                            const ChVector<>& point2,
                                            double rest_length,
                                            std::shared_ptr<ChLinkTSDA::ForceFunctor> force,
                                            std::shared_ptr<ChPointPointShape> geometry) {
    TSDA t;
    t.tsda = nullptr;
    t.body1 = body1;
    t.body2 = body2;
    t.point1 = point1;
    t.point2 = point2;
    t.rest_length = rest_length;
    t.force = force;
    t.geometry = geometry;

    if (!mirrored) {
        m_tsdas.insert({{name, -1}, t});
    } else {
        m_tsdas.insert({{name, 0}, t});
        m_tsdas.insert({{name, 1}, t});
    }
}

// -----------------------------------------------------------------------------

std::string ChGenericWheeledSuspension::Name(const PartKey& id) const {
    if (id.side == -1)
        return m_name + "_" + id.name;

    return m_name + "_" + id.name + (id.side == VehicleSide::LEFT ? "_L" : "_R");
}

ChVector<> ChGenericWheeledSuspension::TransformPosition(const ChVector<>& pos_loc, int side) const {
    ChVector<> pos = pos_loc;
    if (side == VehicleSide::RIGHT)
        pos.y() = -pos.y();
    return m_X_SA.TransformLocalToParent(pos);
}

ChQuaternion<> ChGenericWheeledSuspension::TransformRotation(const ChQuaternion<>& rot_local, int side) const {
    ChQuaternion<> rot = rot_local.GetNormalized();
    if (side == VehicleSide::RIGHT) {
        ChVector<> u = rot.GetXaxis();
        ChVector<> w = rot.GetZaxis();
        u.y() *= -1;
        w.y() *= -1;
        ChVector<> v = Vcross(w, u);
        ChMatrix33<> R(u, v, w);
        rot = R.Get_A_quaternion();
    }
    return m_X_SA.GetRot() * rot;
}

std::shared_ptr<ChBody> ChGenericWheeledSuspension::FindBody(const std::string& name, int side) const {
    // Search for a non-mirrored body
    auto b1 = m_bodies.find({name, -1});
    if (b1 != m_bodies.end())
        return b1->second.body;

    // Search for a mirrored body
    auto b2 = m_bodies.find({name, side});
    assert(b2 != m_bodies.end());
    return b2->second.body;
}

void ChGenericWheeledSuspension::Initialize(std::shared_ptr<ChChassis> chassis,
                                            std::shared_ptr<ChSubchassis> subchassis,
                                            std::shared_ptr<ChSteering> steering,
                                            const ChVector<>& location,
                                            double left_ang_vel,
                                            double right_ang_vel) {
    m_parent = chassis;
    m_rel_loc = location;

    // Chassis orientation (expressed in absolute frame)
    // Recall that the suspension reference frame is aligned with the chassis.
    ChQuaternion<> chassisRot = chassis->GetBody()->GetFrame_REF_to_abs().GetRot();

    // Express the suspension reference frame in the absolute coordinate system.
    m_X_SA.SetPos(location);
    m_X_SA.ConcatenatePreTransformation(chassis->GetBody()->GetFrame_REF_to_abs());

    // Initialize all bodies in the suspension subsystem
    for (auto& item : m_bodies) {
        ChVector<> pos = TransformPosition(item.second.pos, item.first.side);
        ChQuaternion<> rot = TransformRotation(item.second.rot, item.first.side);
        item.second.body = chrono_types::make_shared<ChBody>();
        item.second.body->SetMass(item.second.mass);
        item.second.body->SetInertiaXX(item.second.inertia_moments);
        item.second.body->SetInertiaXY(item.second.inertia_products);
        item.second.body->SetPos(pos);
        item.second.body->SetRot(rot);
        item.second.body->SetNameString(Name(item.first));
        chassis->GetSystem()->AddBody(item.second.body);
    }

    // Create and initialize joints in the suspension subsystem
    for (auto& item : m_joints) {
        // Extract the two bodies connected by this joint
        std::shared_ptr<ChBody> body1;
        if (item.second.body1.chassis)
            body1 = chassis->GetBody();
        else if (item.second.body1.subchassis && subchassis != nullptr)
            body1 = subchassis->GetBeam(item.first.side == 0 ? LEFT : RIGHT);
        else if (item.second.body1.steering && steering != nullptr)
            body1 = steering->GetSteeringLink();
        else
            body1 = FindBody(item.second.body1.name, item.first.side);

        std::shared_ptr<ChBody> body2;
        if (item.second.body2.chassis)
            body2 = chassis->GetBody();
        else if (item.second.body2.subchassis && subchassis != nullptr)
            body2 = subchassis->GetBeam(item.first.side == 0 ? LEFT : RIGHT);
        else if (item.second.body2.steering && steering != nullptr)
            body2 = steering->GetSteeringLink();
        else
            body2 = FindBody(item.second.body2.name, item.first.side);

        // Create joint
        ChVector<> pos = TransformPosition(item.second.pos, item.first.side);
        ChQuaternion<> rot = TransformRotation(item.second.rot, item.first.side);
        item.second.joint = chrono_types::make_shared<ChVehicleJoint>(item.second.type,        //
                                                                      Name(item.first),        //
                                                                      body1,                   //
                                                                      body2,                   //
                                                                      ChCoordsys<>(pos, rot),  //
                                                                      item.second.bdata);
        chassis->AddJoint(item.second.joint);
    }

    // Create and initialize distance constraints in the suspension subsystem
    for (auto& item : m_dists) {
        // Extract the two bodies connected by this constraint
        std::shared_ptr<ChBody> body1;
        if (item.second.body1.chassis)
            body1 = chassis->GetBody();
        else if (item.second.body1.subchassis && subchassis != nullptr)
            body1 = subchassis->GetBeam(item.first.side == 0 ? LEFT : RIGHT);
        else if (item.second.body1.steering && steering != nullptr)
            body1 = steering->GetSteeringLink();
        else
            body1 = FindBody(item.second.body1.name, item.first.side);

        std::shared_ptr<ChBody> body2;
        if (item.second.body2.chassis)
            body2 = chassis->GetBody();
        else if (item.second.body2.subchassis && subchassis != nullptr)
            body2 = subchassis->GetBeam(item.first.side == 0 ? LEFT : RIGHT);
        else if (item.second.body2.steering && steering != nullptr)
            body2 = steering->GetSteeringLink();
        else
            body2 = FindBody(item.second.body2.name, item.first.side);

        // Create distance constraint
        ChVector<> point1 = TransformPosition(item.second.point1, item.first.side);
        ChVector<> point2 = TransformPosition(item.second.point2, item.first.side);
        item.second.dist = chrono_types::make_shared<ChLinkDistance>();
        item.second.dist->SetNameString(Name(item.first));
        item.second.dist->Initialize(body1, body2, false, point1, point2);
        chassis->GetSystem()->AddLink(item.second.dist);
    }

    // Create and initialize the spindle bodies, the spindle revolute joints, and the suspension axles
    for (int side = 0; side < 2; side++) {
        double ang_vel = (side == LEFT) ? left_ang_vel : right_ang_vel;
        double sign = (side == LEFT) ? -1 : +1;

        auto spindlePos = TransformPosition(getSpindlePos(), side);
        auto spindleRot = chassisRot * Q_from_AngZ(sign * getToeAngle()) * Q_from_AngX(sign * getCamberAngle());

        // Spindle body
        m_spindle[side] = std::shared_ptr<ChBody>(chassis->GetSystem()->NewBody());
        m_spindle[side]->SetNameString(Name({"spindle", side}));
        m_spindle[side]->SetPos(spindlePos);
        m_spindle[side]->SetRot(spindleRot);
        m_spindle[side]->SetWvel_loc(ChVector<>(0, ang_vel, 0));
        m_spindle[side]->SetMass(getSpindleMass());
        m_spindle[side]->SetInertiaXX(getSpindleInertia());
        chassis->GetSystem()->AddBody(m_spindle[side]);

        // Spindle revolute joint
        auto abody_id = getSpindleAttachmentBody();
        std::shared_ptr<ChBody> abody;
        if (abody_id.chassis)
            abody = chassis->GetBody();
        else
            abody = FindBody(abody_id.name, side);

        ChCoordsys<> rev_csys(spindlePos, spindleRot * Q_from_AngX(CH_C_PI_2));
        m_revolute[side] = chrono_types::make_shared<ChLinkLockRevolute>();
        m_revolute[side]->SetNameString(Name({"spindle_rev", side}));
        m_revolute[side]->Initialize(m_spindle[side], abody, rev_csys);
        chassis->GetSystem()->AddLink(m_revolute[side]);

        // Axle shaft
        m_axle[side] = chrono_types::make_shared<ChShaft>();
        m_axle[side]->SetNameString(Name({"axle", side}));
        m_axle[side]->SetInertia(getAxleInertia());
        m_axle[side]->SetPos_dt(-ang_vel);
        chassis->GetSystem()->AddShaft(m_axle[side]);

        // Axle connection to spindle
        m_axle_to_spindle[side] = chrono_types::make_shared<ChShaftsBody>();
        m_axle_to_spindle[side]->SetNameString(Name({"axle_to_spindle", side}));
        m_axle_to_spindle[side]->Initialize(m_axle[side], m_spindle[side], ChVector<>(0, -1, 0));
        chassis->GetSystem()->Add(m_axle_to_spindle[side]);
    }

    // Create and initialize TSDAs in the suspension subsystem
    for (auto& item : m_tsdas) {
        // Extract the two bodies connected by this joint
        std::shared_ptr<ChBody> body1;
        if (item.second.body1.chassis)
            body1 = chassis->GetBody();
        else if (item.second.body1.subchassis && subchassis != nullptr)
            body1 = subchassis->GetBeam(item.first.side == 0 ? LEFT : RIGHT);
        else if (item.second.body1.steering && steering != nullptr)
            body1 = steering->GetSteeringLink();
        else
            body1 = FindBody(item.second.body1.name, item.first.side);

        std::shared_ptr<ChBody> body2;
        if (item.second.body2.chassis)
            body2 = chassis->GetBody();
        else if (item.second.body2.subchassis && subchassis != nullptr)
            body2 = subchassis->GetBeam(item.first.side == 0 ? LEFT : RIGHT);
        else if (item.second.body2.steering && steering != nullptr)
            body2 = steering->GetSteeringLink();
        else
            body2 = FindBody(item.second.body2.name, item.first.side);

        // Create TSDA
        ChVector<> point1 = TransformPosition(item.second.point1, item.first.side);
        ChVector<> point2 = TransformPosition(item.second.point2, item.first.side);
        item.second.tsda = chrono_types::make_shared<ChLinkTSDA>();
        item.second.tsda->SetNameString(Name(item.first));
        item.second.tsda->Initialize(body1, body2, false, point1, point2);
        item.second.tsda->SetRestLength(item.second.rest_length);
        item.second.tsda->RegisterForceFunctor(item.second.force);
        chassis->GetSystem()->AddLink(item.second.tsda);
    }
}

void ChGenericWheeledSuspension::InitializeInertiaProperties() {
    m_mass = 2 * getSpindleMass();
    for (const auto& item : m_bodies)
        m_mass += item.second.body->GetMass();
}

void ChGenericWheeledSuspension::UpdateInertiaProperties() {
    m_parent->GetTransform().TransformLocalToParent(ChFrame<>(m_rel_loc, QUNIT), m_xform);

    // Calculate COM and inertia expressed in global frame
    utils::CompositeInertia composite;

    ChMatrix33<> inertiaSpindle(getSpindleInertia());
    composite.AddComponent(m_spindle[LEFT]->GetFrame_COG_to_abs(), getSpindleMass(), inertiaSpindle);
    composite.AddComponent(m_spindle[RIGHT]->GetFrame_COG_to_abs(), getSpindleMass(), inertiaSpindle);
    for (const auto& item : m_bodies)
        composite.AddComponent(item.second.body->GetFrame_COG_to_abs(), item.second.body->GetMass(),
                               item.second.body->GetInertia());

    // Express COM and inertia in subsystem reference frame
    m_com.coord.pos = m_xform.TransformPointParentToLocal(composite.GetCOM());
    m_com.coord.rot = QUNIT;

    m_inertia = m_xform.GetA().transpose() * composite.GetInertia() * m_xform.GetA();
}

// -----------------------------------------------------------------------------

double ChGenericWheeledSuspension::GetTrack() {
    return 2 * getSpindlePos().y();
}

std::shared_ptr<ChBody> ChGenericWheeledSuspension::GetAntirollBody(VehicleSide side) const {
    auto abody_id = getAntirollBody();
    return FindBody(abody_id.name, side);
}

// -----------------------------------------------------------------------------

std::vector<ChSuspension::ForceTSDA> ChGenericWheeledSuspension::ReportSuspensionForce(VehicleSide side) const {
    std::vector<ChSuspension::ForceTSDA> forces;

    for (const auto& item : m_tsdas) {
        if (item.first.side != side)
            continue;
        auto tsda = item.second.tsda;
        ////std::cout << "TSDA " << item.first.name << std::endl;
        ////std::cout << "  force:    " << tsda->GetForce() << std::endl;
        ////std::cout << "  length:   " << tsda->GetLength() << std::endl;
        ////std::cout << "  velocity: " << tsda->GetVelocity() << std::endl;
        forces.push_back(
            ChSuspension::ForceTSDA(item.first.name, tsda->GetForce(), tsda->GetLength(), tsda->GetVelocity()));
    }

    return forces;
}

void ChGenericWheeledSuspension::LogConstraintViolations(VehicleSide side) {
    // Spindle revolute joint
    {
        const auto& C = m_revolute[side]->GetConstraintViolation();
        std::cout << "Spindle revolute " << std::endl;
        std::cout << "   " << C.transpose() << std::endl;
    }

    // Suspension joints
    for (const auto& item : m_joints) {
        if (item.first.side != side)
            continue;
        auto joint = item.second.joint;
        if (!joint->IsKinematic())
            continue;
        auto link = joint->GetAsLink();
        const auto& C = link->GetConstraintViolation();
        assert(C.size() == link->GetDOC_c());
        std::cout << "Joint " << item.first.name << " type: " << ChVehicleJoint::GetTypeString(item.second.type)
                  << std::endl;
        std::cout << "   " << C.transpose() << std::endl;
    }

    // Distance constraints
    for (const auto& item : m_dists) {
        if (item.first.side != side)
            continue;
        auto dist = item.second.dist;
        const auto& C = dist->GetConstraintViolation();
        std::cout << "Distance constraint " << item.first.name << std::endl;
        std::cout << "   " << C.transpose() << std::endl;
    }
}

// -----------------------------------------------------------------------------

void ChGenericWheeledSuspension::AddVisualizationAssets(VisualizationType vis) {
    ChSuspension::AddVisualizationAssets(vis);

    if (vis == VisualizationType::NONE)
        return;

    for (const auto& item : m_bodies)
        if (item.second.geometry)
            item.second.geometry->CreateVisualizationAssets(item.second.body, vis);
    for (const auto& item : m_tsdas)
        if (item.second.geometry)
            item.second.tsda->AddVisualShape(item.second.geometry);
    for (const auto& item : m_dists)
        item.second.dist->AddVisualShape(chrono_types::make_shared<ChSegmentShape>());
}

void ChGenericWheeledSuspension::RemoveVisualizationAssets() {
    for (const auto& item : m_bodies)
        ChPart::RemoveVisualizationAssets(item.second.body);
    for (const auto& item : m_tsdas)
        ChPart::RemoveVisualizationAssets(item.second.tsda);
    for (const auto& item : m_dists)
        ChPart::RemoveVisualizationAssets(item.second.dist);

    ChSuspension::RemoveVisualizationAssets();
}

// -----------------------------------------------------------------------------

void ChGenericWheeledSuspension::ExportComponentList(rapidjson::Document& jsonDocument) const {
    ChPart::ExportComponentList(jsonDocument);

    std::vector<std::shared_ptr<ChBody>> bodies;
    bodies.push_back(m_spindle[0]);
    bodies.push_back(m_spindle[1]);
    for (const auto& item : m_bodies)
        bodies.push_back(item.second.body);
    ExportBodyList(jsonDocument, bodies);

    std::vector<std::shared_ptr<ChShaft>> shafts;
    shafts.push_back(m_axle[0]);
    shafts.push_back(m_axle[1]);
    ExportShaftList(jsonDocument, shafts);

    std::vector<std::shared_ptr<ChLink>> joints;
    std::vector<std::shared_ptr<ChLoadBodyBody>> bushings;
    joints.push_back(m_revolute[0]);
    joints.push_back(m_revolute[1]);
    for (const auto& item : m_joints)
        item.second.joint->IsKinematic() ? joints.push_back(item.second.joint->GetAsLink())
                                         : bushings.push_back(item.second.joint->GetAsBushing());
    for (const auto& item : m_dists)
        joints.push_back(item.second.dist);
    ExportJointList(jsonDocument, joints);
    ExportBodyLoadList(jsonDocument, bushings);

    std::vector<std::shared_ptr<ChLinkTSDA>> springs;
    for (const auto& item : m_tsdas)
        springs.push_back(item.second.tsda);
    ExportLinSpringList(jsonDocument, springs);
}

void ChGenericWheeledSuspension::Output(ChVehicleOutput& database) const {
    if (!m_output)
        return;

    std::vector<std::shared_ptr<ChBody>> bodies;
    bodies.push_back(m_spindle[0]);
    bodies.push_back(m_spindle[1]);
    for (const auto& item : m_bodies)
        bodies.push_back(item.second.body);
    database.WriteBodies(bodies);

    std::vector<std::shared_ptr<ChShaft>> shafts;
    shafts.push_back(m_axle[0]);
    shafts.push_back(m_axle[1]);
    database.WriteShafts(shafts);

    std::vector<std::shared_ptr<ChLink>> joints;
    std::vector<std::shared_ptr<ChLoadBodyBody>> bushings;
    joints.push_back(m_revolute[0]);
    joints.push_back(m_revolute[1]);
    for (const auto& item : m_joints)
        item.second.joint->IsKinematic() ? joints.push_back(item.second.joint->GetAsLink())
                                         : bushings.push_back(item.second.joint->GetAsBushing());
    for (const auto& item : m_dists)
        joints.push_back(item.second.dist);
    database.WriteJoints(joints);
    database.WriteBodyLoads(bushings);

    std::vector<std::shared_ptr<ChLinkTSDA>> springs;
    for (const auto& item : m_tsdas)
        springs.push_back(item.second.tsda);
    database.WriteLinSprings(springs);
}

// -----------------------------------------------------------------------------

ChGenericWheeledSuspension::BodyIdentifier::BodyIdentifier(const std::string& part_name,
                                                           bool is_chassis,
                                                           bool is_subchassis,
                                                           bool is_steering)
    : name(part_name), chassis(is_chassis), subchassis(is_subchassis), steering(is_steering) {}

ChGenericWheeledSuspension::BodyIdentifier::BodyIdentifier()
    : name(""), chassis(false), subchassis(false), steering(false) {}

ChGenericWheeledSuspension::ChassisIdentifier::ChassisIdentifier() : BodyIdentifier("", true, false, false) {}
ChGenericWheeledSuspension::SubchassisIdentifier::SubchassisIdentifier() : BodyIdentifier("", false, true, false) {}
ChGenericWheeledSuspension::SteeringIdentifier::SteeringIdentifier() : BodyIdentifier("", false, false, true) {}

bool ChGenericWheeledSuspension::PartKey::operator==(const PartKey& rhs) const {
    return name == rhs.name && side == rhs.side;
}

std::size_t ChGenericWheeledSuspension::PartKeyHash::operator()(const PartKey& id) const {
    /*
    std::string ext_name = id.name;
    if (id.side != -1)
        ext_name += (id.side == VehicleSide::LEFT ? "_L" : "_R");
    return std::hash<std::string>{}(ext_name);
    */
    std::size_t h1 = std::hash<std::string>{}(id.name);
    if (id.side == -1)
        return h1;
    std::string ext = (id.side == VehicleSide::LEFT ? "_L" : "_R");
    std::size_t h2 = std::hash<std::string>{}(id.name);
    return h1 ^ (h2 << 1);
}

}  // end namespace vehicle
}  // end namespace chrono
