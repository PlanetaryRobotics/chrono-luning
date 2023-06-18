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
// Authors: Alessandro Tasora
// =============================================================================

#include "chrono/assets/ChCylinderShape.h"

namespace chrono {

// Register into the object factory, to enable run-time dynamic creation and persistence
CH_FACTORY_REGISTER(ChCylinderShape)

ChCylinderShape::ChCylinderShape() {
    SetMutable(false);
}

ChCylinderShape::ChCylinderShape(double radius, double height) {
    gcylinder.r = radius;
    gcylinder.h = height;
    SetMutable(false);
}

ChCylinderShape::ChCylinderShape(const geometry::ChCylinder& cyl) : gcylinder(cyl) {
    SetMutable(false);
}

void ChCylinderShape::ArchiveOut(ChArchiveOut& marchive) {
    // version number
    marchive.VersionWrite<ChCylinderShape>();
    // serialize parent class
    ChVisualShape::ArchiveOut(marchive);
    // serialize all member data:
    marchive << CHNVP(gcylinder);
}

void ChCylinderShape::ArchiveIn(ChArchiveIn& marchive) {
    // version number
    /*int version =*/ marchive.VersionRead<ChCylinderShape>();
    // deserialize parent class
    ChVisualShape::ArchiveIn(marchive);
    // stream in all member data:
    marchive >> CHNVP(gcylinder);
}

}  // end namespace chrono
