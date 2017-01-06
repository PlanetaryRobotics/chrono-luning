// =============================================================================
// PROJECT CHRONO - http://projectchrono.org
//
// Copyright (c) 2016 projectchrono.org
// All right reserved.
//
// Use of this source code is governed by a BSD-style license that can be found
// in the LICENSE file at the top level of the distribution and at
// http://projectchrono.org/license-chrono.txt.
//
// =============================================================================
// Authors: Nic Olsen
// =============================================================================

#ifndef CHRONO_DISTRIBUTED_PHYSICS_CHBODYDISTR_H_
#define CHRONO_DISTRIBUTED_PHYSICS_CHBODYDISTR_H_

namespace chrono {

class ChBodyDistr {
public:
	ChBodyDistr();
	virtual ~ChBodyDistr();
	void SetGlobalId(int id) { if (id >= 0) global_id = id; }
	int GetGloablId() {return global_id;}
	int GetPos(int dim) {return pos[dim];}
	int GetVel(int dim) {return vel[dim];}
	void SetPos(double p, int dim) { pos[dim] = p; }
	void SetVel(double v, int dim) { vel[dim] = v; }

protected:
	double pos[3];
	double vel[3];
	double force[3];
	//TODO: Member variables for sphere, but general enough for others
	
	int global_id;
};

} /* namespace chrono */

#endif /* CHRONO_DISTRIBUTED_PHYSICS_CHBODYDISTR_H_ */
