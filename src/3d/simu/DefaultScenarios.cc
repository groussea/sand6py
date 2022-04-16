/*
 * This file is part of Sand6, a C++ continuum-based granular simulator.
 *
 * Copyright 2016 Gilles Daviet <gilles.daviet@inria.fr> (Inria - Université Grenoble Alpes)
 *
 * Sand6 is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.

 * Sand6 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.

 * You should have received a copy of the GNU General Public License
 * along with Sand6.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "Scenario.hh"

#include "Simu.hh"
#include "RigidBody.hh"

#include "geo/LevelSet.hh"

#include "utils/Config.hh"
#include "utils/string.hh"
#include "utils/Log.hh"

namespace d6
{

// Default scenars

struct RayleighScenar : public Scenario
{
	Scalar particle_density(const Vec &x) const override
	{
		return (x[2] > .5 * m_config->box[2]) ? 1. : 0.;
	}
};

struct BedScenar : public Scenario
{
	Scalar particle_density(const Vec &x) const override
	{
		return (x[2] < .5 * m_config->box[2]) ? 1. : 0.;
	}
};

struct CollapseScenar : public Scenario
{
	Scalar particle_density(const Vec &x) const override
	{
		return (x[0] < 0.5 * m_config->box[0] || x[2] < 1. * m_config->box[2]) ? 1. : 0.;
	}

};

struct BridsonScenar : public Scenario
{
	Vec center;
	Scalar radius;

	virtual void init(const Params &) override
	{
		center = Vec(.75 * m_config->box[0], .75 * m_config->box[1], .5 * m_config->box[2]);
		radius = .125 * m_config->box[0];
	}

	Scalar particle_density(const Vec &x) const override
	{
		return (x[0] > center[0] && x[1] > center[1] && (x - center).squaredNorm() > radius * radius) ? 1. : 0.;
	}
};

struct TowerScenar : public Scenario
{
	Vec center;
	Scalar radius;

	Scalar volMass;

	Scalar h0;
	Scalar hvel;
	Scalar zvel;
	Scalar avel;

	Scalar particle_density(const Vec &x) const override
	{
		return (
				   (std::fabs(x[0] - center[0]) < radius / 2 && std::fabs(x[1] - center[1]) < radius) || (std::fabs(x[0] - center[0]) < radius && std::fabs(x[1] - center[1]) < radius / 2)
				   //std::pow( x[1] - center[1], 2 ) + std::pow( x[0] - center[0], 2 ) < radius * radius
				   || std::fabs(x[2]) < .25 * radius)
				   ? 1.
				   : 0.;
	}

	virtual void init(const Params &params) override
	{
		center = Vec(.5 * m_config->box[0], .5 * m_config->box[1], .5 * m_config->box[2]);
		radius = .125 * m_config->box[0];

		volMass = scalar_param(params, "vm", Units::VolumicMass, 1.5 * m_config->units().R);
		h0 = scalar_param(params, "h0", Units::Length, 1.);
		hvel = scalar_param(params, "hvel", Units::Velocity, 1.);
		avel = scalar_param(params, "avel", Units::Frequency, 0.);

		const Scalar t = h0 / hvel;
		zvel = .5 * m_config->gravity.norm() * t;
	}

	void add_rigid_bodies(std::vector<RigidBody> &rbs) const override
	{
		LevelSet::Ptr ls = LevelSet::make_sphere();
		Vec pos = center - h0 * Vec(1, 1, 0) / M_SQRT2;
		pos[2] -= .125 * m_config->box[2];
		ls->scale(.5 * .125 * m_config->box[0]).set_origin(pos);

		rbs.emplace_back(ls, volMass);
		rbs.back().set_velocity(Vec(hvel / M_SQRT2, hvel / M_SQRT2, zvel), Vec(avel, 0, 0));
	}

	void update(Simu &simu, Scalar /*time*/, Scalar dt) const override
	{
		for (RigidBody &rb : simu.rigidBodies())
		{
			rb.integrate_gravity(dt, m_config->gravity);
		}
	}
};

struct RbPlaneTestScenar : public Scenario
{
	// Scalar particle_density(const Vec &x) const override
	// {
	// 	return (x[2] < .5 * m_config->box[2]) ? 1. : 0.;
	// }

		Scalar particle_density(const Vec &x) const override
	{
		return ((x[0] < 0.9 * m_config->box[0] || x[2] < 0. * m_config->box[2]) && (x[2] < 0.5* m_config->box[2])) ? 1. : 0.;
	}

	void init(const Params &params) override
	{

		R = scalar_param(params, "d", Units::None, 0.05) * m_config->box[0];
		H = scalar_param(params, "h", Units::None, 0.5) * m_config->box[2];
		S = scalar_param(params, "s", Units::None, 1.);


	}


	void add_rigid_bodies(std::vector<RigidBody> &rbs) const override
	{
		Vec box = m_config->box;

		const Scalar L = box[2] * (0.9);
		const Scalar W = 2 * m_config->typicalLength();
		Vec doorBox = Vec(W, 0.6 * box[1], 0.5 * L);

		LevelSet::Ptr ls = LevelSet::make_box(doorBox);
		ls->set_origin(Vec(0.9*box[0] + W, box[1] * 0.5, -0.1 * box[2] - L * 0.5));
		rbs.emplace_back(ls, 1.e99);


	}

	void update(Simu &simu, Scalar time, Scalar dt) const override
	{
		Vec box = m_config->box;
		const Scalar L = box[2] * (0.9);
		const Scalar W = 2 * m_config->typicalLength();
		if (time*m_config->units().toSI( Units::Time )>4){
					for (RigidBody &rb : simu.rigidBodies())
		{
			if (time*m_config->units().toSI( Units::Time )<6){
					rb.move_to(Vec(0.9*box[0] + W, box[1] * 0.5, 0.4 * box[2] + L * 0.5));
			} else if(time*m_config->units().toSI( Units::Time )>6 and time*m_config->units().toSI( Units::Time )<6.1) {
				rb.move_to(Vec(1*box[0] - W, box[1] * 0.5, 0.1 * box[2] + L * 0.5));
			}
		rb.set_velocity(Vec(0,0,0.5*m_config->units().toSI( Units::Velocity )), Vec::Zero());
		}

		} 
	}


Eigen::Matrix<Scalar, 2, Eigen::Dynamic> bezier;
	Scalar R;
	Scalar H;
	Scalar S;

};

struct ImpactScenarlhe : public Scenario
{
	virtual void init(const Params &params)
	{
		volMass = scalar_param(params, "vm", Units::VolumicMass, 1.5 * m_config->units().R);
		d = scalar_param(params, "d", Units::Length, 0.25);
		hBall = scalar_param(params, "h_ball", Units::Length, 1. / 3);
		hBed = scalar_param(params, "h_bed", Units::Length, 1. / 3);
	}

	Scalar particle_density(const Vec &x) const override
	{
		Scalar w_sideW=m_config->units().fromSI(Units::Length) * 0.01;
		return (x[2] < hBed && (x[1]<(m_config->box[1]-w_sideW))&& (x[1]>w_sideW) && (x[0]<(m_config->box[0]-w_sideW)) && (x[0]>w_sideW)) ? 1. : 0.;
	}

	void add_rigid_bodies(std::vector<RigidBody> &rbs) const override
	{
		Vec box = m_config->box;
		Scalar w_sideW=m_config->units().fromSI(Units::Length) * 0.01;
		LevelSet::Ptr ls = LevelSet::make_sphere();
		ls->scale(radius()).set_origin(Vec(0.5 * m_config->box[0], 0.5 * m_config->box[1], 1.5*m_config->box[2]));
		// ls->scale(radius()).set_origin(Vec(0.5 * m_config->box[0], 0.5 * m_config->box[1], hBall+hBed));
		// determine the fall velocity when the ball enters in the domain
		rbs.emplace_back(ls, volMass);
		rbs.back().set_velocity(Vec(0, 0, -std::sqrt(-2 * m_config->gravity[2] * (hBall - (1.5*m_config->box[2] - hBed)))), Vec(0, 0, 0));
		// rbs.back().set_velocity(Vec(0, 0,0), Vec(0, 0, 0));
		
		std::ofstream RBout(m_config->base_dir + "/RBinfo.txt");
		d6::dump(RBout, "Nframe,time,X,Y,Z,Ux,Uy,Uz\n");
		// introduce walls

		Vec wallBox = Vec(box[0], w_sideW/2, box[2]);
		auto leftWall = LevelSet::make_box(wallBox);
		leftWall->set_origin(Vec(box[0] * 0.5, w_sideW/2 ,box[2] * 0.5));
		rbs.emplace_back(leftWall, 1.e99);
		auto rightWall = LevelSet::make_box(wallBox);
		rightWall->set_origin(Vec(box[0] * 0.5, box[1]- w_sideW/2  ,box[2] * 0.5));
		rbs.emplace_back(rightWall, 1.e99);

		Vec bWallBox = Vec( w_sideW/2, box[1],box[2]);

		auto backWall = LevelSet::make_box(bWallBox);
		backWall->set_origin(Vec(w_sideW/2 , 0.5*box[1],box[2] * 0.5));
		rbs.emplace_back(backWall, 1.e99);

		auto frontWall = LevelSet::make_box(bWallBox);
		frontWall->set_origin(Vec(box[0]- w_sideW/2 , 0.5*box[1],box[2] * 0.5));
		rbs.emplace_back(frontWall, 1.e99);


	}

	void update(Simu &simu, Scalar time, Scalar dt) const override
	{
		int n =0 ;
		for (RigidBody &rb : simu.rigidBodies())
		{
			if (n==0){
			rb.integrate_gravity(dt, m_config->gravity);
			Vec vel = rb.velocity();
			rb.set_velocity(vel, rb.angularVelocity());
			std::ofstream RBout(m_config->base_dir + "/RBinfo.txt", std::fstream::app);
			RBout << std::round(time * m_config->fps)<<",";
			RBout << time * m_config->units().toSI(Units::Time)<<",";
			RBout << rb.levelSet().origin()[0] * m_config->units().toSI(Units::Length)<<",";
			RBout << rb.levelSet().origin()[1] * m_config->units().toSI(Units::Length)<<",";
			RBout << rb.levelSet().origin()[2] * m_config->units().toSI(Units::Length)<<",";
			RBout << vel[0] * m_config->units().toSI(Units::Velocity)<<",";
			RBout << vel[1] * m_config->units().toSI(Units::Velocity)<<",";
			RBout << vel[2] * m_config->units().toSI(Units::Velocity)<<"\n";
			RBout.close();
			n = n+1;}
		}
	}

	Scalar radius() const
	{
		return d / 2;
	}


private:
	Scalar volMass;
	Scalar d;
	Scalar hBed;
	Scalar hBall;
};

struct ImpactScenar : public Scenario
{

	Scalar particle_density(const Vec &x) const override
	{
		return (x[2] < h * m_config->box[2]) ? 1. : 0.;
	}

	virtual void init(const Params &params)
	{
		volMass = scalar_param(params, "vm", Units::VolumicMass, 1.5 * m_config->units().R);
		zvel = scalar_param(params, "zvel", Units::Velocity, 0.);
		avel = scalar_param(params, "avel", Units::Frequency, 0.);
		d = scalar_param(params, "d", Units::None, 0.25);
		h = scalar_param(params, "h", Units::None, 1. / 3);
	}

	void add_rigid_bodies(std::vector<RigidBody> &rbs) const override
	{
		LevelSet::Ptr ls = LevelSet::make_sphere();
		ls->scale(radius()).set_origin(.5 * m_config->box + Vec(0, 0, .25 * m_config->box[2]));

		rbs.emplace_back(ls, volMass);
		rbs.back().set_velocity(Vec(0, 0, -zvel), Vec(avel, 0, 0));
	}

	void update(Simu &simu, Scalar /*time*/, Scalar dt) const override
	{
		for (RigidBody &rb : simu.rigidBodies())
		{
			rb.integrate_gravity(dt, m_config->gravity);
			Vec vel = rb.velocity();
			vel[2] = std::max(vel[2], (radius() - rb.levelSet().origin()[2]) / dt);
			rb.set_velocity(vel, rb.angularVelocity());
		}
	}

	Scalar radius() const
	{
		return d / 2 * m_config->box[0];
	}

private:
	Scalar volMass;
	Scalar zvel;
	Scalar avel;
	Scalar d;
	Scalar h;
};

struct WheelScenar : public Scenario
{

	Scalar particle_density(const Vec &x) const override
	{
		return (x[2] < 1. / 2. * m_config->box[2]) ? 1. : 0.;
	}

	virtual void init(const Params &params)
	{
		volMass = scalar_param(params, "vm", Units::VolumicMass, 1.5 * m_config->units().R);
		xvel = scalar_param(params, "xvel", Units::Velocity, 0.);
		avel = scalar_param(params, "avel", Units::Frequency, 0.);
		torque = scalar_param(params, "torque", Units::Torque, 0.);
		tilt = scalar_param(params, "tilt", Units::None, 0.);
		R = scalar_param(params, "R", Units::None, 0.125);
		r = scalar_param(params, "r", Units::None, 0.375);
		angFric = scalar_param(params, "af", Units::None, 0) * m_config->units().fromSI(Units::Torque) * m_config->units().fromSI(Units::Time);
	}

	void add_rigid_bodies(std::vector<RigidBody> &rbs) const override
	{
		LevelSet::Ptr ls = LevelSet::make_torus(r);
		ls->scale(radius())
			.rotate(Vec(1, 0, 0), (1 - tilt) * M_PI_2)
			.set_origin(.5 * m_config->box + Vec(.25 * m_config->box[0], -.25 * m_config->box[1], radius() * (1 + 2 * r)));

		rbs.emplace_back(ls, volMass);
		rbs.back().set_velocity(Vec(0, 0, -xvel), Vec(0, -avel, 0));
	}

	void update(Simu &simu, Scalar /*time*/, Scalar dt) const override
	{
		Vec6 forces;
		forces.head<3>().setZero();

		for (RigidBody &rb : simu.rigidBodies())
		{
			forces.tail<3>() = rb.levelSet().rotation() * Vec(0, 0, torque - angFric * rb.angularVelocity().norm());
			rb.integrate_gravity(dt, m_config->gravity);
			rb.integrate_forces(dt, forces);
		}
	}

	Scalar radius() const
	{
		return R * m_config->box[0];
	}

private:
	Scalar volMass;
	Scalar xvel;
	Scalar avel;
	Scalar torque;
	Scalar tilt;
	Scalar R;
	Scalar r;
	Scalar angFric;
};

struct SiloScenar : public Scenario
{

	constexpr static Scalar s_h = 2;

	Scalar particle_density(const Vec &x) const override
	{
		return (x[2] > (.5 * m_config->box[2] + s_h) //(1-Dbar)*R )
													 //	|| x[2] < h*m_config->box[2]
				)
				   ? 1.
				   : 0.;
	}

	virtual void init(const Params &params) override
	{
		W = .5 * m_config->box[0];
		R = scalar_param(params, "d", Units::None, 0.5) * W;
	}

	void add_rigid_bodies(std::vector<RigidBody> &rbs) const override
	{
		LevelSet::Ptr ls = LevelSet::make_hole(R / s_h + 1);
		ls->scale(s_h).set_origin(.5 * m_config->box);
		rbs.emplace_back(ls, 1.e99);
	}

	void update(Simu &simu, Scalar /*time*/, Scalar /*dt*/) const override
	{
		DynParticles &particles = simu.particles();
		const Scalar zmin = m_config->box[2] / 3;

#pragma omp parallel for
		for (size_t i = 0; i < particles.count(); ++i)
		{
			if (particles.geo().centers()(2, i) < zmin)
			{
				particles.remove(i);
			}
		}
	}

private:
	Scalar W;
	Scalar R;
};

struct HourglassScenar : public Scenario
{

	Scalar particle_density(const Vec &x) const override
	{
		return (hg_ls2->eval_at(x) < 0 &&
				x[2] < 0.3 * m_config->box[2] && x[0])
				   ? 1
				   : 0;
	// 		return (
	// 		x[2] > .5 * m_config->box[2] && x[2] < .8 * m_config->box[2])
	// 			? 1
	// 			: 0;
	}

	void init(const Params &params) override
	{
		const Scalar D = scalar_param(params, "L", Units::None, 0.75);
		const Scalar d = scalar_param(params, "d", Units::None, 0.25);
		const Scalar S = .5 * D * m_config->box[0];
		const Scalar H = m_config->box[0] / 2 / S;


		hg_ls = LevelSet::make_hourglass(H, d);

		hg_ls->scale(S).rotate(Vec(1, 0, 0), M_PI_2).set_origin(.5 * m_config->box);

		hg_ls2 = LevelSet::make_hourglass(H, d);

		hg_ls2->scale(S).rotate(Vec(1, 0, 0), M_PI_2).set_origin(.5 * m_config->box);

		avel = scalar_param(params, "avel", Units::Frequency, 0.062);
	}

	void add_rigid_bodies(std::vector<RigidBody> &rbs) const override
	{
		rbs.emplace_back(hg_ls, 1.e99);
		rbs.back().set_velocity(Vec(0, 0, 0), Vec(0,avel , 0));
	}


	Scalar avel;
	mutable LevelSet::Ptr hg_ls;
	mutable LevelSet::Ptr hg_ls2;
};

struct BunnyScenar : public Scenario
{

	Scalar particle_density(const Vec &x) const override
	{
		// std::cout<<"init"<<std::endl;

		return (bunny_ls2->eval_at(x) < 0 &&
				x[2] * 5 < m_config->box[0])
				   ? 1
				   : 0;
	}

	void init(const Params &params) override
	{
		std::cout<<"init"<<std::endl;
		const std::string &meshname = string_param(params, "mesh", "../scenes/bunny.obj");
		// const std::string &meshname = string_param(params, "mesh", "../scenes/testVall.obj");
		bunny_ls = LevelSet::from_mesh(meshname.c_str());


		const Scalar S = m_config->box[0];

		bunny_ls->scale(S )
			.rotate(Vec(0, 1, 0), M_PI / 4)
			.rotate(Vec(1, 0, 0), M_PI / 4)
			.set_origin(S * Vec(.5, .25, -.5));
		bunny_ls->compute();

		const std::string &meshname2 =  "../scenes/bunny.obj";
		// const std::string &meshname = string_param(params, "mesh", "../scenes/testVall.obj");
		bunny_ls2 = LevelSet::from_mesh(meshname2.c_str());


		bunny_ls2->scale(S * 4)
			.rotate(Vec(0, 1, 0), M_PI / 4)
			.rotate(Vec(1, 0, 0), M_PI / 4)
			.set_origin(S * Vec(.5, .25, -.5));
		bunny_ls2->compute();


	}

	void add_rigid_bodies(std::vector<RigidBody> &rbs) const override
	{
		
		rbs.push_back(RigidBody{bunny_ls, 1.e99});
	}

	void update(Simu &simu, Scalar time, Scalar /*dt*/) const override
	{
		// const Scalar speed = m_config->units().toSI(Units::Time);

		// const Scalar tw = time * speed;
		// Vec vel = Vec::Zero();

		// if (tw > .5)
		// {
		// 	const Scalar t = tw - .5;
		// 	if (t < 1)
		// 	{
		// 		vel = Vec(0, 0, 6 * (t - t * t)) * speed * (.3 * m_config->box[0]);
		// 	}
		// }

		for (RigidBody &rb : simu.rigidBodies())
		{
			rb.set_velocity(vel, Vec::Zero());
		}
	}

	mutable LevelSet::Ptr bunny_ls;
};

struct WritingScenar : public Scenario
{

	Scalar particle_density(const Vec &x) const override
	{
		return (x[2] < .5 * m_config->box[2]) ? 1 : 0;
	}

	void init(const Params &params) override
	{

		R = scalar_param(params, "d", Units::None, 0.05) * m_config->box[0];
		H = scalar_param(params, "h", Units::None, 0.5) * m_config->box[2];
		S = scalar_param(params, "s", Units::None, 1.);

		bezier.resize(2, 12);
		bezier << 0.196, 1.046, -0.358, 0.492, 0.609, 0.858, 0.677, 0.737, 1.068, 0.518, 1.630, 0.990,
			0.140, 0.269, 0.450, 0.561, 0.482, 0.686, 0.630, 0.106, 0.574, -0.25, 0.231, 0.364;

		bezier *= m_config->box[0] / 1.4;
	}

	void add_rigid_bodies(std::vector<RigidBody> &rbs) const override
	{
		LevelSet::Ptr cyl_ls;
		cyl_ls = LevelSet::make_cylinder(H / R);
		cyl_ls->scale(R).set_origin(Vec::Zero());
		rbs.emplace_back(cyl_ls, 1.e99);
	}

	void update(Simu &simu, Scalar time, Scalar /*dt*/) const override
	{
		const unsigned nCurves = bezier.cols() / 4;
		const Scalar speed = m_config->units().toSI(Units::Time) * S;
		const Scalar trans = .1;

		const Scalar tw = time * speed;
		const unsigned cid = std::floor(tw / (1. + trans));

		Scalar t = tw - (1. + trans) * cid;
		const bool inbetween = t < trans;
		t -= trans;

		Vec vel = Vec::Zero();

		if (cid < nCurves && !inbetween)
		{
			vel.head<2>() = -3 * (1 - t) * (1 - t) * bezier.col(0 + 4 * cid) + 3 * ((1 - t) - 2 * t) * (1 - t) * bezier.col(1 + 4 * cid) + 3 * (2 * (1 - t) - t) * t * bezier.col(2 + 4 * cid) + 3 * t * t * bezier.col(3 + 4 * cid);
			const Scalar S = 5.e1;
			vel[2] = -(H / 2 - R) * 2 * S * t * std::exp(-S * t * t);
			vel *= speed;
		}

		for (RigidBody &rb : simu.rigidBodies())
		{
			if (cid < nCurves && inbetween)
			{
				Vec origin;
				origin.head<2>() = bezier.col(4 * cid);
				origin[2] = .5 * m_config->box[2] + H / 2;
				rb.move_to(origin);
			}

			rb.set_velocity(vel, Vec::Zero());
		}
	}

	Eigen::Matrix<Scalar, 2, Eigen::Dynamic> bezier;
	Scalar R;
	Scalar H;
	Scalar S;
};

struct DiggingScenar : public Scenario
{

	Scalar particle_density(const Vec &x) const override
	{
		return (x[2] < .5 * m_config->box[2]) ? 1 : 0;
	}

	void init(const Params &params) override
	{
		speedup = scalar_param(params, "s", Units::None, .5);

		std::string meshname = string_param(params, "mesh", "../scenes/Schmuck2l.obj");
		left_hand_ls = LevelSet::from_mesh(meshname.c_str());
		meshname = string_param(params, "mesh", "../scenes/Schmuck2r.obj");
		right_hand_ls = LevelSet::from_mesh(meshname.c_str());

		const Scalar S = m_config->box[0];
		const Scalar H = m_config->box[2];

		left_hand_ls->scale(.5 * S)
			.rotate(Vec(0, 1, 0), M_PI_2)
			.set_origin(Vec(S / 2, 0, H));
		left_hand_ls->compute();

		right_hand_ls->scale(.5 * S)
			.rotate(Vec(0, 1, 0), -M_PI_2)
			.rotate(Vec(0, 0, 1), -M_PI)
			.set_origin(Vec(S / 2, S, H));
		right_hand_ls->compute();
	}

	void add_rigid_bodies(std::vector<RigidBody> &rbs) const override
	{
		rbs.emplace_back(left_hand_ls, 1.e99);
		rbs.emplace_back(right_hand_ls, 1.e99);
	}

	void update(Simu &simu, Scalar time, Scalar /*dt*/) const override
	{
		const Scalar speed = speedup * M_PI_2 * m_config->units().toSI(Units::Time);

		const Scalar tw = time * speed;
		Vec vel = Vec::Zero();
		Vec avel = Vec::Zero();

		if (tw > 2)
		{
			const Scalar t = tw - 2;
			if (t < 1)
			{
				avel = Vec(-M_PI_2, 0, 0);
			}
		}
		else if (tw > 1)
		{
			const Scalar t = tw - 1;
			if (t < 1)
			{
				vel = Vec(0, 0, 6 * (t - t * t)) * .4 * m_config->box[2];
			}
		}
		else if (tw > 0)
		{
			const Scalar t = tw;

			const Eigen::AngleAxis<Scalar> firstRot(
				Eigen::AngleAxis<Scalar>(-M_PI / 6, Vec(0, 1, 0)) *
				Eigen::AngleAxis<Scalar>(.8 * M_PI_2, Vec(1, 0, 0)));

			if (t < 1)
			{
				avel = firstRot.angle() * firstRot.axis();
				vel = Vec(0, std::sin(t) * .75 * m_config->box[1], -std::cos(t) * .7 * m_config->box[2]);
			}
		}

		vel *= speed;
		avel *= speed;

		int i = 1;
		for (RigidBody &rb : simu.rigidBodies())
		{
			rb.set_velocity(Vec(vel[0], i * vel[1], vel[2]), Vec(i * avel[0], avel[1], i * avel[2]));
			i = -i;
		}
	}

	Scalar speedup;
	mutable LevelSet::Ptr left_hand_ls;
	mutable LevelSet::Ptr right_hand_ls;
};

struct SplittingScenar : public Scenario
{
	Vec center;
	Scalar radius;

	Scalar volMass;

	Scalar h0;
	Scalar hvel;
	Scalar zvel;
	Scalar avel;

	Scalar particle_density(const Vec &x) const override
	{
		return (std::fabs(x[0] - center[0]) < radius && x[2] > center[2] + radius)
				   ? 1
				   : 0;
	}

	virtual void init(const Params &) override
	{
		center = Vec(.5 * m_config->box[0], .5 * m_config->box[1], .5 * m_config->box[2]);
		radius = .125 * m_config->box[0];
	}

	void add_rigid_bodies(std::vector<RigidBody> &rbs) const override
	{
		LevelSet::Ptr ls = LevelSet::make_cylinder(m_config->box[1] / radius);

		Vec pos = center;
		ls->scale(radius).set_origin(pos).rotate(Vec(1, 0, 0), M_PI_2);

		rbs.emplace_back(ls, 1.e99);
	}
};

struct CylinderScenar : public Scenario
{
	Vec center;
	Scalar radius;

	Scalar particle_density(const Vec &x) const override
	{

		return (x.head<2>() - center.head<2>()).squaredNorm() < radius * radius ? 1 : 0;
	}

	virtual void init(const Params &params) override
	{
		center = Vec(.5 * m_config->box[0], .5 * m_config->box[1], .5 * m_config->box[2]);
		Scalar r = scalar_param(params, "r", Units::None, .125);
		radius = r * m_config->box[0];
	}
};

struct CollapseScenarLHEDoor : public Scenario
{
	Scalar particle_density(const Vec &x) const override
	{
		return ((x[0] < columnLength || x[2] < .1 * m_config->box[2]) && (x[2] < (fracH + 0.1) * m_config->box[2]) && (x[1]<(m_config->box[1]-w_sideW))&& (x[1]>w_sideW) && (x[0]>w_sideW) ) ? 1. : 0.;
	}

	virtual void init(const Params &params) override
	{
		tauD = scalar_param(params, "taudoor", Units::Time, .06);
		velD = scalar_param(params, "veldoor", Units::Velocity, 1.0);
		muDoor = scalar_param(params, "mud", Units::None, -99);
		ts = scalar_param(params, "ts", Units::Time, 0);
		fracH = scalar_param(params, "frac_h", Units::None, 1.);
		columnLength = scalar_param(params, "column_length", Units::Length, 1.);
		w_sideW = scalar_param(params, "wsw", Units::Length, 0.01);
	}

	void add_rigid_bodies(std::vector<RigidBody> &rbs) const override
	{
		// 		const Scalar a = 0.1 ;
		Vec box = m_config->box;
		const Scalar L = box[2] * (0.9);
		const Scalar W = 2 * m_config->typicalLength();
		Vec doorBox = Vec(W, 0.6 * box[1], 0.5 * L);

		LevelSet::Ptr ls = LevelSet::make_box(doorBox);
		ls->set_origin(Vec(columnLength + W, box[1] * 0.5, .1 * box[2] + L * 0.5));
		rbs.emplace_back(ls, 1.e99);
		rbs.back().set_mu(muDoor);
		d6::Log::Debug() << "mu_door = " << muDoor << std::endl ;
		std::ofstream RBout(m_config->base_dir + "/door.txt");
		d6::dump(RBout, "Nframe,time,X,Y,Z,Ux,Uy,Uz");
		RBout << "\n";
		RBout.close();
		// Scalar w_sideW=m_config->units().fromSI(Units::Length) * 0.01;
		
		if (w_sideW > 0.0) {
		Vec wallBox = Vec(box[0], w_sideW/2, box[2]);
		auto leftWall = LevelSet::make_box(wallBox);
		leftWall->set_origin(Vec(box[0] * 0.5, w_sideW/2 ,box[2] * 0.5));
		rbs.emplace_back(leftWall, 1.e99);
		auto rightWall = LevelSet::make_box(wallBox);
		rightWall->set_origin(Vec(box[0] * 0.5, box[1]- w_sideW/2  ,box[2] * 0.5));
		rbs.emplace_back(rightWall, 1.e99);
		}

		Scalar back_widthWall = 0.01*m_config->units().fromSI(Units::Length);
		Vec bWallBox = Vec(back_widthWall / 2, box[1], box[2]);
		auto backWall = LevelSet::make_box(bWallBox);
		backWall->set_origin(Vec(w_sideW/2 , 0.5*box[1],box[2] * 0.5));
		rbs.emplace_back(backWall, 1.e99);

		// auto rightWall = LevelSet::make_plane();
		// rightWall->rotate(Vec(1, 0, 0), M_PI_2).set_origin(Vec(0, box[1], 0));
		// rbs.emplace_back(rightWall, 1.e99);

		// auto sideWall = LevelSet::make_plane();
		// sideWall->rotate(Vec(0, 1, 0), M_PI_2).set_origin(Vec(0, 0, 0));
		// rbs.emplace_back(sideWall, 1.e99);

	}

	void update(Simu &simu, Scalar time, Scalar /*dt*/) const override
	{
		Scalar iF = time * m_config->fps;

		Scalar speedz = velD * (1 - exp(-(time - ts) / tauD));
		//         Scalar speedy = 0.1*m_config->units().fromSI( Units::Velocity)*(1-exp(-t/(0.05)));
		Vec vel = Vec::Zero();
		if (time > ts)
		{
			vel[2] = speedz;
		}

			simu.rigidBodies()[0].set_velocity(vel, VecR::Zero());
			simu.rigidBodies()[0];
			Vec position = simu.rigidBodies()[0].levelSet().origin();
			if (std::abs(iF - static_cast<int>(std::round(iF))) < 1.e-8)
			{
			std::ofstream RBout(m_config->base_dir + "/door.txt", std::fstream::app);
			RBout << std::round(time * m_config->fps)<<",";
			RBout << time * m_config->units().toSI(Units::Time)<<",";
			RBout << position[0] * m_config->units().toSI(Units::Length)<<",";
			RBout << position[1] * m_config->units().toSI(Units::Length)<<",";
			RBout << position[2] * m_config->units().toSI(Units::Length)<<",";
			RBout << vel[0] * m_config->units().toSI(Units::Velocity)<<",";
			RBout << vel[1] * m_config->units().toSI(Units::Velocity)<<",";
			RBout << vel[2] * m_config->units().toSI(Units::Velocity)<<"\n";
			RBout.close();
				RBout.close();
			}


		}
	

private:
	Scalar tauD;
	Scalar velD;
	Scalar ts;
	Scalar fracH;
	Scalar columnLength;
	Scalar w_sideW;
	Scalar muDoor;
};


struct CollapseScenarLHEDoorH : public Scenario
{
	Scalar particle_density(const Vec &x) const override
	{
		return ((x[0] < columnLength || x[2] < .1 * m_config->box[2]) && (x[2] < (fracH + 0.1) * m_config->box[2]) && (x[1]<(m_config->box[1]-w_sideW))&& (x[1]>w_sideW) && (x[0]>w_sideW) ) ? 1. : 0.;
	}

	virtual void init(const Params &params) override
	{
		tauD = scalar_param(params, "taudoor", Units::Time, .06);
		velD = scalar_param(params, "veldoor", Units::Velocity, 1.0);
		ts = scalar_param(params, "ts", Units::Time, 0);
		fracH = scalar_param(params, "frac_h", Units::None, 1.);
		columnLength = scalar_param(params, "column_length", Units::Length, 1.);
		w_sideW = scalar_param(params, "wsw", Units::Length, 0.01);
	}

	void add_rigid_bodies(std::vector<RigidBody> &rbs) const override
	{
		// 		const Scalar a = 0.1 ;
		Vec box = m_config->box;
		const Scalar L = box[2] * (0.9);
		const Scalar W = 2 * m_config->typicalLength();
		Vec doorBox = Vec(W, 0.6 * box[1], 0.5 * L);

		LevelSet::Ptr ls = LevelSet::make_box(doorBox);
		ls->set_origin(Vec(columnLength + W, box[1] * 0.5, .1 * box[2] + L * 0.5));
		rbs.emplace_back(ls, 1.e99);
		std::ofstream RBout(m_config->base_dir + "/door.txt");
		d6::dump(RBout, "Nframe,time,X,Y,Z,Ux,Uy,Uz");
		RBout << "\n";
		RBout.close();
		// Scalar w_sideW=m_config->units().fromSI(Units::Length) * 0.01;
		
		if (w_sideW > 0.0) {
		Vec wallBox = Vec(box[0], w_sideW/2, box[2]);
		auto leftWall = LevelSet::make_box(wallBox);
		leftWall->set_origin(Vec(box[0] * 0.5, w_sideW/2 ,box[2] * 0.5));
		rbs.emplace_back(leftWall, 1.e99);
		auto rightWall = LevelSet::make_box(wallBox);
		rightWall->set_origin(Vec(box[0] * 0.5, box[1]- w_sideW/2  ,box[2] * 0.5));
		rbs.emplace_back(rightWall, 1.e99);
		}

		Scalar back_widthWall = 0.01*m_config->units().fromSI(Units::Length);
		Vec bWallBox = Vec(back_widthWall / 2, box[1], box[2]);
		auto backWall = LevelSet::make_box(bWallBox);
		backWall->set_origin(Vec(w_sideW/2 , 0.5*box[1],box[2] * 0.5));
		rbs.emplace_back(backWall, 1.e99);

		// auto rightWall = LevelSet::make_plane();
		// rightWall->rotate(Vec(1, 0, 0), M_PI_2).set_origin(Vec(0, box[1], 0));
		// rbs.emplace_back(rightWall, 1.e99);

		// auto sideWall = LevelSet::make_plane();
		// sideWall->rotate(Vec(0, 1, 0), M_PI_2).set_origin(Vec(0, 0, 0));
		// rbs.emplace_back(sideWall, 1.e99);

	}

	void update(Simu &simu, Scalar time, Scalar /*dt*/) const override
	{
		Scalar iF = time * m_config->fps;

		Scalar speedz = velD * (1 - exp(-(time - ts) / tauD));
		//         Scalar speedy = 0.1*m_config->units().fromSI( Units::Velocity)*(1-exp(-t/(0.05)));
		Vec vel = Vec::Zero();
		if (time > ts)
		{
			vel[0] = speedz;
		}

			simu.rigidBodies()[0].set_velocity(vel, VecR::Zero());
			simu.rigidBodies()[0];
			Vec position = simu.rigidBodies()[0].levelSet().origin();
			if (std::abs(iF - static_cast<int>(std::round(iF))) < 1.e-8)
			{
			std::ofstream RBout(m_config->base_dir + "/door.txt", std::fstream::app);
			RBout << std::round(time * m_config->fps)<<",";
			RBout << time * m_config->units().toSI(Units::Time)<<",";
			RBout << position[0] * m_config->units().toSI(Units::Length)<<",";
			RBout << position[1] * m_config->units().toSI(Units::Length)<<",";
			RBout << position[2] * m_config->units().toSI(Units::Length)<<",";
			RBout << vel[0] * m_config->units().toSI(Units::Velocity)<<",";
			RBout << vel[1] * m_config->units().toSI(Units::Velocity)<<",";
			RBout << vel[2] * m_config->units().toSI(Units::Velocity)<<"\n";
			RBout.close();
				RBout.close();
			}


		}
	

private:
	Scalar tauD;
	Scalar velD;
	Scalar ts;
	Scalar fracH;
	Scalar columnLength;
	Scalar w_sideW;
};


struct CollapseIonescu : public Scenario
{
	Scalar particle_density(const Vec &x) const override
	{	Scalar w_sideW=m_config->units().fromSI(Units::Length) * 0.01;
		return ((x[0] < columnLength) && (x[2] > (0.) * m_config->box[2]) && (x[2] < (fracH) * m_config->box[2]) && (x[1]<(m_config->box[1]-w_sideW))&& (x[1]>w_sideW) && (x[0]>w_sideW) ) ? 1. : 0.;
	}

	virtual void init(const Params &params) override
	{
		tauD = scalar_param(params, "taudoor", Units::Time, .06);
		velD = scalar_param(params, "veldoor", Units::Velocity, 1.0);
		ts = scalar_param(params, "ts", Units::Time, 0.5);
		fracH = scalar_param(params, "frac_h", Units::None, 1.);
		columnLength = scalar_param(params, "column_length", Units::Length, 1.);
	}

	void add_rigid_bodies(std::vector<RigidBody> &rbs) const override
	{
		// 		const Scalar a = 0.1 ;
		Vec box = m_config->box;
		const Scalar L = box[2] * (0.9);
		const Scalar W = 2 * m_config->typicalLength();
		Vec doorBox = Vec(W, 0.6 * box[1], 0.5 * L);

		LevelSet::Ptr ls = LevelSet::make_box(doorBox);
		ls->set_origin(Vec(columnLength + W, box[1] * 0.5, L * 0.5));
		rbs.emplace_back(ls, 1.e99);
	
		std::ofstream RBout(m_config->base_dir + "/door.txt");
		d6::dump(RBout, "Nframe,time,X,Y,Z,Ux,Uy,Uz");
		RBout << "\n";
		RBout.close();
		Scalar w_sideW=m_config->units().fromSI(Units::Length) * 0.01;
		Vec wallBox = Vec(box[0], w_sideW/2, box[2]);

		auto leftWall = LevelSet::make_box(wallBox);
		leftWall->set_origin(Vec(box[0] * 0.5, w_sideW/2 ,box[2] * 0.5));
		rbs.emplace_back(leftWall, 1.e99);
		auto rightWall = LevelSet::make_box(wallBox);
		rightWall->set_origin(Vec(box[0] * 0.5, box[1]- w_sideW/2  ,box[2] * 0.5));
		rbs.emplace_back(rightWall, 1.e99);

		Vec bWallBox = Vec( w_sideW/2, box[1],box[2]);
		auto backWall = LevelSet::make_box(bWallBox);
		backWall->set_origin(Vec(w_sideW/2 , 0.5*box[1],box[2] * 0.5));
		rbs.emplace_back(backWall, 1.e99);

		// auto bottom = LevelSet::make_box(Vec( box[0], box[1],0.05*box[2]));
		// bottom->set_origin(Vec(0.5*box[0] , 0.5*box[1],box[2] * 0.0));
		// rbs.emplace_back(bottom, 1.e99);
		// rbs.back().set_mu(1);
		// auto rightWall = LevelSet::make_plane();
		// rightWall->rotate(Vec(1, 0, 0), M_PI_2).set_origin(Vec(0, box[1], 0));
		// rbs.emplace_back(rightWall, 1.e99);

		// auto sideWall = LevelSet::make_plane();
		// sideWall->rotate(Vec(0, 1, 0), M_PI_2).set_origin(Vec(0, 0, 0));
		// rbs.emplace_back(sideWall, 1.e99);

	}

	void update(Simu &simu, Scalar time, Scalar /*dt*/) const override
	{
		Scalar t = time * m_config->fps;
		Scalar tau = tauD * m_config->fps;

		Scalar speedz = m_config->units().fromSI( Units::Velocity)*(2.3);
		//         Scalar speedy = 0.1*m_config->units().fromSI( Units::Velocity)*(1-exp(-t/(0.05)));
		Vec vel = Vec::Zero();
		if (time > ts)
		{
			vel[2] = speedz;
		}
		int n =0 ;
		for (RigidBody &rb : simu.rigidBodies())
		{
			if(n == 0){
			rb.set_velocity(vel, VecR::Zero());

			std::cout << time * m_config->fps << std::endl;
			std::cout << std::abs(t - static_cast<int>(std::round(t))) << std::endl;

			if (std::abs(t - static_cast<int>(std::round(t))) < 1.e-8)
			{
				std::ofstream RBout(m_config->base_dir + "/door.txt", std::fstream::app);
				RBout << std::round(time * m_config->fps);
				RBout << ",";
				RBout << time * m_config->units().toSI(Units::Time);
				RBout << ",";
				RBout << rb.levelSet().origin()[0] * m_config->units().toSI(Units::Length);
				RBout << ",";
				RBout << rb.levelSet().origin()[1] * m_config->units().toSI(Units::Length);
				RBout << ",";
				RBout << rb.levelSet().origin()[2] * m_config->units().toSI(Units::Length);
				RBout << ",";
				RBout << vel[0] * m_config->units().toSI(Units::Velocity);
				RBout << ",";
				RBout << vel[1] * m_config->units().toSI(Units::Velocity);
				RBout << ",";
				RBout << vel[2] * m_config->units().toSI(Units::Velocity);
				RBout << "\n";
				RBout.close();
			}
			n = n + 1;
			}

			break;
		}
	}

private:
	Scalar tauD;
	Scalar velD;
	Scalar ts;
	Scalar fracH;
	Scalar columnLength;

};


struct InclinedPlane : public Scenario
{
	Scalar particle_density(const Vec &x) const override
	{	

		return ((x[0] < columnLength || x[2] < zD*0.0)  && (x[2] < (fracH) * m_config->box[2]) ) ? 1. : 0.;
	}

	virtual void init(const Params &params) override
	{
		tauD = scalar_param(params, "taudoor", Units::Time, .06);
		velD = scalar_param(params, "veldoor", Units::Velocity, 0.0);
		zD = scalar_param(params, "zdoor", Units::Length, 0);
		ts = scalar_param(params, "ts", Units::Time, 0.5);
		fracH = scalar_param(params, "frac_h", Units::None, 1.);
		columnLength = scalar_param(params, "column_length", Units::Length, 1.);

	}

	void add_rigid_bodies(std::vector<RigidBody> &rbs) const override
	{
		// 		const Scalar a = 0.1 ;
		Vec box = m_config->box;
		const Scalar L = box[2] * (0.9);
		const Scalar W = 2 * m_config->typicalLength();
		Vec doorBox = Vec(W, 0.6 * box[1], 0.5 * L);

		LevelSet::Ptr ls = LevelSet::make_box(doorBox);
		ls->set_origin(Vec(columnLength + W, box[1] * 0.5, L * 0.5+zD));
		rbs.emplace_back(ls, 1.e99);
		
		std::ofstream RBout(m_config->base_dir + "/door.txt");
		d6::dump(RBout, "Nframe,time,X,Y,Z,Ux,Uy,Uz");
		RBout << "\n";
		RBout.close();



		// auto bottom = LevelSet::make_box(Vec( box[0], box[1],0.05*box[2]));
		// bottom->set_origin(Vec(0.5*box[0] , 0.5*box[1],box[2] * 0.0));
		// rbs.emplace_back(bottom, 1.e99);
		// rbs.back().set_mu(1);
		// auto rightWall = LevelSet::make_plane();
		// rightWall->rotate(Vec(1, 0, 0), M_PI_2).set_origin(Vec(0, box[1], 0));
		// rbs.emplace_back(rightWall, 1.e99);

		// auto sideWall = LevelSet::make_plane();
		// sideWall->rotate(Vec(0, 1, 0), M_PI_2).set_origin(Vec(0, 0, 0));
		// rbs.emplace_back(sideWall, 1.e99);

	}

	void update(Simu &simu, Scalar time, Scalar /*dt*/) const override
	{
		Scalar t = time * m_config->fps;
		Scalar tau = tauD * m_config->fps;

		Scalar speedz = m_config->units().fromSI( Units::Velocity)*(0.0);
		//         Scalar speedy = 0.1*m_config->units().fromSI( Units::Velocity)*(1-exp(-t/(0.05)));
		Vec vel = Vec::Zero();
		if (time > ts)
		{
			vel[2] = speedz;
		}
		int n =0 ;
		for (RigidBody &rb : simu.rigidBodies())
		{
			if(n == 0){
			rb.set_velocity(vel, VecR::Zero());

			std::cout << time * m_config->fps << std::endl;
			std::cout << std::abs(t - static_cast<int>(std::round(t))) << std::endl;

			if (std::abs(t - static_cast<int>(std::round(t))) < 1.e-8)
			{
				std::ofstream RBout(m_config->base_dir + "/door.txt", std::fstream::app);
				RBout << std::round(time * m_config->fps);
				RBout << ",";
				RBout << time * m_config->units().toSI(Units::Time);
				RBout << ",";
				RBout << rb.levelSet().origin()[0] * m_config->units().toSI(Units::Length);
				RBout << ",";
				RBout << rb.levelSet().origin()[1] * m_config->units().toSI(Units::Length);
				RBout << ",";
				RBout << rb.levelSet().origin()[2] * m_config->units().toSI(Units::Length);
				RBout << ",";
				RBout << vel[0] * m_config->units().toSI(Units::Velocity);
				RBout << ",";
				RBout << vel[1] * m_config->units().toSI(Units::Velocity);
				RBout << ",";
				RBout << vel[2] * m_config->units().toSI(Units::Velocity);
				RBout << "\n";
				RBout.close();
			}
			n = n + 1;
			}

			break;
		}
	}

private:
	Scalar tauD;
	Scalar velD;
	Scalar ts;
	Scalar fracH;
	Scalar columnLength;
	Scalar zD;
};

struct CollapseScenarDoor : public Scenario
{
	Scalar particle_density(const Vec &x) const override
	{
		return ((x[0] < columnLength || x[2] < .1 * m_config->box[2]) && (x[2] < (fracH + 0.1) * m_config->box[2])) ? 1. : 0.;
	}

	virtual void init(const Params &params) override
	{
		tauD = scalar_param(params, "taudoor", Units::Time, .06);
		velD = scalar_param(params, "veldoor", Units::Velocity, 1.0);
		ts = scalar_param(params, "ts", Units::None, 15);
		fracH = scalar_param(params, "frac_h", Units::None, 1.);
		columnLength = scalar_param(params, "column_length", Units::Length, 1.);
	}

	void add_rigid_bodies(std::vector<RigidBody> &rbs) const override
	{
		// 		const Scalar a = 0.1 ;
		Vec box = m_config->box;
		const Scalar L = box[2] * (0.9);
		const Scalar W = 2 * m_config->typicalLength();
		Vec doorBox = Vec(W, 0.6 * box[1], 0.5 * L);

		LevelSet::Ptr ls = LevelSet::make_box(doorBox);
		ls->set_origin(Vec(columnLength + W, box[1] * 0.5, .1 * box[2] + L * 0.5));
		rbs.emplace_back(ls, 1.e99);
		std::ofstream RBout(m_config->base_dir + "/door.txt");
		d6::dump(RBout, "Nframe,time,X,Y,Ux,Uy");
		RBout << "\n";
		RBout.close();
	}

	void update(Simu &simu, Scalar time, Scalar /*dt*/) const override
	{
		Scalar t = time * m_config->fps;
		Scalar tau = tauD * m_config->fps;

		Scalar speedy = velD * (1 - exp(-(t - ts) / tau));
		//         Scalar speedy = 0.1*m_config->units().fromSI( Units::Velocity)*(1-exp(-t/(0.05)));
		Vec vel = Vec::Zero();
		if (t > ts)
		{
			vel[2] = speedy;
		}

		for (RigidBody &rb : simu.rigidBodies())
		{
			rb.set_velocity(vel, VecR::Zero());
			//             std::cout << vel << std::endl;
			// 			std::cout <<  rb.levelSet().origin() << std::endl;

			std::cout << time * m_config->fps << std::endl;
			std::cout << std::abs(t - static_cast<int>(std::round(t))) << std::endl;

			if (std::abs(t - static_cast<int>(std::round(t))) < 1.e-8)
			{

				std::ofstream RBout(m_config->base_dir + "/door.txt", std::fstream::app);
				RBout << std::round(time * m_config->fps);
				RBout << ",";
				RBout << time * m_config->units().toSI(Units::Time);
				RBout << ",";
				RBout << rb.levelSet().origin()[0] * m_config->units().toSI(Units::Length);
				RBout << ",";
				RBout << rb.levelSet().origin()[1] * m_config->units().toSI(Units::Length);
				RBout << ",";
				RBout << vel[0] * m_config->units().toSI(Units::Velocity);
				RBout << ",";
				RBout << vel[1] * m_config->units().toSI(Units::Velocity);
				RBout << "\n";
				RBout.close();
			}

			break;
		}
	}

private:
	Scalar tauD;
	Scalar velD;
	Scalar ts;
	Scalar fracH;
	Scalar columnLength;
};

// Factories & stuff

std::unique_ptr<Scenario> DefaultScenarioFactory::make(const std::string &str) const
{
	if (str == "rayleigh")
		return std::unique_ptr<Scenario>(new RayleighScenar());
	if (str == "collapselhedoor")
		return std::unique_ptr<Scenario>(new CollapseScenarLHEDoor());
	if (str == "collapselhedoorh")
		return std::unique_ptr<Scenario>(new CollapseScenarLHEDoorH());
	if (str == "collapseionescu")
		return std::unique_ptr<Scenario>(new CollapseIonescu());
	if (str == "inclinedplane")
		return std::unique_ptr<Scenario>(new InclinedPlane());
	if (str == "collapsedoor")
		return std::unique_ptr<Scenario>(new CollapseScenarDoor());
	if (str == "collapse")
		return std::unique_ptr<Scenario>(new CollapseScenar());
	if (str == "bridson")
		return std::unique_ptr<Scenario>(new BridsonScenar());
	if (str == "tower")
		return std::unique_ptr<Scenario>(new TowerScenar());
	if (str == "rb_plane_test")
		return std::unique_ptr<Scenario>(new RbPlaneTestScenar());
	if (str == "impact")
		return std::unique_ptr<Scenario>(new ImpactScenar());
	if (str == "impactlhe")
		return std::unique_ptr<Scenario>(new ImpactScenarlhe());
	if (str == "silo")
		return std::unique_ptr<Scenario>(new SiloScenar());
	if (str == "bunny")
		return std::unique_ptr<Scenario>(new BunnyScenar());
	if (str == "writing")
		return std::unique_ptr<Scenario>(new WritingScenar());
	if (str == "hourglass")
		return std::unique_ptr<Scenario>(new HourglassScenar());
	if (str == "wheel")
		return std::unique_ptr<Scenario>(new WheelScenar());
	if (str == "digging")
		return std::unique_ptr<Scenario>(new DiggingScenar());
	if (str == "splitting")
		return std::unique_ptr<Scenario>(new SplittingScenar());
	if (str == "cylinder")
		return std::unique_ptr<Scenario>(new CylinderScenar());

	return std::unique_ptr<Scenario>(new BedScenar());
}

} // namespace d6
