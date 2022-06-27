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
		return (x[1] > .5 * box()[1] &&
				(x - .5 * box()).squaredNorm() > std::pow(box()[0] / 64, 2))
				   ? 1.
				   : 0.;
	}
};

struct Sedimentation : public Scenario
{
	Scalar particle_density(const Vec &x) const override
	{
		return (x[1] < h0 * box()[1]) ? phi0 : 0;
	}
	virtual void init(const Params &params) override
	{
		h0 = scalar_param(params, "h0", Units::None, .9);
		phi0 = scalar_param(params, "phi0", Units::None, .3);
	}

private:
	Scalar h0;
	Scalar phi0;
};

struct FallingBallScenar : public Scenario
{
	Scalar particle_density(const Vec &x) const override
	{
		return (x - (.5 * box() + Vec(0, .375 * box()[1]))).squaredNorm() < std::pow(box()[0] / 4, 2) ? 1. : 0.;
	}
};

struct BedScenar : public Scenario
{
	Scalar particle_density(const Vec &x) const override
	{
		return (x[1] < .5 * m_config->box[1]) ? 1. : 0.;
	}
};

struct CollapseScenarLHEDoor : public Scenario
{
	Scalar particle_density(const Vec &x) const override
	{
		return ((x[0] < columnLength || x[1] < Hbed * m_config->box[1]) & (x[1] < (fracH + Hbed) * m_config->box[1])) ? 1. : 0.;
	}

	virtual void init(const Params &params) override
	{
		tauD = scalar_param(params, "taudoor", Units::Time, .06);
		velD = scalar_param(params, "veldoor", Units::Velocity, 1.0);
		ts = scalar_param(params, "ts", Units::Time, 15);
		fracH = scalar_param(params, "frac_h", Units::None, 1.);
		columnLength= scalar_param(params, "column_length", Units::Length, 0.);
		zD= scalar_param(params, "zdoor", Units::Length, 0.);
		Hbed= scalar_param(params, "hbed", Units::Length, 0.1);
	}

	void add_rigid_bodies(std::vector<RigidBody> &rbs) const override
	{
		// 		const Scalar a = 0.1 ;
		const Scalar L = m_config->box[1] * (1-Hbed);

		LevelSet::Ptr ls = LevelSet::make_cylinder(L);
		ls->set_origin(Vec(columnLength + m_config->typicalLength(), Hbed* m_config->box[1] + L / 2 + m_config->typicalLength()+zD));
		ls->set_rotation(0);
		rbs.emplace_back(ls, 1.);
		rbs.back().set_velocity(Vec(0, 0), 0.);
		std::ofstream RBout(m_config->base_dir + "/door.txt");
		d6::dump(RBout, "Nframe,time,X,Y,Ux,Uy");
		RBout << "\n";
		RBout.close();
	}

	void update(Simu &simu, Scalar time, Scalar /*dt*/) const override
	{
		Scalar iF = time * m_config->fps;
		Scalar speedy = velD * (1 - exp(-(time - ts) / tauD));
		//         Scalar speedy = 0.1*m_config->units().fromSI( Units::Velocity)*(1-exp(-t/(0.05)));
		Vec vel = Vec::Zero();
		if (time > ts)
		{
			vel[1] = speedy;
		}
			if (std::abs(iF - static_cast<int>(std::round(iF))) < 1.e-8)
			{
			simu.rigidBodies()[0].set_velocity(vel, 0);
			Vec position = simu.rigidBodies()[0].levelSet().origin();
			std::ofstream RBout(m_config->base_dir + "/door.txt", std::fstream::app);
			RBout << std::round(time * m_config->fps)<<",";
			RBout << time * m_config->units().toSI(Units::Time)<<",";
			RBout << position[0] * m_config->units().toSI(Units::Length)<<",";
			RBout << position[1] * m_config->units().toSI(Units::Length)<<",";
			RBout << vel[0] * m_config->units().toSI(Units::Velocity)<<",";
			RBout << vel[1] * m_config->units().toSI(Units::Velocity)<<"\n";
			RBout.close();
		}
	}

private:
	Scalar tauD;
	Scalar velD;
	Scalar ts;
	Scalar fracH;
	Scalar columnLength;
	Scalar zD;
	Scalar Hbed;
};

struct CollapseScenarLHEDoorH : public Scenario
{
	Scalar particle_density(const Vec &x) const override
	{
		return ((x[0] < columnLength || x[1] < .1 * m_config->box[1]) & (x[1] < (fracH + 0.1) * m_config->box[1])) ? 1. : 0.;
	}

	virtual void init(const Params &params) override
	{
		tauD = scalar_param(params, "taudoor", Units::Time, .06);
		velD = scalar_param(params, "veldoor", Units::Velocity, 1.0);
		ts = scalar_param(params, "ts", Units::Time, 15);
		fracH = scalar_param(params, "frac_h", Units::None, 1.);
		columnLength= scalar_param(params, "column_length", Units::Length, 0.);
		zD= scalar_param(params, "zdoor", Units::Length, 0.);
	}

	void add_rigid_bodies(std::vector<RigidBody> &rbs) const override
	{
		// 		const Scalar a = 0.1 ;
		const Scalar L = m_config->box[1] * (0.9);

		LevelSet::Ptr ls = LevelSet::make_cylinder(L);
		ls->set_origin(Vec(columnLength + m_config->typicalLength(), .1 * m_config->box[1] + L / 2 + m_config->typicalLength()+zD));
		ls->set_rotation(0);
		rbs.emplace_back(ls, 1.);
		rbs.back().set_velocity(Vec(0, 0), 0.);
		std::ofstream RBout(m_config->base_dir + "/door.txt");
		d6::dump(RBout, "Nframe,time,X,Y,Ux,Uy");
		RBout << "\n";
		RBout.close();
	}

	void update(Simu &simu, Scalar time, Scalar /*dt*/) const override
	{
		Scalar iF = time * m_config->fps;
		Scalar speedx = velD * (1 - exp(-(time - ts) / tauD));
		//         Scalar speedy = 0.1*m_config->units().fromSI( Units::Velocity)*(1-exp(-t/(0.05)));
		Vec vel = Vec::Zero();
		if (time > ts)
		{
			vel[0] = -speedx;
		}
			if (std::abs(iF - static_cast<int>(std::round(iF))) < 1.e-8)
			{
			simu.rigidBodies()[0].set_velocity(vel, 0);
			Vec position = simu.rigidBodies()[0].levelSet().origin();
			std::ofstream RBout(m_config->base_dir + "/door.txt", std::fstream::app);
			RBout << std::round(time * m_config->fps)<<",";
			RBout << time * m_config->units().toSI(Units::Time)<<",";
			RBout << position[0] * m_config->units().toSI(Units::Length)<<",";
			RBout << position[1] * m_config->units().toSI(Units::Length)<<",";
			RBout << vel[0] * m_config->units().toSI(Units::Velocity)<<",";
			RBout << vel[1] * m_config->units().toSI(Units::Velocity)<<"\n";
			RBout.close();
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


struct CollapseScenarLHE : public Scenario
{
	Scalar particle_density(const Vec &x) const override
	{
		return (x[0] < columnLength || x[1] < .1 * m_config->box[1]) ? 1. : 0.;
	}

	virtual void init(const Params &params) override
	{
		columnLength = scalar_param(params, "column_length", Units::Length, .25);
		
	}

private:
	Scalar l0;
	Scalar h0;
	Scalar columnLength;
};

struct CollapseScenarCohesion : public Scenario
{
	Scalar particle_density(const Vec &x) const override
	{
		return (x[0] < columnLength) ? 1. : 0.;
	}

	virtual void init(const Params &params) override
	{
		columnLength = scalar_param(params, "column_length", Units::Length, .25);
	}

private:
	Scalar l0;
	Scalar h0;
	Scalar columnLength;
};

struct CollapseScenar : public Scenario
{
	Scalar particle_density(const Vec &x) const override
	{
		return (x[0] < l0 * box()[0] && x[1] < h0 * box()[1]) ? 1. : 0.;
	}

	virtual void init(const Params &params) override
	{
		l0 = scalar_param(params, "l0", Units::None, .25);
		h0 = scalar_param(params, "h0", Units::None, .75);
	}

private:
	Scalar l0;
	Scalar h0;
};

struct HeapScenar : public Scenario
{
	Scalar particle_density(const Vec &x) const override
	{
		return (x[0] > (.375 - .175 / 2) * box()[0] && x[0] < (.375 + .175 / 2) * box()[0] && x[1] < .75 * box()[1]) ? 1. : 0.;
	}
};

struct PlaneTestScenar : public Scenario
{
	Scalar particle_density(const Vec &x) const override
	{
		return (x[1] < .2 * box(1)) ? 1. : 0.;
	}

	void add_rigid_bodies(std::vector<RigidBody> &rbs) const override
	{
		LevelSet::Ptr ls = LevelSet::make_plane();
		ls->set_origin(.5 * m_config->box - Vec(0, .5 * m_config->box[1]));
		ls->set_rotation(M_PI / 2);

		rbs.emplace_back(ls, 1.);
		rbs.back().set_velocity(Vec(0, 0), 0);
	}
};

struct ImpactScenar : public Scenario
{

	Scalar particle_density(const Vec &x) const override
	{
		return (x[1] < 1. / 3. * m_config->box[1]) ? 1. : 0.;
	}

	virtual void init(const Params &params)
	{
		volMass = scalar_param(params, "vm", Units::VolumicMass, 1.5 * m_config->units().R);
		zvel = scalar_param(params, "zvel", Units::Velocity, 0.);
		avel = scalar_param(params, "avel", Units::Frequency, 0.);
		d = scalar_param(params, "d", Units::None, 0.25);
	}

	void add_rigid_bodies(std::vector<RigidBody> &rbs) const override
	{
		LevelSet::Ptr ls = LevelSet::make_sphere();

		ls->scale(radius()).set_origin(.5 * m_config->box + Vec(0, .25 * m_config->box[1]));

		rbs.emplace_back(ls, volMass);
		rbs.back().set_velocity(Vec(0, -zvel), avel);
	}

	Scalar radius() const
	{
		return d / 2 * m_config->box[0];
	}

	void update(Simu &simu, Scalar /*time*/, Scalar dt) const override
	{
		for (RigidBody &rb : simu.rigidBodies())
		{
			rb.integrate_gravity(dt, m_config->gravity);
		}
	}

private:
	Scalar volMass;
	Scalar zvel;
	Scalar avel;
	Scalar d;
};

struct SiloScenar : public Scenario
{

	virtual void init(const Params &params)
	{
		volMass = scalar_param(params, "vm", Units::VolumicMass, 0.1 * m_config->units().R);
		hvel = scalar_param(params, "hvel", Units::Velocity, 1);
		avel = scalar_param(params, "avel", Units::Frequency, 0.);
		d = scalar_param(params, "d", Units::None, 0.2);
		vpos = scalar_param(params, "vpos", Units::None, 0.5);
	}
	Scalar particle_density(const Vec &x) const override
	{
		return (x[1] - 1 > vpos * m_config->box[1]) ? 1. : 0.;
	}
	void add_rigid_bodies(std::vector<RigidBody> &rbs) const override
	{

		const Scalar L = m_config->box[0] * (1 - d) / 2;

		LevelSet::Ptr ls = LevelSet::make_cylinder(L);
		ls->set_origin(Vec(L / 2, vpos * m_config->box[1]));
		ls->set_rotation(M_PI / 2);

		LevelSet::Ptr ls2 = LevelSet::make_cylinder(L);
		ls2->set_origin(Vec(m_config->box[0] - L / 2, vpos * m_config->box[1]));
		ls2->set_rotation(M_PI / 2);

		rbs.emplace_back(ls2, 1.);
		rbs.emplace_back(ls, 1.);
	}

private:
	Scalar volMass;
	Scalar hvel;
	Scalar avel;
	Scalar d;
	Scalar vpos;
};

struct SiloScenarObstacle : public Scenario
{
	Scalar particle_density(const Vec &x) const override
	{
		return (x[1] - 1 > vpos * m_config->box[1] and
				std::pow(std::pow(x[0] - vpos * box(0), 2.) + std::pow(x[1] - vpos * box(1) - 3. * radius(), 2.), 0.5) > radius())
				   ? 1.
				   : 0.;
	}
	virtual void init(const Params &params)
	{
		volMass = scalar_param(params, "vm", Units::VolumicMass, 0.1 * m_config->units().R);
		hvel = scalar_param(params, "hvel", Units::Velocity, 1);
		avel = scalar_param(params, "avel", Units::Frequency, 0.);
		d = scalar_param(params, "d", Units::None, 0.2);
		vpos = scalar_param(params, "vpos", Units::None, 0.5);
	}
	void add_rigid_bodies(std::vector<RigidBody> &rbs) const override
	{

		const Scalar L = m_config->box[0] * (1 - d) / 2;

		LevelSet::Ptr ls = LevelSet::make_cylinder(L);
		ls->set_origin(Vec(L / 2, vpos * m_config->box[1]));
		ls->set_rotation(M_PI / 2);

		LevelSet::Ptr ls2 = LevelSet::make_cylinder(L);
		ls2->set_origin(Vec(m_config->box[0] - L / 2, vpos * m_config->box[1]));
		ls2->set_rotation(M_PI / 2);

		rbs.emplace_back(ls2, 1.);
		rbs.emplace_back(ls, 1.);

		LevelSet::Ptr ls3 = LevelSet::make_sphere();
		ls3->scale(radius()).set_origin(Vec(.5 * box(0), vpos * box(1) + 3 * radius()));

		// 		const Scalar t = .25*box(0) / hvel ;
		// 		const Scalar zvel = .5 * m_config->gravity.norm() * t ;

		rbs.emplace_back(ls3, std::numeric_limits<Scalar>::infinity());
		rbs.back().set_velocity(Vec(0, 0), 0);
	}
	Scalar radius() const
	{
		return d / 2 * m_config->box[0];
	}

private:
	Scalar volMass;
	Scalar hvel;
	Scalar avel;
	Scalar d;
	Scalar vpos;
};

struct TowerScenar : public Scenario
{

	Scalar particle_density(const Vec &x) const override
	{
		return (x[0] >= .25 * box(0) && x[0] <= .375 * box(0)) ? 1. : 0.;
	}

	virtual void init(const Params &params)
	{
		volMass = scalar_param(params, "vm", Units::VolumicMass, 0.1 * m_config->units().R);
		hvel = scalar_param(params, "hvel", Units::Velocity, 1);
		avel = scalar_param(params, "avel", Units::Frequency, 0.);
		d = scalar_param(params, "d", Units::None, 0.25);
	}

	void add_rigid_bodies(std::vector<RigidBody> &rbs) const override
	{
		LevelSet::Ptr ls = LevelSet::make_sphere();
		ls->scale(radius()).set_origin(Vec(0, .375 * box(1)));

		const Scalar t = .25 * box(0) / hvel;
		const Scalar zvel = .5 * m_config->gravity.norm() * t;

		rbs.emplace_back(ls, volMass);
		rbs.back().set_velocity(Vec(hvel, zvel), avel);
	}

	Scalar radius() const
	{
		return d / 2 * m_config->box[1];
	}

	void update(Simu &simu, Scalar /*time*/, Scalar dt) const override
	{
		for (RigidBody &rb : simu.rigidBodies())
		{
			rb.integrate_gravity(dt, m_config->gravity);
		}
	}

private:
	Scalar volMass;
	Scalar hvel;
	Scalar avel;
	Scalar d;
};

// Factories & stuff

std::unique_ptr<Scenario> DefaultScenarioFactory::make(const std::string &str) const
{
	if (str == "rayleigh")
		return std::unique_ptr<Scenario>(new RayleighScenar());
	if (str == "ball")
		return std::unique_ptr<Scenario>(new FallingBallScenar());
	if (str == "collapselhe")
		return std::unique_ptr<Scenario>(new CollapseScenarLHE());
	if (str == "collapselhedoor")
		return std::unique_ptr<Scenario>(new CollapseScenarLHEDoor());
	if (str == "collapselhedoorh")
		return std::unique_ptr<Scenario>(new CollapseScenarLHEDoorH());
	if (str == "collapsecohesion")
		return std::unique_ptr<Scenario>(new CollapseScenarCohesion());
	if (str == "collapse")
		return std::unique_ptr<Scenario>(new CollapseScenar());
	if (str == "planetest")
		return std::unique_ptr<Scenario>(new PlaneTestScenar());
	if (str == "impact")
		return std::unique_ptr<Scenario>(new ImpactScenar());
	if (str == "silo")
		return std::unique_ptr<Scenario>(new SiloScenar());
	if (str == "silobstacle")
		return std::unique_ptr<Scenario>(new SiloScenarObstacle());
	if (str == "tower")
		return std::unique_ptr<Scenario>(new TowerScenar());
	if (str == "sedim")
		return std::unique_ptr<Scenario>(new Sedimentation());
	if (str == "heap")
		return std::unique_ptr<Scenario>(new HeapScenar());
	return std::unique_ptr<Scenario>(new BedScenar());
}

} // namespace d6
