#include "RigidBody.hh"

#include "LevelSet.hh"

#include <iostream>

namespace d6 {

RigidBody::RigidBody(std::unique_ptr<LevelSet> &ls, Scalar volMass )
		: m_levelSet ( std::move(ls) ), m_volumicMass( volMass )
{
	m_velocity.setZero() ;
}

Vec RigidBody::velocity_at(const Vec &x) const
{
	return velocity() + angularVelocity()(0,0) * Vec( -x[1], x[0] ) ;
}

void RigidBody::integrate_forces(const Scalar dt, const VecS &forces)
{
	MatS Mi ;
	inv_inertia( Mi );

	m_velocity += dt * Mi * forces ;
}

void RigidBody::integrate_gravity(const Scalar dt, const Vec &gravity)
{
	VecS forces ;
	forces.head<WD>() = m_volumicMass * m_levelSet->volume() * gravity ;
	forces.tail<RD>().setZero() ;

	integrate_forces( dt, forces );
}

void RigidBody::move(const Scalar dt) const
{
	m_levelSet->move( dt * velocity(), dt * angularVelocity()(0,0) );
}

void RigidBody::move_to(const Vec &pos) const
{
	m_levelSet->set_origin( pos ) ;
}

void RigidBody::inv_inertia( MatS& Mi ) const
{
	m_levelSet->inv_inertia( Mi ) ;
	Mi /= m_volumicMass ;
}

} //d6
