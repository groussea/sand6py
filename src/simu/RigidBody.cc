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
	return velocity() + angularVelocity().cross( x ) ;
}

void RigidBody::integrate_forces(const Scalar dt, const Vec6 &forces)
{
	Mat66 Mi ;
	inv_inertia( Mi );

	m_velocity += dt * Mi * forces ;
}

void RigidBody::integrate_gravity(const Scalar dt, const Vec &gravity)
{
	Vec6 forces ;
	forces.head<3>() = m_volumicMass * m_levelSet->volume() * gravity ;
	forces.tail<3>().setZero() ;

	integrate_forces( dt, forces );
}

void RigidBody::move(const Scalar dt) const
{

	const Scalar avn = angularVelocity().norm() ;
	Vec axis = avn > 1.e-12
			? Vec(angularVelocity()/avn)
			: Vec(1,0,0) ;
	Eigen::AngleAxis< Scalar > aa( dt * avn, axis ) ;

	m_levelSet->move( dt * velocity(), Quaternion( aa ) );
}

void RigidBody::inv_inertia( Mat66& Mi ) const
{
	m_levelSet->inv_inertia( Mi ) ;
	Mi /= m_volumicMass ;
}

} //d6
