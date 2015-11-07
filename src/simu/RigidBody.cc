#include "RigidBody.hh"

#include "LevelSet.hh"

namespace d6 {

RigidBody::RigidBody(std::unique_ptr<LevelSet> &ls )
		: m_levelSet ( std::move(ls) )
{
	m_velocity.setZero() ;
}

Vec RigidBody::velocity_at(const Vec &x) const
{
	return velocity() + angularVelocity().cross( x ) ;
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


} //d6
