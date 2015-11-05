#ifndef D6_RIGID_BODY_HH
#define D6_RIGID_BODY_HH

#include "utils/alg.hh"
#include <memory>

namespace d6 {

class LevelSet ;

class RigidBody
{
public:
	explicit RigidBody( std::unique_ptr< LevelSet > &ls ) ;

	Vec velocity_at( const Vec& x ) const ;

	const LevelSet& levelSet() const
	{ return *m_levelSet ; }

	void set_velocity( const Vec& velocity, const Vec& angularVelocity )
	{
		m_velocity.head<3>() = velocity ;
		m_velocity.tail<3>() = angularVelocity ;
	}

	void predict_velocity( const Scalar dt, const Vec6& forces ) const ;

	void move( const Scalar dt ) const ;

private:
	std::unique_ptr< LevelSet > m_levelSet ;

	Vec6 m_velocity ;

} ;

} //d6


#endif
