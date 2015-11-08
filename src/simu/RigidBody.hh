#ifndef D6_RIGID_BODY_HH
#define D6_RIGID_BODY_HH

#include "utils/alg.hh"
#include <memory>

namespace d6 {

class LevelSet ;

class RigidBody
{
public:

	EIGEN_MAKE_ALIGNED_OPERATOR_NEW

	typedef Eigen::Matrix< Scalar, 6, 1, Eigen::DontAlign > Vel ;

	typedef typename Vel::     FixedSegmentReturnType< 3 >::Type      VelComp ;
	typedef typename Vel::ConstFixedSegmentReturnType< 3 >::Type ConstVelComp ;

	explicit RigidBody( std::unique_ptr< LevelSet > &ls ) ;

	Vec velocity_at( const Vec& x ) const ;

	const LevelSet& levelSet() const
	{ return *m_levelSet ; }
	const LevelSet* levelSetPtr() const
	{ return m_levelSet.get() ; }

	void set_velocity( const Vec& vel, const Vec& angularVel )
	{
		velocity() = vel ;
		angularVelocity() = angularVel ;
	}

	void predict_velocity( const Scalar dt, const Vec6& forces ) const ;

	void move( const Scalar dt ) const ;

	const Vel &velocities() const {
		return m_velocity ;
	}


	VelComp velocity() { return m_velocity.head<3>() ; }
	VelComp angularVelocity() { return m_velocity.tail<3>() ; }

	ConstVelComp velocity() const { return m_velocity.head<3>() ; }
	ConstVelComp angularVelocity() const { return m_velocity.tail<3>() ; }

private:
	std::unique_ptr< LevelSet > m_levelSet ;

	Vel m_velocity ;

} ;

} //d6


#endif
