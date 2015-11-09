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

	RigidBody( std::unique_ptr< LevelSet > &ls, Scalar volMass ) ;

	// Accessors

	Vec velocity_at( const Vec& x ) const ;

	const Vel &velocities() const {
		return m_velocity ;
	}

	VelComp velocity() { return m_velocity.head<3>() ; }
	VelComp angularVelocity() { return m_velocity.tail<3>() ; }

	ConstVelComp velocity() const { return m_velocity.head<3>() ; }
	ConstVelComp angularVelocity() const { return m_velocity.tail<3>() ; }

	const LevelSet& levelSet() const
	{ return *m_levelSet ; }
	const LevelSet* levelSetPtr() const
	{ return m_levelSet.get() ; }

	Scalar volumicMass() const
	{
		return m_volumicMass ;
	}

	void inv_inertia( Mat66& Mi ) const ;

	// Modifiers

	void set_velocity( const Vec& vel, const Vec& angularVel )
	{
		velocity() = vel ;
		angularVelocity() = angularVel ;
	}

	void integrate_gravity( const Scalar dt, const Vec& gravity ) ;
	void integrate_forces( const Scalar dt, const Vec6& forces ) ;

	void move( const Scalar dt ) const ;


private:
	std::unique_ptr< LevelSet > m_levelSet ;

	Scalar m_volumicMass ;

	Vel m_velocity ;

} ;

} //d6


#endif
