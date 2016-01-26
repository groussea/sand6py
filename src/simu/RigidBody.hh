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

	typedef Eigen::Matrix< Scalar, SD, 1, Eigen::DontAlign > Vel ;

	typedef typename Vel::     FixedSegmentReturnType< WD >::Type         VelComp ;
	typedef typename Vel::ConstFixedSegmentReturnType< WD >::Type    ConstVelComp ;
	typedef typename Vel::     FixedSegmentReturnType< RD >::Type      AngVelComp ;
	typedef typename Vel::ConstFixedSegmentReturnType< RD >::Type ConstAngVelComp ;

	RigidBody( std::unique_ptr< LevelSet > &ls, Scalar volMass ) ;

	// Accessors

	Vec velocity_at( const Vec& x ) const ;

	const Vel &velocities() const {
		return m_velocity ;
	}

	VelComp velocity() { return m_velocity.head<WD>() ; }
	AngVelComp angularVelocity() { return m_velocity.tail<RD>() ; }

	ConstVelComp velocity() const { return m_velocity.head<WD>() ; }
	ConstAngVelComp angularVelocity() const { return m_velocity.tail<RD>() ; }

	const LevelSet& levelSet() const
	{ return *m_levelSet ; }
	const LevelSet* levelSetPtr() const
	{ return m_levelSet.get() ; }

	Scalar volumicMass() const
	{
		return m_volumicMass ;
	}

	void inv_inertia( MatS& Mi ) const ;

	// Modifiers

	void set_velocity( const Vec& vel, const VecR& angularVel )
	{
		velocity() = vel ;
		angularVelocity() = angularVel ;
	}

	void integrate_gravity( const Scalar dt, const Vec& gravity ) ;
	void integrate_forces( const Scalar dt, const VecS& forces ) ;

	void move( const Scalar dt ) const ;

	void move_to( const Vec& pos ) const ;

private:
	std::unique_ptr< LevelSet > m_levelSet ;

	Scalar m_volumicMass ;

	Vel m_velocity ;

} ;

} //d6


#endif
