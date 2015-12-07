#ifndef D6_LEVEL_SET_HH
#define D6_LEVEL_SET_HH

#include "utils/alg.hh"

#include <memory>
#include <Eigen/Geometry>

namespace d6 {

typedef Eigen::Quaternion< Scalar > Quaternion ;

class LevelSet {

public:

	EIGEN_MAKE_ALIGNED_OPERATOR_NEW

	typedef std::unique_ptr< LevelSet > Ptr ;

	// Accessors
	Scalar eval_at( const Vec& x ) const ;
	void   grad_at( const Vec& x, Vec& grad ) const ;

	const Vec& origin() const
	{
		return m_origin ;
	}
	const Quaternion& rotation() const
	{
		return m_frame ;
	}
	Scalar scale() const {
		return m_scale ;
	}

	//Constructors
	static Ptr make_sphere( ) ;
	static Ptr make_plane( ) ;
	static Ptr make_torus( Scalar radius ) ;
	static Ptr make_cylinder( Scalar height ) ;
	static Ptr make_hole( Scalar radius ) ;
	static Ptr from_mesh( const char* objFile ) ;

	virtual bool compute() { return true ; }

	// Absolute positioning
	LevelSet& scale( const Scalar s )
	{ m_scale = s ; return *this ; }

	LevelSet& set_origin( const Vec& pos ) {
		m_origin = pos ;
		return *this ;
	}
	LevelSet& set_rotation( const Quaternion& frame ) {
		m_frame = frame ;
		return *this ;
	}
	LevelSet& set_rotation( const Vec& axis, Scalar angle ) {
		m_frame = Eigen::AngleAxis< Scalar >( angle, axis )  ;
		return *this ;
	}

	void inv_inertia( Mat66 & Mi ) const ;

	Scalar volume() const {
		return local_volume() * std::pow( m_scale, 3. ) ;
	}

	//Movement
	void move( const Vec& depl, const Quaternion& rot ) ;

	template<class Archive>
	static void register_derived(Archive &ar) ;
	template<class Archive>
	void serialize(Archive &ar, const unsigned int version) ;

protected:

	void to_local( const Vec &world, Vec &local) const ;
	void to_local_mat ( Mat &mat ) const ;

	LevelSet() ;

	virtual Scalar eval_local( const Vec&x ) const = 0 ;
	virtual Vec    grad_local( const Vec&x ) const = 0 ;
	virtual Scalar local_volume( ) const = 0 ;
	virtual void local_inv_inertia( Mat &I ) const = 0 ;

private:
	Vec 	   m_origin ;
	Quaternion m_frame ;

	Scalar 	   m_scale ;
};


} //d6


#endif
