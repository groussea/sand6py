#ifndef D6_LEVEL_SET_HH
#define D6_LEVEL_SET_HH

#include "utils/alg.hh"

#include <memory>
#include <Eigen/Geometry>

namespace d6 {

template< Index Dim >
struct RotationTraits
{
	typedef Scalar Type ;
};
template< >
struct RotationTraits<3>
{
	typedef Eigen::Quaternion< Scalar, Eigen::DontAlign > Type ;
};
typedef typename RotationTraits<WD>::Type Rotation ;


class LevelSet {

public:

	typedef std::unique_ptr< LevelSet > Ptr ;

	// Accessors
	Scalar eval_at( const Vec& x ) const ;
	void   grad_at( const Vec& x, Vec& grad ) const ;

	const Vec& origin() const
	{
		return m_origin ;
	}
	const Rotation& rotation() const
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
	static Ptr make_hourglass( Scalar height, Scalar radius ) ;
	static Ptr from_mesh( const char* objFile ) ;

	virtual bool compute() { return true ; }

	// Absolute positioning
	LevelSet& scale( const Scalar s )
	{ m_scale = s ; return *this ; }

	LevelSet& set_origin( const Vec& pos ) {
		m_origin = pos ;
		return *this ;
	}
	LevelSet& set_rotation( const Rotation& frame ) {
		m_frame = frame ;
		return *this ;
	}
	LevelSet& rotate( const Rotation& frame ) {
		m_frame = frame * m_frame ;
		return *this ;
	}

#if D6_DIM==3
	LevelSet& set_rotation( const Vec& axis, Scalar angle ) {
		m_frame = Eigen::AngleAxis< Scalar >( angle, axis )  ;
		return *this ;
	}
	LevelSet& rotate( const Vec& axis, Scalar angle ) {
		m_frame = Eigen::AngleAxis< Scalar >( angle, axis ) * m_frame  ;
		return *this ;
	}
#endif

	void inv_inertia( MatS & Mi ) const ;

	Scalar volume() const {
		return local_volume() * std::pow( m_scale, WD ) ;
	}

	//Movement
	void move( const Vec& depl, const Rotation& rot ) ;

	template<class Archive>
	static void register_derived(Archive &ar) ;
	template<class Archive>
	void serialize(Archive &ar, const unsigned int version) ;

		virtual ~LevelSet() {}

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
	Rotation m_frame ;

	Scalar 	   m_scale ;
};


} //d6


#endif
