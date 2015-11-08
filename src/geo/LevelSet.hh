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

	const Vec& origin() const {
		return m_origin ;
	}

	//Constructors
	static Ptr make_sphere( ) ;
	static Ptr make_plane( ) ;
	static Ptr make_box( const Vec &box ) ;
	static Ptr make_cylinder( Scalar height ) ;

	// Absolute positioning
	LevelSet& scale( const Scalar s )
	{ m_scale = s ; return *this ; }

	LevelSet& set_origin( const Vec& pos ) {
		m_origin = pos ;
		return *this ;
	}
	LevelSet& set_frame( const Quaternion& frame ) {
		m_frame = frame ;
		return *this ;
	}

	//Movement
	void move( const Vec& depl, const Quaternion& rot ) ;

	template<class Archive>
	static void register_derived(Archive &ar) ;
	template<class Archive>
	void serialize(Archive &ar, const unsigned int version) ;

protected:

	void to_local( const Vec &world, Vec &local) const ;

	LevelSet() ;

	virtual Scalar eval_local( const Vec&x ) const = 0 ;
	virtual Vec    grad_local( const Vec&x ) const = 0 ;

private:
	Vec 	   m_origin ;
	Quaternion m_frame ;

	Scalar 	   m_scale ;
};


} //d6


#endif
