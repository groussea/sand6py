#ifndef D6_LEVEL_SET_IMPL_HH
#define D6_LEVEL_SET_IMPL_HH

#include "LevelSet.hh"

#include <limits>

namespace d6 {

struct SphereLevelSet : public LevelSet
{
	virtual Scalar eval_local(const Vec &x) const {
		return 1. - x.norm() ;
	}

	virtual Vec grad_local(const Vec &x) const {
		return -x / ( 1.e-12 + x.norm() ) ;
	}

	virtual void local_inv_inertia( Mat& I ) const {
		I = Mat::Identity() / ( local_volume() * 2./5. ) ;
	}

	virtual Scalar local_volume() const {
		return 4./3 * M_PI ;
	}

	template<class Archive>
	void serialize(Archive &ar, const unsigned int version ) ;
};
struct PlaneLevelSet : public LevelSet
{
	virtual Scalar eval_local(const Vec &x) const {
		return - x[2] ;
	}

	virtual Vec grad_local(const Vec & ) const {
		return Vec(0, 0, -1) ;
	}

	virtual void local_inv_inertia( Mat& I ) const {
		I.setZero() ;
	}

	virtual Scalar local_volume() const {
		return std::numeric_limits<Scalar>::infinity() ;
	}

	template<class Archive>
	void serialize(Archive &ar, const unsigned int version ) ;
};

struct TorusLevelSet : public LevelSet
{

	explicit TorusLevelSet( Scalar rad = 1. )
		: m_radius( rad )
	{}

	Scalar radius() const { return m_radius ; }

	virtual Scalar eval_local(const Vec &x) const
	{
		Vec proj ;
		proj_on_circle( x, proj );

		return m_radius - (x - proj).norm() ;
	}

	virtual Vec grad_local(const Vec &x) const {
		Vec proj ;
		proj_on_circle( x, proj );

		const Scalar n = (x - proj).norm() ;
		return (proj - x) / (1.e-12 + n) ;
	}

	virtual void local_inv_inertia( Mat& I ) const {
		//From http://mathworld.wolfram.com/Torus.html
		I.setZero() ;
		I(0,0) = I(1,1) = 1./( 5./8 * m_radius * m_radius + .5 ) ;
		I(2,2)          = 1./( 3./4 * m_radius * m_radius + 1. ) ;
		I.diagonal() /= local_volume() ;
	}

	virtual Scalar local_volume() const {
		return 2 * M_PI * M_PI * m_radius * m_radius ;
	}

	template<class Archive>
	void serialize(Archive &ar, const unsigned int version ) ;

private:

	static void proj_on_circle( const Vec&x, Vec& proj )
	{
		const Scalar n = x.head<2>().norm() ;

		// Projection on unit circle on XY axis
		if( n < 1.e-12 ) proj = Vec(0,1,0) ;
		else {
			proj.head<2>() = x.head<2>() / n ;
			proj[2] = 0 ;
		}
	}

	Scalar m_radius ;
};

struct CylinderLevelSet : public LevelSet
{

	explicit CylinderLevelSet( Scalar h = 1. )
		: m_height( h )
	{}

	Scalar height() const { return m_height ; }

	virtual Scalar eval_local(const Vec &x) const
	{
		Vec proj ;
		proj_on_axis( x, proj );

		return 1. - (x - proj).norm() ;
	}

	virtual Vec grad_local(const Vec &x) const {
		Vec proj ;
		proj_on_axis( x, proj );

		const Scalar n = (x - proj).norm() ;
		return (proj - x) / (1.e-12 + n) ;
	}

	virtual void local_inv_inertia( Mat& I ) const {
		//From http://mathworld.wolfram.com/Torus.html
		I.setZero() ;
		I(0,0) = I(1,1) = 6./( 3. + m_height * m_height ) ;
		I(2,2)          = 1. ;
		I.diagonal() *= 2./local_volume() ;
	}

	virtual Scalar local_volume() const {
		return m_height * M_PI ;
	}

	template<class Archive>
	void serialize(Archive &ar, const unsigned int version ) ;

private:

	void proj_on_axis( const Vec&x, Vec& proj ) const
	{
		proj[0] = 0. ;
		proj[1] = 0. ;
		proj[2] = std::min(-.5*m_height, std::max(.5*m_height, x[2])) ;
	}

	Scalar m_height ;
};

} //ns d6

#endif
