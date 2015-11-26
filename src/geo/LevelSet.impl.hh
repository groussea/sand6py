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

} //ns d6

#endif
