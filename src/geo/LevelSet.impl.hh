#ifndef D6_LEVEL_SET_IMPL_HH
#define D6_LEVEL_SET_IMPL_HH

#include "LevelSet.hh"

namespace d6 {

struct SphereLevelSet : public LevelSet
{
	virtual Scalar eval_local(const Vec &x) const {
		return 1. - x.norm() ;
	}

	virtual Vec grad_local(const Vec &x) const {
		return -x / ( 1.e-12 + x.norm() ) ;
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

	template<class Archive>
	void serialize(Archive &ar, const unsigned int version ) ;
};

} //ns d6

#endif
