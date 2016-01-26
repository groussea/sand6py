#ifndef D6_LEVEL_SET_IMPL_HH
#define D6_LEVEL_SET_IMPL_HH

#include "LevelSet.hh"

#include "Grid.hh"
#include "ScalarField.hh"

#include <limits>

namespace d6 {

struct SphereLevelSet : public LevelSet
{
	Scalar eval_local(const Vec &x) const override {
		return 1. - x.norm() ;
	}

	Vec grad_local(const Vec &x) const override {
		return -x / ( 1.e-12 + x.norm() ) ;
	}

	void local_inv_inertia( MatR& I ) const override {
		I(0,0) = 2. / local_volume() ;
	}

	Scalar local_volume() const override {
		return M_PI ;
	}

	template<class Archive>
	void serialize(Archive &ar, const unsigned int version ) ;
};
struct PlaneLevelSet : public LevelSet
{
	Scalar eval_local(const Vec &x) const override {
		return - x[1] ;
	}

	Vec grad_local(const Vec & ) const override {
		return Vec(0, -1) ;
	}

	void local_inv_inertia( MatR& I ) const override {
		I.setZero() ;
	}

	Scalar local_volume() const override {
		return std::numeric_limits<Scalar>::infinity() ;
	}

	template<class Archive>
	void serialize(Archive &ar, const unsigned int version ) ;
};


} //ns d6

#endif
