#ifndef D6_LEVEL_SET_IMPL_HH
#define D6_LEVEL_SET_IMPL_HH

#include "LevelSet.hh"

#include <boost/serialization/base_object.hpp>

namespace d6 {

struct SphereLevelSet : public LevelSet
{
	virtual Scalar eval_local(const Vec &x) const {
		return 1. - x.norm() ;
	}

	virtual Vec grad_local(const Vec &x) const {
		return -x / ( 1.e-12 + x.norm() ) ;
	}
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
	void serialize(Archive &ar, const unsigned int )
	{
		// save/load base class information
		ar & boost::serialization::base_object<LevelSet>(*this);
	}
};

template<class Archive>
void LevelSet::register_derived(Archive &ar )
{
	ar.register_type(static_cast<SphereLevelSet *>(NULL));
	ar.register_type(static_cast< PlaneLevelSet *>(NULL));
}

template<class Archive>
void LevelSet::serialize(Archive &ar, const unsigned int )
{
	ar & m_origin ;
	ar & m_scale ;
	ar & m_frame.w() & m_frame.x() & m_frame.y() & m_frame.z() ;
}

} //ns d6

#endif
