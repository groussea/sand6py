#ifndef D6_LEVEL_SET_IO_HH
#define D6_LEVEL_SET_IO_HH

#include "LevelSet.impl.hh"

#include <boost/serialization/base_object.hpp>

namespace d6 {

// save/load base class information
template<class Archive>
void PlaneLevelSet::serialize(Archive &ar, const unsigned int )
{
	ar & boost::serialization::base_object<LevelSet>(*this);
}
template<class Archive>
void SphereLevelSet::serialize(Archive &ar, const unsigned int )
{
	ar & boost::serialization::base_object<LevelSet>(*this);
}

// register derived class ptrs
template<class Archive>
void LevelSet::register_derived(Archive &ar )
{
	ar.register_type(static_cast<   SphereLevelSet *>(NULL));
	ar.register_type(static_cast<    PlaneLevelSet *>(NULL));
}

//base class serilization
template<class Archive>
void LevelSet::serialize(Archive &ar, const unsigned int )
{
	ar & m_origin ;
	ar & m_scale ;
	ar & m_frame ;
}


} //d6

#endif
