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
template<class Archive>
void TorusLevelSet::serialize(Archive &ar, const unsigned int )
{
	ar & boost::serialization::base_object<LevelSet>(*this);
	ar & m_radius ;
}
template<class Archive>
void HourglassLevelSet::serialize(Archive &ar, const unsigned int )
{
	ar & boost::serialization::base_object<LevelSet>(*this);
	ar & m_height ;
	ar & m_radius ;
}
template<class Archive>
void HoleLevelSet::serialize(Archive &ar, const unsigned int )
{
	ar & boost::serialization::base_object<LevelSet>(*this);
	ar & m_radius ;
}
template<class Archive>
void CylinderLevelSet::serialize(Archive &ar, const unsigned int )
{
	ar & boost::serialization::base_object<LevelSet>(*this);
	ar & m_height ;
}
template<class Archive>
void MeshLevelSet::serialize(Archive &ar, const unsigned int )
{
	ar & boost::serialization::base_object<LevelSet>(*this);
	ar & m_objFile ;
}


// register derived class ptrs
template<class Archive>
void LevelSet::register_derived(Archive &ar )
{
	ar.register_type(static_cast<   SphereLevelSet *>(NULL));
	ar.register_type(static_cast<    PlaneLevelSet *>(NULL));
	ar.register_type(static_cast<    TorusLevelSet *>(NULL));
	ar.register_type(static_cast< CylinderLevelSet *>(NULL));
	ar.register_type(static_cast<HourglassLevelSet *>(NULL));
	ar.register_type(static_cast<     MeshLevelSet *>(NULL));
	ar.register_type(static_cast<     HoleLevelSet *>(NULL));
}

//base class serilization
template<class Archive>
void LevelSet::serialize(Archive &ar, const unsigned int )
{
	ar & m_origin ;
	ar & m_scale ;
	ar & m_frame.w() & m_frame.x() & m_frame.y() & m_frame.z() ;
}


} //d6

#endif
