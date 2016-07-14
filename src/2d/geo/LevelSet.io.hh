/*
 * This file is part of Sand6, a C++ continuum-based granular simulator.
 *
 * Copyright 2016 Gilles Daviet <gilles.daviet@inria.fr> (Inria - Université Grenoble Alpes)
 *
 * Sand6 is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.

 * Sand6 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.

 * You should have received a copy of the GNU General Public License
 * along with Sand6.  If not, see <http://www.gnu.org/licenses/>.
*/

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
void SegLevelSet::serialize(Archive &ar, const unsigned int )
{
	ar & boost::serialization::base_object<LevelSet>(*this);
	ar & m_len ;
}

// register derived class ptrs
template<class Archive>
void LevelSet::register_derived(Archive &ar )
{
	ar.register_type(static_cast<   SphereLevelSet *>(NULL));
	ar.register_type(static_cast<    PlaneLevelSet *>(NULL));
	ar.register_type(static_cast<      SegLevelSet *>(NULL));
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
