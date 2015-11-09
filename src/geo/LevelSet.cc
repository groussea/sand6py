#include "LevelSet.hh"

#include "LevelSet.impl.hh"

#include <iostream>

namespace d6 {

LevelSet::LevelSet()
	: m_origin(Vec::Zero()), m_frame( Quaternion::Identity() )
{
	scale(1.) ;
}

void LevelSet::move(const Vec &depl, const Quaternion &rot)
{
	m_origin += depl ;
	m_frame   = rot * m_frame ;
}

Scalar LevelSet::eval_at(const Vec &x) const
{
	Vec loc;
	to_local(x, loc);
	return m_scale * eval_local( loc ) ;
}

void LevelSet::grad_at(const Vec &x, Vec &grad) const
{
	Vec loc;
	to_local(x, loc);
	grad = m_frame * ( grad_local( loc ).array() ) ;
}

void LevelSet::to_local(const Vec &world, Vec &local) const
{
	const Quaternion &fi = m_frame.inverse() ;
//	Quaternion fi = m_frame ; fi.w() = -fi.w() ;
	local = ( fi * ( world - m_origin ) ) / m_scale  ;
}

//LevelSet::Ptr LevelSet::make_cube() { return Ptr( new CubeLevelSet() ) ; }
LevelSet::Ptr LevelSet::make_sphere() { return Ptr( new SphereLevelSet() ) ; }
//LevelSet::Ptr LevelSet::make_cylinder() { return Ptr( new CylinderLevelSet() ) ; }
LevelSet::Ptr LevelSet::make_plane() { return Ptr( new PlaneLevelSet() ) ; }

} //d6
