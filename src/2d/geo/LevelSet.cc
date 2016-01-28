#include "LevelSet.hh"

#include "LevelSet.impl.hh"

#include <iostream>

namespace d6 {

LevelSet::LevelSet()
	: m_origin(Vec::Zero()), m_frame( 0 )
{
	scale(1.) ;
}

void LevelSet::move(const Vec &depl, const Rotation &rot)
{
	m_origin += depl ;
	m_frame   = rot + m_frame ;
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
	grad = d6::rotate( m_frame, ( grad_local( loc ).array() ) ) ;
}

void LevelSet::to_local(const Vec &world, Vec &local) const
{
//	Quaternion fi = m_frame ; fi.w() = -fi.w() ;
	local = d6::rotate( -m_frame, ( world - m_origin ) ) / m_scale  ;
}
void LevelSet::to_local_mat(MatR &mat) const
{
	mat(0,0) = 1. / m_scale ;
}

void LevelSet::inv_inertia(MatS &Mi) const
{
	Mi.setIdentity() ;
	Mi.diagonal().head<WD>() /= local_volume() ;

	MatR I ;
	local_inv_inertia( I );

	Mi.block<RD,RD>(WD,WD) = I/m_scale/m_scale ;

	Mi *= std::pow( m_scale, -WD ) ;
}

LevelSet::Ptr LevelSet::make_sphere()   { return Ptr( new SphereLevelSet() ) ; }
LevelSet::Ptr LevelSet::make_plane()    { return Ptr( new PlaneLevelSet() ) ; }
LevelSet::Ptr LevelSet::make_cylinder(const Scalar len) { return Ptr( new SegLevelSet(len) ) ; }

} //d6
