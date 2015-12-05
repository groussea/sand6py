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
//	Quaternion fi = m_frame ; fi.w() = -fi.w() ;
	local = ( m_frame.inverse() * ( world - m_origin ) ) / m_scale  ;
}
void LevelSet::to_local_mat( Mat &mat) const
{
	mat = m_frame.inverse().matrix() / m_scale ;
}

void LevelSet::inv_inertia(Mat66 &Mi) const
{
	Mi.setIdentity() ;
	Mi.diagonal().head<3>() /= local_volume() ;

	Mat I ;
	local_inv_inertia( I );
	Mat w2l ;
	to_local_mat( w2l );

	Mi.block<3,3>(3,3) = w2l.transpose() * I * w2l ;


	Mi *= std::pow( m_scale, -3. ) ;
}

void MeshLevelSet::compute() 
{
	//FIXME
}

LevelSet::Ptr LevelSet::make_sphere() { return Ptr( new SphereLevelSet() ) ; }
LevelSet::Ptr LevelSet::make_cylinder(Scalar height) { return Ptr( new CylinderLevelSet( height ) ) ; }
LevelSet::Ptr LevelSet::make_plane() { return Ptr( new PlaneLevelSet() ) ; }
LevelSet::Ptr LevelSet::make_torus(Scalar radius) { return Ptr( new TorusLevelSet( radius ) ) ; }
LevelSet::Ptr LevelSet::from_mesh(const char* objFile) { return Ptr( new MeshLevelSet( objFile ) ) ; }

} //d6
