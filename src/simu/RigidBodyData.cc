#include "RigidBodyData.hh"

#include "RigidBody.hh"

#include "geo/LevelSet.hh"

#include "geo/TensorField.hh"

#include "geo/Grid.hh"
#include "geo/Voxel.hh"

namespace d6 {

RigidBodyData::RigidBodyData( RigidBody& rb_, TensorField &s )
	: rb(rb_), stresses(s)
{
}

Scalar RigidBodyData::phi( const Vec &x ) const
{
	return 1. + rb.levelSet().eval_at( x ) ;
}
void RigidBodyData::grad_phi(const Vec &x, Vec &grad) const
{
	rb.levelSet().grad_at( x, grad );
}

void RigidBodyData::compute_active( const Active& phaseNodes, BoundaryConditions& bc )
{
	const MeshType &mesh = stresses.mesh() ;

	typename MeshType::NodeList nodelist ;
	typename MeshType::CellGeo geo ;

	nodes.reset( mesh.nNodes() );

	for( const typename MeshType::Cell& cell : phaseNodes.cells )
	{
		mesh.get_geo( cell, geo );

		bool boundary = false ;

		for( Index k = 0 ; k < MeshType::NV ; ++k  ) {
			if( phi( geo.vertex(k) ) >= 1. ) {
//			if( phi( geo.vertex(k) ) >  0. ) {
				boundary = true ;
				break ;
			}
		}

		if( boundary ) {
			nodes.cells.push_back( cell ) ;
			mesh.list_nodes( cell, nodelist );

			for( Index k = 0 ; k < MeshType::NV ; ++k  ) {
				if( nodes.indices[ nodelist[k] ] == Active::s_Inactive ) {
					nodes.indices[ nodelist[k] ] = nodes.nNodes++ ;
					bc[ nodelist[k] ].bc = BoundaryInfo::Interior ;
				}
			}
		}
	}
}

void RigidBodyData::compute_fraction( )
{
	const MeshType &mesh = stresses.mesh() ;

	typename MeshType::NodeList nodelist ;
	typename MeshType::CellGeo geo ;

	fraction.resize( nodes.count() );

	for( const typename MeshType::Cell& cell : nodes.cells )
	{
		mesh.get_geo( cell, geo );
		mesh.list_nodes( cell, nodelist );

		for( Index k = 0 ; k < MeshType::NV ; ++k  ) {
			fraction( nodes.indices[nodelist[k]] ) = phi( geo.vertex(k) ) ;
		}
	}
}

} //d6
