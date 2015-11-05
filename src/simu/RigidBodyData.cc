#include "RigidBodyData.hh"

#include "RigidBody.hh"

#include "geo/LevelSet.hh"

#include "geo/TensorField.hh"

#include "geo/Grid.hh"
#include "geo/Voxel.hh"

namespace d6 {


Scalar RigidBodyData::phi( const Vec &x ) const
{
	return 1. + rb.levelSet().eval_at( x ) ;
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


} //d6
