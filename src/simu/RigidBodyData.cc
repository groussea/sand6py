#include "RigidBodyData.hh"

#include "RigidBody.hh"

#include "FormBuilder.hh"
#include "FormBuilder.impl.hh"

#include "geo/LevelSet.hh"

#include "geo/TensorField.hh"
#include "geo/Tensor.hh"

#include "geo/Grid.hh"
#include "geo/Voxel.hh"

#include <bogus/Core/Block.impl.hpp>

namespace d6 {

RigidBodyData::RigidBodyData( RigidBody& rb_, TensorField &s )
	: rb(rb_), stresses(s)
{
}

Scalar RigidBodyData::phi( const Vec &x ) const
{
	return std::min( 1., std::max(0., 1. + rb.levelSet().eval_at( x ) ) );
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

		bool occupied = false ;
		bool boundary = false ;
		bool interior = true  ;

		for( Index k = 0 ; k < MeshType::NV ; ++k  ) {
			if( phi( geo.vertex(k) ) >  0. ) {
				occupied = true ;

				if( phi( geo.vertex(k) ) >= 1. ) {
					boundary = true ;
				} else {
					interior = false ;
				}
			}
		}

		if( occupied ) {
			occupiedCells.push_back( cell ) ;
		}

		if( boundary && !interior ) {
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

void RigidBodyData::integrate(const Active &phaseNodes, Index totNodes)
{
	typedef const typename MeshType::Location& Loc ;
	typedef const typename MeshType::Interpolation& Itp ;
	typedef const typename MeshType::Derivatives& Dcdx ;

	const Index m = phaseNodes.count() ;

	const MeshType &mesh = stresses.mesh() ;
	FormBuilder builder( mesh ) ;
	builder.reset( totNodes );
	builder.addToIndex( occupiedCells, phaseNodes.indices, phaseNodes.indices );
	builder.addToIndex(   nodes.cells,      nodes.indices, phaseNodes.indices );
	builder.makeCompressed();

	jacobian.clear() ;
	jacobian.setRows( totNodes );
	jacobian.setCols( m );
	jacobian.cloneIndex( builder.index() ) ;
	jacobian.setBlocksToZero() ;

	builder.integrate_qp( occupiedCells, [&]( Scalar w, Loc loc, Itp itp, Dcdx )
	{
		Vec dphi_dx ;
		grad_phi( mesh.pos( loc ), dphi_dx ) ;

		FormBuilder:: addUTauGphi( jacobian, w, itp, dphi_dx, phaseNodes.indices, phaseNodes.indices ) ;
	}
	);
	builder.integrate_qp(  nodes.cells, [&]( Scalar w, Loc loc, Itp itp, Dcdx )
	{
		Vec dphi_dx ;
		grad_phi( mesh.pos( loc ), dphi_dx ) ;

		FormBuilder:: addUTauGphi( jacobian, w, itp, dphi_dx,      nodes.indices, phaseNodes.indices ) ;
	}
	);

}

void RigidBodyData::assemble_matrices(const Active &phaseNodes, Index totNodes)

{
	const MeshType &mesh = stresses.mesh() ;

	typename MeshType::NodeList nodelist ;
	typename MeshType::CellGeo geo ;

	fraction.resize( nodes.count() );
	fraction.setConstant( -1 ) ;

	typename FormMat<3,6>::UncompressedType	proj ;
	proj.clear() ;
	proj.setRows( phaseNodes.count() );
	proj.setCols( 1 );

	typename FormMat<3,6>::BlockT P ;
	P.block<3,3>(0,0).setIdentity() ;

	for( const typename MeshType::Cell& cell : nodes.cells )
	{
		mesh.get_geo( cell, geo );
		mesh.list_nodes( cell, nodelist );

		for( Index k = 0 ; k < MeshType::NV ; ++k  ) {
			const Index glb_idx = nodelist[k] ;
			const Index loc_idx = nodes.indices[glb_idx] - nodes.offset ;
			const Vec pos = geo.vertex(k) ;

			if( fraction( loc_idx ) < 0. ) {
				fraction( loc_idx ) = phi( pos ) ;
				const Vec dx = pos - rb.levelSet().origin() ;
				make_cross_mat( dx, P.block<3,3>(0,3) ) ;
				proj.insert( phaseNodes.indices[ glb_idx ], 0  ) = P ;
			}
		}
	}

	// FIXME bogus add out-of-bounds insert assert
	// FIXME bogus copy with fixed block type

	proj.finalize() ;

	projection = proj.transpose() ;
//	projection.setZero() ;
//	projection.setCols( phaseNodes.count() );
//	projection.setRows( 1 );
//	projection.finalize() ;
//	projection += proj ;

	integrate( phaseNodes, totNodes );
}


} //d6
