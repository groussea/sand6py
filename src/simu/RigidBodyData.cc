#include "RigidBodyData.hh"

#include "RigidBody.hh"

#include "FormBuilder.hh"
#include "FormBuilder.impl.hh"

#include "geo/LevelSet.hh"

#include "geo/TensorField.hh"
#include "geo/Tensor.hh"

#include "geo/MeshImpl.hh"
#include "geo/MeshCell.hh"

#include <bogus/Core/Block.impl.hpp>

namespace d6 {

const Scalar RigidBodyData::s_splatRad = 1 ;

RigidBodyData::RigidBodyData( RigidBody& rb_, TensorField &s )
	: rb(rb_), stresses(s)
{
}

Scalar RigidBodyData::phi( const Vec &x ) const
{
	return std::min( 1., std::max(0., 1. + rb.levelSet().eval_at( x ) / s_splatRad ) );
}
void RigidBodyData::grad_phi(const Vec &x, Vec &grad) const
{
	rb.levelSet().grad_at( x, grad );
	grad /= s_splatRad ;
}

void RigidBodyData::compute_active( const Active& phaseNodes, BoundaryConditions& bc )
{
	const TensorField::ShapeFunc &shape = stresses.shape() ;
	const MeshType &mesh = shape.derived().mesh() ;

	typename Linear<MeshImpl>::NodeList nodelist ;
	typename MeshType::CellGeo geo ;

	nodes.reset( mesh.nNodes() );

	for( const typename MeshType::Cell& cell : phaseNodes.cells )
	{
		mesh.get_geo( cell, geo );

		bool occupied = false ;
		bool boundary = false ;
		bool interior = true  ;

		for( Index k = 0 ; k < MeshType::NV ; ++k  ) {
			const Scalar phi_k = phi( geo.vertex(k) ) ;

			if( phi_k >= 1. ) {
				boundary = true ;
			} else {
				if( phi_k >  0. ) {
					occupied = true ;
				}
				interior = false ;
			}
		}

		if( occupied || boundary ) {
			occupiedCells.push_back( cell ) ;
		}

		if( boundary && !interior )
		{
			nodes.cells.push_back( cell ) ;
			typename TensorField::ShapeFunc::Location loc ; loc.cell = cell ;
			shape.list_nodes( loc, nodelist );

			for( Index k = 0 ; k < Linear<MeshImpl>::NI ; ++k  ) {
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
	typedef const typename Linear<MeshImpl>::Location& Loc ;
	typedef const typename Linear<MeshImpl>::Interpolation& Itp ;
	typedef const typename Linear<MeshImpl>::Derivatives& Dcdx ;

	const Index m = phaseNodes.count() ;

	const TensorField::ShapeFunc &shape = stresses.shape() ;
	const MeshType &mesh = shape.derived().mesh() ;

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

		FormBuilder:: addUTaunGphi( jacobian, w, itp, dphi_dx, phaseNodes.indices, phaseNodes.indices ) ;
	}
	);
	builder.integrate_node(  nodes.cells, [&]( Scalar w, Loc loc, Itp itp )
	{
		Vec dphi_dx ;
		grad_phi( mesh.pos( loc ), dphi_dx ) ;

		FormBuilder:: addUTauGphi ( jacobian, w, itp, dphi_dx,      nodes.indices, phaseNodes.indices ) ;
	}
	);

}

void RigidBodyData::assemble_matrices(const Active &phaseNodes, Index totNodes)

{
	const TensorField::ShapeFunc &shape = stresses.shape() ;
	const MeshType &mesh = shape.derived().mesh() ;

	typename Linear<MeshImpl>::NodeList nodelist ;
	typename MeshType::CellGeo geo ;

	fraction.resize( nodes.count() );
	fraction.setConstant( -1 ) ;

	typename FormMat<WD,SD>::UncompressedType	proj ;
	proj.clear() ;
	proj.setRows( phaseNodes.count() );
	proj.setCols( 1 );

	typename FormMat<WD,SD>::BlockT P ;
	P.block<WD,WD>(0,0).setIdentity() ;

	for( const typename MeshType::Cell& cell : nodes.cells )
	{
		mesh.get_geo( cell, geo );
		typename TensorField::ShapeFunc::Location loc ; loc.cell = cell ;
		shape.list_nodes( loc, nodelist );

		for( Index k = 0 ; k < Linear<MeshImpl>::NI ; ++k  ) {
			const Index glb_idx = nodelist[k] ;
			const Index loc_idx = nodes.indices[glb_idx] - nodes.offset ;
			const Vec pos = geo.vertex(k) ;

			if( fraction( loc_idx ) < 0. ) {
				fraction( loc_idx ) = phi( pos ) ;
				const Vec dx = pos - rb.levelSet().origin() ;
				make_cross_mat( dx, P.block<WD,RD>(0,WD) ) ;
				proj.insert( phaseNodes.indices[ glb_idx ], 0  ) = P ;
			}
		}
	}


	proj.finalize() ;

	projection = proj.transpose() ;

	integrate( phaseNodes, totNodes );
}


} //d6
