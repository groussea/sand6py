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

void RigidBodyData::compute_active( const Active& phaseNodes )
{
	const TensorField::ShapeFuncImpl &shape = stresses.shape() ;
	const MeshType &mesh = shape.mesh() ;

	typename TensorField::ShapeFuncType::NodeList nodelist ;
	typename MeshType::CellGeo geo ;

	nodes.reset( shape.nDOF() );

	for( const typename MeshType::Cell& cell : phaseNodes.cells )
	{
		mesh.get_geo( cell, geo );

		bool boundary = false ;
		bool interior = true  ;

		for( Index k = 0 ; k < MeshType::NV ; ++k  ) {
			const Scalar phi_k = phi( geo.vertex(k) ) ;

			if( phi_k >= 1. ) {
				boundary = true ;
			} else {
				interior = false ;
			}
		}

		if( boundary && !interior )
		{
			nodes.cells.push_back( cell ) ;
			shape.list_nodes( cell, nodelist );

			for( Index k = 0 ; k < nodelist.rows() ; ++k  ) {
				if( nodes.indices[ nodelist[k] ] == Active::s_Inactive ) {
					nodes.indices[ nodelist[k] ] = nodes.nNodes++ ;
				}
			}
		}
	}
}

void RigidBodyData::integrate(const PrimalShape& primalShape, const DualShape& dualShape,
							  const Active &primalNodes, const Active& dualNodes,
							  Index totNodes)
{
	//FIXME other approxes
	typedef typename RBStresses::ShapeFuncImpl RBShapeFunc ;

	static_assert( std::is_same< PrimalShape, RBShapeFunc >::value,
			"Different RB and primal shape func not allowed yet" ) ;

	typedef const typename PrimalShape::Interpolation& P_Itp ;
	typedef const typename PrimalShape::Derivatives&   P_Dcdx ;
	typedef const typename DualShape  ::Interpolation& D_Itp ;
	typedef const typename DualShape  ::Derivatives&   D_Dcdx ;

	const Index m = primalNodes.count() ;

	{
		typedef FormBuilder< DualShape, PrimalShape > Builder ;
		Builder builder( dualShape, primalShape ) ;

		builder.reset( totNodes );
		builder.addToIndexIf( dualNodes.indices, primalNodes.indices,
							  [&](const Vec& pos){ return phi( pos ) > 0 ; } );
		builder.makeCompressed();

		jacobian.clear() ;
		jacobian.setRows( totNodes );
		jacobian.setCols( m );
		jacobian.cloneIndex( builder.index() ) ;
		jacobian.setBlocksToZero() ;

		builder.integrate_qp( [&]( Scalar w, const Vec& pos, D_Itp l_itp, D_Dcdx, P_Itp r_itp, P_Dcdx )
		{
			if( phi( pos ) > 0 ) {
				Vec dphi_dx ;
				grad_phi( pos, dphi_dx ) ;

				Builder:: addUTaunGphi( jacobian, w, l_itp, r_itp, dphi_dx, dualNodes.indices, primalNodes.indices ) ;
			}
		}
		);

		// Fraction
		intFraction.resize( dualNodes.count() ) ;
		intFraction.setZero() ;
		// TODO do not process all cells -- or factorize for all rbs
		for( auto qpIt = dualShape.qpBegin() ; qpIt != dualShape.qpEnd() ; ++qpIt  ) {
			typename DualShape::Location loc ;
			typename DualShape::Interpolation itp ;
			qpIt.locate( loc ) ;
			dualShape.interpolate( loc, itp ) ;

			const Scalar phiRb = std::max( 0., std::min(1., rb.levelSet().eval_at( qpIt.pos() ) ) ) ;

			for( Index k = 0 ; k < itp.nodes.rows() ; ++ k ) {
				const Index idx = dualNodes.indices[ itp.nodes[k] ] ;
				if( idx != Active::s_Inactive ) {
					intFraction[ idx ] += qpIt.weight() * itp.coeffs[k] * phiRb ;
				}
			}
		}
	}

	{
		FormMat<SD,WD>::Type jacobian_2 ;

		typedef FormBuilder< RBShapeFunc, PrimalShape > Builder ;
		Builder builder( primalShape, primalShape ) ;

		builder.reset( totNodes );
		builder.addToIndex(   nodes.cells.begin(),  nodes.cells.end(),      nodes.indices, primalNodes.indices );
		builder.makeCompressed();

		jacobian_2.clear() ;
		jacobian_2.setRows( totNodes );
		jacobian_2.setCols( m );
		jacobian_2.cloneIndex( builder.index() ) ;
		jacobian_2.setBlocksToZero() ;

		builder.integrate_node(  nodes.cells.begin(), nodes.cells.end(), [&]( Scalar w, const Vec& pos, P_Itp l_itp, P_Itp r_itp )
		{
			Vec dphi_dx ;
			grad_phi( pos, dphi_dx ) ;

			Builder:: addUTauGphi ( jacobian_2, w, l_itp, r_itp, dphi_dx,      nodes.indices, primalNodes.indices ) ;
		}
		);

		jacobian += jacobian_2 ;

	}

}

void RigidBodyData::assemble_matrices( const PrimalShape& primalShape, const DualShape& dualShape,
									   const Active &primalNodes, const Active& dualNodes,
									   Index totNodes)

{
	typedef typename PrimalShape::MeshType MeshType ;
	const MeshType &mesh = primalShape.mesh() ;

	typename TensorField::ShapeFuncType::NodeList nodelist ;

	typename FormMat<WD,SD>::UncompressedType	proj ;
	proj.clear() ;
	proj.setRows( primalNodes.count() );
	proj.setCols( 1 );

	typename FormMat<WD,SD>::BlockT P ;
	P.topLeftCorner<WD,WD>().setIdentity() ;

	for( const typename MeshType::Cell& cell : nodes.cells )
	{
		typename PrimalShape::Location ploc ;
		ploc.cell = cell ;
		primalShape.list_nodes( ploc, nodelist );

		for( Index k = 0 ; k < nodelist.rows() ; ++k  ) {

			const Index glb_idx = nodelist[k] ;
			const Index row = primalNodes.indices[ glb_idx ] ;

			if( proj.majorIndex().size(row) == 0 ) {

				primalShape.locate_dof( ploc, k ) ;
				const Vec pos = mesh.pos( ploc ) ;
				const Vec dx = pos - rb.levelSet().origin() ;
				make_cross_mat( dx, P.block<WD,RD>(0,WD) ) ;
				proj.insert( row, 0  ) = P ;
			}
		}
	}
	proj.finalize() ;

	projection = proj.transpose() ;

	integrate( primalShape, dualShape, primalNodes, dualNodes, totNodes );
}


} //d6
