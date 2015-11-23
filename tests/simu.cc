
#include <gtest/gtest.h>

#include "geo/Grid.hh"
#include "simu/FormBuilder.hh"
#include "simu/FormBuilder.impl.hh"

#include <bogus/Core/Block.impl.hpp>

using namespace d6 ;



TEST( simu, quad ) {

	Vec3i dim( 10, 5, 15 ) ;
	Vec   box( 12, 17, 2 ) ;
	Grid g( box, dim ) ;

	Grid::Cells cells ;
	for( Grid::CellIterator it = g.cellBegin() ; it != g.cellEnd() ; ++it )
		cells.push_back( * it ) ;

	ASSERT_EQ( (Index) cells.size(), g.nCells() ) ;

//	std::vector< Index > indices ( g.nNodes() ) ;
//	std::iota( indices.begin(), indices.end(), 0 );
//	for( unsigned i = 0 ; i < indices.size() ; ++i ) std::cout << indices[i] << " " ;

	typedef const typename Grid::Location& Loc ;
	typedef const typename Grid::Interpolation& Itp ;
	typedef const typename Grid::Derivatives& Dcdx ;

	Scalar f_cst  = 0. ;
	Scalar f_lin  = 0. ;
	Scalar f_quad = 0. ;

	FormBuilder builder( g ) ;
	builder.integrate_qp( cells, [&]( Scalar w, Loc loc, Itp , Dcdx ) {
		Grid::CellGeo geo ;
		g.get_geo( loc.cell, geo );
		f_cst += w ;
		f_lin += w * geo.pos( loc.coords ).prod() ;
		f_quad += w * geo.pos( loc.coords ).prod()  * geo.pos( loc.coords ).prod() ;
	} ) ;

	ASSERT_FLOAT_EQ( g.box().prod() , f_cst  ) ;
	ASSERT_FLOAT_EQ( g.box().prod() * g.box().prod() / 8 , f_lin  ) ;
	ASSERT_FLOAT_EQ( g.box().prod() * g.box().prod() * g.box().prod() / 27,  f_quad ) ;
}

TEST( simu, DuDv ) {
	FormMat<3,3>::Type A ;
	A.setRows(2); A.setCols(2);
	A.insertBack(0,0).setZero() ;
	A.insertBack(0,1).setZero() ;
	A.insertBack(1,0).setZero() ;
	A.insertBack(1,1).setZero() ;
	A.finalize();

	std::vector< Index > indices( 8, 1 ) ;

	typename Grid::Location loc ;
	typename Grid::Interpolation itp ;
	typename Grid::Derivatives dc_dx ;

	Grid g( Vec(1,1,1), Vec3i(1,1,1) ) ;
	Voxel vx ;
	g.get_geo( *g.cellBegin(), vx );
	loc.cell = *g.cellBegin() ;

	typename Voxel::QuadPoints  qp ;
	typename Voxel::QuadWeights qp_w ;

	vx.get_qp( qp, qp_w );

	for( int q = 0 ; q < Voxel::NQ ; ++q ) {
		loc.coords = qp.col(q) ;

		g.interpolate( loc, itp );
		g.get_derivatives( loc, dc_dx );

		indices[ itp.nodes[q] ] = 0 ;
		FormBuilder::addDuDv( A, qp_w[q], itp, dc_dx, indices, indices ) ;
		indices[ itp.nodes[q] ] = 1 ;
	}

	ASSERT_TRUE( ( A.block(0) - A.block(0).trace()/3 * Mat::Identity() ).isZero() ) ;
}
