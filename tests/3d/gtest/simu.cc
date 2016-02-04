
#include <gtest/gtest.h>

#include "geo/MeshImpl.hh"
#include "simu/FormBuilder.hh"
#include "simu/FormBuilder.impl.hh"

#include <bogus/Core/Block.impl.hpp>

using namespace d6 ;



TEST( simu, quad ) {

	Vec3i dim( 1, 1, 1 ) ;
	Vec   box( 2, 6, .3 ) ;
	MeshImpl g( box, dim ) ;
	Linear<MeshImpl> shape(g) ;

	MeshType::Cells cells ;
	for( MeshType::CellIterator it = g.cellBegin() ; it != g.cellEnd() ; ++it )
		cells.push_back( * it ) ;

	ASSERT_EQ( (Index) cells.size(), g.nCells() ) ;

//	std::vector< Index > indices ( g.nNodes() ) ;
//	std::iota( indices.begin(), indices.end(), 0 );
//	for( unsigned i = 0 ; i < indices.size() ; ++i ) std::cout << indices[i] << " " ;

	typedef const typename Linear<MeshImpl>::Interpolation& Itp ;
	typedef const typename Linear<MeshImpl>::Derivatives& Dcdx ;

	Scalar f_cst   = 0. ;
	Scalar f_lin   = 0. ;
	Scalar f_quad  = 0. ;
	Scalar f_quad2 = 0. ;

	typedef FormBuilder< Linear<MeshImpl>, Linear<MeshImpl> > Builder ;
	Builder builder( g.shaped<Linear>(), g.shaped<Linear>() ) ;
	builder.integrate_cell<form::Left>( cells.begin(), cells.end(), [&]( Scalar w, const Vec&pos, Itp, Dcdx, Itp, Dcdx ) {
		f_cst += w ;
		f_lin +=   w * pos.prod() ;
		f_quad  += w * pos.prod()  * pos.prod() ;
		f_quad2 += w * pos.dot( pos ) ;

	} ) ;

	ASSERT_FLOAT_EQ( g.box().prod() , f_cst  ) ;
	ASSERT_FLOAT_EQ( g.box().prod() * g.box().prod() / 8 , f_lin  ) ;
	EXPECT_FLOAT_EQ( g.box().prod() * g.box().prod() * g.box().prod() / 27,  f_quad ) ;
	EXPECT_FLOAT_EQ( (std::pow(g.box()[0],3)*g.box()[2]*g.box()[1]
		+ std::pow(g.box()[1],3)*g.box()[2]*g.box()[0]
		+ std::pow(g.box()[2],3)*g.box()[0]*g.box()[1]) / 3,  f_quad2 ) ;
}

#if (D6_MESH_IMPL == D6_MESH_GRID)
TEST( simu, DuDv ) {
	FormMat<3,3>::Type A ;
	A.setRows(2); A.setCols(2);
	A.insertBack(0,0).setZero() ;
	A.insertBack(0,1).setZero() ;
	A.insertBack(1,0).setZero() ;
	A.insertBack(1,1).setZero() ;
	A.finalize();

	std::vector< Index > indices( 8, 1 ) ;

	typename MeshType::Location loc ;
	typename Linear<MeshImpl>::Interpolation itp ;
	typename Linear<MeshImpl>::Derivatives dc_dx ;

	MeshImpl g( Vec(1,1,1), Vec3i(1,1,1) ) ;
	Linear<MeshImpl> shape(g) ;
	typedef FormBuilder< Linear<MeshImpl>, Linear<MeshImpl> > Builder ;

	MeshType::CellGeo vx ;
	g.get_geo( *g.cellBegin(), vx );
	loc.cell = *g.cellBegin() ;

	typename MeshType::CellGeo::QuadPoints  qp ;
	typename MeshType::CellGeo::QuadWeights qp_w ;

	vx.get_qp( qp, qp_w );

	for( int q = 0 ; q < MeshType::CellGeo::NQ ; ++q ) {
		loc.coords = qp.col(q) ;

		shape.interpolate( loc, itp );
		shape.get_derivatives( loc, dc_dx );

		indices[ itp.nodes[q] ] = 0 ;
		Builder::addDuDv( A, qp_w[q], itp, dc_dx, itp, dc_dx, indices, indices ) ;
		indices[ itp.nodes[q] ] = 1 ;
	}

	ASSERT_TRUE( ( A.block(0) - A.block(0).trace()/3 * Mat::Identity() ).isZero() ) ;
}
#endif
