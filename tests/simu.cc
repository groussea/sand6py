
#include <gtest/gtest.h>

#include "utils/string.hh"
#include "utils/File.hh"

#include "simu/Config.hh"

#include "geo/Grid.hh"
#include "simu/FormBuilder.hh"
#include "simu/FormBuilder.impl.hh"

using namespace d6 ;

TEST( simu, config )
{
	const std::string &test_scene =
			FileInfo(__FILE__).parentDirectory().filePath("../scenes/test.conf") ;

	Config c ;
	ASSERT_TRUE( c.from_file( test_scene ) ) ;

	ASSERT_EQ( c.res, Vec3i(1, 5, 10) ) ;
	ASSERT_TRUE( c.gravity.isApprox( Vec(40,0,0) ) ) ;
	ASSERT_TRUE( c.box.isApprox( Vec(10,20,30) ) ) ;
	ASSERT_DOUBLE_EQ( c.dt(), 1.e-2 ) ;

	const Scalar Re = ( c.box.minCoeff() * std::sqrt( c.box.minCoeff() * c.gravity.norm()) * c.volMass ) / c.viscosity ;
	ASSERT_DOUBLE_EQ( 2.e8, Re ) ;

	c.internalize();

	ASSERT_TRUE( c.box.isApprox( Vec(1,2,3) ) ) ;
	ASSERT_TRUE( c.gravity.isApprox( Vec(1,0,0) ) ) ;
	ASSERT_DOUBLE_EQ( c.units().U, 20. ) ;
	ASSERT_DOUBLE_EQ( c.units().T, .5 ) ;
	ASSERT_DOUBLE_EQ( c.units().P, 4.e5 ) ;
	ASSERT_DOUBLE_EQ( c.units().M, 2.e5 ) ;
	ASSERT_DOUBLE_EQ( c.units().M, c.units().R * c.units().U * c.units().L ) ;

	ASSERT_DOUBLE_EQ( c.viscosity, 1./Re ) ;

}


TEST( simu, quad ) {

	Vec3i dim( 10, 5, 15 ) ;
	Vec   box( 1, 1, 2 ) ;
	Grid g( box, dim ) ;

	Grid::Cells cells ;
	for( Grid::CellIterator it = g.cellBegin() ; it != g.cellEnd() ; ++it )
		cells.push_back( * it ) ;

	ASSERT_EQ( cells.size(), g.nCells() ) ;

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
