
#include "utils/Config.hh"

#include "geo/MeshShapeFunction.hh"
#include "geo/P2ShapeFunction.hh"
#include "geo/TetGrid.hh"

#include "ScalarField.hh"

#include <gtest/gtest.h>

using namespace d6 ;

TEST( geo, p2 )
{
	VecWi dim(2,2) ;
	Vec   box(2.,2.) ;

	TetGrid g ( box,   dim ) ;
	TetGrid gd( box, 2*dim ) ;

	ASSERT_EQ(  8, g.nCells() ) ;
	ASSERT_EQ( 16, g.nEdges() ) ;

	P2< TetGrid > p2( g ) ;
	ASSERT_EQ( 25, p2.nDOF() ) ;
	ASSERT_EQ( gd.nNodes(), p2.nDOF() ) ;

	AbstractScalarField< P2< TetGrid > > f( p2 ) ;
	f.set_constant( 1. ) ;
	ASSERT_DOUBLE_EQ( 1., f(Vec(.3,.3)) ) ;
	ASSERT_DOUBLE_EQ( 1., f(Vec(1.3,.3)) ) ;
	ASSERT_DOUBLE_EQ( 1., f(Vec(.3,1.3)) ) ;
	ASSERT_DOUBLE_EQ( 1., f(Vec(1.3,1.3)) ) ;


	f.add_at( g.locate(Vec(.5,.5)), 1. );
	ASSERT_DOUBLE_EQ( 2., f(Vec(.5,.5)) ) ;
	ASSERT_DOUBLE_EQ( 1., f(Vec(1.,1.)) ) ;

	for( auto qpIt = p2.qpBegin() ; qpIt != p2.qpEnd() ; ++qpIt ) {
		TetGrid::Location loc, loc_bis, loc_ter ;
		qpIt.locate(loc) ;
		p2.locate( qpIt.pos(), loc_bis ) ;

		ASSERT_TRUE( loc.cell.isApprox( loc_bis.cell ) ) ;
		ASSERT_TRUE( loc.coords.isApprox( loc_bis.coords ) ) ;

		typename P2<TetGrid>::Interpolation itp;
		p2.interpolate( loc, itp );
		ASSERT_FLOAT_EQ( 1, itp.coeffs.sum() ) ;
		p2.interpolate_tpz( loc, itp );
		ASSERT_FLOAT_EQ( 1, itp.coeffs.sum() ) ;
		ASSERT_LE( 0, itp.coeffs.minCoeff() ) ;
		ASSERT_GE( 1, itp.coeffs.maxCoeff() ) ;

		loc_ter.coords.setZero() ;
		for( Index k = 0 ; k < itp.nodes.size() ; ++k ) {
			p2.locate_dof( loc_bis, k );
			loc_ter.coords += loc_bis.coords * itp.coeffs[k] ;
		}
		ASSERT_TRUE( loc.coords.isApprox( loc_ter.coords ) ) ;

	}

	p2.compute_lumped_mass( f.flatten() );
	ASSERT_DOUBLE_EQ( g.box().prod(), f.flatten().sum() );

	p2.compute_tpz_mass( f.flatten() );
	ASSERT_DOUBLE_EQ( g.box().prod(), f.flatten().sum() );

}
