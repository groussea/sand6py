
#include <gtest/gtest.h>

#include "geo/Grid.hh"

#include "geo/ScalarField.hh"
#include "geo/VectorField.hh"
#include "geo/TensorField.hh"

#include "geo/Tensor.hh"
#include "geo/Voxel.hh"

#include "geo/TetGrid.hh"
#include "geo/Tet.hh"

#include "geo/MeshShapeFunction.hh"

#include "geo/BoundaryInfo.hh"

#include <Eigen/Geometry>

using namespace d6 ;

TEST( geo, grid )
{
	Vec3i dim( 10, 5, 4 ) ;
	Vec   box( 1, 1, 1 ) ;
	Grid g( box, dim ) ;
	Linear<Grid> shape(g) ;

	ASSERT_EQ( 200, g.nCells() ) ;
	ASSERT_EQ( 11*6*5, g.nNodes() ) ;
	ASSERT_TRUE( box.isApprox( g.box() ) ) ;

	Linear<Grid>::Interpolation itp ;
	shape.interpolate( g.locate(Vec(0,0,0)), itp );

	ASSERT_EQ( 0, itp.nodes[0] ) ;
	ASSERT_EQ( dim[2]+1, itp.nodes[2] ) ;
	ASSERT_EQ( (dim[1]+1)*(dim[2]+1)+(dim[2]+1)+1, itp.nodes[7] ) ;
	ASSERT_DOUBLE_EQ( 1, itp.coeffs[0] ) ;

	shape.interpolate( g.locate(Vec(1,1,1)), itp );

	ASSERT_EQ( g.nNodes()-1, itp.nodes[7] ) ;
	ASSERT_DOUBLE_EQ( 1, itp.coeffs[7] ) ;

	Index i = 0 ;
	for( GridIterator it( g.cellBegin() ) ;  it != g.cellEnd() ; ++it, ++i )
	{
		ASSERT_EQ( i, it.index() ) ;
	}
	ASSERT_EQ( i, g.nCells() ) ;

	Voxel vx ;
	g.get_geo( Vec3i(1,2,3), vx );
	const Vec x1 = vx.center() + 1.e-2 * Vec(1,-2,3) ;

	Linear<Grid>::Derivatives dc_dx ;
	Grid::Location loc ;
	g.locate( x1, loc );
	shape.interpolate( loc, itp ) ;
	shape.get_derivatives( loc, dc_dx );

	for( unsigned i = 0 ; i < 3 ; ++i ) {
		Vec dx ( Vec::Zero() ) ;
		dx[i] = 1.e-2 ;

		Linear<Grid>::CoefList coeffs_pred = itp.coeffs + ( dc_dx * dx ) ;

		const Vec x2 = x1 + dx ;
		Linear<Grid>::Interpolation itp2 ;
		shape.interpolate( g.locate(x2), itp2 ) ;

		ASSERT_TRUE( coeffs_pred.isApprox( itp2.coeffs ) ) ;
	}

	Vec red = dc_dx.colwise().sum() ;
	ASSERT_TRUE( red.isZero() ) ;

}

TEST( geo, field )
{
	Vec3i dim( 10, 10, 10 ) ;
	Vec   box( 1, 1, 1 ) ;
	Grid g( box, dim ) ;

	Linear<Grid> shape(g) ;
	AbstractScalarField< Linear<Grid> > phi( shape ) ;
	ASSERT_EQ( g.nNodes(), phi.flatten().rows() ) ;

	phi.set_constant( 3 );
	ASSERT_DOUBLE_EQ( 3., phi( (Vec( 0.2, 0.7, 0.5 ) ) )) ;

	phi.set_zero();
	ASSERT_DOUBLE_EQ( 0., phi( (Vec( 0.3, 0.4, 0.1 ) ) )) ;
	phi.add_at( g.locate(Vec( 0.3, 0.4, 0.1 )), 1 );
	ASSERT_FLOAT_EQ( 1,    phi( (Vec( 0.3, 0.4, 0.1 ) ) )) ;
	ASSERT_FLOAT_EQ( 0.3,  phi( (Vec( 0.37, 0.4, 0.1 ) ) )) ;
	ASSERT_FLOAT_EQ( 0.7,  phi( (Vec( 0.3, 0.37, 0.1 ) ) )) ;
	ASSERT_FLOAT_EQ( 0.5,  phi( (Vec( 0.3, 0.4, 0.15 ) ) )) ;

	AbstractVectorField< Linear<Grid> > u( shape ) ;
	u.set_constant( Vec(0,1,0) ) ;
	ASSERT_TRUE( Vec(0,1,0).isApprox( u( (Vec( 0.2, 0.7, 0.5 ) ) ) )) ;

	u.add_at( g.locate(Vec( 0.35, 0.45, 0.15 )), Vec(1,0,0) );
	ASSERT_TRUE( Vec(0.125,1,0).isApprox( u( (Vec( 0.3, 0.4, 0.2 ) ) )) ) ;
}

TEST( geo, tensor )
{

	{
		Vec6 vec ;
		tensor_view( vec ).set_diag( Vec(1,1,1) ) ;
		ASSERT_DOUBLE_EQ( std::sqrt(6.)/2., vec[0] ) ;
		ASSERT_DOUBLE_EQ( 1.5, vec.squaredNorm() ) ;

		const Vec6 vec2 = vec ;
		Mat mat ;
//		tensor_view( vec2 ).set_diag( Vec(1,1,1) ) ; //Compile error
		tensor_view( vec2 ).get( mat ) ;
	}

	{
		DynMat mat(9,9) ;
		mat.setZero() ;
		tensor_view( mat.block<1,6>(0,0) ).set_diag( Vec(1,-1,0) ) ;
		ASSERT_DOUBLE_EQ( 1., mat(0,1) ) ;
		ASSERT_DOUBLE_EQ( 1., mat.squaredNorm() ) ;
	}

	{
		Mat mat, mat2 ;
		mat.setRandom() ;

		Eigen::Matrix< Scalar, 9, 1 > vec;
		tensor_view( vec ).set( mat ) ;
		tensor_view( vec ).get( mat2 ) ;

		ASSERT_TRUE( mat.isApprox(mat2) ) ;
		ASSERT_DOUBLE_EQ( .5 * mat.squaredNorm(), vec.squaredNorm() ) ;

		Mat sym, sym2 ;
		sym = .5 * (mat + mat.transpose() ) ;
		tensor_view( vec.head<6>() ).get( sym2 ) ;
		ASSERT_TRUE( sym.isApprox(sym2) ) ;

		Mat spi = .5 * ( mat - mat.transpose() ) ;
		ASSERT_TRUE( spi.isApprox( tensor_view( vec.segment<3>(6) ).as_mat() ) ) ;
	}


}

TEST( geo, quad ) {

	Voxel vx ;
	vx.origin.setZero() ;
	vx.box = Vec( 2, 3, .5 ) ;

	typename Voxel::QuadPoints  qp ;
	typename Voxel::QuadWeights qp_w ;

	Scalar f_cst  = 0 ;
	Scalar f_lin  = 0 ;
	Scalar f_quad = 0 ;

	vx.get_qp( qp, qp_w );

	for( int k = 0 ; k < Voxel::NQ ; ++k ) {
		f_cst  += qp_w[k] ;
		f_lin  += qp_w[k] * qp.col(k).prod() ;
		f_quad += qp_w[k] * qp.col(k).prod() * qp.col(k).prod() ;
	}

	ASSERT_DOUBLE_EQ(   vx.box.prod()      , f_cst  ) ;
	ASSERT_DOUBLE_EQ( ( vx.box / 2 ).prod(), f_lin  ) ;
	ASSERT_DOUBLE_EQ( ( vx.box / 3 ).prod(), f_quad ) ;

}

TEST( geo, field_func )
{
	Vec3i dim( 10, 10, 10 ) ;
	Vec   box( 1, 1, 1 ) ;
	MeshImpl g( box, dim ) ;

	Linear<MeshImpl> shape(g) ;
	AbstractTensorField<Linear<MeshImpl>> tf( shape ) ;

	Vec6 sym ( Vec6::Zero() ) ;
	sym[0] = 3. ;
	sym[3] = 4. ;

	tf.set_constant( sym );

	AbstractScalarField<Linear<MeshImpl>> sf ( tf.trace() );

	for( Index i = 0 ; i < sf.size() ; ++i ) {
		ASSERT_DOUBLE_EQ( 3., sf[i] ) ;
	}

	sf = tf.norm() ;
	for( Index i = 0 ; i < sf.size() ; ++i ) {
		ASSERT_DOUBLE_EQ( 5., sf[i] ) ;
	}

}

TEST( geo, bc )
{
	{
		StrBoundaryMapper mapper(" left:Slip  right:norMal \t bottom:free ") ;
		EXPECT_EQ( BoundaryInfo::Stick, mapper("top") ) ;
		EXPECT_EQ( BoundaryInfo::Slip, mapper("left") ) ;
		EXPECT_EQ( BoundaryInfo::Normal, mapper("right") ) ;
		EXPECT_EQ( BoundaryInfo::Free, mapper("bottom") ) ;
	}
	{
		StrBoundaryMapper mapper("cuve") ;
		EXPECT_EQ( BoundaryInfo::Normal, mapper("top") ) ;
		EXPECT_EQ( BoundaryInfo::Slip, mapper("left") ) ;
		EXPECT_EQ( BoundaryInfo::Slip, mapper("right") ) ;
		EXPECT_EQ( BoundaryInfo::Slip, mapper("front") ) ;
		EXPECT_EQ( BoundaryInfo::Slip, mapper("back") ) ;
		EXPECT_EQ( BoundaryInfo::Stick, mapper("bottom") ) ;
	}
}

TEST( geo, cross )
{
	const Vec x ( 1, 2, 3 ) ;
	Mat mat ;
	make_cross_mat( x, mat ) ;

	ASSERT_DOUBLE_EQ( (mat*x).norm(), 0. ) ;

	Vec test( -5, 12, 4) ;
	ASSERT_TRUE( (mat * test).isApprox( test.cross(x) ) ) ;
}


TEST( geo, tet )
{
	Tet tet ;
	tet.box = Vec(.5, 1, 1.5 ) ;
	tet.origin = Vec(1.,1.,1.) ;

	Vec point( 1.3, 1.2, 1.4 ) ;
	for( tet.orientation = 0 ; tet.orientation < 16 ; ++tet.orientation ) {
		Tet::Coords c ;
		tet.compute_coords( point, c );
		ASSERT_DOUBLE_EQ(c.sum(), 1) ;
		ASSERT_TRUE( point.isApprox( tet.pos( c ) )) ;
	}

}

TEST( geo, tetGrid )
{
	Vec3i dim( 1, 1, 1 ) ;
	Vec   box( .75, 1.5, 2 ) ;

	TetGrid g( box, dim ) ;

	ASSERT_EQ( 6, g.nCells() ) ;
	ASSERT_EQ( 8, g.nNodes() ) ;

	Tet tet ;
	for( TetGrid::CellIterator it = g.cellBegin() ; it != g.cellEnd() ; ++it )
	{
		//std::cout << (*it).transpose() << std::endl ;
		g.get_geo( *it, tet );

		//std::cout << "origin " << tet.origin.transpose() << std::endl ;
		for( int k = 0 ; k < 4 ; ++ k ) {
			//std::cout << tet.vertex(k).transpose() << std::endl ;
			ASSERT_TRUE( tet.vertex(k).minCoeff() >= 0 ) ;
			ASSERT_TRUE( ( box - tet.vertex(k)).minCoeff() >= 0 ) ;
		}
		ASSERT_TRUE( Vec(0,0,box[2]).isApprox( tet.vertex(3) ) ) ;
	}

	TetGrid::Location loc ;

	for( int i = 0 ; i < 10 ; ++i ) {
		for( int j = 0 ; j < 10 ; ++j ) {
			for( int k = 0 ; k < 10 ; ++k ) {
				Vec p ( box[0]/10*i, box[1]/10*j, box[2]/10*k ) ;
				g.locate( p, loc ) ;
				ASSERT_TRUE( loc.coords.minCoeff() >= -1.e-12 ) ;
				ASSERT_DOUBLE_EQ( 1, loc.coords.sum() ) ;
				g.get_geo(loc.cell, tet) ;
				ASSERT_TRUE( p.isApprox( tet.pos( loc.coords ) ) ) ;
			}
		}
	}

	Linear<TetGrid> shape (g)  ;
	AbstractScalarField< Linear<TetGrid> > field( shape ) ;
	field.set_zero() ;
	field[7] = 1. ;

	ASSERT_DOUBLE_EQ( 1,   field( (box) ) ) ;
	ASSERT_DOUBLE_EQ( 0,   field( (Vec(box[0],0,box[2]) )) ) ;
	ASSERT_DOUBLE_EQ( 0.5, field( (Vec(box[0],0.5*box[1],box[2]) )) ) ;


	for( unsigned k = 0 ; k < 6 ; ++k ) {
		TetGrid::Cell cell = *g.cellBegin() ;
		cell[3] = k ;

		g.get_geo( cell, tet );
		const Vec x1 = tet.center() + 1.e-2 * Vec(1,-2,3) ;

		Linear<TetGrid>::Derivatives dc_dx ;
		Linear<TetGrid>::Interpolation itp ;
		g.locate( x1, loc );
		shape.interpolate( loc, itp ) ;
		shape.get_derivatives( loc, dc_dx );

		for( unsigned i = 0 ; i < 3 ; ++i ) {
			Vec dx ( Vec::Zero() ) ;
			dx[i] = 1.e-2 ;

			Linear<TetGrid>::CoefList coeffs_pred = itp.coeffs + ( dc_dx * dx ) ;

			TetGrid::Location loc2 ;
			Linear<TetGrid>::Interpolation itp2 ;

			const Vec x2 = x1 + dx ;
			g.locate( x2, loc2 );
			shape.interpolate( loc2, itp2 ) ;

			ASSERT_TRUE( coeffs_pred.isApprox( itp2.coeffs ) ) ;

		}
	}

}

TEST( geo, adjacency )
{
	Vec3i dim( 10, 15, 9 ) ;
	Vec   box( 1, 1, 1 ) ;

	MeshImpl g( box, dim ) ;
	Linear<MeshImpl> shape(g) ;
	Eigen::VectorXi nAdj ( g.nNodes() ) ;
	nAdj.setZero() ;

	for( typename MeshType::CellIterator it = g.cellBegin() ; it != g.cellEnd() ; ++it ) {

		typename Linear<MeshImpl>::NodeList nodes ;
		typename MeshImpl::Location loc ; loc.cell = *it ;
		shape.list_nodes( loc, nodes );

		for( Index k = 0 ; k < Linear<MeshImpl>::NI ; ++k ) {
			nAdj[nodes[k]]++ ;
		}
	}

	StrBoundaryMapper mapper("diri") ;
	std::vector< BoundaryInfo > bc ( nAdj.rows() ) ;
	g.make_bc( mapper, bc ) ;

	for( Index i = 0 ; i < nAdj.rows() ; ++i ) {
		if( bc[i].bc == BoundaryInfo::Interior ) {
			ASSERT_EQ( nAdj[i], g.nAdjacent( i ) ) ;
		}
	}

}

TEST( geo, aniso_matrix )
{
	Mat orir ;
	orir.setRandom() ;
	Mat ori = .5 * ( orir + orir.transpose() ) ;

	Mat66 Abar ;
	compute_anisotropy_matrix( ori, Abar );

	ASSERT_TRUE( Abar.isApprox( Abar.transpose() ) ) ;

	Mat taur ;
	taur.setRandom() ;
	Mat tau = .5 * ( taur + taur.transpose() );

	Mat dev = tau -  1./3 * tau.trace() * Mat::Identity() ;
	Mat ani_tau = ori * dev * ori ;
	Mat ani_tau_dev = ani_tau -  1./3 * ani_tau.trace() * Mat::Identity() ;

	ASSERT_NEAR( 0, dev.trace(), 1.e-12 ) ;
	ASSERT_NEAR( 0, ani_tau_dev.trace(), 1.e-12 ) ;

	Vec6 taubar, Abartau, tauref ;
	tensor_view( taubar ).set( tau ) ;
	tensor_view( tauref ).set( ani_tau_dev ) ;
	tauref[0] = taubar[0] ;

	Abartau = Abar * taubar ;

	ASSERT_DOUBLE_EQ( Abartau[0], taubar[0] ) ;
	ASSERT_TRUE( tauref.isApprox( Abartau )) ;
}
