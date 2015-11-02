
#include <gtest/gtest.h>

#include "geo/Grid.hh"

#include "geo/ScalarField.hh"
#include "geo/VectorField.hh"
#include "geo/TensorField.hh"

#include "geo/Tensor.hh"
#include "geo/Voxel.hh"

using namespace d6 ;

TEST( geo, grid )
{
	Vec3i dim( 10, 5, 4 ) ;
	Vec   box( 1, 1, 1 ) ;
	Grid g( box, dim ) ;

	ASSERT_EQ( 200, g.nCells() ) ;
	ASSERT_EQ( 11*6*5, g.nNodes() ) ;
	ASSERT_TRUE( box.isApprox( g.box() ) ) ;

	Grid::Interpolation itp ;
	g.interpolate( Vec(0,0,0), itp );

	ASSERT_EQ( 0, itp.nodes[0] ) ;
	ASSERT_EQ( dim[2]+1, itp.nodes[2] ) ;
	ASSERT_EQ( (dim[1]+1)*(dim[2]+1)+(dim[2]+1)+1, itp.nodes[7] ) ;
	ASSERT_DOUBLE_EQ( 1, itp.coeffs[0] ) ;

	g.interpolate( Vec(1,1,1), itp );

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

	Grid::Derivatives dc_dx ;
	Grid::Location loc ;
	g.locate( x1, loc );
	g.interpolate( loc, itp ) ;
	g.get_derivatives( loc, dc_dx );

	for( unsigned i = 0 ; i < 3 ; ++i ) {
		Vec dx ( Vec::Zero() ) ;
		dx[i] = 1.e-2 ;

		Grid::CoefList coeffs_pred = itp.coeffs + ( dc_dx * dx ) ;

		const Vec x2 = x1 + dx ;
		Grid::Interpolation itp2 ;
		g.interpolate( x2, itp2 ) ;

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

	AbstractScalarField< Grid > phi( g ) ;
	ASSERT_EQ( g.nNodes(), phi.flatten().rows() ) ;

	phi.set_constant( 3 );
	ASSERT_DOUBLE_EQ( 3., phi( Vec( 0.2, 0.7, 0.5 ) ) ) ;

	phi.set_zero();
	ASSERT_DOUBLE_EQ( 0., phi( Vec( 0.3, 0.4, 0.1 ) ) ) ;
	phi.add_at( Vec( 0.3, 0.4, 0.1 ), 1 );
	ASSERT_FLOAT_EQ( 1,  phi( Vec( 0.3, 0.4, 0.1 ) ) ) ;
	ASSERT_FLOAT_EQ( 0.3,  phi( Vec( 0.37, 0.4, 0.1 ) ) ) ;
	ASSERT_FLOAT_EQ( 0.7,  phi( Vec( 0.3, 0.37, 0.1 ) ) ) ;
	ASSERT_FLOAT_EQ( 0.5,  phi( Vec( 0.3, 0.4, 0.15 ) ) ) ;

	AbstractVectorField< Grid > u( g ) ;
	u.set_constant( Vec(0,1,0) ) ;
	ASSERT_TRUE( Vec(0,1,0).isApprox( u( Vec( 0.2, 0.7, 0.5 ) ) ) ) ;

	u.add_at( Vec( 0.35, 0.45, 0.15 ), Vec(1,0,0) );
	ASSERT_TRUE( Vec(0.125,1,0).isApprox( u( Vec( 0.3, 0.4, 0.2 ) ) ) ) ;
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
	Grid g( box, dim ) ;

	TensorField tf( g ) ;

	Vec6 sym ( Vec6::Zero() ) ;
	sym[0] = 3. ;
	sym[3] = 4. ;

	tf.set_constant( sym );

	ScalarField sf ( tf.trace() );

	for( Index i = 0 ; i < sf.size() ; ++i ) {
		ASSERT_DOUBLE_EQ( 3., sf[i] ) ;
	}

	sf = tf.norm() ;
	for( Index i = 0 ; i < sf.size() ; ++i ) {
		ASSERT_DOUBLE_EQ( 5., sf[i] ) ;
	}

}
