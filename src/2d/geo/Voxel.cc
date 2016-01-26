#include "Voxel.hh"

#include "Tensor.hh"

namespace d6 {

Index Voxel::sample_uniform(const unsigned N, const Index start, Points &points, Frames &frames) const
{
	Scalar min = box.minCoeff() ;

	VecWi Nsub ;
	for( int k = 0 ; k < WD ; ++ k)
		Nsub[k] = N * std::round( box[k] / min ) ;

	const Vec subBox = box.array() / Nsub.array().cast< Scalar >() ;

	VecS frame ;
	tensor_view( frame ).set_diag( Vec( .25 * subBox.array() * subBox.array() ) ) ;

	Index p = start ;
	for( int i = 0 ; i < Nsub[0] ; ++i )
		for( int j = 0 ; j < Nsub[1] ; ++j ) {
				points.col(p) = origin + (Vec(i+.5,j+.5).array() * subBox.array()).matrix() ;
				frames.col(p) = frame ;
				++p ;
			}

	return p - start ;
}

Voxel::QuadPoints Voxel::Qps()
{
	// .5 * ( 1 +- 1./sqrt(3) )
	const Vec dqp = Vec::Constant( 1./sqrt(3.) );
	//		const Vec dqp = Vec::Constant( 1. );
	const Vec qp0 = .5 * ( Vec::Ones() - dqp );

	QuadPoints qps ;
	for( int i = 0 ; i < 2 ; ++i ) {
		for( int j = 0 ; j < 2 ; ++j ) {
				VecWi corner ( i, j ) ;
				qps.col( cornerIndex(corner) ) = qp0.array() + corner.cast< Scalar >().array()*dqp.array() ;
		}
	}
	return qps ;
}

} //d6
