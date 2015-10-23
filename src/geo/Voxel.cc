#include "Voxel.hh"

#include "Tensor.hh"

namespace d6 {

Index Voxel::sample_uniform(const unsigned N, const Index start, Points &points, Frames &frames) const
{
	Scalar min = box.minCoeff() ;

	Vec3i Nsub ;
	for( int k = 0 ; k < 3 ; ++ k)
		Nsub[k] = N * std::round( box[k] / min ) ;

	const Vec subBox = box.array() / Nsub.array().cast< Scalar >() ;

	Vec6 frame ;
	tensor_view( frame ).set_diag( Vec( .25 * subBox.array() * subBox.array() ) ) ;

	Index p = start ;
	for( int i = 0 ; i < Nsub[0] ; ++i )
		for( int j = 0 ; j < Nsub[1] ; ++j )
			for( int k = 0 ; k < Nsub[2] ; ++k ) {
				points.col(p) = corner + (Vec(i+.5,j+.5,k+.5).array() * subBox.array()).matrix() ;
				frames.col(p) = frame ;
				++p ;
			}

	return p - start ;
}

}
