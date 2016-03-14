#include "Tensor.hh"

#include <unsupported/Eigen/MatrixFunctions>

namespace d6
{

void compute_anisotropy_matrix( const Mat& A, MatS & Abar )
{

	Abar.setIdentity() ;

	VecS taubar ;
	Mat tau, ani_tau ;

	for( int k = 1 ; k < SD ; ++k ) {
		taubar.setZero() ;
		taubar[k] = 1. ;

		tensor_view( taubar ).get( tau ) ;
		ani_tau = A * tau * A ;

		tensor_view( Abar.col(k) ).set( ani_tau );
		Abar(0,k) = 0 ;
	}

}

void compute_convection_matrix( const Mat& A, const Scalar dt, MatS & Aexp )
{
	MatS Abar ;

	VecS taubar ;
	Mat tau, conv_tau ;

	for( int k = 0 ; k < SD ; ++k ) {
		taubar.setZero() ;
		taubar[k] = 1. ;

		tensor_view( taubar ).get( tau ) ;
		conv_tau = A * tau + tau * A.transpose() ;

		tensor_view( Abar.col(k) ).set( conv_tau );
	}

	Aexp = (dt * Abar).exp() ;
}

} //d6
