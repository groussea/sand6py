#include "FormBuilder.hh"

#include "geo/MeshImpl.hh"
#include "geo/Tensor.hh"

#include <bogus/Core/Block.impl.hpp>

namespace d6 {


// FIXME differente approx for

/*
 *
 * D( u )
 * a = Tr D(u) / sqrt(6) = ( dux_dx + duy_dy + duz_dz ) / sqrt(6)
 * b = .5 ( dux_dx - duy_dy )
 * c = ( 2*duz_dz - duy_dy - dux_dx ) / (2 * s_sqrt_3)
 * d = .5 (dux_dy + duy_dx)
 * e = .5 (dux_dz + duz_dx)
 * f = .5 (duz_dy + duy_dz)
 */

void FormBuilder::addDuDv( FormMat<3,3>::Type& A, Scalar w,
						   Index rowIndex, const DcdxRow& row_dx,
						   Itp itp, Dcdx dc_dx, Indices colIndices )
{
	typedef FormMat<3,3>::Type::BlockType Block ;

	for( int j = 0 ; j < Linear<MeshImpl>::NI ; ++j ) {
			Block &b = A.block( rowIndex, colIndices[itp.nodes[j]] ) ;

			// dux_dx, duy_dy, duz_dz
			b(0,0) +=      w * row_dx(0) * dc_dx(j, 0)
				   +  .5 * w * row_dx(1) * dc_dx(j, 1)
				   +  .5 * w * row_dx(2) * dc_dx(j, 2) ;
			b(1,1) +=      w * row_dx(1) * dc_dx(j, 1)
				   +  .5 * w * row_dx(0) * dc_dx(j, 0)
				   +  .5 * w * row_dx(2) * dc_dx(j, 2) ;
			b(2,2) +=      w * row_dx(2) * dc_dx(j, 2)
				   +  .5 * w * row_dx(0) * dc_dx(j, 0)
				   +  .5 * w * row_dx(1) * dc_dx(j, 1) ;

			// 2 * .25 * (dux_dy + duy_dx) * (dvx_dy + dvy_dx)
			b(0,1) += .5 * w * row_dx(1) * dc_dx(j, 0) ;
			b(1,0) += .5 * w * row_dx(0) * dc_dx(j, 1) ;

			// 2 * .25 * (dux_dz + duz_dx) * (dvx_dz + dvz_dx)
			b(0,2) += .5 * w * row_dx(2) * dc_dx(j, 0) ;
			b(2,0) += .5 * w * row_dx(0) * dc_dx(j, 2) ;

			// 2 * .25 * (duz_dy + duy_dz) * (dvz_dy + dvy_dz)
			b(2,1) += .5 * w * row_dx(1) * dc_dx(j, 2) ;
			b(1,2) += .5 * w * row_dx(2) * dc_dx(j, 1) ;
	}
}

void FormBuilder::addTauDu( FormMat<6,3>::Type& A, Scalar m, Index rowIndex, Itp itp, Dcdx dc_dx, Indices colIndices )
{
	typedef FormMat<6,3>::Type::BlockType Block ;

	for( int j = 0 ; j < Linear<MeshImpl>::NI ; ++j ) {
		Block &b = A.block( rowIndex, colIndices[itp.nodes[j]] ) ;

		// a * sqrt2_3 * (dux_dx + duy_dy + duz_dz)
		b(0,0) += m * s_sqrt_23 * dc_dx(j, 0) ;
		b(0,1) += m * s_sqrt_23 * dc_dx(j, 1) ;
		b(0,2) += m * s_sqrt_23 * dc_dx(j, 2) ;

		// b * (dux_dx - duy_dy )
		b(1,0) += m * dc_dx(j, 0) ;
		b(1,1) -= m * dc_dx(j, 1) ;

		// c * isqrt_3 * ( -dux_dx - duy_dy + 2*duz_dz)
		b(2,0) -= m * s_isqrt_3 * dc_dx(j, 0) ;
		b(2,1) -= m * s_isqrt_3 * dc_dx(j, 1) ;
		b(2,2) += m * s_isqrt_3 * dc_dx(j, 2) * 2;

		// d * ( dux_dy + duy_dx )
		b(3,0) += m * dc_dx(j, 1) ;
		b(3,1) += m * dc_dx(j, 0) ;

		// e * ( dux_dz + duz_dx )
		b(4,0) += m * dc_dx(j, 2) ;
		b(4,2) += m * dc_dx(j, 0) ;

		// f * ( duz_dy + duy_dz )
		b(5,2) += m * dc_dx(j, 1) ;
		b(5,1) += m * dc_dx(j, 2) ;
	}
}

void FormBuilder::addVDp ( FormMat<3,1>::Type& A, Scalar m, Index rowIndex, Itp itp, Dcdx dc_dx, Indices colIndices )
{
	typedef FormMat<3,1>::Type::BlockType Block ;

	for( int j = 0 ; j < Linear<MeshImpl>::NI ; ++j ) {
		Block &b = A.block( rowIndex, colIndices[itp.nodes[j]] ) ;

		// ( vx da_dx + vy da_dy + vz da_dz)
		b += m * dc_dx.row(j) ;
	}
}

void FormBuilder::addTauWu( FormMat<3,3>::Type& A, Scalar m, Index rowIndex, Itp itp, Dcdx dc_dx, Indices colIndices )
{
	typedef FormMat<3,3>::Type::BlockType Block ;

	for( int j = 0 ; j < Linear<MeshImpl>::NI ; ++j ) {
		Block &b = A.block( rowIndex, colIndices[itp.nodes[j]] ) ;

		// i * ( dux_dy - duy_dx )
		b(0,0) += m * dc_dx(j, 1) ;
		b(0,1) -= m * dc_dx(j, 0) ;

		// j * ( dux_dz - duz_dx )
		b(1,0) += m * dc_dx(j, 2) ;
		b(1,2) -= m * dc_dx(j, 0) ;

		// k * ( duy_dz - duz_dy )
		b(2,1) += m * dc_dx(j, 2) ;
		b(2,2) -= m * dc_dx(j, 1) ;
	}
}

void FormBuilder::addUTauGphi( FormMat<6,3>::Type& A, Scalar w, Itp itp, const Vec& dphi_dx, Indices rowIndices, Indices colIndices )
{
	typedef FormMat<6,3>::Type::BlockType Block ;

//#pragma omp parallel for
	for( int k = 0 ; k < Linear<MeshImpl>::NI ; ++k ) {
		for( int j = 0 ; j < Linear<MeshImpl>::NI ; ++j ) {
			Block &b = A.block( rowIndices[itp.nodes[k]], colIndices[itp.nodes[j]] ) ;
			const Scalar m = w * itp.coeffs[k] * itp.coeffs[j] ;

			// a * sqrt2_3 * (dphi_dx ux + dphi_dy uy + dphi_dz uz)
			b(0,0) += m * s_sqrt_23 * dphi_dx(0) ;
			b(0,1) += m * s_sqrt_23 * dphi_dx(1) ;
			b(0,2) += m * s_sqrt_23 * dphi_dx(2) ;

			// b * (dphi_dx ux - dphi_dy uy )
			b(1,0) += m * dphi_dx(0) ;
			b(1,1) -= m * dphi_dx(1) ;

			// c * isqrt_3 * ( -dux_dx - duy_dy + 2*duz_dz)
			b(2,0) -= m * s_isqrt_3 * dphi_dx(0) ;
			b(2,1) -= m * s_isqrt_3 * dphi_dx(1) ;
			b(2,2) += m * s_isqrt_3 * dphi_dx(2) * 2;

			// d * ( dux_dy + duy_dx )
			b(3,0) += m * dphi_dx(1) ;
			b(3,1) += m * dphi_dx(0) ;

			// e * ( dux_dz + duz_dx )
			b(4,0) += m * dphi_dx(2) ;
			b(4,2) += m * dphi_dx(0) ;

			// f * ( duz_dy + duy_dz )
			b(5,2) += m * dphi_dx(1) ;
			b(5,1) += m * dphi_dx(2) ;
		}
	}

}

void FormBuilder::addUTaunGphi( FormMat<6,3>::Type& A, Scalar w, Itp itp, const Vec& dphi_dx, Indices rowIndices, Indices colIndices )
{
	typedef FormMat<6,3>::Type::BlockType Block ;

//#pragma omp parallel for
	for( int k = 0 ; k < Linear<MeshImpl>::NI ; ++k ) {
		for( int j = 0 ; j < Linear<MeshImpl>::NI ; ++j ) {
			Block &b = A.block( rowIndices[itp.nodes[k]], colIndices[itp.nodes[j]] ) ;
			const Scalar m = w * itp.coeffs[k] * itp.coeffs[j] ;

			// a * sqrt2_3 * (dphi_dx ux + dphi_dy uy + dphi_dz uz)
			b.row(0) += m * s_sqrt_23 * dphi_dx ;
		}
	}
}

} //d6
