#include "FormBuilder.hh"

#include "geo/MeshImpl.hh"
#include "geo/Tensor.hh"

#include <bogus/Core/Block.impl.hpp>

namespace d6 {


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

void FormBuilder::addDuDv( FormMat<2,2>::Type& A, Scalar w,
						   Index rowIndex, const typename MeshType::Derivatives::ConstRowXpr& row_dx,
						   Itp itp, Dcdx dc_dx, Indices colIndices )
{
	typedef FormMat<2,2>::Type::BlockType Block ;

	for( int j = 0 ; j < MeshType::NV ; ++j ) {
			Block &b = A.block( rowIndex, colIndices[itp.nodes[j]] ) ;

			// dux_dx, duy_dy, duz_dz
			b(0,0) +=      w * row_dx(0) * dc_dx(j, 0)
				   +  .5 * w * row_dx(1) * dc_dx(j, 1) ;
			b(1,1) +=      w * row_dx(1) * dc_dx(j, 1)
				   +  .5 * w * row_dx(0) * dc_dx(j, 0) ;

			// 2 * .25 * (dux_dy + duy_dx) * (dvx_dy + dvy_dx)
			b(0,1) += .5 * w * row_dx(1) * dc_dx(j, 0) ;
			b(1,0) += .5 * w * row_dx(0) * dc_dx(j, 1) ;
	}
}

void FormBuilder::addTauDu( FormMat<3,2>::Type& A, Scalar m, Index rowIndex, Itp itp, Dcdx dc_dx, Indices colIndices )
{
	typedef FormMat<3,2>::Type::BlockType Block ;

	for( int j = 0 ; j < MeshType::NV ; ++j ) {
		Block &b = A.block( rowIndex, colIndices[itp.nodes[j]] ) ;

		// a * (dux_dx + duy_dy )
		b(0,0) += m * dc_dx(j, 0) ;
		b(0,1) += m * dc_dx(j, 1) ;

		// b * (dux_dx - duy_dy )
		b(1,0) += m * dc_dx(j, 0) ;
		b(1,1) -= m * dc_dx(j, 1) ;

		// c * ( dux_dy + duy_dx )
		b(2,0) += m * dc_dx(j, 1) ;
		b(2,1) += m * dc_dx(j, 0) ;
	}
}

void FormBuilder::addVDp ( FormMat<2,1>::Type& A, Scalar m, Index rowIndex, Itp itp, Dcdx dc_dx, Indices colIndices )
{
	typedef FormMat<2,1>::Type::BlockType Block ;

	for( int j = 0 ; j < MeshType::NV ; ++j ) {
		Block &b = A.block( rowIndex, colIndices[itp.nodes[j]] ) ;

		// ( vx da_dx + vy da_dy + vz da_dz)
		b += m * dc_dx.row(j) ;
	}
}

void FormBuilder::addTauWu( FormMat<1,2>::Type& A, Scalar m, Index rowIndex, Itp itp, Dcdx dc_dx, Indices colIndices )
{
	typedef FormMat<1,2>::Type::BlockType Block ;

	for( int j = 0 ; j < MeshType::NV ; ++j ) {
		Block &b = A.block( rowIndex, colIndices[itp.nodes[j]] ) ;

		// d * ( dux_dy - duy_dx )
		b(0,0) += m * dc_dx(j, 1) ;
		b(0,1) -= m * dc_dx(j, 0) ;
	}
}

void FormBuilder::addUTauGphi( FormMat<3,2>::Type& A, Scalar w, Itp itp, const Vec& dphi_dx, Indices rowIndices, Indices colIndices )
{
	typedef FormMat<3,2>::Type::BlockType Block ;

//#pragma omp parallel for
	for( int k = 0 ; k < MeshType::NV ; ++k ) {
		for( int j = 0 ; j < MeshType::NV ; ++j ) {
			Block &b = A.block( rowIndices[itp.nodes[k]], colIndices[itp.nodes[j]] ) ;
			const Scalar m = w * itp.coeffs[k] * itp.coeffs[j] ;

			// a * (dphi_dx ux + dphi_dy uy + dphi_dz uz)
			b(0,0) += m * dphi_dx(0) ;
			b(0,1) += m * dphi_dx(1) ;

			// b * (dphi_dx ux - dphi_dy uy )
			b(1,0) += m * dphi_dx(0) ;
			b(1,1) -= m * dphi_dx(1) ;

			// c * ( dux_dy + duy_dx )
			b(2,0) += m * dphi_dx(1) ;
			b(2,1) += m * dphi_dx(0) ;
		}
	}

}

void FormBuilder::addUTaunGphi( FormMat<3,2>::Type& A, Scalar w, Itp itp, const Vec& dphi_dx, Indices rowIndices, Indices colIndices )
{
	typedef FormMat<3,2>::Type::BlockType Block ;

//#pragma omp parallel for
	for( int k = 0 ; k < MeshType::NV ; ++k ) {
		for( int j = 0 ; j < MeshType::NV ; ++j ) {
			Block &b = A.block( rowIndices[itp.nodes[k]], colIndices[itp.nodes[j]] ) ;
			const Scalar m = w * itp.coeffs[k] * itp.coeffs[j] ;

			// a * (dphi_dx ux + dphi_dy uy + dphi_dz uz)
			b.row(0) += m * dphi_dx ;
		}
	}
}

} //d6
