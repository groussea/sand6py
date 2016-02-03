
#include "geo/Tensor.hh"

namespace d6 {

struct FormBlocks {

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

    template < typename RowDerived, typename ColDerived >
    static inline void addDuDv( Eigen::Matrix<Scalar, WD, WD>& b,
                                const Scalar w,
                                const Eigen::MatrixBase< RowDerived >& row_dx,
                                const Eigen::MatrixBase< ColDerived >& col_dx )
    {

        // dux_dx, duy_dy, duz_dz
        b(0,0) +=      w * row_dx(0) * col_dx(0)
                + .5 * w * row_dx(1) * col_dx(1) ;
        b(1,1) +=      w * row_dx(1) * col_dx(1)
                + .5 * w * row_dx(0) * col_dx(0) ;

        // 2 * .25 * (dux_dy + duy_dx) * (dvx_dy + dvy_dx)
        b(0,1) += .5 * w * row_dx(1) * col_dx(0) ;
        b(1,0) += .5 * w * row_dx(0) * col_dx(1) ;
    }

    template < typename ColDerived >
    static inline void addTauDu( Eigen::Matrix<Scalar, SD, WD>& b,
                                 const Scalar m,
                                 const Eigen::MatrixBase< ColDerived >& col_dx )
    {
        // a * (dux_dx + duy_dy )
        b(0,0) += m * col_dx(0) ;
        b(0,1) += m * col_dx(1) ;

        // b * (dux_dx - duy_dy )
        b(1,0) += m * col_dx(0) ;
        b(1,1) -= m * col_dx(1) ;

        // c * ( dux_dy + duy_dx )
        b(2,0) += m * col_dx(1) ;
        b(2,1) += m * col_dx(0) ;
    }

    template < typename ColDerived >
    static inline void addVDp( Eigen::Matrix<Scalar, WD, 1>& b,
                               const Scalar m,
                               const Eigen::MatrixBase< ColDerived >& col_dx )
    {
        // ( vx da_dx + vy da_dy + vz da_dz)
        b += m * col_dx ;
    }

    template < typename ColDerived >
    static inline void addTauWu( Eigen::Matrix<Scalar, RD, WD>& b,
                                 const Scalar m,
                                 const Eigen::MatrixBase< ColDerived >& col_dx )
    {
        // d * ( dux_dy - duy_dx )
        b(0,0) += m * col_dx(1) ;
        b(0,1) -= m * col_dx(0) ;
    }

    static inline void addUTauGphi( Eigen::Matrix<Scalar, SD, WD>& b,
                                    const Scalar m,
                                    const Vec& dphi_dx )
    {
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

    static inline void addUTaunGphi( Eigen::Matrix<Scalar, SD, WD>& b,
                                     const Scalar m,
                                     const Vec& dphi_dx )
    {
        // a * (dphi_dx ux + dphi_dy uy + dphi_dz uz)
        b.row(0) += m * dphi_dx ;
    }

} ;

} //d6
