/*
 * This file is part of Sand6, a C++ continuum-based granular simulator.
 *
 * Copyright 2016 Gilles Daviet <gilles.daviet@inria.fr> (Inria - Universit√© Grenoble Alpes)
 *
 * Sand6 is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.

 * Sand6 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.

 * You should have received a copy of the GNU General Public License
 * along with Sand6.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef D6_TENSOR_HH
#define D6_TENSOR_HH

#include "utils/alg.hh"

namespace d6 {

static constexpr Scalar s_sqrt_2  = M_SQRT2 ;
static constexpr Scalar s_sqrt_3  = 1.73205080756887729352744634151 ;
static constexpr Scalar s_sqrt_2_d = 1 ;

template< typename Derived >
struct TensorView
{
    static constexpr Index Size = Derived::SizeAtCompileTime ;
    static constexpr bool HasSymPart = ( Size == 3 || Size == 4 ) ;
    static constexpr bool HasSpiPart = ( Size == 1 || Size == 4 ) ;
    static constexpr Index SpiOff = HasSymPart ? 3 : 0 ;

    explicit TensorView( Derived map )
        : m_map( map )
    {
        EIGEN_STATIC_ASSERT_VECTOR_ONLY(Derived);
    }

    TensorView& set_diag( const Vec& diag ) {
        if ( HasSymPart ) {
            m_map.setZero() ;

            m_map[0] = diag.sum() / 2 ;
            m_map[1] = ( diag[0] - diag[1] ) / 2. ;
        }
        return *this ;
    }

    TensorView& set( const Mat& mat ) {

        if ( HasSymPart ) {
            m_map[0] = ( mat(0,0) + mat(1,1) ) / 2 ;
            m_map[1] = ( mat(0,0) - mat(1,1) ) / 2. ;

            m_map[2] = .5*( mat(0,1) + mat(1,0) ) ;
        }

        if( HasSpiPart ) {
            m_map[SpiOff+0] = .5*( mat(0,1) - mat(1,0) ) ;
        }
        return *this ;
    }

    void get( Mat& mat ) const {
        mat.setZero() ;
        add( mat ) ;
    }
    void add( Mat& mat ) const {

        if ( HasSymPart ) {
            mat(0,0) += m_map[0] + m_map[1] ;
            mat(1,1) += m_map[0] - m_map[1] ;
            mat(0,1) += m_map[2] ;
            mat(1,0) += m_map[2] ;
        }

        if( HasSpiPart ) {
            mat(0,1) += m_map[SpiOff+0] ;
            mat(1,0) -= m_map[SpiOff+0] ;
        }
    }

    Mat as_mat() const {
        Mat mat ;
        get( mat ) ;
        return mat ;
    }

    const Derived& vec() const {
        return m_map ;
    }
    Derived& vec() {
        return m_map ;
    }

private:

    Derived m_map ;

};

template< typename Derived >
TensorView< Derived > tensor_view( const Eigen::MapBase< Derived, Eigen::ReadOnlyAccessors >& map )
{
    return TensorView< Derived > ( map.derived() ) ;
}

template< typename Derived >
TensorView< Derived > tensor_view( const Eigen::MapBase< Derived, Eigen::WriteAccessors >& map )
{
    return TensorView< Derived > ( map.derived() ) ;
}

template< typename Derived >
TensorView< typename Derived::MapType > tensor_view( Eigen::MatrixBase< Derived >& vec )
{
    return tensor_view( Derived::Map(vec.derived().data(), vec.size()) ) ;
}

template< typename Derived >
TensorView< typename Derived::ConstMapType > tensor_view( const Eigen::MatrixBase< Derived >& vec )
{
    return TensorView< typename Derived::ConstMapType >( Derived::Map(vec.derived().data(), vec.size()) ) ;
}

inline TensorView< typename Eigen::Matrix<Scalar,1,1>::MapType > tensor_view( Scalar &s ) {
    return tensor_view( Eigen::Matrix<Scalar,1,1>::Map( &s ) ) ;
}
inline TensorView< typename Eigen::Matrix<Scalar,1,1>::ConstMapType > tensor_view( const Scalar &s ) {
    return TensorView< typename Eigen::Matrix<Scalar,1,1>::ConstMapType >( Eigen::Matrix<Scalar,1,1>::Map( &s ) ) ;
}

// mat * a == a ^ x
template <typename Derived>
void make_cross_mat( const Vec& x, const Eigen::MapBase< Derived, Eigen::WriteAccessors >& map )
{
    Derived mat = map.derived() ;
    mat(0) = -x[1] ;
    mat(1) =  x[0] ;
}
inline void make_cross_mat( const Vec& x, Vec& mat ) {
    return make_cross_mat( x, mat.block<WD,RD>(0,0) ) ;
}

// such that Abar bar(tau) = (tauN ; bar( A Dev(tau) A ) )
void compute_anisotropy_matrix( const Mat& A, MatS & Abar ) ;

// exp(dt * Abar) such that Abar bar(tau) = A tau + tau A'
void compute_convection_matrix( const Mat& A, const Scalar dt, MatS & Aexp ) ;

template <typename Derived>
void convect( const Scalar dt, const Mat&A, TensorView< Derived > tau ) {
    MatS Aexp ;
    compute_convection_matrix( A, dt, Aexp );
    tau.vec() = Aexp * tau.vec() ;
}

} //d6

#endif
