#ifndef D6_TENSOR_HH
#define D6_TENSOR_HH

#include "utils/alg.hh"

namespace d6 {

static constexpr Scalar s_sqrt_2 = std::sqrt(2.) ;
static constexpr Scalar s_sqrt_3 = std::sqrt(3.) ;
static constexpr Scalar s_sqrt_6 = std::sqrt(6.) ;

template< typename Derived >
struct TensorView
{
	static constexpr Index Size = Derived::SizeAtCompileTime ;
	static constexpr bool HasSymPart = ( Size == 6 || Size == 9 ) ;
	static constexpr bool HasSpiPart = ( Size == 3 || Size == 9 ) ;
	static constexpr Index SpiOff = HasSymPart ? 6 : 0 ;

	explicit TensorView( Derived map )
		: m_map( map )
	{
		EIGEN_STATIC_ASSERT_VECTOR_ONLY(Derived);
	}

	TensorView& set_diag( const Vec& diag ) {
		if ( HasSymPart ) {
			m_map.setZero() ;

			m_map[0] = diag.sum() / s_sqrt_6 ;
			m_map[1] = ( diag[0] - diag[1] ) / 2. ;
			m_map[2] = ( 2*diag[2] - diag[0] - diag[1] ) / (2 * s_sqrt_3) ;
		}
		return *this ;
	}

	TensorView& set( const Mat& mat ) {

		if ( HasSymPart ) {
			m_map[0] = ( mat(0,0)+mat(1,1)+mat(2,2) ) / s_sqrt_6 ;
			m_map[1] = ( mat(0,0) - mat(1,1) ) / 2. ;
			m_map[2] = ( 2*mat(2,2) - mat(1,1) - mat(0,0) ) / (2 * s_sqrt_3) ;

			m_map[3] = .5*( mat(0,1) + mat(1,0) ) ;
			m_map[4] = .5*( mat(0,2) + mat(2,0) ) ;
			m_map[5] = .5*( mat(1,2) + mat(2,1) ) ;
		}

		if( HasSpiPart ) {
			m_map[SpiOff+0] = .5*( mat(0,1) - mat(1,0) ) ;
			m_map[SpiOff+1] = .5*( mat(0,2) - mat(2,0) ) ;
			m_map[SpiOff+2] = .5*( mat(1,2) - mat(2,1) ) ;
		}
		return *this ;
	}

	void get( Mat& mat ) const {
		if ( HasSymPart ) {
			mat(0,0) = 2./s_sqrt_6 * m_map[0] + m_map[1] - 1./s_sqrt_3 * m_map[2] ;
			mat(1,1) = 2./s_sqrt_6 * m_map[0] - m_map[1] - 1./s_sqrt_3 * m_map[2] ;
			mat(2,2) = 2./s_sqrt_6 * m_map[0]            + 2./s_sqrt_3 * m_map[2] ;
			mat(0,1) = mat(1,0) = m_map[3] ;
			mat(0,2) = mat(2,0) = m_map[4] ;
			mat(1,2) = mat(2,1) = m_map[5] ;
		} else mat.setZero() ;

		if( HasSpiPart ) {
			mat(0,1) += m_map[SpiOff+0] ;
			mat(1,0) -= m_map[SpiOff+0] ;
			mat(0,2) += m_map[SpiOff+1] ;
			mat(2,0) -= m_map[SpiOff+1] ;
			mat(1,2) += m_map[SpiOff+2] ;
			mat(2,1) -= m_map[SpiOff+2] ;
		}
	}

	Mat as_mat() const {
		Mat mat ;
		get( mat ) ;
		return mat ;
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


} //d6

#endif