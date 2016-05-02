#ifndef D6_SERIALIZATION_HH
#define D6_SERIALIZATION_HH

#include <boost/serialization/array.hpp>
#include <boost/serialization/split_free.hpp>
#include <bogus/Core/Eigen/EigenSerialization.hpp>

namespace boost
{
namespace serialization
{

template<typename Archive, typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols>
inline void load(
	   Archive & ar,
	   Eigen::Array<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols> & matrix,
	   const unsigned int file_version
   )
{
	static_assert( _Rows != Eigen::Dynamic, "Dynamic array serial not implemented" ) ;
	static_assert( _Cols != Eigen::Dynamic, "Dynamic array serial not implemented" ) ;
	(void) file_version ;
	ar & make_array( matrix.data(), matrix.size() ) ;
}

template<typename Archive, typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols>
inline void save(
	   Archive & ar,
	   const Eigen::Array<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols> & matrix,
	   const unsigned int file_version
   )
{
	static_assert( _Rows != Eigen::Dynamic, "Dynamic array serial not implemented" ) ;
	static_assert( _Cols != Eigen::Dynamic, "Dynamic array serial not implemented" ) ;
	(void) file_version ;
	ar & make_array( matrix.data(), matrix.size() ) ;
}

template<typename Archive, typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols>
inline void serialize(
	   Archive & ar,
	   Eigen::Array<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols> & matrix,
	   const unsigned int file_version
   )
{
	split_free( ar, matrix, file_version ) ;
}

} //serialization
} //boost

#endif
