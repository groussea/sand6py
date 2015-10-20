#ifndef D6_STRING_HH
#define D6_STRING_HH

#include <string>
#include <vector>
#include <sstream>

#include "alg.hh"

namespace d6 {

//! Left- and Right-trim \p res, removing all the characters in \p chars
std::string trim ( const std::string& str, const std::string& chars ) ;
//! Same as trim( str, all_white_space_characters )
std::string trim ( const std::string& str ) ;

std::string to_lower ( const std::string& str ) ;
std::string to_upper ( const std::string& str ) ;

float       to_float  ( const std::string& str ) ;
double      to_double ( const std::string& str ) ;
int         to_int    ( const std::string& str ) ;
unsigned    to_uint   ( const std::string& str ) ;
std::size_t to_size_t ( const std::string& str ) ;
bool        to_bool   ( const std::string& str ) ;


//! Trim a string, replace all whitespace sequences with a single space, and puts the result in lowercase
std::string canonify( const std::string& str ) ;
//! Same thing, with less frenglish
std::string canonicalize( const std::string& str ) ;

//! Explodes \p src into an array \p res, splitting on any character that is present in \p separators
void split ( const std::string& str, const std::string &separators, std::vector< std::string > &res ) ;
//! Same as split( src, all_white_space_characters )
void split ( const std::string& str, std::vector< std::string > &res ) ;

//! \return a string which is the concatenation of all elements of \p src separated by \p separator
std::string join ( const std::vector< std::string > &src, const std::string &separator ) ;
//! Same as join( src, ' ' )
std::string join ( const std::vector< std::string > &src ) ;

//! \return whether src starts with \p start
bool starts_with( const std::string &src, const std::string &start ) ;
//! \return whether src ends with \p end
bool ends_with( const std::string &src, const std::string &end ) ;

// Poor-man version of QString::arg
bool split_on_next_marker ( const std::string &src, std::string &first_part, std::string &second_part ) ;

//!
template< typename T >
std::string arg( const std::string &src, const T &replacement )
{
	std::string p1, p2 ;
	if( split_on_next_marker( src, p1, p2 ) )
	{
		std::stringstream os( std::stringstream::out ) ;
		os << p1 << replacement << p2 ;
		return os.str() ;
	}
	return src ;
}
//! Same with two arguments
template< typename T1, typename T2 >
std::string arg( const std::string &src, const T1 &r1, const T2& r2 )
{
  return arg< T2 >(arg< T1 >( src, r1 ), r2 ) ;
}

//! Same with three arguments
template< typename T1, typename T2, typename T3 >
std::string arg3( const std::string &src, const T1 &r1, const T2& r2, const T3& r3 )
{
  return arg< T3 >(arg< T1 >( src, r1, r2 ), r3 ) ;
}

// Templated comversion functions using operator>>

template< typename T >
bool cast ( std::istringstream& stream, T& res )
{
	return ( stream >> res ) ;
}

template< typename Scalar, int Rows, int Cols >
bool cast ( std::istringstream& stream, Eigen::Matrix<Scalar, Rows, Cols>& res )
{
	for( size_t k = 0 ; stream && k < res.size() ; ++ k) {
		stream >> res.data()[k] ;
	}
	return stream ;
}


template< typename T >
bool cast ( const std::string& str, T& res )
{
	std::istringstream stream( str ) ;
	return cast( stream, res ) ;
}

template< typename NumType >
NumType to_num ( const std::string& str )
{
	NumType n ;
	return cast( str, n ) ? n : 0 ;
}

} //namespace d6


#define D6_stringify(s) D6_preproc_str(s)
#define D6_preproc_str(s) #s

#endif
