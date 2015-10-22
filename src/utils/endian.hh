#ifndef D6_ENDIAN_HH
#define D6_ENDIAN_HH

#include <cstdint>

namespace d6 {

// Reading/Writing  Big Endian streams

template < size_t >
struct byte_order {} ;

template <>
struct byte_order< 4 >
{
	typedef union {
		std::uint32_t i;
		unsigned char c[4];
	} word ;

	static const unsigned char* order()
	{
		static const word s_order = {0x00010203};
		return s_order.c ;
	}
} ;

template <>
struct byte_order< 8 >
{
	typedef union {
		std::uint64_t i;
		unsigned char c[8];
	} word ;

	static const unsigned  char* order()
	{
		static const word s_order = {0x0001020304050607};
		return s_order.c ;
	}
} ;

template < typename T >
T from_big_endian( const unsigned char* data )
{
	typedef byte_order< sizeof( T ) > BO ;
	union {
		T i;
		unsigned char c[8];
	} reordered ;

	for( unsigned i = 0 ; i < sizeof( T ) ; ++i )
	{
		reordered.c[ i ] = data[ BO::order()[i] ] ;
	}

	return  reordered.i ;
}

template < typename T >
void to_big_endian( const T& num, unsigned char* data )
{
	typedef byte_order< sizeof( T ) > BO ;
	union {
		T i;
		unsigned char c[8];
	} reordered ;

	reordered.i = num ;
	for( unsigned i = 0 ; i < sizeof( T ) ; ++i )
	{
		data[ BO::order()[i] ] = reordered.c[ i ] ;
	}
}

} //d6

#endif
