#ifndef D6_THREADLOCAL_HH
#define D6_THREADLOCAL_HH

/*!
 * \file ThreadLocal.hh
 * \brief Simple utility struct to use TLS with non-basic types
 *
 * Uses placement new operator on thread-local buffer ( so no memleak will occur ), BUT:
 * /!\ \b Warning: The object's destructor will \b NEVER by called. Don't rely on it.
 *
 * Usage:
  \code
	//Declare a thread local VectorXd name x
	D6_TLS( Eigen::VectorXd, x ) = Eigen::VectorXd::Zero( 10 ) ;
	// Access via operator-> or operator*
	*x += Eigen::VectorXd::Ones( 10 ) ;
	std::cout << x->transpose() << std::endl ;

  \endcode
 *
 */

#ifdef _MSC_VER
#define D6_TLS_SPEC __declspec( thread )
#else
#define D6_TLS_SPEC __thread
#endif

#define D6_TLS( Type, Name ) \
  static D6_TLS_SPEC char s_data_##Name[ sizeof( Type ) ] ; \
  static D6_TLS_SPEC char s_state_##Name  = 0 ; \
  tls_impl::ThreadLocal< Type > Name( s_data_##Name, s_state_##Name )

namespace d6 {

namespace tls_impl
{

template < typename T >
class ThreadLocal
{

  enum State
  {
	UNINITIALIZED = 0,
	CONSTRUCTED,
	SET
  } ;

public:

  ThreadLocal( char * data, char &state )
	: _data( data ), _state( state )
  {
	if( ((char) UNINITIALIZED ) == _state )
	{
	  new ( _data ) T() ;
	  _state = (char) CONSTRUCTED ;
	}
  }

  //! Assignement operator that will ony be exectued once
  T& operator=( const T& rhs )
  {
	T& lhs = *get() ;
	if( _state == (char) CONSTRUCTED ){
	  lhs = rhs ;
	  _state = (char) SET ;
	}
	return lhs ;
  }

  T& operator * ( ) {
	return *get();
  }
  T const& operator * ( ) const {
	return *get();
  }
  T* operator -> ( ) {
	return get();
  }
  T const* operator -> ( ) const {
	return get();
  }

private:
  T* get() const
  {
	return reinterpret_cast< T* > ( _data ) ;
  }

  void* raw() const
  {
	return reinterpret_cast< void* > ( _data ) ;
  }

  char *_data ;
  char &_state ;

} ;

}

} //namespace d6

#endif // THREADLOCAL_D6_HPP
