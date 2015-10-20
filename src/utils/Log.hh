#ifndef LOG_D6_HPP
#define LOG_D6_HPP

/*!
 * \file Log.hpp
 * \brief Filtered, thread-safe, textual output mechanism
 * The stream will be exclusively owned by one thread until
 * it is flushed, for example by appending std::endl.
 *
 * /!\ Not flushing the streams will cause deadlocks !
 *
 */

//! If defined, do not ensure thread-safety of textual output
//#define D6_LOG_NO_THREAD_SAFETY

#include "LogLevels.hh"

#ifndef D6_LOG_NO_THREAD_SAFETY
#include "Mutex.hh"
#endif

#include "ThreadLocal.hh"

#include <iostream>
#include <sstream>

namespace d6 {

//! Log-related definitions
/*! Usage
  \code
  Log::Info() << " This is message is relatively important" << std::endl ;
  Log::Debug()::as_std_ostream() << " This is much less so " << std::endl ;
  \endcode
  */
namespace Log
{

  //! Config : Defines verbosity level and output streams
  /*! Should be accessed only through Config::get() */
  struct Config
  {
	//! Normal output stream ( defaults to Log::Debug() )
	std::ostream* out ;
	//! Error output stream ( defaults to Log::Error() )
	std::ostream* err ;
	//! Output level ( defaults to L_All )
	Level level ;
	//! If true, do not prefix output with '[StreamLetter]'
	bool no_prefix ;

	//! Returns the static global instance
	static Config& get()
	{
	  static Config s_config ;
	  return s_config ;
	}

  private:
	//! Forbid construction outside of class
	Config() : out( NULL ), err(NULL), level( L_All ), no_prefix( false )
	{}

  } ;

  namespace log_impl {
  template < Level L > struct Log ;
  //! Should accessed directly with the typedefs below
  template < Level L > struct Stream {

	//! Filtered, thread safe operator<<
	template< typename T >
	Log< L >& operator<< (  const T& o ) const ;

	//! Filtered, thread safe operator<< for stream manipulators
	Log< L >& operator<< ( std::ostream& modifier( std::ostream& )  ) const ;

	//! Returns an almost thread safe std::ostream
	/*! Actually, the rarely used sputc() method of the streambuf won't be thread safe */
	static std::ostream& as_std_ostream( ) ;

	//! Thread-unsafe stream, with should not perform any blocking call to malloc
	static std::ostream& unsafe( ) ;

  private:
	static Log< L >& impl() ;
  } ;
  }

  //! Convenience typedef for each log level
#define D6_LOG_LEVEL( level ) typedef log_impl::Stream< L_##level > level ;
	D6_LOG_LEVELS
	typedef log_impl::Stream< L_Nothing > Nothing ;
#undef D6_LOG_LEVEL

  //! Full and abbreviated names
  template< Level L > struct  LevelNames {
	static const char* full() { return "Nothing" ; }
	static char abbrev() { return  'N' ; }
  } ;

  // Specializations for each Level
#define D6_LOG_LEVEL( level ) \
template< > struct  LevelNames< L_##level > { \
  static const char* full() { return #level ; } \
  static char abbrev() { return  #level[0] ; } \
} ;
  D6_LOG_LEVELS
#undef D6_LOG_LEVEL

  namespace log_impl
  {

   struct LogBase
   {

	 explicit LogBase( std::stringbuf *buf )
	   : stream( buf )
	 {
	   Config &c = Config::get() ;
	   if( !c.out )
	   {
		 c.out = &std::cout ;
	   }
	   if( !c.err )
	   {
		 c.err = &std::cerr ;
	   }
	 }

	template< Level L >
	 static const char * prefix( )
	 {
	  static const char s_prefix[5] = { '[', LevelNames<L>::abbrev(), ']', ' ', '\0' };

	  return ( Config::get().no_prefix ) ? "" : s_prefix ;
	 }

#ifndef D6_LOG_NO_THREAD_SAFETY
	 static Mutex& global_mutex()
	 {
	   static Mutex s_global( Mutex::Recursive ) ;
	   return s_global ;
	 }
#endif

	std::ostream stream ;

  } ;

  template < Level L   >
   struct Buf : public std::stringbuf {
	 std::ostream* &stream ;
#ifndef D6_LOG_NO_THREAD_SAFETY
	 Mutex mutex ;
	 volatile bool locked ;
#endif

	 Buf( std::ostream* &stream_ptr )
	   : stream( stream_ptr )
#ifndef D6_LOG_NO_THREAD_SAFETY
	   , mutex( Mutex::Recursive )
  #endif
	 {}

	 //! Called each time the stream is flushed ( e.g. for each std::endl )
	 virtual int sync ( )
	 {
	   {
		 //Output prefix + buffered string
#ifndef D6_LOG_NO_THREAD_SAFETY
		 LockGuard lock( LogBase::global_mutex() ) ;
#endif
		 *stream << LogBase::prefix<L>( ) << str();
		 stream->flush();
	   }

	   str("");

		// Release mutex
#ifndef D6_LOG_NO_THREAD_SAFETY
	   if( locked )
	   {
		 locked = false ;
		 mutex.unlock() ;
	   }
#endif

	   return 0;
	 }

#ifndef D6_LOG_NO_THREAD_SAFETY
	 virtual std::streamsize xsputn (const char* s, std::streamsize n)
	 {
	   acquire_mutex();
	   return std::stringbuf::xsputn( s, n ) ;
	 }
#endif


	 void acquire_mutex()
	 {
#ifndef D6_LOG_NO_THREAD_SAFETY
	   mutex.lock() ;
	   if( locked ) mutex.unlock() ; //This thread was already owning the recursive mutex, release it once

	   locked = true ;
#endif
	 }

   } ;

  template < Level L   >
  struct Log : public LogBase {
	Buf< L > buffer ;
	Log() : LogBase( &buffer ), buffer ( Config::get().out ) {}
	Log& thread_safe() { buffer.acquire_mutex() ; return *this ; }
	std::ostream& unsafe() { return  *buffer.stream << LogBase::prefix< L >() ; }
  } ;

  template < >
  struct Log< L_Error > : public LogBase {
	Buf< L_Error > buffer ;
	Log() : LogBase( &buffer ), buffer ( Config::get().err ) {}
	Log& thread_safe() { buffer.acquire_mutex() ; return *this ; }
	std::ostream& unsafe() { return *buffer.stream << LogBase::prefix< L_Error >() ; }
  } ;

  template < Level L, typename T >
  Log< L > & operator << ( Log< L > & s, const T& o )
  {
	if( L < Config::get().level ) return s ;
	s.thread_safe().stream << o ;
	return s ;
  }

  template < Level L >
  Log< L > & operator << ( Log< L > & s, std::ostream& modifier( std::ostream& )  )
  {
	if( L < Config::get().level ) return s ;
	s.thread_safe().stream << modifier ;
	return s ;
  }


  template < > struct Stream< L_Nothing > : public std::ostream {

	Stream(): std::ostream(0)
	{
	}

	static std::ostream& as_std_ostream( )
	{
	  D6_TLS( Stream< L_Nothing >, black_hole ) ;
	  return *black_hole ;
	}

	static std::ostream& unsafe( )
	{
	  return as_std_ostream() ;
	}
  };


  //! Filtered, thread safe operator<<
  template< Level L >
  template< typename T >
  Log< L >& Stream< L >::operator<< (  const T& o ) const
  {
	return ( impl() << o ) ;
  }

  //! Filtered, thread safe operator<< for stream manipulators
  template< Level L >
  Log< L >& Stream< L >::operator<< ( std::ostream& modifier( std::ostream& )  ) const
  {
	return ( impl() << modifier ) ;
  }

  //! Returns an almost thread safe std::ostream
  /*! Actually, the rarely used sputc() method of the streambuf won't be thread safe */
  template< Level L >
  std::ostream& Stream< L >::as_std_ostream( )
  {
	if ( L < Config::get().level )
	  return Nothing::as_std_ostream() ;

	return impl().thread_safe().stream ;
  }

  //! Thread-unsafe stream, with should not perform any blocking call to malloc
  template< Level L >
  std::ostream& Stream< L >::unsafe( )
  {
	if ( L < Config::get().level )
	  return Nothing::as_std_ostream() ;

	return impl().unsafe() ;
  }

  template< Level L >
  Log< L >& Stream< L >::impl()
  {
	static Log< L > s_impl ;
	return s_impl ;
  }


  } // log_impl

} // Log

} // namespace d6



#endif // LOG_D6_HPP
