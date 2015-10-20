#ifndef D6_MUTEX_HH
#define D6_MUTEX_HH

/*!
 * \file Mutex.hpp
 * \brief RAII wrapper around pthread's mutexes
 *
 */

namespace d6 {

class MutexImpl ;

//! Simple wrapper around pthread/OpenMP mutexes
class Mutex
{

public:

  enum Type  {
	  Default,
	  Recursive
  }  ;

  explicit Mutex( Type t = Default ) ;

  ~Mutex() ;

  void lock() const ;

  bool try_lock() const ;

  void unlock() const ;

private:

  Mutex( const Mutex& ) ;
  Mutex& operator= ( const Mutex& ) ;

  MutexImpl* const _impl ;
} ;

//! Scoped lock
class LockGuard
{
public:
  LockGuard()
	: _owned( new Mutex() ), _mutex( *_owned )
  {
	_mutex.lock() ;
  }

  explicit LockGuard( const Mutex& mutex )
	:  _owned( 0 ), _mutex( mutex )
  {
	_mutex.lock() ;
  }

  ~LockGuard()
  {
	_mutex.unlock() ;
	delete _owned ;
  }

private:
  LockGuard( const LockGuard& ) ;
  LockGuard& operator= ( const LockGuard& ) ;

  Mutex* _owned ;
  const Mutex &_mutex ;
} ;

} //namespace d6

#endif // MUTEX_D6_HPP
