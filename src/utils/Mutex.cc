#ifndef MUTEXIMPL_HPP
#define MUTEXIMPL_HPP

#include "Mutex.hh"

#ifdef _OPENMP
#include <omp.h>
#else
#include <pthread.h>
#endif

namespace d6 {

class MutexImpl
{

public:

  MutexImpl( Mutex::Type t )
  {
#ifdef _OPENMP
	_type = t ;
	if( _type == Mutex::Recursive )
	{
	  omp_init_nest_lock( &_mutex.nest_lock ) ;
	} else {
	  omp_init_lock( &_mutex.lock ) ;
	}
#else
	pthread_mutexattr_init( &_attr ) ;
	if( t == Mutex::Recursive )
	{
	  pthread_mutexattr_settype(&_attr, PTHREAD_MUTEX_RECURSIVE );
	} else {
	  pthread_mutexattr_settype(&_attr, PTHREAD_MUTEX_NORMAL );
	}
	pthread_mutex_init( &_mutex, &_attr ) ;
#endif

  }

  ~MutexImpl()
  {
#ifdef _OPENMP
	if( _type == Mutex::Recursive )
	{
	  omp_destroy_nest_lock( &_mutex.nest_lock ) ;
	} else {
	  omp_destroy_lock( &_mutex.lock ) ;
	}
#else
	pthread_mutex_destroy( &_mutex ) ;
	pthread_mutexattr_destroy( &_attr ) ;
#endif
  }

  void lock()
  {
#ifdef _OPENMP
	if( _type == Mutex::Recursive )
	{
	  omp_set_nest_lock( &_mutex.nest_lock ) ;
	} else {
	  omp_set_lock( &_mutex.lock ) ;
	}
#else
	pthread_mutex_lock( &_mutex ) ;
#endif
  }

  bool try_lock()
  {
#ifdef _OPENMP
	if( _type == Mutex::Recursive )
	{
	  return 0 != omp_test_nest_lock( &_mutex.nest_lock ) ;
	} else {
	  return 0 != omp_test_lock( &_mutex.lock ) ;
	}
#else
	return !pthread_mutex_trylock( &_mutex ) ;
#endif
  }

  void unlock()
  {
#ifdef _OPENMP
	if( _type == Mutex::Recursive )
	{
	  omp_unset_nest_lock( &_mutex.nest_lock ) ;
	} else {
	  omp_unset_lock( &_mutex.lock ) ;
	}
#else
	pthread_mutex_unlock( &_mutex ) ;
#endif
  }

private:

  MutexImpl( const Mutex& ) ;
  MutexImpl& operator= ( const Mutex& ) ;

#ifdef _OPENMP
  Mutex::Type _type ;
  union {
	omp_lock_t lock ;
	omp_nest_lock_t nest_lock ;
  } _mutex ;
#else
  pthread_mutex_t _mutex ;
  pthread_mutexattr_t _attr ;
#endif

} ;

Mutex::Mutex( Mutex::Type t )
  : _impl( new MutexImpl(t) )
{
}

Mutex::~Mutex( )
{
  delete _impl ;
}

void Mutex::lock() const
{
  _impl->lock();
}

bool Mutex::try_lock() const
{
  return _impl->try_lock();
}

void Mutex::unlock() const
{
  _impl->unlock();
}


} //namespace d6

#endif // MUTEXIMPL_HPP
