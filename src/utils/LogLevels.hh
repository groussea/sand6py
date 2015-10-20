#ifndef D6_LOGLEVELS_HH
#define D6_LOGLEVELS_HH
/*!
 * \file Loglevels.hpp
 * \brief Definition of the output verbosity levels
 *
 * \sa Log.hh
 *
 */


//! Logs levels, ordered from more verbose to less verbose
/*! For each of those levels, one typedef will be defined in the namespace Log.
  It will be then possible to use it as an out stream. \sa Log */
#define D6_LOG_LEVELS \
  D6_LOG_LEVEL( All     ) \
  D6_LOG_LEVEL( Debug   ) \
  D6_LOG_LEVEL( Verbose ) \
  D6_LOG_LEVEL( Info    ) \
  D6_LOG_LEVEL( Warning ) \
  D6_LOG_LEVEL( Error   )

namespace d6 {
namespace Log {

//! Enum containing all log levels, prefixed by 'L_' ( e.g L_Debug )
enum Level {
#define D6_LOG_LEVEL( level ) L_##level,
  D6_LOG_LEVELS
#undef D6_LOG_LEVEL
  L_Nothing
} ;

} // Log

} // namespace d6

#endif // LOGLEVELS_DEF_D6_H

