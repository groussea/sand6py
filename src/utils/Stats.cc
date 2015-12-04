#include "Stats.hh"

#include <iomanip>

namespace d6 {

template < typename T >
struct FieldWidth {
	static constexpr int width = 7 ;
};
template < >
struct FieldWidth< Scalar > {
	static constexpr int width = 11 ;
};

Stats::Stats(const char *base_dir)
	: m_out( FileInfo(base_dir).filePath("stats.txt"), std::ios_base::out )
#define STAT_FIELD( name, type, abv ) \
	, name(0)
	EXPAND_STAT
#undef STAT_FIELD
{

#define STAT_FIELD( name, type, abv ) \
	m_out << std::setw( FieldWidth<type>::width ) << abv << " " ;
	EXPAND_STAT
#undef STAT_FIELD

	m_out << std::endl ;

}

void Stats::dump() {

#define STAT_FIELD( name, type, abv ) \
	m_out << std::setw( FieldWidth<type>::width ) << name << " " ;
	EXPAND_STAT
#undef STAT_FIELD

   m_out << std::endl ;

	++ stepId ;
}

} //d6
