#ifndef D6_VTK_WRITER_HH
#define D6_VTK_WRITER_HH

#include <string>
#include "utils/File.hh"

namespace d6 {

class VTKWriter {

public:
	enum Mode { Ascii, Binary } ;

	void setMode( Mode mode ) { m_mode = mode ; }
	Mode mode() const { return m_mode ; }

protected:

	explicit VTKWriter( const char* base_dir ) ;

	std::string fileName(  const unsigned frame, const char* dataName ) const ;
	bool open( const unsigned frame, const char* dataName, File& file ) const ;

	void writeHeader( File &file, const char *title ) const ;

	void writeDataHeader( File &file, const int Dim, const char* name ) const ;

	template< typename Scalar >
	void write( File &file, const Scalar* data, const int Dim, const size_t size ) const ;

	const char* m_base_dir ;
	Mode m_mode ;
};

} //d6

#endif
