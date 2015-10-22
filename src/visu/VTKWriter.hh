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

	bool startFile( const char* name, unsigned frame ) ;

protected:

	explicit VTKWriter( const char* base_dir ) ;
	virtual ~VTKWriter() {}

	virtual void writeMesh( File &file ) const = 0 ;
	virtual size_t nDataPoints() const = 0 ;

	//! Write attribute header + attribute data
	template< typename Scalar >
	void writeAttribute( const char* name, const Scalar* data, const int Dim ) ;

	//! Writes raw data
	template< typename Scalar >
	void write( File &file, const Scalar* data, const int Dim, const size_t size ) const ;

	File m_file ;

private:

	std::string fileName(  const unsigned frame, const char* dataName ) const ;
	bool open( const unsigned frame, const char* dataName ) ;

	void writeHeader( File &file, const char *title ) const ;
	void writeAttributeHeader( File &file, const int Dim, const char* name ) const ;

	const char* m_base_dir ;
	Mode m_mode ;

};

} //d6

#endif
