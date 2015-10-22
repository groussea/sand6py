#include "VTKWriter.hh"

#include "utils/File.hh"
#include "utils/Log.hh"
#include "utils/endian.hh"

#include "geo/Tensor.hh"

namespace d6 {

VTKWriter::VTKWriter(const char *base_dir)
	: m_base_dir( base_dir ), m_mode( Binary )
{
//	if( !FileInfo( m_base_dir ).exists() ) {
//		Log::Error
//	}
}

std::string VTKWriter::fileName(const unsigned frame, const char *dataName) const
{
	return FileInfo( FileInfo( m_base_dir ).filePath("vtk") ).filePath( arg("%1-%2.vtk", dataName, frame) ) ;
}

bool VTKWriter::open(const unsigned frame, const char *dataName, File &file) const
{
	const std::string fn = fileName(frame, dataName) ;
	FileInfo(fn).makePath() ;
	if( !file.open(fn , std::ios_base::out ) ) {
		Log::Error() << "Could not write into " << file.name() <<std::endl ;
		return false ;
	}

	Log::Info() << "Writing " << fn << std::endl ;

	return true ;
}

void VTKWriter::writeHeader( File &file, const char *title ) const
{
	file << "# vtk DataFile Version 2.0\n" ;
	file << title << "\n" ;
	if( Ascii )
		file << "ASCII\n" ;
	else
		file << "BINARY\n" ;
}

void VTKWriter::writeDataHeader( File& file, const int Dim, const char* name) const
{
	switch( Dim ) {
		case 1:
			file << "SCALARS " << name << " float\n" ;
			file << "LOOKUP_TABLE default\n" ;
			break ;
		case 3:
			file << "VECTORS " << name << " float\n" ;
			break ;
		case 6:
			file << "TENSORS " << name << " float\n" ;
			break ;
	}

}

template< typename T >
struct binary_type
{
	typedef T type ;
} ;
template< >
struct binary_type< double >
{
	typedef float type ;
} ;
template< >
struct binary_type< Index >
{
	typedef int type ;
} ;


template< typename Scalar >
static void write_scalar_ascii( File& file, const Scalar* data, const size_t size )
{
	for( size_t i = 0 ; i < size ; ++ i ) {
		file << data[i] << " " ;
	}
}

template< typename Scalar >
static void write_scalar_binary( File& file, const Scalar* data, const size_t size )
{
	typedef typename binary_type< Scalar >::type BinScalar ;
	unsigned char bin[ sizeof(BinScalar) ] ;

	for( size_t i = 0 ; i < size ; ++ i ) {
		to_big_endian( BinScalar(data[i]), bin ) ;
		for( unsigned j=0 ; j < sizeof(BinScalar) ; ++j )
			file << bin[j] ;
	}
}

template< typename Scalar >
static void write_tensor_ascii( File& file, const Scalar* data, const size_t size )
{
	Mat mat ;
	for( size_t i = 0 ; i < size ; ++ i ) {
		tensor_view( Vec6::Map( data + 6*i ) ).get( mat ) ;
		write_scalar_ascii( file, mat.data(), mat.size() ) ;
	}
}

template< typename Scalar >
static void write_tensor_binary( File& file, const Scalar* data, const size_t size )
{
	Mat mat ;
	for( size_t i = 0 ; i < size ; ++ i ) {
		tensor_view( Vec6::Map( data + 6*i ) ).get( mat ) ;
		write_scalar_binary( file, mat.data(), mat.size() ) ;
	}
}

template< typename Scalar >
void VTKWriter::write( File& file, const Scalar* data, int Dim, const size_t size ) const
{
	if( m_mode == Ascii ) {
		if( Dim == 6 ) {
			write_tensor_ascii( file, data, size ) ;
		} else {
			write_scalar_ascii( file, data, size*Dim ) ;
		}
	} else {
		if( Dim == 6 ) {
			write_tensor_binary( file, data, size ) ;
		} else {
			write_scalar_binary( file, data, size*Dim ) ;
		}
	}

	file << "\n" ;
}

template void VTKWriter::write( File&, const double*, int, size_t) const;

}
