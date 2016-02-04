#ifndef D6_VTK_FIELD_WRITER_HH
#define D6_VTK_FIELD_WRITER_HH

#include "VTKWriter.hh"

#include "utils/alg.hh"
#include "geo/geo.fwd.hh"

namespace d6 {

template < typename ShapeFuncT >
class VTKFieldWriter : public VTKWriter
{

public:

	VTKFieldWriter( const char* base_dir, const ShapeFuncT& shape ) ;

	template< typename Derived >
	bool dump( const char* name, const FieldBase< Derived >& field ) ;

protected:
	void writeMesh(File &file) const ;
	size_t nDataPoints( ) const ;

private:
	const ShapeFuncT& m_shape ;

};

} //d6

#endif

