#ifndef D6_VTK_FIELD_WRITER_HH
#define D6_VTK_FIELD_WRITER_HH

#include "VTKWriter.hh"

#include "utils/alg.hh"
#include "geo/geo.fwd.hh"

namespace d6 {


class VTKFieldWriter : public VTKWriter
{

public:

	VTKFieldWriter( const char* base_dir, const MeshType& mesh ) ;

	template< typename Derived >
	bool dump( const char* name, const FieldBase< Derived >& field ) ;

protected:
	void writeMesh(File &file) const ;
	size_t nDataPoints( ) const ;

private:

	const MeshType& m_mesh ;

};

} //d6

#endif

