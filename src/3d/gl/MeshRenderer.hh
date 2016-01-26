#ifndef D6_MESH_RENDERER_HH
#define D6_MESH_RENDERER_HH

#include "utils/alg.hh"

#include "VertexBuffer.hh"

namespace d6 {

class TriangularMesh ;
class Shader ;

class MeshRenderer 
{
	public:
	
	void reset( const TriangularMesh& mesh ) ;

	void draw( const Shader &shader ) const ;

	bool ok() const 
	{
		return m_vertices.valid() ;
	}

private:
	gl::VertexBuffer3f m_vertices ;
	gl::VertexBuffer3f m_normals ;
	gl::VertexBuffer3f m_uvs ;

} ;

} //d6


#endif

