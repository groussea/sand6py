#ifndef D6_SHAPE_RENDERER_HH
#define D6_SHAPE_RENDERER_HH

#include "utils/alg.hh"

#include "VertexBuffer.hh"
#include "Shader.hh"

namespace d6 {

class LevelSet ;

class ShapeRenderer
{

public:
	void init() ;

	void draw(const LevelSet &ls, const Vec &box, const Eigen::Vector3f &lightPos) const ;


	const gl::VertexBuffer3f& sphereVertices() const
	{ return m_sphereVertices ; }
	const gl::IndexBuffer& sphereQuadIndices() const
	{ return m_sphereQuadIndices ; }

	const gl::VertexBuffer3f& squareVertices() const
	{ return m_squareVertices ; }

private:

	gl::VertexBuffer3f m_sphereVertices ;
	gl::IndexBuffer	   m_sphereQuadIndices ;

	gl::VertexBuffer3f m_squareVertices ;

	Shader m_ballShader ;

};

} //d6

#endif
