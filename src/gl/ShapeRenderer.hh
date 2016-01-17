#ifndef D6_SHAPE_RENDERER_HH
#define D6_SHAPE_RENDERER_HH

#include "utils/alg.hh"

#include "VertexBuffer.hh"
#include "Shader.hh"

#include <string>
#include <unordered_map>

namespace d6 {

class LevelSet ;
class MeshRenderer ;
class TriangularMesh ;
struct Texture ;

class ShapeRenderer
{

public:
	void init() ;

	void compute_shadow( const LevelSet &ls, 
		const Eigen::Matrix4f& depthModelView, const Eigen::Matrix4f& depthProjection ) const ;
	void draw(const LevelSet &ls, const Vec &box, const Eigen::Vector3f &lightPos,
		bool shadowed, const Texture& depthTexture,
		const Eigen::Matrix4f& depthModelView, const Eigen::Matrix4f& depthProjection ) const ;


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
	Shader m_ballDepthShader ;
	Shader m_solidShader ;
	Shader m_solidDepthShader ;
};

} //d6

#endif
