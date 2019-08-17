/*
 * This file is part of Sand6, a C++ continuum-based granular simulator.
 *
 * Copyright 2016 Gilles Daviet <gilles.daviet@inria.fr> (Inria - Université Grenoble Alpes)
 *
 * Sand6 is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.

 * Sand6 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.

 * You should have received a copy of the GNU General Public License
 * along with Sand6.  If not, see <http://www.gnu.org/licenses/>.
*/

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
	void config_shaders() ;
	void setup_vaos() ;

	void compute_shadow( const LevelSet &ls, 
		const Eigen::Matrix4f& depthModelView, const Eigen::Matrix4f& depthProjection ) const ;
	void draw(const LevelSet &ls, const Vec &box, const Eigen::Vector3f &lightPos,
		bool shadowed, const Texture& depthTexture,
		const Eigen::Matrix4f& modelView, const Eigen::Matrix4f& projection, 
		const Eigen::Matrix4f& depthModelView, const Eigen::Matrix4f& depthProjection ) const ;


	const gl::VertexBuffer3f& sphereVertices() const
	{ return m_sphereVertices ; }
	const gl::IndexBuffer& sphereTriIndices() const
	{ return m_sphereTriIndices ; }
	const gl::ArrayObject& sphereVertexArrays() const
	{ return m_sphereVertexArrays ; }

	const gl::VertexBuffer3f& squareVertices() const
	{ return m_squareVertices ; }

private:

	gl::VertexBuffer3f m_sphereVertices ;
	gl::IndexBuffer	   m_sphereTriIndices ;
	gl::ArrayObject    m_sphereVertexArrays;

	gl::VertexBuffer3f m_squareVertices ;
	gl::ArrayObject    m_billboardArrays;

	Shader m_ballShader ;
	Shader m_ballDepthShader ;
	Shader m_solidShader ;
	Shader m_solidDepthShader ;
};

} //d6

#endif
