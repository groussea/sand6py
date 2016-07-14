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

#ifndef D6_GLVIEWER_HH
#define D6_GLVIEWER_HH

#include "mono/PhaseFields.hh"

#include "gl/VertexBuffer.hh"
#include "gl/Shader.hh"
#include "gl/ShapeRenderer.hh"

#include <vector>

namespace d6 {

class Offline ;

class GLViewer
{
public:

	enum Element
	{
		eGrid,
		eColors,
		eVectors,
		eTensors,
		eStream,
		eParticles,
		nElements
	};

	enum ScalarEntity {
		sePressure,
		seFraction,
		seSpin,
		seDivergence,
		nSE
	} ;

	enum VectorEntity {
		veVelocity,
		veGradPhi,
		veProj,
		veForce,
		nVE
	} ;

	enum TensorEntity {
		teStress,
		teDu,
		nTE
	} ;

	GLViewer( const Offline& offline,
			const int width, const int height )
		: m_offline( offline ), m_width( width ), m_height( height ),
		  m_xOffset( 0 ), m_yOffset( 0 ), m_zoom( 1 ), m_snapId(0),
		  m_drawOrientations( false ),
		  m_scalarEntity( seFraction ), m_vectorEntity( veVelocity ),
		  m_tensorEntity( teDu )
	{

	}

	int width() const { return m_width ; }
	int height() const { return m_height ; }

	void init( ) ;

	void update_buffers( ) ;
	void update_particle_buffers( ) ;
	void update_color_buffers( ) ;
	void update_vector_buffers( ) ;
	void update_tensor_buffers( ) ;
	void update_texture( ) ;

	void draw( ) const ;

	void toggleRendering( Element elt )
	{
		m_shouldRender[ elt ] = !m_shouldRender[ elt ] ;
	}

	void toggleParticleRepr()
	{
		m_drawOrientations = !m_drawOrientations ;
	}

	void toggleScalarEntity()
	{
		m_scalarEntity = (ScalarEntity)( (((unsigned)m_scalarEntity) + 1 ) % nSE ) ;
	}
	void toggleVectorEntity()
	{
		m_vectorEntity = (VectorEntity)( (((unsigned)m_vectorEntity) + 1 ) % nVE ) ;
	}
	void toggleTensorEntity()
	{
		m_tensorEntity = (TensorEntity)( (((unsigned)m_tensorEntity) + 1 ) % nVE ) ;
	}

	void move( double x, double y ) ;

	void zoom( double factor ) ;

	void updateViewport() ;

	void snap() ;

private:
	typedef PrimalShape Shape ;
	typedef AbstractScalarField< Shape > ScalarField ;
	typedef AbstractVectorField< Shape > VectorField ;
	typedef AbstractTensorField< Shape > TensorField ;

	const VectorField &getVectorEntity() const ;
	ScalarField getScalarEntity() const ;
	TensorField getTensorEntity() const ;

	const Offline& m_offline ;

	int m_width ;
	int m_height ;

	double m_xOffset ;
	double m_yOffset ;
	double m_zoom ;

	unsigned m_snapId ;

	bool m_shouldRender[ nElements ] ;
	bool m_drawOrientations ;

	ScalarEntity m_scalarEntity ;
	VectorEntity m_vectorEntity ;
	TensorEntity m_tensorEntity ;

	ShapeRenderer m_shapeRenderer ;

	Shader	m_particleShader ;
	Shader	m_vectorShader ;
	Shader	m_tensorShader ;

	//Grid
	gl::VertexBuffer2f m_gridVertices ;
	gl::IndexBuffer    m_gridQuadIndices ;
	int				   m_gridPrimitive ;

	//Color field
	gl::VertexBuffer3f m_gridColors ;

	//Vectors
	gl::VertexBuffer2f m_arrows ;

	//Streamlines
	GLuint m_texId ;
	std::vector< float > m_texData ;
	std::vector< float > m_rndData ;

	//Tensor
	gl::VertexBuffer3f m_tensors ;

	//particles
	gl::VertexBuffer2d m_particles ;
	gl::VertexBuffer4f m_particleFrames ;
	gl::ArrayBufferf   m_particleAlpha ;
	gl::ArrayBufferf   m_particleColors ;
};

}

#endif
