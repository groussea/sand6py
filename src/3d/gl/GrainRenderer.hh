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

#ifndef D6_GRAIN_RENDERER_HH
#define D6_GRAIN_RENDERER_HH

#include "visu/Sampler.hh"

#include "VertexBuffer.hh"
#include "Shader.hh"

namespace d6 {

struct Texture ;
class ShapeRenderer ;

class GrainRenderer
{

public:
	GrainRenderer( const Offline & offline, const ShapeRenderer &shapeRenderer, unsigned nSamples )
		: m_offline( offline ), m_shapeRenderer( shapeRenderer ), m_sampler(offline),
		  m_nSamples( nSamples ),  m_grainSizeFactor( 1 )
	{

	}

	void setGrainSizeFactor( float factor )
	{
		m_grainSizeFactor = factor ;
	}
	void cutAndColorVelocities() {
		m_sampler.setMode( Sampler::VelocityCut );
	}
	void useDiscs() {
		m_sampler.setMode( Sampler::Discs );
	}

	const Sampler& sampler() const { return m_sampler ; }

	bool valid() const {
		return m_nSamples > 0 ;
	}

	void compute_shadow( const float baseGrainSize, const Eigen::Matrix4f &depthMVP ) ;
	void draw(const Texture &depthTexture, const Eigen::Vector3f& lightPosition,
			  const float baseGrainSize, const Eigen::Matrix4f &depthMVP ) ;

	void init() ;
	void update_buffers() ;

	void move() ;

private:

	void draw_grains ( const Shader &shader, const float pixelSize,
					   const Eigen::Matrix4f &depthMVP, bool instanced ) const ;

	const Offline& m_offline ;
	const ShapeRenderer &m_shapeRenderer ;

	Sampler  m_sampler ;

	unsigned m_nSamples  ;

	float 	 m_grainSizeFactor ;

	gl::VertexBuffer3f m_grainVertices ;
	gl::VertexBuffer3f m_grainNormals ;
	gl::ArrayBufferf m_grainVisibility ;
	gl::ArrayBufferf m_grainNoise ;

	Shader m_grainsShader ;
	Shader m_depthShader ;
} ;



} //d6

#endif
