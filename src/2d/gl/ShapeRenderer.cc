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

#include "ShapeRenderer.hh"

#include "geo/LevelSet.impl.hh"

#include "opengl.hh"

namespace d6 {



void ShapeRenderer::init()
{

}


void ShapeRenderer::draw( const LevelSet &ls, const Vec &box ) const
{
    const Eigen::Vector2f translation = ls.origin().cast < GLfloat >() ;
    Eigen::Matrix2f rotation ;
    {
        const Scalar c = std::cos( ls.rotation() ) ;
        const Scalar s = std::sin( ls.rotation()  ) ;
        rotation << c, -s,
                    s,  c;
    }

    Eigen::Matrix4f mat = Eigen::Matrix4f::Identity() ;
    mat.block<2,2>(0,0) = rotation * ls.scale() ;
    mat.block<2,1>(0,3) = translation ;

    glColor4f(1., 0., .8, 1);

    glPushMatrix();
    glMultMatrixf( mat.data() );

    const SegLevelSet* seg = dynamic_cast< const SegLevelSet* >(&ls) ;

    if( dynamic_cast<const SphereLevelSet*>(&ls) )  {
        const unsigned res = 64 ;

        glBegin( GL_LINE_STRIP ) ;
        for( unsigned k = 0 ; k < res ; ++k )
        {
            const float t = 2. * k * M_PI / (res - 1) ;
            glVertex2d( std::cos(t), std::sin(t) ) ;
        }
        glEnd() ;

        glBegin( GL_LINES ) ;
        glVertex2d( 1, 0 ) ;
        glVertex2d( 0, 0 ) ;
        glEnd() ;
    } else if ( dynamic_cast<const PlaneLevelSet*>(&ls) ) {

        glBegin( GL_LINES );
        glVertex3f( -box[0], 0, 0 );
        glVertex3f(  box[0], 0, 0 );
        glEnd( ) ;

    } else if ( seg ) {

        glBegin( GL_LINES );
        glVertex3f( 0, -.5*seg->len(), 0 );
        glVertex3f( 0,  .5*seg->len(), 0 );
        glEnd( ) ;

    }


    glPopMatrix();


}

} // d6
