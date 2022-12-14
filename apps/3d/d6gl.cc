
/*
 * This file is part of Sand6, a C++ continuum-based granular simulator.
 *
 * Copyright 2019 Gilles Daviet <gilles.daviet@inria.fr> (Inria - Universit√© Grenoble Alpes)
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

#include "GLViewer.hh"

#include "utils/string.hh"
#include "utils/Config.hh"

#include "visu/Offline.hh"

#include <GLFW/glfw3.h>

#include <cstdlib>

namespace d6
{

class Appli
{

public:
	Appli(
      d6::Offline& offline, unsigned frame,
	  	const int nSamples,
	       const int width, const int height , bool snaps, bool run, bool anim) :
	    m_pWindow( NULL ),
	    m_offline(offline), m_viewer( offline, nSamples, width, height ),
	    m_currentFrame( frame ),
	    m_running( run ), m_snapshotting( snaps ), m_anim(anim),
	    m_mouseX( 0 ), m_mouseY( 0 )
	{
		assert( ! s_instance ) ;
		s_instance = this ;
	}

	~Appli()
	{
		quit(0) ;
	}

	GLViewer& viewer() { return m_viewer; }

	void set_frame( unsigned frame ) {
		if( m_offline.load_frame( frame ) ) {
			std::cout << m_offline.particles().count() << std::endl;
			m_currentFrame = frame ;
			m_viewer.update_buffers();

			}
		}
	


	bool next_frame() {
		unsigned nextFrame = m_currentFrame + 1  ;
		set_frame( nextFrame ) ;
		return nextFrame == m_currentFrame ;
	}
	bool prev_frame() {
		unsigned nextFrame = m_currentFrame - 1  ;
		if( m_currentFrame > 0 )
			set_frame( nextFrame ) ;
		return nextFrame == m_currentFrame ;
	}

	int run( )
	{ 
		m_viewer.frameAll();
		Eigen::Vector3f lookPos;
		Eigen::Vector3f box = m_offline.box().cast<float>();
        lookPos[0]=0.4 * box[0];
        lookPos[1]=0.5 * box[1];
        lookPos[2]=0.5* box[2];
		m_viewer.look_at(lookPos);
        Eigen::Vector3f position;
        position[0] = 0.7 * box[0];
        position[1] = -3. * box[1];
        position[2] = 0.4 * box[2];
		m_viewer.set_cam_pos(position);
		m_viewer.zoom(0.7) ;

		do
		{
			if( m_running && !next_frame() ) {
				m_running = false ;
			}
			if (m_anim ) {	
				m_viewer.rotate(0.1,0.05 ) ;
				m_viewer.zoom(0.999) ;}
			m_viewer.draw( ) ;
		if(  m_snapshotting ) {
				m_viewer.snap(m_currentFrame);}

			glfwSwapBuffers( m_pWindow);
			glfwPollEvents();
		}
		while( !glfwWindowShouldClose(m_pWindow)
		       && ( glfwGetKey(m_pWindow, GLFW_KEY_ESCAPE ) != GLFW_PRESS )
		       && ( glfwGetKey(m_pWindow, GLFW_KEY_Q) != GLFW_PRESS )) ;
		return 0 ;
	}

	static Appli& instance()
	{
		return *s_instance ;
	}

	static void key_callback(GLFWwindow *, int key, int /*scancode*/, int action, int mods )
	{
		instance().process_key( key, action, mods ) ;
	}

	static void mouse_click_callback(GLFWwindow *, int button, int action, int /*mods*/ )
	{
		instance().process_mouse_click( button, action );
	}

	static void mouse_motion_callback(GLFWwindow *, double x, double y )
	{
		instance().process_mouse_motion( x, y );
	}

	static void mouse_scroll_callback(GLFWwindow *, double, double y_offset )
	{
		instance().process_mouse_scroll( y_offset );
	}

	void init()
	{
		if ( !glfwInit() )
			quit(1);

		glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
		glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
		glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
		glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);

		if(!( m_pWindow = glfwCreateWindow(m_viewer.width(), m_viewer.height(), "Sand6GL", NULL, NULL) ))
		{
			quit(1);
		}

		glfwMakeContextCurrent(m_pWindow);

		glfwSetInputMode( m_pWindow, GLFW_STICKY_KEYS, GL_TRUE ) ;
		glfwSetKeyCallback( m_pWindow, &Appli::key_callback ) ;
		glfwSetMouseButtonCallback( m_pWindow, &Appli::mouse_click_callback ) ;
		glfwSetScrollCallback( m_pWindow, &Appli::mouse_scroll_callback ) ;

		std::cout << "Using OpenGL " << glGetString(GL_VERSION)
		          << " rendered on " << glGetString(GL_RENDERER)
		          << " from " << glGetString(GL_VENDOR)
		          << std::endl  ;

		set_frame( m_currentFrame );
		m_viewer.init( ) ;

	}

private:

	void quit(int rc)
	{
		if( m_pWindow) glfwDestroyWindow(m_pWindow);

		std::clog << "Bye. [code " << rc << "]" << std::endl ;
		glfwTerminate();
		std::exit(rc);
	}

	void process_key( int key, int action, int mods )
	{
		if(action == GLFW_PRESS) m_keyMods |= mods;
		if(action == GLFW_RELEASE) m_keyMods = 0;

		if( action != GLFW_RELEASE && action != GLFW_REPEAT )
			return ;

		switch( key )
		{
		case GLFW_KEY_ENTER:
			m_running = !m_running ;
			break;
		case 'A':
			m_fastDrawForCameraMotion = !m_fastDrawForCameraMotion;
			break;
		case 'C':
			break;
		case 'D':
			break;
		case 'F':
			m_viewer.frameAll();
			break;
		case 'G':
			break;
		case 'I':
			if( ! m_running) next_frame();
			break ;
		case 'M':
			break;
		case 'O':
			break;
		case 'P':
			if( ! m_running) prev_frame();
			break;
		case 'R':
			m_snapshotting = ! m_snapshotting ;
			break;
		case 'S':
			break;
		case 'T':
			break;
		case 'V':
			break;
		case GLFW_KEY_HOME:
			set_frame(0);
			break;
		case GLFW_KEY_PERIOD:
			m_viewer.snap(m_currentFrame ) ;
			break ;
		case GLFW_KEY_LEFT:
			m_viewer.translate(10., 0.);
			break ;	
		case GLFW_KEY_RIGHT:
			m_viewer.translate(-10., 0.);
			break ;	
		case GLFW_KEY_UP:
			m_viewer.translate(0, 10.);
			break ;	
		case GLFW_KEY_DOWN:
			m_viewer.translate(0, -10.);
			break ;	

		}

	}

	void process_mouse_motion( double x, double y )
	{
		bool translate = m_keyMods & GLFW_MOD_SUPER;
		if (translate)
			m_viewer.translate(x - m_mouseX, y - m_mouseY);
		else
			m_viewer.rotate(x - m_mouseX, y - m_mouseY);
		m_mouseX = x ; m_mouseY = y ;
		Eigen::Vector3f const  pos = m_viewer.get_cam_pos();
		std::cout << "pos_cam: "<< pos << std::endl;
	}

	void process_mouse_click( int button, int action )
	{
		if( button != GLFW_MOUSE_BUTTON_1 )
			return ;

		double x, y ;
		glfwGetCursorPos( m_pWindow, &x, &y ) ;
		m_mouseX = x ; m_mouseY = y ;
		if( action == GLFW_PRESS )
		{
			glfwSetCursorPosCallback( m_pWindow, &Appli::mouse_motion_callback ) ;
			if(m_fastDrawForCameraMotion) m_viewer.setFastDraw(true);
		} else {
			glfwSetCursorPosCallback( m_pWindow, NULL ) ;
			m_viewer.setFastDraw(false);
		}



	}

	void process_mouse_scroll( double offset )
	{
		if( offset > 0 )
			m_viewer.zoom( 0.9 ) ;
		else if( offset < 0 )
			m_viewer.zoom( 1.0/0.9 ) ;
	}


	GLFWwindow* m_pWindow ;

	Offline& m_offline ;
	GLViewer m_viewer ;

	unsigned m_currentFrame ;
	bool m_running ;
	bool m_snapshotting ;
	bool m_anim ;

	static Appli* s_instance ;

	int m_keyMods = 0;
	double m_mouseX, m_mouseY ;
	bool m_fastDrawForCameraMotion = false;
};

Appli* Appli::s_instance = NULL ;

} //d6

static void usage( const char *name )
{
	std::cout << "Usage: " << name
	          << " [sim_dir=out] [options] "
	          << "\n3D OpenGL viewer."
	          << "\n\n" ;

	std::cout << "Options:\n"
	          << "-? \t Display this help message and exit\n"
	          << "-n frame_id \t Jump to frame frame_id \n"
	          << "-w width \t Specify viewport width\n"
	          << "-h height \t Specify viewport height\n"
	          << std::endl ;
}

int main( int argc, const char * argv[] )
{


	const char * base_dir = "out" ;
	unsigned frame = 0 ;

	unsigned width  = 1600 ;
	unsigned height = 900 ;
	unsigned nSamples = 0;

	bool discs  = false ;
	bool colVel = false;
	bool sn = false;
	bool run = false;
	bool anim = false;
	float grainSizeFactor = 1;

	for( int i = 1 ; i < argc ; ++i )
	{
		if( argv[i][0] == '-' ){
			switch(argv[i][1]) {
			case '?':
				usage( argv[0]) ;
				return 0;
			case 'n':
				if( ++i == argc ) break ;
				frame = d6::to_uint( argv[i] ) ;
				break ;
			case 'w':
				if( ++i == argc ) break ;
				width = d6::to_uint( argv[i] ) ;
				break;
			case 'h':
				if( ++i == argc ) break ;
				height = d6::to_uint( argv[i] ) ;
				break;
			case 's':
				if( ++i == argc ) break ;
				nSamples = d6::to_uint( argv[i] ) ;
				break;
			case 'd':
				discs = true ;
				break;
			case 'v':
				colVel=true;
				break;
			case 'g':
				if( ++i == argc ) break ;
				grainSizeFactor = d6::to_float( argv[i] ) ;
				break;
			case 'S':
				sn=true;
				break;
			case 'r':
				run=true;
				break;
			case 'A':
				anim = true;
				break;
			
			}
		} else {
			base_dir = argv[i] ;
		}
	}

	d6::Offline offline( base_dir ) ;

	if( width == 0 || height == 0 ) {
		const float a = offline.config().box[0]/offline.config().box[1] ;
		height = std::sqrt( 5.e5 / a ) ;
		width = height * a ;
	}

	d6::Appli appli( offline, frame, nSamples, width, height , sn, run,anim);
	if (discs)
		appli.viewer().grainsRenderer().useDiscs();
	if (colVel)
		appli.viewer().grainsRenderer().cutAndColorVelocities();
	appli.viewer().grainsRenderer().setGrainSizeFactor( grainSizeFactor );

	appli.init();
	return appli.run( ) ;

}

