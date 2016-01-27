
#include "GLViewer.hh"

#include "utils/string.hh"
#include "utils/Config.hh"

#include "visu/Offline.hh"

#ifdef GLFW3
#include <GLFW/glfw3.h>
#else
#include <GL/glfw.h>
#endif

#include <cstdlib>

namespace d6
{


class Appli
{

public:
	Appli( d6::Offline& offline, unsigned frame,
		   const int width, const int height ) :
	#ifdef GLFW3
		m_pWindow( NULL ),
	#else
		m_mouseWheelPos(0),
	#endif
		m_offline(offline), m_viewer( offline, width, height ),
		m_currentFrame( frame ),
		m_running( false ), m_snapshotting( false ),
		m_mouseX( 0 ), m_mouseY( 0 )
	{
		assert( ! s_instance ) ;
		s_instance = this ;

		init() ;
	}

	~Appli()
	{
		quit(0) ;
	}

	void set_frame( unsigned frame ) {
		if( m_offline.load_frame( frame ) ) {
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

		do {
			if( m_running ) {
				next_frame() ;
			}

			m_viewer.draw( ) ;

			if( m_running && m_snapshotting ) {
				m_viewer.snap();
			}

#ifdef GLFW3
			glfwSwapBuffers( m_pWindow);
			glfwPollEvents();
		}
		while( !glfwWindowShouldClose(m_pWindow)
			   && ( glfwGetKey(m_pWindow, GLFW_KEY_ESCAPE ) != GLFW_PRESS )
			   && ( glfwGetKey(m_pWindow, GLFW_KEY_Q) != GLFW_PRESS )) ;
#else
			glfwSwapBuffers();
			glfwPollEvents();

		}
		while( glfwGetWindowParam( GLFW_OPENED )
			   && ( glfwGetKey(GLFW_KEY_ESC) != GLFW_PRESS )
			   && ( glfwGetKey('Q') != GLFW_PRESS )) ;
#endif
		return 0 ;
	}

	static Appli& instance()
	{
		return *s_instance ;
	}

#ifdef GLFW3
	static void key_callback(GLFWwindow *, int key, int /*scancode*/, int action, int /*mods*/ )
#else
	static void key_callback( int key, int action )
#endif
	{
		instance().process_key( key, action ) ;
	}

#ifdef GLFW3
	static void mouse_click_callback(GLFWwindow *, int button, int action, int /*mods*/ )
#else
	static void mouse_click_callback(int button, int action )
#endif
	{
		instance().process_mouse_click( button, action );
	}

#ifdef GLFW3
	static void mouse_motion_callback(GLFWwindow *, double x, double y )
#else
	static void mouse_motion_callback( int x, int y )
#endif
	{
		instance().process_mouse_motion( x, y );
	}

#ifdef GLFW3
	static void mouse_scroll_callback(GLFWwindow *, double, double y_offset )
#else
	static void mouse_scroll_callback( int y_offset )
#endif
	{
		instance().process_mouse_scroll( y_offset );
	}


private:

	void quit(int rc)
	{
#ifdef GLFW3
		if( m_pWindow) glfwDestroyWindow(m_pWindow);
#endif
		std::cout << "Bye. [code " << rc << "]" << std::endl ;
		glfwTerminate();
		std::exit(rc);
	}

	void init()
	{
		if ( !glfwInit() )
			quit(1);
#ifdef GLFW3
		if(!( m_pWindow = glfwCreateWindow(m_viewer.width(), m_viewer.height(), "Hyb2D", NULL, NULL) ))
#else
		if (glfwOpenWindow(m_viewer.width(), m_viewer.height(), 5, 6, 5,
						   0, 0, 0, GLFW_WINDOW) != GL_TRUE)
#endif
		{
			quit(1);
		}

#ifdef GLFW3
		glfwMakeContextCurrent(m_pWindow);

		glfwSetInputMode( m_pWindow, GLFW_STICKY_KEYS, GL_TRUE ) ;
		glfwSetKeyCallback( m_pWindow, &Appli::key_callback ) ;
		glfwSetMouseButtonCallback( m_pWindow, &Appli::mouse_click_callback ) ;
		glfwSetScrollCallback( m_pWindow, &Appli::mouse_scroll_callback ) ;
#else
		glfwSetWindowTitle("Hyb2D");
		glfwEnable( GLFW_STICKY_KEYS );
		glfwSetKeyCallback( &Appli::key_callback ) ;
		glfwSetMouseButtonCallback( &Appli::mouse_click_callback ) ;
		glfwSetMouseWheelCallback( &Appli::mouse_scroll_callback ) ;
#endif

		std::cout << "Using OpenGL " << glGetString(GL_VERSION)
				  << " rendered on " << glGetString(GL_RENDERER)
				  << " from " << glGetString(GL_VENDOR)
				  << std::endl  ;

		m_viewer.init( ) ;
		set_frame( m_currentFrame );

	}

	void process_key( int key, int action )
	{
		if( action != GLFW_RELEASE )
			return ;

		switch( key )
		{
		case GLFW_KEY_ENTER:
			m_running = !m_running ;
			break;
		case 'A':
			m_viewer.toggleScalarEntity();
			m_viewer.update_color_buffers();
			break;
		case 'C':
			m_viewer.toggleRendering( GLViewer::eColors );
			m_viewer.update_color_buffers();
			break;
		case 'D':
			m_viewer.toggleTensorEntity();
			m_viewer.update_tensor_buffers();
			break;
		case 'F':
			m_viewer.toggleVectorEntity();
			m_viewer.update_vector_buffers();
			break;
		case 'G':
			m_viewer.toggleRendering( GLViewer::eGrid );
			break;
		case 'I':
			if( ! m_running) next_frame();
			break ;
		case 'M':
			m_viewer.toggleRendering( GLViewer::eParticles );
			m_viewer.update_particle_buffers();
			break;
		case 'O':
			m_viewer.toggleParticleRepr();
			m_viewer.update_particle_buffers();
			break;
		case 'P':
			if( ! m_running) prev_frame();
			break;
		case 'R':
			m_snapshotting = ! m_snapshotting ;
			break;
		case 'S':
			m_viewer.toggleRendering( GLViewer::eStream );
			m_viewer.update_texture();
			break;
		case 'T':
			m_viewer.toggleRendering( GLViewer::eTensors );
			m_viewer.update_tensor_buffers();
			break;
		case 'V':
			m_viewer.toggleRendering( GLViewer::eVectors );
			m_viewer.update_vector_buffers();
			break;
		case GLFW_KEY_HOME:
			set_frame(0);
			break;
#ifdef GLFW3
		case GLFW_KEY_PERIOD:
#else
		case '.':
#endif
			m_viewer.snap() ;
			break ;
		}

	}

	void process_mouse_motion( double x, double y )
	{
		m_viewer.move( x - m_mouseX, y - m_mouseY );
		m_mouseX = x ; m_mouseY = y ;
	}

	void process_mouse_click( int button, int action )
	{
		if( button != GLFW_MOUSE_BUTTON_1 )
			return ;

#ifdef GLFW3
		double x, y ;
		glfwGetCursorPos( m_pWindow, &x, &y ) ;
		m_mouseX = x ; m_mouseY = y ;
		if( action == GLFW_PRESS )
		{
			glfwSetCursorPosCallback( m_pWindow, &Appli::mouse_motion_callback ) ;
		} else {
			glfwSetCursorPosCallback( m_pWindow, NULL ) ;
		}
#else
		int x, y ;
		glfwGetMousePos( &x, &y ) ;
		m_mouseX = x ; m_mouseY = y ;
		if( action == GLFW_PRESS )
		{
			glfwSetMousePosCallback( &Appli::mouse_motion_callback ) ;
		} else {
			glfwSetMousePosCallback( NULL ) ;
		}
#endif

	}

	void process_mouse_scroll( double offset )
	{
#ifndef GLFW3
		offset -= m_mouseWheelPos ;
		m_mouseWheelPos += offset ;
#endif

		if( offset > 0 )
			m_viewer.zoom( .8 ) ;
		else if( offset < 0 )
			m_viewer.zoom( 1.25 ) ;
	}


#ifdef GLFW3
	GLFWwindow* m_pWindow ;
#else
	int m_mouseWheelPos ;
#endif

	Offline& m_offline ;
	GLViewer m_viewer ;

	unsigned m_currentFrame ;
	bool m_running ;
	bool m_snapshotting ;

	static Appli* s_instance ;

	double m_mouseX, m_mouseY ;
};

Appli* Appli::s_instance = NULL ;

} //d6

static void usage( const char *name )
{
	std::cout << "Usage: " << name
			  << " [sim_dir=out] [options] "
			  << "\n2D OpenGL viewer."
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

	unsigned width  = 1280 ;
	unsigned height = 0 ;

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
			}
		} else {
			base_dir = argv[i] ;
		}
	}

	d6::Offline offline( base_dir ) ;

	if( height == 0 ) {
		height = width * (offline.config().box[1]/offline.config().box[0]) ;
	}

	return d6::Appli( offline, frame, width, height ).run( ) ;

}

