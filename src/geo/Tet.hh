#ifndef D6_TET_HH
#define D6_TET_HH

#include "utils/alg.hh"

namespace d6 {

struct Tet {

	static constexpr Index NV = 4 ;
	static constexpr Index NC = 4 ;
	static constexpr Index NQ = 4 ;

	typedef Eigen::Matrix< Scalar, NC, 1 > Coords ;
	typedef Eigen::Matrix< Scalar, NV, 3 > Derivatives ;
	
	typedef Eigen::Matrix< Scalar, NC, NQ> QuadPoints ;
	typedef Eigen::Matrix< Scalar,  1, NQ> QuadWeights ;

	typedef Eigen::Matrix< Scalar, 3, 4 > Vertices ;
	typedef Eigen::Matrix< Scalar, 3, Eigen::Dynamic > Points ;
	typedef Eigen::Matrix< Scalar, 6, Eigen::Dynamic > Frames ;

	struct {
		unsigned char part : 1 ;
		unsigned char rot  : 3 ;
		unsigned char sym  : 3 ;

		inline char rotate ( int d ) const { return rot&(1<<d) ; }
		inline char sign   ( int d ) const { 
			return 1 - 2*( ( sym&(1<<d) ) >> d )  ; 
		}

	} geometry ;

	enum Part {
		Left  = 0,
		Right = 1
	};

	Vec origin ;  //!< 3D coords of first corner
	Arr box    ;  //!< Dimensions of cell

	int orientation ; //!< 1 bit part + 3 bits symmetry

	Tet()
	: geometry{0,0,0}
	{}

	Vec pos( const Coords& coords ) const {

		Vertices vtx  ;
		compute_vertices( vtx ) ;

		Vec p = vtx*coords ;

		to_world( p ) ;
		return p + origin ;
	}

	Vec center() const
	{
		//TODO optimize
		return pos( Coords::Constant(1./NC) ) ;
	}

	Vec vertex( int cornerIndex ) const
	{
		Vec v ;
		offset( cornerIndex, v ) ;

		to_world(v) ;
		return v + origin ;
	}

	Scalar volume() const {
		return box.prod() / 6 ;
	}

	void compute_coords( const Vec& pos, Coords & coords ) const {
		Vec local = pos - origin ;
		to_local( local );
		local_coords( local, coords );
	}
	void compute_derivatives( const Coords & coords, Derivatives& dc_dx ) const ;

	Index sample_uniform( const unsigned N, const Index start, Points &points, Frames &frames ) const ;

	static QuadPoints Qps() ;

	void get_qp( QuadPoints& qp, QuadWeights& weights ) const {
		static const QuadPoints s_qps = Qps() ;
		qp = s_qps ;
		weights.setConstant( volume() / NQ ) ;
	}

	void update_geometry( unsigned char rotation, int num ) ;

private:

	void compute_vertices( Vertices& vertices ) const ;

	void offset( int cornerIndex, Vec &v ) const ;

	void to_world( Vec& pos ) const ;
	void to_local( Vec& pos ) const ;

	void local_coords( const Vec& pos, Coords& coords ) const ;

} ;

} //d6

#endif
