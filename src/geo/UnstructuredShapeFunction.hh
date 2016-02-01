#ifndef D6_UNSTRUCTURED_SHAPE_FUNCTION
#define D6_UNSTRUCTURED_SHAPE_FUNCTION

#include "ShapeFunctionBase.hh"

namespace d6 {

struct UnstructuredShapeFunc ;

template<>
struct ShapeFuncTraits< UnstructuredShapeFunc >
{
    typedef DynMatW Vertices ;
    typedef Index Location ;
    typedef Vertices DOFDefinition ;

    enum{ NI = 1, NQ = 1 } ;

    static bool constexpr is_mesh_based = false ;
} ;

struct UnstructuredShapeFunc : public ShapeFuncBase< UnstructuredShapeFunc >
{
    typedef ShapeFuncBase< UnstructuredShapeFunc > Base ;
    typedef ShapeFuncTraits< UnstructuredShapeFunc > Traits ;
    typedef typename Traits::Vertices Vertices ;

    UnstructuredShapeFunc( Vertices &v ) : vertices(v)
    {}

    void compute_volumes( DynVec& volumes ) const
    {
        // FIXME
        // TODO compute half-distance to nearest-neighbour
        volumes.setZero( nDoF() ) ;
    }

    Index nDoF() const { return vertices.cols() ; }

    void interpolate( const Location& loc, typename Base::Interpolation& itp ) const {
        itp.nodes[0] = loc ;
        itp.coeffs[0] = 1 ;
    }
    void get_derivatives( const Location&, typename Base::Derivatives& dc_dx ) const {
        dc_dx.setZero() ;
    }

    void list_nodes( const Location& loc, typename Base::NodeList& list ) const {
        list[0] = loc ;
    }

    void all_dof_positions( DynMatW& v  ) const
    {
        v = vertices ;
    }


    const Vertices& vertices ;
};

} //d6

#endif

