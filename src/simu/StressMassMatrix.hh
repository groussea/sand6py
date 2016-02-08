#ifndef D6_STRESS_MASS_MATRIX_HH
#define D6_STRESS_MASS_MATRIX_HH

#include "PhaseFields.hh"

#include "geo/MeshShapeFunction.hh"

#include "utils/block_mat.hh"

namespace d6 {

struct Active ;

template< typename Shape >
struct AbstractStressMassMatrix
{
    static constexpr bool use_identity = true ;
    typename FormMat<SD, SD>::SymType  inv_sqrt ;

    void compute( const Shape& shape , const Active& nodes, const Index totNodes )  ;
};

template < typename MeshT>
struct AbstractStressMassMatrix<DGLinear<MeshT>>
{
    typedef DGLinear<MeshT> Shape ;

    static constexpr bool use_identity = false ;
    typename FormMat<SD, SD>::SymType  inv_sqrt ;

    void compute( const Shape& shape, const Active& nodes, const Index totNodes )  ;

} ;


typedef AbstractStressMassMatrix< DualShape > StressMassMatrix ;

} //d6

#endif
