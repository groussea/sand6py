#include "StressMassMatrix.hh"

#include "ActiveIndices.hh"
#include "MeshImpl.hh"

#include <Eigen/Eigenvalues>

#include <bogus/Core/Block.impl.hpp>


namespace d6 {

template < typename Shape>
void AbstractStressMassMatrix<Shape>::compute( const Shape& shape, const Active& nodes, const Index totNodes )
{
    inv_sqrt.setRows( totNodes ) ;
    inv_sqrt.setIdentity() ;
}

template < typename MeshT>
void AbstractStressMassMatrix<DGLinear<MeshT>>::compute( const Shape& shape, const Active& nodes, const Index totNodes )
{

    assert( 0 == ( nodes.count() % Shape::NI ) ) ;
    const Index nCells = nodes.count() / Shape::NI ;

    typedef Eigen::Matrix < Scalar, Shape::NI, Eigen::Dynamic > BlockAgg ;
    BlockAgg diagBlocks ;

    diagBlocks.resize( Shape::NI, Shape::NI*nCells );
    diagBlocks.setZero() ;

    typename Shape::Location loc ;
    typename Shape::Interpolation itp ;

    for( auto qpIt = shape.qpBegin() ; qpIt != shape.qpEnd() ; ++qpIt ) {
        qpIt.locate( loc ) ;
        shape.interpolate( loc, itp ) ;
        for( Index k = 0 ; k < itp.nodes.rows() ; ++k ) {
            const Index rowIndex = nodes.indices[itp.nodes[k]] ;
            if( rowIndex == Active::s_Inactive ) continue ;

            const Index cellIndex = rowIndex / Shape::NI ;
            const Index rowInCell = rowIndex - cellIndex * Shape::NI ;

            for( Index j = 0 ; j < itp.nodes.rows() ; ++j ) {
                const Index colIndex = nodes.indices[itp.nodes[j]] ;
                assert( cellIndex == (colIndex / Shape::NI ) ) ;
                if( colIndex == Active::s_Inactive ) continue ;

                diagBlocks( rowInCell, colIndex ) += itp.coeffs[k]*itp.coeffs[j] ;
            }
        }
    }

    typedef Eigen::Matrix< Scalar, Shape::NI, Shape::NI > MatDG ;

#pragma omp parallel for
    for( Index i = 0 ; i < nCells ; ++i) {
        typedef Eigen::Block< BlockAgg, Shape::NI, Shape::NI > BlockDG ;
        BlockDG B = diagBlocks.template block< Shape::NI, Shape::NI >( 0, i * Shape::NI  ) ;
        MatDG   M = B ;

        Eigen::SelfAdjointEigenSolver< MatDG > es( M ) ;
        B = es.eigenvectors() * (1. / es.eigenvalues().array().sqrt() ).matrix().asDiagonal() * es.eigenvectors().transpose() ;
//        B.setIdentity() ;

//        std::cout << (B*B*M) << std::endl ;
    }


    inv_sqrt.setRows( totNodes ) ;

    for( Index i = 0 ; i < nCells ; ++i ) {

        for( Index k = 0 ; k < Shape::NI ; ++k ) {
            for( Index j = 0 ; j <= k ; ++j ) {
                inv_sqrt.insertBack( i*Shape::NI + k, i*Shape::NI + j )
                        = diagBlocks( k, i*Shape::NI + j ) * MatS::Identity() ;
            }

        }
    }
    for( Index i = nodes.count() ; i < totNodes ; ++i ) {
        inv_sqrt.insertBack( i,i ).setIdentity() ;
    }

    inv_sqrt.finalize() ;

}


template struct AbstractStressMassMatrix< DualShape > ;

} //d6

