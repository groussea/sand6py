#ifndef D6_FORM_BUILDER_HH
#define D6_FORM_BUILDER_HH

#include "geo/geo.fwd.hh"
#include "utils/alg.hh"

#include "geo/MeshBase.hh"

#include <bogus/Core/Block.hpp>

namespace d6 {

class Particles ;

template < Index Rows, Index Cols >
struct FormMat {
	typedef Eigen::Matrix< Scalar, Rows, Cols > BlockT ;
	typedef bogus::SparseBlockMatrix< BlockT > Type ;
	typedef bogus::SparseBlockMatrix< BlockT, bogus::SYMMETRIC > SymType ;
};


class FormBuilder {

	typedef typename FormMat< 3,3 >::Type::MajorIndexType          CompressedIndexType ;
	typedef typename CompressedIndexType::Index BgIndex ;

public:

	FormBuilder( const MeshType& mesh )
		: m_mesh(mesh)
	{}

	void reset( Index rows ) { m_data.clear(); m_data.resize( rows ); }

	void addToIndex(
			const typename MeshType::Cells& cells,
			const std::vector< Index > &rowIndices,
			const std::vector< Index > &colIndices
					 ) ;

	void makeCompressed() ;


	template < typename Func >
	void integrate_qp( const typename MeshType::Cells& cells, Func func	) const ;

	template < typename Func >
	void integrate_particle( const Particles& particles, Func func ) const  ;

	const CompressedIndexType& index() { return m_compressed ; }

private:
	const MeshType &m_mesh ;

	CompressedIndexType m_compressed ;

	std::vector< std::vector< BgIndex > > m_data ;

};


} //d6

#endif
