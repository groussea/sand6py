#ifndef D6_FORM_BUILDER_HH
#define D6_FORM_BUILDER_HH

#include "geo/geo.fwd.hh"
#include "geo/MeshBase.hh"

#include "utils/block_mat.hh"

namespace d6 {

class Particles ;


class FormBuilder {

	typedef typename FormMat< 3,3 >::Type::MajorIndexType          CompressedIndexType ;
	typedef typename CompressedIndexType::Index BgIndex ;

	typedef const std::vector< Index > &Indices ;
	typedef const typename MeshType::Interpolation& Itp ;
	typedef const typename MeshType::Derivatives& Dcdx ;


public:

	FormBuilder( const MeshType& mesh )
		: m_mesh(mesh)
	{}

	void reset( Index rows ) ;

	void addToIndex(
			const typename MeshType::Cells& cells,
			Indices rowIndices, Indices colIndices	 ) ;

	void addRows( Index rows ) ;

	void makeCompressed() ;


	template < typename Func >
	void integrate_qp( const typename MeshType::Cells& cells, Func func	) const ;

	template < typename Func >
	void integrate_particle( const Particles& particles, Func func ) const  ;

	static void addDuDv ( FormMat<3,3>::Type& A, Scalar w, Itp itp, Dcdx dc_dx, Indices rowIndices, Indices colIndices ) ;
	static void addVDsig( FormMat<3,6>::Type& A, Scalar w, Itp itp, Dcdx dc_dx, Indices rowIndices, Indices colIndices ) ;

	static void addTauDu( FormMat<6,3>::Type& A, Scalar w, Itp itp, Dcdx dc_dx, Indices rowIndices, Indices colIndices ) ;
	static void addTauWu( FormMat<3,3>::Type& A, Scalar w, Itp itp, Dcdx dc_dx, Indices rowIndices, Indices colIndices ) ;

	static void addUTauGphi( FormMat<6,3>::Type& A, Scalar w, Itp itp, const Vec& dphi_dx, Indices rowIndices, Indices colIndices ) ;

	const CompressedIndexType& index() { return m_compressed ; }

private:

	const MeshType &m_mesh ;

	CompressedIndexType m_compressed ;

	std::vector< std::vector< BgIndex > > m_data ;

};


} //d6

#endif
