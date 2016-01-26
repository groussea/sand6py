#ifndef D6_FORM_BUILDER_HH
#define D6_FORM_BUILDER_HH

#include "geo/geo.fwd.hh"
#include "geo/MeshBase.hh"

#include "utils/block_mat.hh"

namespace d6 {

class Particles ;

//! Utility class for building matrices of bilinear forms
class FormBuilder {

	typedef typename FormMat< WD,WD >::Type::MajorIndexType          CompressedIndexType ;
	typedef typename CompressedIndexType::Index BgIndex ;

	typedef const std::vector< Index > &Indices ;
	typedef const typename MeshType::Interpolation& Itp ;
	typedef const typename MeshType::Derivatives& Dcdx ;


public:

	FormBuilder( const MeshType& mesh )
		: m_mesh(mesh)
	{}

	void reset( Index rows ) ;
	Index rows() const { return m_data.size() ; }

	//! Computes matrices non-zero blocks (nodes sharing a cell) from list of active cells
	void addToIndex(
			const typename MeshType::Cells& cells,
			Indices rowIndices, Indices colIndices	 ) ;

	void addRows( Index rows ) ;

	//! Compute compressed sparse block matrix index
	void makeCompressed() ;


	//! Integrate over quadrature points
	template < typename Func >
	void integrate_qp( const typename MeshType::Cells& cells, Func func	) const ;
	template < typename CellIterator, typename Func >
	void integrate_qp( const CellIterator& cellBegin, const CellIterator& cellEnd, Func func ) const ;

	//! Integrate over nodes ( trapezoidal approx )
	template < typename Func >
	void integrate_node( const typename MeshType::Cells& cells, Func func	) const ;
	template < typename CellIterator, typename Func >
	void integrate_node( const CellIterator& cellBegin, const CellIterator& cellEnd, Func func ) const ;

	//! Integrate over particles (MPM)
	template < typename Func >
	void integrate_particle( const Particles& particles, Func func ) const  ;

	// Building blocks

	static void addDuDv     ( FormMat<WD,WD>::Type& A, Scalar w,
							  Index rowIndex, const typename MeshType::Derivatives::ConstRowXpr& row_dx,
							   Itp itp, Dcdx dc_dx, Indices colIndices ) ;

	static void addVDp      ( FormMat<WD, 1>::Type& A, Scalar w, Index rowIndex, Itp itp, Dcdx dc_dx, Indices colIndices ) ;
	static void addTauDu    ( FormMat<SD,WD>::Type& A, Scalar w, Index rowIndex, Itp itp, Dcdx dc_dx, Indices colIndices ) ;
	static void addTauWu    ( FormMat<RD,WD>::Type& A, Scalar w, Index rowIndex, Itp itp, Dcdx dc_dx, Indices colIndices ) ;

	static void addDuDv     ( FormMat<WD,WD>::Type& A, Scalar w, Itp itp, Dcdx dc_dx, Indices rowIndices, Indices colIndices ) ;
	static void addVDp      ( FormMat<WD, 1>::Type& A, Scalar w, Itp itp, Dcdx dc_dx, Indices rowIndices, Indices colIndices ) ;
	static void addTauDu    ( FormMat<SD,WD>::Type& A, Scalar w, Itp itp, Dcdx dc_dx, Indices rowIndices, Indices colIndices ) ;
	static void addTauWu    ( FormMat<RD,WD>::Type& A, Scalar w, Itp itp, Dcdx dc_dx, Indices rowIndices, Indices colIndices ) ;

	static void addUTaunGphi( FormMat<SD,WD>::Type& A, Scalar w, Itp itp, const Vec& dphi_dx, Indices rowIndices, Indices colIndices ) ;
	static void addUTauGphi ( FormMat<SD,WD>::Type& A, Scalar w, Itp itp, const Vec& dphi_dx, Indices rowIndices, Indices colIndices ) ;

	const CompressedIndexType& index() { return m_compressed ; }

private:

	const MeshType &m_mesh ;

	CompressedIndexType m_compressed ;

	std::vector< std::vector< BgIndex > > m_data ;

};


} //d6

#endif
