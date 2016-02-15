#ifndef D6_FORM_BUILDER_HH
#define D6_FORM_BUILDER_HH

#include "geo/geo.fwd.hh"
#include "geo/MeshBase.hh"

#include "utils/block_mat.hh"

// FIXME includes
#include "geo/MeshImpl.hh"

namespace d6 {

class Particles ;


//! Utility class for building matrices of bilinear forms
//! a(u,v) = u' A v, i.e. row == left, col == right
template < typename LhsShape, typename RhsShape >
class FormBuilder {

	typedef typename FormMat< WD,WD >::Type::MajorIndexType          CompressedIndexType ;
	typedef typename CompressedIndexType::Index BgIndex ;

	typedef const std::vector< Index > &Indices ;

	typedef const typename LhsShape::Interpolation& LhsItp ;
	typedef const typename LhsShape::Derivatives& LhsDcdx ;
	typedef const typename LhsShape::Derivatives::ConstRowXpr& LhsDcdxRow ;

	typedef const typename RhsShape::Interpolation& RhsItp ;
	typedef const typename RhsShape::Derivatives& RhsDcdx ;

public:

	FormBuilder( const LhsShape& lhsShape, const RhsShape& rhsShape )
		: m_lhsShape( lhsShape ), m_rhsShape( rhsShape )
	{}

	void reset( Index rows ) ;
	Index rows() const { return m_data.size() ; }

	//! Computes matrices non-zero blocks (nodes sharing a cell) from list of active cells
	template < typename CellIterator>
	void addToIndex(
			const CellIterator& cellBegin, const CellIterator& cellEnd,
			Indices rowIndices, Indices colIndices	 ) ;


	void addRows( Index rows ) ;

	//! Compute compressed sparse block matrix index
	void makeCompressed() ;


	//! Integrate over quadrature points
	template < typename QPIterator, typename Func >
	void integrate_qp( const QPIterator& qpBegin, const QPIterator& qpEnd, Func func ) const ;

	template < typename Func >
	void integrate_qp( Func func ) const ;
	template < typename CellIterator, typename Func >
	void integrate_cell( const CellIterator& cellBegin, const CellIterator& cellEnd, Func func ) const ;

	//! Integrate over nodes ( trapezoidal approx )
	template < typename CellIterator, typename Func >
	void integrate_node( const CellIterator& cellBegin, const CellIterator& cellEnd, Func func ) const ;

	//! Integrate over particles (MPM)
	template < typename Func >
	void integrate_particle( const Particles& particles, Func func ) const  ;

	// Building blocks

	static void addDuDv     ( FormMat<WD,WD>::Type& A, Scalar w, Index rowIndex, LhsDcdxRow row_dx,
							  RhsItp itp, RhsDcdx dc_dx, Indices colIndices ) ;
	static void addDpV      ( FormMat< 1,WD>::Type& A, Scalar w, Index colIndex, LhsItp itp, LhsDcdx dc_dx, Indices rowIndices ) ;
	static void addTauDu    ( FormMat<SD,WD>::Type& A, Scalar w, Index rowIndex, RhsItp itp, RhsDcdx dc_dx, Indices colIndices ) ;
	static void addTauWu    ( FormMat<RD,WD>::Type& A, Scalar w, Index rowIndex, RhsItp itp, RhsDcdx dc_dx, Indices colIndices ) ;

	static void addDuDv     ( FormMat<WD,WD>::Type& A, Scalar w, LhsItp lhs_itp, LhsDcdx lhs_dc_dx, RhsItp rhs_itp, RhsDcdx rhs_dc_dx, Indices rowIndices, Indices colIndices ) ;
	static void addDpV      ( FormMat<1, WD>::Type& A, Scalar w, LhsItp lhs_itp, LhsDcdx dc_dx, RhsItp rhs_itp, Indices rowIndices, Indices colIndices ) ;
	static void addTauDu    ( FormMat<SD,WD>::Type& A, Scalar w, LhsItp lhs_itp, RhsItp rhs_itp, RhsDcdx dc_dx, Indices rowIndices, Indices colIndices ) ;
	static void addTauWu    ( FormMat<RD,WD>::Type& A, Scalar w, LhsItp lhs_itp, RhsItp rhs_itp, RhsDcdx dc_dx, Indices rowIndices, Indices colIndices ) ;

	static void addUTaunGphi( FormMat<SD,WD>::Type& A, Scalar w, LhsItp lhs_itp, RhsItp rhs_itp, const Vec& dphi_dx, Indices rowIndices, Indices colIndices ) ;
	static void addUTauGphi ( FormMat<SD,WD>::Type& A, Scalar w, LhsItp lhs_itp, RhsItp rhs_itp, const Vec& dphi_dx, Indices rowIndices, Indices colIndices ) ;

	const CompressedIndexType& index() { return m_compressed ; }


	//FIXME make private
	LhsShape m_lhsShape ;
	RhsShape m_rhsShape ;
private:

	CompressedIndexType m_compressed ;

	std::vector< std::vector< BgIndex > > m_data ;

};


} //d6

#endif
