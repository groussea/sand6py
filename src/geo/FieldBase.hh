#ifndef D6_FIELD_BASE_HH
#define D6_FIELD_BASE_HH

#include "MeshBase.hh"

namespace d6 {

template< typename Derived >
struct FieldTraits
{

};

template< typename Derived >
struct FieldBase
{

public:

	typedef FieldTraits< Derived > Traits ;
	typedef MeshBase< typename Traits::MeshType > MeshType ;

	static constexpr unsigned D = Traits::Dimension ;
	typedef typename Traits::ValueType ValueType ;
	typedef Eigen::Matrix< Scalar, D, 1 > Seg ;

	FieldBase( const MeshType& mesh )
		: m_mesh( mesh )
	{
		m_data.resize( D * m_mesh.nNodes() );
	}

	const DynVec& flatten() const { return m_data ; }
	DynVec& flatten() { return m_data ; }

	void eval_at( const Vec& x, Seg& res ) const ;

	ValueType eval_at( const Vec& x ) const {
		Seg seg ; eval_at( x, seg );
		return seg ;
	}

protected:
	const MeshType & m_mesh ;
	DynVec   m_data ;

} ;

} //d6

#endif

