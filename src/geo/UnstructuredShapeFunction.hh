#ifndef D6_UNSTRUCTURED_SHAPE_FUNCTION
#define D6_UNSTRUCTURED_SHAPE_FUNCTION

#include "ShapeFunctionBase.hh"

namespace d6 {

struct UnstructuredShapeFunc ;
struct Config ;

struct UnstructuredDOFs {
	typedef Eigen::Matrix< Scalar, WD, Eigen::Dynamic > Vertices ;
	typedef Eigen::Matrix< Scalar,  1, Eigen::Dynamic > Weights ;

	const Vertices& vertices ;
	Weights weights ;

	explicit UnstructuredDOFs( const Vertices& v )
		: vertices(v), m_count(0)
	{}

	Index count() const { return m_count ; }

	void resize( Index n ) { m_count = n ; }
	void compute_weights_from_vertices( const Config &c) ;

	template <typename Ar>
	void serialize( Ar& ar, const unsigned int ) { ar&m_count ; }

private:
	Index m_count ;
};

struct UnstructuredDOFIterator
{

	explicit UnstructuredDOFIterator( const UnstructuredDOFs& def, const Index i )
		: dofDef(def), index(i)
	{}

	bool operator==( const UnstructuredDOFIterator& o ) const
	{ return o.index == index ;}
	bool operator!=( const UnstructuredDOFIterator& o ) const
	{ return o.index != index ; }

	UnstructuredDOFIterator& operator++() {
		++index ;
		return *this ;
	}

	Vec pos() const {
		return dofDef.vertices.col( index ) ;
	}
	Scalar weight() const {
		return dofDef.weights[index] ;
	}
	void locate( Index &loc ) const {
		loc = index ;
	}

private:
	UnstructuredDOFs dofDef ;
	Index index ;
};

template<>
struct ShapeFuncTraits< UnstructuredShapeFunc >
{

	typedef Index Location ;
	typedef UnstructuredDOFs DOFDefinition ;

	enum{ NI = 1, NQ = 1 } ;

	static bool constexpr is_mesh_based = false ;

	template <typename CellIterator = Index >
	struct QPIterator {
		typedef UnstructuredDOFIterator Type ;
	};
} ;

struct UnstructuredShapeFunc : public ShapeFuncBase< UnstructuredShapeFunc >
{

	typedef ShapeFuncBase< UnstructuredShapeFunc > Base ;
	typedef ShapeFuncTraits< UnstructuredShapeFunc > Traits ;
	typedef typename Traits::DOFDefinition DOFDefinition ;

	UnstructuredShapeFunc( const DOFDefinition &v ) : m_dofDef(v)
	{}

	void compute_lumped_mass( DynVec& volumes ) const
	{
		volumes = m_dofDef.weights ;
	}

	Index nDOF() const { return m_dofDef.count() ; }

	void locate_by_pos_or_id( const Vec&, const Index id, typename Base::Location & loc ) const {
		loc = id ;
	}

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

	UnstructuredDOFIterator qpBegin() const {
		return UnstructuredDOFIterator( m_dofDef, 0 ) ;
	}

	UnstructuredDOFIterator qpEnd() const {
		return UnstructuredDOFIterator( m_dofDef, nDOF() ) ;
	}

	template <typename CellIterator>
	UnstructuredDOFIterator qpIterator( const Index &it ) const {
		return UnstructuredDOFIterator( m_dofDef, it ) ;
	}

protected:
	const DOFDefinition &m_dofDef ;
};

} //d6

#endif
