#ifndef D6_EXPR_HH
#define D6_EXPR_HH

#include "utils/alg.hh"
#include "geo.fwd.hh"

namespace d6 {

template< typename ValueType >
struct Expr {
	virtual ValueType operator() ( const Vec&  x ) const = 0 ;
};

typedef Expr< Scalar > ScalarExpr ;

template < typename Functor, typename ValueType = typename std::result_of< Functor( const Vec& ) >::type >
struct LambdaExpr : public Expr< ValueType > {
	const Functor& expr ;

	explicit LambdaExpr( const Functor& xpr ) : expr(xpr) {}
	ValueType operator() ( const Vec&  x ) const  { return expr(x) ; }
};

template < typename Functor >
LambdaExpr< Functor > make_expr( const Functor& expr )
{
	return LambdaExpr< Functor >{ expr } ;
}


} //d6

#endif
