#include "ActiveIndices.hh"

#include "PhaseFields.hh"

#include "geo/ScalarField.hh"
#include "geo/VectorField.hh"
#include "geo/TensorField.hh"

#include "geo/MeshImpl.hh"

namespace d6 {

const Index Active::s_Inactive = -1 ;

void Active::computeRevIndices()
{
	revIndices.resize( nNodes );

#pragma omp parallel for
	for( size_t i = 0 ; i < indices.size() ; ++i  ) {
		const Index idx = indices[ i ] ;
		if( idx != Active::s_Inactive ) {
			revIndices[idx-offset] = i ;
		}
	}
}

void Active::setOffset(const Index o)
{
	offset = o ;

#pragma omp parallel for
	for( size_t i = 0 ; i < indices.size() ; ++i  ) {
		if( indices[i] != Active::s_Inactive )
			indices[i] += o ;
	}
}

template < typename Derived >
void Active::field2var( const FieldBase<Derived> &field, DynVec& var, bool resize ) const
{
	constexpr Index D = FieldBase<Derived>::D ;

	if( resize )
		var.resize( (offset + D) * count() );

#pragma omp parallel for
	for( Index i = 0 ; i < nNodes ; ++ i) {
		const Index idx = revIndices[ i ] ;
		Segmenter<D>::segment( var, offset+i ) = field[ idx ] ;
	}
}

template < typename Derived >
void Active::field2var( const FieldBase<Derived> &field, DynArr& var, bool resize ) const
{
	constexpr Index D = FieldBase<Derived>::D ;

	if( resize )
		var.resize( (offset + D) * count() );

#pragma omp parallel for
	for( Index i = 0 ; i < nNodes ; ++ i) {
		const Index idx = revIndices[ i ] ;
		Segmenter<D, DynArr>::segment( var, offset+i ) = field[ idx ] ;
	}
}

template < typename Derived >
void Active::var2field( const DynVec& var,  FieldBase<Derived> &field ) const
{
	constexpr Index D = FieldBase<Derived>::D ;

	field.set_zero();

#pragma omp parallel for
	for( Index i = 0 ; i < nNodes ; ++ i) {
		const Index idx = revIndices[ i ] ;
		field[ idx ] = Segmenter<D>::segment( var, offset+i ) ;
	}
}


template void Active::var2field( const DynVec & var, FieldBase<AbstractScalarField<Linear<MeshImpl>>> &field ) const;
template void Active::var2field( const DynVec & var, FieldBase<AbstractVectorField<Linear<MeshImpl>>> &field ) const;
template void Active::var2field( const DynVec & var, FieldBase<AbstractTensorField<Linear<MeshImpl>>> &field ) const;
template void Active::var2field( const DynVec & var, FieldBase<AbstractSkewTsField<Linear<MeshImpl>>> &field ) const;
template void Active::var2field( const DynVec & var, FieldBase<AbstractScalarField<DGLinear<MeshImpl>>> &field ) const;
template void Active::var2field( const DynVec & var, FieldBase<AbstractTensorField<DGLinear<MeshImpl>>> &field ) const;
template void Active::var2field( const DynVec & var, FieldBase<AbstractSkewTsField<DGLinear<MeshImpl>>> &field ) const;

template void Active::field2var( const FieldBase<AbstractScalarField<Linear<MeshImpl>>> &field, DynArr & var, bool ) const;
template void Active::field2var( const FieldBase<AbstractScalarField<Linear<MeshImpl>>> &field, DynVec & var, bool ) const;
template void Active::field2var( const FieldBase<AbstractVectorField<Linear<MeshImpl>>> &field, DynVec & var, bool ) const;
template void Active::field2var( const FieldBase<AbstractTensorField<Linear<MeshImpl>>> &field, DynVec & var, bool ) const;
template void Active::field2var( const FieldBase<AbstractSkewTsField<Linear<MeshImpl>>> &field, DynVec & var, bool ) const;

template void Active::field2var( const FieldBase<AbstractScalarField<DGLinear<MeshImpl>>> &field, DynArr & var, bool ) const;
template void Active::field2var( const FieldBase<AbstractTensorField<DGLinear<MeshImpl>>> &field, DynVec & var, bool ) const;
template void Active::field2var( const FieldBase<AbstractSkewTsField<DGLinear<MeshImpl>>> &field, DynVec & var, bool ) const;

} //d6
