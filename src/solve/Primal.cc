#include "Primal.hh"

namespace d6 {

bool PrimalData::load(const char *file)
{
	return false ;
}

bool PrimalData::dump(const char *file) const
{
	return false ;
}


Primal::Primal(const PrimalData &data)
	: m_data( data )
{}


Scalar Primal::solve(DynVec &lambda, DynVec &gamma) const
{
	lambda.setZero() ;
	gamma.setZero() ;

	return 0. ;
}


}
