/*
 * This file is part of Sand6, a C++ continuum-based granular simulator.
 *
 * Copyright 2016 Gilles Daviet <gilles.daviet@inria.fr> (Inria - Université Grenoble Alpes)
 *
 * Sand6 is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.

 * Sand6 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.

 * You should have received a copy of the GNU General Public License
 * along with Sand6.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef D6_LEVEL_SET_IMPL_HH
#define D6_LEVEL_SET_IMPL_HH

#include "LevelSet.hh"

#include "Grid.hh"
#include "ScalarField.hh"

#include <limits>

namespace d6 {

struct SphereLevelSet : public LevelSet
{
	Scalar eval_local(const Vec &x) const override {
		return 1. - x.norm() ;
	}

	Vec grad_local(const Vec &x) const override {
		return -x / ( 1.e-12 + x.norm() ) ;
	}

	void local_inv_inertia( MatR& I ) const override {
		I(0,0) = 2. / local_volume() ;
	}

	Scalar local_volume() const override {
		return M_PI ;
	}

	template<class Archive>
	void serialize(Archive &ar, const unsigned int version ) ;
};
struct PlaneLevelSet : public LevelSet
{
	Scalar eval_local(const Vec &x) const override {
		return - x[1] ;
	}

	Vec grad_local(const Vec & ) const override {
		return Vec(0, -1) ;
	}

	void local_inv_inertia( MatR& I ) const override {
		I.setZero() ;
	}

	Scalar local_volume() const override {
		return std::numeric_limits<Scalar>::infinity() ;
	}

	template<class Archive>
	void serialize(Archive &ar, const unsigned int version ) ;
};
struct SegLevelSet : public LevelSet
{
	SegLevelSet( const Scalar len = 1)
		: m_len(len)
	{}

	Scalar eval_local(const Vec &x) const override {
		Vec p ;
		proj_on_seg( x, p );
		return 1 - (x-p).norm() ;
	}

	Vec grad_local(const Vec &x ) const override {
		Vec p ;
		proj_on_seg( x, p );
		const Scalar n = (x-p).norm() ;
		return (p-x)/(n+1.e-12) ;
	}

	void local_inv_inertia( MatR& I ) const override {
		I.setZero() ;
	}

	Scalar local_volume() const override {
		return std::numeric_limits<Scalar>::infinity() ;
	}

	template<class Archive>
	void serialize(Archive &ar, const unsigned int version ) ;

	Scalar len() const { return m_len ; }

private:

	void proj_on_seg( const Vec& x, Vec &p ) const {
		p[1] = std::max( -.5*m_len, std::min( .5*m_len, x[1] ) ) ;
		p[0] = 0 ;
	}

	Scalar m_len ;
};


} //ns d6

#endif
