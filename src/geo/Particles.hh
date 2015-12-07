#ifndef D6_PARTICLES_HH
#define D6_PARTICLES_HH

#include "Expr.hh"
#include "utils/Mutex.hh"

namespace d6 {

class Particles {

public:

	template< int D >
	struct Data {
		typedef Eigen::Matrix< Scalar, D, Eigen::Dynamic > Type ;
	};

	struct Event {
		enum Type { Split, Merge, Remove } ;
		Type   type   ;
		size_t first  ;
		size_t second ;

		Vec    dx ;

		template < typename Archive >
		void serialize( Archive &ar, unsigned int ) ;

		static Event split( size_t first, size_t second, Vec dx ) {
			return Event{ Split, first, second, dx } ;
		}
		static Event merge( size_t first, size_t second, Vec dx ) {
			return Event{ Merge, first, second, dx } ;
		}
		static Event remove( size_t removed, size_t reloc  ) {
			return Event{ Remove, removed, reloc, Vec::Zero() } ;
		}
	};
	struct EventLog {
		void log( const Event& event )	;

		void start() {
			m_log.emplace_back() ;
		}
		void clear() {
			m_log.clear() ;
		}
		const std::vector< std::vector< Event > >& events() const 	{
			return m_log ;
		}

		template < typename Archive >
		void serialize( Archive &ar, unsigned int ) {
			ar & m_log ;
		}

		const std::vector< std::vector< Event > >& log() const {
			return m_log ;
		}

	private:
		std::vector< std::vector< Event > > m_log ;
		Mutex m_log_mutex ;

	};

	static const size_t s_MAX ;

	Particles() ;

	void generate(const ScalarExpr &expr, const unsigned nSamples, const MeshType& mesh ) ;

	size_t count() const { return m_count ; }

	template < typename Archive >
	void serialize( Archive &ar, unsigned int ) ;

	const typename Data< 1 >::Type&    volumes() const { return    m_volumes ; }
	const typename Data< 3 >::Type&    centers() const { return    m_centers ; }
	const typename Data< 3 >::Type& velocities() const { return m_velocities ; }
	const typename Data< 6 >::Type&     frames() const { return     m_frames ; }
	const typename Data< 6 >::Type&     orient() const { return     m_orient ; }


private:

	std::size_t m_count ;

	typename Data< 1 >::Type m_volumes ;

	typename Data< 3 >::Type m_centers ;
	typename Data< 3 >::Type m_velocities ;

	typename Data< 6 >::Type m_frames ;
	typename Data< 6 >::Type m_orient ; // Aniso

	void resize( size_t n ) ;

	friend class DynParticles ;

};

} //d6

#endif

