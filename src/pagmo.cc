// Copyright (c) 2014 CNRS
// Authors: Benjamin Chretien

// This file is part of roboptim-core-plugin-pagmo
// roboptim-core-plugin-pagmo is free software: you can redistribute it
// and/or modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation, either version
// 3 of the License, or (at your option) any later version.

// roboptim-core-plugin-pagmo is distributed in the hope that it will be
// useful, but WITHOUT ANY WARRANTY; without even the implied warranty
// of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Lesser Public License for more details.  You should have
// received a copy of the GNU Lesser General Public License along with
// roboptim-core-plugin-pagmo  If not, see
// <http://www.gnu.org/licenses/>.

#include <cstring>
#include <time.h>

#include <boost/shared_ptr.hpp>

#include <boost/assign/list_of.hpp>
#include <boost/preprocessor/array/elem.hpp>
#include <boost/preprocessor/iteration/local.hpp>
#include <boost/preprocessor/stringize.hpp>

#include <roboptim/core/function.hh>
#include <roboptim/core/linear-function.hh>
#include <roboptim/core/problem.hh>
#include <roboptim/core/solver-error.hh>

#include "roboptim/core/plugin/pagmo/pagmo.hh"
#include "roboptim/core/plugin/pagmo/log.hh"

#include <pagmo/src/pagmo.h>
#include <pagmo/src/archipelago.h>
#include <pagmo/src/topologies.h>

#include <pagmo/src/algorithm/bee_colony.h>
#include <pagmo/src/algorithm/cs.h>
#include <pagmo/src/algorithm/cstrs_co_evolution.h>
#include <pagmo/src/algorithm/cmaes.h>
#include <pagmo/src/algorithm/de.h>
#include <pagmo/src/algorithm/de_1220.h>
#include <pagmo/src/algorithm/ihs.h>
#include <pagmo/src/algorithm/mbh.h>
#include <pagmo/src/algorithm/mde_pbx.h>
#include <pagmo/src/algorithm/ms.h>
#include <pagmo/src/algorithm/nsga2.h>
#include <pagmo/src/algorithm/sa_corana.h>

// TODO: add Ipopt support if Ipopt has been found
//#include <pagmo/src/algorithm/ipopt.h>

namespace roboptim
{
  namespace pagmo
  {
    namespace detail
    {
      template <typename T>
      boost::shared_ptr< ::pagmo::algorithm::base> createInstance ()
      {
	return boost::shared_ptr< ::pagmo::algorithm::base> (new T);
      }

    } // namespace detail

    SolverNlp::SolverNlp (const problem_t& problem) :
      parent_t (problem),
      n_ (problem.function ().inputSize ()),
      m_ (problem.function ().outputSize ()),
      x_ (n_),
      solverState_ (problem),
      wrapper_ (problem),
      algo_map_ ()
    {
      // Initialize x
      x_.setZero ();

      // Load <algo string, algo> map
      algo_map_ = boost::assign::map_list_of
#define N_ALGO 12
#define ALGO_LIST (N_ALGO, (bee_colony,cmaes,cs,cstrs_co_evolution,de,de_1220, \
                            ihs,mbh,mde_pbx,ms,nsga2,sa_corana))
#define GET_ALGO(n) BOOST_PP_ARRAY_ELEM(n,ALGO_LIST)
#define BOOST_PP_LOCAL_MACRO(n)				\
	(std::string (BOOST_PP_STRINGIZE(GET_ALGO(n))), \
	 GET_ALGO(n))
#define BOOST_PP_LOCAL_LIMITS (0,N_ALGO-1)
#include BOOST_PP_LOCAL_ITERATE()
	;
#undef ALGO_LIST
#undef N_ALGO

      // Initialize parameters
      initializeParameters ();

    }

    SolverNlp::~SolverNlp ()
    {
    }

#define DEFINE_PARAMETER(KEY, DESCRIPTION, VALUE)	\
    do {						\
      parameters ()[KEY].description = DESCRIPTION;	\
      parameters ()[KEY].value = VALUE;			\
    } while (0)

    void SolverNlp::initializeParameters ()
    {
      // Clear parameters
      parameters ().clear ();

      // Shared parameters
      DEFINE_PARAMETER ("max-iterations", "number of iterations", 10000);

      // PaGMO-specific parameters
      DEFINE_PARAMETER ("pagmo.candidates", "number of candidates", 30);
      DEFINE_PARAMETER ("pagmo.seed", "random seed", static_cast<int> (time (0)));
      DEFINE_PARAMETER ("pagmo.algorithm", "algorithm", "mde_pbx");
      DEFINE_PARAMETER ("pagmo.local_algorithm", "local algorithm", "de_1220");
      DEFINE_PARAMETER ("pagmo.output_file", "output file", "");
      DEFINE_PARAMETER ("pagmo.penalty_weight",
                        "penalty weight for constrained problems", 1e5);
      DEFINE_PARAMETER ("pagmo.threads", "number of threads", 1);
      // FIXME: use max-iterations instead?
      DEFINE_PARAMETER ("pagmo.generations", "number of generations", 3000);
    }

    void SolverNlp::solve ()
    {
      using namespace Eigen;
      using namespace boost;

      // Load optional starting point
      if (problem ().startingPoint ())
	{
	  x_ = *(problem ().startingPoint ());
	}
      // TODO: find what can be done with given starting point

      if (parameters ().find ("pagmo.algorithm") == parameters ().end ())
	{
          result_ = SolverError ("Undefined PaGMO algorithm.");
          return;
	}

      // Use callback (even if mostly empty)
      if (!callback_.empty ())
	{
          callback_ (problem (), solverState_);
	}

      unsigned int seed = static_cast<unsigned int> (get<int> (parameters ()["pagmo.seed"].value));
      int n_candidates = get<int> (parameters ()["pagmo.candidates"].value);
      int generations = get<int> (parameters ()["pagmo.generations"].value);
      double penalty_weight = get<double> (parameters ()["pagmo.penalty_weight"].value);
      int n_threads = get<int> (parameters ()["pagmo.threads"].value);

      bool is_constrained = (wrapper_.get_c_dimension () > 0);

      std::string algo_str = get<std::string>
	(parameters ()["pagmo.algorithm"].value);


      // 1 - Create population.
      // Here, evolution will take place on the same thread as main.
      typedef ::pagmo::population pop_t;
      typedef shared_ptr<pop_t> pop_ptr;
      pop_ptr pop;

      typedef ::pagmo::archipelago archi_t;
      typedef shared_ptr<archi_t> archi_ptr;
      archi_ptr archi;

      // For NSGA-II at least 5 individuals in the population are needed
      // and the population size must be a multiple of 4.
      if (algo_map_[algo_str] == nsga2)
	n_candidates -= n_candidates%4;

      // If the problem is unconstrained, simply use the given problem
      if (!is_constrained)
	pop = pop_ptr (new pop_t (wrapper_, n_candidates, seed));
      // Else, if the problem is constrained, we rely on death penalty
      // to transform the problem into an unconstrained problem with the
      // weighted sum of the violations added to the objective function.
      else
	{
          std::vector<double> penalties (wrapper_.get_c_dimension (),
                                         penalty_weight);
          pop = pop_ptr (new pop_t
			 (::pagmo::problem::death_penalty
			  (wrapper_,
			   ::pagmo::problem::death_penalty::WEIGHTED,
			   penalties),
			  n_candidates, seed));
	}

      std::string output_file = get<std::string>
	(parameters ()["pagmo.output_file"].value);
      bool use_log = (output_file.compare ("") != 0);
      LogRedirector redirector;

      // 2 - We instantiate the algorithm with a given number of generations.
      // Note: this is where we can change the algorithm that will be used,
      // e.g Differential Evolution, Ipopt etc.

#define SOLVE(ALGO,POP)							\
      if (use_log) {							\
	ALGO.set_screen_output (true);					\
	redirector.start ();						\
      }									\
      if (n_threads > 1) {						\
	/* TODO: handle user-defined topology */			\
	::pagmo::topology::one_way_ring topo;				\
	::pagmo::archipelago::distribution_type				\
	    distrib = ::pagmo::archipelago::point_to_point;		\
	::pagmo::archipelago::migration_direction			\
	    direction = ::pagmo::archipelago::destination;		\
	if (!is_constrained) {						\
	  archi = archi_ptr (new archi_t (ALGO, wrapper_, n_threads,	\
					  n_candidates, topo, distrib, direction)); \
	  archi->evolve ();						\
	} else {							\
          std::vector<double> penalties (wrapper_.get_c_dimension (),	\
                                         penalty_weight);		\
          ::pagmo::problem::death_penalty				\
	      pb = ::pagmo::problem::death_penalty			\
	      (wrapper_,						\
	       ::pagmo::problem::death_penalty::WEIGHTED,		\
	       penalties);						\
	  archi = archi_ptr (new archi_t (ALGO, pb, n_threads,		\
					  n_candidates, topo, distrib, direction)); \
	  archi->evolve ();						\
	}								\
      }									\
      else ALGO.evolve (*POP);						\
      if (use_log) {							\
	redirector.dump (output_file);					\
	redirector.stop ();						\
      }

      // Choose appropriate algorithm
      switch (algo_map_[algo_str])
	{
	case bee_colony:
	  {
	    // Instantiate algorithm
	    ::pagmo::algorithm::bee_colony algo (generations);

	    // Evolve population with the algorithm and solve the problem
	    SOLVE(algo,pop);
	  }
	  break;

	case cmaes:
	  {
	    // Instantiate algorithm
	    ::pagmo::algorithm::cmaes algo (generations);

	    // Evolve population with the algorithm and solve the problem
	    SOLVE(algo,pop);
	  }
	  break;

	case cs:
	  {
	    // Instantiate algorithm
	    ::pagmo::algorithm::cs algo (generations);

	    // Evolve population with the algorithm and solve the problem
	    SOLVE(algo,pop);
	  }
	  break;

	case cstrs_co_evolution:
	  {
	    // Instantiate algorithm
	    ::pagmo::algorithm::de_1220 original (generations);
	    ::pagmo::algorithm::cstrs_co_evolution algo (original);

	    // Evolve population with the algorithm and solve the problem
	    SOLVE(algo,pop);
	  }
	  break;

	case de:
	  {
	    // Instantiate algorithm
	    ::pagmo::algorithm::de algo (generations);

	    // Evolve population with the algorithm and solve the problem
	    SOLVE(algo,pop);
	  }
	  break;

	default:
	case de_1220:
	  {
	    // Instantiate algorithm
	    ::pagmo::algorithm::de_1220 algo (generations);

	    // Evolve population with the algorithm and solve the problem
	    SOLVE(algo,pop);
	  }
	  break;

	case ihs:
	  {
	    // Instantiate algorithm
	    ::pagmo::algorithm::ihs algo (generations);

	    // Evolve population with the algorithm and solve the problem
	    SOLVE(algo,pop);
	  }
	  break;

	case mbh:
	  {
	    // Instantiate algorithm
	    ::pagmo::algorithm::de_1220 local (generations);
	    ::pagmo::algorithm::mbh algo (local);

	    // Evolve population with the algorithm and solve the problem
	    SOLVE(algo,pop);
	  }
	  break;

	case mde_pbx:
	  {
	    // Instantiate algorithm
	    ::pagmo::algorithm::mde_pbx algo (generations);

	    // Evolve population with the algorithm and solve the problem
	    SOLVE(algo,pop);
	  }
	  break;

	case ms:
	  {
	    // Instantiate algorithm
	    ::pagmo::algorithm::de_1220 local (generations);
	    ::pagmo::algorithm::ms algo (local, n_candidates);

	    // Evolve population with the algorithm and solve the problem
	    SOLVE(algo,pop);
	  }
	  break;

	case nsga2:
	  {
	    // Instantiate algorithm
	    ::pagmo::algorithm::nsga2 algo (generations);

	    // Evolve population with the algorithm and solve the problem
	    SOLVE(algo,pop);
	  }
	  break;

	case sa_corana:
	  {
	    // Instantiate algorithm
	    ::pagmo::algorithm::sa_corana algo (generations);

	    // Evolve population with the algorithm and solve the problem
	    SOLVE(algo,pop);
	  }
	  break;

	}

#undef SOLVE

      Result result (n_, m_);

      // multiple threads (archipelago)
      if (n_threads > 1)
	{
	  // TODO: find best comparison strategy for multiobjective functions
	  archi_t::size_type best_i = 0;
          double best_min = archi->get_island (0)->get_population ().champion ().f[0];

          for (archi_t::size_type i = 1; i < archi->get_size (); ++i)
	    {
              if (archi->get_island (i)->get_population ().champion ().f[0] < best_min)
		{
                  best_i = i;
                  best_min = archi->get_island (i)->get_population ().champion ().f[0];
		}
	    }

	  // Note: apparently, using a reference or a map here fails.
	  const ::pagmo::population::champion_type
	    champion = archi->get_island (best_i)->get_population ().champion ();

          result.x.resize (n_);
          result.value.resize (m_);
          for (size_t i = 0; i < static_cast<size_t> (n_); ++i)
	    result.x[static_cast<vector_t::Index> (i)] = champion.x[i];
          for (size_t i = 0; i < static_cast<size_t> (m_); ++i)
	    result.value[static_cast<vector_t::Index> (i)] = champion.f[i];
	}
      else // 1 thread
	{
          Map<const argument_t> map_x (pop->champion ().x.data (), n_);
          Map<const result_t> map_obj (pop->champion ().f.data (), m_);

          result.x = map_x;
          result.value = map_obj;
	}
      function_t::result_t constraints;
      constraints.resize (wrapper_.get_c_dimension ()/2);

      typedef problem_t::constraints_t::const_iterator
	citer_t;
      typedef wrapper_t::linearFunction_t linearFunction_t;
      typedef wrapper_t::nonlinearFunction_t nonlinearFunction_t;

      function_t::size_type idx = 0;
      for (citer_t it = problem ().constraints ().begin ();
           it != problem ().constraints ().end (); ++it)
	{
          shared_ptr<DifferentiableFunction> g;
          if (it->which () == linearFunctionId)
	    g = get<shared_ptr<linearFunction_t> > (*it);
          else
	    g = get<shared_ptr<nonlinearFunction_t> > (*it);

          constraints.segment (idx, g->outputSize ()) = (*g) (result.x);
          idx += g->outputSize ();
	}

      // Fill result structure
      result.constraints = constraints;
      result_ = result;

      // Use callback for last iteration
      if (!callback_.empty ())
	{
          solverState_.x () = result.x;
          solverState_.cost () = result.value[0];
          callback_ (problem (), solverState_);
	}

    }
  } // namespace pagmo
} // end of namespace roboptim

extern "C"
{
  using namespace roboptim::pagmo;
  typedef SolverNlp::parent_t solver_t;

  ROBOPTIM_DLLEXPORT unsigned getSizeOfProblem ();
  ROBOPTIM_DLLEXPORT const char* getTypeIdOfConstraintsList ();
  ROBOPTIM_DLLEXPORT solver_t* create (const SolverNlp::problem_t& pb);
  ROBOPTIM_DLLEXPORT void destroy (solver_t* p);

  ROBOPTIM_DLLEXPORT unsigned getSizeOfProblem ()
  {
    return sizeof (solver_t::problem_t);
  }

  ROBOPTIM_DLLEXPORT const char* getTypeIdOfConstraintsList ()
  {
    return typeid (solver_t::problem_t::constraintsList_t).name ();
  }

  ROBOPTIM_DLLEXPORT solver_t* create (const SolverNlp::problem_t& pb)
  {
    return new SolverNlp (pb);
  }

  ROBOPTIM_DLLEXPORT void destroy (solver_t* p)
  {
    delete p;
  }
}
