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

#include <pagmo/src/algorithm/bee_colony.h>
#include <pagmo/src/algorithm/cs.h>
#include <pagmo/src/algorithm/cstrs_co_evolution.h>
#include <pagmo/src/algorithm/cmaes.h>
#include <pagmo/src/algorithm/de.h>
#include <pagmo/src/algorithm/de_1220.h>
#include <pagmo/src/algorithm/ihs.h>
#include <pagmo/src/algorithm/mbh.h>

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
#define N_ALGO 8
#define ALGO_LIST (N_ALGO, (bee_colony,cmaes,cs,cstrs_co_evolution,de,de_1220, \
                            ihs,mbh))
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

    SolverNlp::~SolverNlp () throw ()
    {
    }

#define DEFINE_PARAMETER(KEY, DESCRIPTION, VALUE)	\
    do {						\
      parameters ()[KEY].description = DESCRIPTION;	\
      parameters ()[KEY].value = VALUE;			\
    } while (0)

    void SolverNlp::initializeParameters () throw ()
    {
      // Clear parameters
      parameters ().clear ();

      // Shared parameters
      DEFINE_PARAMETER ("max-iterations", "number of iterations", 10000);

      // PaGMO-specific parameters
      DEFINE_PARAMETER ("pagmo.candidates", "number of candidates", 20);
      DEFINE_PARAMETER ("pagmo.seed", "random seed", 123);
      DEFINE_PARAMETER ("pagmo.algorithm", "algorithm", "de_1220");
      DEFINE_PARAMETER ("pagmo.algorithm", "algorithm", "mbh");
      DEFINE_PARAMETER ("pagmo.output_file", "output file", "");
      // FIXME: use max-iterations instead?
      DEFINE_PARAMETER ("pagmo.generations", "number of generations", 3000);
    }

    void SolverNlp::solve () throw ()
    {
      using namespace Eigen;

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

      int seed = boost::get<int> (parameters ()["pagmo.seed"].value);
      int n_c = boost::get<int> (parameters ()["pagmo.candidates"].value);
      int generations = boost::get<int> (parameters ()["pagmo.generations"].value);

      // 1 - Create population.
      // Here, evolution will take place on the same thread as main.
      ::pagmo::population pop (wrapper_, n_c, seed);

      // TODO: redirect verbose to log file
      std::string output_file = boost::get<std::string>
          (parameters ()["pagmo.output_file"].value);
      bool use_log = (output_file.compare ("") != 0);
      LogRedirector redirector;

      // 2 - We instantiate the algorithm with a given number of generations.
      // Note: this is where we can change the algorithm that will be used,
      // e.g Differential Evolution, Ipopt etc.
      std::string algo_str = boost::get<std::string>
	(parameters ()["pagmo.algorithm"].value);

#define SOLVE(ALGO,POP) \
      if (use_log) { \
          ALGO.set_screen_output (true); \
          redirector.start (); \
      } \
      ALGO.evolve (POP); \
      if (use_log) redirector.dump (output_file);

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
	}

    #undef SOLVE

      // Note: for parallel processing, evolutions can take place in parallel
      // on multiple separate islands containing, each several candidate
      // solutions to the problem.
      //::pagmo::archipelago archi (algo, wrapper_, 8, 20);
      //archi.evolve ();

      Map<const argument_t> map_x (pop.champion ().x.data (), n_);
      Map<const result_t> map_obj (pop.champion ().f.data (), m_);
      Map<const vector_t> map_cstr (pop.champion ().c.data (),
                                    pop.champion ().c.size ());

      // Fill result structure
      Result result (n_, 1);
      result.x = map_x;
      result.value = map_obj;
      result.constraints = map_cstr;
      result_ = result;
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
