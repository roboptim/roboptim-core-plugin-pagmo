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
#include <map>
#include <limits> // epsilon

#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>

#include <roboptim/core/function.hh>
#include <roboptim/core/linear-function.hh>
#include <roboptim/core/problem.hh>
#include <roboptim/core/solver-error.hh>

#include "roboptim/core/plugin/pagmo/pagmo.hh"

namespace roboptim
{
  namespace pagmo
  {
    SolverNlp::SolverNlp (const problem_t& problem) :
      parent_t (problem),
      n_ (problem.function ().inputSize ()),
      m_ (problem.function ().outputSize ()),
      x_ (n_),
      solverState_ (problem),
      algo_map_ (),
      global_algos_ ()
    {
      // Initialize x
      x_.setZero ();
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
    }

    void SolverNlp::solve () throw ()
    {
      using namespace Eigen;
      typedef VectorXd::Index index_t;

      // Load optional starting point
      if (problem ().startingPoint ())
	{
	  x_ = *(problem ().startingPoint ());
	}

      result_ = SolverError ("Error");
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
