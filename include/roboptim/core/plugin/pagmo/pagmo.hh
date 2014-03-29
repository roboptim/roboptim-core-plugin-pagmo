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

#ifndef ROBOPTIM_CORE_PLUGIN_PAGMO_PAGMO_HH
# define ROBOPTIM_CORE_PLUGIN_PAGMO_PAGMO_HH

# include <map>

# include <roboptim/core/solver.hh>
# include <roboptim/core/solver-state.hh>

# include "roboptim/core/plugin/pagmo/problem_wrapper.hh"

namespace roboptim {
  namespace pagmo {

    // Available PaGMO algorithms
    enum Algorithm
      {
	bee_colony,         // Artificial Bee Colony
	cmaes,              // Covariance Matrix Adaptation Evolutionary Strategy
	cs,                 // Compass Search Solver
	cstrs_co_evolution, // Co-Evolution constraints handling
	de,                 // Differential Evolution
	de_1220,            // Differential Evolution - 1220 variant
	ihs,                // Improved Harmony Search
	ipopt,              // Ipopt
	mbh,                // Monotonic Basin Hopping
	mde_pbx,            // MDE_pBX - Differential Evolution variant
	ms,                 // Multi-start
	nsga2,              // Nondominated Sorting genetic algorithm II
	sa_corana           // Simulated Annealing with adaptive neighbourhood
      };

    /// \brief Solver interfacing with the PaGMO library.
    class SolverNlp :
      public Solver<DifferentiableFunction,
                    boost::mpl::vector<LinearFunction,
                                       DifferentiableFunction> >
    {
    public:
      /// \brief Parent type
      typedef Solver<DifferentiableFunction,
                     boost::mpl::vector<LinearFunction,
                                        DifferentiableFunction> > parent_t;
      /// \brief Cost function type
      typedef problem_t::function_t function_t;
      /// \brief Argument type
      typedef function_t::argument_t argument_t;
      /// \brief type of result
      typedef function_t::result_t result_t;
      /// \brief type of gradient
      typedef DifferentiableFunction::gradient_t gradient_t;
      /// \brief Size type
      typedef Function::size_type size_type;
      /// \brief Constraints type
      typedef problem_t::constraints_t constraints_t;
      /// \brief Constraint type
      typedef problem_t::constraint_t constraint_t;
      /// \brief Intervals type
      typedef problem_t::intervals_t intervals_t;
      /// \brief Interval type
      typedef problem_t::interval_t interval_t;

      /// \brief Solver state
      typedef SolverState<parent_t::problem_t> solverState_t;

      /// \brief RobOptim callback
      typedef parent_t::callback_t callback_t;

      /// \brief PaGMO problem wrapper
      typedef ProblemWrapper<parent_t> wrapper_t;

      /// \brief Constructor by problem
      explicit SolverNlp (const problem_t& problem);
      virtual ~SolverNlp () throw ();
      /// \brief Solve the optimization problem
      virtual void solve () throw ();

      /// \brief Return the number of variables.
      size_type n () const
      {
	return n_;
      }

      /// \brief Return the number of functions.
      size_type m () const
      {
	return m_;
      }

      /// \brief Get the optimization parameters.
      Function::argument_t& parameter ()
      {
	return x_;
      }

      /// \brief Get the optimization parameters.
      const Function::argument_t& parameter () const
      {
	return x_;
      }

      /// \brief Set the callback called at each iteration.
      virtual void
      setIterationCallback (callback_t callback) throw (std::runtime_error)
      {
        callback_ = callback;
      }

      /// \brief Get the callback called at each iteration.
      const callback_t& callback () const throw ()
      {
        return callback_;
      }

    private:
      void initializeParameters () throw ();

    public:
      static const int linearFunctionId = 0;
      static const int nonlinearFunctionId = 1;

    private:
      /// \brief Number of variables
      size_type n_;
      /// \brief Dimension of the cost function
      size_type m_;

      /// \brief Parameter of the function
      Function::argument_t x_;

      /// \brief State of the solver at each iteration
      solverState_t solverState_;

      /// \brief Intermediate callback (called at each end of iteration).
      callback_t callback_;

      /// \brief PaGMO problem wrapper.
      wrapper_t wrapper_;

      /// \brief Map string to PaGMO algorithm
      std::map<std::string, Algorithm> algo_map_;
    }; // class SolverNlp
  } // namespace pagmo
} // namespace roboptim
#endif // ROBOPTIM_CORE_PLUGIN_PAGMO_PAGMO_HH
