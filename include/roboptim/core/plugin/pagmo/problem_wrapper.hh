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

#ifndef ROBOPTIM_CORE_PLUGIN_PAGMO_PROBLEM_WRAPPER_HH
# define ROBOPTIM_CORE_PLUGIN_PAGMO_PROBLEM_WRAPPER_HH

# include <map>
# include <set>

# include <boost/mpl/vector.hpp>

# include <roboptim/core/solver.hh>
# include <roboptim/core/solver-state.hh>

# include <pagmo/src/algorithms.h>
# include <pagmo/src/problem/base.h>

namespace roboptim {
  namespace pagmo {
    /// \brief Class that makes a bridge to PaGMO's problem interface.
    template <typename T>
    class ProblemWrapper : public ::pagmo::problem::base
    {
    public:
      typedef ::pagmo::problem::base base_t;
      typedef ::pagmo::problem::base_ptr base_ptr;
      typedef ::pagmo::decision_vector decision_vector_t;
      typedef ::pagmo::fitness_vector fitness_vector_t;
      typedef ::pagmo::constraint_vector constraint_vector_t;

      typedef typename T::callback_t callback_t;
      typedef typename T::problem_t problem_t;
      typedef typename T::result_t result_t;
      typedef typename problem_t::function_t function_t;
      typedef typename problem_t::interval_t interval_t;
      typedef typename problem_t::intervals_t intervals_t;

      ProblemWrapper (const problem_t& pb);
      virtual ~ProblemWrapper ();

      virtual base_ptr clone () const;
      std::string get_name () const;

    protected:
      virtual void objfun_impl
      (fitness_vector_t&, const decision_vector_t&) const;

      virtual void compute_constraints_impl
      (constraint_vector_t&, const decision_vector_t&) const;

    private:
      const problem_t& pb_;

    }; // class ProblemWrapper
  } // namespace pagmo
} // namespace roboptim
# include "roboptim/core/plugin/pagmo/problem_wrapper.hxx"
#endif // ROBOPTIM_CORE_PLUGIN_PAGMO_PROBLEM_WRAPPER_HH
