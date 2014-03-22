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

#ifndef ROBOPTIM_CORE_PLUGIN_PAGMO_PROBLEM_WRAPPER_HXX
# define ROBOPTIM_CORE_PLUGIN_PAGMO_PROBLEM_WRAPPER_HXX

# include <boost/shared_ptr.hpp>

# include <roboptim/core/function.hh>
# include <roboptim/core/linear-function.hh>
# include <roboptim/core/problem.hh>
# include <roboptim/core/solver-error.hh>

namespace roboptim
{
  namespace pagmo
  {
    template<typename T>
    ProblemWrapper<T>::ProblemWrapper (const problem_t& pb)
      : base_t (static_cast<int> (pb.function ().inputSize ()),0,1,0,
		// TODO: handle multidimensional constraints
		static_cast<int> (pb.constraints ().size ())),
	pb_ (pb)
    {
      std::vector<double> lb;
      std::vector<double> ub;

      for (typename intervals_t::const_iterator
             iter  = pb.argumentBounds ().begin ();
	   iter != pb.argumentBounds ().end ();
	   ++iter)
        {
	  lb.push_back (iter->first);
	  ub.push_back (iter->second);
        }

      // Set bounds
      set_bounds (lb,ub);
    }

    template <typename T>
    ProblemWrapper<T>::~ProblemWrapper ()
    {
    }

    template <typename T>
    typename ProblemWrapper<T>::base_ptr ProblemWrapper<T>::clone () const
    {
      return base_ptr (new ProblemWrapper<T> (*this));
    }

    template <typename T>
    std::string ProblemWrapper<T>::get_name () const
    {
      return "RobOptim problem";
    }

    template <typename T>
    void ProblemWrapper<T>::objfun_impl
    (fitness_vector_t& f, const decision_vector_t& x) const
    {
      using namespace Eigen;

      // map decision vector (std::vector<double>) to Eigen
      Map<const VectorXd> map_x (x.data (), pb_.function ().inputSize ());

      // objective function
      f[0] = pb_.function () (map_x)[0];
    }

    template <typename T>
    void ProblemWrapper<T>::compute_constraints_impl
    (constraint_vector_t& g, const decision_vector_t& x) const
    {
      using namespace Eigen;

      // map decision vector (std::vector<double>) to Eigen
      Map<const VectorXd> map_x (x.data (), pb_.function ().inputSize ());

      double infinity = function_t::infinity ();

      // PaGMO expects constraints in the form: g <= 0

      typedef typename problem_t::constraints_t::const_iterator cstr_iter_t;

      size_t cstrs_id = 0;
      size_t iter = 0;

      // For each (multidimensional) constraint
      for (cstr_iter_t it = pb_.constraints ().begin ();
	   it != pb_.constraints ().end (); ++it)
        {
	  boost::shared_ptr<Function> cstr;
	  cstr = boost::get<boost::shared_ptr<Function> > (*it);

	  typename function_t::result_t res = (*cstr)(map_x);

	  // For each constraint dimension
	  for (size_t i = 0; i < cstr->outputSize (); ++i)
            {
	      // for each constraint
	      const interval_t&
                interval = pb_.boundsVector ()[cstrs_id][i];

	      double lb = function_t::getLowerBound (interval);
	      double ub = function_t::getUpperBound (interval);
	      bool has_lb = (lb > -infinity);
	      bool has_ub = (ub < infinity);

	      if (has_lb)
                {
		  // if one lower bound, i.e. lb <= g
		  g[iter] = lb - res[iter];
		  iter++;
                }
	      if (has_ub)
                {
		  // if one upper bound, i.e. g <= ub
		  g[iter] = res[iter] - ub;
		  iter++;
                }
            }
	  cstrs_id++;
        }
    }
  } // namespace pagmo
}


#endif // ROBOPTIM_CORE_PLUGIN_PAGMO_PROBLEM_WRAPPER_HXX
