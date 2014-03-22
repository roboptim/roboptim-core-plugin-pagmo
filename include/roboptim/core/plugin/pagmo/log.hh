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

#ifndef ROBOPTIM_CORE_PLUGIN_PAGMO_LOG_HH
# define ROBOPTIM_CORE_PLUGIN_PAGMO_LOG_HH

# include <iostream>
# include <sstream>

namespace roboptim {
  namespace pagmo {

    /// \brief Allows to redirect standard output to a log file.
    class LogRedirector
    {
    public:
      LogRedirector ();
      ~LogRedirector ();

      /// \brief Start redirection.
      void start ();

      /// \brief Reset redirection.
      void stop ();

      /// \brief Dump redirected data to file.
      void dump (std::string filename);

    private:
      std::streambuf* initial_cout_;
      std::ostringstream str_cout_;
    };

  } // namespace pagmo
} // namespace roboptim
#endif // ROBOPTIM_CORE_PLUGIN_PAGMO_LOG_HH
