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

#include "roboptim/core/plugin/pagmo/log.hh"

#include <iostream>
#include <fstream>

namespace roboptim
{
  namespace pagmo
  {
    LogRedirector::LogRedirector ()
      : initial_cout_ (0x0)
    {
    }

    LogRedirector::~LogRedirector ()
    {
      stop ();
    }

    void LogRedirector::start ()
    {
      // Get initial standard output
      initial_cout_ = std::cout.rdbuf ();

      // Redirect stringstream instead
      std::cout.rdbuf (this->str_cout_.rdbuf ());
    }

    void LogRedirector::stop ()
    {
      if (initial_cout_ != 0x0)
	{
	  std::cout.rdbuf (initial_cout_);
	  initial_cout_ = 0x0;
	}
    }

    void LogRedirector::dump (std::string filename)
    {
      std::ofstream log_file;
      log_file.open (filename.c_str (), std::ofstream::out);
      log_file << str_cout_.str ();
      log_file.close ();
    }
  } // namespace pagmo
} // end of namespace roboptim
