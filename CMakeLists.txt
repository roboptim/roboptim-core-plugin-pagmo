# Copyright (c) 2014 CNRS
# Authors: Benjamin Chretien
#
#
# This file is part of roboptim-core-plugin-pagmo
# roboptim-core-plugin-pagmo is free software: you can redistribute it
# and/or modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation, either version
# 3 of the License, or (at your option) any later version.
#
# roboptim-core-plugin-pagmo is distributed in the hope that it will be
# useful, but WITHOUT ANY WARRANTY; without even the implied warranty
# of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Lesser Public License for more details.  You should have
# received a copy of the GNU Lesser General Public License along with
# roboptim-core-plugin-pagmo  If not, see
# <http://www.gnu.org/licenses/>.

CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

#SET(CXX_DISABLE_WERROR TRUE)
INCLUDE(cmake/base.cmake)
INCLUDE(cmake/boost.cmake)
INCLUDE(cmake/eigen.cmake)
INCLUDE(cmake/cpack.cmake)

SET(PROJECT_NAME roboptim-core-plugin-pagmo)
SET(PROJECT_DESCRIPTION "Nonlinear optimization using PaGMO.")
SET(PROJECT_URL "http://github.com/roboptim/roboptim-core-plugin-pagmo")

SET(PKG_CONFIG_ADDITIONAL_VARIABLES plugindir ${PKG_CONFIG_ADDITIONAL_VARIABLES})
SETUP_PROJECT()

# Declare headers
SET(${PROJECT_NAME}_HEADERS
  include/roboptim/core/plugin/pagmo/pagmo.hh
  include/roboptim/core/plugin/pagmo/log.hh
  include/roboptim/core/plugin/pagmo/problem_wrapper.hh
  include/roboptim/core/plugin/pagmo/problem_wrapper.hxx
  )

# Search for dependencies.
SET(BOOST_COMPONENTS date_time filesystem system unit_test_framework)
SEARCH_FOR_BOOST()

ADD_REQUIRED_DEPENDENCY("eigen3 >= 3.2.0")
ADD_REQUIRED_DEPENDENCY("roboptim-core >= 3.0")
ADD_REQUIRED_DEPENDENCY("pagmo")

MESSAGE(STATUS "ROBOPTIM_CORE_PLUGINDIR: " ${ROBOPTIM_CORE_PLUGINDIR})

OPTION(DISABLE_TESTS "Disable test programs" OFF)

ADD_SUBDIRECTORY(src)

IF(NOT DISABLE_TESTS)
  ADD_SUBDIRECTORY(tests)
ELSE()
  MESSAGE(WARNING
    "Tests should only be disabled for speficic cases. Do it at your own risk.")
ENDIF()

HEADER_INSTALL("${HEADERS}")

SETUP_PROJECT_FINALIZE()
SETUP_PROJECT_CPACK()
