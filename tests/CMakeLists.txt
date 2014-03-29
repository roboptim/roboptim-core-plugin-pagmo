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


# Prepare test suite
SET(PLUGIN_PATH "${CMAKE_BINARY_DIR}/src")
INCLUDE(shared-tests/tests.cmake)

# Schittkowski problems that are expected to fail
SET(SCHITTKOWSKI_PROBLEMS
      1  2  3  4  5  6  7  8  9 # Missing problems in PDF:
  10 11 12 13 14 15 16 17 18 19 # ------------------------
  20 21 22 23 24 25 26 27 28 29 #
  30 31 32 33 34 35 36 37 38 39 #
  40 41    43 44 45    47 48 49 # 46
  50 51 52 53    55 56          # 58
  60 61    63 64 65 66          #
     71 72          76    78 79 #
  80 81                         #

  # Particular variants of some problems.
  71b)

SET(SCHITTKOWSKI_PROBLEMS_FAIL
    6 7 10 13 15 24 26 27 31 39 40 47 56 60 61 63 71 78 79 80 81 71b)

SET(SOLVER_NAME "pagmo")
SET(FUNCTION_TYPE ::roboptim::EigenMatrixDense)
SET(PROGRAM_SUFFIX "")
SET(COST_FUNCTION_TYPE ::roboptim::GenericDifferentiableFunction)
SET(CONSTRAINT_TYPE_1 ::roboptim::GenericLinearFunction)
SET(CONSTRAINT_TYPE_2 ::roboptim::GenericDifferentiableFunction)
BUILD_SCHITTKOWSKI_PROBLEMS()

# RobOptim test problems
SET(SOLVER_NAME "pagmo")
SET(FUNCTION_TYPE ::roboptim::EigenMatrixDense)
SET(PROGRAM_SUFFIX "")
SET(COST_FUNCTION_TYPE ::roboptim::GenericDifferentiableFunction)
BUILD_ROBOPTIM_PROBLEMS()
