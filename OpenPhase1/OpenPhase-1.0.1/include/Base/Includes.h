/*
 *   This file is part of the OpenPhase (R) software library.
 *
 *   Copyright (c) 2009-2022 Ruhr-Universitaet Bochum,
 *                 Universitaetsstrasse 150, D-44801 Bochum, Germany
 *             AND 2018-2022 OpenPhase Solutions GmbH,
 *                 Universitaetsstrasse 136, D-44799 Bochum, Germany.
 *
 *   This program is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 *   File created :   2011
 *   Main contributors :   Oleg Shchyglo; Efim Borukhovich; Dmitry Medvedev
 *
 */

// Standard C++ headers

#ifndef INCLUDES_H
#define INCLUDES_H

#include <vector>
#include <array>
#include <algorithm>
#include <typeinfo>
#include <type_traits>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <ctime>
#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <sys/stat.h>
#include <sys/types.h>
#include <cmath>
#include <complex>
#include <map>
#include <cfloat>
#include <limits>
#include <errno.h>
#include <signal.h>
#include <stdexcept>
#include <random>
#include <ciso646>
#include <functional>
#include <exception>
#include <algorithm>
#include <cassert>
#include <any>

#ifndef SERIAL
#include <omp.h>
#endif

#ifdef MPI_PARALLEL
#include "mpi.h"
#endif

#include "Base/Definitions.h"
#include "Base/Macros.h"
#include "Base/Storages.h"
#include "Base/OPObject.h"
#include "GradientStencil.h"
#include "LaplacianStencil.h"
#include "PhysicalConstants.h"

#endif

