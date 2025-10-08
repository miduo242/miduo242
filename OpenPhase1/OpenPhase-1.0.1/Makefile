#   This file is part of the OpenPhase (R) software library.
#
#   Copyright (c) 2009-2022 Ruhr-Universitaet Bochum,
#                 Universitaetsstrasse 150, D-44801 Bochum, Germany
#             AND 2018-2022 OpenPhase Solutions GmbH,
#                 Universitaetsstrasse 136, D-44799 Bochum, Germany.
#
#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#   
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program.  If not, see <http://www.gnu.org/licenses/>.

# Root make file for OpenPhase project.

.PHONY : all clean cleanall clean_doc examples benchmarks docs
SETTINGS = default
include Makefile.defs

all: Makefile
	$(call print_header,Release)
	@echo "Compile main library."
	$(MAKE) --no-print-directory -C src
	$(call print_ending,$(COMPILE_CURRENT_TIME),$(COMPILE_START_TIME))

examples:
	make -C examples

exercises:
	make -C exercises

benchmarks:
	make -C benchmarks

docs:
	doxygen

clean:
	rm -rf hdf5 HighFive tinyxml2
	make --no-print-directory -C src clean
	make --no-print-directory -C benchmarks clean
	make --no-print-directory -C examples clean

cleanall:
	rm -rf hdf5 HighFive tinyxml2
	make --no-print-directory -C src clean
	make --no-print-directory -C benchmarks cleanall
	make --no-print-directory -C examples cleanall

clean_doc:
	-rm -r ./documentation/doxygen
