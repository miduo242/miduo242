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

# Makefile for subfolder DamageModels

# Depth relative to Openphase main folder
DEPTH = ../../..

# Include makefile definitions
include $(DEPTH)/Makefile.defs

# Object SUBDIRS
OBJDIR := $(DEPTH)$(OBJDIR)
# Finds cpp files in current folder and remove path
SRC := $(wildcard *.cpp)
SRC := $(notdir $(SRC))
# Exlude entries from EXCLUDE list in src/Makefile
SRC := $(filter-out $(EXCLUDE), $(SRC))

# List of object files to create
OBJS = $(SRC:%.cpp=$(OBJDIR)/%.o)
# List of header files in include and its subfolders
HEADERS = $(wildcard $(DEPTH)/include/*.h) $(wildcard $(DEPTH)/include/**/*.h)
# Dependency files
DEPS = $(patsubst %.o,%.d,$(OBJS))

INCLUDES = -I $(DEPTH)/include/
LIB = $(DEPTH)/lib/$(LIBNAME)

.PHONY : all clean

all: $(LIB)

$(LIB): $(OBJS)
	ar $(AROPTIONS) $(LIB) $(OBJS)
	ranlib $(LIB)

$(OBJDIR)/%.o: %.cpp Makefile
	$(CXX) -MMD -MP $(CXXFLAGS) $(INCLUDES) -c $< -o $@

clean:

# Include List of dependency files in /obj folder. Will not complain
# if not present at first compilation.
-include $(DEPS)
