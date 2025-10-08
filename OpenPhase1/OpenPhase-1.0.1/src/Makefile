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

# Makefile for src/ SUBDIRS.

# ---------------------------------------------------------------------
# Exclude list, i.e. a.cpp b.cpp x.cpp. Will be passed to submakefiles
export EXCLUDE =
# ---------------------------------------------------------------------
include ../Makefile.defs

DEPTH = ..

# Include makefile definitions

# Object SUBDIRS
OBJDIR := ..$(OBJDIR)
# Finds cpp files in current folder and remove path
SRC := $(wildcard *.cpp)
SRC := $(notdir $(SRC))
# Exlude entries from list above
SRC := $(filter-out $(EXCLUDE), $(SRC))

# List of object files to create
OBJS = $(SRC:%.cpp=$(OBJDIR)/%.o)
# List of directories in src/
SUBDIRS = $(sort $(dir $(wildcard */)))
# List of header files in include and its subfolders
HEADERS = $(wildcard ../include/*.h) $(wildcard ../include/**/*.h)
# Dependency files
DEPS = $(patsubst %.o,%.d,$(OBJS))

ALLSRC = $(wildcard *.cpp) $(wildcard */*.cpp) $(wildcard */*/*.cpp)
ALLSRC := $(notdir $(ALLSRC))
# Exlude entries from list above
ALLSRC := $(filter-out $(EXCLUDE), $(ALLSRC))
# List of object files to create
OFILES = $(ALLSRC:%.cpp=$(OBJDIR)/%.o)

# pull in dependency info for *existing* .o files
#-include $(OBJS:.o=.d)

LIB = ../lib/$(LIBNAME)
ifeq ($(findstring OPS, $(SETTINGS)),)
     INCLUDES += -I../include
     LDFLAGS += -L../lib
endif
.PHONY : all clean $(SUBDIRS) subdirs

all: subdirs $(OBJS) $(LIB)

subdirs: $(SUBDIRS)

ifneq ($(findstring shared, $(SETTINGS)),)
$(LIB): $(OBJS) $(SUBDIRS)
	$(CXX) $(SOOPTIONS) -Wl,-soname,$(LIBNAME) -o $(LIB) $(OFILES) $(STDLIBS)
else
$(LIB): $(OBJS)
	ar $(AROPTIONS) $(LIB) $(OBJS)
	ranlib $(LIB)
endif

$(OBJDIR)/%.o: %.cpp Makefile ../Makefile.defs
	$(CXX) -MMD -MP $(CXXFLAGS) $(INCLUDES) $ -c $< -o $@

$(SUBDIRS):
	$(MAKE) -C $@

clean:
	rm -f $(OBJDIR)/*.bak $(OBJDIR)/*.o $(OBJDIR)/*.g* $(OBJDIR)/*.d* $(LIB)

# Include List of dependency files in /obj folder. Will not complain
# if not present at first compilation.
-include $(DEPS)
