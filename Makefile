# Copyright (c) 2017, 
# Eyal Shalev (eyal@gsi.gov.il)
# Vladimir Lyakhovsky
# Harel Levin (harellevin@gmail.com)
# Gal Oren (galoren.com@gmail.com)
# All rights reserved to:
# Geological Survey of Israel (GSI) &
# Nuclear Research Center - Negev (NRCN).
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#    * Redistributions of source code must retain the above copyright
#      notice, this list of conditions and the following disclaimer.
#    * Redistributions in binary form must reproduce the above copyright
#      notice, this list of conditions and the following disclaimer in the
#      documentation and/or other materials provided with the distribution.
#    * Neither the name of Harel Levin or Gal Oren, nor the
#      names of its contributors may be used to endorse or promote products
#      derived from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL Harel Levin & Gal Oren BE LIABLE FOR ANY
# DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

OS := $(shell uname)

include make.inc.$(OS)
 
.SUFFIXES: 
.SUFFIXES: .o .c .cpp .f90 .f .mod
HEADER  = -I$(MPI_DIR)/include -I$(TRILINOS_INC) 
LIBS	= $(MATHLIB) $(BLASLIB) $(TRILINOSLIBS) $(MPLIB)

OBJ = solver.o modules.o main_3d.o input.o functions.o  node_mass.o\
efdlm.o elastic.o move_grid.o output.o plastic.o damage.o drop.o\
pformpar.o tetra4.o \
diffusion.o bc_stress.o derivation.o

all: poro.exe dat2tec.exe

dat2tec.exe: 
	$(FORTRAN) dat2tec.f90 $(FFLAG) $(LOADOPTS) -o $@
 
poro.exe: $(OBJ)
	$(LOADER) $(FFLAG) $(LOADOPTS) $(OBJ) $(LIBS) $(CUDALIBS) -o $@
	
	
.c.o:
	$(CC) $(CFLAGS) $(CDEFS) $(HEADER) -c $< $(VERBOSE)

.cpp.o:
	$(CXX) $(CXXFLAGS) $(CDEFS) $(HEADER) -c $< $(VERBOSE)

solver.cpp: solver_$(PETRA).cpp
	touch solver.cpp; \
	rm solver.cpp; \
	ln -s solver_$(PETRA).cpp $@
.f90.o:
	$(FORTRAN) $(FFLAGS) -c $< $(VERBOSE)
	
clean:	
	rm -f *.o *.exe solver.cpp

doc:
	doxygen

docsave:
	git commit -a
	git push
