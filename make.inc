PLAT = _OPENMP

BLASDEF   =
BLASLIB   = -mkl 
MATHLIB   = -lm
MPLIB     =

#  The archiver and the flag(s) to use when building archive (library)
#  If your system has no ranlib, set RANLIB = echo.
#
ARCH         = ar
ARCHFLAGS    = cr
RANLIB       = ranlib
 
# Definitions used by CPP:
PREDEFS         = -D__OPENMP
#

CC           = icc -qopenmp
CFLAGS       = -O3 $(PREDEFS) -DPRNTlevel=0 -DDEBUGlevel=0 -qopt-report-phase=vec,loop -qopt-report-file=stderr -lpthread
#
# There are a few files that should NOT be optimized, which appear in the
# end of SRC/Makefile, taking the flag $(NOOPTS)
# NOOPTS should be defined to explicitly turn off any optimization.
NOOPTS	 = -O0

############################################################################
FORTRAN      = ifort -qopenmp
#FFLAGS       = -O3 -align array64byte -qopt-report-phase=offload,vec,loop -qopt-report-file=stderr -lpthread -mkl -xAVX -cpp
FFLAGS       = -O3 -align array64byte -lpthread -mkl -xAVX -cpp
LOADER       = icc -qopenmp
LOADOPTS     =

#
#  C preprocessor defs for compilation for the Fortran interface
#  (-DNoChange, -DAdd_, -DAdd__, or -DUpCase)
#
CDEFS        = -DAdd_
