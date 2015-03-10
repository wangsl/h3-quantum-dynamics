
# $Id$

.DEFAULT: .F .F90 .c .C
.SUFFIXES: .F .F90 .c .C

O = .

F77 = ifort
F90 = ifort
CC = icc
CXX = icpc

CFLAGS =  

CXXFLAGS = -std=c++0x -I$(MATLAB_ROOT)/extern/include

FFLAGS = 

Link = $(CXX)

LIBS = -lfftw3_omp -lifcoremt

MEXA64Files = \
	TimeEvolutionMex.mexa64 \
	BKMP2Mex.mexa64 \
	GaussLegendreMex.mexa64

OBJS = $(O)/matutils.o $(O)/indent.o $(O)/die.o $(O)/out.o \
	$(O)/MatlabStructures.o $(O)/MatlabStructuresio.o \
	$(O)/$(O)/bkmp2.o \
	$(O)/GaussLegendre.o $(O)/sortcpp.o \
	$(O)/psitest.o \
	$(O)/fftwinterface.o $(O)/timeEvol.o  $(O)/LegTransform.o \
	$(O)/zeros.o $(O)/FortAux.o

QMLibs = libqmdyn.a

.DEFAULT_GOAL := TimeEvolutionMex.mexa64

all: $(MEXA64Files)

#$(EXE) : $(OBJS)
#	$(Link) $(CXXFLAGS) -o $(EXE) $(OBJS) $(LIBS)

$(O)/%.o: %.c
	cd $(O) ; $(CC)  $(cFLAGS) -c $<
$(O)/%.o: %.C
	cd $(O) ; $(CXX) $(CXXFLAGS) -c $<
$(O)/%.o: %.F
	cd $(O) ; $(F77) $(FFLAGS) -c $<
$(O)/%.o: %.F90
	cd $(O) ; $(F90) $(FFLAGS) -c $<
%io.C: %.h
	perl io.pl $<

$(QMLibs): $(OBJS)
	ar -crusv $(QMLibs) $(OBJS)

%.mexa64: %.o $(QMLibs)
	$(Link) -shared $(CXXFLAGS) -o $@ $^ $(LIBS)

clean:
	rm -f *.o *~ *.mod $(EXE) depend $(MEXA64Files) $(QMLibs)

depend :
	$(CXX) $(CXXFLAGS) -MM *.[cC] | perl dep.pl | tee $@

include depend
