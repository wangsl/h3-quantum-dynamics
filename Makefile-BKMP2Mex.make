
# $Id$

.DEFAULT: .f .for .c .f90
.SUFFIXES: .f .for .c .f90

O = .

F77 = ifort
F90 = ifort
CC = icc
CXX = icpc

CFLAGS =  

CXXFLAGS = -I$(MATLAB_ROOT)/extern/include

FFLAGS = 

Link = $(CXX) -shared-intel -shared

LIBS = -lifcoremt

EXE = BKMP2Mex.mexa64

OBJS = $(O)/bkmp2.o  $(O)/BKMP2Mex.o

$(EXE) : $(OBJS)
	$(Link) $(CXXFLAGS) -o $(EXE) $(OBJS) $(LIBS)

$(O)/%.o: %.c
	cd $(O) ; $(CC)  $(cFLAGS) -c $<
$(O)/%.o: %.cc
	cd $(O) ; $(CXX) $(CXXFLAGS) -c $<
$(O)/%.o: %.cpp
	cd $(O) ; $(CXX) $(CXXFLAGS) -c $<
$(O)/%.o: %.C
	cd $(O) ; $(CXX) $(CXXFLAGS) -c $<
$(O)/%.o: %.f
	cd $(O) ; $(F77) $(FFLAGS) -c $<
$(O)/%.o: %.for
	cd $(O) ; $(F77) $(FFLAGS) -c $<
$(O)/%.o: %.f90
	cd $(O) ; $(F90) $(FFLAGS) -c $<

clean:
	rm -f *.o *.dat *~ *.mod *.ti *.ii $(EXE) depend

depend :
	$(CXX) $(CXXFLAGS) -MM *.[cC] | perl dep.pl | tee $@

include depend

