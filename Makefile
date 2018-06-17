MKLROOT=/progs/intel/composerxe/mkl/
MPIROOT=/progs/intel/impi/4.0.3.008/


PROGRAM=APSOgenEAM.x
OBJS = APSOgenEAM.o

VARIABLE_FLAGS=
debug := VARIABLE_FLAGS=-g
d64 := VARIABLE_FLAGS=-double-size 64 -real-size 64
d128 := VARIABLE_FLAGS=-double-size 128 -real-size 128

FFLAGS = $(FLAGS) $(VARIABLE_FLAGS)

all debug d64 d128:	$(PROGRAM)
#ifort -c $< -o $@ -fp-model precise -real-size 128 -double-size 128 -i4 -shared-intel $(FFLAGS) -prec-div -g

%.o: %.f
	${MPIROOT}/bin64/mpiifort -c $< -o $@ -i8 -I$(MKLROOT)/include/intel64/ilp64  -I$(MKLROOT)/include -fp-model precise -real-size 64 -double-size 64 $(FFLAGS) -prec-div

$(PROGRAM) :	$(OBJS)
	${MPIROOT}/bin64/mpiifort -o $(PROGRAM) $(FFLAGS) $(OBJS)  $(MKLROOT)/lib/intel64/libmkl_lapack95_lp64.a  -Wl,--start-group  $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a $(MKLROOT)/lib/intel64/libmkl_sequential.a $(MKLROOT)/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm

clean: 
	-rm $(OBJS) APSOgenEAM.x
