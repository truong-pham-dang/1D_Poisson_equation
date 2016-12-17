CF         = gfortran
FFLAGS     = -O3 -g
LD         = gfortran
LDFLAGS    = 
PREPROC    = 

OBJS =  1D_Poisson.o \
        matinv.o  \

.SUFFIXES: .o .f90 .f
.f90.o:
	$(LD) -c $(FFLAGS) $<
.f.o:
	$(LD) -c $(FFLAGS) $<

code :$(OBJS) 
	$(LD) $(LDFLAGS) -o $@ $(OBJS) $(MPILIBDIR) $(MPILIBS)

clean :
	rm -f code *.o core *.mod

