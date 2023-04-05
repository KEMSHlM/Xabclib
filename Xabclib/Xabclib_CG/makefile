ifeq ($(OPENATI_COMP),HITACHI)
  F90 = /opt/ofort90/bin/f90 -Oss -64 -omp -cpp -i,L -fixed -DHF_COMPILER
else
#  F90 = ifort -O3 -m64 -openmp -shared-intel -cpp -fixed
#  F90 = ifort -O3 -m64 -parallel -shared-intel -cpp -fixed
  F90 = ifort -O3 -m64 -fopenmp -shared-intel -cpp -fixed -fPIC
endif
#F90   = gfortran -O3 -m64 -fopenmp -shared-intel -cpp -ffixed-form -fPIC
CC  = cc

SRCDIR    = .


IFILE = Xabclib_CG.o Xabclib_CG_PCHK.o Xabclib_PCGENS.o Xabclib_PCSLV_SYM.o Xabclib_EVAL_SYM.o Xabclib_IC0D.o Xabclib_DSRMV_ForJacobiIter.o


all : $(IFILE)

.c.o	:
	$(CC) -c $(INCLUDE) $<

Xabclib_CG.o : $(SRCDIR)/Xabclib_CG.f90
	$(F90) -c $(SRCDIR)/Xabclib_CG.f90
Xabclib_CG_PCHK.o : $(SRCDIR)/Xabclib_CG_PCHK.f90
	$(F90) -c $(SRCDIR)/Xabclib_CG_PCHK.f90
Xabclib_PCGENS.o : $(SRCDIR)/Xabclib_PCGENS.f90
	$(F90) -c $(SRCDIR)/Xabclib_PCGENS.f90
Xabclib_PCSLV_SYM.o : $(SRCDIR)/Xabclib_PCSLV_SYM.f90
	$(F90) -c $(SRCDIR)/Xabclib_PCSLV_SYM.f90
Xabclib_EVAL_SYM.o : $(SRCDIR)/Xabclib_EVAL_SYM.f90
	$(F90) -c $(SRCDIR)/Xabclib_EVAL_SYM.f90
Xabclib_IC0D.o : $(SRCDIR)/Xabclib_IC0D.f90
	$(F90) -c $(SRCDIR)/Xabclib_IC0D.f90
Xabclib_DSRMV_ForJacobiIter.o : $(SRCDIR)/Xabclib_DSRMV_ForJacobiIter.f90
	$(F90) -c $(SRCDIR)/Xabclib_DSRMV_ForJacobiIter.f90

clean :
	rm -f *.o
