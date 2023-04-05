ifeq ($(OPENATI_COMP),HITACHI)
  F90   = /opt/ofort90/bin/f90 -Oss -64 -omp -cpp -i,L -fixed -DHF_COMPILER
  F90O3 = /opt/ofort90/bin/f90 -O3 -64 -omp -cpp -i,L -fixed -DHF_COMPILER
else
  F90   = ifort -O3 -m64 -fopenmp -shared-intel -cpp -fixed
  F90O3 = ifort -O0 -m64 -fopenmp -shared-intel -cpp -fixed
endif
CC  = cc

SRCDIR       = .

IFILE = Xabclib_LANCZOS.o  Xabclib_MGramSch.o Xabclib_QSORTD.o LAPACK_dstevd.o BLAS_dstevd.o

all : $(IFILE)

.c.o	:
	$(CC) -c $(INCLUDE) $<


Xabclib_LANCZOS.o : $(SRCDIR)/Xabclib_LANCZOS.f90
	$(F90) -c $(SRCDIR)/Xabclib_LANCZOS.f90
Xabclib_MGramSch.o : $(SRCDIR)/Xabclib_MGramSch.f90
	$(F90) -c $(SRCDIR)/Xabclib_MGramSch.f90
Xabclib_QSORTD.o : $(SRCDIR)/Xabclib_QSORTD.f90
	$(F90) -c $(SRCDIR)/Xabclib_QSORTD.f90

LAPACK_dstevd.o : $(SRCDIR)/LAPACK_dstevd.f90
	$(F90O3) -c $(SRCDIR)/LAPACK_dstevd.f90
BLAS_dstevd.o : $(SRCDIR)/BLAS_dstevd.f90
	$(F90) -c $(SRCDIR)/BLAS_dstevd.f90

.PRECIOUS : $(LOADNAME) ;

clean :
	rm -f $(LOADNAME)
	rm -f *.o

