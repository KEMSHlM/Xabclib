ifeq ($(OPENATI_COMP),HITACHI)
  F90   = /opt/ofort90/bin/f90 -Oss -64 -omp -cpp -i,L -fixed -DHF_COMPILER
  F90O3 = /opt/ofort90/bin/f90 -O3 -64 -omp -cpp -i,L -fixed -DHF_COMPILER
else
  F90   = ifort -O3 -m64 -fopenmp -shared-intel -cpp -fixed
  F90O3 = ifort -O0 -m64 -fopenmp -shared-intel -cpp -fixed
endif
#F90O3 = $(F90)
CC  = cc

SRCDIR       = .

IFILE = Xabclib_ARNOLDI.o LAPACK_dgeev.o BLAS_dgeev.o

all : $(IFILE)

.c.o	:
	$(CC) -c $(INCLUDE) $<


Xabclib_ARNOLDI.o : $(SRCDIR)/Xabclib_ARNOLDI.f90
	$(F90) -c $(SRCDIR)/Xabclib_ARNOLDI.f90

LAPACK_dgeev.o : $(SRCDIR)/LAPACK_dgeev.f90
	$(F90O3) -c $(SRCDIR)/LAPACK_dgeev.f90
BLAS_dgeev.o : $(SRCDIR)/BLAS_dgeev.f90
	$(F90) -c $(SRCDIR)/BLAS_dgeev.f90

.PRECIOUS : $(LOADNAME) ;

clean :
	rm -f $(LOADNAME)
	rm -f *.o

