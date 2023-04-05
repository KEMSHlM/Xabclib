ifeq ($(OPENATI_COMP),HITACHI)
  F90 = /opt/ofort90/bin/f90 -Oss -64 -omp -cpp -i,L -fixed -DHF_COMPILER
else
  F90 = ifort -O3 -m64 -fopenmp -shared-intel -cpp -fixed
#  F90 = ifort -O3 -m64 -parallel -mcmodel=medium -cpp -fixed
endif
CC  = cc

SRCDIR   = .

IFILE = OpenATI_EIGENSOLVE.o OpenATI_LINEARSOLVE.o OpenATI_UTIL.o OpenATI_DAFRT.o OpenATI_DSRMV.o OpenATI_DURMV.o OpenATI_DSRMV_Setup.o OpenATI_DURMV_Setup.o OpenATI_DAFGS.o OpenATI_NZDIV.o OpenATI_init.o OpenATI_GET_VER.o OpenATI_DAFMC_CCS2CRS.o OpenATI_DAFSTG.o OpenATI_QSORTI.o

all : $(IFILE)

.c.o	:
	$(CC) -c $(INCLUDE) $<

OpenATI_EIGENSOLVE.o : $(SRCDIR)/OpenATI_EIGENSOLVE.f90
	$(F90) -c $(SRCDIR)/OpenATI_EIGENSOLVE.f90
OpenATI_LINEARSOLVE.o : $(SRCDIR)/OpenATI_LINEARSOLVE.f90
	$(F90) -c $(SRCDIR)/OpenATI_LINEARSOLVE.f90
OpenATI_UTIL.o : $(SRCDIR)/OpenATI_UTIL.f90
	$(F90) -c $(SRCDIR)/OpenATI_UTIL.f90

OpenATI_DSRMV.o : $(SRCDIR)/OpenATI_DSRMV.f90
	$(F90) -c $(SRCDIR)/OpenATI_DSRMV.f90
OpenATI_DURMV.o : $(SRCDIR)/OpenATI_DURMV.f90
	$(F90) -c $(SRCDIR)/OpenATI_DURMV.f90
OpenATI_DAFRT.o : $(SRCDIR)/OpenATI_DAFRT.f90
	$(F90) -c $(SRCDIR)/OpenATI_DAFRT.f90

OpenATI_DAFGS.o : $(SRCDIR)/OpenATI_DAFGS.f90
	$(F90) -c $(SRCDIR)/OpenATI_DAFGS.f90

OpenATI_DSRMV_Setup.o : $(SRCDIR)/OpenATI_DSRMV_Setup.f90
	$(F90) -c $(SRCDIR)/OpenATI_DSRMV_Setup.f90
OpenATI_DURMV_Setup.o : $(SRCDIR)/OpenATI_DURMV_Setup.f90
	$(F90) -c $(SRCDIR)/OpenATI_DURMV_Setup.f90

OpenATI_NZDIV.o : $(SRCDIR)/OpenATI_NZDIV.f90
	$(F90) -c $(SRCDIR)/OpenATI_NZDIV.f90

OpenATI_init.o : $(SRCDIR)/OpenATI_init.f90
	$(F90) -c $(SRCDIR)/OpenATI_init.f90
OpenATI_GET_VER.o : $(SRCDIR)/OpenATI_GET_VER.f90
	$(F90) -c $(SRCDIR)/OpenATI_GET_VER.f90
OpenATI_DAFMC_CCS2CRS.o : $(SRCDIR)/OpenATI_DAFMC_CCS2CRS.f90
	$(F90) -c $(SRCDIR)/OpenATI_DAFMC_CCS2CRS.f90

OpenATI_DAFSTG.o : $(SRCDIR)/OpenATI_DAFSTG.f90
	$(F90) -c $(SRCDIR)/OpenATI_DAFSTG.f90

OpenATI_QSORTI.o : $(SRCDIR)/OpenATI_QSORTI.f90
	$(F90) -c $(SRCDIR)/OpenATI_QSORTI.f90

.PRECIOUS : $(LOADNAME) ;

clean :
	rm -f *.o

