ifeq ($(OPENATI_COMP),HITACHI)
  F90 = /opt/ofort90/bin/f90 -Os -64 -omp -cpp -i,L -fixed -DHF_COMPILER
else
  F90 = ifort -O3 -m64 -fopenmp -shared-intel -cpp -fixed
endif

SRCDIR       = .
OpenATLIB = ../../libOpenAT.a
LOADNAME = gload
LD    = $(F90) -fopenmp

IFILE = sampleLINEARSOLVE_TS.o


$(LOADNAME) : $(IFILE)
	-rm -f $(LOADNAME)
	$(LD) -o $(LOADNAME) $(IFILE) $(OpenATLIB)
	@echo "gload is now up-to-date"

.c.o	:
	$(CC) -c $(INCLUDE) $<

sampleLINEARSOLVE_TS.o : $(SRCDIR)/sampleLINEARSOLVE_TS.f90
	$(F90) -c $(SRCDIR)/sampleLINEARSOLVE_TS.f90

.PRECIOUS : gload ;

clean :
	rm -f $(LOADNAME)
	rm -f *.o
