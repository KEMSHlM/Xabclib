ifeq ($(OPENATI_COMP),HITACHI)
  F90 = /opt/ofort90/bin/f90 -Os -64 -omp -cpp -i,L -fixed -DHF_COMPILER
else
  F90 = ifort -O3 -m64 -fopenmp -shared-intel -cpp -fixed
endif

SRCDIR       = .
OpenATLIB = ../../libOpenAT.a
LOADNAME = lload
LD    = $(F90)


IFILE = sampleEIGENSOLVE_TS.o


$(LOADNAME) : $(IFILE)
	-rm -f $(LOADNAME)
	$(LD) -o $(LOADNAME) $(IFILE) $(OpenATLIB) 
	@echo "lload is now up-to-date"

.c.o	:
	$(CC) -c $(INCLUDE) $<

sampleEIGENSOLVE_TS.o : $(SRCDIR)/sampleEIGENSOLVE_TS.f90
	$(F90) -c $(SRCDIR)/sampleEIGENSOLVE_TS.f90

.PRECIOUS : lload ;

clean :
	rm -f $(LOADNAME)
	rm -f *.o
