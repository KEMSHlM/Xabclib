ifeq ($(OPENATI_COMP),HITACHI)
  F90 = /opt/ofort90/bin/f90 -Os -64 -omp -cpp -i,L -fixed -DHF_COMPILER
else
  F90 = ifort -O3 -m64 -fopenmp -shared-intel -cpp -fixed
endif
CC  = cc

SRCDIR    = .
OpenATLIB = ../../libOpenAT.a
LOADNAME  = aload
LD        = $(F90)

IFILE = sampleARNOLDI.o



$(LOADNAME) : $(IFILE)
	-rm -f $(LOADNAME)
	$(LD) -o $(LOADNAME) $(IFILE) $(OpenATLIB)
	@echo "aload is now up-to-date"

.c.o	:
	$(CC) -c $(INCLUDE) $<

sampleARNOLDI.o : $(SRCDIR)/sampleARNOLDI.f90
	$(F90) -c $(SRCDIR)/sampleARNOLDI.f90

.PRECIOUS : aload ;

clean :
	rm -f $(LOADNAME)
	rm -f *.o
