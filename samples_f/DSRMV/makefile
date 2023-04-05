ifeq ($(OPENATI_COMP),HITACHI)
  F90 = /opt/ofort90/bin/f90 -Os -64 -omp -cpp -i,L -fixed -DHF_COMPILER
else
  F90 = ifort -O3 -m64 -fopenmp -shared-intel -cpp -fixed
endif
SRCDIR       = .
LIB = ../../libOpenAT.a

LOADNAME = lload
LD    = $(F90)

IFILE = sampleDSRMV.o


$(LOADNAME) : $(IFILE)
	-rm -f $(LOADNAME)
	$(LD) -o $(LOADNAME) $(IFILE) $(LIB)
	@echo "$(LOADNAME) is now up-to-date"

sampleDSRMV.o : $(SRCDIR)/sampleDSRMV.f90
	$(F90) -c $(SRCDIR)/sampleDSRMV.f90

.PRECIOUS : $(LOADNAME) ;

clean :
	rm -f $(LOADNAME)
	rm -f *.o
