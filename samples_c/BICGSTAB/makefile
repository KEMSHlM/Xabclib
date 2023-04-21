CC = icc -O3 -m64 -fopenmp -shared-intel -lifport -lifcore -shared-intel 

SRCDIR = .
# OpenATLIB = ../../build/libXabClib.a
OpenATLIB = ../../libOpenAT.a
LOADNAME = gload
LD = $(CC)

IFILE = sampleBICGSTAB.o

$(LOADNAME) : $(IFILE)
	-rm -f $(LOADNAME)
	$(LD) -o $(LOADNAME) $(IFILE) $(OpenATLIB)
	@echo "gload is now up-to-date"

.c.o    :
	$(CC) -c $(INCLUDE) $<

sampleBICGSTAB.o : $(SRCDIR)/sampleBICGSTAB.c
	$(CC) -c $(SRCDIR)/sampleBICGSTAB.c

.PRECIOUS : gload ;

clean :
	rm -f $(LOADNAME)
	rm -f *.o

