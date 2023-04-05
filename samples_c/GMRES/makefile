CC = icc -O3 -m64 -shared-intel -fopenmp -lifport -lifcore -shared-intel

SRCDIR = .
OpenATLIB = ../../libOpenAT.a
LOADNAME = gload
LD = $(CC)

IFILE = sampleGMRES.o

$(LOADNAME) : $(IFILE)
	-rm -f $(LOADNAME)
	$(LD) -o $(LOADNAME) $(IFILE) $(OpenATLIB)
	@echo "gload is now up-to-date"

.c.o    :
	$(CC) -c $(INCLUDE) $<

sampleGMRES.o : $(SRCDIR)/sampleGMRES.c
	$(CC) -c $(SRCDIR)/sampleGMRES.c

.PRECIOUS : gload ;

clean :
	rm -f $(LOADNAME)
	rm -f *.o

