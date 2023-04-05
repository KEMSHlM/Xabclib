CC = icc -O3 -m64 -fopenmp -shared-intel -lifport -lifcore -shared-intel

SRCDIR = .
OpenATLIB = ../../libOpenAT.a
LOADNAME = gload
LD = $(CC)

IFILE = sampleLINEARSOLVE.o

$(LOADNAME) : $(IFILE)
	-rm -f $(LOADNAME)
	$(LD) -o $(LOADNAME) $(IFILE) $(OpenATLIB)
	@echo "gload is now up-to-date"

.c.o    :
	$(CC) -c $(INCLUDE) $<

sampleLINEARSOLVE.o : $(SRCDIR)/sampleLINEARSOLVE.c
	$(CC) -c $(SRCDIR)/sampleLINEARSOLVE.c

.PRECIOUS : gload ;

clean :
	rm -f $(LOADNAME)
	rm -f *.o

