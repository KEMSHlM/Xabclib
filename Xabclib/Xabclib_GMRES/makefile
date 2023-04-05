ifeq ($(OPENATI_COMP),HITACHI)
  F90 = /opt/ofort90/bin/f90 -Oss -64 -omp -cpp -i,L -fixed -DHF_COMPILER
else
  F90 = ifort -O3 -m64 -fopenmp -shared-intel -cpp -fixed
endif
CC  = cc

SRCDIR   = .

IFILE =  Xabclib_GMRES.o Xabclib_GMRES2.o Xabclib_GMRES_INIT.o Xabclib_GMRES_GS_Policy.o Xabclib_GMRES_MGS_MIN.o  Xabclib_GMRES_MIN.o Xabclib_GMRES_PCHK.o

all : $(IFILE)

.c.o:
	$(CC) -c $<

Xabclib_IDIAGPGEN.o : $(SRCDIR)/Xabclib_IDIAGPGEN.f90
	$(F90) -c $(SRCDIR)/Xabclib_IDIAGPGEN.f90
Xabclib_PCGEN.o : $(SRCDIR)/Xabclib_PCGEN.f90
	$(F90) -c $(SRCDIR)/Xabclib_PCGEN.f90
Xabclib_PCSLV.o : $(SRCDIR)/Xabclib_PCSLV.f90
	$(F90) -c $(SRCDIR)/Xabclib_PCSLV.f90
Xabclib_GMRES.o : $(SRCDIR)/Xabclib_GMRES.f90
	$(F90) -c $(SRCDIR)/Xabclib_GMRES.f90
Xabclib_GMRES2.o : $(SRCDIR)/Xabclib_GMRES2.f90
	$(F90) -c $(SRCDIR)/Xabclib_GMRES2.f90
Xabclib_GMRES_INIT.o : $(SRCDIR)/Xabclib_GMRES_INIT.f90
	$(F90) -c $(SRCDIR)/Xabclib_GMRES_INIT.f90
Xabclib_GMRES_GS_Policy.o : $(SRCDIR)/Xabclib_GMRES_GS_Policy.f90
	$(F90) -c $(SRCDIR)/Xabclib_GMRES_GS_Policy.f90
Xabclib_GMRES_DOT.o : $(SRCDIR)/Xabclib_GMRES_DOT.f90
	$(F90) -c $(SRCDIR)/Xabclib_GMRES_DOT.f90
Xabclib_GMRES_MGS_MIN.o : $(SRCDIR)/Xabclib_GMRES_MGS_MIN.f90
	$(F90) -c $(SRCDIR)/Xabclib_GMRES_MGS_MIN.f90
Xabclib_GMRES_MIN.o : $(SRCDIR)/Xabclib_GMRES_MIN.f90
	$(F90) -c $(SRCDIR)/Xabclib_GMRES_MIN.f90
Xabclib_GMRES_PCHK.o : $(SRCDIR)/Xabclib_GMRES_PCHK.f90
	$(F90) -c $(SRCDIR)/Xabclib_GMRES_PCHK.f90
Xabclib_GMRES_EVAL.o : $(SRCDIR)/Xabclib_GMRES_EVAL.f90
	$(F90) -c $(SRCDIR)/Xabclib_GMRES_EVAL.f90

clean :
	rm -f *.o
