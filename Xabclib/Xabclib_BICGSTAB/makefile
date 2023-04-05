ifeq ($(OPENATI_COMP),HITACHI)
  F90 = /opt/ofort90/bin/f90 -Oss -64 -omp -cpp -i,L -fixed -DHF_COMPILER
else
	# -shared-intel ってのは，静的に確保されるメモリの上限を大きくするものらしい．
	# intel コンパイラにはなかったから，代わりが-shared-intelみたい．
  F90 = ifort -O3 -m64 -fopenmp -shared-intel -cpp -fixed
#  F90 = ifort -O3 -m64 -openmp -shared-intel -cpp -fixed
#  F90 = ifort -O3 -m64 -parallel -shared-intel -cpp -fixed
endif
CC  = cc

SRCDIR    = .


IFILE = Xabclib_BLAS1.o Xabclib_BICGSTAB.o Xabclib_BICGSTABV.o Xabclib_PCGENE.o Xabclib_BICGSTAB_PCHK.o ILUT.o ILUTSOL.o Xabclib_PCSLV.o Xabclib_EVAL.o Xabclib_IDIAGPGEN.o


all : $(IFILE)

.c.o	:
	$(CC) -c $(INCLUDE) $<

Xabclib_BICGSTABV.o : $(SRCDIR)/Xabclib_BICGSTABV.f90
	$(F90) -c $(SRCDIR)/Xabclib_BICGSTABV.f90

Xabclib_BICGSTAB.o : $(SRCDIR)/Xabclib_BICGSTAB.f90
	$(F90) -c $(SRCDIR)/Xabclib_BICGSTAB.f90

Xabclib_BLAS1.o : $(SRCDIR)/Xabclib_BLAS1.f90
	$(F90) -c $(SRCDIR)/Xabclib_BLAS1.f90
Xabclib_BICGSTAB_PCHK.o : $(SRCDIR)/Xabclib_BICGSTAB_PCHK.f90
	$(F90) -c $(SRCDIR)/Xabclib_BICGSTAB_PCHK.f90
Xabclib_PCGENE.o : $(SRCDIR)/Xabclib_PCGENE.f90
	$(F90) -c $(SRCDIR)/Xabclib_PCGENE.f90
Xabclib_PCSLV.o : $(SRCDIR)/Xabclib_PCSLV.f90
	$(F90) -c $(SRCDIR)/Xabclib_PCSLV.f90
Xabclib_EVAL.o : $(SRCDIR)/Xabclib_EVAL.f90
	$(F90) -c $(SRCDIR)/Xabclib_EVAL.f90
Xabclib_IDIAGPGEN.o : $(SRCDIR)/Xabclib_IDIAGPGEN.f90
	$(F90) -c $(SRCDIR)/Xabclib_IDIAGPGEN.f90

ILUT.o : $(SRCDIR)/ILUT.f90
	$(F90) -c $(SRCDIR)/ILUT.f90
ILUTSOL.o : $(SRCDIR)/ILUTSOL.f90
	$(F90) -c $(SRCDIR)/ILUTSOL.f90



clean :
	rm -f *.o
