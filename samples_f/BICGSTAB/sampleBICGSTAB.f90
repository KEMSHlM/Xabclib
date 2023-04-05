      PROGRAM MAIN
*
*-----------------------------------------------------------------------
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
      INTEGER NZMAX, NMAX
      PARAMETER (NZMAX = 10000000, NMAX = 300000)
      INTEGER  N, NZ
      INTEGER IRP (NMAX), ICOL (NZMAX)
      DOUBLE PRECISION VAL (NZMAX)
*
      PARAMETER (MGMAX=100,NIPC=1,NRPC=NMAX,NIWK=1)
      DIMENSION B(NMAX),X(NMAX),EPS(2),IOPT(3)
      DIMENSION IDP(NMAX),
     $          IWK(NMAX),RWK(9*NMAX+(NMAX-1)/2+1)
*
      PARAMETER (NZMAX2=NMAX*(20*2+1),IWPC=NZMAX2*3/2+NMAX*3+20)
      DOUBLE PRECISION RGRPARM(10),RPARM(10),PRECOND(IWPC)
      INTEGER IGRPARM(10),IPARM(10),INFO,IDIAGP(NMAX)
      INTEGER IAT(10),IATPARAM(50)
      DOUBLE PRECISION RATPARAM(50)
*
      INTEGER IPCPARM(10)
      DOUBLE PRECISION RPCPARM(10)
      DOUBLE PRECISION Y(NMAX),Z(NMAX)
*
      REAL ETIME,TARRAY(2),TIME0,TIME1
      CHARACTER*60 FNAME
      CHARACTER*20 CVER
*
      external omp_get_num_threads,omp_get_max_threads
      integer  omp_get_num_threads,omp_get_max_threads
*
      LDWORK=9*NMAX+(NMAX-1)/2+1+1
*
      MAXP=omp_get_max_threads()
      WRITE(6,*) ' OpenMP Num. of Max Threads =',MAXP
*
      FNAME="../MatrixData/ecl32.dat"
*-----------------------------------------------------------------------
      OPEN(5,FILE=FNAME)
      READ(5,*) N,NZ
      READ(5,*) (IRP (I),I=1,N+1)
      READ(5,*) (ICOL(I),I=1,NZ)
      READ(5,*) (VAL (I),I=1,NZ)
      CLOSE(5)
      WRITE(6,*) ' FILE=',FNAME
      WRITE(6,*) ' N=',N,', NZ=',NZ
      WRITE(6,*) '---------'
*----------------------------------------------------------------------
      DO I=1,N
         X(I)=1.0D0
      END DO
*----------------------------------------------------------------------
      ICASE=11
      CALL OpenATI_DURMV_11(n,nz,IRP,ICOL,VAL,X,B)
*----------------------------------------------------------------------
      DO I=1,N
         X(I)=0.0D0
      END DO
*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      CALL OPENATI_INIT(IATPARAM,RATPARAM,INFO)
*-----------ILU(0)
      IATPARAM(25)=5
*
      IDEBUG=0
      IATPARAM(50) = IDEBUG
      IATPARAM(22) = 10000
      RATPARAM(22) = 600.0D0
*
      WRITE(6,*) ''
      WRITE(6,*) '=================================================='
      CALL OpenATI_GET_VER(CVER)
      WRITE(6,*) '    OpenATI_VER IS ','" ',CVER,' "'
      WRITE(6,*) '=================================================='
*
*------------Shoji Itoh version
      CALL Xabclib_BICGSTAB(N,NZ,IRP,ICOL,VAL,B,X,
     $                      PRECOND,IWPC,IATPARAM,RATPARAM,
     $                      RWK,LDWORK,INFO)
      IF (INFO .EQ. 0) THEN
         WRITE(6,*) ' ## BiCGSTAB PROCESS SUCCESSFULLY ENDED.'
      ELSE
         WRITE(6,*) ' ## BiCGSTAB PROCESS FAILED. INFO=',INFO
      END IF
*================
*---------------------------------------------------------------------
      WRITE(6,*)' >>>>> Selected Mat*Vec Impl.=',IATPARAM(10)
      IF(IATPARAM(4).EQ. 1) THEN
         WRITE(6,*)' START KRYLOV SUBSPACE DIMENSION =',IATPARAM(28)
         WRITE(6,*)' FINAL KRYLOV SUBSPACE DIMENSION =',IATPARAM(29)
      END IF
      WRITE(6,6020) IATPARAM(23),RATPARAM(29)
 6020 FORMAT('  NUM OF ITERS =',I4,', RESID. =',1PD22.15)
      WRITE(6,6021) RATPARAM(31),RATPARAM(32)
 6021 FORMAT('  PRECONDITIONING TIME =',1PD12.5,' (SEC) , ETIME =',
     $       1PD12.5,'(SEC)')
      WRITE(6,6022) RATPARAM(28)
 6022 FORMAT('  RHS norm =',1PD12.5)
      WRITE(6,6023) RATPARAM(30)
 6023 FORMAT('  Floating Ope. counts =',1PD12.5)
      WRITE(6,*) '===================================================='
*
      CALL Xabclib_EVAL(N,NZ,IRP,ICOL,VAL,X,B,RWK,RERR,INFO)
      WRITE(6,*) ' ||b-Ax_i||_2 / ||b||_2    = ',RERR
*
*---------------------------------------------------------------------
 9999 CONTINUE
C
      STOP
      END
