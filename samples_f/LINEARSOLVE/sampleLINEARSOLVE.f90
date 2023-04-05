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
      PARAMETER (MGMAX=100)
      DIMENSION B(NMAX),X(NMAX),EPS(2),IOPT(3)
      DIMENSION IDP(NMAX),
     $          IWK(NMAX),RWK((MGMAX+3)*NMAX+2*(MGMAX+1)*(MGMAX+1))
*
      PARAMETER (NZMAX2=MAX(NZMAX,NMAX*(20*2+1)),
     $           IWPC=NZMAX2*3/2+NMAX*3+50)
      DOUBLE PRECISION RGRPARM(10),RPARM(10),PRECOND(IWPC)
      INTEGER IGRPARM(10),IPARM(10),INFO,IDIAGP(NMAX)
      INTEGER IAT(10)
      INTEGER IATPARAM(50)
      DOUBLE PRECISION RATPARAM(50)
*
      INTEGER IPCPARM(10)
      DOUBLE PRECISION RPCPARM(10)
      DOUBLE PRECISION Y(NMAX),Z(NMAX)
*
      REAL ETIME,TARRAY(2),TIME0,TIME1
      CHARACTER FNAME32 *32
*
      external omp_get_num_threads,omp_get_max_threads
      integer  omp_get_num_threads,omp_get_max_threads
      external OpenATI_DURMV,Xabclib_GMRES,Xabclib_GMRES_EVAL
      external OpenATI_LINEARSOLVE
*
      LDWORK=(MGMAX+3)*NMAX+2*(MGMAX+1)*(MGMAX+1)
*
      MAXP=omp_get_max_threads()
      WRITE(6,*) ' OpenMP Num. of Max Threads =',MAXP
*
      OPEN(5,FILE="../MatrixData/ecl32.dat")
*-----------------------------------------------------------------------
      READ(5,*) N,NZ
      if (n.gt.nmax .or. nz.gt.nzmax) then
          write(6,*) ' not enough array dim.'
          stop
      endif
      READ(5,*) (IRP (I),I=1,N+1)
      READ(5,*) (ICOL(I),I=1,NZ)
      READ(5,*) (VAL (I),I=1,NZ)
      WRITE(6,*) ' N=',N,', NZ=',NZ
      WRITE(6,*) '---------'
*----------------------------------------------------------------------
      DO I=1,N
         X(I)=1.0D0
      END DO
*----------------------------------------------------------------------
      CALL OpenATI_DURMV_11(n,nz,IRP,ICOL,VAL,X,B)
*----------------------------------------------------------------------
      DO I=1,N
         X(I)=0.0D0
      END DO
*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      CALL OpenATI_INIT(IATPARAM, RATPARAM, INFO)
      IATPARAM(50)=0
      CALL OpenATI_LINEARSOLVE(N,NZ,IRP,ICOL,VAL,B,X,
     $                         IATPARAM,RATPARAM,INFO)
*
      ICASE=11
      CALL Xabclib_EVAL(N,NZ,IRP,ICOL,VAL,X,B,RWK,RERR,INFO)
      WRITE(6,*) ' [Policy ] ||b-Ax_i||_2 / ||b||_2    = ',RERR
*
 9999 CONTINUE
C
      STOP
      END
