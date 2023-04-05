      PROGRAM MAIN
*
*-----------------------------------------------------------------------
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
      INTEGER NTMP, NZTMP
*
      INTEGER IRPTMP, ICOLTMP
      ALLOCATABLE :: IRPTMP(:), ICOLTMP(:)
      DOUBLE PRECISION VALTMP
      ALLOCATABLE :: VALTMP(:)
*
      INTEGER N, NZ, INFO
*
      INTEGER IRP, ICOL, IATPARAM
      ALLOCATABLE :: IRP(:), ICOL(:), IATPARAM(:)
*
      DOUBLE PRECISION VAL, B, X, RATPARAM, RWK, RERR
      ALLOCATABLE :: VAL(:), B(:), X(:), RATPARAM(:), RWK(:)
*
      CHARACTER FNAME32 *32
*
      external omp_get_num_threads,omp_get_max_threads
      integer  omp_get_num_threads,omp_get_max_threads
      external omp_get_thread_num,omp_set_thread_num
      integer  omp_get_thread_num,omp_set_thread_num
*
      MAXP=omp_get_max_threads()
      WRITE(6,*) ' OpenMP Num. of Max Threads =',MAXP
*
      OPEN(5,FILE="../MatrixData/ecl32.dat")
*-----------------------------------------------------------------------
      READ(5,*) NTMP, NZTMP
      ALLOCATE(IRPTMP(NTMP+1),ICOLTMP(NZTMP), VALTMP(NZTMP))
      READ(5,*) (IRPTMP (I),I=1,NTMP+1)
      READ(5,*) (ICOLTMP(I),I=1,NZTMP)
      READ(5,*) (VALTMP (I),I=1,NZTMP)
      CLOSE(5)
      WRITE(6,*) ' N=',NTMP,', NZ=',NZTMP
      WRITE(6,*) '---------'
*----------------------------------------------------------------------
*
!$omp parallel default(none)
!$omp+ private(N, NZ, IRP, ICOL, VAL, B, X, INFO)
!$omp+ private(IATPARAM, RATPARAM, IP, IO)
!$omp+ private(RWK, RERR)
!$omp+ shared(NTMP, NZTMP, IRPTMP, ICOLTMP, VALTMP)
*
      IP=OMP_GET_THREAD_NUM()
*
      N=NTMP
      NZ=NZTMP
*
      ALLOCATE(IRP(N+1),ICOL(NZ), VAL(NZ))
      ALLOCATE(B(N), X(N), RWK(N))
      ALLOCATE(IATPARAM(50), RATPARAM(50))
*
      DO I=1,N+1
         IRP(I)=IRPTMP(I)
      END DO
      DO I=1,NZ
         ICOL(I)=ICOLTMP(I)
         VAL(I)=VALTMP(I)
      END DO
      DO I=1,N
         X(I)=1.0D0
      END DO
*----------------------------------------------------------------------
      CALL OpenATI_DURMV_11(N,NZ,IRP,ICOL,VAL,X,B)
*----------------------------------------------------------------------
      DO I=1,N
         X(I)=0.0D0
      END DO
*
      CALL OpenATI_INIT(IATPARAM, RATPARAM, INFO)
!$omp barrier
      write(6,*)'*** OpenATI_LINEARSOLVE THREAD-SAFE TEST ***',IP
C      IATPARAM(50)=1
      CALL OpenATI_LINEARSOLVE(N,NZ,IRP,ICOL,VAL,B,X,
     $                         IATPARAM,RATPARAM,INFO)
*
      CALL Xabclib_EVAL(N,NZ,IRP,ICOL,VAL,X,B,RWK,RERR,INFO)
      WRITE(6,*) ' IP = ',IP, ' [Policy ] ||b-Ax_i||_2 / ||b||_2    = '
     $           ,RERR
*
      DEALLOCATE(IRP, ICOL, VAL, B, X, RWK, IATPARAM, RATPARAM)
!$omp barrier
!$omp end parallel
      DEALLOCATE(IRPTMP, ICOLTMP, VALTMP)
 9999 CONTINUE
C
      STOP
      END
