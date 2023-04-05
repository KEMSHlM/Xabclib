      SUBROUTINE Xabclib_PCSLV_SYM(N,NNZ,IRP,ICOL,VAL,
     $                             X,Y,PRECOND,IATPARAM,RATPARAM,WK,
     $                             SINF,LSINF,INFO)
      IMPLICIT DOUBLE PRECISION ( A-H , O-Z )
*     ..
*     .. Scalar Auruments ..
      INTEGER           N, NNZ, LSINF, INFO
*     ..
*     .. Array Arguments
      INTEGER           IRP(N+1), ICOL(NNZ)
      DOUBLE PRECISION  VAL(NNZ), X(N), Y(N), PRECOND(*), SINF(LSINF)
      DOUBLE PRECISION  WK(N,*)
*
      INTEGER           IATPARAM(50)
      DOUBLE PRECISION  RATPARAM(50)
*     ..
*
*  =====================================================================
*  Purpose
*  =======
*
*  Xabclib_PCSLV_SYM solves Mx = y for x. (x = inv(M)*y)
*     M : the preconditioner used.
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          The order of the matrix.
*
*  NNZ     (input) INTEGER
*          Number of nonzero elements in the matrix.
*
*  IRP     (input) INTEGER array, dimension (N+1)
*          The index of the elements in VAL which start a row of the
*          upper traiangler matrix.
*
*  ICOL    (input) INTEGER array, dimension (NNZ)
*          The column indices of each nonzero.
*
*  VAL     (input) DOUBLE PRECISION array, dimension (NNZ)
*          Nonzero elements of the matrix in Compressed Row Storage
*          format.
*
*  X       (output) DOUBLE PRECISION array, dimension (N)
*          The solution vector of Mx = y.
*
*  Y       (input) DOUBLE PRECISION array, dimension (N)
*          The right hand side vector of Mx = y.
*
*  PRECOND (input) DOUBLE PRECISION array
*          Preconditioner data.
*
*  IATPARAM (input) Integer array, dimension (50)
*          Integer parameter set for OpenATI and Xabclib.
*           IATPARAM(25)    : Preconditioner type.
*
*  RATPARAM (input) DOUBLE PRECISION array, dimension (50)
*          Double precision parameter set for OpenATI and Xabclib.
*           RATPARAM(25)    : Preconditioner parameter.
*
*  WK(N,   (Workspase) DOUBLE PRECISION array
*  NUM_SMP)Workspace for parallel vector reduction array.
*          (NUM_SMP = IATPARAM(3))
*
*  SINF    (input/output) DOUBLE PRECISION array
*  (LSINF) Setup Information for Matrix-Vector product Impl.
*
*  LSINF   (input) INTEGER
*          The size of SINF.
*          case IATPARAM(7)=1 :
*             if IATPARAM(8) = 11 :
*                LSINF >= 0
*             else if IATPARAM(8) = 12 :
*                LSINF >= INT(0.5*NUM_SMP)+1
*             else if IATPARAM(8) = 13 :
*                LSINF >= N+NUM_SMP+3
*          case IATPARAM(7)=2 :
*             LSINF >= INT(0.5*NUM_SMP)+1
*          case IATPARAM(7)=3 :
*             LSINF >= N+NUM_SMP+3
*          (NUM_SMP = IATPARAM(3))
*
*  INFO    (output) INTEGER
*          Return code.
*
*  =====================================================================
*     ..
*     .. Local variables
*
      INTEGER          I, KIND_PRECOND
      DOUBLE PRECISION OMEGA
*
      INFO=0
*
      KIND_PRECOND = IATPARAM(25)
*
      IF (KIND_PRECOND .EQ. 1 .OR. KIND_PRECOND .EQ. 2) THEN
*
*-- PRECONDITIONER NOT APPLIED
*
!$omp parallel do
         DO I=1,N
            X(I)=Y(I)
         END DO
!$omp end parallel do
*
      ELSE IF (KIND_PRECOND .EQ. 3) THEN
*
*-- JACOBI ITERATED PRECONDITIONER
*
         ITER = IATPARAM(26)
         CALL Xabclib_PCSLV1_SYM(N,NNZ,IRP,ICOL,VAL,X,Y,PRECOND,ITER,
     &                           IATPARAM,RATPARAM,WK,SINF,LSINF,INFO)
         RATPARAM(30)=RATPARAM(30)+ITER*2.0D0*(2.0D0*NNZ-N)
      ELSE IF (KIND_PRECOND .EQ. 4) THEN
*
*-- ICCG PRECONDITIONER
*
         CALL Xabclib_PCSLV4_SYM(N,NNZ,IRP,ICOL,VAL,X,Y,PRECOND,INFO)
         RATPARAM(30)=RATPARAM(30)+2.0D0*(2.0D0*NNZ-N)
      ELSE 
*
*-- KIND_PRECOND > 4 NOT SUPPORTED
*
         INFO=-8
      END IF
*
      RETURN
      END
*-----------------------------------------------------------------------
      SUBROUTINE Xabclib_PCSLV1_SYM(N,NNZ,IRP,ICOL,VAL,X,Y,PIVOTS,ITER,
     &                              IATPARAM,RATPARAM,WK,SINF,LSINF,
     &                              INFO)
      IMPLICIT DOUBLE PRECISION ( A-H , O-Z )
*     ..
*     .. Scalar Auruments ..
      INTEGER           N, NNZ, LSINF, INFO
*     ..
*     .. Array Arguments
      INTEGER           IRP(N+1), ICOL(NNZ)
      DOUBLE PRECISION  VAL(NNZ), X(N), Y(N), PIVOTS(N), SINF(LSINF)
      DOUBLE PRECISION  WK(N,*)
*
      INTEGER           IATPARAM(50)
      DOUBLE PRECISION  RATPARAM(50)
*     ..
*
*  =====================================================================
*  Purpose
*  =======
*
*  Xabclib_PCSLV1 solves Mx = y for x. (x = inv(M)*y)
*     M : the Jacobi preconditioner.
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          The order of the matrix.
*
*  NNZ     (input) INTEGER
*          Number of nonzero elements in the matrix.
*
*  IRP     (input) INTEGER array, dimension (N+1)
*          The index of the elements in VAL which start a row of the
*          matrix.
*
*  ICOL    (input) INTEGER array, dimension (NNZ)
*          The column indices of each nonzero.
*
*  VAL     (input) DOUBLE PRECISION array, dimension (NNZ)
*          Nonzero elements of the matrix in Compressed Row Storage
*          format.
*
*  X       (output) DOUBLE PRECISION array, dimension (N)
*          The solution vector of Mx = y.
*
*  Y       (input) DOUBLE PRECISION array, dimension (N)
*          The right hand side vector of Mx = y.
*
*  PIVOTS  (input) DOUBLE PRECISION array, dimension (N)
*          Jacobi preconditioner data (RECIPROCAL DIAGONAL).
*
*  ITER    (input) INTEGER
*           Iterations of the Jacobi refinement.
*
*  IATPARAM (input/output) INTEGER array, dimension (50)
*          Integer parameter set for OpenATI and Xabclib.
*          On entry,
*           IATPARAM(3)     : Number of threads.
*           IATPARAM(8)     : OpenATI_DSRMV impl. method.
*           IATPARAM(50)    : Debug print control flag.
*          On exit,
*
*  RATPARAM (input/output) DOUBLE PRECISION array, dimension (50)
*          Double precision parameter set for OpenATI and Xabclib.
*
*  WK(N,   (Workspase) DOUBLE PRECISION array
*  NUM_SMP)Workspace for parallel vector reduction array.
*          (NUM_SMP = IATPARAM(3))
*
*  SINF    (input/output) DOUBLE PRECISION array
*  (LSINF) Setup Information for Matrix-Vector product Impl.
*
*  LSINF   (input) INTEGER
*          The size of SINF.
*          case IATPARAM(7)=1 :
*             if IATPARAM(8) = 11 :
*                LSINF >= 0
*             else if IATPARAM(8) = 12 :
*                LSINF >= INT(0.5*NUM_SMP)+1
*             else if IATPARAM(8) = 13 :
*                LSINF >= N+NUM_SMP+3
*          case IATPARAM(7)=2 :
*             LSINF >= INT(0.5*NUM_SMP)+1
*          case IATPARAM(7)=3 :
*             LSINF >= N+NUM_SMP+3
*          (NUM_SMP = IATPARAM(3))
*
*  INFO    (output) INTEGER
*          Return code.
*
*  =====================================================================
*     ..
*     .. Local variables
*
      INTEGER I
      DOUBLE PRECISION ,ALLOCATABLE :: LWK(:)
*
      INFO=0
      ALLOCATE(LWK(N),STAT=IERR_ALLOC)
      IF (IERR_ALLOC .NE.0 ) THEN
         INFO=-500
         WRITE(IO,*) '   <<<< ALLOCATBLE ERROR CODE >>>> = ',INFO
         RETURN
      END IF
*
      DO I=1,N
         X(I)=PIVOTS(I)*Y(I)
      ENDDO
*
      DO K=1,ITER
*
*     ----Mat*Vec except diagonal
      CALL Xabclib_DSRMV_ForJacobiIter(N,NNZ,IRP,ICOL,VAL,X,LWK,
     $                                 IATPARAM,RATPARAM,WK,SINF,
     $                                 LSINF,INFO)
**
!$omp parallel do
      DO I=1,N
         X(I)=PIVOTS(I)*(Y(I)-LWK(I))
      END DO
!$omp end parallel do

*
      ENDDO
*
*
      DEALLOCATE(LWK,STAT=IERR_DEALLOC)
      IF (IERR_DEALLOC.NE.0) THEN
         INFO=500
         WRITE(6,*) '[Error] IERR_DEALLOC=',IERR_DEALLOC
      END IF
*
      RETURN
      END
*-----------------------------------------------------------------------
      SUBROUTINE Xabclib_PCSLV4_SYM(N,NNZ,IRP,ICOL,VAL,
     $                              X,Y,PIVOTS,INFO)
      IMPLICIT DOUBLE PRECISION ( A-H , O-Z )
*     ..
*     .. Scalar Auruments ..
      INTEGER           N, NNZ, INFO
*     ..
*     .. Array Arguments
      INTEGER           IRP(N+1), ICOL(NNZ)
      DOUBLE PRECISION  VAL(NNZ), X(N), Y(N), PIVOTS(N)
*     ..
*
*  =====================================================================
*  Purpose
*  =======
*
*  Xabclib_PCSLV4 solves Mx = y for x. (x = inv(M)*y)
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          The order of the matrix.
*
*  NNZ     (input) INTEGER
*          Number of nonzero elements in the matrix.
*
*  IRP     (input) INTEGER array, dimension (N+1)
*          The index of the elements in VAL which start a row of the
*          upper traiangler matrix.
*
*  ICOL    (input) INTEGER array, dimension (NNZ)
*          The column indices of each nonzero.
*
*  VAL     (input) DOUBLE PRECISION array, dimension (NNZ)
*          Nonzero elements of the matrix in Compressed Row Storage
*          format.
*
*  X       (output) DOUBLE PRECISION array, dimension (N)
*          The solution vector of Mx = y.
*
*  Y       (input) DOUBLE PRECISION array, dimension (N)
*          The right hand side vector of Mx = y.
*
*  PIVOTS  (input) DOUBLE PRECISION array, dimension (N)
*          ILU(0) preconditioner data.
*
*  INFO    (output) INTEGER
*          Return code.
*
*  =====================================================================
*     ..
*     .. Local variables
*
      INTEGER          I,J,IP
      DOUBLE PRECISION SUM,SUM_ARRAY(N)
*
      INFO=0
*
*     forward solve
*
      SUM_ARRAY(1:N) = 0.0D0
      DO I=1,N
         X(I) = PIVOTS(I) * (Y(I) - SUM_ARRAY(I))
         DO J = IRP(I)+1, IRP(I+1)-1
            IP = ICOL(J)
            SUM_ARRAY(IP) = SUM_ARRAY(IP) + VAL(J) * X(I)
         END DO
      END DO
*
*     backward solve
*
      DO I=N,1,-1
         SUM = 0.0D0
         DO J = IRP(I)+1, IRP(I+1)-1
            SUM = SUM + VAL(J) * X(ICOL(J))
         END DO
         X(I) = X(I) - PIVOTS(I) * SUM
      END DO
*
      RETURN
      END
