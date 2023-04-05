      SUBROUTINE Xabclib_GMRES_INIT (N,NNZ,IRP,ICOL,VAL,X,B,R0,V1,BETA,
     &              BETA0,IDIAGP,KIND_PRECOND,MOPE,PRECOND,
     &              IATPARAM,RATPARAM,
     &              ITRCNT,RERR,UINF,LUINF,INFO)
      IMPLICIT DOUBLE PRECISION ( A-H , O-Z )
*     ..
*     .. Scalar Auruments ..
      INTEGER           N, NNZ, KIND_PRECOND, ITRCNT, INFO
      INTEGER           LUINF, NUM_SMP
      DOUBLE PRECISION  BETA, BETA0, RERR
*     ..
*     .. Array Arguments
      INTEGER           IRP(N+1), ICOL(NNZ), IDIAGP(N)
      INTEGER           ICASE
      DOUBLE PRECISION  VAL(NNZ), X(N), B(N), R0(N), V1(N), PRECOND(*)
*
      INTEGER           IATPARAM(50)
      DOUBLE PRECISION  RATPARAM(50)
*
      DOUBLE PRECISION  UINF(LUINF)
*     ..
*
*  =====================================================================
*  Purpose
*  =======
*
*  Xabclib_GMRES_INIT computes the first Arnoldi vector v_1 and 
*  estimates the error of the approximate solution vector.
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
*  X       (input) DOUBLE PRECISION array, dimension (N)
*          The approximate solution vector.
*
*  B       (input) DOUBLE PRECISION array, dimension (N)
*          The right hand side vector.
*
*  R0      (output) DOUBLE PRECISION array, dimension (N)
*          The residual vector. ( R0 = inv(M) (b - A x_0) )
*
*  V1      (output) DOUBLE PRECISION array, dimension (N)
*          The first Arnoldi vector v_1.
*
*  BETA    (output) DOUBLE PRECISION
*          The square norm of the residual vector.
*          ( BETA = || inv(M) R0 ||_2 )
*
*  BETA0   (input/output) DOUBLE PRECISION
*          The square norm of the right hand side vector.
*          ( BETA0 = || inv(M) b ||_2 )
*
*  IDIAGP  (input) INTEGER array, dimension (N)
*          The index of the diagonal elements in VAL.
*
*  KIND_PRECOND  (input) INTEGER
*          ID of the preconditioner used.
*          KIND_PRECOND = 1 : None
*          KIND_PRECOND = 2 : Jacobi preconditioner
*          KIND_PRECOND = 3 : SSOR   preconditioner
*          KIND_PRECOND = 4 : ILU(0) preconditioner
*
*  PRECOND (input) DOUBLE PRECISION array
*          Preconditioner data.
*
*  IPCPARM (input) INTEGER array, dimension (10)
*          INTEGER parameter set for the preconditioner.
*          IPCPARM(1) = 1 : The preconditioner not generated yet.
*          IPCPARM(1) = 2 : The preconditioner already generated.
*          IPCPARM(2:10)  : Not used.
*
*  RPCPARM (input) DOUBLE PRECISION array, dimension (10)
*          DOUBLE PRECISION parameter set for the preconditioner.
*          if KIND_PRECOND = 3,
*          RPCPARM(1) : The omega parameter for SSOR preconditioning.
*          if KIND_PRECOND = 4,
*          RPCPARM(1) : The breakdown criterion of incomplete LU
*                        decomposition process for ILU(0) precondition-
*                        ing.
*
*  ICASE   (input) INTEGER
*          Selected Algorithm NO. of matrix-vector product.
*
*  ITRCNT  (input) INTEGER
*          The iteration count of restarted GMRES.
*
*  RERR    (output) DOUBLE PRECISION
*          The error of the approximate solution vector.
*          ( RERR = || inv(M) R0 ||_2 / || inv(M) b ||_2 )
*
*  UINF    (input) DOUBLE PRECISION array, dimension (LUINF)
*          Setup Information for OpenATI_DURMV
*
*  LUINF   (input) INTEGER
*          The dimension of UINF.
*          if ICASE = 11 :
*             LUINF >= 0
*          else if ICASE = 12 :
*             LUINF >= INT(0.5*NUM_SMP)+1
*          else if ICASE = 13 :
*             LUINF >= INT(1.5*N)+INT(4.25*JL)+2
*          else if ICASE = 21 :
*             LUINF >= INT(1.125*NNZ)+INT(2.125*JL)+1
*
*  NUM_SMP (input) INTEGER
*          number of threads
*
*  INFO    (output) INTEGER
*          Return code.
*
*  =====================================================================
*     ..
*     .. Local variables
*
      INTEGER          I
      DOUBLE PRECISION RBETA
*
*     beta_0 = ||inv(M) b||_2
*     r_0 = inv(M) (b - A x_0); beta = ||r_0||_2; v_1=r_0/beta
*
*
      INFO=0
*   
*  beta_0 = || inv(M) b ||_2
*
      IF (ITRCNT .EQ. 1) THEN
         IF (MOPE.EQ.1) THEN
            IF (KIND_PRECOND.GE.5) then
               CALL ILUTSOL('N',N,NNZ,IRP,ICOL,VAL,R0,B,IDIAGP,
     $                 KIND_PRECOND,PRECOND,IATPARAM,RATPARAM,INFO)
            ELSE
               CALL Xabclib_PCSLV(N,NNZ,IRP,ICOL,VAL,R0,B,IDIAGP,
     $                  PRECOND,IATPARAM,RATPARAM,INFO)
            END IF
         IF (INFO .NE. 0) RETURN
         ELSE
!$omp parallel do
            DO I=1,N
               R0(I)=B(I)
            END DO
!$omp end parallel do
         ENDIF
         BETA0=0.0D0
!$omp parallel do reduction(+:BETA0)
         DO I=1,N
            BETA0=BETA0+R0(I)*R0(I)
         END DO
!$omp end parallel do
         IF (BETA0 .EQ. 0.0D0) THEN
            INFO=-6
            RETURN
         END IF
         BETA0=DSQRT(BETA0)
      END IF
*   
*  r_0 = A x_0
*
      CALL OpenATI_DURMV(N,NNZ,IRP,ICOL,VAL,X,R0,IATPARAM,RATPARAM,
     $        UINF,LUINF,INFO_A)
*   
*  r_0 = b - A x_0 
*
!$omp parallel do
      DO I=1,N
         R0(I)=B(I)-R0(I)
      END DO
!$omp end parallel do
*   
*  r_0 = inv(M) r_0 = inv(M)(b - A x_0)
*
      IF (MOPE.EQ.1) THEN
            IF (KIND_PRECOND.GE.5) then
               CALL ILUTSOL('N',N,NNZ,IRP,ICOL,VAL,R0,R0,IDIAGP,
     $                 KIND_PRECOND,PRECOND,IATPARAM,RATPARAM,INFO)
            ELSE
               CALL Xabclib_PCSLV(N,NNZ,IRP,ICOL,VAL,R0,R0,IDIAGP,
     &                      PRECOND,IATPARAM,RATPARAM,INFO)
            END IF
      IF (INFO .NE. 0) RETURN
      ENDIF
*
*  beta = || r_0 ||_2
*
      BETA=0.0D0
!$omp parallel do reduction(+:BETA)
      DO I=1,N
         BETA=BETA+R0(I)*R0(I)
      END DO
!$omp end parallel do
      IF (INFO .NE. 0) RETURN
      BETA=DSQRT(BETA)
*
*  CALCULATE RELATIVE ERROR
*
      RERR=BETA/BETA0 
*
*  v_1 = r_0 / beta
*
      IF (BETA .NE. 0.0D0) THEN
         RBETA=1.0D0/BETA
!$omp parallel do
         DO I=1,N
            V1(I)=R0(I)*RBETA
         END DO
!$omp end parallel do
      END IF
*
      RETURN
      END
