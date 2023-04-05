      SUBROUTINE Xabclib_PCSLV (N,NNZ,IRP,ICOL,VAL,X,Y,IDIAGP,
     $                          PRECOND,IATPARAM,RATPARAM,INFO)
      IMPLICIT DOUBLE PRECISION ( A-H , O-Z )
*     ..
*     .. Scalar Auruments ..
      INTEGER           N, NNZ, KIND_PRECOND, INFO
*     ..
*     .. Array Arguments
      INTEGER           IRP(N+1), ICOL(NNZ), IDIAGP(N)
      DOUBLE PRECISION  VAL(NNZ), X(N), Y(N), PRECOND(*)
*
      INTEGER           IATPARAM(50)
      DOUBLE PRECISION  RATPARAM(50)
*     ..
*
*  =====================================================================
*  Purpose
*  =======
*
*  Xabclib_PCSLV solves Mx = y for x. (x = inv(M)*y)
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
*  IDIAGP  (input) INTEGER array, dimension (N)
*          The index of the diagonal elements in VAL.
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
*  INFO    (output) INTEGER
*          Return code.
*
*  =====================================================================
*     ..
*     .. Local variables
*
      INTEGER          I
      DOUBLE PRECISION OMEGA
*
      INFO=0
*
      KIND_PRECOND = IATPARAM(25)
*
      IF (KIND_PRECOND .EQ. 1) THEN
*
*-- PRECONDITIONER NOT APPLIED
*
!$omp parallel do
         DO I=1,N
            X(I)=Y(I)
         END DO
!$omp end parallel do
*
      ELSE IF (KIND_PRECOND .EQ. 2) THEN
*
*-- JACOBI PRECONDITIONER
*
         CALL Xabclib_PCSLV1(N,NNZ,IRP,ICOL,VAL,X,Y,IDIAGP,PRECOND,INFO)
      ELSE IF (KIND_PRECOND .EQ. 3) THEN
*
*-- SSOR PRECONDITIONER
*
         OMEGA = RATPARAM(25)
         CALL Xabclib_PCSLV2(N,NNZ,IRP,ICOL,VAL,X,Y,IDIAGP,PRECOND,
     $           OMEGA,INFO)
      ELSE IF (KIND_PRECOND .EQ. 4) THEN
*
*-- ILU(0) PRECONDITIONER
*
         CALL Xabclib_PCSLV3(N,NNZ,IRP,ICOL,VAL,X,Y,IDIAGP,PRECOND,
     $           INFO)
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
      SUBROUTINE Xabclib_PCSLV1(N,NNZ,IRP,ICOL,VAL,X,Y,IDIAGP,PIVOTS,
     $              INFO)
      IMPLICIT DOUBLE PRECISION ( A-H , O-Z )
*     ..
*     .. Scalar Auruments ..
      INTEGER           N, NNZ, INFO
*     ..
*     .. Array Arguments
      INTEGER           IRP(N+1), ICOL(NNZ), IDIAGP(N)
      DOUBLE PRECISION  VAL(NNZ), X(N), Y(N), PIVOTS(N)
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
*  IDIAGP  (input) INTEGER array, dimension (N)
*          The index of the diagonal elements in VAL.
*
*  PIVOTS  (input) DOUBLE PRECISION array, dimension (N)
*          Jacobi preconditioner data.
*
*  INFO    (output) INTEGER
*          Return code.
*
*  =====================================================================
*     ..
*     .. Local variables
*
      INTEGER I
*
      INFO=0
*
!$omp parallel do
      DO I=1,N
         X(I)=PIVOTS(I)*Y(I)
      END DO
!$omp end parallel do
*
      RETURN
      END
*-----------------------------------------------------------------------
      SUBROUTINE Xabclib_PCSLV2(N,NNZ,IRP,ICOL,VAL,X,Y,IDIAGP,PIVOTS,
     $              OMEGA,INFO)
      IMPLICIT DOUBLE PRECISION ( A-H , O-Z )
*     ..
*     .. Scalar Auruments ..
      INTEGER           N, NNZ, INFO
      DOUBLE PRECISION  OMEGA
*     ..
*     .. Array Arguments
      INTEGER           IRP(N+1), ICOL(NNZ), IDIAGP(N)
      DOUBLE PRECISION  VAL(NNZ), X(N), Y(N), PIVOTS(N)
*     ..
*
*  =====================================================================
*  Purpose
*  =======
*
*  Xabclib_PCSLV2 solves Mx = y for x. (x = inv(M)*y)
*     M : the SSOR preconditioner.
*     M = w/(2-w) (1/w D + L) inv(D) (1/w D + U), w : omega
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
*  IDIAGP  (input) INTEGER array, dimension (N)
*          The index of the diagonal elements in VAL.
*
*  PIVOTS  (input) DOUBLE PRECISION array, dimension (N)
*          SSOR preconditioner data.
*
*  OMEGA   (input) DOUBLE PRECISION
*          The omega parameter for SSOR preconditioning.
*
*  INFO    (output) INTEGER
*          Return code.
*
*  =====================================================================
*     ..
*     .. Local variables
*
      INTEGER          I,J,K,J_PTR
      DOUBLE PRECISION W2,SUM
*
      INFO=0
      W2=(2.0D0-OMEGA)/OMEGA
*
*     SOLVE w/(2-w) (1/w D + L) z = y FOR z.
*
      DO I=1,N
         SUM=0.0D0
         DO J_PTR=IRP(I),IDIAGP(I)-1
            SUM=SUM+VAL(J_PTR)*X(ICOL(J_PTR))
         END DO
         X(I)=OMEGA*PIVOTS(I)*(W2*Y(I)-SUM)
      END DO
*
*     SOLVE (1/w I + inv(D) U) x = z FOR x.
*
      DO I=N,1,-1
         SUM=0.0D0
         DO J_PTR=IDIAGP(I)+1,IRP(I+1)-1
            SUM=SUM+VAL(J_PTR)*X(ICOL(J_PTR))
         END DO
         X(I)=OMEGA*(X(I)-PIVOTS(I)*SUM)
      END DO
*
      RETURN
      END
*-----------------------------------------------------------------------
      SUBROUTINE Xabclib_PCSLV3(N,NNZ,IRP,ICOL,VAL,X,Y,IDIAGP,PIVOTS,
     $              INFO)
      IMPLICIT DOUBLE PRECISION ( A-H , O-Z )
*     ..
*     .. Scalar Auruments ..
      INTEGER           N, NNZ, INFO
*     ..
*     .. Array Arguments
      INTEGER           IRP(N+1), ICOL(NNZ), IDIAGP(N)
      DOUBLE PRECISION  VAL(NNZ), X(N), Y(N), PIVOTS(N)
*     ..
*
*  =====================================================================
*  Purpose
*  =======
*
*  Xabclib_PCSLV3 solves Mx = y for x. (x = inv(M)*y)
*     M : the ILU(0) preconditioner.
*     M = (D + L) inv(D) (D + U)
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
*  IDIAGP  (input) INTEGER array, dimension (N)
*          The index of the diagonal elements in VAL.
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
      INTEGER          I,J,K,J_PTR
      DOUBLE PRECISION SUM
*
      INFO=0
*
*     SOLVE (D + L_A) z = y FOR z.
*
      DO I=1,N
         SUM=0.0D0
         DO J_PTR=IRP(I),IDIAGP(I)-1
            SUM=SUM+VAL(J_PTR)*X(ICOL(J_PTR))
         END DO
         X(I)=PIVOTS(I)*(Y(I)-SUM)
      END DO
*
*     SOLVE (I + inv(D) U_A) x = z FOR x.
*
      DO I=N,1,-1
         SUM=0.0D0
         DO J_PTR=IDIAGP(I)+1,IRP(I+1)-1
            SUM=SUM+VAL(J_PTR)*X(ICOL(J_PTR))
         END DO
         X(I)=X(I)-PIVOTS(I)*SUM
      END DO
*
      RETURN
      END
