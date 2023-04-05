      SUBROUTINE Xabclib_CG_PCHK
     $           (N,NNZ,IRP,ICOL,VAL,IATPARAM,RATPARAM,NPRE,LWORK,INFO)
      IMPLICIT DOUBLE PRECISION ( A-H , O-Z )
*     ..
*     .. Scalar Auruments ..
      INTEGER           N, NNZ, NPRE, LWORK, INFO
*     ..
*     .. Array Arguments
      INTEGER           IRP(N+1), ICOL(NNZ)
      INTEGER           IATPARAM(50)
      DOUBLE PRECISION  VAL(NNZ), RATPARAM(50)
*     ..
*
*  =====================================================================
*  Purpose
*  =======
*
*  Xabclib_CG_PCHK check parameters of ICCG.
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
*  IATPARAM (input) Integer array, dimension (50)
*          Integer parameter set for OpenATI and Xabclib.
*           IATPARAM(24)    : Preconditioner operations flag.
*           IATPARAM(25)    : Preconditioner type.
*
*  RATPARAM (input) DOUBLE PRECISION array, dimension (50)
*          Double precision parameter set for OpenATI and Xabclib.
*
*  NPRE    (input) INTEGER
*          The dimension of PRECOND.
*          if KIND_PRECOND = 1,      NPRE >= 0.
*          if KIND_PRECOND = 2,      NPRE >= N.
*          if KIND_PRECOND = 4,      NPRE >= N+(N+1)/2+1+NZ/2+1+NZ.
*
*  LWORK   (input) INTEGER
*          The dimension of WORK.
*          LWORK >= (4)*N
*
*  INFO    (output) INTEGER
*          Return code.
*
*  Local valiables
*  ===============
*
*  IOPTPC     INTEGER
*             Control flag for generating preconditioner.
*  MAX_ITER   INTEGER
*             Maximum number of restart iterations.
*  STOP_TOL   DOUBLE PRECISION
*             Convergence criterion.
*  MAX_ETIME  DOUBLE PRECISION
*             Maximum elapsed time for the solver.
*
*  =====================================================================
*     ..
*     .. Local variables
*
      INTEGER          I,J,I_PTR,IDPCNT,IOPTPC,ITER_MAX
      INTEGER          KIND_PRECOND
      DOUBLE PRECISION STOP_TOL,MAX_ETIME
*
      INFO=0
*
      KIND_PRECOND = IATPARAM(25)
      IOPTPC       = IATPARAM(24)
*
*-- N >= 2
      IF (N .LT. 2) THEN
         INFO=-1
         GO TO 9999
      END IF
*-- NNZ >= 2, NNZ=IRP(N+1)-1
      IF (N.LT.2 .OR. NNZ.NE.IRP(N+1)-1) THEN
         INFO=-2
         GO TO 9999
      END IF
*-- IRP
*
*-- NPRE VALIDITY CHECK
      IF (KIND_PRECOND .EQ. 1) THEN
         IF (NPRE .LT. 0) THEN
            INFO=-9
            GO TO 9999
         END IF
      ELSE IF (KIND_PRECOND.GE.2) THEN
         IF (NPRE .LT. N) THEN
            INFO=-9
            GO TO 9999
         END IF
      ELSE IF (KIND_PRECOND.LE.4) THEN
         IF (NPRE .LT. N+(N+1)/2+1+NZ/2+1+NZ) THEN
            INFO=-9
            GO TO 9999
         END IF
      ELSE
*
      END IF
*
*-- LWORK >= 4*N
      IF (LWORK .LT. 4*N) THEN
         INFO=-13
         GO TO 9999
      END IF
*
*-- DIAG CHEACK
      DO I = 1, N
         DO J = IRP(I), IRP(I+1)-1
            IF (I .EQ. ICOL(J)) THEN
               GO TO 100
            END IF
         END DO
         INFO = 100
         GO TO 9999
  100    CONTINUE
      END DO
*-- RATPARAM
      STOP_TOL = RATPARAM(23)
      IF (STOP_TOL .LE. 0.0D0) THEN
         RATPARAM(23) = 1.0D-8
      END IF
*---------
 9999 CONTINUE
      RETURN
      END
