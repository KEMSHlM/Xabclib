      SUBROUTINE Xabclib_BICGSTAB_PCHK
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
*  Xabclib_BICGSTAB_PCHK check parameters of BICGSTAB.
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
*           IATPARAM(10)    : Fastest OpenATI_DURMV impl. method.
*           IATPARAM(24)    : Preconditioner operations flag.
*           IATPARAM(25)    : Preconditioner type.
*
*  RATPARAM (input) DOUBLE PRECISION array, dimension (50)
*          Double precision parameter set for OpenATI and Xabclib.
*
*  NPRE    (input) INTEGER
*          The dimension of PRECOND.
*          if KIND_PRECOND = 1,        NPRE >= 0.
*          if KIND_PRECOND = 2,3 or 4, NPRE >= N.          
*          if KIND_PRECOND = 5,        NPRE >= 3*NNZ/2+2*N+50.
*          if KIND_PRECOND = 6,        
*                               NPRE >= 3*(2.0*IFILL+1)*N/2+3*N+50.
*                               (IFILL=IATPARAM(26))
*
*  LWORK   (input) INTEGER
*          The dimension of WORK.
*          LWORK >= (9)*N + (N-1)/2+1
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
      INTEGER          I,J,I_PTR,IDPCNT,IOPTPC,ITER_MAX,NNZLMT
      INTEGER          MVCASE_USYM,KIND_PRECOND
      PARAMETER        (NNZLMT=1908874111)
      DOUBLE PRECISION STOP_TOL,MAX_ETIME
*
      INFO=0
*
      MVCASE_USYM  = IATPARAM(10)
      KIND_PRECOND = IATPARAM(25)
      IFILL        = IATPARAM(26)
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
*-- NNZ LIMIT CHECK ( for original Seg. Scan MV)
      IF (MVCASE_USYM.EQ.21 .AND. NNZ.GT.NNZLMT) THEN
         INFO=700
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
      ELSE IF (KIND_PRECOND.GE.2 .AND. KIND_PRECOND.LE.4) THEN
         IF (NPRE .LT. N) THEN
            INFO=-9
            GO TO 9999
         END IF
      ELSE IF (KIND_PRECOND.EQ.5) THEN
         IF (NPRE .LT. 3*NNZ/2+2*N+50) THEN
            INFO=-9
            GO TO 9999
         END IF
      ELSE IF (KIND_PRECOND.EQ.6) THEN
         IF (NPRE .LT. 3*(2.0D0*IFILL+1)*N/2+3*N+50) THEN
            INFO=-9
            GO TO 9999
         END IF
      ELSE IF (KIND_PRECOND.EQ.7) THEN
C        IF (NPRE .LT. 3*(2.0D0*IFILL+1)*N/2+3*N+50) THEN
C           INFO=-9
C           GO TO 9999
C        END IF
      ELSE
*
      END IF
*
*-- Convergence criterion CHECK
      IF (RATPARAM(23) .LT. 1.0D-15 ) THEN
       WRITE(6,*) "Warning by Xabclib: Convergence criterion
     $ (RATPARAM(23)) is too small! RATPARAM(23) = ", RATPARAM(23)
      END IF
*
*-- LWORK >= 8*N+(N-1)/2+1
      IF (LWORK .LT. 9*N+(N-1)/2+1) THEN
         INFO=-13
         GO TO 9999
      END IF 
*---------
 9999 CONTINUE
      RETURN
      END
