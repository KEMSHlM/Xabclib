      SUBROUTINE Xabclib_EVAL_SYM (N,NNZ,IRP,ICOL,VAL,X,B,RI,RERR,INFO)
      IMPLICIT DOUBLE PRECISION ( A-H , O-Z )
*     ..
*     .. Scalar Auruments ..
      INTEGER           N, NNZ, INFO
      DOUBLE PRECISION  RERR
*     ..
*     .. Array Arguments
      INTEGER           IRP(N+1), ICOL(NNZ)
      DOUBLE PRECISION  VAL(NNZ), X(N), B(N), RI(N)
*     ..
*
*  =====================================================================
*  Purpose
*  =======
*
*  Xabclib_EVAL_SYM estimates the error of the approximate solution 
*  vector. The preconditioner is not applied to the estimation.
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
*  RI      (workspace) DOUBLE PRECISION array, dimension (N)
*
*  RERR    (output) DOUBLE PRECISION
*          The error of the approximate solution vector.
*          ( RERR = || b - A x_i ||_2 / ||  b ||_2 )
*
*  INFO    (output) INTEGER
*          Return code.
*
*  Local valirables
*  ================
*
*  BETA0   DOUBLE PRECISION 
*          The square norm of the right hand sede vector.
*          ( BETA0 = || b ||_2 )
*  BETA    DOUBLE PRECISION
*          The square norm of the residual vector.
*          ( BETA  = || b - A x_i ||_2 )
*
*  =====================================================================
*     ..
*     .. Local variables
*
      INTEGER          I
      DOUBLE PRECISION BETA0,BETA
*
      INFO=0
*
*  beta_0 = || b ||_2
*
      BETA0=0.0D0
      DO I=1,N
         BETA0=BETA0+B(I)*B(I)
      END DO
      BETA0=DSQRT(BETA0)
*
*  beta = || b - A x_i ||_2
*
      CALL OpenATI_DSRMV_11(N,NNZ,IRP,ICOL,VAL,X,RI)
      DO I=1,N
         RI(I)=RI(I)-B(I)
      END DO
      BETA=0.0D0
      DO I=1,N
         BETA=BETA+RI(I)*RI(I)
      END DO
      BETA=DSQRT(BETA)
*
*  error = beta / beta_0
*
      RERR=BETA/BETA0
*
      RETURN
      END
