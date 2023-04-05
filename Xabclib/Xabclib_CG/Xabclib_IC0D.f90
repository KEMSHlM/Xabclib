      SUBROUTINE Xabclib_IC0D(N,NNZ,IRP,ICOL,VAL,PIVOTS,FLOP,INFO)
      IMPLICIT DOUBLE PRECISION ( A-H , O-Z )
*     ..
*     .. Scalar Auruments ..
      INTEGER           N, NNZ, INFO
      DOUBLE PRECISION  FLOP
*     ..
*     .. Array Arguments
      INTEGER           IRP(N+1), ICOL(NNZ)
      DOUBLE PRECISION  VAL(NNZ), PIVOTS(N)
*     ..
*
*  =====================================================================
*  Purpose
*  =======
*
*   Xabclib_IC0D generates the Incomplete Cholesky diagonal 
*   preconditioner.
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
*          upper triangular matrix.
*
*  ICOL    (input) INTEGER array, dimension (NNZ)
*          The column indices of each nonzero.
*
*  VAL     (input) DOUBLE PRECISION array, dimension (NNZ)
*          Nonzero elements of the matrix in Compressed Row Storage
*          format.
*
*  PIVOTS  (output) DOUBLE PRECISION array, dimension (N)
*          Generated ICOD preconditioner data.
*
*  FLOP    (output) DOUBLE PRECISION
*          Floating operations(*10^9 operations).
*
*  INFO    (output) INTEGER
*          Return code.
*
*  =====================================================================
*     ..
*     .. Local variables
      INTEGER          I,IC,J
      DOUBLE PRECISION PIV
*
      INFO = 0
*
      FLOP = 0.0D0
*
      DO I = 1, N
         PIVOTS(I) = VAL(IRP(I))
      END DO
      DO I = 1, N
         PIV = PIVOTS(I)
         IF (DABS(PIV) .LT. 1.0D-13) THEN
            IF (PIV .GT. 0.0D0) THEN
               PIV =  1.0D-13
            ELSE
               PIV = -1.0D-13
            END IF
         END IF
         PIVOTS(I) = 1.0D0 / PIV
         DO J = IRP(I)+1, IRP(I+1)-1
            IC = ICOL(J)
            PIVOTS(IC) = PIVOTS(IC) - VAL(J) * VAL(J) * PIVOTS(I)
            FLOP = FLOP + 3.0D0
         END DO
      END DO
      FLOP = FLOP + 2.0D0 * N
*
      RETURN
      END
