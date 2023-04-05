      SUBROUTINE Xabclib_IDIAGPGEN (N,NNZ,IRP,ICOL,IDIAGP,INFO)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
*     ..
*     .. Scalar Auruments ..
      INTEGER           N, NNZ, INFO
*     ..
*     .. Array Arguments
      INTEGER           IRP(N+1), ICOL(NNZ), IDIAGP(N)
*     ..
*
*  =====================================================================
*  Purpose
*  =======
*
*  Xabclib_IDIAGPGEN generates the index list of the diagonal elements
*  in VAL.
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
*  IDIAGP  (output) INTEGER array, dimension (N)
*          The index of the diagonal elements in VAL.
*
*  INFO    (output) INTEGER
*          Return code.
*
*  =====================================================================
*     ..
*     .. Local variables
*
      INTEGER I,J,IIND
*
      INFO=0
*
      DO J=1,N
         IDIAGP(J)=0
         DO IIND=IRP(J),IRP(J+1)-1
            IF (J .EQ. ICOL(IIND)) THEN
               IDIAGP(J)=IIND
               GO TO 100
            END IF
         END DO
         INFO=100
         RETURN
  100    CONTINUE
      END DO
*
      RETURN
      END
