      SUBROUTINE XNORMD(N,X,S)
      IMPLICIT DOUBLE PRECISION ( A-H , O-Z )
*     ..
*     .. Scalar Auruments ..
      INTEGER          N
      DOUBLE PRECISION S
*     ..
*     .. Array Arguments
      DOUBLE PRECISION X(N)
*     ..
*
*  =====================================================================
*  Purpose
*  =======
*
*  XNORMD computes vector-vector products.
*
*  S = S + X * X
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          The order of the matrix.
*
*  X       (input) DOUBLE PRECISION array, dimension (N)
*          Input vector.
*
*  S       (output) DOUBLE PRECISION
*          Result.
*
*  =====================================================================
*     ..
*     .. Local variables
*
      INTEGER I
*
      S=0.0D0
      DO I=1,N
         S = S + X(I) * X(I)
      ENDDO
*
      RETURN
      END
*
*
      SUBROUTINE XDAXPY(N,ALF,X,Y)
      IMPLICIT DOUBLE PRECISION ( A-H , O-Z )
*     ..
*     .. Scalar Auruments ..
      INTEGER          N
      DOUBLE PRECISION ALF
*     ..
*     .. Array Arguments
      DOUBLE PRECISION X(N), Y(N)
*     ..
*
*  =====================================================================
*  Purpose
*  =======
*
*  XDAXPY computes vector-vector products.
*
*  Y = ALFA * X + Y
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          The order of the matrix.
*
*  ALF     (input) DOUBLE PRECISION
*          Coefficient value.
*
*  X       (input) DOUBLE PRECISION array, dimension (N)
*          Input vector.
*
*  Y       (input/output) DOUBLE PRECISION array, dimension (N)
*          On entry,
*           Input vector.
*          On exit,
*           Output vector.
*
*  =====================================================================
*     ..
*     .. Local variables
*
      INTEGER I
*
      DO I=1,N
         Y(I) = ALF * X(I) + Y(I)
      ENDDO
*
      RETURN
      END
*
*
      SUBROUTINE XDXPAY(N,ALF,X,Y)
      IMPLICIT DOUBLE PRECISION ( A-H , O-Z )
*     ..
*     .. Scalar Auruments ..
      INTEGER          N
      DOUBLE PRECISION ALF
*     ..
*     .. Array Arguments
      DOUBLE PRECISION X(N), Y(N)
*     ..
*
*  =====================================================================
*  Purpose
*  =======
*
*  XDXPAY computes vector-vector products.
*
*  Y = X + ALFA * Y
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          The order of the matrix.
*
*  ALF     (input) DOUBLE PRECISION
*          Coefficient value.
*
*  X       (input) DOUBLE PRECISION array, dimension (N)
*          Input vector.
*
*  Y       (input/output) DOUBLE PRECISION array, dimension (N)
*          On entry,
*           Input vector.
*          On exit,
*           Output vector.
*
*  =====================================================================
*     ..
*     .. Local variables
*
      INTEGER I
*
      DO I=1,N
         Y(I) = X(I) + ALF * Y(I)
      ENDDO
*
      RETURN
      END
*
*
      SUBROUTINE XDDOT(N,X,Y,S)
      IMPLICIT DOUBLE PRECISION ( A-H , O-Z )
*     ..
*     .. Scalar Auruments ..
      INTEGER          N
      DOUBLE PRECISION S
*     ..
*     .. Array Arguments
      DOUBLE PRECISION X(N), Y(N)
*     ..
*
*  =====================================================================
*  Purpose
*  =======
*
*  XDDOT computes vector-vector products.
*
*  Y = X + ALFA * Y
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          The order of the matrix.
*
*  X       (input) DOUBLE PRECISION array, dimension (N)
*          Input vector.
*
*  Y       (input) DOUBLE PRECISION array, dimension (N)
*          Input vector.
*
*  S       (output) DOUBLE PRECISION
*          Result.
*
*  =====================================================================
*     ..
*     .. Local variables
*
      INTEGER I
*
      S=0.0D0
      DO I=1,N
         S = S + X(I) * Y(I)
      ENDDO
*
      RETURN
      END
*
*
      SUBROUTINE XNORM2(N,X,S)
      IMPLICIT DOUBLE PRECISION ( A-H , O-Z )
*     ..
*     .. Scalar Auruments ..
      INTEGER          N
      DOUBLE PRECISION S
*     ..
*     .. Array Arguments
      DOUBLE PRECISION X(N)
*     ..
*
*  =====================================================================
*  Purpose
*  =======
*
*  XNORM2 computes vector-vector products.
*
*  S = S + X * X
*  S = SQRT(S)
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          The order of the matrix.
*
*  X       (input) DOUBLE PRECISION array, dimension (N)
*          Input vector.
*
*  S       (output) DOUBLE PRECISION
*          Result.
*
*  =====================================================================
*     ..
*     .. Local variables
*
      INTEGER I
*
      S=0.0D0
      DO I=1,N
         S = S + X(I) * X(I)
      ENDDO
      S = DSQRT(S)
*
      RETURN
      END
