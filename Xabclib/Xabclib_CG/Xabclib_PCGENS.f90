      SUBROUTINE Xabclib_PCGENS(N,NNZ,IRP,ICOL,VAL,
     $                          PRECOND,IATPARAM,RATPARAM,INFO)
      IMPLICIT DOUBLE PRECISION ( A-H , O-Z )
*     ..
*     .. Scalar Auruments ..
      INTEGER           N, NNZ, KIND_PRECOND, INFO
*     ..
*     .. Array Arguments
      INTEGER           IRP(N+1), ICOL(NNZ)
      DOUBLE PRECISION  VAL(NNZ), PRECOND(*)
*
      INTEGER           IATPARAM(50)
      DOUBLE PRECISION  RATPARAM(50)
*     ..
*
*  =====================================================================
*  Purpose
*  =======
*
*  Xabclib_PCGENS generates the preconditioner used for the r
*  ICCG algorithm.
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
*  PRECOND (output) DOUBLE PRECISION array
*          Generated preconditioner data.
*
*  IATPARAM (input) INTEGER array, dimension (50)
*          Integer parameter set for OpenATI and Xabclib.
*           IATPARAM(25)    : Preconditioner type.
*
*  RATPARAM (input/output) DOUBLE PRECISION array, dimension (50)
*          Double precision parameter set for OpenATI and Xabclib.
*          On entry,
*          if IATPARAM(25) = 3,
*          RATPARAM(25) : The omega parameter for SSOR preconditioning.
*          if IATPARAM(25) = 4,
*          RATPARAM(25) : The breakdown criterion of incomplete LU
*                         decomposition process for ILU(0) precondition-
*                         ing.
*          On exit,
*          RATPARAM(30) : Floating operations(*10^9 operations).
*
*  INFO    (output) INTEGER
*          Return code.
*
*  =====================================================================
*     ..
*     .. Local variables
*
      DOUBLE PRECISION FLOP
*
      INFO=0
*
      KIND_PRECOND = IATPARAM(25)
*
      IF (KIND_PRECOND.EQ.1) THEN
*
         FLOP=0.0D0
*-- PRECONDITIONER NOT USED
*
      ELSE IF (KIND_PRECOND.EQ.2) THEN
*
*-- JACOBI PRECONDITIONER
*
         CALL Xabclib_PCGENS2(N,NNZ,IRP,ICOL,VAL,PRECOND,INFO)
         FLOP=2.0D0*NZ+10.0D0*N
*
      ELSE IF (KIND_PRECOND.EQ.3) THEN
*
*-- JACOBI ITERATIVE REFINEMENT
*
         CALL Xabclib_PCGENS3(N,NNZ,IRP,ICOL,VAL,PRECOND,INFO)
         FLOP=1.0D0*N
*
      ELSE IF (KIND_PRECOND .EQ. 4) THEN
*
*-- ICOD PRECONDITIONER
*
         CALL Xabclib_IC0D(N,NNZ,IRP,ICOL,VAL,PRECOND,FLOP,INFO)
         IF (INFO .NE. 0) THEN
            WRITE(6,*)" DIAGONAL ELEMENTS WARNING "
            INFO = 0
         END IF
         FLOP=2.0D0*NZ
*
      ELSE 
*
*-- KIND_PRECOND > 4 NOT SUPPORTED
*
         INFO=-8
      END IF
*
      RATPARAM(30)=RATPARAM(30)+FLOP
      RETURN
      END
*-----------------------------------------------------------------------
      SUBROUTINE Xabclib_PCGENS2 (N,NNZ,IRP,ICOL,VAL,PIVOTS,INFO)
      IMPLICIT DOUBLE PRECISION ( A-H , O-Z )
*     ..
*     .. Scalar Auruments ..
      INTEGER           N, NNZ, INFO
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
*  Xabclib_PCGENS2 generates the Jacobi preconditioner.
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
*  PIVOTS  (output) DOUBLE PRECISION array, dimension (N)
*          Generated Jacobi preconditioner data.
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
      DO I=1,N
        IF (ABS(VAL(IRP(I))) .LT. 2.2D-16) THEN
            INFO=100
            RETURN
         ELSE IF (VAL(IRP(I)) .LT. 0.0D0) THEN
            write(6,*) "DIAG<0.0D0",I,VAL(IRP(I))
            INFO=110
            RETURN
         ELSE
            PIVOTS(I)=1.0D0/DSQRT(VAL(IRP(I)))
            VAL(IRP(I))=1.0D0
         END IF
      END DO
*
*    ------------ Jacobi scaling : A(I,J) := D(I) * A(I,J) * D(J)
*    ------------ Matrix 'A' is Update.
*
!$omp parallel do schedule(static,4) private(I,DI,JC,J,DJ)
      DO I=1,N
        DI=PIVOTS(I)
        DO JC=IRP(I)+1,IRP(I+1)-1
          J=ICOL(JC) 
          DJ=PIVOTS(J)
          VAL(JC)=VAL(JC)*DI*DJ 
        ENDDO
      ENDDO
!$omp end parallel do
*
      RETURN
      END
*-----------------------------------------------------------------------
      SUBROUTINE Xabclib_PCGENS3 (N,NNZ,IRP,ICOL,VAL,PIVOTS,INFO)
      IMPLICIT DOUBLE PRECISION ( A-H , O-Z )
*     ..
*     .. Scalar Auruments ..
      INTEGER           N, NNZ, INFO
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
*  Xabclib_PCGENS3 generates the Jacobi Iterative preconditioner.
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
*  PIVOTS  (output) DOUBLE PRECISION array, dimension (N)
*          Generated Jacobi preconditioner data.
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
      DO I=1,N
        IF (ABS(VAL(IRP(I))) .LT. 2.2D-16) THEN
            INFO=100
            RETURN
        END IF
      END DO
*
*    ------------ diagonal dominant check
      DO I=1,N
         PIVOTS(I)=0.0D0
      END DO
*
      DO I=1,N
        DI=DABS(VAL(IRP(I)))
        SUM=0.0D0
        DO JC=IRP(I)+1,IRP(I+1)-1
          SUM=SUM+DABS(VAL(JC))
          PIVOTS(ICOL(JC))=PIVOTS(ICOL(JC))+DABS(VAL(JC))
        ENDDO
        SUM=SUM+PIVOTS(I)
        IF (DI*1.00001D0 .LT. SUM) THEN
          INFO=120
          write(6,*) 'diag. domi. check I=',i,DI,SUM
          RETURN
        ENDIF
      ENDDO
*
*    ------------ reciprocal diagonal
      DO I=1,N
         PIVOTS(I)=1.0D0/VAL(IRP(I))
      END DO
*
      RETURN
      END
