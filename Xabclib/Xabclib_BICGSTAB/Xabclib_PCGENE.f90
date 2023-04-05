      SUBROUTINE Xabclib_PCGENE(N,NNZ,IRP,ICOL,VAL,IDIAGP,
     $                          PRECOND,IATPARAM,RATPARAM,INFO)
      IMPLICIT DOUBLE PRECISION ( A-H , O-Z )
*     ..
*     .. Scalar Auruments ..
      INTEGER           N, NNZ, KIND_PRECOND, INFO
*     ..
*     .. Array Arguments
      INTEGER           IRP(N+1), ICOL(NNZ), IDIAGP(N)
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
*  Xabclib_PCGEN generates the preconditioner used for the restarted
*  GMRES algorithm.
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
*  IDIAGP  (input) INTEGER array, dimension (N)
*          The index of the diagonal elements in VAL.
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
      INTEGER          NZMAX,IP_ALU,IP_JLU,IP_JU,IP_IW
      DOUBLE PRECISION OMEGA,FLOP,BD_TOL
*
      INFO=0
*
      KIND_PRECOND = IATPARAM(25)
*
      IF (KIND_PRECOND.EQ.1) THEN
*
*-- PRECONDITIONER NOT USED
*
      ELSE IF (KIND_PRECOND.EQ.2 .OR. KIND_PRECOND.EQ.3) THEN
*
*-- JACOBI/SOR PRECONDITIONER
*
         IF (KIND_PRECOND .EQ. 3) THEN
            OMEGA = RATPARAM(25)
            IF (OMEGA.LE.0.0D0 .OR. OMEGA.GE.2.0D0) THEN
               INFO=-12
               RETURN
            END IF
            IF (OMEGA.LT.1.0D0) THEN
               RATPARAM(25) = 1.0D0
               OMEGA = RATPARAM(25)
            END IF
         END IF
         CALL Xabclib_PCGEN1(N,NNZ,IRP,ICOL,VAL,IDIAGP,PRECOND,INFO)
         FLOP=1.0D0*N
      ELSE IF (KIND_PRECOND .EQ. 4) THEN
*
*-- ILU(0) PRECONDITIONER
*
         BD_TOL = RATPARAM(25)
         CALL Xabclib_PCGEN2(N,NNZ,IRP,ICOL,VAL,IDIAGP,PRECOND,BD_TOL,
     $                       FLOP,INFO)
      ELSE IF (KIND_PRECOND .EQ. 5) THEN
*-- KIND_PRECOND = 5 ILU0 with LU-update
*
         NZMAX=NNZ
         IP_ALU=1
         IP_JLU=IP_ALU + NZMAX+10
         IP_JU =IP_JLU + NZMAX/2+1 +10
         IP_IW =IP_JU  + N/2+1
         CALL ILU0(N,VAL,ICOL,IRP,PRECOND(IP_ALU),
     $             PRECOND(IP_JLU),PRECOND(IP_JU),PRECOND(IP_IW),
     $             INFO)
         IF (INFO.NE.0) THEN
            INFO=100
         ENDIF
      ELSE IF (KIND_PRECOND .EQ. 6) THEN
*
*-- KIND_PRECOND = 6 ILUT ()
         LFIL=IATPARAM(26)
         TOL =RATPARAM(25)
*
         NZMAX=N*(LFIL*2+1)
         IP_ALU=1
         IP_JLU=IP_ALU + NZMAX +10
         IP_JU =IP_JLU + (NZMAX-1)/2+1 +10
         IP_W  =IP_JU  + N/2+1
         IP_JW =IP_W   + N+1
         CALL ILUT(N,VAL,ICOL,IRP,LFIL,TOL,PRECOND(IP_ALU),
     $             PRECOND(IP_JLU),PRECOND(IP_JU),NZMAX,PRECOND(IP_W),
     $             PRECOND(IP_JW),INFO)
         IF (INFO.NE.0) THEN
            INFO=100
         ENDIF
*
C===========ILU(K) implemented 2012.08
*     ELSE IF (KIND_PRECOND .EQ. 7) THEN
*        LFIL=IATPARAM(26)
*
*        NZMAX=20*NNZ
*        write(6,*) ' ILUK NZMAX=',NZMAX
*        IP_ALU  =1
*        IP_JLU  =IP_ALU + NZMAX     +10
*        IP_JU   =IP_JLU + NZMAX/2+1 +10
*        IP_LEVS =IP_JU  + N/2+1     +10
*        IP_W    =IP_LEVS+ NZMAX/2+1
*        IP_JW   =IP_W   + N
*        CALL iluk(N,VAL,ICOL,IRP,LFIL,PRECOND(IP_ALU),
*    $            PRECOND(IP_JLU),PRECOND(IP_JU),PRECOND(IP_LEVS),
*    $            NZMAX,
*    $            PRECOND(IP_W),PRECOND(IP_JW),INFO)
*        write(6,*) '*ILU(K) K,INFO=',LFIL,INFO
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
      SUBROUTINE Xabclib_PCGEN1 (N,NNZ,IRP,ICOL,VAL,IDIAGP,PIVOTS,INFO)
      IMPLICIT DOUBLE PRECISION ( A-H , O-Z )
*     ..
*     .. Scalar Auruments ..
      INTEGER           N, NNZ, INFO
*     ..
*     .. Array Arguments
      INTEGER           IRP(N+1), ICOL(NNZ), IDIAGP(N)
      DOUBLE PRECISION  VAL(NNZ), PIVOTS(N)
*     ..
*
*  =====================================================================
*  Purpose
*  =======
*
*  Xabclib_PCGEN1 generates the Jacobi/SSOR preconditioner.
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
*  IDIAGP  (input) INTEGER array, dimension (N)
*          The index of the diagonal elements in VAL.
*
*  PIVOTS  (output) DOUBLE PRECISION array, dimension (N)
*          Generated Jacobi/SSOR preconditioner data.
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
         IF (DABS(VAL(IDIAGP(I))) .LT. 2.2D-16) THEN
            INFO=100
            RETURN
         ELSE
            PIVOTS(I)=1.0D0/VAL(IDIAGP(I))
         END IF
      END DO
*
      RETURN
      END
*-----------------------------------------------------------------------
      SUBROUTINE Xabclib_PCGEN2(N,NNZ,IRP,ICOL,VAL,IDIAGP,PIVOTS,
     &              BD_TOL,FLOP,INFO)
      IMPLICIT DOUBLE PRECISION ( A-H , O-Z )
*     ..
*     .. Scalar Auruments ..
      INTEGER           N, NNZ, INFO
      DOUBLE PRECISION  BD_TOL, FLOP
*     ..
*     .. Array Arguments
      INTEGER           IRP(N+1), ICOL(NNZ), IDIAGP(N), IPCPARM(10)
      DOUBLE PRECISION  VAL(NNZ), PIVOTS(N), RPCPARM(10)
*     ..
*
*  =====================================================================
*  Purpose
*  =======
*
*  Xabclib_PCGEN2 generates the ILU(0) preconditioner.
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
*  IDIAGP  (input) INTEGER array, dimension (N)
*          The index of the diagonal elements in VAL.
*
*  PIVOTS  (output) DOUBLE PRECISION array, dimension (N)
*          Generated ILU(0) preconditioner data.
*
*  BD_TOL  (input) DOUBLE PRECISION
*          The breakdown criterion of incomplete LU decomposition
*          process for ILU(0) preconditioning.
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
*
      INTEGER          I,J,K,J_PTR,K_PTR
      INTEGER          IFOUND
      DOUBLE PRECISION ELEMENT
*
      INFO=0
*
      FLOP=0.0D0
*
      IF (BD_TOL .LE. 0.0D0) THEN
         INFO=-12
         RETURN
      END IF
*
      DO I=1,N
         PIVOTS(I)=VAL(IDIAGP(I))
      ENDDO
*
      DO I=1,N
         PIVOTS(I)=1.0D0/PIVOTS(I)
         DO J=IDIAGP(I)+1,IRP(I+1)-1
            IFOUND=0
            DO K=IRP(ICOL(J)),IDIAGP(ICOL(J))-1
               IF (ICOL(K) .EQ. I) THEN
                  IFOUND=1
                  ELEMENT=VAL(K)
               END IF
            END DO
            IF (IFOUND .EQ. 1) THEN
               PIVOTS(ICOL(J))=PIVOTS(ICOL(J))-ELEMENT*PIVOTS(I)*VAL(J)
               FLOP=FLOP+3.0D0
            END IF
         END DO
      END DO
*
      DO I=1,N
         ELEMENT=1.0D0/(BD_TOL*DABS(VAL(IDIAGP(I))))
         IF (DABS(PIVOTS(I)) .GT. ELEMENT) THEN
            INFO=100
            RETURN
         END IF
      END DO
      FLOP=FLOP+2.0D0*N
*
      RETURN
      END
