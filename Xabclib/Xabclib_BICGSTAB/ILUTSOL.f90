      SUBROUTINE ILUTSOL (POS,N,NNZ,IRP,ICOL,VAL,X,Y,IDIAGP,
     $                    KIND_PRECOND,PRECOND,IATPARAM,RATPARAM,INFO)
      IMPLICIT DOUBLE PRECISION ( A-H , O-Z )
*     ..
*     .. Scalar Auruments ..
      CHARACTER*1 POS
      INTEGER           N, NNZ, KIND_PRECOND, INFO
*     ..
*     .. Array Arguments
      INTEGER           IRP(N+1), ICOL(NNZ), IDIAGP(N), IATPARAM(50)
      DOUBLE PRECISION  VAL(NNZ), X(N), Y(N), PRECOND(*)
      DOUBLE PRECISION  RATPARAM(50)
*     ..
*
*  =====================================================================
*  Purpose
*  =======
*
*  ILUTSOL 
*     IF POS='N'THEN, solves        Mx = y for x. (x = inv(M)*y)
*     IF POS='T'THEN, solves trans(M)x = y for x. (x = inv(M)*y)
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
*  IATPARAM (input) INTEGER array, dimension (50)
*
*  RATPARAM (input) DOUBLE PRECISION array, dimension (50)
*          if KIND_PRECOND = 6,
*
*  INFO    (output) INTEGER
*          Return code.
*
*  =====================================================================
*
      INFO=0
      IF (KIND_PRECOND .EQ. 1) THEN
*
*-- PRECONDITIONER NOT APPLIED
*
         DO I=1,N
            X(I)=Y(I)
         END DO
*
      ELSE IF (KIND_PRECOND .EQ. 5 ) THEN
         NZMAX=NNZ
         IP_ALU=1
         IP_JLU=IP_ALU + NZMAX +10
         IP_JU =IP_JLU + NZMAX/2+1 +10
*
         IF (POS.EQ.'N') THEN
            CALL LUSOL (N,Y,X,PRECOND(IP_ALU),PRECOND(IP_JLU),
     $                  PRECOND(IP_JU))
         ELSE IF (POS.EQ.'T') THEN
            CALL LUTSOL(N,Y,X,PRECOND(IP_ALU),PRECOND(IP_JLU),
     $                  PRECOND(IP_JU))
         END IF
*
      ELSE IF (KIND_PRECOND .EQ. 6 ) THEN
         LFIL=IATPARAM(26)
         NZMAX=N*(LFIL*2+1)
         IP_ALU=1
         IP_JLU=IP_ALU + NZMAX  +10
         IP_JU =IP_JLU + (NZMAX-1)/2+1 +10
*
         IF (POS.EQ.'N') THEN
            CALL LUSOL (N,Y,X,PRECOND(IP_ALU),PRECOND(IP_JLU),
     $                  PRECOND(IP_JU))
         ELSE IF (POS.EQ.'T') THEN
            CALL LUTSOL(N,Y,X,PRECOND(IP_ALU),PRECOND(IP_JLU),
     $                  PRECOND(IP_JU))
         END IF
*
      ELSE IF (KIND_PRECOND .EQ. 7 ) THEN
         LFIL=IATPARAM(26)
         NZMAX=20*NNZ
         IP_ALU=1
         IP_JLU=IP_ALU + NZMAX         +10
         IP_JU =IP_JLU + (NZMAX)/2+1 +10
*
         IF (POS.EQ.'N') THEN
            CALL LUSOL (N,Y,X,PRECOND(IP_ALU),PRECOND(IP_JLU),
     $                  PRECOND(IP_JU))
         ELSE IF (POS.EQ.'T') THEN
            CALL LUTSOL(N,Y,X,PRECOND(IP_ALU),PRECOND(IP_JLU),
     $                  PRECOND(IP_JU))
         END IF
      ELSE 
*-- KIND_PRECOND > 6 NOT SUPPORTED
         INFO=-8
      END IF
*
      RETURN
      END
