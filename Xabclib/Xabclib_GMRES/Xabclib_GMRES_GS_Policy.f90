      SUBROUTINE Xabclib_GMRES_GSPOL (N,NNZ,IRP,ICOL,VAL,IDIAGP,
     $              KIND_PRECOND,MOPE,PRECOND,IATPARAM,RATPARAM,
     $              H,M,M0,V,W_J,Z,UINF,LUINF,INFO)
      IMPLICIT DOUBLE PRECISION ( A-H , O-Z )
*     ..
*     .. Scalar Auruments ..
      INTEGER           N, NNZ, KIND_PRECOND, M, M0, INFO
      INTEGER           LUINF, NUM_SMP
*     ..
*     .. Array Arguments
      INTEGER           IRP(N+1), ICOL(NNZ), IDIAGP(N)
      INTEGER           ICASE
      DOUBLE PRECISION  VAL(NNZ), PRECOND(*), RPCPARM(10), H(M+1,M)
      DOUBLE PRECISION  V(N,M), W_J(N), Z(N)
      DOUBLE PRECISION  UINF(LUINF)
*
      INTEGER           IATPARAM(50)
      DOUBLE PRECISION  RATPARAM(50)
*     ..
*
*  =====================================================================
*  Purpose
*  =======
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
*  IPCPARM (input) INTEGER array, dimension (10)
*          INTEGER parameter set for the preconditioner.
*          IPCPARM(1) = 1 : The preconditioner not generated yet.
*          IPCPARM(1) = 2 : The preconditioner already generated.
*          IPCPARM(2:10)  : Not used.
*
*  RPCPARM (input) DOUBLE PRECISION array, dimension (10)
*          DOUBLE PRECISION parameter set for the preconditioner.
*          if KIND_PRECOND = 3,
*          RPCPARM(1) : The omega parameter for SSOR preconditioning.
*          if KIND_PRECOND = 4,
*          RPCPARM(1) : The breakdown criterion of incomplete LU
*                        decomposition process for ILU(0) precondition-
*                        ing.
*
*  H       (output) DOUBLE PRECISION array, dimension (M+1,M)
*          The Hessenberg matrix H_m.
*
*  M       (input) INTEGER
*          The size of Krylov subspace.
*
*  M0      (output) INTEGER
*          The step count j of the Arnoldi loop at which the GMRES
*          algorithm breaks down.
*          if the GMRES algorithm does not break down, M0 = M.
*
*  V       (output) DOUBLE PRECISION array, dimension (N,M)
*          The Arnoldi vectors (v_1,v_2, ..., v_m) which form an ortho-
*          normal basis of the Krylov subspace
*          K_m = span{v_1,A v_1, ... ,A^(m-1) v_1}.
*
*  W_J     (workspace) DOUBLE PRECISION array, dimension (N)
*
*  Z       (workspace) DOUBLE PRECISION array, dimension (N)
*
*  ICASE   (input) INTEGER
*          Selected Algorithm NO. of matrix-vector product.
*
*  UINF    (input) DOUBLE PRECISION array, dimension (LUINF)
*          Setup Information for OpenATI_DURMV
*
*  LUINF   (input) INTEGER
*          The dimension of UINF.
*          if ICASE = 11 :
*             LUINF >= 0
*          else if ICASE = 12 :
*             LUINF >= INT(0.5*NUM_SMP)+1
*          else if ICASE = 13 :
*             LUINF >= INT(1.5*N)+INT(4.25*JL)+2
*          else if ICASE = 21 :
*             LUINF >= INT(1.125*NNZ)+INT(2.125*JL)+1
*
*  NUM_SMP (input) INTEGER
*          number of threads
*
*  INFO    (output) INTEGER
*          Return code.
*
*  =====================================================================
*     ..
*     .. Local variables
*
      INTEGER          I,J,K
      DOUBLE PRECISION RH
*
      INFO=0
*
      NORMALFLG=0
*
      DO J=1,M
*
         IF (MOPE.EQ.1) THEN
*  COMPUTE w_j = inv(M) A v_j
*          
            CALL OpenATI_DURMV(N,NNZ,IRP,ICOL,VAL,V(1,J),W_J,
     $              IATPARAM,RATPARAM,UINF,LUINF,INFO_A)
            IF (KIND_PRECOND.GE.5) then
               CALL ILUTSOL('N',N,NNZ,IRP,ICOL,VAL,W_J,W_J,IDIAGP,
     $                 KIND_PRECOND,PRECOND,IATPARAM,RATPARAM,INFO)
            ELSE
               CALL Xabclib_PCSLV(N,NNZ,IRP,ICOL,VAL,W_J,W_J,IDIAGP,
     $              PRECOND,IATPARAM,RATPARAM,INFO)
            END IF
         ELSE
*-------------right preconditioner
            IF (KIND_PRECOND.GE.5) then
               CALL ILUTSOL('N',N,NNZ,IRP,ICOL,VAL,Z,V(1,J),IDIAGP,
     $                 KIND_PRECOND,PRECOND,IATPARAM,RATPARAM,INFO)
            ELSE
               CALL Xabclib_PCSLV(N,NNZ,IRP,ICOL,VAL,Z,V(1,J),IDIAGP,
     $              PRECOND,IATPARAM,RATPARAM,INFO)
            END IF
            CALL OpenATI_DURMV(N,NNZ,IRP,ICOL,VAL,Z,W_J,
     $              IATPARAM,RATPARAM,UINF,LUINF,INFO_A)
         ENDIF
*
*====================================================
         CALL OpenATI_DAFGS
     &          (NORMALFLG,N,W_J,V,N,J,H(1,J),IATPARAM,RATPARAM,INFO_A)
*====================================================
*
*  h_ij = ||w_j||^2
*
         RH=0.0D0
!$omp parallel do reduction(+:RH)
         DO K=1,N
            RH=RH+W_J(K)*W_J(K)
         END DO
!$omp end parallel do
         H(J+1,J)=DSQRT(RH)
*
         IF (H(J+1,J) .EQ. 0.0D0) THEN
            M0=J
            GO TO 100
         END IF
*
*  v_j+1 = w_j/h_j+1,j
*
         IF (J .LT. M) THEN
            RH=1.0D0/H(J+1,J)
!$omp parallel do
            DO K=1,N
               V(K,J+1)=W_J(K)*RH
            END DO
!$omp end parallel do
         END IF
      END DO
*
      M0=M
  100 CONTINUE
*
      RETURN
      END
