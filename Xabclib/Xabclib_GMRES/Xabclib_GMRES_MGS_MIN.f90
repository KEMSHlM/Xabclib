      SUBROUTINE Xabclib_GMRES_MGS_MIN (N,NNZ,IRP,ICOL,VAL,X,B,IDIAGP,
     $              KIND_PRECOND,MOPE,PRECOND,IATPARAM,RATPARAM,
     $              H,M,M0,V,W_J,Z,G,R,TMP,BETA,RERR,STOP_TOL,UINF,
     $              LUINF,INFO)
      IMPLICIT DOUBLE PRECISION ( A-H , O-Z )
*     ..
*     .. Scalar Auruments ..
      INTEGER           N, NNZ, KIND_PRECOND, M, M0, INFO
      INTEGER           LUINF, NUM_SMP
*     ..
*     .. Array Arguments
      INTEGER           IRP(N+1), ICOL(NNZ), IDIAGP(N)
      INTEGER           ICASE, MFLG
      DOUBLE PRECISION  BETA, BETA0, RERR, STOP_TOL
      DOUBLE PRECISION  VAL(NNZ), PRECOND(*), RPCPARM(10), H(M+1,M)
      DOUBLE PRECISION  R(M+1,M),TMP(2,M)
      DOUBLE PRECISION  V(N,M), W_J(N), Z(N), X(N), G(M+1), B(N)
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
*  Xabclib_GMRES_MGS_MIN generates the Arnoldi basis and the Hessenberg
*  matirx H_m using the Arnoldi algorithm with the modified Gram-Schmidt
*  algorithm.
*  Xabclib_GMRES_MGS_MIN computes the minimizer y_m of
*  ||beta e_1 - H_m y||_2 and updates the approximate solution
*  vector x (x_m = x_0 + V_m y_m).
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
*          format.a
*
*  X       (input/output) DOUBLE PRECISION array, dimension (N)
*          On entry, The approximate solution vector not updated.
*          On exit,  The approximate solution vector updated.
*
*  B       (input) DOUBLE PRECISION array, dimension (N)
*          The right hand side vector.
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
*  G       (workspace) DOUBLE PRECISION array, dimension (M+1)
*
*  R       (workspace) DOUBLE PRECISION array, dimension (M+1,M)
*
*  TMP     (workspace) DOUBLE PRECISION array, dimension (2,M)
*
*  BETA    (output) DOUBLE PRECISION
*          The square norm of the residual vector.
*          ( BETA = || inv(M) R0 ||_2 )
*
*  RERR    (output) DOUBLE PRECISION
*          The error of the approximate solution vector.
*          ( RERR = ||G(J+1)|| / || inv(M) b ||_2 )
*
*  STOP_TOL   DOUBLE PRECISION
*             Convergence criterion.
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
*  INFO    (output) INTEGER
*          Return code.
*
*  =====================================================================
*     ..
*     .. Local variables
*
      INTEGER          I,J,K,RESTART
      DOUBLE PRECISION RH, CS, R1, R2, ETIME, PRERR
*
      INFO=0
*
      IDBG = IATPARAM(50)
*
      NORMALFLG=0
*
      BETA0 = 0.0D0
      MFLG = 0
      IF (RATPARAM(35) .NE. 0.0D0) THEN
         PRERR = RERR*RATPARAM(35)
         RESTART = 1
      END IF
!$omp parallel do reduction(+:BETA0)
      DO I=1,N
         BETA0=BETA0+B(I)*B(I)
      END DO
!$omp end parallel do
      BETA0 = DSQRT(BETA0)
      G(1) = BETA
!$omp parallel do
      DO I = 2, M+1
         G(I) = 0.0D0
      END DO
!$omp end parallel do
      DO J = 1, M
*
         IF (MOPE .EQ. 1) THEN
*  COMPUTE w_j = inv(M) A v_j
*          
            CALL OpenATI_DURMV(N,NNZ,IRP,ICOL,VAL,V(1,J),W_J,
     $              IATPARAM,RATPARAM,UINF,LUINF,INFO_A)
            IF (KIND_PRECOND .GE. 5) then
               CALL ILUTSOL('N',N,NNZ,IRP,ICOL,VAL,W_J,W_J,IDIAGP,
     $                 KIND_PRECOND,PRECOND,IATPARAM,RATPARAM,INFO)
            ELSE
               CALL Xabclib_PCSLV(N,NNZ,IRP,ICOL,VAL,W_J,W_J,IDIAGP,
     $              PRECOND,IATPARAM,RATPARAM,INFO)
            END IF
         ELSE
*-------------right preconditioner
            IF (KIND_PRECOND .GE. 5) then
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
         H(J+1,J) = 0.0D0
!$omp parallel do reduction(+:H)
         DO K = 1, N
            H(J+1,J) = H(J+1,J) + W_J(K) * W_J(K)
         END DO
!$omp end parallel do
         H(J+1,J) = DSQRT(H(J+1,J))
*
         IF (H(J+1,J) .EQ. 0.0D0) THEN
            MFLG = 1
         ELSE
*
*  v_j+1 = w_j/h_j+1,j
*
            IF (J .LT. M) THEN
               RH = 1.0D0 / H(J+1,J)
!$omp parallel do
               DO K = 1, N
                  V(K,J+1) = W_J(K) * RH
               END DO
!$omp end parallel do
            END IF
         END IF
************************************************
         R(1,J) = H(1,J)
         DO K = 1, J-1
            R1 =  TMP(1,K) * R(K,J) + TMP(2,K) * H(K+1,J)
            R2 = -TMP(2,K) * R(K,J) + TMP(1,K) * H(K+1,J)
            R(K  ,J) = R1
            R(K+1,J) = R2
         END DO
         CS = DSQRT(R(J,J)**2 + H(J+1,J)**2)
         TMP(1,J) = R(J  ,J) / CS
         TMP(2,J) = H(J+1,J) / CS
         G(J+1) = -TMP(2,J) * G(J)
         G(J  ) =  TMP(1,J) * G(J)
         R(J  ,J) = CS
         R(J+1,J) = 0.0D0
         RERR = dabs(G(J+1)/BETA0)
         IF (IDBG .EQ. 1) then
            WRITE(6,*) ' ==serial update GMRES resid=',J,RERR
         END IF
         IF (RERR .LE. STOP_TOL .OR. MFLG .EQ. 1) THEN
            M0 = J
            go to 200
         END IF
         IF (RESTART .EQ. 1) THEN
            IF (RERR .LT. PRERR) THEN
               M0 = J
               go to 200
            END IF
         END IF
************************************************
      END DO
      M0 = M
  200 CONTINUE
************************************************
*
*  SOLVE  R_m y_m = g_m for y_m
*
      DO I=M0,1,-1
         DO J=I+1,M0
            G(I)=G(I)-R(I,J)*G(J)
         END DO
         IF (R(I,I) .EQ. 0.0D0) THEN
            INFO=200
            write(6,*) ' I,Rii=',I,R(I,I)
            GO TO 9000
         END IF
         G(I)=G(I)/R(I,I)
      END DO
*
*  UPDATE SOLUTION VECTOR X
*
      IF (MOPE.EQ.1) THEN
!$omp parallel do
         DO I=1,N
            DO J=1,M0
               X(I)=X(I)+V(I,J)*G(J)
            END DO
         END DO
!$omp end parallel do
      ELSE
!$omp parallel do
         DO I=1,N
            SUM=0.0D0
            DO J=1,M0
               SUM=SUM+V(I,J)*G(J)
            END DO
            Z(I)=SUM
         END DO
!$omp end parallel do
         IF (KIND_PRECOND.GE.5) then
            CALL ILUTSOL('N',N,NNZ,IRP,ICOL,VAL,Z,Z,IDIAGP,
     $                 KIND_PRECOND,PRECOND,IATPARAM,RATPARAM,INFO)
         ELSE
            CALL Xabclib_PCSLV(N,NNZ,IRP,ICOL,VAL,Z,Z,IDIAGP,
     $          PRECOND,IATPARAM,RATPARAM,INFO)
         END IF
!$omp parallel do
         DO I=1,N
            X(I)=X(I)+Z(I)
         END DO
!$omp end parallel do
      ENDIF
*****************************************************
*
 9000 CONTINUE
*
      RETURN
      END
