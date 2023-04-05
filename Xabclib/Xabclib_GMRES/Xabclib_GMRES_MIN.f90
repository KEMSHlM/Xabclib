      SUBROUTINE Xabclib_GMRES_MIN (N,H,M,M0,G,V,BETA,X,
     $              NNZ,IRP,ICOL,VAL,IDIAGP,
     $              KIND_PRECOND,MOPE,PRECOND,IATPARAM,RATPARAM,Z,
     $              INFO)
      IMPLICIT DOUBLE PRECISION ( A-H , O-Z )
*     ..
*     .. Scalar Auruments ..
      INTEGER           N, M, M0, INFO
      DOUBLE PRECISION  BETA
*     ..
*     .. Array Arguments
      DOUBLE PRECISION  H(M+1,M), G(M+1), V(N,M), X(N)
*
      INTEGER           IRP(N+1), ICOL(NNZ), IDIAGP(N)
      INTEGER           KIND_PRECOND
      DOUBLE PRECISION  VAL(NNZ), Z(N)
*
      INTEGER           IATPARAM(50)
      DOUBLE PRECISION  RATPARAM(50)
*     ..
*
*  =====================================================================
*  Purpose
*  =======
*
*  Xabclib_GMRES_MIN computes the minimizer y_m of
*  ||beta e_1 - H_m y||_2 and updates the approximate solution
*  vector x (x_m = x_0 + V_m y_m).
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          The order of the matrix.
*
*  H       (input) DOUBLE PRECISION array, dimension (M+1,M)
*          The Hessenberg matrix H_m.
*
*  M       (input) INTEGER
*          The size of Krylov subspace.
*
*  M0      (input) INTEGER
*          The step count j of the Arnoldi loop at which the GMRES
*          algorithm breaks down.
*          if the GMRES algorithm does not break down, M0 = M.
*
*  G       (workspace) DOUBLE PRECISION array, dimension (M+1)
*          
*  V       (input) DOUBLE PRECISION array, dimension (N,M)
*          The Arnoldi vectors (v_1,v_2, ..., v_m) which form an ortho-
*          normal basis of the Krylov subspace
*          K_m = span{v_1,A v_1, ... ,A^(m-1) v_1}.
*
*  BETA    (input) DOUBLE PRECISION
*          The square norm of the residual vector.
*          ( BETA = || inv(M) (b - A x_0) ||_2 )
*
*  X       (input/output) DOUBLE PRECISION array, dimension (N)
*          On entry, The approximate solution vector not updated.
*          On exit,  The approximate solution vector updated.
*
*  INFO    (output) INTEGER
*          Return code.
*
*  =====================================================================
*     ..
*     .. Local variables
*
      INTEGER          I,J
      DOUBLE PRECISION WK,SI,CI,H1,H2,G1,G2
*
*
      INFO=0
      MMAX=M
      M=M0
*
*  TRANSFORM THE HESSENBERG MATRIX INTO UPPER TRIANGULAR FORM
*
      G(1)=BETA
      DO I=2,M+1      
         G(I)=0.0D0
      END DO
      DO I=1,M
         IF (H(I+1,I) .NE. 0.0D0) THEN  
            WK=1.0D0/DSQRT(H(I,I)*H(I,I)+H(I+1,I)*H(I+1,I))
            SI=H(I+1,I)*WK
            CI=H(I  ,I)*WK
            DO J=I,M
               H1=H(I  ,J)
               H2=H(I+1,J)
               H(I  ,J)= CI*H1+SI*H2
               H(I+1,J)=-SI*H1+CI*H2
            END DO
            G1=G(I  )
            G2=G(I+1)
            G(I  )= CI*G1+SI*G2
            G(I+1)=-SI*G1+CI*G2
         END IF
      END DO
*
*  SOLVE  R_m y_m = g_m for y_m
*
      DO I=M,1,-1
         DO J=I+1,M
            G(I)=G(I)-H(I,J)*G(J)
         END DO
         IF (H(I,I) .EQ. 0.0D0) THEN
            INFO=200
            write(6,*) ' I,Hii=',I,H(I,I)
            GO TO 9000
         END IF
         G(I)=G(I)/H(I,I)
      END DO
*
*  UPDATE SOLUTION VECTOR X   
*
      IF (MOPE.EQ.1) THEN
!$omp parallel do private(J)
         DO I=1,N
            DO J=1,M
               X(I)=X(I)+V(I,J)*G(J)
            END DO
         END DO
!$omp end parallel do
      ELSE
!$omp parallel do private(SUM,J)
         DO I=1,N
            SUM=V(I,1)*G(1)
            DO J=2,M
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
*
 9000 CONTINUE
      RETURN
      END
