      SUBROUTINE OpenATI_DAFGS
     &                 (NORMALFLG,N,X,Q,LQ,MM,HR,IATPARAM,RATPARAM,INFO)
*-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
      DOUBLE PRECISION X(N),Q(LQ,MM)
      DOUBLE PRECISION HR(MM)
      DOUBLE PRECISION , ALLOCATABLE :: WK(:)
*
      INTEGER           IATPARAM(50)
      DOUBLE PRECISION  RATPARAM(50)
*
*  =====================================================================
*  Purpose
*  =======
*
*  OpenATI_DAFGS computes a vector "X(N)" orthogonalize 
*  to vectors "Q(1:N,MM)" by various Gram-Schmidt method.
*
*
*      X <- X  |  Q(1:N,MM)
*             --- 
*
*  GS procedure select below;
*      IATPARAM(12) = 0 : CGS (Crassical GS)
*      IATPARAM(12) = 1 : DGKS (Double iterated CGS)
*                         (Daniel-Gragg-Kaufman-Stewart)
*      IATPARAM(12) = 2 : MGS (Modified GS) **default**
*      IATPARAM(12) = 3 : BCGS (Blocked CGS)
*                         Inter-block with MGS,
*                         Intra-block with CGS, block=4
*
*
*  Arguments
*  =========
*
*  NORMALFLG (input) INTEGER
*             Whether output vector normalize or not.
*             NORMALFLG = 0 : not normalize
*                       = 1 : normalize
*
*  N         (input) INTEGER
*             vector length.
*
*  X(N)      (input/output) Double precision
*             input vector.
*
*  Q(LQ,MM)  (input) Double precision
*             Orthonomal vectors in Q(1:N,MM).
*
*  LQ        (input) INTEGER
*             Leading dimension of Q.
*
*  MM        (input) INTEGER
*             the number of vectors in Q.
*
*  HR(MM)    (output) Double precision
*             inner product <X,Q(1:N,M)> in HR(M)
*
*  IDGKS     INTEGER
*             Whether 2nd DGKS iterate running or not.
*             
*  =====================================================================
*
      IDGKS=0
*
      IGSPOL = IATPARAM(12)
*----------------------------------------CRASSICAL GRAM-SCHMIDT
      IF (IGSPOL.EQ.0) THEN
          CALL OpenATI_CGS(N,X,Q,LQ,MM,HR)
          IF (NORMALFLG.EQ.1) THEN
             WNORM=0.0D0
             DO 10 I=1,N
                WNORM=WNORM+X(I)*X(I)
   10        CONTINUE
             WNORM=DSQRT(WNORM)
          END IF
*----------------------------------------DGKS
      ELSE IF (IGSPOL.EQ.1) THEN
          ETA=1.0D0/DSQRT(2.0D0)
*--------------CRASSICAL GRAM-SCHMIDT
          CALL OpenATI_CGS(N,X,Q,LQ,MM,HR)
*
          WNORM=0.0D0
!$omp parallel do reduction(+:WNORM)
          DO 20 I=1,N
             WNORM=WNORM+X(I)*X(I)
   20     CONTINUE
!$omp end parallel do
          WNORM=DSQRT(WNORM)
*
          HNORM=0.0D0
          DO 30 K=1,MM
             HNORM=HNORM+HR(K)*HR(K)
   30     CONTINUE
          HNORM=DSQRT(HNORM)
*------------------------------DGKS-STEP WITH CGS
          IF ( WNORM .LT. HNORM*ETA) THEN
             ALLOCATE(WK(MM),STAT=IERR_DGKS)
             IF (IERR_DGKS .NE.0) THEN
                WRITE(6,*) 
     &           '[OpenATI_DAFGS] DGKS Alloc ERROR'
                GOTO 2000
             END IF
             IDGKS=1
             CALL OpenATI_CGS(N,X,Q,LQ,MM,WK)
 2000        CONTINUE
             IF (NORMALFLG.EQ.1) THEN
                WNORM=0.0D0
                DO 40 I=1,N
                   WNORM=WNORM+X(I)*X(I)
   40           CONTINUE
                WNORM=DSQRT(WNORM)
             END IF
             DEALLOCATE(WK)
          END IF
          IATPARAM(13) = IDGKS
*-------------------------------------Blocked CRASSICAL GRAM-SCHMIDT
      ELSE IF (IGSPOL.EQ.3) THEN
          KMOD=MOD(MM,4)
*--------------Blocked CRASSICAL GRAM-SCHMIDT
          DO 130 K=1,MM-KMOD,4
             S0=0.0D0
             S1=0.0D0
             S2=0.0D0
             S3=0.0D0
!$omp parallel 
!$omp do reduction(+:S0,S1,S2,S3)
             DO 110 I=1,N
                S0=S0+Q(I,K  )*X(I)
                S1=S1+Q(I,K+1)*X(I)
                S2=S2+Q(I,K+2)*X(I)
                S3=S3+Q(I,K+3)*X(I)
  110        CONTINUE
!$omp end do
!$omp do
             DO 120 I=1,N
                X(I)=X(I)-S0*Q(I,K  )-S1*Q(I,K+1)
     &                   -S2*Q(I,K+2)-S3*Q(I,K+3)
  120        CONTINUE
!$omp end do
             HR(K  )=S0
             HR(K+1)=S1
             HR(K+2)=S2
             HR(K+3)=S3
!$omp end parallel
  130     CONTINUE
*
          DO 160 K=MM-KMOD+1,MM
               S0=0.0D0
!$omp parallel
!$omp do reduction(+:S0)
               DO 140 I=1,N
                  S0=S0+Q(I,K  )*X(I)
  140          CONTINUE
!$omp end do
!$omp do
               DO 150 I=1,N
                  X(I)=X(I)-S0*Q(I,K)
  150          CONTINUE
!$omp end do
               HR(K  )=S0
!$omp end parallel
  160     CONTINUE
          IF (NORMALFLG.EQ.1) THEN
             WNORM=0.0D0
             DO 170 I=1,N
                WNORM=WNORM+X(I)*X(I)
  170        CONTINUE
             WNORM=DSQRT(WNORM)
          END IF
*----------------------------------------MGS
      ELSE
*-------------------------ORIGINAL MGS
            DO 220 K=1,MM
               S=0.0D0
!$omp parallel
!$omp do reduction(+:S)
               DO 200 I=1,N
                  S=S+Q(I,K)*X(I)
  200          CONTINUE
!$omp end do
!$omp do
               DO 210 I=1,N
                  X(I)=X(I)-S*Q(I,K)
  210          CONTINUE
!$omp end do
               HR(K)=S
!$omp end parallel
  220       CONTINUE
*
            IF (NORMALFLG.EQ.1) THEN
               WNORM=0.0D0
               DO 230 I=1,N
                  WNORM=WNORM+X(I)*X(I)
  230          CONTINUE
               WNORM=DSQRT(WNORM)
            END IF
*
      END IF
*
      IF (NORMALFLG.EQ.1) THEN
         DO 500 I=1,N
            X(I)=X(I)/WNORM
  500    CONTINUE
      END IF
*
      RETURN
      END
********************************************************************
      SUBROUTINE OpenATI_CGS(N,X,Q,LQ,MM,HR)
*-------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
      DOUBLE PRECISION X(N),Q(LQ,MM)
      DOUBLE PRECISION HR(MM)
*
*  =====================================================================
*  Purpose
*  =======
*
*  OpenATI_CGS computes a vector "X(N)" orthogonalize
*  to vectors "Q(1:N,MM)".
*
*  Arguments
*  =========
*
*  N         (input) INTEGER
*             vector length.
*
*  X(N)      (input/output) Double precision
*             input vector.
*
*  Q(LQ,MM)  (input) Double precision
*             Orthonomal vectors in Q(1:N,MM).
*
*  LQ        (input) INTEGER
*             Leading dimension of Q.
*
*  MM        (input) INTEGER
*             the number of vectors in Q.
*
*  HR(MM)    (output) Double precision
*             inner product <X,Q(1:N,M)> in HR(M)
*
*  =====================================================================
*
      IF(MOD(MM,16).EQ.0) THEN
      KMOD=MOD(MM,4)
      NMOD=MOD(N ,2)

!$omp parallel do
      DO 50 K=1,MM-KMOD,4
         S0=0.0D0
         S1=0.0D0
         S2=0.0D0
         S3=0.0D0
         DO 40 I=1,N
            S0=S0+Q(I,K  )*X(I)
            S1=S1+Q(I,K+1)*X(I)
            S2=S2+Q(I,K+2)*X(I)
            S3=S3+Q(I,K+3)*X(I)
   40    CONTINUE
         HR(K  )=S0
         HR(K+1)=S1
         HR(K+2)=S2
         HR(K+3)=S3
   50 CONTINUE
!$omp end parallel do

      DO 70 K=MM-KMOD+1,MM
         S0=0.0D0
!$omp parallel do reduction(+:S0)
         DO 60 I=1,N
            S0=S0+Q(I,K  )*X(I)
   60    CONTINUE
!$omp end parallel do
         HR(K  )=S0
   70 CONTINUE
*
      DO 120 K=1,MM-KMOD,4
         S0=HR(K)
         S1=HR(K+1)
         S2=HR(K+2)
         S3=HR(K+3)
!$omp parallel do
         DO 110 I=1,N
            X(I)=X(I)-Q(I,K  )*S0-Q(I,K+1)*S1
     &               -Q(I,K+2)*S2-Q(I,K+3)*S3
  110    CONTINUE
!$omp end parallel do
  120 CONTINUE
*
      DO 140 K=MM-KMOD+1,MM
         S0=HR(K)
!$omp parallel do
         DO 130 I=1,N
            X(I)=X(I)-Q(I,K  )*S0
  130    CONTINUE
!$omp end parallel do
  140 CONTINUE
*
      ELSE
*
      DO 75 K=1,MM
         S0=0.0D0
!$omp parallel do reduction(+:S0)
         DO 65 I=1,N
            S0=S0+Q(I,K  )*X(I)
  65    CONTINUE
!$omp end parallel do
         HR(K  )=S0
!$omp parallel do
         DO 135 I=1,N
            X(I)=X(I)-Q(I,K  )*S0
  135    CONTINUE
!$omp end parallel do
  75  CONTINUE
*
      END IF
*
*
*
      RETURN
      END
