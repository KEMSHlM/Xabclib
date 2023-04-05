      SUBROUTINE Xabclib_CG
     $                     (N, NNZ, IRP, ICOL, VAL, B, X,
     $                      PRECOND, NPRE, IATPARAM, RATPARAM,
     $                      WORK, LWORK, INFO)
*
*     Xabclib_CG : 
*                   
      IMPLICIT DOUBLE PRECISION ( A-H , O-Z )
*     ..
*     .. Scalar Auruments ..
      INTEGER           N, NNZ, NPRE, LWORK, INFO
*     ..
*     .. Array Arguments
      INTEGER           IRP(N+1), ICOL(NNZ)
      DOUBLE PRECISION  VAL(NNZ), B(N), X(N), PRECOND(NPRE)
      DOUBLE PRECISION  WORK(LWORK)
*
      INTEGER           IATPARAM(50)
      DOUBLE PRECISION  RATPARAM(50)
*     ..
*
*  =====================================================================
*  Purpose
*  =======
*
*  Xabclib_CG computes the solution to a linear equation system 
*  by Itoh's right-preconditioned ICCG algorithm.
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
*  B       (input) DOUBLE PRECISION array, dimension (N)
*          The right hand side vector.
*
*  X       (input/output) DOUBLE PRECISION array, dimension (N)
*          On entry, An initial guess to the solution vector.
*          On exit,  The approximate solution vector.
*
*  PRECOND (input/output) DOUBLE PRECISION array, dimension (NPRE)
*          On entry,
*          if IATPARAM(24) = 1, Workspace for the preconditioner used.
*          if IATPARAM(24) = 2, Already defined preconditioner data.
*          On exit,
*          Defined preconditioner data.
*
*  NPRE    (input) INTEGER
*          The dimension of PRECOND.
*          if IATPARAM(25) = 1,      NPRE >= 0.
*          if IATPARAM(25) = 2,      NPRE >= N.
*          if IATPARAM(25) = 4,      NPRE >= N+(N+1)/2+1+NZ/2+1+NZ.
*
*  IATPARAM (input/output) Integer array, dimension (50)
*          IATPARAM(3:20) : OpenATI parameter set.
*           IATPARAM(3)     : Number of threads.
*           IATPARAM(7)     : OpenATI_DURMV auto-tuned on/off.
*           IATPARAM(8)     : Fastest OpenATI_DURMV impl. method.
*           IATPARAM(11)    : Columns of Segmented Scan's algorithm.
*           IATPARAM(12)    : Type of Gram-Schmidt procedure.
*           IATPARAM(13)    : DGKS refinement done or not.
*
*          IATPARAM(21:50) : Xabclib parameter set.
*           IATPARAM(21)    : Number of threads.
*           IATPARAM(22)    : Max. iterations.
*           IATPARAM(23)    : Number of iterations.
*           IATPARAM(24)    : Preconditioner operations flag.
*                             1 : not generated yet
*                             2 : already generated
*           IATPARAM(25)    : Preconditioner type.
*                             1 : None
*                             2 : Jacobi preconditioner
*                             4 : ICCG preconditioner
*           IATPARAM(31)    : Mat-Vec times.
*           IATPARAM(32)    : Krylov subspace iteration times.
*           IATPARAM(33)    : The flag of detectiong stagnation.
*           IATPARAM(34)    : Minimmum running iterations.
*           IATPARAM(50)    : Debug print control flag.
*
*  RATPARAM (input/output) DOUBLE PRECISION array, dimension (50)
*           RATPARAM(21:50) : Xabclib parameter set.
*           RATPARAM(22)    : Max. elapsed time(limit time)
*           RATPARAM(23)    : Convergence criterion.
*           RATPARAM(25)    : Preconditioner parameter.
*           RATPARAM(28)    : 2-norm of RHS.
*           RATPARAM(29)    : 2-norm of max. residual.
*           RATPARAM(30)    : Floating operations(*10^9 operations).
*           RATPARAM(31)    : Preconditioner time.
*           RATPARAM(32)    : Total solve time(elapsed).
*           RATPARAM(34)    : Minimum running time.
*
*  WORK    (workspace) DOUBLE PRECISION array, dimension (LWORK)
*
*  LWORK   (input) INTEGER
*          The dimension of WORK.
*          LWORK >= (4)*N
*          
*  INFO    (output) INTEGER
*          Return code.
*          =   0:  Successful exit.
*          <   0:  if INFO = -i, the i-th argument had an illegal value.
*          = 100:  The preconditioner generation process failed.
*          = 200:  The ICCG algorithm breaks down.
*          = 400:  CPU time exceeds to the limit.
*          = 500:  Iteration count exceeds to the limit.
*          = 600:  Stagnation of relative residual.
*          = 700:  Upper limit of NNZ.(if IATPARAM(10)=21)
*          =1000:  Stagnation of relative residual.
*
*  Local variables
*  ===============
*
*  IOPTPC     INTEGER
*             Control flag for generating preconditioner.
*  IDBG       INTEGER
*             Debug print control flag.
*
*  =====================================================================
*     ..
*     .. Local variables
*
      INTEGER           IOPTPC,M
      DOUBLE PRECISION  ETIME1,ETIME2,OMP_GET_WTIME
*
      INTEGER           IDBG, INFO_A
*
      INTEGER           NUM_SMP,IERR_ALLOC
*
      INTEGER           IP_Z,IP_R,IP_P,IP_Q
      INTEGER           IP_IRP2,IP_ICOL2,IP_VAL2,IP_PIV
      INTEGER           DUMMY1,DUMMY2
      DOUBLE PRECISION  DUMMY3
*
      INTEGER           MV_AT_SYM,MVCASE_SYM
      INTEGER           LSINF,IP_WK,IP_SI,LWK
      DOUBLE PRECISION  WK
      ALLOCATABLE::     WK(:)
*
      DOUBLE PRECISION  ETIME0,RERR
*
      EXTERNAL OMP_GET_MAX_THREADS
      INTEGER  OMP_GET_MAX_THREADS
*
      IDBG = IATPARAM(50)
*
*+++++ INITIALIZATION FOR RESTARTED ICCG  ++++++++++++++++++++++++++++++
*
      INFO=0
      IERR_ALLOC=-1
      IATPARAM(31)=0
      IATPARAM(32)=0
*--------------------------Floating Ope. counts
      RATPARAM(30)=0.0D0
*
*  CHECK INPUT PARAMETERS FOR ICCG
*
      CALL Xabclib_CG_PCHK(N,NNZ,IRP,ICOL,VAL,IATPARAM,RATPARAM,
     $                     NPRE,LWORK,INFO)
      IF (INFO .NE. 0) THEN
         IF (IDBG .EQ. 1) THEN
            IF (INFO .EQ. 100)THEN
               WRITE(6,*) ' DIAGONAL ELEMENTS NOT FINED'
            END IF
            WRITE(6,*) ' # PCHK, INFO=',INFO
         END IF
         GO TO 9999
      END IF
      IF (IDBG .EQ. 1) THEN
         WRITE(6,*) ' <<< INPUT PARAMTER CHECK SUCCESSFULLY ENDED.>>>'
         WRITE(6,*) ' '
      END IF
*
*  INITIALIZE PARAMETERS & POINTERS FOR ARRAYS HAVING M-INDEPENDENT SIZE
*
*
      IP_R     = 1
      IP_Z     = IP_R    +  N
      IP_P     = IP_Z    +  N
      IP_Q     = IP_P    +  N
*
      NUM_SMP   = IATPARAM( 3)
*
      MV_AT_SYM = IATPARAM( 7)
*
      IF (IDBG .EQ. 1) THEN
         WRITE(6,*) 'NUM_SMP,MV_AT_SYM=',NUM_SMP,MV_AT_SYM
      ENDIF
************************************************************************
* ---------------------------------------------------------------
* ----- Matrix-vector product auto-tuning.(symmetric) -----------
* ---------------------------------------------------------------
      IF (MV_AT_SYM .EQ. 2 .OR. MV_AT_SYM .EQ. 3) THEN
*
         DO I=1,N
            WORK(I)=1.0D0
         ENDDO
*
         IF (MV_AT_SYM .EQ. 2) THEN
            IATPARAM(7)=2
            LWK=0
            LSINF=INT(0.5*NUM_SMP)+1
            ALLOCATE(WK(LWK+LSINF),STAT=IERR_ALLOC)
         ELSE
            IATPARAM(7)=3
            LWK  =N*NUM_SMP
            LSINF=N+NUM_SMP+3
            ALLOCATE(WK(LWK+LSINF),STAT=IERR_ALLOC)
            IF (IERR_ALLOC .NE. 0) THEN
               IATPARAM(7)=2
               LWK=0
               LSINF=INT(0.5*NUM_SMP)+1
               ALLOCATE(WK(LWK+LSINF),STAT=IERR_ALLOC)
            END IF
         END IF
*
         IP_WK=1
         IP_SI=IP_WK+LWK
*
         CALL OpenATI_DSRMV(N,NNZ,IRP,ICOL,VAL,X,WORK(IP_R),
     &                      IATPARAM,RATPARAM,
     &                      WK(IP_WK),WK(IP_SI),LSINF,INFO_A)
*
         IATPARAM(7)=0
      ELSE
         IATPARAM(7)=0
         MVCASE_SYM=IATPARAM(8)
*
         IF (MVCASE_SYM .EQ. 11) THEN
            LWK=0
            LSINF=0
         ELSE IF (MVCASE_SYM .EQ. 12) THEN
            LWK=0
            LSINF=INT(0.5*NUM_SMP)+1
         ELSE IF (MVCASE_SYM .EQ. 13) THEN
            LWK=N*NUM_SMP
            LSINF=N+NUM_SMP+3
         END IF
         ALLOCATE(WK(LWK+LSINF),STAT=IERR_ALLOC)
         IF (IERR_ALLOC .NE. 0) THEN
            INFO=500
            IF (IDBG .EQ. 1) WRITE(6,*) ' ALLOC ERROR, INFO=',INFO
            GO TO 9999
         END IF
*
         IP_WK=1
         IP_SI=IP_WK+LWK
*
         CALL OpenATI_DSRMV_Setup(N,NNZ,IRP,ICOL,
     &        IATPARAM,RATPARAM,WK(IP_SI),LSINF,INFO_A)
      ENDIF
      IF (IDBG .EQ. 1) THEN
         WRITE(6,*) ' Matrix-vector product auto-tuning finished.'
      END IF
************************************************************************
*
*     IOPCTPC=1 : PRECONDITIONER NOT DEFINDED YET
*     IOPCTPC=0 : PRECONDITIONER ALREADY DEFINDED
*
      IOPTPC       = IATPARAM(24)
      KIND_PRECOND = IATPARAM(25)
*
*  GENERATE PRECONDITIONER FOR CG
*
      ETIME0= OMP_GET_WTIME()
*
      RATPARAM(31)=0.0D0
      IF (IOPTPC .EQ. 1) THEN
         ETIME1= OMP_GET_WTIME()
         IP_PIV = 1
         CALL Xabclib_PCGENS(N,NNZ,IRP,ICOL,VAL,
     $                       PRECOND(IP_PIV),IATPARAM,
     $                       RATPARAM,INFO)
         IF (IDBG .EQ. 1) THEN
            WRITE(6,*) ' ## PCGENS SUCCESSFULLY ENDED.'
         END IF
         ETIME2= OMP_GET_WTIME()
         IF (INFO .NE. 0) THEN
            IF (IDBG .EQ. 1) WRITE(6,*) ' # PCGEN, INFO=',INFO
            GO TO 9999
         END IF
*
* ---------------------preconditioning time
         RATPARAM(31)=ETIME2-ETIME1
      END IF
*
*+++ RESTARTED ICCG MAIN LOOP +++++++++++++++++++++++++++++++++++++++++
*
*
      CALL ICCG (N,NNZ,IRP,ICOL,VAL,B,X,KIND_PRECOND,PRECOND(IP_PIV),
     $           IATPARAM,RATPARAM,WORK(IP_Z),WORK(IP_R),WORK(IP_P),
     $           WORK(IP_Q),WK(IP_WK),WK(IP_SI),LSINF,LWK,NUM_SMP,
     $           IDBG,INFO)
*
      ETIME2= OMP_GET_WTIME()
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*
      RATPARAM(32)=ETIME2-ETIME0
*
      IF (INFO .NE. 0) GO TO 9999
*-----------
      IF (IDBG .EQ. 1) THEN
         WRITE(6,*) '=================================================='
         IF (INFO .EQ. 0) THEN
            WRITE(6,*) ' ## CG PROCESS SUCCESSFULLY ENDED.'
         ELSE
            WRITE(6,*) ' ## CG PROCESS FAILED.'
         END IF
         WRITE(6,6001) IATPARAM(23),RATPARAM(29),RATPARAM(32)
         WRITE(6,*) '=================================================='
      END IF
*
*  EVALUATE ACCURACY OF OBTAINED SOLUTION VECTOR X
*
*      CALL Xabclib_EVAL_SYM(N,NNZ,IRP,ICOL,VAL,X,B,WORK(IP_R),RERR,INFO)
*
*      IF (IDBG .EQ. 1) THEN
*         WRITE(6,*) ' '
*         WRITE(6,*) ' ||b-Ax_i||_2 / ||b||_2    = ',RERR
*      END IF
*
 6000 FORMAT('  ITR =',I4,', ERR =',1PD22.15)
 6001 FORMAT('  NUM OF ITERS =',I4,', ERR =',1PD22.15,', ETIME =',
     $       D10.3,'(SEC)')
 9999 CONTINUE
*
      IF (IERR_ALLOC .EQ. 0) THEN
         DEALLOCATE(WK)
      END IF
*
      RETURN
      END
***********************************************************
      SUBROUTINE ICCG(N,NZ,IRP,ICOL,VAL,B,X,
     $                KIND_PRECOND,PRECOND,IATPARAM,RATPARAM,
     $                Z,R,P,Q,WK,SINF,LSINF,LWK,NUM_SMP,IDBG,INFO)
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*     ..Scalar Arguments
      INTEGER N,NZ,KIND_PRECOND,LSINF,LWK,NUMSMP,IDBG,INFO
*
*     ..Array Arguments
      INTEGER          IRP(N+1),ICOL(NZ)
      DOUBLE PRECISION VAL(NZ),B(N),X(N),PRECOND(N)
      DOUBLE PRECISION Z(N),R(N),P(N),Q(N),WK(LWK),SINF(LSINF)
      INTEGER          IATPARAM(50)
      DOUBLE PRECISION RATPARAM(50)
*
* ======================================================================
* Purpose
* =======
*
* ICCG
*
* Arguments
* =========
*
*  N       (input) INTEGER
*          The order of the matrix.
*
*  NZ      (input) INTEGER
*          Number of nonzero elements in the matrix.
*
*  IRP     (input) INTEGER array, dimension (N+1)
*          The index of the elements in VAL which start a row of the
*          upper traianglar matrix.
*
*  ICOL    (input) INTEGER array, dimension (NNZ)
*          The column indices of each nonzero.
*
*  VAL     (input) DOUBLE PRECISION array, dimension (NNZ)
*          Nonzero elements of the matrix in Compressed Row Storage
*          format.
*
*  B       (input) DOUBLE PRECISION array, dimension (N)
*          The right hand side vector.
*
*  X       (output) DOUBLE PRECISION array, dimension (N)
*          The solution vector of Mx = y.
*
*  KIND_PRECOND (input) INTEGER
*
*  PRECOND (input) DOUBLE PRECISION array, dimension (N)
*          Preconditioner data.
*
*  IATPARAM (input/output) INTEGER array, dimension (50)
*
*  RATPARAM (input/output) DOUBLE PRECISION array, dimension(50)
*
*  Z       (workspace) DOUBLE PRECISION array, dimension (N)
*
*  R       (workspace) DOUBLE PRECISION array, dimension (N)
*
*  P       (workspace) DOUBLE PRECISION array, dimension (N)
*
*  Q       (workspace) DOUBLE PRECISION array, dimension (N)
*
*  WK      (workspase) DOUBLE PRECISION array, dimension (LWK)
*
*  SINF    (input) DOUBLE PRECISION array, dimension (LSINF)
*          Setup Information for Matrix-Vector product Impl.
*
*  LSINF   (input) INTEGER
*          The size of SINF.
*          case IATPARAM(7)=1 :
*             if IATPARAM(8) = 11 :
*                LSINF >= 0
*             else if IATPARAM(8) = 12 :
*                LSINF >= INT(0.5*NUM_SMP)+1
*             else if IATPARAM(8) = 13 :
*                LSINF >= N+NUM_SMP+3
*          case IATPARAM(7)=2 :
*             LSINF >= INT(0.5*NUM_SMP)+1
*          case IATPARAM(7)=3 :
*             LSINF >= N+NUM_SMP+3
*          (NUM_SMP = IATPARAM(3))
*
*  LWK     (input) INTEGER
*          The size of WK.
*
*  NUM_SMP (input) INTEGER
*          Number of threads.
*
*  IDBG    (input) INTEGER
*          Debug print control flag.
*
*  INFO    (output) INTEGER
*
* =====================================================================
*
*
*     ..Local variables
      INTEGER          I,J,ITR,MAX_ITR
      DOUBLE PRECISION rho,rho0,BETA,ALF,td,STOP_TOL
      DOUBLE PRECISION RESID,rnorm,bnorm
      DOUBLE PRECISION MAX_ETIME,ETIME1,ETIME2,ETIME3
*
      DOUBLE PRECISION  ETIME4,ETIME5,MTVC_TIME
      MTVC_TIME = 0.0D0
*
      FLOP=0.0D0
      STOP_TOL  = RATPARAM(23)
      MAX_ETIME = RATPARAM(22)
      MAX_ITER  = IATPARAM(22)
*
      IF (MAX_ITER .LT. 0) MAX_ITER = N
*
      ETIME1 = OMP_GET_WTIME()
      ETIME2 = ETIME1
*
      IF (KIND_PRECOND.EQ.2) THEN
        DO I=1,N
          B(I)=B(I)*PRECOND(I)
        ENDDO
      ENDIF
* r = b - Ax
*
      ETIME4 = OMP_GET_WTIME()
      CALL OpenATI_DSRMV(N,NZ,IRP,ICOL,VAL,X,R,IATPARAM,RATPARAM,
     $                   WK,SINF,LSINF,INFO_A)
      ETIME5 = OMP_GET_WTIME()
      MTVC_TIME = MTVC_TIME + ETIME5 - ETIME4
      DO I = 1, N
         R(I) = B(I) - R(I)
      END DO
*
      CALL XNORM2(N,R,rnorm)
      CALL XNORM2(N,B,bnorm)
      RESID = rnorm /bnorm
      IF (IDBG .EQ. 1) THEN
         WRITE(6,*) 'initial bnorm,rnorm,residual=',bnorm,rnorm,RESID
      END IF
      FLOP=FLOP+5.0D0*N+2.0D0*(2.0D0*NZ-N)
*
* ICCG MAIN LOOP
      DO ITR = 1, MAX_ITER
*
* solve M(Z) = R
         CALL Xabclib_PCSLV_SYM(N,NZ,IRP,ICOL,VAL,
     $                          Z,R,PRECOND,IATPARAM,RATPARAM,WK,
     $                          SINF,LSINF,INFO)
*
         CALL XDDOT(N,R,Z,rho)
         IF (rho .EQ. 0.0D0) THEN
            IF (IDBG .EQ. 1) WRITE(6,*) ' ICCG break down <r,z>' ,rho
            GO TO 9000
         END IF
*
         IF (ITR .EQ. 1) THEN
            DO I = 1, N
               P(I) = Z(I)
            ENd DO
         ELSE
            BETA = rho / rho0
            DO I = 1, N
               P(I) = Z(I) + BETA * P(I)
            END DO
         END IF
*
* Q = AP
         ETIME4 = OMP_GET_WTIME()
         CALL OpenATI_DSRMV(N,NZ,IRP,ICOL,VAL,P,Q,IATPARAM,RATPARAM,
     $                      WK,SINF,LSINF,INFO_A)
*
         ETIME5 = OMP_GET_WTIME()
         MTVC_TIME = MTVC_TIME + ETIME5 - ETIME4
         CALL XDDOT(N,P,Q,td)
         IF (td .EQ. 0.0D0) THEN
            IF (IDBG .EQ. 1) WRITE(6,*) ' ICCG break down <p,q>' ,td
            GO TO 9000
         END IF
         ALF = rho / td
         DO I = 1, N
            X(I) = X(I) + ALF * P(I)
            R(I) = R(I) - ALF * Q(I)
         END DO
         rho0 = rho
*
         CALL XNORM2(N,R,rnorm)
         RESID = rnorm / bnorm
         IF (IDBG .EQ. 1) THEN
            WRITE(6,*) ' ITR=',ITR,' RESID=',RESID
         END IF
*
         ETIME3 = ETIME2
         ETIME2 = OMP_GET_WTIME()
         ETIME3 = ETIME2 - ETIME3
         IF (MAX_ETIME .GT. 0.0D0) THEN
            IF (ETIME2-ETIME1+RATPARAM(31) .GT. MAX_ETIME) THEN
               IF (IDBG .EQ. 1) WRITE(6,*) ' MAX_ETIME OVER'
               INFO = 400
               GO TO 9000
            END IF
         END IF
         IF (RESID .LT. STOP_TOL) THEN
            GO TO 9000
         END IF
      END DO
*
 9000 CONTINUE
      FLOP=FLOP+12.0D0*ITR*N+2.0D0*ITR*(2.0D0*NZ-N)
      IF (RESID .GT. STOP_TOL .AND. ITR .GE. MAX_ITER) THEN
         INFO = 500
      END IF
*****
      IF (KIND_PRECOND.EQ.2) THEN
        DO I=1,N
          X(I)=X(I)*PRECOND(I)
        ENDDO
      ENDIF
      RATPARAM(36)=MTVC_TIME
      RATPARAM(28) = bnorm
      RATPARAM(29) = rnorm / bnorm
      RATPARAM(30) = RATPARAM(30) +FLOP
      IATPARAM(23) = ITR
*
      RETURN
      END
