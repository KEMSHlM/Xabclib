      SUBROUTINE Xabclib_GMRES2 (N, NNZ, IRP, ICOL, VAL, B, X,
     $                          PRECOND, NPRE, IATPARAM,RATPARAM,
     $                          WORK, LWORK,INFO)
*
*     Xabclib_GMRES : Compute the solution to a linear equation system
*                     by preconditioned restarted GMRES algorithm.
*     January 2012
*                   
      IMPLICIT DOUBLE PRECISION ( A-H , O-Z )
*     ..
*     .. Scalar Auruments ..
      INTEGER           N, NNZ, KIND_PRECOND, NPRE, MSIZE, LWORK, INFO
*     ..
*     .. Array Arguments
      INTEGER           IRP(N+1), ICOL(NNZ)
      INTEGER           IATPARAM(50)
      DOUBLE PRECISION  RATPARAM(50)
      DOUBLE PRECISION  VAL(NNZ), B(N), X(N), PRECOND(NPRE)
      DOUBLE PRECISION  WORK(LWORK)
*     ..
*
*  =====================================================================
*  Purpose
*  =======
*
*  Xabclib_GMRES computes the solution to a linear equation system 
*  by preconditioned restarted GMRES algorithm.
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
*          if IATPARAM(25) = 1,        NPRE >= 0.
*          if IATPARAM(25) = 2,3 or 4, NPRE >= N.
*          if IATPARAM(25) = 5,    NPRE >= 3*NNZ/2+2*N+50.
*          if IATPARAM(25) = 6,    NPRE >= 3*(2.0*IFILL+1)*N/2+3*N+50.
*                                  (IFILL=IATPARAM(26))
*
*  IATPARAM (input/output) Integer array, dimension (50)
*          IATPARAM(3:20) : OpenATI parameter set.
*           IATPARAM(3)     : Number of threads.
*           IATPARAM(4)     : Flag of Krylov subspace expand by
*                             MM-ratio
*           IATPARAM(5)     : Value of Krylov subspace expand.
*           IATPARAM(9)     : OpenATI_DURMV auto-tuned on/off.
*           IATPARAM(10)    : Fastest OpenATI_DURMV impl. method.
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
*                             3 : SSOR   preconditioner
*                             4 : ILU(0)diag  preconditioner
*                             5 : ILU(0) preconditioner
*                             6 : ILUT(p,t) preconditioner
*           IATPARAM(26)    : Fill-in of ILUT.
*           IATPARAM(27)    : Input size of Krylov subspace.
*           IATPARAM(28)    : Start size of Krylov subspace at
*                             subspace expand AT-on.
*           IATPARAM(29)    : Final size of Krylov subspace.
*           IATPARAM(31)    : Mat-Vec times.
*           IATPARAM(32)    : Kryrov subspace iteration times.
*           IATPARAM(33)    : The flag of detectiong stagnation.
*           IATPARAM(34)    : Minimmum running iterations.
*           IATPARAM(50)    : Debug print control flag.
*
*  RATPARAM (input/output) DOUBLE PRECISION array, dimension (50)
*          RATPARAM(3:20) : OpenATI parameter set.
*           RATPARAM(4)     : Threshhold of MM-ratio.
*           RATPARAM(5)     : Value of MM-ratio.
*
*          RATPARAM(21:50) : Xabclib parameter set.
*           RATPARAM(22)    : Max. elapsed time(limit time)
*           RATPARAM(23)    : Convergence criterion.
*           RATPARAM(25)    : Preconditioner parameter.
*           RATPARAM(28)    : 2-norm of RHS.
*           RATPARAM(29)    : 2-norm of max. residual.
*           RATPARAM(30)    : Floating operations(*10^9 operations).
*           RATPARAM(31)    : Preconditioner time.
*           RATPARAM(32)    : Total solve time(elapsed).
*           RATPARAM(34)    : Minimum running time.
*           RATPARAM(35)    : Restart condition.
*            
*  MSIZE   (input) INTEGER
*          The maximum size of Krylov subspace.
*  
*  WORK    (workspace) DOUBLE PRECISION array, dimension (LWORK)
*
*  LWORK   (input) INTEGER
*          The dimension of WORK.
*          LWORK>=(MSIZE+2)*N + 2*(MSIZE+1)*(MSIZE+1) + 3*(MSIZE+1) + (N-1)/2+1
*          
*  INFO    (output) INTEGER
*          Return code.
*          =   0:  Successful exit.
*          <   0:  if INFO = -i, the i-th argument had an illegal value.
*          = 100:  The preconditioner generation process failed.
*          = 200:  The GMRES algorithm breaks down.
*          = 300:  OpenATI_DAFRT input data invalid.
*          = 400:  CPU time exceeds to the limit.
*          = 500:  Iteration count exceeds to the limit.
*          = 600:  Not Enough free memory.(if IATPARAM(10)=12,13,21)
*          = 700:  Upper limit of NNZ.(if IATPARAM(10)=21)
*          =1000:  Stagnation of relative residual.
*
*  Local variables
*  ===============
*
*  IOPTPC     INTEGER
*             Control flag for generating preconditioner.
*  MAX_ITER   INTEGER
*             Maximum number of restart iterations.
*  M          INTEGER
*             The size of Krylov subspace.
*  STOP_TOL   DOUBLE PRECISION
*             Convergence criterion.
*  MAX_ETIME  DOUBLE PRECISION
*             Maximum elapsed time for the solver.
*  BETA       DOUBLE PRECISION
*             The square norm of the residual vector.
*             ( BETA = || inv(M) R0 ||_2 )
*  BETA0      DOUBLE PRECISION
*             The square norm of the right hand side vector.
*             ( BETA0 = || inv(M) b ||_2 )
*  ICASE      INTEGER
*             Selected Algorithm NO. of matrix-vector product.
*  IDBG       INTEGER
*             Debug print control flag.
*  SAMP       DOUBLE PRECISION array, dimension (5)
*             Sampling accuracy of approximate solution vectors.
*  ISTGCNT    INTEGER
*             The counter for detcting stagnation.
*
*  =====================================================================
*     ..
*     .. Local variables
*
      INTEGER           IOPTPC,MAX_ITER,M
      DOUBLE PRECISION  STOP_TOL,MAX_ETIME,BETA,BETA0
      DOUBLE PRECISION  ETIME1,ETIME2,OMP_GET_WTIME
*
      INTEGER           ICASE,IDBG
*
      DOUBLE PRECISION  SAMP(5)
*
      DOUBLE PRECISION  UINF
      ALLOCATABLE ::    UINF(:)
      INTEGER           LUINF,NUM_SMP,JL,IERR_ALLOC
*
      EXTERNAL OMP_GET_MAX_THREADS
      INTEGER  OMP_GET_MAX_THREADS
*
      INTEGER           ISTG,ISTGCNT
      DOUBLE PRECISION  ETIME3,ETIME4,PERR,EMA,MIN_ETIME
*
      IDBG = IATPARAM(50)
*
*+++++ INITIALIZATION FOR RESTARTED GMRES ++++++++++++++++++++++++++++++
*
      INFO=0
      IERR_ALLOC=-1
      IRT=0
      IATPARAM(31)=0
      IATPARAM(32)=0
      ISTGCNT=0
      MPTH=IATPARAM(6)
      ISTG=IATPARAM(33)
      PERR=1.0D0
      EMA=0.0D0
*--------------------------Floating Ope. counts
      RATPARAM(30)=0.0D0
*
*  CHECK INPUT PARAMETERS FOR GMRES(M)
*
      CALL Xabclib_GMRES_PCHK(N,NNZ,IRP,ICOL,VAL,IATPARAM,RATPARAM,
     $                        NPRE,LWORK,INFO)
      IF (INFO .NE. 0) THEN
         IF (IDBG .EQ. 1) WRITE(6,*) ' # PCHK, INFO=',INFO
         GO TO 9999
      END IF
      IF (IDBG .EQ. 1) THEN
         WRITE(6,*) ' <<< INPUT PARAMTER CHECK SUCCESSFULLY ENDED.>>>'
         WRITE(6,*) ' '
      END IF
*
*  INITIALIZE PARAMETERS & POINTERS FOR ARRAYS HAVING M-INDEPENDENT SIZE
*
      STOP_TOL   = RATPARAM(23)
      MAX_ITER   = IATPARAM(22)
          IF (MAX_ITER .LT. 0) MAX_ITER = N
      MAX_ETIME  = RATPARAM(22)
      NUM_SMP    = IATPARAM( 3)
*
      IP_IDIAGP=1
      IP_Z     =IP_IDIAGP+ (N-1)/2+1
      IP_WJ    =IP_Z     +  N
      IP_V     =IP_WJ    +  N
*
*  SET PERFORMANCE TUNING PARAMTERS
*
*
      MSIZE_AT  = IATPARAM( 4)
      MSIZE     = IATPARAM(27)
      MSIZE_INI = IATPARAM(28)
*------------------------OpenATI_DAFRT option
      DO I=1,5
         SAMP(I)=1.0D0
      ENDDO
*+++++++++++++++++++++++++++++ Msize AT
      IF (MSIZE_AT .EQ. 1) THEN
*-- Tuning Point ---
         IF (MSIZE_INI .LE. 0) THEN
            M=2
         ELSE
            M=MSIZE_INI
         END IF
*-------------------
         IF (M .GT. MSIZE) M=MSIZE
         MOLD=M
      ELSE
         M=MSIZE
      END IF
*
*+++++++++++++++++++++++++++++ MatVec AT
      MV_AT_USYM = IATPARAM( 9)
      JL         = IATPARAM(11)
      IF (MV_AT_USYM .EQ. 2 .OR. MV_AT_USYM .EQ. 3) THEN
         IF (MV_AT_USYM .EQ. 2) THEN
            IATPARAM( 9) = 2
            LUINF=INT(0.5*NUM_SMP)+1
            ALLOCATE(UINF(LUINF),STAT=IERR_ALLOC)
         ELSE
            IATPARAM( 9) = 3
*
            LUINF=INT(1.5*N)+INT(4.25*JL)+10
            ALLOCATE(UINF(LUINF),STAT=IERR_ALLOC)
            IF (IERR_ALLOC .NE. 0) THEN
               IATPARAM( 9) = 2
               LUINF=INT(0.5*NUM_SMP)+1
               ALLOCATE(UINF(LUINF),STAT=IERR_ALLOC)
            END IF
         END IF
*
         CALL OpenATI_DURMV(N,NNZ,IRP,ICOL,VAL,X,WORK(IP_Z),
     $           IATPARAM,RATPARAM,UINF,LUINF,INFO_A)
         IATPARAM( 9) = 0
      ELSE
         IATPARAM( 9) = 0
         MVCASE_USYM  = IATPARAM(10)
*
         IF (MVCASE_USYM .EQ. 11) THEN
            LUINF=0
            ICASE=11
         ELSE IF (MVCASE_USYM .EQ. 12) THEN
            LUINF=INT(0.5*NUM_SMP)+1
            ICASE=12
         ELSE IF (MVCASE_USYM .EQ. 13) THEN
            LUINF=INT(1.5*N)+INT(4.25*JL)+10
            ICASE=13
         ELSE IF (MVCASE_USYM .EQ. 21) THEN
            LUINF=INT(1.125*NNZ)+INT(2.125*JL)+10
            ICASE=21
         END IF
         ALLOCATE(UINF(LUINF),STAT=IERR_ALLOC)
         IF (IERR_ALLOC .NE. 0) THEN
            INFO=600
            GO TO 9999
         END IF
         CALL OpenATI_DURMV_Setup(N,NNZ,IRP,IATPARAM,RATPARAM,
     &           UINF,LUINF,INFO_A)
      END IF
*
*  GENERATE POINTER LIST FOR DIAGONAL ELEMENTS OF MATRIX A
*
      CALL Xabclib_IDIAGPGEN(N,NNZ,IRP,ICOL,WORK(IP_IDIAGP),INFO)
      IF (INFO .NE. 0) THEN
        IF (IATPARAM(25).NE.6 .AND. IATPARAM(25).NE.1) THEN
         IF (IDBG .EQ. 1) WRITE(6,*) ' # IDIAGPGEN, INFO=',INFO
         GO TO 9999
        ELSE
         INFO=0
        END IF
      END IF
*
*     IOPCTPC=0 : PRECONDITIONER NOT DEFINDED YET
*     IOPCTPC=1 : PRECONDITIONER ALREADY DEFINDED
*
      IOPTPC       = IATPARAM(24)
      KIND_PRECOND = IATPARAM(25)
* ----------------------precond operation ( 1: left, else: right)
      MOPE = 2
*
*  GENERATE PRECONDITIONER FOR GMRES(M)
*
      ETIME0= OMP_GET_WTIME()
*
      RATPARAM(31)=0.0D0
      IF (IOPTPC .EQ. 1) THEN
         ETIME1= OMP_GET_WTIME()
         CALL Xabclib_PCGENE(N,NNZ,IRP,ICOL,VAL,WORK(IP_IDIAGP),
     $                       PRECOND,IATPARAM,RATPARAM,INFO)
         ETIME2= OMP_GET_WTIME()
         IF (INFO .NE. 0) THEN
            IF (IDBG .EQ. 1) WRITE(6,*) ' # PCGEN, INFO=',INFO
            GO TO 9999
         END IF
* ---------------------preconditioning time
         RATPARAM(31)=ETIME2-ETIME1
      END IF
*
      FLOP=0.0D0
*+++ RESTARTED GMRES MAIN LOOP +++++++++++++++++++++++++++++++++++++++++
*
      ETIME1= OMP_GET_WTIME()
      ETIME2= ETIME1
*
      ISAMP = 0
*
*  beta_0 = || inv(M) b ||_2
*
      IF (MOPE.EQ.1) THEN
         IF (KIND_PRECOND.GE.5) then
            CALL ILUTSOL('N',N,NNZ,IRP,ICOL,VAL,WORK(IP_Z),B,
     $                   WORK(IP_IDIAGP),KIND_PRECOND,PRECOND,
     $                   IATPARAM,RATPARAM,INFO)
         ELSE
            CALL Xabclib_PCSLV(N,NNZ,IRP,ICOL,VAL,WORK(IP_Z),B,
     $                         WORK(IP_IDIAGP),PRECOND,IATPARAM,
     $                         RATPARAM,INFO)
         END IF
      ELSE
         DO I=1,N
            WORK(IP_Z+I-1)=B(I)
         END DO
      ENDIF
      BETA0=0.0D0
!$omp parallel do reduction(+:BETA0)
      DO I=1,N
         BETA0=BETA0+WORK(IP_Z+I-1)*WORK(IP_Z+I-1)
      END DO
!$omp end parallel do
      IF (BETA0 .EQ. 0.0D0) THEN
         INFO=-6
      END IF
      BETA0=DSQRT(BETA0)
         
      DO ITRCNT=1,MAX_ITER+1
*
*  r_0 = A x_0
*
         CALL OpenATI_DURMV(N,NNZ,IRP,ICOL,VAL,X,WORK(IP_Z),IATPARAM,
     $                      RATPARAM,UINF,LUINF,INFO_A)
*
*  r_0 = b - A x_0
*
         DO I=1,N
            WORK(IP_Z+I-1)=B(I)-WORK(IP_Z+I-1)
         END DO
*
*  r_0 = inv(M) r_0 = inv(M)(b - A x_0)
*
         IF (MOPE.EQ.1) THEN
               IF (KIND_PRECOND.GE.5) then
                  CALL ILUTSOL('N',N,NNZ,IRP,ICOL,VAL,WORK(IP_Z),
     $                         WORK(IP_Z),WORK(IP_IDIAGP),KIND_PRECOND,
     $                         PRECOND,IATPARAM,RATPARAM,INFO)
               ELSE
                  CALL Xabclib_PCSLV(N,NNZ,IRP,ICOL,VAL,WORK(IP_Z),
     $                               WORK(IP_Z),WORK(IP_IDIAGP),PRECOND,
     $                               IATPARAM,RATPARAM,INFO)
               END IF
         ENDIF
*
*  beta = || r_0 ||_2
*
         BETA=0.0D0
!$omp parallel do reduction(+:BETA)
         DO I=1,N
            BETA=BETA+WORK(IP_Z+I-1)*WORK(IP_Z+I-1)
         END DO
!$omp end parallel do
         BETA=DSQRT(BETA)
*
*  CALCULATE RELATIVE ERROR
*
         IF (ITRCNT .EQ. 1) THEN
            RERR = BETA/BETA0
            IF (RERR .LE. STOP_TOL) GOTO 1000
         END IF
*
*  v_1 = r_0 / beta
*
         IF (BETA .NE. 0.0D0) THEN
            RBETA=1.0D0/BETA
            DO I=1,N
               WORK(IP_V+I-1)=WORK(IP_Z+I-1)*RBETA
            END DO
         END IF
*
         FLOP=FLOP+2.0D0*NNZ+6.0D0*N
         IATPARAM(31)=IATPARAM(31)+1
         IF (INFO .NE. 0) THEN 
            IF (IDBG .EQ. 1) WRITE(6,*) ' # INIT, INFO=',INFO
            GO TO 9999
         END IF
*
         IF (ISAMP .EQ. 5) ISAMP=0
         ISAMP=ISAMP+1
         SAMP(ISAMP)=RERR
*
         IF (MSIZE_AT .EQ. 1) THEN
            CALL OpenATI_DAFRT(5,SAMP,IRT,IATPARAM,RATPARAM,INFO_A)
            IF (IRT .EQ. 1) THEN
*-- Tuning Point ---
               MOLD=M
               M=M+IATPARAM(5)
*-------------------
               IF (M .GT. MSIZE) THEN
                M=MSIZE
               ENDIF
               IF (IDBG .EQ. 1) THEN
                  WRITE(6,*) ' OpenATI_DAFRT Subspace expand',M
               ENDIF
            ENDIF
         END IF
*
         ETIME3= ETIME2
         ETIME2= OMP_GET_WTIME()
         ETIME3= ETIME2 - ETIME3
*-- Detect Stagnation ---
*
        IF(ISTG.EQ.1) THEN
         CALL OpenATI_DAFSTG(ISTGCNT,EMA,RERR,PERR,STOP_TOL,
     $    ITRCNT-1,MAX_ITER,ETIME2-ETIME1+RATPARAM(31),ETIME3,MAX_ETIME,
     $    IATPARAM,RATPARAM,INFO_A)
         IF ((ISTGCNT.GE.MPTH .AND.
     $       ETIME2-ETIME1+RATPARAM(31) .GT. RATPARAM(34) .AND.
     $       ITRCNT-1 .GT. IATPARAM(34)) .OR.
     $        ISTGCNT.GT.MPTH*10) THEN
          INFO=1000
          IF (IDBG.EQ.1) WRITE(6,*)
     $    '    Stagnation of rerative residual occurred.'
          GO TO 1000
         ENDIF
         PERR=RERR
        ENDIF
*--------------------------
*
*  SET POINTERS TO WORK ARRAYS
*
         IP_H   = IP_V   +  N*M
         IP_G   = IP_H   + (M+1)*M
         IP_R   = IP_G   +  M+1
         IP_TMP = IP_R   + (M+1)*M
         IWRKSZ = IP_TMP +  M*2    -1
*
*  ALNOLDI - MODIFIED GRAM-SCHMIDT
*  COMPUTATE THE MINIMIZER OF ||beta e_1 - H_m y||2 AND
*  UPDATE SOLUTION VECTOR X
*
         CALL Xabclib_GMRES_MGS_MIN(N,NNZ,IRP,ICOL,VAL,X,B,
     $           WORK(IP_IDIAGP),KIND_PRECOND,MOPE,PRECOND,IATPARAM,
     $           RATPARAM,WORK(IP_H),M,M0,WORK(IP_V),WORK(IP_WJ),
     $           WORK(IP_Z),WORK(IP_G),WORK(IP_R),WORK(IP_TMP),BETA,
     $           RERR,STOP_TOL,UINF,LUINF,INFO)
         IATPARAM(31)=IATPARAM(31)+M0
         IATPARAM(32)=IATPARAM(32)+M0
         FLOP=FLOP+4.0D0*NNZ*M0+4.0D0*N*M0+2.0D0*M0*M0*N
         FLOP=FLOP+2.0D0*M0*N+2.0D0*NNZ
         IF (INFO .NE. 0) THEN
            IF (IDBG .EQ. 1) WRITE(6,*) ' # MGS_MIN,  INFO=',INFO
            GO TO 9999
         END IF
         IF (IDBG .EQ. 1) THEN
            WRITE(6,6000) ITRCNT,RERR
         END IF
         IF (RERR.LE.STOP_TOL .OR. ITRCNT.EQ.MAX_ITER) GO TO 1000
         ETIME4= OMP_GET_WTIME()
         IF (MAX_ETIME .GT. 0.0D0) THEN
            IF (ETIME4-ETIME1+RATPARAM(31) .GT. MAX_ETIME) THEN
               INFO=400
               GO TO 1000
            END IF
         END IF
      END DO
*
 1000 CONTINUE
*
      IF (IDBG .EQ. 1) THEN
         WRITE(6,*)"last M:",M0
      END IF
      ETIME2= OMP_GET_WTIME()
      IF (RERR.GT.STOP_TOL .AND. ITRCNT.EQ.MAX_ITER) THEN
         INFO=500
      END IF
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*
      IATPARAM(23)=ITRCNT
      IATPARAM(29)=MSIZE
      IF (MSIZE_AT .EQ. 1) THEN
         IF (IRT .EQ. 1) THEN
            IATPARAM(29)=MOLD
         ELSE
            IATPARAM(29)=M
         END IF
      END IF
      RATPARAM(32)=ETIME2-ETIME0
      RATPARAM(28)=BETA0
      RATPARAM(29)=RERR
      FLOP=FLOP+(ITRCNT)*(24.0D0*N+8.0D0*NNZ)
      RATPARAM(30)=1.0D-9*(RATPARAM(30)+FLOP)
*
      IF (INFO .NE. 0) GO TO 9999
*-----------
      IF (IDBG .EQ. 1) THEN
         WRITE(6,*) '=================================================='
         IF (INFO .EQ. 0) THEN
            WRITE(6,*) ' ## GMRES2 PROCESS SUCCESSFULLY ENDED.'
         ELSE
            WRITE(6,*) ' ## GMRES2 PROCESS FAILED.'
         END IF
         WRITE(6,6001) ITRCNT,RERR,ETIME2-ETIME1
         WRITE(6,*) '=================================================='
      END IF
*
*  EVALUATE ACCURACY OF OBTAINED SOLUTION VECTOR X
*
      CALL Xabclib_EVAL(N,NNZ,IRP,ICOL,VAL,X,B,WORK(IP_Z),RERR,INFO)
*
      IF (IDBG .EQ. 1) THEN
         WRITE(6,*) ' '
         WRITE(6,*) ' ||b-Ax_i||_2 / ||b||_2    = ',RERR
      END IF
*
 6000 FORMAT('  ITR =',I4,', ERR =',1PD22.15)
 6001 FORMAT('  NUM OF ITERS =',I4,', ERR =',1PD22.15,', ETIME =',
     $       D10.3,'(SEC)')
 9999 CONTINUE
*
      IF (IERR_ALLOC .EQ. 0) THEN
         DEALLOCATE(UINF)
      END IF
*
      RETURN
      END
