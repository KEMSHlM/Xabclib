      SUBROUTINE OpenATI_EIGENSOLVE(N,NNZ,IRP,ICOL,VAL,NEV,EV,EVEC,
     $                              IATPARAM,RATPARAM,INFO)
*
*     OpenATI_EIGENSOLVE : Compute the solution to a eigen system
*                          with policy controller.
*
*     January 2010
*     December 2011
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*     ..
*     .. Scalar Auruments ..
      INTEGER           N, NNZ, NEV, INFO
*     ..
*     .. Array Arguments
      INTEGER           IRP(N+1), ICOL(NNZ)
      DOUBLE PRECISION  VAL(NNZ), EV(*), EVEC(*)
*
      INTEGER           IATPARAM(50)
      DOUBLE PRECISION  RATPARAM(50)
*     ..
*  =====================================================================
*  Purpose
*  =======
*              Compute the solution to a eigen system with policy 
*              controller.
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
*  NEV     (input) INTEGER
*          The number of requirement eigen-pairs.
*
*  EV      (output) DOUBLE PRECISION array, dimension (NEV)
*  (NEV)   Eigenvalues.
*
*  EVEC    (output) DOUBLE PRECISION array, dimension (N,NEV)
*  (N,NEV) Eigenvectors.
*
*  IATPARAM (input/output) Integer array, dimension (50)
*          Integer parameter set for OpenATI and Xabclib.
*           IATPARAM(1:2)   : Mandatory.
*
*          IATPARAM(3:20) : OpenATI parameter set.
*           IATPARAM(3)     : Number of threads.
*           IATPARAM(4)     : Flag of Krylov subspace expand by
*                             MM-ratio
*           IATPARAM(5)     : Value of Krylov subspace expand.
*           IATPARAM(6)     : Reserved.
*           IATPARAM(7)     : OpenATI_DSRMV auto-tuned on/off.
*           IATPARAM(8)     : Fastest OpenATI_DSRMV impl. method.
*           IATPARAM(9)     : OpenATI_DURMV auto-tuned on/off.
*           IATPARAM(10)    : Fastest OpenATI_DURMV impl. method.
*           IATPARAM(11)    : Columns of Segmented Scan's algorithm.
*           IATPARAM(12)    : Type of Gram-Schmidt procedure.
*           IATPARAM(13)    : DGKS refinement done or not.
*           IATPARAM(14)    : Access to meminfo.
*           IATPARAM(15)    : Number of retried solver.
*           IATPARAM(16)    : Total restart of solver.
*           IATPARAM(17)    : Total Mat-Vec times.
*           IATPARAM(18:20) : Reserved.
*
*          IATPARAM(21:50) : Xabclib parameter set.
*           IATPARAM(21)    : Number of threads.
*           IATPARAM(22)    : Max. iterations.
*           IATPARAM(23)    : Number of iterations.
*           IATPARAM(24)    : Preconditioner operations flag.
*           IATPARAM(25)    : Preconditioner type.
*           IATPARAM(26)    : Fill-in of ILUT.
*           IATPARAM(27)    : Input size of Krylov subspace.
*           IATPARAM(28)    : Start size of Krylov subspace at
*                             subspace expand AT-on.
*           IATPARAM(29)    : Final size of Krylov subspace.
*           IATPARAM(30)    : Eigenvalue order option
*           IATPARAM(31)    : Mat-Vec times.
*           IATPARAM(32)    : The flag of detectiong stagnation.
*           IATPARAM(33)    : Minimmum running iterations.
*           IATPARAM(34:49) : Reserved.
*           IATPARAM(50)    : Debug print control flag.
*
*  RATPARAM (input/output) DOUBLE PRECISION array, dimension (50)
*          Double precision parameter set for OpenATI and Xabclib.
*           RATPARAM(1:2)   : Mandatory.
*
*          RATPARAM(3:20) : OpenATI parameter set.
*           RATPARAM(3)     : Reserved.
*           RATPARAM(4)     : Threshhold of MM-ratio.
*           RATPARAM(5)     : Value of MM-ratio.
*           RATPARAM(6:13)  : Reserved.
*           RATPARAM(14)    : Residual norm.
*           RATPARAM(15)    : Set-up time.
*           RATPARAM(16)    : Preconditioner time.
*           RATPARAM(17)    : Solver time.
*           RATPARAM(18)    : Total time.
*           RATPARAM(19:20) : Reserved.
*
*          RATPARAM(21:50) : Xabclib parameter set.
*           RATPARAM(21)    : Reserved.
*           RATPARAM(22)    : Max. elapsed time(limit time)
*           RATPARAM(23)    : Convergence criterion.
*           RATPARAM(24)    : Reserved.
*           RATPARAM(25)    : Preconditioner parameter.
*           RATPARAM(26:27) : Reserved.
*           RATPARAM(28)    : 2-norm of RHS.
*           RATPARAM(29)    : 2-norm of max. residual.
*           RATPARAM(30)    : Floating operations(*10^9 operations).
*           RATPARAM(31)    : Preconditioner time.
*           RATPARAM(32)    : Total solve time(elapsed).
*           RATPARAM(33)    : Minimum running time.
*           RATPARAM(34:50) : Reserved.
*
*  INFO    (output) INTEGER
*          Return code.
*          =    0:  Successful exit
*          = -100:  Invalid separate word('=') in keyword.
*          = -200:  Invalid the number of threads value at POLICY FILE.
*          = -300:  Invalid in keyword 'POLICY' in POLICY FILE.
*          = -400:  Not enough space.
*          = -500:  Cannot allocate working-set space.
*  =====================================================================
*     ..
*     .. Local variables
*
      INTEGER  NUM_POL, LP_SMPS, LP_POLICY, LP_SOLVER
      INTEGER  MP_MAXMEM, MP_RESID, MP_MAXTIME
      INTEGER  I, K
      INTEGER  MM, ISIZE, MSIZE, LIP, LOOP, IREST, IORDER
*
      INTEGER  ITYPE, NBACK, NUMCPU, INFO_A
      INTEGER  IO, IDBG, IDUMMY, IERR_ALLOC, IERR_DEALLOC
*
      INTEGER  TIME, MEMORY, ACCURACY, STABLE
      INTEGER  SOLVER_XABCLIB_LANCZOS
      INTEGER  SOLVER_XABCLIB_ARNOLDI
*
      DOUBLE PRECISION USEMEM
      DOUBLE PRECISION ETIMES, TIMEMAX, ETIME0, ETIME1
      DOUBLE PRECISION S, S2, SS, SS2, X, Z, RES
      DOUBLE PRECISION DUMMY1, DUMMY2, DUMMY3
*
      CHARACTER*8   POLICY(5),GSPOL(0:4)
      CHARACTER*10  DAY(2)
      CHARACTER*32  SOLVER
*
      PARAMETER (NUM_POL=20)
*
      INTEGER           THREAD_NUM
      CHARACTER*16      THREAD_NUM_CHA
*
      INTEGER          IPOLICY(NUM_POL)
      DOUBLE PRECISION RPOLICY(NUM_POL)
      CHARACTER*32     SOLVERLIST(NUM_POL)
*
      INTEGER           OMP_GET_NUM_THREADS, OMP_GET_THREAD_NUM
      DOUBLE PRECISION  OMP_GET_WTIME
*
      PARAMETER (LP_SMPS   = 1,
     $           LP_POLICY = 2, TIME=1, MEMORY=2, ACCURACY=3, STABLE=4,
     $           LP_SOLVER = 4, SOLVER_XABCLIB_LANCZOS=11,
     $                          SOLVER_XABCLIB_ARNOLDI=12,
     $           MP_MAXMEM = 1,
     $           MP_RESID  = 2,
     $           MP_MAXTIME= 3)
*
      DOUBLE PRECISION ,ALLOCATABLE :: WK1(:)
*
      POLICY(1)='TIME'
      POLICY(2)='MEMORY'
      POLICY(3)='ACCURACY'
      POLICY(4)='STABLE'
*
      GSPOL(0)='CGS'
      GSPOL(1)='DGKS'
      GSPOL(2)='MGS'
      GSPOL(3)='BCGS'
*
      SOLVERLIST(11)='XABCLIB_LANCZOS'
      SOLVERLIST(12)='XABCLIB_ARNOLDI'
*
*=======================================================================
*
      CALL DATE_AND_TIME(DAY(1),DAY(2))
      ETIMES= OMP_GET_WTIME()
*
      THREAD_NUM= OMP_GET_THREAD_NUM()
      WRITE(THREAD_NUM_CHA,*) THREAD_NUM
*
      IO=16+THREAD_NUM
      OPEN(IO,FILE='OPENATI_POLICY_REPORT.'//ADJUSTL(THREAD_NUM_CHA))
      WRITE(IO,7000)
      WRITE(IO,7010)
      WRITE(IO,7011) DAY(1)(1:4),DAY(1)(5:8),
     $               DAY(2)(1:2),DAY(2)(3:4)
      WRITE(IO,7000)
      WRITE(IO,*)
 7000 FORMAT(1H ,'****************************************************')
 7010 FORMAT(1H ,'*****  OpenATI  EIGEN  SOLVER POLICY REPORT    *****')
 7011 FORMAT(1H ,'*****',21X,A4,'.',A4,'   ',A2,':',A2,'    *****')
*
      WRITE(IO,*) ' [Environment variables]'
*----------------------Debug print control
      IDBG=IATPARAM(50)
      WRITE(IO,*) '    OPENATI_DEBUG  = ',IDBG
*
*-----------------------------------------------------------------------
*         READ POLICY DEFINITION FILE
*-----------------------------------------------------------------------
      ITYPE=10
      CALL OPENATI_READ_POLICY(ITYPE,IPOLICY,RPOLICY,IDBG,INFO)
      IF (INFO.NE.0) THEN
         WRITE(IO,*) '   <<<< ERROR CODE >>>> = ',INFO
         GOTO 9900
      END IF
*-----------------------------------------------------------------------
*         SIZE FOR ALLOCATE
*-----------------------------------------------------------------------
      SOLVER = SOLVERLIST(IPOLICY(LP_SOLVER))
*
      CALL OPENATI_PRESET(IATPARAM,SOLVER,IPOLICY,RPOLICY,N,NNZ,NEV,
     $                    MM,ISIZE,IDBG,INFO)
*
      WRITE(IO,*) ' [Policy Definitions]'
      WRITE(IO,*) '   POLICY         = ',POLICY(IPOLICY(LP_POLICY)) 
      WRITE(IO,*) '   SMPs           = ',IPOLICY(LP_SMPS)
      WRITE(IO,*) '   SOLVER         = ',SOLVER
      WRITE(IO,*) '   REQUIREMENT WORKING MEMORY = ',RPOLICY(MP_MAXMEM)
      WRITE(IO,*) '     <<< Upper Bound 16GBYTE >>>'
      WRITE(IO,*) '   REQUIREMENT RESIDUAL       = ',RPOLICY(MP_RESID)
      WRITE(IO,*) '   REQUIREMENT MAX. TIME      = ',RPOLICY(MP_MAXTIME)
      WRITE(IO,*)
      WRITE(IO,*) '   MAX. SUBSPACE SIZE   =',MM
*
      IF (INFO.NE.0) THEN
         WRITE(IO,*) '   <<<< ERROR CODE >>>> = ',INFO
         GOTO 9900
      END IF

      ALLOCATE(WK1(ISIZE),STAT=IERR_ALLOC)
      IF (IERR_ALLOC .NE.0 ) THEN
         INFO=-500
         WRITE(IO,*) '   <<<< ERROR CODE >>>> = ',INFO
         GOTO 9900
      END IF
*
      USEMEM = 1.0D0*ISIZE
      USEMEM = USEMEM*8.0D-9
      WRITE(IO,7020) USEMEM
 7020 FORMAT(1H ,'   RUNTIME MEMORY USE   =',F10.2,' [GBYTE]')
*
*-----------------------------------------------------------------------
*         PARAMETER SETTING
*-----------------------------------------------------------------------
      NBACK = OMP_GET_NUM_THREADS()
      NUMCPU = IPOLICY(LP_SMPS   )
      IATPARAM(3)=NUMCPU
      IATPARAM(21)=NUMCPU
      CALL OMP_SET_NUM_THREADS(NUMCPU)
*
      TIMEMAX = RPOLICY(MP_MAXTIME)
*---------------------MAX. RESTARTS
      IATPARAM(22) = MIN(10000,N)
*---------------------REQUIREMENT ACCURACY
      RATPARAM(23) = RPOLICY(MP_RESID  )
*---------------------Maximum elapsed time for the solver
      RATPARAM(22) = TIMEMAX
*
*---------------------KRYLOV SUBSPACE SIZE FROM OPENATI_PRESET
      MSIZE=MM
*---------------------POLICY = STABLE
      IF (IPOLICY(LP_POLICY) .EQ. STABLE ) THEN
         IF (MSIZE.GE.5*NEV) THEN
            IATPARAM(27) = 5*NEV
         ELSE
            INFO=-400
            GOTO 9700
         ENDIF 
         IATPARAM(28) = IATPARAM(27)
         IATPARAM(4) = 0
         IATPARAM(12) = 2
         IF      (SOLVER .EQ.'XABCLIB_LANCZOS') THEN
            IATPARAM(7) = 0
            IATPARAM(8) = 12
         ELSE IF (SOLVER .EQ.'XABCLIB_ARNOLDI') THEN
            IATPARAM(9) = 0
            IATPARAM(10)= 12
         END IF
      ELSE
         IATPARAM(27)=MM
*---------------------INIT VALUE OF M
         IATPARAM(28) = MIN(2*NEV,IATPARAM(27))
*---------------------KRYLOV SUBSPACE EXPAND AT ON
         IATPARAM(4) = 1
*---------------------FASTEST MATRIX-VECTOR PRODUCT AT ON
         IF      (SOLVER .EQ.'XABCLIB_LANCZOS') THEN
            IF (IPOLICY(LP_POLICY) .EQ. MEMORY ) THEN
               IATPARAM(7) = 2
            ELSE
               IATPARAM(7) = 3
            END IF
         ELSE IF (SOLVER .EQ.'XABCLIB_ARNOLDI') THEN
            IF (IPOLICY(LP_POLICY) .EQ. MEMORY ) THEN
               IATPARAM(9) = 2
            ELSE
               IATPARAM(9) = 3
            END IF
         END IF
*=======================================================================
*            POLICY CONTROL
*=======================================================================
*-------------<<POLICY 1>>: GramSchmidt STRATEGY
         IF (IPOLICY(LP_POLICY   ).EQ. TIME) THEN
*----------------------------------BCGS
            IATPARAM(12) = 3
         ELSE
            IF (RPOLICY(MP_RESID  ) .LE. 1.0D-10) THEN
*----------------------------------MGS
               IATPARAM(12) = 2
            ELSE
*----------------------------------BCGS
               IATPARAM(12) = 3
            END IF
         END IF
*
      END IF
*
      WRITE(IO,*)
      IF      (SOLVER .EQ.'XABCLIB_LANCZOS') THEN
         WRITE(IO,7030) IATPARAM(4),IATPARAM(7)
      ELSE IF (SOLVER .EQ.'XABCLIB_ARNOLDI') THEN
         WRITE(IO,7030) IATPARAM(4),IATPARAM(9)
      END IF
 7030 FORMAT(1H ,'   KRYLOV SUBSPACE EXPAND AT = ',I2,
     $           '   ,MATVEC AT = ',I2)
*
      WRITE(IO,*) '   Initial Gram-Schmidt Strategy = ',
     $            GSPOL(IATPARAM(12))
      WRITE(IO,*)
*------------Integer WORK POINTER
      IF      (SOLVER .EQ.'XABCLIB_LANCZOS') THEN
         LIP=(MSIZE+1)*N+2*MSIZE*MSIZE+7*MSIZE+5*NEV+3
      ELSE IF (SOLVER .EQ.'XABCLIB_ARNOLDI') THEN
         LIP=(MSIZE+5)*N+5*MSIZE*MSIZE+9*MSIZE+6*NEV+1
      END IF
*
      ETIME0= OMP_GET_WTIME()
      ETIME0= ETIME0 - ETIMES
*
      RES=-1.0D0
      ETIME1= 0.0D0
      LOOP=1
      IREST=0
*
      IATPARAM(17)=0
*
 3000 CONTINUE
*
      IF      (SOLVER .EQ.'XABCLIB_LANCZOS') THEN
         CALL Xabclib_LANCZOS( N, NNZ, IRP, ICOL, VAL,
     $                         NEV, EV, EVEC, N,
     $                         IATPARAM, RATPARAM,
     $                         WK1, LIP-1, WK1(LIP), 5*MSIZE+3, INFO )
      ELSE IF (SOLVER .EQ.'XABCLIB_ARNOLDI') THEN
         CALL Xabclib_ARNOLDI( N, NNZ, IRP, ICOL, VAL,
     $                         NEV, EV, EVEC, N,
     $                         IATPARAM, RATPARAM,
     $                         WK1, LIP-1, WK1(LIP), MSIZE, INFO )
      END IF
*
      ETIME1= ETIME1 + RATPARAM(32)
      IREST=IREST+IATPARAM(23)
*
      IATPARAM(17)=IATPARAM(17) + IATPARAM(31)
*
      RATPARAM(14)=RATPARAM(29)
*
      IF (INFO.NE.0 ) THEN
         GOTO 9700
      END IF
*
      IF (TIMEMAX .GT. 0  .AND. ETIME0+ETIME1 .GT. TIMEMAX) THEN
         GOTO 9700
      END IF
*-------------<<POLICY 2>>: ACCURACY
      IF (IPOLICY(LP_POLICY   ).EQ. ACCURACY) THEN
         RES=-1.0D0
*-------------------RESIDUAL CHECK
         IF      (SOLVER .EQ.'XABCLIB_LANCZOS') THEN
            CALL OpenATI_RESID(N,NNZ,IRP,ICOL,VAL,EV,NEV,EVEC,
     $                                WK1,RES)
         ELSE IF (SOLVER .EQ.'XABCLIB_ARNOLDI') THEN
            CALL OpenATI_RESIDZ(N,NNZ,IRP,ICOL,VAL,EV,NEV,EVEC,
     $                                WK1,RES)
         END IF
         RATPARAM(14)=RES
         IF (IDBG.EQ.1)
     $     WRITE(6,*) '[POLICY:ACCURACY] RESIDUAL=',RES
*
         IF (RES .LT. RPOLICY(MP_RESID  ) ) THEN
            GOTO 9000
         ELSE
*----------------------------------DGKS
            IATPARAM(12) = 1
            RATPARAM(22) = TIMEMAX -ETIME1 -ETIME0
            RATPARAM(23) = RATPARAM(23)*0.1D0
            LOOP=LOOP+1
            GOTO 3000
         END IF
      ELSE
         GOTO 9000
      END IF
*
*--------------------------------------SUCCESSFULY ENDING
 9000 CONTINUE
      WRITE(IO,*) ' ====== OPENATI_EIGENSOLVE SUCCESSFULY ENDED ======'
      GOTO 9800
*--------------------------------------ERROR AFTER ALLOCATE
 9700 CONTINUE
      WRITE(IO,*) ' ****** OPENATI_EIGENSOLVE ERROR DETECTED *********'
      WRITE(IO,*) '   <<<< ERROR CODE >>>> = ',INFO
*
 9800 CONTINUE
      IATPARAM(15)=LOOP
      IATPARAM(16)=IREST
*
      RATPARAM(15)=ETIME0
      RATPARAM(16)=0.0D0
      RATPARAM(17)=ETIME1
      RATPARAM(18)=ETIME0 + ETIME1
*
      WRITE(IO,*)
      WRITE(IO,*) ' [OPENATI_EIGENSOLVE RESULT]'
      WRITE(IO,*) '   MATRIX DATA : N= ',N,' NNZ= ',NNZ
      IF      (SOLVER .EQ.'XABCLIB_LANCZOS') THEN
      WRITE(IO,*) '   FASTEST MATVEC NO. = ',IATPARAM(8)
      ELSE IF (SOLVER .EQ.'XABCLIB_ARNOLDI') THEN
      WRITE(IO,*) '   FASTEST MATVEC NO. = ',IATPARAM(10)
      END IF
      WRITE(IO,*) '   FINAL KRYLOV SUBSPACE SIZE =',IATPARAM(29)
      WRITE(IO,*) '   FINAL Gram-Schmidt Strategy = ',
     $                     GSPOL(IATPARAM(12))
      IF      (SOLVER .EQ.'XABCLIB_LANCZOS') THEN
      WRITE(IO,*) '   NUMBER OF RETRYED LANCZOS=',LOOP
      WRITE(IO,*) '   TOTAL RESTARTS of LANCZOS=',IREST
      ELSE IF (SOLVER .EQ.'XABCLIB_ARNOLDI') THEN
      WRITE(IO,*) '   NUMBER OF RETRYED ARNOLDI=',LOOP
      WRITE(IO,*) '   TOTAL RESTARTS of ARNOLDI=',IREST
      END IF
      IF (IPOLICY(LP_POLICY   ).EQ. ACCURACY)
     $WRITE(IO,*) '   MAXIMUM RESIDUAL NORM    =',RES
      WRITE(IO,*) '   SET-UP TIME              =',ETIME0,' [SEC]'
      WRITE(IO,*) '   SOLVER TIME              =',ETIME1,' [SEC]'
      WRITE(IO,*)
      WRITE(IO,*) '   TOTAL  TIME              =',ETIME0+ETIME1,' [SEC]'
*
      DEALLOCATE(WK1,STAT=IERR_DEALLOC)
      IF (IERR_DEALLOC.NE.0) THEN
         WRITE(6,*) '[Error] IERR_DEALLOC=',IERR_DEALLOC
      END IF
*--------------------------------------PARAMETER ERROR
 9900 CONTINUE
      RETURN
      END
*
*
*
      SUBROUTINE OpenATI_RESID(N,NNZ,IRP,ICOL,VAL,EV,NEV,EVEC,
     $                                WK1,RES)
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*     ..
*     ..Scalar Argument..
      INTEGER           N,NNZ,NEV
      DOUBLE PRECISION  RES
*     ..
*     ..Array Argument
      INTEGER           IRP(N+1),ICOL(NNZ)
      DOUBLE PRECISION  VAL(NNZ)
*
      DOUBLE PRECISION  EV(NEV), EVEC(N,NEV), WK1(N)
*
*  =====================================================================
*
      INTEGER           I,J_PTR,K
      DOUBLE PRECISION  S, S2, SS, SS2, X, Z
*
      DO K=1,NEV
         CALL OpenATI_DSRMV_11(N,NNZ,IRP,ICOL,VAL,EVEC(1,K),WK1)
         S =0.0D0
         S2=0.0D0
         DO I=1,N
            X =WK1(I)-EV(K)*EVEC(I,K)
            S =S +X*X
            S2=S2+WK1(I)*WK1(I)
         ENDDO
         SS =DSQRT(S )
         SS2=DSQRT(S2)
         Z  =SS/SS2
         RES=DMAX1(RES,Z)
      END DO
*
      RETURN
      END
*
*
      SUBROUTINE OpenATI_RESIDZ(N,NNZ,IRP,ICOL,VAL,EV,NEV,EVEC,
     $                                WK1,RES)
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*     ..
*     ..Scalar Argument..
      INTEGER           N,NNZ,NEV
      DOUBLE PRECISION  RES
*     ..
*     ..Array Argument
      INTEGER           IRP(N+1),ICOL(NNZ)
      DOUBLE PRECISION  VAL(NNZ)
*
      COMPLEX*16        EV(NEV), EVEC(N,NEV), WK1(N)
*
*  =====================================================================
*
      INTEGER           I,J_PTR,K
      DOUBLE PRECISION  S1, S2, SS1, SS2,Z
      COMPLEX*16        S, X
*
      DO K=1,NEV
!$omp parallel
!$omp do private(S,J_PTR,I)
         DO I=1,N
            S=DCMPLX(0.0D0,0.0D0)
            DO J_PTR=IRP(I),IRP(I+1)-1
               S=S+VAL(J_PTR)*EVEC(ICOL(J_PTR),K)
            END DO
            WK1(I)=S
         END DO
!$omp end do
!$omp end parallel
*
         S1=0.0D0
         S2=0.0D0
         DO I=1,N
            X =WK1(I)-EV(K)*EVEC(I,K)
            S1=S1+X*DCONJG(X)
            S2=S2+WK1(I)*DCONJG(WK1(I))
         ENDDO
         SS1=DSQRT(S1)
         SS2=DSQRT(S2)
         Z  =SS1/SS2
         RES=DMAX1(RES,Z)
      END DO
*
      RETURN
      END
