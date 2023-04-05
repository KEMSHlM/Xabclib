      SUBROUTINE OpenATI_LINEARSOLVE(N,NNZ,IRP,ICOL,VAL,B,X,IATPARAM,
     $           RATPARAM,INFO)
*
*     OpenATI_LINEARSOLVE : Compute the solution to a linear equation
*                           system with policy controller.
*
*     January 2010
*     December 2011
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*     ..
*     .. Scalar Auruments ..
      INTEGER           N, NNZ, INFO
*     ..
*     .. Array Arguments
      INTEGER           IRP(N+1), ICOL(NNZ)
      DOUBLE PRECISION  VAL(NNZ), B(N), X(N)
*
      INTEGER           IATPARAM(50)
      DOUBLE PRECISION  RATPARAM(50)
*     ..
*  =====================================================================
*  Purpose
*  =======
*              Compute the solution to a linear equation system with
*              policy controller.
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
*  IATPARAM (input/output) Integer array, dimension (50)
*          Integer parameter set for OpenATI and Xabclib.
*           IATPARAM(1:2)   : Mandatory.
*
*          IATPARAM(3:20) : OpenATI parameter set.
*           IATPARAM(3)     : Number of threads.
*           IATPARAM(4)     : Flag of Krylov subspace expand by
*                             MM-ratio
*           IATPARAM(5)     : Value of Krylov subspace expand.
*           IATPARAM(6)     : Threshold for judging stagnation.
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
*           IATPARAM(32)    : Krylov subspace iteration times.
*           IATPARAM(33)    : The flag of detectiong stagnation.
*           IATPARAM(34)    : Minimmum running iterations.
*           IATPARAM(35:49) : Reserved.
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
*           RATPARAM(6)     : Exponen.
*           RATPARAM(7:13)  : Reserved.
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
*           RATPARAM(34)    : Minimum running time.
*           RATPARAM(34:50) : Reserved.
*
*  INFO    (output) INTEGER
*          Return code.
*          =    0:  Successful exit
*          = -100:  Invalid separate word('=') in keyword.
*          = -200:  Invalid the number of threads value at POLICY FILE.
*          = -300:  Invalid in keyword 'POLICY' in POLICY FILE.
*          = -310:  Invalid in keyword 'PRECONDITIONER' in POLICY FILE.
*          = -400:  Not enough space.
*          = -500:  Cannot allocate working-set space.
*  =====================================================================
*     ..
*     .. Local variables
*
      INTEGER  NUM_POL, LP_SMPS, LP_POLICY, LP_PRECON, LP_SOLVER
      INTEGER  MP_MAXMEM, MP_RESID, MP_MAXTIME
      INTEGER  TIME,MEMORY,ACCURACY,STABLE
      INTEGER  NO, JACOBI, SSOR, ILU0D, ILU0, ILUT, AUTO
      INTEGER  SOLVER_XABCLIB_GMRES
      INTEGER  SOLVER_XABCLIB_BICGSTAB
      INTEGER  SOLVER_XABCLIB_BICGSTAB_ORG
*
      INTEGER  I, LOOP, IREST, ITYPE
      INTEGER  IO, IDBG
      INTEGER  NEV, MM, ISIZE, NBACK, NUMCPU
      INTEGER  KIND_PRECOND, MSIZE, NPRE
      INTEGER  INFO_A, IERR_ALLOC, IERR_DEALLOC
*
      DOUBLE PRECISION ETIMES, ETIME0, ETIME1, TIMEMAX, TIMELIM
      DOUBLE PRECISION ETIMES2,ETIMEP
      DOUBLE PRECISION ETIMEA, ETIMEB, ETIMEC
      DOUBLE PRECISION USEMEM
      DOUBLE PRECISION BNORM, RES, S, SS

      CHARACTER*8   POLICY(5),PRECON(10),GSPOL(0:4)
      CHARACTER*10  DAY(2)
      CHARACTER*32  SOLVER
      CHARACTER*20  CVER
*
      PARAMETER (NUM_POL=20)
*
      INTEGER           THREAD_NUM
      CHARACTER*16      THREAD_NUM_CHA
*
      INTEGER           IPOLICY(NUM_POL)
      DOUBLE PRECISION  RPOLICY(NUM_POL)
      CHARACTER*32      SOLVERLIST(NUM_POL)
*
      INTEGER           OMP_GET_NUM_THREADS
      INTEGER           OMP_GET_THREAD_NUM
      DOUBLE PRECISION  OMP_GET_WTIME
*
      DOUBLE PRECISION ,ALLOCATABLE :: WK1(:), WK2(:), ANS(:)
*
      INTEGER           PS_PAIR(2,100), IPRELIST(100)
      DOUBLE PRECISION  RPRELIST(100)
      INTEGER           TOVER
*
      PARAMETER (LP_SMPS   = 1,
     $           LP_POLICY = 2, TIME=1, MEMORY=2, ACCURACY=3, STABLE=4,
     $           LP_PRECON = 3, NO=1, JACOBI=2, SSOR=3, ILU0D=4,
     $                          ILU0=5, ILUT=6, AUTO=10,
     $           LP_SOLVER = 4, SOLVER_XABCLIB_GMRES=1,
     $                          SOLVER_XABCLIB_BICGSTAB=2,
     $                          SOLVER_XABCLIB_BICGSTAB_ORG=3,
     $           MP_MAXMEM = 1,
     $           MP_RESID  = 2,
     $           MP_MAXTIME= 3)
*
      POLICY(1)='TIME'
      POLICY(2)='MEMORY'
      POLICY(3)='ACCURACY'
      POLICY(4)='STABLE'
*
      PRECON(1) ='NONE'
      PRECON(2) ='JACOBI'
      PRECON(3) ='SSOR'
      PRECON(4) ='ILU0D'
      PRECON(5) ='ILU0'
      PRECON(6) ='ILUT'
      PRECON(10)='AUTO'
*
      GSPOL(0)='CGS'
      GSPOL(1)='DGKS'
      GSPOL(2)='MGS'
      GSPOL(3)='BCGS'
*
      SOLVERLIST(1) ='XABCLIB_GMRES'
      SOLVERLIST(2) ='XABCLIB_BICGSTAB'
      SOLVERLIST(3) ='XABCLIB_BICGSTAB_ORG'
      SOLVERLIST(10)='AUTO'
*
*=======================================================================
*---------------------get OpenATI version info.
      CALL OpenATI_GET_VER(CVER)
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
      WRITE(IO,*) '  <<<  OpenATI Version  ','" ',CVER,' " >>>'
      WRITE(IO,*)
 7000 FORMAT(1H ,'****************************************************')
 7010 FORMAT(1H ,'*****  OpenATI  LINEAR SOLVER POLICY REPORT    *****')
 7011 FORMAT(1H ,'*****',21X,A4,'.',A4,'   ',A2,':',A2,'    *****')
*
      WRITE(IO,*) ' [Environment variables]'
*----------------------Debug print control
      IDBG= IATPARAM(50)
      WRITE(IO,*) '    OPENATI_DEBUG  = ',IDBG
*
*-----------------------------------------------------------------------
*         READ POLICY DEFINITION FILE
*-----------------------------------------------------------------------
      ITYPE=1
      CALL OPENATI_READ_POLICY(ITYPE,IPOLICY,RPOLICY,IDBG,INFO)
      IF (INFO.NE.0) THEN
         WRITE(IO,*) '   <<<< ERROR CODE >>>> = ',INFO
         GOTO 9900
      END IF
*
      KIND_SOLVER  = IPOLICY(LP_SOLVER)
      SOLVER       = SOLVERLIST(KIND_SOLVER)
      KIND_PRECOND = IPOLICY(LP_PRECON)
*
      WRITE(IO,*) ' [Policy Definitions]'
      WRITE(IO,*) '   POLICY         = ',POLICY(IPOLICY(LP_POLICY))
      WRITE(IO,*) '   SMPs           = ',IPOLICY(LP_SMPS)
      WRITE(IO,*) '   SOLVER         = ',SOLVER
      WRITE(IO,*) '   PRECONDITIONER = ',PRECON(IPOLICY(LP_PRECON))
      WRITE(IO,*) '   REQUIREMENT WORKING MEMORY = ',RPOLICY(MP_MAXMEM)
      WRITE(IO,*) '     <<< Upper Bound 16GBYTE >>>'
      WRITE(IO,*) '   REQUIREMENT RESIDUAL       = ',RPOLICY(MP_RESID)
      WRITE(IO,*) '   REQUIREMENT MAX. TIME      = ',RPOLICY(MP_MAXTIME)
      WRITE(IO,*)
*
      IF (KIND_SOLVER.EQ.10 .OR. KIND_PRECOND.EQ.10) THEN
         IFILL=IATPARAM(26)
         TAU=RATPARAM(25)
         CALL OPENATI_SET_STRATEGY(KIND_SOLVER,KIND_PRECOND,IFILL,TAU,
     $                             NTRY,PS_PAIR,IPRELIST,RPRELIST)
*----------------------------STAGNATION DETECT ON
         IATPARAM(33)=1
         WRITE(IO,*) '   STAGNATION DETECT          = ',IATPARAM(33)
         WRITE(IO,*) ' ++++++++++++++++++++++++++++++++++++++++++++++++
     $++++++'
         WRITE(IO,*) '   NUMBER OF AVAILABLE CHOICE = ',NTRY
         DO I=1,NTRY
            I1=PS_PAIR(1,I)
            I2=PS_PAIR(2,I)
            I3=IPRELIST(I)
            R1=RPRELIST(I)
            IF (I1.EQ.6) THEN
               WRITE(IO,7013) I,I3,R1,SOLVERLIST(I2)
 7013       FORMAT(1H ,'   TRIAL-',I2,' ILUT(',I3,',',D10.2') - ',
     $A32)
            ELSE
               WRITE(IO,7014) I,PRECON(I1),SOLVERLIST(I2)
 7014       FORMAT(1H ,'   TRIAL-',I2,5X,A8,'-',A32)
            END IF
         ENDDO
      ELSE
         NTRY=1
         PS_PAIR(1,NTRY)=IPOLICY(LP_PRECON)
         PS_PAIR(2,NTRY)=IPOLICY(LP_SOLVER)
         IPRELIST(NTRY) =IATPARAM(26)
         RPRELIST(NTRY) =RATPARAM(25)
      END IF
*
      TIMEMAX = RPOLICY(MP_MAXTIME)
      TIMELIM = TIMEMAX
      IF ( TIMEMAX .GT. 0.0D0 ) THEN
         TRYTIME = TIMEMAX/NTRY/4
      ELSE
         TRYTIME = -1.0D0
      END IF
      RATPARAM(34) = TRYTIME
      TOVER=0
*
      NUMCPU = IPOLICY(LP_SMPS   )
      IATPARAM(3) = NUMCPU
      IATPARAM(21) = NUMCPU
      CALL OMP_SET_NUM_THREADS(NUMCPU)
*
*------------BNORM
      ALLOCATE(ANS(N),STAT=IERR_ALLOC)
      IF (IERR_ALLOC .NE.0 ) THEN
         INFO=-500
         WRITE(IO,*) '   <<<< ERROR CODE >>>> = ',INFO
         GOTO 9900
      END IF
      CALL OpenATI_DURMV_11(N,NNZ,IRP,ICOL,VAL,X,ANS)
      BNORM=0.0D0
      S=0.0D0
      DO I=1,N
         BNORM=BNORM+B(I)*B(I)
         S=S+(ANS(I)-B(I))*(ANS(I)-B(I))
      ENDDO
      BNORM=DSQRT(BNORM)
      IF (BNORM.LT.1.0D-15) BNORM=1.0D0
      SS=DSQRT(S)
      GOODANS=SS/BNORM
*
      ANS(1:N)=X(1:N)
*
      ETIME0= OMP_GET_WTIME()
      ETIMES= ETIME0 - ETIMES
      ETIMEP= 0.0D0
      ETIME0= 0.0D0
      ETIME1= 0.0D0
      TIMELIM = TIMELIM - ETIMES
      IRETRY=0
      IRETRYMM=IATPARAM(28)
*[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[
      DO 100 ITRY=1,NTRY+1
*[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[
      IF (ITRY.EQ.NTRY+1) THEN
         IF (IRETRY.NE.0) THEN
            IATPARAM(33)=0
         ELSE
            GOTO 200   
         ENDIF
      ENDIF
*
      ETIMEA= OMP_GET_WTIME()
*
      KIND_PRECOND = PS_PAIR(1,ITRY)
      IFILL=IPRELIST(ITRY)
      TAU  =RPRELIST(ITRY)
      KIND_SOLVER  = PS_PAIR(2,ITRY)
      SOLVER       = SOLVERLIST(KIND_SOLVER)
      IATPARAM(25)=KIND_PRECOND
*
      IF ( NTRY.NE.1) THEN
         WRITE(IO,*)
         WRITE(IO,*) ' ++++++++++++++++++++++++++++++++++++++++++++++++
     $++++++'
         WRITE(IO,*) ' << SELECT PRECOND. & SOLVER  >> TRIAL=',ITRY
         WRITE(IO,*) '  *SOLVER         = ',SOLVER
         WRITE(IO,*) '  *PRECONDITIONER = ',PRECON(KIND_PRECOND)
         IF ( KIND_PRECOND.EQ.3 ) THEN
            WRITE(IO,*) '  * Relaxation Parameter= ',RPRELIST(ITRY)
         END IF
         IF ( KIND_PRECOND.EQ.4 ) THEN
            WRITE(IO,*) '  * Break Down Threshold= ',RPRELIST(ITRY)
         END IF
         IF ( KIND_PRECOND.EQ.6 ) THEN
            WRITE(IO,7015) IPRELIST(ITRY),RPRELIST(ITRY)
 7015       FORMAT(1H ,'  * ILUT(p,tau) ,p=',I3,', tau=',D12.5)
         END IF
         WRITE(IO,*) ' ++++++++++++++++++++++++++++++++++++++++++++++++
     $++++++'
         IF (IDBG.EQ.1) THEN
         WRITE(6,*)
         WRITE(6,*) ' ++++++++++++++++++++++++++++++++++++++++++++++++
     $++++++'
         WRITE(6,*) ' << SELECT PRECOND. & SOLVER  >> TRIAL=',ITRY
         WRITE(6,*) '  *SOLVER         = ',SOLVER
         WRITE(6,*) '  *PRECONDITIONER = ',PRECON(KIND_PRECOND)
         IF ( KIND_PRECOND.EQ.3 ) THEN
            WRITE(6,*) '  * Relaxation Parameter= ',RPRELIST(ITRY)
         END IF
         IF ( KIND_PRECOND.EQ.4 ) THEN
            WRITE(6,*) '  * Break Down Threshold= ',RPRELIST(ITRY)
         END IF
         IF ( KIND_PRECOND.EQ.6 ) THEN
            WRITE(6,7015) IPRELIST(ITRY),RPRELIST(ITRY)
         END IF
         WRITE(6,*) ' ++++++++++++++++++++++++++++++++++++++++++++++++
     $++++++'
         END IF
      END IF
*
*-----------------------------------------------------------------------
*         SIZE FOR ALLOCATE
*-----------------------------------------------------------------------
      CALL OPENATI_PRESET2(IATPARAM,SOLVER,KIND_PRECOND,IFILL,TAU,
     $                     IPOLICY,RPOLICY,N,NNZ,1,MM,ISIZE,IADDSIZE,
     $                     IDBG,INFO)
*
      WRITE(IO,*)
      WRITE(IO,*) '   MAX. SUBSPACE SIZE   =',MM
*
      IF (INFO.NE.0) THEN
         WRITE(IO,*) '   <<<< ERROR CODE >>>> = ',INFO
         GOTO 9900
      END IF
*
*-------------------------------MEMORY ALLOC CONTROL
      IF (ALLOCATED(WK1)) DEALLOCATE(WK1)
       ALLOCATE(WK1(ISIZE),STAT=IERR_ALLOC)
      IF (IERR_ALLOC .NE.0 ) THEN
         INFO=-500
         WRITE(IO,*) '   <<<< ERROR CODE >>>> = ',INFO
         GOTO 9900
      END IF
*
*-----------------------------------------------------------------------
*         PARAMETER SETTING
*-----------------------------------------------------------------------
*----------------SSOR's omega
      IF (KIND_PRECOND .EQ.3) THEN
         IF (RPRELIST(ITRY).LT.1.0D0 .OR. RPRELIST(ITRY).GE.2.0D0) THEN
            RATPARAM(25)=1.0D0
         ELSE
            RATPARAM(25)=RPRELIST(ITRY)
         ENDIF
*----------------ILU(0)_diag,ILU(0),ILUT(p,TAU)'s threshold of breakdown
      ELSE IF (KIND_PRECOND .GE.4) THEN
         RATPARAM(25)=RPRELIST(ITRY)
      END IF
      IF (KIND_PRECOND .EQ.6) THEN
         IATPARAM(26)=IPRELIST(ITRY)
      END IF
*---------------------Preconditioner create or re-use
      IPALLOC=0
      IF (ITRY.EQ.1) THEN
         IATPARAM(24) = 1
      ELSE
         IF (KIND_PRECOND.NE.PS_PAIR(1,ITRY-1) ) THEN
            IATPARAM(24) = 1
         ELSE
            IF (IFILL.NE.IPRELIST(ITRY-1) ) THEN
               IATPARAM(24) = 1
            ELSE
               IPALLOC=1
               IATPARAM(24) = 2
            END IF
         END IF
      END IF
*
      IF (IPALLOC.EQ.0) THEN
         IF (IDBG.EQ.1)
     $      WRITE(6,*) '[POLICY:PRECONDITIONER] ALLOCATE=',IADDSIZE
         IF (ALLOCATED(WK2)) DEALLOCATE(WK2)
         ALLOCATE(WK2(IADDSIZE),STAT=IERR_ALLOC)
         IF (IERR_ALLOC .NE.0 ) THEN
            INFO=-500
            WRITE(IO,*) '   <<<< ERROR CODE >>>> = ',INFO
            GOTO 9900
         END IF
      END IF
*
      USEMEM = 1.0D0*(ISIZE+IADDSIZE)
      USEMEM = USEMEM*8.0D-9
      WRITE(IO,7020) USEMEM
 7020 FORMAT(1H ,'   RUNTIME MEMORY USE   =',F10.3,' [GBYTE]')
*
*---------------------MAX. RESTARTS
      IATPARAM(22) = MIN(10000,N)
*---------------------REQUIREMENT ACCURACY
      RATPARAM(23) = RPOLICY(MP_RESID  )
*
*---------------------POLICY = STABLE
      IF (IPOLICY(LP_POLICY) .EQ. STABLE ) THEN
         IF (SOLVER.EQ.SOLVERLIST(1)) THEN
            MSIZE=MM
            IF (MSIZE.GE.30) THEN
               MSIZE=30
               IATPARAM(27) = MSIZE
            ELSE
               INFO=-400
               GOTO 9700
            ENDIF
            IATPARAM(28) = IATPARAM(27)
         END IF
         IATPARAM(4) = 0
         IATPARAM(9) = 0
         IATPARAM(10) = 12
         IATPARAM(12) = 2
      ELSE
*---------------------KRYLOV SUBSPACE MAX. SIZE FROM OPENATI_PRESET2
         MSIZE=MM
         IATPARAM(27) = MM
*---------------------INIT VALUE OF M
         IATPARAM(28) = 2 
         IATPARAM(28) = IRETRYMM
*
*---------------------KRYLOV SUBSPACE EXPAND AT ON
         IATPARAM(4)=1
*---------------------FASTEST MATRIX-VECTOR PRODUCT AT ON
         IF (IPOLICY(LP_POLICY) .EQ. MEMORY ) THEN
            IATPARAM(9) = 2
         ELSE
            IATPARAM(9) = 3
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
      WRITE(IO,7030) IATPARAM(4),IATPARAM(9)
 7030 FORMAT(1H ,'   KRYLOV SUBSPACE EXPAND AT = ',I2,
     $           '   ,MATVEC AT = ',I2)
*
      WRITE(IO,*) '   Initial Gram-Schmidt Strategy = ',
     $            GSPOL(IATPARAM(12))
      WRITE(IO,*)
*------------WORK POINTER
*
      IF      (SOLVER.EQ.'XABCLIB_GMRES') THEN
*
         LDWK=(MSIZE+2)*N+(MSIZE+1)*(MSIZE+1)+(N-1)/2+2
*
      ELSE IF (SOLVER.EQ.'XABCLIB_BICGSTAB') THEN
*
         LDWK=10*N
*
*     ELSE IF (SOLVER.EQ.'XABCLIB_BICGSTAB_ORG') THEN
*
*        NPRE=10*N
*
      END IF
*
      RES=-1.0D0
      LOOP=1
      IREST=0
*
      IATPARAM(17)=0
*
      IF (ITRY.NE.1) THEN
         X(1:N) = ANS(1:N)
      END IF
*
      ETIMEB  = OMP_GET_WTIME()
      ETIMEP  = ETIMEP  + (ETIMEB - ETIMEA)
      TIMELIM = TIMELIM - (ETIMEB - ETIMEA)
*
 3000 CONTINUE
*---------------------Maximum elapsed time for the solver
      IF (TIMEMAX.GT.0.0D0) THEN
         RATPARAM(22) = TIMELIM
      ELSE
         RATPARAM(22) = TIMEMAX
      END IF
      IF (IDBG.EQ.1) WRITE(6,*) '[POLICY:TIMELIMIT] =',TIMELIM,' [SEC]'
*
      IF      (SOLVER.EQ.'XABCLIB_GMRES') THEN
*
         CALL Xabclib_GMRES (N,NNZ,IRP,ICOL,VAL,B,X,WK2,
     $                       IADDSIZE,IATPARAM,RATPARAM,WK1,LDWK,INFO)
*
      ELSE IF (SOLVER.EQ.'XABCLIB_BICGSTAB') THEN
*
         CALL Xabclib_BICGSTAB (N,NNZ,IRP,ICOL,VAL,B,X,WK2,
     $                       IADDSIZE,IATPARAM,RATPARAM,WK1,LDWK,INFO)
*
*     ELSE IF (SOLVER.EQ.'XABCLIB_BICGSTAB_ORG') THEN
*
*        CALL Xabclib_BICGSTABV (N,NNZ,IRP,ICOL,VAL,B,X,WK1(NPRE),
*    $                       IADDSIZE,IATPARAM,RATPARAM,WK1,LDWK,INFO)
*
      END IF
*
*----------- ETIME0 : SUM OF PCOND, ETIME1 : SUM OF TOTAL XABCLIB SOLVE
      ETIME0= ETIME0 + RATPARAM(31)
      ETIME1= ETIME1 + RATPARAM(32)
      TIMELIM = TIMELIM - RATPARAM(32)
      IREST=IREST+IATPARAM(23)
*
      IATPARAM(17)=IATPARAM(17) + IATPARAM(31)
*
      RATPARAM(14)=RATPARAM(29)
*
         IRETRYMM=IATPARAM(29)
      IF (INFO .EQ. 1000 .AND. RATPARAM(29) .LT. GOODANS) THEN
         GOODANS = RATPARAM(29)
         ANS(1:N)=X(1:N)
         IRETRY=ITRY
         IRETRYMM=IATPARAM(29)
         PS_PAIR(1,NTRY+1)=PS_PAIR(1,ITRY)
         PS_PAIR(2,NTRY+1)=PS_PAIR(2,ITRY)
         IPRELIST(NTRY+1) =IPRELIST(ITRY)
         RPRELIST(NTRY+1) =RPRELIST(ITRY)
      END IF
*
      IF (INFO.NE.0 ) THEN
         IF (TIMEMAX .GT. 0  .AND.  TIMELIM .LE. 0.0D0) THEN
            TOVER=1
         END IF
         GOTO 9700
      END IF
      IF (TIMEMAX .GT. 0  .AND.  TIMELIM .LE. 0.0D0) THEN
         TOVER=1
         GOTO 9700
      END IF
*
*-------------<<POLICY 2>>: ACCURACY
      IF (IPOLICY(LP_POLICY   ).EQ. ACCURACY) THEN
*-------------------RESIDUAL CHECK
         CALL OpenATI_DURMV_11(N,NNZ,IRP,ICOL,VAL,X,WK1)
         S=0.0D0
         DO I=1,N
            S=S+(WK1(I)-B(I))*(WK1(I)-B(I))
         ENDDO
         SS=DSQRT(S)
         RES=SS/BNORM
         RATPARAM(14)=RES
         IF (IDBG.EQ.1)
     $      WRITE(6,*) '[POLICY:ACCURACY] RESIDUAL=',RES
*
         IF (RES .LT. RPOLICY(MP_RESID  ) ) THEN
            GOTO 9000
         ELSE
*----------------------------------DGKS
            IATPARAM(12) = 1
            IATPARAM(25) = 1
            RATPARAM(23) = RATPARAM(23)*0.1D0
            LOOP=LOOP+1
            GOTO 3000
         END IF
      ELSE
         RES=RATPARAM(29)
         GOTO 9000
      END IF
*
*--------------------------------------SUCCESSFULY ENDING
 9000 CONTINUE
      WRITE(IO,*) ' ====== OPENATI_LINEARSOLVE SUCCESSFULY ENDED ======'
      GOTO 9800
*--------------------------------------ERROR AFTER ALLOCATE
 9700 CONTINUE
      WRITE(IO,*) ' ****** OPENATI_LINEARSOLVE ERROR DETECTED *********'
      WRITE(IO,*) '   <<<< ERROR CODE >>>> = ',INFO
*
 9800 CONTINUE
      IATPARAM(15)=LOOP
      IATPARAM(16)=IREST
      IATPARAM(18)=KIND_PRECOND
      IATPARAM(19)=IATPARAM(26)
      IATPARAM(20)=KIND_SOLVER
*
      WRITE(IO,*)
      WRITE(IO,*) ' [OPENATI_LINEARSOLVE RESULT]'
      WRITE(IO,*) '   MATRIX DATA : N= ',N,' NNZ= ',NNZ
      WRITE(IO,*) '   FASTEST MATVEC NO. = ',IATPARAM(10)
      WRITE(IO,*) '   2-Norm of RHS =',RATPARAM(28)
      IF      (SOLVER.EQ.'XABCLIB_GMRES') THEN
      WRITE(IO,*) '   FINAL KRYLOV SUBSPACE SIZE =',IATPARAM(29)
      WRITE(IO,*) '   FINAL Gram-Schmidt Strategy = ',
     $                     GSPOL(IATPARAM(12))
      WRITE(IO,*) '   NUMBER OF RETRYED GMRES =',LOOP
      WRITE(IO,*) '   TOTAL RESTARTS of GMRES =',IREST
      ELSE IF (SOLVER.EQ.'XABCLIB_BICGSTAB') THEN
      WRITE(IO,*) '   NUMBER OF RETRYED BICGSTAB =',LOOP
      WRITE(IO,*) '   TOTAL RESTARTS of BICGSTAB =',IREST
      END IF
      WRITE(IO,*) '   NUMBER OF MAT-VEC       =',IATPARAM(31)
      WRITE(IO,*) '   RESIDUAL NORM           =',RATPARAM(14)
*
      IF (INFO.EQ.0) THEN
         WRITE(IO,*) '   <<<< THIS LINEAR SYSTEM IS SOLVED >>>> '
         GOTO 200
      END IF
      IF (TOVER.EQ.1) THEN
         WRITE(IO,*) '   <<<< UPPER LIMIT OF TIME >>>> '
         GOTO 200
      END IF
*
*[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[
  100 CONTINUE
*[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[
  200 CONTINUE
*
*
      RATPARAM(15)=ETIMES + ETIMEP
      RATPARAM(16)=ETIME0
      RATPARAM(17)=ETIME1 - ETIME0
      RATPARAM(18)=ETIME1 + RATPARAM(15)
      RATPARAM(19)=RATPARAM(25)
*
      WRITE(IO,*) '   SET-UP   TIME           =',RATPARAM(15),' [SEC]'
      WRITE(IO,*) '   PRECOND. TIME           =',RATPARAM(16),' [SEC]'
      WRITE(IO,*) '   SOLVER   TIME           =',RATPARAM(17),' [SEC]'
      WRITE(IO,*)
      WRITE(IO,*) '   TOTAL    TIME           =',RATPARAM(18),' [SEC]'
*
      IF (ALLOCATED(ANS)) THEN
         DEALLOCATE(ANS,STAT=IERR_DEALLOC)
         IF (IERR_DEALLOC.NE.0) THEN
         WRITE(6,*) '[Error] IERR_DEALLOC=',IERR_DEALLOC
         END IF
      END IF
      DEALLOCATE(WK1,STAT=IERR_DEALLOC)
      IF (IERR_DEALLOC.NE.0) THEN
         WRITE(6,*) '[Error] IERR_DEALLOC=',IERR_DEALLOC
      END IF
      DEALLOCATE(WK2,STAT=IERR_DEALLOC)
      IF (IERR_DEALLOC.NE.0) THEN
         WRITE(6,*) '[Error] IERR_DEALLOC=',IERR_DEALLOC
      END IF
*--------------------------------------PARAMETER ERROR
 9900 CONTINUE
      RETURN
      END
