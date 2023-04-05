      SUBROUTINE OPENATI_PRESET2(IATPARAM,SOLVER,KIND_PRECOND,IFILL,TAU,
     $                           IPOLICY,RPOLICY,N,NNZ,NEV,
     $                           MM,ISIZE,IADDSIZE,IDBG,INFO)
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
      INTEGER           IATPARAM(50)
*
      PARAMETER (NUM_POL=20)
*
      INTEGER  TIME,MEMORY,ACCURACY
      INTEGER  NO, JACOBI, SSOR, ILU0D, ILU0, ILUT, AUTO
      PARAMETER (LP_SMPS   = 1,
     $           LP_POLICY = 2, TIME=1, MEMORY=2, ACCURACY=3,
     $           LP_PRECON = 3, NO=1, JACOBI=2, SSOR=3, ILU0D=4,
     $                          ILU0=5, ILUT=6, AUTO=10,
     $           LP_SOLVER = 4,
     $           MP_MAXMEM = 1,
     $           MP_RESID  = 2,
     $           MP_MAXTIME= 3)
*
      CHARACTER*2   DIM
      CHARACTER*15  KEY
      CHARACTER*32  SOLVER
      CHARACTER*80  C80,C80W
*     ..
*     ..Array Arguments
      INTEGER IPOLICY(NUM_POL)
      DOUBLE PRECISION RPOLICY(NUM_POL)
*
      DOUBLE PRECISION ,ALLOCATABLE :: TMP(:)
*
*  =====================================================================
*  Purpose
*  =======
*
*  OpenATI_PRESET2 sets up MAX. KRYLOV SUBSPACE SIZE 'MM' 
*  for OPENATI_LINEARSOLVE.
*
*  Arguments
*  =========
*  SOLVER  (input) CHARACTER*15
*          target of calculation.
*          The content is XABCLIB_GMRES or XABCLIB_LANCZOS
*
*  IPOLICY (output) INTEGER array, dimension (NUM_POL)
*  (NUM    Integer parameter set for policy.
*    _POL)  IPOLICY(1)      : The number of threads.
*           IPOLICY(2)      : Set User policy kinds.
*                      =  1 : Execute time.
*                      =  2 : Run time's memory size.
*                      =  3 : Accuracy.
*           IPOLICY(3)      : Set preconditioner kinds.
*                      =  0 : None.
*                      =  1 : Use Jacobi.
*                      =  2 : Use SSOR.
*                      =  3 : Use ILU(0).
*           IPOLICY(4)      : Value of solver type.
*           IPOLICY(5:20)   : Don't use.
*
*
*  RPOLICY (output) DOUBLE PRECISION array, dimension (NUM_POL)
*  (NUM    Double precision parameter set for policy.
*    _POL)  RPOLICY(1)      : The threshold value of memory size.
*           RPOLICY(2)      : The threshold value of convergence.
*           RPOLICY(3)      : The threshold value of execute time.
*           RPOLICY(4:20)   : Don't use.
*
*  N       (input) INTEGER
*          The order of the matrix.
*
*  NNZ     (input) INTEGER
*          Number of nonzero elements in the matrix.
*
*  NEV     (input) INTEGER
*          request number of eigenvalue
*
*  MM      (output) INTEGER
*          restart revolution
*
*  ISIZE    (output) INTEGER
*          USING WORK AREA's SIZE
*
*  IDBG    (input) INTEGER
*          Debug print status.
*
*  INFO    (output) INTEGER
*          Return code.
*          =    0:  Successful exit
*          = -400:  cannot allocate minimum require working set
*
*  =====================================================================
*
      INTEGER           THREAD_NUM, THREAD_MAX
      CHARACTER*16      THREAD_NUM_CHA
      CHARACTER*20      FNAME
*
      INTEGER           OMP_GET_THREAD_NUM, OMP_GET_MAX_THREADS
      EXTERNAL          OMP_GET_THREAD_NUM, OMP_GET_MAX_THREADS

      INFO=0


*-----------------------------------------------------------------------
*         set max memory size 
*-----------------------------------------------------------------------
      XMAXMEM=RPOLICY(MP_MAXMEM)
      IF (XMAXMEM.GT.0 .OR. IATPARAM(14).EQ.0 ) THEN
         XMEMINFO=XMAXMEM
         GOTO 2100
      END IF
*-----------------------------------------------------------------------
*         get free memory by Linux '/proc/meminfo' 
*-----------------------------------------------------------------------
      THREAD_MAX= OMP_GET_MAX_THREADS()
      THREAD_NUM= OMP_GET_THREAD_NUM()
      WRITE(THREAD_NUM_CHA,*) THREAD_NUM
      FNAME='OpenATI_MemFree.'//ADJUSTL(THREAD_NUM_CHA)
#ifdef HF_COMPILER
      CALL HF_SH('cat /proc/meminfo |grep -i memfree
     $             > '//FNAME,idum)
#else
      CALL SYSTEM('cat /proc/meminfo |grep -i memfree
     $             >'//FNAME)
#endif
*
      LU=50+THREAD_NUM
      OPEN(LU,FILE=FNAME)
      READ(LU,1000,END=999) C80
 1000 FORMAT(A80)
      IEQ=INDEX(C80,':')
      C80W=C80(IEQ+1:80)
      KEY=ADJUSTL(C80W)
      READ(KEY,*) MEM,DIM
      IF (DIM.EQ.'kB' .OR. DIM.EQ.'KB' .OR. DIM.EQ.'kb' ) THEN
         XMEMINFO=MEM/1.0D+6
      ELSE
         GOTO 999
      END IF
      GOTO 2000
*
  999 CONTINUE
      WRITE(6,*) '[OpenATI Warnig] not found command "/proc/meminfo"
     $ on your system. OpenATI search free memory.'
      XMEMINFO=-1.0D0
*
 2000 CONTINUE
      CLOSE(LU)
#ifdef HF_COMPILER
      CALL HF_SH('rm '//FNAME,idum)
#else
      CALL SYSTEM('rm '//FNAME)
#endif
*
      IF (IDBG.EQ.1) WRITE(6,*) 'XMEMINFO=',XMEMINFO,' [GB]'
*
      IF (XMEMINFO.GE.0) THEN
            XMAXMEM = DMIN1(XMEMINFO,16.0D0)
      END IF
*
*-----------------------------------------------------------------------
 2100 CONTINUE
      IF (IDBG.EQ.1) WRITE(6,*) 'XMEMINFO=',XMEMINFO,' [GB]'
      IF (IDBG.EQ.1) WRITE(6,*) 'XMAXMEM =',XMAXMEM,' [GB]'
      RPOLICY(MP_MAXMEM) = XMAXMEM
*----------------------------------------------------------------------
      IF      (KIND_PRECOND.LE.4) THEN
         IADDSIZE=N
         ADDSIZE=8.0D0*1.0D-9*IADDSIZE
      ELSE IF      (KIND_PRECOND.EQ.5) THEN
         IADDSIZE=3*NNZ/2+2*N+50
         ADDSIZE=8.0D0*1.0D-9*IADDSIZE
      ELSE IF (KIND_PRECOND.EQ.6) THEN
         IADDSIZE=3*(2.0D0*IFILL+1)*N/2+3*N+100
         ADDSIZE=8.0D0*1.0D-9*IADDSIZE
      END IF
      IF (IDBG.EQ.1) WRITE(6,*) 'PRECOND SIZE =',ADDSIZE,' [GB]'
*----------------------------------------------------------------------
      IF      (SOLVER .EQ. 'XABCLIB_GMRES') THEN
         M=2
      ELSE IF (SOLVER .EQ. 'XABCLIB_BICGSTAB') THEN
         M=0
*     ELSE IF (SOLVER .EQ. 'XABCLIB_BICGSTAB_ORG') THEN
*        M=0
      ELSE
         WRITE(6,*) '[Error] SOLVER NOT SELECT  ',SOLVER
      END IF
      CALL OPENATI_GET_WORKSIZE(SOLVER,M,N,NEV,SIZEA)
      SIZEA =SIZEA+ADDSIZE
*
      IF      (SOLVER .EQ. 'XABCLIB_GMRES') THEN
         M=100
      ELSE IF (SOLVER .EQ. 'XABCLIB_BICGSTAB') THEN
         M=0
*     ELSE IF (SOLVER .EQ. 'XABCLIB_BICGSTAB_ORG') THEN
*        M=0
      ELSE
         WRITE(6,*) '[Error] SOLVER NOT SELECT  ',SOLVER
      END IF
      CALL OPENATI_GET_WORKSIZE(SOLVER,M,N,NEV,SIZEB)
      SIZEB =SIZEB+ADDSIZE
*
      IF      (SOLVER .EQ.'XABCLIB_GMRES') THEN
         SIZEC=SIZEB + 1.0D-9*9*NNZ
      ELSE IF (SOLVER .EQ. 'XABCLIB_BICGSTAB') THEN
         SIZEC=SIZEB + 1.0D-9*9*NNZ
*     ELSE IF (SOLVER .EQ. 'XABCLIB_BICGSTAB_ORG') THEN
*        SIZEC=SIZEB + 1.0D-9*9*NNZ
      ELSE
         WRITE(6,*) '[Error] SOLVER NOT SELECT  ',SOLVER
      END IF
*
      IF (IDBG.EQ.1) WRITE(6,*) 'SIZEA,SIZEB,SIZEC,=',SIZEA,SIZEB,SIZEC,
     $                          ' [GB]'
*
      IF (XMAXMEM.GT.0  .AND. XMEMINFO.GT.0 ) THEN
         IF (SIZEA.GT.XMAXMEM) THEN
            INFO=-400
            WRITE(6,*) '[Not enough space]'
         END IF
         IF (SIZEB.LT.XMAXMEM) THEN
            MM=100
         ELSE
            CALL OPENATI_GET_MAXMSIZE(SOLVER,XMAXMEM-ADDSIZE,N,NEV,MM)
            GOTO 8000
         END IF
         IF (IPOLICY(LP_POLICY).EQ.MEMORY) THEN
            GOTO 8000
         END IF
         IF (SIZEC.LT.XMAXMEM) THEN
            CALL OPENATI_GET_MAXMSIZE(SOLVER,XMAXMEM-ADDSIZE,N,NEV,MM)
            IF (IDBG.EQ.1) WRITE(6,*) 'MAX MSIZE=',MM
         ELSE
            GOTO 8000
         END IF
      ELSE
*-----------------------------------------------------------------------
*       UNKNOWN XMEMINFO OR UNKNOWN XMAXMEM     
*-----------------------------------------------------------------------
         IF (IDBG.EQ.1) WRITE(6,*) '!!! MEMORY SEARCH SEQUENCE START'
         IF (IPOLICY(LP_POLICY).EQ.MEMORY) THEN 
*------------------------------------------TRY TO M=100
            M=100
            CALL OPENATI_GET_ALLOCSIZE(SOLVER,M,N,NEV,JSIZE)
            JSIZE=JSIZE+ADDSIZE
            ALLOCATE(TMP(JSIZE),STAT=IERR_ALLOC)
            IF (IERR_ALLOC .EQ.0 ) THEN
               DEALLOCATE(TMP,STAT=IERR_DALLOC)
               MM=100
               GOTO 8000
            ELSE
*------------------------could't allocate M=100
               MSEARCHMAX=90
            END IF
         ELSE
            IF (XMAXMEM.GT.0) THEN
              CALL OPENATI_GET_MAXMSIZE(SOLVER,XMAXMEM-ADDSIZE,
     $                                  N,NEV,MSEARCHMAX)
            ELSE
              XTMP=16.0D0
              CALL OPENATI_GET_MAXMSIZE(SOLVER,XTMP-ADDSIZE,
     $                                  N,NEV,MSEARCHMAX)
            END IF
         END IF
         MSEARCHMIN=2
*
         ITR=0
 3100    CONTINUE
         ITR=ITR+1
         KM4 = (MSEARCHMAX-MSEARCHMIN)/4
         MF=MSEARCHMIN
         DO I=1,5
            M=MF+KM4*(I-1)
            IF (M.GT.N) M=N
            CALL OPENATI_GET_ALLOCSIZE(SOLVER,M,N,NEV,JSIZE)
            JSIZE=JSIZE+ADDSIZE
            ALLOCATE(TMP(JSIZE),STAT=IERR_ALLOC)
            IF (IERR_ALLOC .EQ.0 ) THEN
*--------------------------------------ALLOC SUCCESS
               DEALLOCATE(TMP,STAT=IERR_DALLOC)
               IF (I.EQ.5) THEN
                  MM=M
                  ISIZE=JSIZE
                  GO TO 9000
               ELSE
                  MSEARCHMIN=M
               END IF
            ELSE
               IF (I.EQ.1) THEN
                  INFO=-400
                  WRITE(6,*) '[Error : cannot allocate minimum require
     $working set]'
               ELSE
                  MSEARCHMAX=M
                  GO TO 3000
               END IF
            END IF
         ENDDO
 3000    CONTINUE
         IF (ITR.LE.1) THEN
            GOTO 3100
         ELSE
            MM=MSEARCHMIN
         END IF
*
      END IF
*-----------------------------------------------------------------------
*
*-----------------------------------------------------------------------
*          SIZE FOR ALLOCATE KRYLOV SUBSPACE
*-----------------------------------------------------------------------
 8000 CONTINUE
*-----------------------------------------------------------------------
      CALL OPENATI_GET_ALLOCSIZE(SOLVER,MM,N,NEV,ISIZE)
      IF (IDBG.EQ.1) WRITE(6,*) 'POLICY MSIZE=',MM,' ISIZE=',ISIZE,
     $' IADDSIZE=',IADDSIZE
*
 9000 CONTINUE
      RETURN
      END
      SUBROUTINE OPENATI_PRESET(IATPARAM,SOLVER,IPOLICY,RPOLICY,N,NNZ,
     $                          NEV,MM,ISIZE,IDBG,INFO)
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
      INTEGER           IATPARAM(50)
*
      PARAMETER (NUM_POL=20)
*
      INTEGER  TIME,MEMORY,ACCURACY
      INTEGER  NO, JACOBI, SSOR, ILU0
      PARAMETER (LP_SMPS   = 1,
     $           LP_POLICY = 2, TIME=1, MEMORY=2, ACCURACY=3,
     $           LP_PRECON = 3, NO=1, JACOBI=2, SSOR=3, ILU0=4,
     $           LP_SOLVER = 4,
     $           MP_MAXMEM = 1,
     $           MP_RESID  = 2,
     $           MP_MAXTIME= 3)
*
      CHARACTER*2   DIM
      CHARACTER*15  KEY
      CHARACTER*32  SOLVER
      CHARACTER*80  C80,C80W
*     ..
*     ..Array Arguments
      INTEGER IPOLICY(NUM_POL)
      DOUBLE PRECISION RPOLICY(NUM_POL)
*
      DOUBLE PRECISION ,ALLOCATABLE :: TMP(:)
*
*  =====================================================================
*  Purpose
*  =======
*
*  OpenATI_PRESET sets up MAX. KRYLOV SUBSPACE SIZE 'MM'.
*
*  Arguments
*  =========
*  SOLVER  (input) CHARACTER*15
*          target of calculation.
*          The content is XABCLIB_GMRES or XABCLIB_LANCZOS
*
*  IPOLICY (output) INTEGER array, dimension (NUM_POL)
*  (NUM    Integer parameter set for policy.
*    _POL)  IPOLICY(1)      : The number of threads.
*           IPOLICY(2)      : Set User policy kinds.
*                      =  1 : Execute time.
*                      =  2 : Run time's memory size.
*                      =  3 : Accuracy.
*           IPOLICY(3)      : Set preconditioner kinds.
*                      =  0 : None.
*                      =  1 : Use Jacobi.
*                      =  2 : Use SSOR.
*                      =  3 : Use ILU(0).
*           IPOLICY(4)      : Value of solver type.
*           IPOLICY(5:20)   : Don't use.
*
*
*  RPOLICY (output) DOUBLE PRECISION array, dimension (NUM_POL)
*  (NUM    Double precision parameter set for policy.
*    _POL)  RPOLICY(1)      : The threshold value of memory size.
*           RPOLICY(2)      : The threshold value of convergence.
*           RPOLICY(3)      : The threshold value of execute time.
*           RPOLICY(4:20)   : Don't use.
*
*  N       (input) INTEGER
*          The order of the matrix.
*
*  NNZ     (input) INTEGER
*          Number of nonzero elements in the matrix.
*
*  NEV     (input) INTEGER
*          request number of eigenvalue
*
*  MM      (output) INTEGER
*          restart revolution
*
*  ISIZE    (output) INTEGER
*          USING WORK AREA's SIZE
*
*  IDBG    (input) INTEGER
*          Debug print status.
*
*  INFO    (output) INTEGER
*          Return code.
*          =    0:  Successful exit
*          = -400:  cannot allocate minimum require working set
*
*  =====================================================================
*
      INTEGER           THREAD_NUM, THREAD_MAX
      CHARACTER*16      THREAD_NUM_CHA
      CHARACTER*20      FNAME
*
      INTEGER           OMP_GET_THREAD_NUM, OMP_GET_MAX_THREADS
      EXTERNAL          OMP_GET_THREAD_NUM, OMP_GET_MAX_THREADS

      INFO=0


*-----------------------------------------------------------------------
*         set max memory size 
*-----------------------------------------------------------------------
      XMAXMEM=RPOLICY(MP_MAXMEM)
      IF (XMAXMEM.GT.0 .OR. IATPARAM(14).EQ.0 ) THEN
         XMEMINFO=XMAXMEM
         GOTO 2100
      END IF
*-----------------------------------------------------------------------
*         get free memory by Linux '/proc/meminfo' 
*-----------------------------------------------------------------------
      THREAD_MAX= OMP_GET_MAX_THREADS()
      THREAD_NUM= OMP_GET_THREAD_NUM()
      WRITE(THREAD_NUM_CHA,*) THREAD_NUM
      FNAME='OpenATI_MemFree.'//ADJUSTL(THREAD_NUM_CHA)
#ifdef HF_COMPILER
      CALL HF_SH('cat /proc/meminfo |grep -i memfree
     $             > '//FNAME,idum)
#else
      CALL SYSTEM('cat /proc/meminfo |grep -i memfree
     $             >'//FNAME)
#endif
*
      LU=50+THREAD_NUM
      OPEN(LU,FILE=FNAME)
      READ(LU,1000,END=999) C80
 1000 FORMAT(A80)
      IEQ=INDEX(C80,':')
      C80W=C80(IEQ+1:80)
      KEY=ADJUSTL(C80W)
      READ(KEY,*) MEM,DIM
      IF (DIM.EQ.'kB' .OR. DIM.EQ.'KB' .OR. DIM.EQ.'kb' ) THEN
         XMEMINFO=MEM/1.0D+6
      ELSE
         GOTO 999
      END IF
      GOTO 2000
*
  999 CONTINUE
      WRITE(6,*) '[OpenATI Warnig] not found command "/proc/meminfo"
     $ on your system. OpenATI search free memory.'
      XMEMINFO=-1.0D0
*
 2000 CONTINUE
      CLOSE(LU)
#ifdef HF_COMPILER
      CALL HF_SH('rm '//FNAME,idum)
#else
      CALL SYSTEM('rm '//FNAME)
#endif
*
      IF (IDBG.EQ.1) WRITE(6,*) 'XMEMINFO=',XMEMINFO,' [GB]'
*
      IF (XMEMINFO.GE.0) THEN
            XMAXMEM = DMIN1(XMEMINFO,16.0D0)
      END IF
*
*-----------------------------------------------------------------------
 2100 CONTINUE
      IF (IDBG.EQ.1) WRITE(6,*) 'XMAXMEM =',XMAXMEM,' [GB]'
      RPOLICY(MP_MAXMEM) = XMAXMEM
*----------------------------------------------------------------------
      IF      (SOLVER .EQ. 'XABCLIB_GMRES') THEN
         M=2
      ELSE IF (SOLVER .EQ. 'XABCLIB_BICGSTAB') THEN
         M=0
*     ELSE IF (SOLVER .EQ. 'XABCLIB_BICGSTAB_ORG') THEN
*        M=0
      ELSE IF (SOLVER .EQ. 'XABCLIB_LANCZOS') THEN
         M=NEV
      ELSE IF (SOLVER .EQ. 'XABCLIB_ARNOLDI') THEN
         M=NEV
      ELSE
         WRITE(6,*) '[Error] SOLVER NOT SELECT  ',SOLVER
      END IF
      CALL OPENATI_GET_WORKSIZE(SOLVER,M,N,NEV,SIZEA)
*
      IF      (SOLVER .EQ. 'XABCLIB_GMRES') THEN
         M=100
      ELSE IF (SOLVER .EQ. 'XABCLIB_BICGSTAB') THEN
         M=0
*     ELSE IF (SOLVER .EQ. 'XABCLIB_BICGSTAB_ORG') THEN
*        M=0
      ELSE IF (SOLVER .EQ. 'XABCLIB_LANCZOS') THEN
         M=MAX(2*NEV,100)
      ELSE IF (SOLVER .EQ. 'XABCLIB_ARNOLDI') THEN
         M=MAX(2*NEV,100)
      ELSE
         WRITE(6,*) '[Error] SOLVER NOT SELECT  ',SOLVER
      END IF
      CALL OPENATI_GET_WORKSIZE(SOLVER,M,N,NEV,SIZEB)
*
      IF      (SOLVER .EQ.'XABCLIB_GMRES') THEN
         SIZEC=SIZEB + 1.0D-9*9*NNZ
      ELSE IF (SOLVER .EQ. 'XABCLIB_BICGSTAB') THEN
         SIZEC=SIZEB + 1.0D-9*9*NNZ
*     ELSE IF (SOLVER .EQ. 'XABCLIB_BICGSTAB_ORG') THEN
*        SIZEC=SIZEB + 1.0D-9*9*NNZ
      ELSE IF (SOLVER .EQ. 'XABCLIB_LANCZOS') THEN
         SIZEC=SIZEB
      ELSE IF (SOLVER .EQ. 'XABCLIB_ARNOLDI') THEN
         SIZEC=SIZEB + 1.0D-9*9*NNZ
      ELSE
         WRITE(6,*) '[Error] SOLVER NOT SELECT  ',SOLVER
      END IF
*
      IF (IDBG.EQ.1) WRITE(6,*) 'SIZEA,SIZEB,SIZEC,=',SIZEA,SIZEB,SIZEC,
     $                          ' [GB]'
*
      IF (XMAXMEM.GT.0  .AND. XMEMINFO.GT.0 ) THEN
         IF (SIZEA.GT.XMAXMEM) THEN
            INFO=-400
            WRITE(6,*) '[Not enough space]'
         END IF
         IF (SIZEB.LT.XMAXMEM) THEN
            MM=100
         ELSE
            CALL OPENATI_GET_MAXMSIZE(SOLVER,XMAXMEM,N,NEV,MM)
            GOTO 8000
         END IF
         IF (IPOLICY(LP_POLICY).EQ.MEMORY) THEN
            GOTO 8000
         END IF
         IF (SIZEC.LT.XMAXMEM) THEN
            CALL OPENATI_GET_MAXMSIZE(SOLVER,XMAXMEM,N,NEV,MM)
            IF (IDBG.EQ.1) WRITE(6,*) 'MAX MSIZE=',MM
         ELSE
            GOTO 8000
         END IF
      ELSE
*-----------------------------------------------------------------------
*       UNKNOWN XMEMINFO OR UNKNOWN XMAXMEM     
*-----------------------------------------------------------------------
         IF (IDBG.EQ.1) WRITE(6,*) '!!! MEMORY SEARCH SEQUENCE START'
         IF (IPOLICY(LP_POLICY).EQ.MEMORY) THEN 
*------------------------------------------TRY TO M=100
            M=100
            CALL OPENATI_GET_ALLOCSIZE(SOLVER,M,N,NEV,JSIZE)
            ALLOCATE(TMP(JSIZE),STAT=IERR_ALLOC)
            IF (IERR_ALLOC .EQ.0 ) THEN
               DEALLOCATE(TMP,STAT=IERR_DALLOC)
               MM=100
               GOTO 8000
            ELSE
*------------------------could't allocate M=100
               MSEARCHMAX=90
            END IF
         ELSE
            IF (XMAXMEM.GT.0) THEN
              CALL OPENATI_GET_MAXMSIZE(SOLVER,XMAXMEM,N,NEV,MSEARCHMAX)
            ELSE
              XTMP=16.0D0
              CALL OPENATI_GET_MAXMSIZE(SOLVER,XTMP,   N,NEV,MSEARCHMAX)
            END IF
         END IF
         MSEARCHMIN=NEV
*
         ITR=0
 3100    CONTINUE
         ITR=ITR+1
         KM4 = (MSEARCHMAX-MSEARCHMIN)/4
         MF=MSEARCHMIN
         DO I=1,5
            M=MF+KM4*(I-1)
            IF (M.GT.N) M=N
            CALL OPENATI_GET_ALLOCSIZE(SOLVER,M,N,NEV,JSIZE)
            ALLOCATE(TMP(JSIZE),STAT=IERR_ALLOC)
            IF (IERR_ALLOC .EQ.0 ) THEN
*--------------------------------------ALLOC SUCCESS
               DEALLOCATE(TMP,STAT=IERR_DALLOC)
               IF (I.EQ.5) THEN
                  MM=M
                  ISIZE=JSIZE
                  GO TO 9000
               ELSE
                  MSEARCHMIN=M
               END IF
            ELSE
               IF (I.EQ.1) THEN
                  INFO=-400
                  WRITE(6,*) '[Error : cannot allocate minimum require
     $working set]'
               ELSE
                  MSEARCHMAX=M
                  GO TO 3000
               END IF
            END IF
         ENDDO
 3000    CONTINUE
         IF (ITR.LE.1) THEN
            GOTO 3100
         ELSE
            MM=MSEARCHMIN
         END IF
*
      END IF
*-----------------------------------------------------------------------
*
*-----------------------------------------------------------------------
*          SIZE FOR ALLOCATE KRYLOV SUBSPACE
*-----------------------------------------------------------------------
 8000 CONTINUE
*-----------------------------------------------------------------------
      CALL OPENATI_GET_ALLOCSIZE(SOLVER,MM,N,NEV,ISIZE)
      IF (IDBG.EQ.1) WRITE(6,*) 'POLICY MSIZE=',MM,' ISIZE=',ISIZE
*
 9000 CONTINUE
      RETURN
      END
*
*
      SUBROUTINE OPENATI_GET_MAXMSIZE(SOLVER,XMAXMEM,N,NEV,MM)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*32  SOLVER
*
*
*  =====================================================================
*  Purpose
*  =======
*
*  Calculation for MAXMSIZE
*
*  Arguments
*  =========
*  SOLVER  (input) CHARACTER*32
*          target of calculation.
*          The content is XABCLIB_GMRES or XABCLIB_LANCZOS
*
*  XMAXMEM (input) DOUBLE PRECISION
*          The threshold value of memory size.
*
*  N       (input) INTEGER
*          The order of the matrix.
*
*  NEV     (input) INTEGER
*          request number of eigenvalue
*
*  MM      (output) INTEGER
*          restart revolution
*
* ==================================================================
*
      DN=N*1.0D0
      IF      (SOLVER.EQ.'XABCLIB_GMRES') THEN
         ALF=XMAXMEM/8*1.0D+9  - 3.0D0*N
         XM=0.5D0*(-DN+ SQRT(DN*DN+4.0D0*ALF))
      ELSE IF (SOLVER.EQ.'XABCLIB_BICGSTAB') THEN
         XM=0.0D0
*     ELSE IF (SOLVER.EQ.'XABCLIB_BICGSTAB_ORG') THEN
*        XM=0.0D0
      ELSE IF (SOLVER.EQ.'XABCLIB_LANCZOS') THEN
         ALF=XMAXMEM/8*1.0D+9  - 1.0D0*N - 5.0D0*NEV
         XM=0.25D0*(-DN+ SQRT(DN*DN+8.0D0*ALF))
      ELSE IF (SOLVER.EQ.'XABCLIB_ARNOLDI') THEN
         ALF=XMAXMEM/8*1.0D+9  - 5.0D0*N - 6.0D0*NEV
         XM=0.1D0*(-DN+ SQRT(DN*DN+20.0D0*ALF))
      END IF
*
      MM=INT(XM/10)*10
      MM=MIN(MM,N)

      RETURN
      END

      SUBROUTINE OPENATI_GET_ALLOCSIZE(SOLVER,M,N,NEV,ISIZE)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*32  SOLVER
*
*
*  =====================================================================
*  Purpose
*  =======
*
*  Calculation for Allocate size in double precision.
*
*  Arguments
*  =========
*  SOLVER  (input) CHARACTER*15
*          target of calculation.
*          The content is XABCLIB_GMRES or XABCLIB_LANCZOS
*
*  M       (output) INTEGER
*          restart revolution
*
*  N       (input) INTEGER
*          The order of the matrix.
*
*  NEV     (input) INTEGER
*          request number of eigenvalue
*
*  ISIZE    (output) INTEGER
*          USING WORK AREA's SIZE
*
* ==================================================================
*
*
      IF      (SOLVER.EQ.'XABCLIB_GMRES') THEN
*
         ISIZE=N*(M+2)+(M+1)*(M+1)+(N-1)/2+1
         ISIZE=ISIZE+N
*
      ELSE IF (SOLVER.EQ.'XABCLIB_BICGSTAB') THEN
*
         ISIZE=11*N+1
*
*     ELSE IF (SOLVER.EQ.'XABCLIB_BICGSTAB_ORG') THEN
*
*        ISIZE=11*N+1
*
      ELSE IF (SOLVER.EQ.'XABCLIB_LANCZOS') THEN
*
         ISIZE=N*(M+1)+2*M*M+7*M+5*NEV+2
         ISIZE=ISIZE+5*M+3
*
      ELSE IF (SOLVER.EQ.'XABCLIB_ARNOLDI') THEN
*
         ISIZE=N*(M+5)+5*M*M+9*M+6*NEV
         ISIZE=ISIZE+M
*
      END IF
*---------------------------[GB]
      SIZE=ISIZE*8.0D-9
*
      RETURN
      END
*
*
      SUBROUTINE OPENATI_GET_WORKSIZE(SOLVER,M,N,NEV,SIZE)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*32  SOLVER
*
*
*  =====================================================================
*  Purpose
*  =======
*
*  Calculation for WORK ARIA's SIZE
*
*  Arguments
*  =========
*  SOLVER  (input) CHARACTER*15
*          target of calculation.
*          The content is XABCLIB_GMRES or XABCLIB_LANCZOS
*
*  M       (input) INTEGER
*          restart revolution
*
*  N       (input) INTEGER
*          The order of the matrix.
*
*  NEV     (input) INTEGER
*          request number of eigenvalue
*
*  SIZE    (output) INTEGER
*          USING WORK AREA's SIZE
*
* ================================================================
*
      IF      (SOLVER.EQ.'XABCLIB_GMRES') THEN
*
         SIZE=1.0D0*N*(M+2)+(M+1.0D0)*(M+1.0D0)+1.0D0*(N-1)/2+1.0D0
         SIZE=SIZE+1.0D0*N
*
      ELSE IF (SOLVER.EQ.'XABCLIB_BICGSTAB') THEN
*
         SIZE=11.0D0*N+1.0D0
*
*     ELSE IF (SOLVER.EQ.'XABCLIB_BICGSTAB_ORG') THEN
*
*        SIZE=11.0D0*N+1.0D0
*
      ELSE IF (SOLVER.EQ.'XABCLIB_LANCZOS') THEN
*
         SIZE=1.0D0*N*(M+1)+2.0D0*M*M+7.0D0*M+5.0D0*NEV+2.0D0
         SIZE=SIZE+1.0D0*5*M+3.0D0
*
      ELSE IF (SOLVER.EQ.'XABCLIB_ARNOLDI') THEN
*
         SIZE=1.0D0*N*(M+5)+5.0D0*M*M+9.0D0*M+6.0D0*NEV
         SIZE=SIZE+1.0D0*M
*
      END IF
*---------------------------[GB]
      SIZE=SIZE*8.0D-9
      IF (SIZE.GT.16.0D0) THEN
        SIZE=16.0D0
      END IF
*
      RETURN
      END
*
*
      SUBROUTINE OPENATI_READ_POLICY(ITYPE,IPOLICY,RPOLICY,IDBG,INFO)
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*     ..
*     .. Array Arguments
      PARAMETER (NUM_POL=20)
      DOUBLE PRECISION  RPOLICY(NUM_POL)
      INTEGER           IPOLICY(NUM_POL)
      INTEGER           ITYPE
*     ..
*     .. Scalar Auruments ..
      INTEGER           IDBG, INFO
*
*
*  =====================================================================
*  Purpose
*  =======
*
*  OpenATI_READ_POLICY sets up the Policy of IPOLICY and RPOLICY.
*
*  Arguments
*  =========
*
*  ITYPE   (input) INTERGER
*           The type of problem.
*                      =  1 : Linear problem.
*                      = 10 : Eigensolve problem.
*  IPOLICY (output) INTEGER array, dimension (NUM_POL)
*  (NUM    Integer parameter set for policy.
*    _POL)  IPOLICY(1)      : The number of threads.
*           IPOLICY(2)      : Set User policy kinds.
*                      =  1 : Execute time.
*                      =  2 : Run time's memory size.
*                      =  3 : Accuracy.
*           IPOLICY(3)      : Set preconditioner kinds.
*                      =  1 : None.
*                      =  2 : Use Jacobi.
*                      =  3 : Use SSOR.
*                      =  4 : Use ILU(0)_diag (ILU0D).
*                      =  5 : Use ILU(0).
*                      =  6 : Use ILUT(P,tau).
*                      = 10 : Use AUTO.
*           IPOLICY(4)      : Value of solver type.
*                      =  1 : Use GMRES(m).
*                      =  2 : Use BiCGStab(By Itoh).
*                      =  3 : Use BiCGStab(By Van der Vorst).
*                      = 10 : Use AUTO.
*                      = 11 : Use Lanczos.
*                      = 12 : Use Arnoldi.
*           IPOLICY(5:20)   : Don't use.
*
*
*  RPOLICY (output) DOUBLE PRECISION array, dimension (NUM_POL)
*  (NUM    Double precision parameter set for policy.
*    _POL)  RPOLICY(1)      : The threshold value of memory size.
*           RPOLICY(2)      : The threshold value of convergence.
*           RPOLICY(3)      : The threshold value of execute time.
*           RPOLICY(4:20)   : Don't use.
*
*  IDBG    (input) INTEGER
*          Debug print status.
*
*  INFO    (output) INTEGER
*          =    0:  Successful exit
*          = -100:  No exist separate word('=') in keyword.
*          = -200:  Invalid the number of threads value at POLICY FILE.
*          = -300:  Illeagal value of key=POLICY
*          = -310:  Illeagal value of key=PRECONDITIONER
*          = -320:  Illeagal value of Key=POLICY
*          = -400:  Invalid the number of ITYPE.
*          Return code.
*
*  =====================================================================
*
      INTEGER  TIME,MEMORY,ACCURACY
      INTEGER  NO, JACOBI, SSOR, ILU0
      INTEGER  SOLVER_XABCLIB_GMRES
      INTEGER  SOLVER_XABCLIB_BICGSTAB
      INTEGER  SOLVER_XABCLIB_BICGSTAB_ORG
      INTEGER  SOLVER_XABCLIB_LANCZOS
      INTEGER  SOLVER_XABCLIB_ARNOLDI
*
      PARAMETER (LP_SMPS   = 1,
     $           LP_POLICY = 2, TIME=1, MEMORY=2, ACCURACY=3, STABLE=4,
     $           LP_PRECON = 3, NO=1, JACOBI=2, SSOR=3, ILU0D=4,
     $                          ILU0=5, ILUT=6, AUTO=10,
     $           LP_SOLVER = 4, SOLVER_XABCLIB_GMRES=1,
     $                          SOLVER_XABCLIB_BICGSTAB=2,
     $                          SOLVER_XABCLIB_BICGSTAB_ORG=3,
     $                          SOLVER_XABCLIB_LANCZOS=11,
     $                          SOLVER_XABCLIB_ARNOLDI=12,
     $           MP_MAXMEM = 1,
     $           MP_RESID  = 2,
     $           MP_MAXTIME= 3)
*
      CHARACTER*80  C80,C80W
      CHARACTER*15  KEY
      CHARACTER*120 WFILE,POLICY_FILE
*
      INTEGER OMP_GET_NUM_THREADS,OMP_GET_MAX_THREADS,
     $        OMP_GET_THREAD_NUM
      INTEGER THREAD_NUM
      THREAD_NUM= OMP_GET_THREAD_NUM()
*
*--------------POLICY REPORT FILE
      IO=16+THREAD_NUM
*
      INFO=0
*++++++ default set
      DO I=1,NUM_POL
         IPOLICY(I)=-1
         RPOLICY(I)=-1.0D0
      ENDDO
*
      IPOLICY(LP_SMPS   ) = OMP_GET_NUM_THREADS()
      IPOLICY(LP_POLICY ) = TIME
      IPOLICY(LP_PRECON ) = ILU0D
*
      IF (ITYPE .EQ. 1) THEN
         IPOLICY(LP_SOLVER ) = SOLVER_XABCLIB_GMRES
      ELSE IF (ITYPE .EQ. 10) THEN
      ELSE
         INFO= -320
         GOTO 9999
      END IF
*
      RPOLICY(MP_MAXMEM ) = -1.0D0
      RPOLICY(MP_RESID  ) = 1.0D-8
      RPOLICY(MP_MAXTIME) = 0.0D0
*
*++++++ POLICY description file
*
      WRITE(WFILE,*) THREAD_NUM
      WFILE='OPENATI_POLICY_INPUT.'//ADJUSTL(WFILE)
      LENG=LEN_TRIM(WFILE)
      WRITE(IO,*) '    OPENATI_POLICY = ',ADJUSTL(WFILE(1:LENG))
      IF (ADJUSTL(WFILE) .EQ. '    ') THEN
         WRITE(6,*) ' [ OpenATI_POLICY ] POLICY FILE NOT DEFINED;
     $  CANNOT GETENV "OPENATI_POLICY" ,SET DEFAULT POLICY'
         GOTO 9999
      ELSE
         POLICY_FILE=ADJUSTL(WFILE)
      END IF
*
      LU=50+THREAD_NUM
      OPEN(LU,FILE=POLICY_FILE,STATUS='OLD')
*
*++++++ READ POLICY DESCRIPTION FILE
*
  100 CONTINUE
      READ(LU,1000,END=9000) C80
 1000 FORMAT(A80)

      C80W=ADJUSTL(C80)
      IF (C80W(1:1) .EQ. '#' .OR. C80W(1:1) .EQ. '') THEN
         GOTO 100
      END IF
      IEQ=INDEX(C80,'=')
      IF (IEQ.EQ.0) THEN
         IF (IDBG.EQ.1) WRITE(6,*) 'KEY NOT SEPARATE "="',C80
         INFO=-100
         GOTO 9000
      END IF
      C80W(1:IEQ-1)=C80(1:IEQ-1)
      KEY(1:15)='               '
      KEY=ADJUSTL(C80W(1:IEQ-1))
*
      IF      (KEY .EQ. 'CPU') THEN
         C80W=ADJUSTL(C80(IEQ+1:80))
         READ(C80W,*) NUMCPU
         MAXTH=OMP_GET_MAX_THREADS()
         IF (NUMCPU.GE.1  .AND. NUMCPU.LE.MAXTH) THEN
            IPOLICY(LP_SMPS   ) = NUMCPU
         ELSE
            INFO=-200
            GOTO 9000
         END IF
      ELSE IF (KEY .EQ. 'POLICY') THEN
         C80W=ADJUSTL(C80(IEQ+1:80))
         IF      (C80W(1:8) .EQ. 'TIME    ') THEN
            IPOLICY(LP_POLICY ) = TIME
         ELSE IF (C80W(1:8) .EQ. 'MEMORY  ') THEN
            IPOLICY(LP_POLICY ) = MEMORY
         ELSE IF (C80W(1:8) .EQ. 'ACCURACY') THEN
            IPOLICY(LP_POLICY ) = ACCURACY
         ELSE IF (C80W(1:8) .EQ. 'STABLE  ') THEN
            IPOLICY(LP_POLICY ) = STABLE
         ELSE
            INFO=-300
            write(6,*) '[ERROR] illeagal value of key=POLICY',C80W
            GOTO 9000
         END IF
      ELSE IF (KEY .EQ. 'PRECONDITIONER') THEN
         C80W=ADJUSTL(C80(IEQ+1:80))
         IF      (C80W(1:8) .EQ. 'NO      ') THEN
            IPOLICY(LP_PRECON ) = NO
         ELSE IF (C80W(1:8) .EQ. 'JACOBI  ') THEN
            IPOLICY(LP_PRECON ) = JACOBI
         ELSE IF (C80W(1:8) .EQ. 'SSOR    ') THEN
            IPOLICY(LP_PRECON ) = SSOR
         ELSE IF (C80W(1:8) .EQ. 'ILU0D   ') THEN
            IPOLICY(LP_PRECON ) = ILU0D
         ELSE IF (C80W(1:8) .EQ. 'ILU0    ') THEN
            IPOLICY(LP_PRECON ) = ILU0
         ELSE IF (C80W(1:8) .EQ. 'ILUT    ') THEN
            IPOLICY(LP_PRECON ) = ILUT
         ELSE IF (C80W(1:8) .EQ. 'AUTO    ') THEN
            IPOLICY(LP_PRECON ) = AUTO
         ELSE
            INFO=-310
            WRITE(6,*) '[ERROR] illeagal value of key=PRECONDITIONER',
     $                 C80W
            GOTO 9000
         END IF
      ELSE IF (KEY .EQ. 'SOLVER') THEN
         C80W=ADJUSTL(C80(IEQ+1:80))
         IF (ITYPE .EQ. 1) THEN
            IF      (C80W(1:32) .EQ. 'XABCLIB_GMRES') THEN
               IPOLICY(LP_SOLVER ) = SOLVER_XABCLIB_GMRES
            ELSE IF (C80W(1:32) .EQ. 'XABCLIB_BICGSTAB') THEN
               IPOLICY(LP_SOLVER ) = SOLVER_XABCLIB_BICGSTAB
*           ELSE IF (C80W(1:32) .EQ. 'XABCLIB_BICGSTAB_ORG') THEN
*              IPOLICY(LP_SOLVER ) = SOLVER_XABCLIB_BICGSTAB_ORG
            ELSE IF (C80W(1:32) .EQ. 'AUTO') THEN
               IPOLICY(LP_SOLVER ) = AUTO
            ELSE
               INFO=-310
               WRITE(6,*) '[ERROR] illeagal value of key=SOLVER',
     $                  C80W
               GOTO 9000
            ENDIF
         ELSE IF (ITYPE .EQ. 10) THEN
            IF      (C80W(1:32) .EQ. 'XABCLIB_LANCZOS') THEN
               IPOLICY(LP_SOLVER ) = SOLVER_XABCLIB_LANCZOS
            ELSE IF (C80W(1:32) .EQ. 'XABCLIB_ARNOLDI') THEN
               IPOLICY(LP_SOLVER ) = SOLVER_XABCLIB_ARNOLDI
            ELSE
               INFO=-310
               WRITE(6,*) '[ERROR] illeagal value of key=SOLVER',
     $                  C80W
               GOTO 9000
            ENDIF
            IPOLICY(LP_PRECON ) = ILU0
         ELSE
            INFO=-400
            GOTO 9000
         END IF
      ELSE IF (KEY .EQ. 'MAXMEMORY') THEN
         C80W=ADJUSTL(C80(IEQ+1:80))
         READ(C80W,*) XMAXMEM
         RPOLICY(MP_MAXMEM) = XMAXMEM
      ELSE IF (KEY .EQ. 'MAXTIME') THEN
         C80W=ADJUSTL(C80(IEQ+1:80))
         READ(C80W,*) XMAXTIME
         RPOLICY(MP_MAXTIME) = XMAXTIME
      ELSE IF (KEY .EQ. 'RESIDUAL') THEN
         C80W=ADJUSTL(C80(IEQ+1:80))
         READ(C80W,*) RESID
         RPOLICY(MP_RESID  ) = RESID
      END IF
      GOTO 100
*
 9000 CONTINUE
      CLOSE(LU)
 9999 CONTINUE
      RETURN
      END
      SUBROUTINE OPENATI_SET_STRATEGY(KIND_SOLVER,KIND_PRECOND,IFILL,
     $                                TAU,
     $                                NTRY,PS_PAIR,IPRELIST,RPRELIST)
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*     ..
*     .. Array Arguments
      DOUBLE PRECISION  RPRELIST(100)
      INTEGER           PS_PAIR(2,100), IPRELIST(100)
      INTEGER           KIND_PRECOND,KIND_SOLVER
*     ..
*
*
*  =====================================================================
*  Purpose
*  =======
*
*  OpenATI_SET_STRATEGY sets Solver and Preconditioner pair for 
*  automatical tuning.
*
*  Arguments
*  =========
*
*  KIND_SOLVER   (input) INTEGER
*           SOLVER name.
*                    1 ='XABCLIB_GMRES'
*                    2 ='XABCLIB_BICGSTAB'
*                   10 ='AUTO'
*
*  =====================================================================
*
      IF      (KIND_SOLVER.EQ.10 .AND. KIND_PRECOND.NE.10) THEN
         NTRY=2
         PS_PAIR(1, 1:2) = KIND_PRECOND  ; IPRELIST( 1:2)= IFILL 
                                           RPRELIST( 1:2)= TAU
         PS_PAIR(2, 1) = 2
         PS_PAIR(2, 2) = 1
      ELSE IF (KIND_SOLVER.NE.10 .AND. KIND_PRECOND.EQ.10) THEN
         NTRY=4
         PS_PAIR(1, 1) = 3  ; IPRELIST( 1)= 0 ; RPRELIST( 1)=1.0D0
         PS_PAIR(1, 2) = 4  ; IPRELIST( 2)= 0 ; RPRELIST( 2)=1.0D-8
         PS_PAIR(1, 3) = 5  ; IPRELIST( 3)= 0 ; RPRELIST( 3)=0.0D0
         PS_PAIR(1, 4) = 6  ; IPRELIST( 4)=10 ; RPRELIST( 4)=1.0D-8
*
         PS_PAIR(2, 1:4) = KIND_SOLVER
      ELSE
         NTRY=8
         PS_PAIR(1, 1) = 3  ; IPRELIST( 1)= 0 ; RPRELIST( 1)=1.0D0
         PS_PAIR(1, 2) = 3  ; IPRELIST( 2)= 0 ; RPRELIST( 2)=1.0D0
         PS_PAIR(1, 3) = 4  ; IPRELIST( 3)= 0 ; RPRELIST( 3)=1.0D-8
         PS_PAIR(1, 4) = 4  ; IPRELIST( 4)= 0 ; RPRELIST( 4)=1.0D-8
         PS_PAIR(1, 5) = 5  ; IPRELIST( 5)= 0 ; RPRELIST( 5)=0.0D0
         PS_PAIR(1, 6) = 5  ; IPRELIST( 6)= 0 ; RPRELIST( 6)=0.0D0
         PS_PAIR(1, 7) = 6  ; IPRELIST( 7)=10 ; RPRELIST( 7)=1.0D-8
         PS_PAIR(1, 8) = 6  ; IPRELIST( 8)=10 ; RPRELIST( 8)=1.0D-8
*
         PS_PAIR(2, 1) = 2
         PS_PAIR(2, 2) = 1
         PS_PAIR(2, 3) = 2
         PS_PAIR(2, 4) = 1
         PS_PAIR(2, 5) = 2
         PS_PAIR(2, 6) = 1
         PS_PAIR(2, 7) = 2
         PS_PAIR(2, 8) = 1
*
      END IF
*
*
*
*
 9000 CONTINUE
      RETURN
      END
