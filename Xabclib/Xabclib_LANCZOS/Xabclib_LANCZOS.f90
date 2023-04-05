      SUBROUTINE Xabclib_LANCZOS( N, NNZ, IRP, ICOL, VAL,
     $                            NEV, EV, EVEC, LDE,
     $                            IATPARAM, RATPARAM,
     $                            WORK, LWORK, IWORK, LIWORK, INFO )
*
*     Xabclib_LANCZOS : Compute eigenpairs for symmetric sparse matrix
*                       by Explicitly Restart Deflated LANCZOS algorithm.
*     December 2008
*
      IMPLICIT DOUBLE PRECISION ( A-H, O-Z )
*
*     .. Scalar Arguments ..
      INTEGER            N, NNZ, NEV, LDE, LWORK, LIWORK, INFO
*     ..
*     .. Array Arguments ..
      INTEGER            IRP( N+1 ), ICOL( NNZ ), IWORK( LIWORK )
      INTEGER            IATPARAM(50)
*
      DOUBLE PRECISION   VAL( NNZ ), EV( NEV ), EVEC( LDE, NEV )
      DOUBLE PRECISION   WORK( LWORK )
      DOUBLE PRECISION   RATPARAM(50)
*
*     ..
*
*  =====================================================================
*  Purpose
*  =======
*
*  Xabclib_LANCZOS computes NEV eigenpairs of a real symmetric sparse 
*  matrix by Expilictly Restart Deflated LANCZOS algorithm.
*
*  if      IATPARAM(30) = 1 : compute eigenvalues are descending order,
*  else if IATPARAM(30) = 2 : compute absolute eigenvalue are
*                         descending order.
*
*  The Explicitly Restart Deflated Lanczos algorithm is based on 
*  Krylov subspace method. 
*
*
*  Arguments
*  =========
*
*
*  N       (input) INTEGER
*          The order of the matrix.  N >= 0.
*
*  NNZ     (input) INTEGER
*          Non-Zeros of the matrix.  NNZ >= N.
*
*  IRP     (input) INTEGER array
*  (N+1)   Diagonal pointer of the matrix in Harwell-Boeing symmetric
*          format. 
*
*  ICOL    (input) INTEGER array
*  (NNZ)   Column indecies of the matrix in Harwell-Boeing symmetric
*          format. 
*
*  VAL     (input) DOUBLE PRECISION array 
*  (NNZ)   Value of the matrix in Harwell-Boeing symmetric format. 
*
*  NEV     (input) INTEGER
*          The number of desired eigenpairs.
*
*  EV(NEV) (output) DOUBLE PRECISION array
*          Eigenvalues,
*            if iatparam(30) = 1 : Eigenvalues are Descending order.
*            if iatparam(30) = 2 : Absolute Eigenvalues are Descending
*                              order.
*
*  EVEC    (output) DOUBLE PRECISION array, dimension (LDE, NEV)
*  (LDE,   EVEC contains the orthonormal eigenvectors of the matrix VAL,
*  NEV)    with the i-th column of EVEC holding the eigenvector 
*          associated with EV(i).
*
*  LDE     (input) INTEGER
*          The leading dimension of the array EVEC. LDE >= N.
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
*           IATPARAM(7)     : OpenATI_DSRMV auto-tuned on/off.
*           IATPARAM(8)     : Fastest OpenATI_DSRMV impl. method.
*           IATPARAM(12)    : Type of Gram-Schmidt procedure.
*           IATPARAM(13)    : DGKS refinement done or not.
*
*          IATPARAM(21:40) : Xabclib parameter set.
*           IATPARAM(21)    : Number of threads.
*           IATPARAM(22)    : Max. iterations.
*           IATPARAM(23)    : Number of iterations.
*           IATPARAM(27)    : Input size of Krylov subspace.
*           IATPARAM(28)    : Start size of Krylov subspace at
*                             subspace expand AT-on.
*           IATPARAM(29)    : Final size of Krylov subspace.
*           IATPARAM(30)    : Eigenvalue choice option.
*           IATPARAM(31)    : Mat-Vec times.
*           IATPARAM(50)    : Debug print control flag.
*
*  RATPARAM (input/output) DOUBLE PRECISION array, dimension (50)
*          Double precision parameter set for OpenATI and Xabclib.
*           RATPARAM(1:2)   : Mandatory.
*
*          RATPARAM(3:20) : OpenATI parameter set.
*           RATPARAM(3)     : Reserved.
*           RATPARAM(4)     : Threshhold of MM-ratio.
*
*          RATPARAM(21:40) : Xabclib parameter set.
*           RATPARAM(22)    : Max. elapsed time(limit time)
*           RATPARAM(23)    : Convergence criterion.
*           RATPARAM(29)    : 2-norm of max. residual.
*           RATPARAM(30)    : Floating operations(*10^9 operations).
*           RATPARAM(32)    : Total solve time(elapsed).
*
*  WORK    (workspace/output) DOUBLE PRECISION array,
*  (LWORK)   LWORK >= (1+MSIZE)*N + 2*MSIZE*MSIZE + 7*MSIZE + 5*NEV +2.
*
*  IWORK   (workspace/output) INTEGER array,
*  (LIWORK)  LIWORK >= 5*MSIZE + 3.
*
*  INFO    (output) INTEGER
*          =   0:  Successful exit
*          <   0:  If INFO = -i, the i-th argument had an illegal value
*          = 100:  Breakdown occured because of Zero vector generated.
*          = 200:  The algorithm failed to converge in DSTEVD(LAPCAK).
*          = 300:  The algorithm failed to converge in Max restart.
*          = 400:  Come to Upper Limit of CPU time
*          = 500:  Not Enough free memory.(IATPARAM(8)=12,13)
*
*  Local variables
*  ===============
*  IDBG      INTEGER
*            Debug print flag.
*             IDBG=1 : Print debug info.
*  NUM_SMP   INTEGER
*            OMP_GET_MAX_THREADS() 
*  ICASE     INTEGER
*            Selected Algorithm NO. of matrix-vector product
*  WK        DOUBLE PRECISION dynamic array
*            workspace for parallel vector reduction
*            array size : (N,NUM_SMP)
*  =====================================================================
*     ..
*     .. Local variables
*
      INTEGER            LSINF, LWK, IERR_ALLOC
*
      INTEGER            IDBG, MSIZE, NUM_SMP, MV_AT_SYM, MVCASE_SYM
      INTEGER            I, ICASE, JCASE, ICK, M, INFO_A
      INTEGER            IP1, IP2, IP3, IP4, IP5, IP6, IP7, IP8
      INTEGER            IP_WK, IP_SI
      INTEGER            LRWK, LIWK
      DOUBLE PRECISION   WK
      ALLOCATABLE:: WK(:)
*
      EXTERNAL OMP_GET_MAX_THREADS
      INTEGER  OMP_GET_MAX_THREADS
*
      IDBG     =IATPARAM(50)
      MSIZE    =IATPARAM(27)
      NUM_SMP  =IATPARAM(3)
      MV_AT_SYM=IATPARAM(7)
*
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      IERR_ALLOC = -1
      IATPARAM(31) = 0
*
      RATPARAM(29) = 0.0D0
*--------------------------Floating Ope. counts
      RATPARAM(30)=0.0D0
*
*-------------arg-1
      IF (N .LE. 0 ) THEN
         INFO=-1
         GOTO 100
      ENDIF
*-------------arg-2
      IF (NNZ .NE. IRP(N+1)-1 ) THEN
         INFO=-2
         GOTO 100
      ENDIF
*-------------arg-3
*-------------arg-4
      DO I=1,N
         IF (ICOL(IRP(I)) .NE. I) THEN
            WRITE(6,*) '[Error] ERROR IRP DATA I,IRP(I),ICOL(IRP(I))=',
     &                  I,IRP(I),ICOL(IRP(I))
            INFO=-4
            GOTO 100
         ENDIF
      ENDDO
*-------------arg-5
*-------------arg-6
      IF (NEV .LE. 0 .OR. NEV .GE.N ) THEN
         INFO=-6
         GOTO 100
      ENDIF
*-------------arg-7
*-------------arg-8
*-------------arg-9
      IF (LDE .LT. N ) THEN
         INFO=-9
         GOTO 100
      ENDIF
*-------------arg-10
*-------------arg-11
*-------------arg-12
*-------------arg-13
      ICK=(1+MSIZE)*N + 2*MSIZE*MSIZE + 7*MSIZE + 5*NEV +2
      IF (LWORK .LT. ICK ) THEN
         INFO=-13
         WRITE(6,*) '[Error] PARAMETER ERROR =',LWORK,ICK
         GOTO 100
      ENDIF
*-------------arg-14
*-------------arg-15
      ICK=5*MSIZE + 3
      IF (LIWORK .LT. ICK ) THEN
         INFO=-15
         GOTO 100
      ENDIF
*
  100 CONTINUE
      IF (INFO .NE. 0 ) THEN
         WRITE(6,*) '[Error] PARAMETER ERROR =',INFO
         GOTO 3000
      END IF
*
*-------------------------quick return in N=1
*
      IF( N.EQ.1 ) THEN
         EV(1)=VAL(1)
         EVEC( 1, 1 ) = 1.0D0
         RETURN
      END IF
*
      M=MSIZE
*
*     separate work
*
      LRWK=1+4*M+M*M
      LIWK=3+5*M
*
      IP1=1
      IP2=IP1+N
      IP3=IP2+N*M
      IP4=IP3+M
      IP5=IP4+M+1
      IP6=IP5+M
      IP7=IP6+M*M
      IP8=IP7+LRWK
*
*
************************************************************************
*  If IATPARAM(7) = 2 or 3 then search the fastest Matrix-Vector product
************************************************************************
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
         CALL OpenATI_DSRMV(N,NNZ,IRP,ICOL,VAL,WORK,WORK(IP2),
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
            GO TO 3000
         END IF
*
         IP_WK=1
         IP_SI=IP_WK+LWK
*
         CALL OpenATI_DSRMV_Setup(N,NNZ,IRP,ICOL,
     &        IATPARAM,RATPARAM,WK(IP_SI),LSINF,INFO_A)
      ENDIF
*
      CALL Xabclib_LANZ01(N,NNZ,IRP,ICOL,VAL,
     &                    NEV,EV,EVEC,LDE,MSIZE,
     &                    IATPARAM,RATPARAM,
     &                    WORK(IP1),WORK(IP2),WORK(IP3),WORK(IP4),
     &                    WORK(IP5),WORK(IP6),WORK(IP7),LRWK,IWORK,
     &                    LIWK,INFO,
     &                    WK(IP_WK),WORK(IP8),
     &                    WK(IP_SI),LSINF)
*
*     End of Xabclib_LANCZOS
      IF (IERR_ALLOC.EQ.0) THEN
         DEALLOCATE(WK)
      ENDIF
*
 3000 CONTINUE
      RETURN
      END
*
      SUBROUTINE Xabclib_LANZ01(N,NNZ,IRP,ICOL,VAL,
     &                          NEV,EV,EVEC,LDE,MSIZE,
     &                          IATPARAM,RATPARAM,
     &                          RVEC,Q,ALF,BET,RESID,Z,RWK,LRWK,IWK,
     &                          LIWK,INFO,
     &                          WK,SAMP,SINF,LSINF)
*
      IMPLICIT DOUBLE PRECISION ( A-H, O-Z )
*
*     .. Scalar Arguments ..
      INTEGER            N, NNZ, NEV, LDE, MSIZE
      INTEGER            LRWK, LIWK, INFO, LSINF
*     ..
*     .. Array Arguments ..
      INTEGER            IRP( N+1 ), ICOL( NNZ )
      INTEGER            IWK( LIWK )
*
      DOUBLE PRECISION   VAL( NNZ ), EV( NEV ), EVEC( LDE, NEV )
      DOUBLE PRECISION   RVEC(N), Q(N,MSIZE)
      DOUBLE PRECISION   ALF(MSIZE), BET(0:MSIZE)
      DOUBLE PRECISION   RESID(MSIZE), Z(MSIZE,MSIZE)
      DOUBLE PRECISION   RWK( LRWK ), WK(N,*)
      DOUBLE PRECISION   SAMP(5,NEV)
      DOUBLE PRECISION   SINF(LSINF)
*
      INTEGER            IATPARAM(50)
      DOUBLE PRECISION   RATPARAM(50)
*
*     ..
*
*  =====================================================================
*  Purpose
*  =======
*
*  Xabclib_LANCZ01: Controller of The Explicitly Restart Deflated
*                   Lanczos algorithm.
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          The order of the matrix.  N >= 0.
*
*  NNZ     (input) INTEGER
*          Non-Zeros of the matrix.  NNZ >= N.
*
*  IRP     (input) INTEGER array
*  (N+1)   Diagonal pointer of the matrix in Harwell-Boeing symmetric
*          format.
*
*  ICOL    (input) INTEGER array
*  (NNZ)   Column indecies of the matrix in Harwell-Boeing symmetric
*          format.
*
*  VAL     (input) DOUBLE PRECISION array
*  (NNZ)   Value of the matrix in Harwell-Boeing symmetric format.
*
*  NEV     (input) INTEGER
*          The number of desired eigenpairs.
*
*  EV(NEV) (output) DOUBLE PRECISION array
*          Eigenvalues,
*            if iparm(1) = 1 : Eigenvalues are Descending order.
*            if iparm(1) = 2 : Absolute Eigenvalues are Descending
*                              order.
*
*  EVEC    (output) DOUBLE PRECISION array, dimension (LDE, NEV)
*  (LDE,   EVEC contains the orthonormal eigenvectors of the matrix VAL,
*  NEV)    with the i-th column of EVEC holding the eigenvector
*          associated with EV(i).
*
*  LDE     (input) INTEGER
*          The leading dimension of the array EVEC. LDE >= N.
*
*  MSIZE   (input) INTEGER
*          Krylov subspace size. at least MSIZE >=NEV.
*
*  IPARM   (input/output) INTEGER array
*  (10)    Integer parameter set.
*           IPARM(1) =  1 : Eigenvales are Descending order in EV.
*                    =  2 : Absolute Eigenvales are Descending order
*                           in EV.
*           IPARM(2)      : Max iterations of Restarts. >=1
*           IPARM(3)      : (return) Iterations of Restarts.
*           IPARM(4)      : (return) Krylov subspace dimension.
*           IPARM(5:10)   : Don't use.
*
*  RPARM   (input) DOUBLE PRECISION array
*  (10)    DOUBLE PRECISION parameter set.
*           RPARM(1)      : Convergence criteria.
*           RPARM(2)      : Upper Limit of CPU time
*                           if RPARM(2) <= 0.0e0, then no limit of
*                           CPU time. 
*                         : (return) Elapsed time for solver.
*           RPARM(3)      : OpenATI_DAFRT's RPARAM(1) value.
*                           criterion of MM-ratio.
*           RPARM(4:10)   : Don't use.
*
*  IAT     (input) INTEGER array
*  (10)    Automatically Tuning parameter set.
*           IAT(1)   =  1 : Automatically Restart by tuned MSIZE.
*           IAT(2)   =  1 : Automatically tuned Matrix-vector product.
*                    =  2 : Automatically tuned Matrix-vector product 
*                           to detect a residual memory.
*           IAT(3)        : (return) Algorithm NO. of
*                           matrix-vector product.
*
*  RVEC(N) (Workspace) DOUBLE PRECISION array
*          Working vector for subspace.
*
*  Q(N,    (Workspace) DOUBLE PRECISION array
*  MSIZE)  Krylov subspace vectors.
*
*  ALF     (Workspace) DOUBLE PRECISION array
*  (MSIZE) Diagonal elements of Tridiagonal system.
*          On exit of DSTEVD, ALF is Ritz values.
*
*  BET(0:  (Workspace) DOUBLE PRECISION array
*  (MSIZE) Sub-Diag. elements of Tridiagonal system in BET(1:MSIZE-1).
*
*  RESID   (Workspace) DOUBLE PRECISION array
*  (MSIZE) Residual norm with normalized eigenvalue.
*
*  Z(MSIZE (Workspace) DOUBLE PRECISION array
*  ,MSIZE) Eigenvectors of Tridaiagonal system.
*
*  RWK     (Workspace) DOUBLE PRECISION array,
*  (LRWK)  RWK is workspace for DSTEVD. LRWK >=MSIZE*MSIZE+4*MSIZE+1.
*
*  LRWK    (input) INTEGER
*          Length of RWK.  LRWK >= MSIZE*MSIZE + 4*MSIZE + 1.
*
*  IWK     (workspace) INTEGER array,
*  (LIWK)  IWK is workspace for DSTEVD. LIWORK >= 5*MSIZE + 3.
*          Another purpus, 
*               IWK for Sorted Ritz value list for IPARM(1)=2.
*
*  LIWK    (input) INTEGER
*          Length of IWK. LIWORK >= 5*MSIZE + 3. 
*
*  INFO    (output) INTEGER
*          =   0:  successful exit
*          = 100:  Breakdown occured because of Zero vector generated.
*          = 200:  the algorithm failed to converge in DSTEVD(LAPCAK).
*          = 300:  the algorithm failed to converge in Max restart.
*          = 400:  Come to Upper Limit of CPU time
*
*  WK      (Workspace)  DOUBLE PRECISION dynamic array
*            workspace for parallel vector reduction
*            array size : (N,NUM_SMP)
*           (NUM_SMP = IATPARAM(3))
*
*  SAMP    (Workspace)  DOUBLE PRECISION array
*  (5,NEV)   Latest 5-restart history of residual.
*
*  SINF    (input) DOUBLE PRECISION array, dimension (LSINF)
*          Setup Information for OpenATI_DSRMV
*
*  LSINF   (input) INTEGER
*          The dimension of SINF.
*          if ICASE = 11 :
*             LSINF >= 0
*          else if ICASE = 12 :
*             LSINF >= INT(0.5*NUM_SMP)+1
*          else if ICASE = 13 :
*             LSINF >= N+NUM_SMP+3
*           (NUM_SMP = IATPARAM(3), ICASE = IATPARAM(8))
*
*  Local variables
*  ===============
*
*  LOCK    INTEGER
*          The number of locked Ritz pair.
*
*  IS_CONV INTEGER
*          The number of Converged Ritz pair in current Lanczos step.
*
*  KSIZE   INTEGER
*          The number of Redundant Krylov subspace.
*
*  MSIZE2  INTEGER
*          Automatically tuned Krylov subspace.
*
*  =====================================================================
*     ..
*     .. Local variables
*
      INTEGER            I, J, IR, M, KSIZE
      INTEGER            INFO_A, INFO_M, LAINFO
      INTEGER            MAXRESTART, MDEL, ISAMP, MSIZE2
      INTEGER            ITCHK, NORMALFLG, IS_CONV, LOCK, IDGKS
      INTEGER            IACCEL, IRT, MSIZE_INI, IGSPOL, MSIZE_AT
      INTEGER            IORDER, IDBG
*
      DOUBLE PRECISION   S, BINV, BETM
      DOUBLE PRECISION   TTIME0, TTIME1, CTIME1
      DOUBLE PRECISION   CPUTIME, OMP_GET_WTIME
      DOUBLE PRECISION   MAX_ETIME, MM_RATIO
*
*     CALL CPU_TIME(TTIME0)
      TTIME0=OMP_GET_WTIME()
*
      MAXRESTART=IATPARAM(22)
*
      MAX_ETIME=RATPARAM(22)
      MSIZE_INI=IATPARAM(28)
      IGSPOL   =IATPARAM(12)
      MSIZE_AT =IATPARAM( 4)
      MM_RATIO =RATPARAM( 4)
      IORDER   =IATPARAM(30)
*
      IDBG=IATPARAM(50)
*
      FLOP=0.0D0
*
      ITCHK=1
      IF (MAX_ETIME.LE.0.0D0) THEN
         ITCHK=0
      ENDIF
*
*---------covergence criterion check
*
      IF (RATPARAM(23) .LT. 1.0D-15 ) THEN
       WRITE(6,*) "Warning by Xabclib: Convergence criterion
     $ (RATPARAM(23)) is too small! RATPARAM(23) = ", RATPARAM(23)
      END IF  
*
*---------starting vector
*
      DO I=1,N-1,2
        RVEC(I  )= 1.0D0
        RVEC(I+1)=-1.0D0
      ENDDO
      IF (MOD(N,2).EQ.1) RVEC(N)=1.0D0
*
      S=0.0D0
      DO I=1,N
        S=S+RVEC(I)*RVEC(I)
      ENDDO
      BET(0)=SQRT(S)
      FLOP=FLOP+2.0D0*N
************************************************************************
*        Msize Autotuning
************************************************************************
      IF (MSIZE_AT.EQ.1) THEN
*------------------------initial Krylov size
         ISAMP=0
         IF (MSIZE_INI .LT. NEV) THEN
            MSIZE2=NEV
         ELSE
            MSIZE2=MSIZE_INI
         END IF
         IF (MSIZE2 .GT. MSIZE) MSIZE2=MSIZE
*------------------------OpenATI_DAFRT option
         IF (MM_RATIO.LE.0.0D0) THEN
            RATPARAM(4)=100.0D0
         ENDIF
*------------------------History of residual
         DO I=1,NEV
            DO J=1,5
               SAMP(J,I)=0.0D0
            ENDDO
         ENDDO
      ELSE
         MSIZE2=MSIZE
      ENDIF
*
*......................................................
*------------------------Gram-Schmidt Policy Scenario
      IF (IDBG.EQ.1) THEN
         WRITE(6,*) ' Gram Schmidt Policy=',IGSPOL
      ENDIF
      NORMALFLG=0
*......................................................
*
      IS_CONV=0
      LOCK=0
*
      DO IR=1,MAXRESTART
*
         IF (IDBG.EQ.1) THEN
            WRITE(6,*) ' >>>>> IR= ',IR,' MSIZE2=',MSIZE2
         ENDIF
*
         DO M=LOCK+1,MSIZE2
*>>>>breakdown
            IF (DABS(BET(M-1)) .LE. 1.0D-15) THEN
               INFO=100
               RETURN
            ENDIF
            BINV=1.0D0/BET(M-1)
            DO I=1,N
               Q(I,M)=RVEC(I)*BINV
            ENDDO
*
            CALL OpenATI_DSRMV(N,NNZ,IRP,ICOL,VAL,Q(1,M),RVEC,
     &                         IATPARAM,RATPARAM,WK,SINF,LSINF,INFO_M)
*
            S=0.0D0
            DO I=1,N
               S=S+Q(I,M)*RVEC(I)
            ENDDO
            ALF(M)=S
*
            IF (M.EQ.1) THEN
               DO I=1,N
                  RVEC(I)=RVEC(I)-ALF(M)*Q(I,M)
               ENDDO
            ELSE
               DO I=1,N
                  RVEC(I)=RVEC(I)-ALF(M)*Q(I,M)-BET(M-1)*Q(I,M-1)
               ENDDO
            ENDIF
*
            CALL OpenATI_DAFGS
     &               (NORMALFLG,N,RVEC,Q,N,M-1,RWK,
     &                IATPARAM,RATPARAM,INFO_M)
*
            S=0.0D0
            DO I=1,N
              S=S+RVEC(I)*RVEC(I)
            ENDDO
            BET(M)=SQRT(S)
*
         ENDDO              !M-LOOP(expand subspace loop)
         BETM=BET(MSIZE2)
         FLOP=FLOP+(MSIZE2-LOCK)*(4.0D0*NNZ+8.0D0*N+4.0D0*N*(M-1))
         IATPARAM(31)=IATPARAM(31)+MSIZE2-LOCK
*
*-----------------------------------------------------------------------
*        solve tridiagonal system by lapack(devide & conquer)
*        ( on exit ALF(LOCK+1:MSIZE2) are eigenvalue with acending order
*-----------------------------------------------------------------------
         KSIZE=MSIZE2-LOCK
         CALL DSTEVD('V',KSIZE,ALF(LOCK+1),BET(LOCK+1),Z(LOCK+1,LOCK+1),
     &               MSIZE,RWK,LRWK,IWK,LIWK,LAINFO)
         IF (LAINFO.NE.0) THEN
            INFO=200
            RETURN
         ENDIF
*-----------------------------------------------------------------------
*        permutation list of absolute eigenvalue at IPARM(1)=2
*-----------------------------------------------------------------------
         IF (IORDER.EQ.2) THEN
            DO I=LOCK+1,MSIZE2
               RWK(I) = DABS(ALF(I))
            ENDDO
            CALL Xabclib_QSORTD(IWK(LOCK+1),KSIZE,RWK(LOCK+1))
         ENDIF
*-----------------------------------------------------------------------
*        (1) residual estimate
*        (2) convergence checking
*        (3) create Ritz pair
*        (4) creat Shur vector for restart
*        (5) deflation
*-----------------------------------------------------------------------
         CALL Xabclib_LANZ03(LOCK,NEV,KSIZE,Z(LOCK+1,LOCK+1),MSIZE,
     &                       ALF(LOCK+1),Q,N,BETM,IATPARAM,RATPARAM,
     &                       IWK(LOCK+1),RVEC,EV,EVEC,LDE,RESID,IS_CONV,
     &                       INFO,IDBG,FLOP)
*
*        CALL CPU_TIME(TTIME1)
         TTIME1=OMP_GET_WTIME()
         CPUTIME=TTIME1-TTIME0
         IF ( ITCHK.EQ.1 .AND. CPUTIME.GT.MAX_ETIME ) THEN
            INFO=400
            GOTO 1000
         ENDIF
         IF (INFO.NE.0) THEN
            GOTO 1000
         ENDIF
         IF (IDBG.EQ.1) WRITE(6,*) ' AFTER Deflation LOCK=',LOCK
         IF (LOCK .EQ. NEV) THEN
            RESMAX=0.0D0
            DO I=1,NEV
               RESMAX=MAX(RESMAX,RESID(I))
            ENDDO
            RATPARAM(29)=RESMAX
            GO TO 1000
         ELSE
*------------------Shur vector orthonomalize to locked subspace vectors
            CALL OpenATI_DAFGS
     &               (NORMALFLG,N,RVEC,Q,N,LOCK,RWK,
     &                IATPARAM,RATPARAM,INFO_M)
            S=0.0D0
            DO I=1,N
               S=S+RVEC(I)*RVEC(I)
            ENDDO
            BET(LOCK)=SQRT(S)
            FLOP=FLOP+4.0D0*N*LOCK+2.0D0*N
************************************************************************
*     Automatically subspace expand with Max-Min ratio
************************************************************************
            IF (MSIZE_AT.EQ.1) THEN
               IF (ISAMP.EQ.5) ISAMP=0
               ISAMP=ISAMP+1
               DO I=1,NEV
                  SAMP(ISAMP,I)=RESID(I)
               ENDDO
*
               IACCEL=0
               DO I=LOCK+1,NEV
                  IF (IDBG.EQ.1) THEN
                     WRITE(6,*) ' OpenATI_DAFRT Eigenvalue NO.',I
                  ENDIF
                  CALL OpenATI_DAFRT(5,SAMP(1,I),IRT,
     &                               IATPARAM,RATPARAM,INFO_A)
                  IF (IRT .EQ. 0) THEN
                     IACCEL=IACCEL+1
                  ENDIF
               ENDDO
*--------------------------ALL MM-ratio are stagnation.
               IF ( IACCEL .EQ. 0 ) THEN
                  MSIZE2=MSIZE2+IATPARAM(5)
                  IF (MSIZE2.GT.MSIZE) MSIZE2=MSIZE
               ENDIF
            ENDIF
*
         ENDIF
*
*
      ENDDO              !IR-LOOP(restart loop)
*
*-----------------------------------------------------------------------
*     IR-LOOP finished => no convergence
*-----------------------------------------------------------------------
      INFO=300
*
*
*----------------- lanczos successfully ended --------------------------
 1000 CONTINUE
*     CALL CPU_TIME(CTIME1)
      CTIME1=OMP_GET_WTIME()
      IF (IDBG.EQ.1) 
     &    WRITE(6,*) ' Xabclib_LANCZOS Time=',ctime1-ttime0,'[sec]'
*-----------------------------------------------------------------------
      IATPARAM(23)=IR
      IATPARAM(29)=MSIZE
      IF (MSIZE_AT .EQ. 1) THEN
         IATPARAM(29)=MSIZE2
      END IF
      RATPARAM(30)=FLOP*1.0D-9
      RATPARAM(32)=CTIME1-TTIME0
*-----------------------------------------------------------------------
*
*
      RETURN
      END
      SUBROUTINE Xabclib_LANZ03(LOCK,NEV,KSIZE,Z,MSIZE,
     &                          THET,Q,N,BETM,IATPARAM,RATPARAM,
     &                          IND,RVEC,EV,EVEC,LDE,RESID,IS_CONV,
     &                          INFO,IDBG,FLOP)
*
      IMPLICIT DOUBLE PRECISION ( A-H, O-Z )
*
*     .. Scalar Arguments ..
      INTEGER            LOCK, NEV, KSIZE, MSIZE, N, LDE, IS_CONV
      INTEGER            INFO
      DOUBLE PRECISION   BETM
*
*     ..
*     .. Array Arguments ..
      INTEGER            IND( KSIZE )
*
      DOUBLE PRECISION   Z(MSIZE,KSIZE)
      DOUBLE PRECISION   THET( KSIZE ), Q( N, MSIZE )
      DOUBLE PRECISION   RVEC( N ), EV( NEV ), EVEC( LDE, NEV )
      DOUBLE PRECISION   RESID( MSIZE )
*
      INTEGER            IATPARAM(50)
      DOUBLE PRECISION   RATPARAM(50)
*     ..
*
*  =====================================================================
*  Purpose
*  =======
*
*  Xabclib_LANCZ03:  (1) residual estimate 
*                    (2) convergence checking 
*                    (3) create Ritz pair
*                    (4) creat Shur vector for restart
*                    (5) deflation
*
*  Arguments
*  =========
*
*  LOCK    (input/output) INTEGER
*          The number of locked Ritz pair.
*
*  NEV     (input) INTEGER
*          The number of desired eigenpairs.
*
*  KSIZE   (input) INTEGER
*          The number of Redundant Krylov subspace.
*
*  Z(MSIZE (input) DOUBLE PRECISION array
*  ,KSIZE) Eigenvectors of Deflated Tridaiagonal system.
*
*  MSIZE   (input) INTEGER
*          The leading dimension of the array Z.
*
*  THET    (input) DOUBLE PRECISION array
*  (KSIZE) Ritz values.
*
*  Q(N,    (input) DOUBLE PRECISION array
*  MSIZE)  Krylov subspace vectors. Q(*,1:LOCK) are Locked Shur vectors.
*
*  N       (input) INTEGER
*          The order of the matrix.
*
*  BETM    (input) DOUBLE PRECISION
*          The 2-norm of last subspace vector.
*
*  IPARM   (input) INTEGER array
*  (10)    Integer parameter set.
*           IPARM(1) =  1 : Eigenvales are Descending order in EV.
*                    =  2 : Absolute Eigenvales are Descending order
*                           in EV.
*
*  RPARM   (input) DOUBLE PRECISION array
*  (10)    DOUBLE PRECISION parameter set.
*           RPARM(1)      : Convergence criteria.
*
*  IND     (input) INTEGER array
*  (KSIZE) Sorted Absolute Ritz value list in IPARM(1)=2.
*
*  RVEC(N) (output) DOUBLE PRECISION array
*          Shur vector for Initial vecor of next Lanczos step.
*
*  EV(NEV) (input) DOUBLE PRECISION array
*          Eigenvalues,
*            if iparm(1) = 1 : Eigenvalues are Descending order.
*            if iparm(1) = 2 : Absolute Eigenvalues are Descending
*                              order.
*
*  EVEC    (input) DOUBLE PRECISION array, dimension (LDE, NEV)
*  (LDE,   EVEC contains the orthonormal eigenvectors of the matrix VAL,
*  NEV)    with the i-th column of EVEC holding the eigenvector
*          associated with EV(i).
*
*  LDE     (input) INTEGER
*          The leading dimension of the array EVEC. LDE >= N.
*
*  RESID   (Workspace) DOUBLE PRECISION array
*  (MSIZE) Residual norm with normalized eigenvalue.
*
*  MSIZE   (input) INTEGER
*          Krylov subspace size. at least MSIZE >=NEV.
*
*  IS_CONV INTEGER
*          The number of Converged Ritz pair in current Lanczos step.
*
*  INFO    (output) INTEGER
*          =   0:  successful exit
*          = 100:  Breakdown occured because of Zero vector generated.
*
*  Local variables
*  ===============
*
*  NEW_LOCK  INTEGER
*            The number of new locked Ritz pair.
*  IDBG      INTEGER
*            If ENV. Var. 'ABCLIB_DEBUG=1' then
*               IDBG = 1 and Debug print
*  =====================================================================
*     ..
*     .. Local variables
      INTEGER            I, J, K, L
      INTEGER            IC, NEW_LOCK
      INTEGER            IDBG, IORDER
*
      DOUBLE PRECISION   EPS, EE
      DOUBLE PRECISION   S, S2, SINV
      DOUBLE PRECISION   BET( MSIZE )
*
      EPS    =RATPARAM(23)
      IORDER =IATPARAM(30)
*     IF (EPS .LT. 1.0D-15 )  EPS=1.0D-8
*
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*     (1) residual estimate with Lanczos relation
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*
      IF (IORDER.EQ.1) THEN
*-----------------------------------------------------------------------
*            eigenvalue ascending order
*-----------------------------------------------------------------------
         DO K=1,NEV-LOCK
            EE=DABS(THET(KSIZE-K+1))
            IF (EE .LT. 1.0D-15) THEN
               EE=1.0D0
            END IF
            RESID(LOCK+K)=DABS(BETM*Z(KSIZE,KSIZE-K+1))/EE
            IF (IDBG.EQ.1) 
     $          WRITE(6,*) 'K=',LOCK+K,'EE,RESID=',EE,RESID(LOCK+K)
         ENDDO
      ELSE IF (IORDER.EQ.2) THEN
*-----------------------------------------------------------------------
*            absolute eigenvalue ascending order
*-----------------------------------------------------------------------
         DO K=1,NEV-LOCK
            IC=IND(KSIZE-K+1)
            EE=DABS(THET(IC))
            IF (EE .LT. 1.0D-15) THEN
               EE=1.0D0
            END IF
            RESID(LOCK+K)=DABS(BETM*Z(KSIZE,IC))/EE
            IF (IDBG.EQ.1) 
     $          WRITE(6,*) 'K=',LOCK+K,'EE,RESID=',EE,RESID(LOCK+K)
         ENDDO
      ENDIF
*
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*     (2) convergence checking with contiguous mode
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      IS_CONV = 0
      DO K=1,NEV-LOCK
         IF (RESID(LOCK+K) .LE. EPS) THEN
            IS_CONV = IS_CONV + 1
         ELSE
            GOTO 1000
         ENDIF
      ENDDO
 1000 CONTINUE
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*     (3) create Ritz pair (new locked vectors)
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      NEW_LOCK = MIN(IS_CONV,NEV-LOCK)
      IF (IDBG.EQ.1) WRITE(6,*) 'LOCK,NEW_LOCK=',LOCK,NEW_LOCK
      DO L=1,NEW_LOCK
         IF (IORDER.EQ.1) THEN
            IC=KSIZE-L+1
         ELSE
            IC=IND(KSIZE-L+1)
         ENDIF
*>>>>>>>>>>>>>>>>>>Ritz pair
         EV(LOCK+L) = THET(IC)
         S2=0.0D0
!$omp parallel
!$omp do reduction(+:S2),private(s)
         DO I=1,N
            S=0.0D0
            DO J=1,KSIZE
               S=S+Q(I,LOCK+J)*Z(J,IC)
            ENDDO
            EVEC(I,LOCK+L)=S
            S2=S2+S*S
         ENDDO
!$omp enddo
!$omp end parallel
*>>>>breakdown
         IF (SQRT(S2) .LE. 1.0D-15) THEN
            INFO=100
            RETURN
         ENDIF
!$omp parallel
         SINV=1.0D0/SQRT(S2)
!$omp do
         DO I=1,N
            EVEC(I,LOCK+L)=EVEC(I,LOCK+L)*SINV
         ENDDO
!$omp enddo
!$omp end parallel
*
      ENDDO
      FLOP=FLOP+NEW_LOCK*2.0D0*N*(KSIZE+2)
 3000 continue
*
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*     (4) create Shur vector for restart
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      IF (LOCK+NEW_LOCK .LT. NEV) THEN
         IF (IORDER.EQ.1) THEN
            IC=KSIZE-NEW_LOCK
         ELSE
            IC=IND(KSIZE-NEW_LOCK)
         ENDIF
!$omp parallel
!$omp do private(s)
         DO I=1,N
            S=0.0D0
            DO J=1,KSIZE
               S=S+Q(I,LOCK+J)*Z(J,IC)
            ENDDO
            RVEC(I)=S
         ENDDO
!$omp enddo
!$omp end parallel
         FLOP=FLOP+2.0D0*N*KSIZE
*
      END IF
*
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*     (5) deflation
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      DO L=1,NEW_LOCK
!$omp parallel do
         DO I=1,N
            Q(I,LOCK+L) = EVEC(I,LOCK+L)
         ENDDO
!$omp end parallel do
      ENDDO
      LOCK=LOCK+NEW_LOCK
*
      RETURN
      END
