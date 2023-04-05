      SUBROUTINE Xabclib_ARNOLDI( N, NNZ, IRP, ICOL, VAL,
     $                            NEV, EV, EVEC, LDE,
     $                            IATPARAM, RATPARAM,
     $                            WORK, LWORK, IWORK, LIWORK, INFO )
*
*     Xabclib_ARNOLDI : Compute eigenpairs for unsymmetric sparse matrix
*                       by Explicitly Restart Deflated ARNOLDI algorithm
*     January 2011
*
      IMPLICIT DOUBLE PRECISION ( A-H, O-Z )
*
*     .. Scalar Arguments ..
      INTEGER            INFO, LDE, LIWORK, LWORK, N, NNZ, NEV
*     ..
*     .. Array Arguments ..
      INTEGER            IRP( N+1 ), ICOL( NNZ )
      INTEGER            IWORK( LIWORK )
*
      DOUBLE PRECISION   VAL( NNZ ), WORK( LWORK )
      COMPLEX*16  EV( NEV ), EVEC( LDE, NEV )
*
      INTEGER            IATPARAM( 50 )
      DOUBLE PRECISION   RATPARAM( 50 )
*
      INTEGER            IERR_ALLOC
*
      DOUBLE PRECISION  UINF
      ALLOCATABLE ::    UINF(:)
      INTEGER           LUINF,NUM_SMP,JL
*
      EXTERNAL OMP_GET_MAX_THREADS
      INTEGER  OMP_GET_MAX_THREADS
*
*     ..
*
*  =====================================================================
*  Purpose
*  =======
*
*  Xabclib_ARNOLDI computes NEV eigenpairs of a real unsymmetric sparse 
*  matrix by Expilictly Restart Deflated ARNOLDI algorithm.
*
*
*  The Explicitly Restart Deflated ARNOLDI algorithm is based on 
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
*  (N+1)   Diagonal pointer of the matrix in CRS format. 
*
*  ICOL    (input) INTEGER array
*  (NNZ)   Column indecies of the matrix in CRS format. 
*
*  VAL     (input) DOUBLE PRECISION array 
*  (NNZ)   Value of the matrix in CRS format. 
*
*  NEV     (input) INTEGER
*          The number of desired eigenpairs.
*
*  EV(NEV) (output) COMPLEX*16 array
*          Eigenvalues,
*            if iatparam(30) = 1 : Eigenvalues are largest real-part 
*                                  order.
*            if iatparam(30) = 2 : Eigenvalues are largest magnitude
*                                  order.
*            if iatparam(30) = 3 : Eigenvalues are largest imag.-part
*                                  order.
*
*  EVEC    (output) COMPLEX*16 array, dimension (LDE, NEV)
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
*           IATPARAM(9)     : OpenATI_DURMV auto-tuned on/off.
*           IATPARAM(10)    : Fastest OpenATI_DURMV impl. method.
*           IATPARAM(11)    : Columns of Segmented Scan's algorithm.
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
*  (LWORK)   LWORK >= (5+MSIZE)*N + 5*MSIZE*MSIZE + 9*MSIZE + 6*NEV.
*
*  IWORK   (workspace/output) INTEGER array,
*  (LIWORK)  LIWORK >= MSIZE .
*
*  INFO    (output) INTEGER
*          =   0:  Successful exit
*          <   0:  If INFO = -i, the i-th argument had an illegal value
*          = 100:  Breakdown occured because of Zero vector generated.
*          = 200:  The algorithm failed to converge in DSTEVD(LAPCAK).
*          = 300:  The algorithm failed to converge in Max restart.
*          = 400:  Come to Upper Limit of CPU time
*          = 500:  Not Enough free memory.(IATPARAM(10)=12,13,21)
*
*  Local variables
*  ===============
*  IDBG      INTEGER
*            Debug print flag.
*             IATPARAM(50)!=0 : Print debug info.
*  NUM_SMP   INTEGER
*            OMP_GET_MAX_THREADS() 
*  =====================================================================
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
      NUM_SMP    = IATPARAM( 3)
      MV_AT_USYM = IATPARAM( 9)
      JL         = IATPARAM(11)
      MSIZE      = IATPARAM(27)
      IDBG       = IATPARAM(50)
*
      IF (IDBG .EQ. 1) THEN
         write(6,*) 'MV_AT_USYM,JL,MSIZE=',MV_AT_USYM,JL,MSIZE
      ENDIF
*
      IF (RATPARAM(23) .LT. 1.0D-15 ) THEN
       write(6,*) "Warning by Xabclib: Convergence criterion
     $ (RATPARAM(23)) is too small! RATPARAM(23) = ", RATPARAM(23)
      ENDIF
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
      ICK=(5+MSIZE)*N + 5*MSIZE*MSIZE + 9*MSIZE + 6*NEV
      IF (IDBG .EQ. 1) THEN
         WRITE(6,*) 'ICK=',ICK,' INPUT LWORK=',LWORK
      ENDIF
      IF (LWORK .LT. ICK ) THEN
         INFO=-13
         WRITE(6,*) '[Error] PARAMETER ERROR =',LWORK,ICK
         GOTO 100
      ENDIF
*-------------arg-14
*-------------arg-15
      ICK=MSIZE
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
         EV(1)=DCMPLX(VAL(1),0.0D0)
         EVEC( 1, 1 ) = DCMPLX(1.0D0,0.0D0)
         RETURN
      END IF
*
*
* ---------------------------------------------------------------
* -----       Matrix-vector product auto-tuning.   --------------
* ---------------------------------------------------------------
      IF (MV_AT_USYM .EQ. 2 .OR. MV_AT_USYM .EQ. 3) THEN
* ----------------- In this case AT(MV) ON.
         IF(MV_AT_USYM .EQ. 2) THEN
            IATPARAM( 9) = 2
            LUINF=INT(0.5*NUM_SMP)+1
            ALLOCATE(UINF(LUINF),STAT=IERR_ALLOC)
         ELSE
            IATPARAM( 9) = 3
            LUINF=INT(1.5*N)+INT(4.25*JL)+10
            ALLOCATE(UINF(LUINF),STAT=IERR_ALLOC)
            IF (IERR_ALLOC .NE. 0) THEN
               IATPARAM( 9) = 2
               LUINF=INT(0.5*NUM_SMP)+1
               ALLOCATE(UINF(LUINF),STAT=IERR_ALLOC)
            END IF
         END IF
*
         DO I=1,N
            WORK(I)=1.0D0
         ENDDO
         CALL OpenATI_DURMV(N,NNZ,IRP,ICOL,VAL,WORK,WORK(N+1),
     $           IATPARAM,RATPARAM,UINF,LUINF,INFO_A)
         IATPARAM( 9) = 0
      ELSE
* ----------------- In this case AT(MV) OFF.
         IATPARAM( 9) = 0
         MVCASE_USYM  = IATPARAM(10)
         IF(IDBG .EQ. 1) THEN
            write(6,*) 'MVCASE_USYM=',MVCASE_USYM
         ENDIF
         IF      (MVCASE_USYM .EQ. 11) THEN
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
         ELSE
            WRITE(6,*) '[Error] IATPARAM(10) error',MVCASE_USYM
            goto 9999
         END IF
         ALLOCATE(UINF(LUINF),STAT=IERR_ALLOC)
         IF (IERR_ALLOC .NE. 0) THEN
            INFO=500
            GO TO 9999
         END IF
         CALL OpenATI_DURMV_Setup(N,NNZ,IRP,IATPARAM,RATPARAM,
     &                            UINF,LUINF,INFO_A)
      END IF
*
      M=MSIZE
*
*     separate work
*
      LIWK=M
*
      IP1 =1
      IP2 =IP1 +N*(M+1)
      IP3 =IP2 +N
      IP4 =IP3 +N
      IP5 =IP4 +N
      IP6 =IP5 +N
      IP7 =IP6 +M*M
      IP8 =IP7 +M*M
      IP9 =IP8 +M*M
      IP10=IP9 +M*M
      IP11=IP10+M*M
      IP12=IP11+M
      IP13=IP12+M
      IP14=IP13+M
      IP15=IP14+M
      IP16=IP15+4*M
      IP17=IP16+M
      IP18=IP17+NEV
*
      CALL Xabclib_ARNOL01(N,NNZ,IRP,ICOL,VAL,
     &                     NEV,EV,EVEC,LDE,MSIZE,IATPARAM,RATPARAM,
     &                     WORK(IP1),WORK(IP2),WORK(IP3),WORK(IP4),
     &                     WORK(IP5),WORK(IP6),WORK(IP7),WORK(IP8),
     &                     WORK(IP9),WORK(IP10),WORK(IP11),WORK(IP12),
     &                     WORK(IP13),WORK(IP14),WORK(IP15),WORK(IP16),
     &                     WORK(IP17),WORK(IP18),
     &                     IWORK,LIWK,
     &                     UINF,LUINF,NUM_SMP,
     &                     IDBG,INFO)
*
*     End of Xabclib_ARNOLDI
      IF (IERR_ALLOC.EQ.0) THEN
         DEALLOCATE(UINF)
      ENDIF
*
 3000 CONTINUE
 9999 CONTINUE
      RETURN
      END
*
      SUBROUTINE Xabclib_ARNOL01(N,NNZ,IRP,ICOL,VAL,
     &                    NEV,EV,EVEC,LDE,MSIZE,IATPARAM,RATPARAM,
     &                    Q,WV1,WV2,WV3,
     &                    WV4,H,HB,VH,VHR,VHI,
     &                    HV,DD,EHR,EHI,
     &                    HWORK,RHWK,RES,SAMP,
     &                    IWORK,LIWK,
     &                    UINF,LUINF,NUM_SMP,
     &                    IDBG,INFO)
*
      IMPLICIT DOUBLE PRECISION ( A-H, O-Z )
      INTEGER            IRP( N+1 ), ICOL( NNZ )
      INTEGER            IWORK( LIWK )
*
      COMPLEX*16  EV( NEV ), EVEC( LDE, NEV )
      DOUBLE PRECISION   VAL( NNZ )
*
      DOUBLE PRECISION   Q( N,MSIZE+1 )
      DOUBLE PRECISION   WV1(N),WV2(N),WV3(N),WV4(N)
      DOUBLE PRECISION   H(MSIZE,MSIZE),HB(MSIZE,MSIZE),VH(MSIZE,MSIZE)
      DOUBLE PRECISION   VHR(MSIZE,MSIZE),VHI(MSIZE,MSIZE)
      DOUBLE PRECISION   HV(MSIZE),DD(MSIZE)
      DOUBLE PRECISION   EHR(MSIZE),EHI(MSIZE)
      DOUBLE PRECISION   HWORK(4*MSIZE),RHWK(MSIZE)
*
*
      DOUBLE PRECISION   RES(NEV),SAMP(5,NEV)
      DOUBLE PRECISION   UINF(LUINF)
*
      INTEGER            IATPARAM( 50 )
      DOUBLE PRECISION   RATPARAM( 50 )
      DOUBLE PRECISION   MAX_ETIME
*
      INTEGER            IERR_ALLOC
*
*
*  =====================================================================
*
*
      TTIME0=OMP_GET_WTIME()
*
      EPSZ      = 1.0D-15
      STOP_TOL  = RATPARAM(23)
      MAX_ITER  = IATPARAM(22)
          IF (MAX_ITER .LT. 0) MAX_ITER = N
      MAX_ETIME = RATPARAM(22)
*
      MSIZE_AT  = IATPARAM( 4)
      MSIZE_INI = IATPARAM(28)
*
      RERR =-1.0D0
      RERRD=-1.0D0
*
      FLOP=0.0D0
*+++++++++++++++++++++++++++++ Msize AT
      IF (MSIZE_AT .EQ. 1) THEN
         ISAMP=0
*------------------------OpenATI_DAFRT option
         DO J=1,NEV
         DO I=1,5
            SAMP(I,J)=0.0D0
         ENDDO
         ENDDO
*-- Tuning Point ---
         IF (MSIZE_INI .LT. NEV) THEN
            MNOW=NEV
         ELSE
            MNOW=MSIZE_INI
         END IF
*-------------------
         IF (MNOW .GT. MSIZE) MNOW=MSIZE
      ELSE
         MNOW=MSIZE
      END IF
*
*
*---------starting vector
*
      DO I=1,N-1,2
        WV1(I  )= 1.0D0
        WV1(I+1)=-1.0D0
      ENDDO
      IF (MOD(N,2).EQ.1) WV1(N)=1.0D0
*
*......................................................
*------------------------Gram-Schmidt Policy Scenario
      IF (IDBG.EQ.1) THEN
         WRITE(6,*) ' Gram Schmidt Policy=',IATPARAM(12)
      ENDIF
      NORMALFLG=0
*......................................................
*
      LOCK=0
*
      DO IR=1,MAX_ITER
*
         IF (IDBG.EQ.1) THEN
            WRITE(6,*) '============================='
            WRITE(6,*) ' >>>>> ITER= ',IR,' Num. of Subspace=',MNOW
            WRITE(6,*) '============================='
         ENDIF
*
*
         DO J=LOCK+1,MNOW
            DO I=J+2,MNOW
               H(I,J) = 0.0D0
            ENDDO
         ENDDO
*
         CALL OpenATI_DAFGS
     &           (1,N,WV1,Q,N,LOCK,DD,IATPARAM,RATPARAM,INFO)
*
!$omp parallel do default(none)
!$omp+ shared(N,Q,LOCK,WV1)
!$omp+ private(I)
         DO I=1,N
           Q(I,LOCK+1)=WV1(I)
         ENDDO
!$omp end parallel do
*
         CALL OpenATI_DURMV(N,NNZ,IRP,ICOL,VAL,Q(1,LOCK+1),WV1,
     $                         IATPARAM,RATPARAM,UINF,LUINF,INFO_A)
*
         DO K=1,LOCK+1
            S=0.0D0
!$omp parallel do default(none)
!$omp+ reduction(+:S)
!$omp+ shared(N,K,Q,WV1)
!$omp+ private(I)
            DO I=1,N
              S=S+Q(I,K)*WV1(I)
            ENDDO
!$omp end parallel do
            H(K,LOCK+1)=S
         ENDDO
*
         DO I=1,N
            S=0.0D0
            DO K=1,LOCK+1
               S=S+Q(I,K)*H(K,LOCK+1)
            ENDDO
            WV1(I)=WV1(I)-S
         ENDDO
         FLOP=FLOP+4.0D0*N*LOCK+2.0D0*NNZ+4.0D0*N*(LOCK+1)
         IATPARAM(31)=IATPARAM(31)+1
*
*------------------------------Arnoldi Decompose.
         DO M=LOCK+1,MNOW-1
            BETA=0.0D0
!$omp parallel do default(none)
!$omp+ reduction(+:BETA)
!$omp+ shared(N,WV1)
!$omp+ private(I)
            DO I=1,N
               BETA=BETA+WV1(I)*WV1(I)
            ENDDO
!$omp end parallel do
            BETA=DSQRT(BETA)
            IF (BETA .LT. EPSZ) THEN
               IF (IDBG.EQ.1) THEN
                  WRITE(6,*) ' NO MORE ARNOLDI DECOMPOSE BETA=0'  
               ENDIF
               MFINAL=M
               INFO=100
               GOTO 5000
            ENDIF
            DO I=1,N
               Q(I,M+1) = WV1(I)/BETA
            ENDDO
            H(M+1,M) = BETA
*
            CALL OpenATI_DURMV(N,NNZ,IRP,ICOL,VAL,Q(1,M+1),WV1,
     $                         IATPARAM,RATPARAM,UINF,LUINF,INFO_A)
*
            CALL OpenATI_DAFGS
     $              (0,N,WV1,Q,N,M+1,H(1,M+1),IATPARAM,RATPARAM,INFO)
*
*
         ENDDO              !M-LOOP(expand subspace loop)
         FLOP=FLOP+(MNOW-LOCK-2)*
     $             (3.0D0*N+2.0D0*NNZ+4.0D0*N*(MNOW-LOCK-2)/2)
         IATPARAM(31)=IATPARAM(31)+MNOW-LOCK-1
*
         MFINAL= M
         BETA=0.0D0
!$omp parallel do default(none)
!$omp+ reduction(+:BETA)
!$omp+ shared(N,WV1)
!$omp+ private(I)
         DO I=1,N
            BETA=BETA+WV1(I)*WV1(I)
         ENDDO
!$omp end parallel do
         BETA=DSQRT(BETA)
         IF (BETA .LT. EPSZ) THEN
            IF (IDBG.EQ.1) THEN
               WRITE(6,*) ' BETA=0 BUT MAY BE CONVERGED'
            ENDIF
            INFO=100
         ELSE
!$omp parallel do default(none)
!$omp+ shared(N,Q,MFINAL,WV1,BETA)
!$omp+ private(I)
            DO I=1,N
               Q(I,MFINAL+1) = WV1(I)/BETA
            ENDDO
!$omp end parallel do
         ENDIF
         FLOP=FLOP+3.0D0*N
 5000    CONTINUE
*
*-----------------------------------------------------------------------
*        Q(*,1:LOCK)      : CONVERGED SHUR VECTOR
*        Q(*,LOCK+1:MNOW) : NEW DECOMPOSED SUBSPACE
*        BETA             : FINAL VECTOR'S NORM
*
*    FINISH ARNOLDI DECOMPOSE
*      A * Q(*,MNOW) = Q(*,MNOW) * H(MNOW,MNOW) + BETA*Q(*,MNOW+1)*em^T
*                 em is m'th canonical vector
*-----------------------------------------------------------------------
*               (0) solve Hessenberg eigen system by lapack(DGEEV)
*               (1) residual estimate
*               (2) convergence checking
*               (3) create Ritz pair
*               (4) create Shur vector for restart
*-----------------------------------------------------------------------
         CALL Xabclib_ARN_SOLVE_H(H,HB,MSIZE,MNOW,LOCK,NEV,
     &                            EHR,EHI,VH,VHR,VHI,HWORK,IWORK,
     &                            RES,BETA,NCONV,
     &                            IATPARAM,RATPARAM,RERR,
     &                            INFO_H,IDBG)
         IF (INFO_H.NE.0) THEN
            INFO=200
            RETURN
         ENDIF
*
         IF (NCONV .EQ. 0) THEN
*-----------------------------------------------------------------------
*               (4) create new restart shur vector
*-----------------------------------------------------------------------
            ILIST=IWORK(1)
            CALL Xabclib_ARN_NEW_VEC(
     &                       N,MSIZE,MNOW,LOCK,ILIST,
     &                       Q,EHR,EHI,VH,VHR,VHI,
     &                       WV1,
     &                       IATPARAM,RATPARAM,
     &                       INFO,IDBG,FLOP)
* ------------------------WV1 IS NEW STARTING VECTOR
         ELSE
*-----------------------------------------------------------------------
*               (5) deflation one or two Ritz pair
*-----------------------------------------------------------------------
            CALL Xabclib_ARN_DEFLATION(
     &                       N,NNZ,IRP,ICOL,VAL,NEV,
     &                       MSIZE,MNOW,LOCK,NEW_LOCK,IWORK,
     &                       Q,EHR,EHI,VH,VHR,VHI,H,RHWK,DD,
     &                       WV1,WV2,WV3,WV4,
     &                       IATPARAM,RATPARAM,
     &                       UINF,LUINF,
     &                       INFO,IDBG,FLOP)
* ------------------------WV1 IS NEW STARTING VECTOR
            LOCK=LOCK+NEW_LOCK
            IF (IDBG.EQ.1) WRITE(6,*) ' AFTER Deflation LOCK=',LOCK
         END IF
*
         TTIME1=OMP_GET_WTIME()
         CPUTIME=TTIME1-TTIME0
*
         IF (INFO.NE.0) THEN
            GOTO 1000
         ENDIF
*
         IF ( MAX_ETIME.GT.0.0D0 .AND. CPUTIME.GT.MAX_ETIME ) THEN
            INFO=400
            GOTO 1000
         ENDIF
*
         IF (LOCK .GE. NEV) THEN
*-----------------------------------------------------------------------
*                 CONVERGED ALL EIGENPAIRS
*-----------------------------------------------------------------------
            GO TO 1000
         ELSE
************************************************************************
*     Automatically subspace expand with Max-Min ratio
************************************************************************
            IF (MSIZE_AT.EQ.1) THEN
               IF (ISAMP.EQ.5) ISAMP=0
               ISAMP=ISAMP+1
               DO I=LOCK+1,NEV
                  SAMP(ISAMP,I)=RES(I)
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
                  MNOW=MIN(MSIZE,MNOW+IATPARAM(5))
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
      GOTO 2000
*
*
*----------------- Arnoldi successfully ended --------------------------
 1000 CONTINUE
*------------------------------------------SOLVE ONCE MORE
      IF (IDBG.EQ.1) WRITE(6,*) '<<<< ARNOLDI FINAL SOLVE >>>>'
*
      KLOCK=0
      CALL Xabclib_ARN_SOLVE_H(H,HB,MSIZE,LOCK,KLOCK,NEV,
     &                            EHR,EHI,VH,VHR,VHI,HWORK,IWORK,
     &                            RES,BETA,NCONV,
     &                            IATPARAM,RATPARAM,RERRD,
     &                            INFO_H,IDBG)
*----------------------create Ritz-pair
      DO JJ=1,NEV
         ILIST=IWORK(JJ)
         EV(JJ)=DCMPLX(EHR(ILIST),EHI(ILIST))
      ENDDO
*
      DO JJ=1,NEV
         ILIST=IWORK(JJ)
         EI=EHI(ILIST)
         IF (EI.EQ.0.0D0) THEN
!$omp parallel do default(none)
!$omp+ shared(N,LOCK,Q,VHR,ILIST,JJ,EVEC)
!$omp+ private(I,X,K)
            DO I=1,N
               X=0.0D0
               DO K=1,LOCK
                  X=X+Q(I,K)*VHR(K,ILIST)
               ENDDO
               EVEC(I,JJ)=DCMPLX(X,0.0D0)
            ENDDO
!$omp end parallel do
            FLOP=FLOP+2.0D0*N*LOCK
         ELSE
!$omp parallel do default(none)
!$omp+ shared(N,LOCK,Q,VHR,VHI,ILIST,JJ,EVEC)
!$omp+ private(I,X,Y,K)
            DO I=1,N
               X=0.0D0
               Y=0.0D0
               DO K=1,LOCK
                  X=X+Q(I,K)*VHR(K,ILIST)
                  Y=Y+Q(I,K)*VHI(K,ILIST)
               ENDDO
               EVEC(I,JJ  )=DCMPLX(X, Y)
            ENDDO
!$omp end parallel do
            FLOP=FLOP+4.0D0*N*LOCK
         ENDIF
      ENDDO
*------------------Arnoldi finalize ------------------------------------
 2000 CONTINUE
*     CALL CPU_TIME(CTIME1)
      CTIME1=OMP_GET_WTIME()
      IF (IDBG.EQ.1) 
     &    WRITE(6,*) ' Xabclib_ARNOLDI Time=',ctime1-ttime0,'[sec]'
*-----------------------------------------------------------------------
*
      IATPARAM(23)=IR
      IATPARAM(29)=MFINAL
      RATPARAM(32)=CTIME1-TTIME0
      RATPARAM(29)=RERR
      RATPARAM(30)=1.0D-9*FLOP
*-----------------------------------------------------------------------
*
*
      RETURN
      END
*=======================================================================
      SUBROUTINE Xabclib_ARN_SOLVE_H(H,HB,MSIZE,MNOW,LOCK,NEV,
     &                            EHR,EHI,VH,VHR,VHI,HWORK,IWORK,
     &                            RES,BETA,NCONV,
     &                            IATPARAM,RATPARAM,RERR,
     &                            INFO,IDBG)
      IMPLICIT DOUBLE PRECISION ( A-H, O-Z )
*
      DOUBLE PRECISION   H(MSIZE,MSIZE),HB(MSIZE,MSIZE),VH(MSIZE,MSIZE)
      DOUBLE PRECISION   VHR(MSIZE,MSIZE),VHI(MSIZE,MSIZE)
      DOUBLE PRECISION   EHR(MSIZE),EHI(MSIZE)
      DOUBLE PRECISION   HWORK(4*MSIZE),RHWK(MSIZE)
      DOUBLE PRECISION   RES(NEV)
*
      INTEGER            IWORK( MSIZE )
*
      INTEGER            IATPARAM( 50 )
      DOUBLE PRECISION   RATPARAM( 50 )
*
      COMPLEX*16  CEIG,CV
*
*  =====================================================================
*  Purpose
*  =======
*
*  Xabclib_ARN_SOLVE_H:  (1) solve Hessenberg system
*                        (2) residual estimate 
*                        (3) convergence checking 
*                        (4) create Ritz pair
*
*  =====================================================================
      INFO=0
*
      EPSZ= 1.0D-15
*
      STOP_TOL = RATPARAM(23)
      IOPT     = IATPARAM(30)
      IDBG     = IATPARAM(50)
*
      MDIM = MNOW-LOCK
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*     (1) solve Hessenberg system
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      DO J=1,MDIM
         DO I=1,MIN(J+1,MDIM)
            HB(I,J) =H(I+LOCK,J+LOCK)
         ENDDO
         DO I=J+2,MDIM
            HB(I,J) = 0.0D0
         ENDDO
      ENDDO
*---------------------------call lapack
      LWORK=4*MSIZE
*
      MDUMMY=1
      CALL DGEEV('N','V',MDIM,HB,MSIZE,EHR,EHI,DUMMY,MDUMMY,VH,MSIZE,
     *           HWORK,LWORK,INFO)
      IF(IDBG.EQ.1) WRITE(6,*)'<<LAPACK DGEEV END>> INFO=',INFO
*
*----------------------eigenvector
      DO JJ=1,MDIM
         IF (EHI(JJ).EQ.0.0D0) THEN
            DO I=1,MDIM
               VHR(I,JJ)=VH(I,JJ)
               VHI(I,JJ)=0.0D0
            ENDDO
         ELSE IF (EHI(JJ).GT.0.0D0) THEN
            IF (JJ.EQ.MDIM) THEN
               IF (IDBG.EQ.1) WRITE(6,*) 'DGEEV ILLEAGALE EIGEN PAIR'
            END IF
            DO I=1,MDIM
               VHR(I,JJ  )= VH(I,JJ)
               VHI(I,JJ  )= VH(I,JJ+1)
               VHR(I,JJ+1)= VH(I,JJ)
               VHI(I,JJ+1)=-VH(I,JJ+1)
            ENDDO
         ENDIF
      ENDDO
*
      IF (IOPT.EQ.1) THEN
*-----------------------------------------------------------------------
*            largest real part eigenvalue
*-----------------------------------------------------------------------
         DO I=1,MDIM
            RHWK(I)=-EHR(I)
         ENDDO
      ELSE IF (IOPT.EQ.2) THEN
*-----------------------------------------------------------------------
*            largest magnitude
*-----------------------------------------------------------------------
         DO I=1,MDIM
            CEIG=DCMPLX(EHR(I),EHI(I))
            RHWK(I)=-CDABS(CEIG)
         ENDDO
      ELSE IF (IOPT.EQ.3) THEN
*-----------------------------------------------------------------------
*            largest imaginary part
*-----------------------------------------------------------------------
         DO I=1,MDIM
            RHWK(I)=-EHI(I)
         ENDDO
      ENDIF
*
      CALL Xabclib_QSORTD(IWORK,MDIM,RHWK)
*
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*      (2) residual estimate
*      (3) convergence checking
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      DO I=1,NEV
         ILIST=IWORK(I)
         IF (EHI(ILIST).LT.0.0D0) THEN
         ELSE IF (EHI(ILIST).GT.0.0D0) THEN
            IF (ILIST.EQ.MDIM) THEN
               INFO = 200
               IF (IDBG.EQ.1) WRITE(6,*) 'DGEEV ILLEAGALE EIGEN PAIR'
               GOTO 1000
            ENDIF
         ENDIF
         CEIG = DCMPLX(EHR(ILIST),EHI(ILIST))
         AAA  = CDABS(CEIG)
         IF (AAA.LE.EPSZ) AAA=1.0D0
         CV   = DCMPLX(VHR(MDIM,ILIST),VHI(MDIM,ILIST))
         RES(I)=CDABS(BETA*CV)/AAA
         IF (IDBG.EQ.1) THEN
            WRITE(6,3010) i,res(i),dreal(Ceig),dimag(Ceig)
 3010    FORMAT(1H ,'<<<resid =',i3,D12.5,
     &              ' RitzValue=',D19.12,'+i ',D19.12)
         ENDIF
      ENDDO
*
      NCONV=0
      IFIRST=0
      DO I=1,NEV
         IF (RES(I).LE.STOP_TOL) THEN
            NCONV=NCONV+1
            IF (IFIRST.EQ.0 .AND. RES(I) .GT. RERR) THEN
               RERR=RES(I)
               IFIRST=1
            ENDIF
         ELSE
            GOTO 1000
         ENDIF
      ENDDO
 1000 CONTINUE
      IF (IDBG .EQ. 1) THEN
         write(6,*) 
         write(6,*) ' SOLVE H. Num. of CONVERGED Ritz-pair=',NCONV
         write(6,*) '============'
      ENDIF
*
      RETURN
      END
*=======================================================================
      SUBROUTINE Xabclib_ARN_DEFLATION(
     &                       N,NNZ,IRP,ICOL,VAL,NEV,
     &                       MSIZE,MNOW,LOCK,NEW_LOCK,IWORK,
     &                       Q,EHR,EHI,VH,VHR,VHI,H,RHWK,DD,
     &                       WV1,WV2,WV3,WV4,
     &                       IATPARAM,RATPARAM,
     &                       UINF,LUINF,
     &                       INFO,IDBG,FLOP)
*  =====================================================================
*  Purpose
*  =======
*
*  Xabclib_ARN_DEFLATION:  deflation and create new starting vector.
*
*  =====================================================================
*
      IMPLICIT DOUBLE PRECISION ( A-H, O-Z )
      INTEGER            IRP( N+1 ), ICOL( NNZ )
      INTEGER            IWORK( MSIZE )
*
      DOUBLE PRECISION   VAL( NNZ )
*
      DOUBLE PRECISION   Q( N,MSIZE+1 )
      DOUBLE PRECISION   WV1(N),WV2(N),WV3(N),WV4(N)
      DOUBLE PRECISION   H(MSIZE,MSIZE),VH(MSIZE,MSIZE)
      DOUBLE PRECISION   VHR(MSIZE,MSIZE),VHI(MSIZE,MSIZE)
      DOUBLE PRECISION   HV(MSIZE),DD(MSIZE)
      DOUBLE PRECISION   EHR(MSIZE),EHI(MSIZE)
      DOUBLE PRECISION   RHWK(MSIZE)
*
      DOUBLE PRECISION   UINF(LUINF)
*
      INTEGER            IATPARAM( 50 )
      DOUBLE PRECISION   RATPARAM( 50 )
*
*
*
*  =====================================================================
      INFO=0
*
      MDIM=MNOW-LOCK
      ILIST=IWORK(1)
*
      EI=EHI(ILIST)
*
      IF (EI.EQ.0.0D0) THEN
*-------------------------------CONVERGED REAL RITZ-PAIR
         IF (IDBG.EQ.1) THEN
            WRITE(6,*) ' !REAL DEFLATION START',ILIST
            WRITE(6,*) ' !MDIM,MNOW,LOCK=',MDIM,MNOW,LOCK
         END IF
!$omp parallel do default(none)
!$omp+ shared(N,MDIM,Q,LOCK,VHR,ILIST,WV1)
!$omp+ private(I,X,J)
         DO I=1,N
            X=0.0D0
            DO J=1,MDIM
               X=X+Q(I,J+LOCK)*VHR(J,ILIST)
            ENDDO
            WV1(I)=X
         ENDDO
!$omp end parallel do
*-------------------------------ORTHONORMALIZE TO LOCKED VECTORS
         CALL OpenATI_DAFGS
     &           (1,N,WV1,Q,N,LOCK,DD,IATPARAM,RATPARAM,INFO)
*-------------------------------DEFLATION
         CALL OpenATI_DURMV(N,NNZ,IRP,ICOL,VAL,WV1,WV2,
     $                         IATPARAM,RATPARAM,UINF,LUINF,INFO_A)
         FLOP=FLOP+2.0D0*N*MDIM+4.0D0*N*LOCK+2.0D0*NNZ
         IATPARAM(31)=IATPARAM(31)+1
*
         DO I=1,LOCK
            X=0.0D0
!$omp parallel do default(none)
!$omp+ reduction(+:X)
!$omp+ shared(N,Q,I,WV2)
!$omp+ private(K)
            DO K=1,N
               X=X+Q(K,I)*WV2(K)
            ENDDO
!$omp end parallel do
            H(I,LOCK+1)=X
         ENDDO
         X=0.0D0
!$omp parallel do default(none)
!$omp+ reduction(+:X)
!$omp+ shared(N,WV1,WV2)
!$omp+ private(I)
         DO I=1,N
           X=X+WV1(I)*WV2(I)
         ENDDO
!$omp end parallel do
         H(LOCK+1,LOCK+1)=X
*
         DO I=LOCK+2,MSIZE
           H(I,LOCK+1) = 0.0D0
         ENDDO
*
!$omp parallel do default(none)
!$omp+ shared(N,Q,LOCK,WV1)
!$omp+ private(I)
         DO I=1,N
            Q(I,LOCK+1)=WV1(I)
         ENDDO
!$omp end parallel do
         FLOP=FLOP+2.0D0*N*LOCK+2.0D0*N
*
         NEXT = IWORK(2)
         CALL Xabclib_ARN_NEW_VEC(
     &                       N,MSIZE,MNOW,LOCK,ILIST,
     &                       Q,EHR,EHI,VH,VHR,VHI,
     &                       WV1,
     &                       IATPARAM,RATPARAM,
     &                       INFO,IDBG,FLOP)
*
         NEW_LOCK=1
*
*
      ELSE
         IF (IDBG.EQ.1) THEN
            WRITE(6,*) ' !CONJUGATE DEFLATION START',ILIST
            WRITE(6,*) ' !MDIM,MNOW,LOCK=',MDIM,MNOW,LOCK
         ENDIF
*-------------------------------CONVERGED COMPLEX-CONJGATE RITZ-PAIR
         XX=0.0D0
         YY=0.0D0
!$omp parallel do default(none)
!$omp+ reduction(+:XX)
!$omp+ reduction(+:YY)
!$omp+ shared(N,MDIM,Q,LOCK,VHR,ILIST,VHI,WV1,WV2)
!$omp+ private(I,J,X,Y)
         DO I=1,N
            X=0.0D0
            Y=0.0D0
            DO J=1,MDIM
               X=X+Q(I,J+LOCK)*VHR(J,ILIST)
               Y=Y+Q(I,J+LOCK)*VHI(J,ILIST)
            ENDDO
            WV1(I)=X
            WV2(I)=Y
            XX=XX+X*X
            YY=YY+Y*Y
         ENDDO
!$omp end parallel do
         XX=DSQRT(XX)
         YY=DSQRT(YY)
         IF (XX.NE.0.0D0) THEN
            DO I=1,N
               WV1(I)=WV1(I)/XX
            ENDDO
         END IF
         IF (YY.NE.0.0D0) THEN
            DO I=1,N
               WV2(I)=WV2(I)/YY
            ENDDO
         END IF
         FLOP=FLOP+4.0D0*N*MDIM+2.0D0*N
*
*-------------------------------ORTHONORMALIZE TO LOCKED VECTORS
         CALL OpenATI_DAFGS
     &           (1,N,WV1,Q,N,LOCK,DD,IATPARAM,RATPARAM,INFO)
         CALL OpenATI_DAFGS
     &           (0,N,WV2,Q,N,LOCK,DD,IATPARAM,RATPARAM,INFO)
*-------------------------------WV2 ORTHOGONAILZE TO WV1
         CALL OpenATI_DAFGS
     &           (1,N,WV2,WV1,N,1,DD,IATPARAM,RATPARAM,INFO)
*
*-------------------------------DEFLATION
         CALL OpenATI_DURMV(N,NNZ,IRP,ICOL,VAL,WV1,WV3,
     $                         IATPARAM,RATPARAM,UINF,LUINF,INFO_A)
         CALL OpenATI_DURMV(N,NNZ,IRP,ICOL,VAL,WV2,WV4,
     $                         IATPARAM,RATPARAM,UINF,LUINF,INFO_A)
         IATPARAM(31)=IATPARAM(31)+2
*
         DO I=1,LOCK
            X=0.0D0
            Y=0.0D0
!$omp parallel do default(none)
!$omp+ reduction(+:X)
!$omp+ reduction(+:Y)
!$omp+ shared(N,I,Q,WV3,WV4)
!$omp+ private(K)
            DO K=1,N
               X=X+Q(K,I)*WV3(K)
               Y=Y+Q(K,I)*WV4(K)
            ENDDO
!$omp end parallel do
            H(I,LOCK+1)=X
            H(I,LOCK+2)=Y
         ENDDO
         X1=0.0D0
         X2=0.0D0
         Y1=0.0D0
         Y2=0.0D0
!$omp parallel do default(none)
!$omp+ reduction(+:X1)
!$omp+ reduction(+:X2)
!$omp+ reduction(+:Y1)
!$omp+ reduction(+:Y2)
!$omp+ shared(N,WV1,WV2,WV3,WV4)
!$omp+ private(I)
         DO I=1,N
           X1=X1+WV1(I)*WV3(I)
           X2=X2+WV1(I)*WV4(I)
           Y1=Y1+WV2(I)*WV3(I)
           Y2=Y2+WV2(I)*WV4(I)
         ENDDO
!$omp end parallel do
         H(LOCK+1,LOCK+1)=X1
         H(LOCK+1,LOCK+2)=X2
         H(LOCK+2,LOCK+1)=Y1
         H(LOCK+2,LOCK+2)=Y2
*
         DO I=LOCK+3,MSIZE
           H(I,LOCK+1) = 0.0D0
           H(I,LOCK+2) = 0.0D0
         ENDDO
*
!$omp parallel do default(none)
!$omp+ shared(N,Q,LOCK,WV1,WV2)
!$omp+ private(I)
         DO I=1,N
            Q(I,LOCK+1)=WV1(I)
            Q(I,LOCK+2)=WV2(I)
         ENDDO
!$omp end parallel do
         FLOP=FLOP+8.0D0*N*LOCK+12.0D0*N+4.0D0*NNZ
*
         IF (MDIM.GE.3) THEN
            NEXT = IWORK(3)
            CALL Xabclib_ARN_NEW_VEC(
     &                       N,MSIZE,MNOW,LOCK,ILIST,
     &                       Q,EHR,EHI,VH,VHR,VHI,
     &                       WV1,
     &                       IATPARAM,RATPARAM,
     &                       INFO,IDBG,FLOP)
*
         END IF
*
         NEW_LOCK=2
*
      ENDIF
*
      RETURN
      END
*=======================================================================
      SUBROUTINE Xabclib_ARN_NEW_VEC(
     &                       N,MSIZE,MNOW,LOCK,ILIST,
     &                       Q,EHR,EHI,VH,VHR,VHI,
     &                       WV1,
     &                       IATPARAM,RATPARAM,
     &                       INFO,IDBG,FLOP)
*
      IMPLICIT DOUBLE PRECISION ( A-H, O-Z )
*
      DOUBLE PRECISION   Q( N,MSIZE+1 )
      DOUBLE PRECISION   WV1(N)
      DOUBLE PRECISION   VH(MSIZE,MSIZE)
      DOUBLE PRECISION   VHR(MSIZE,MSIZE), VHI(MSIZE,MSIZE)
      DOUBLE PRECISION   HV(MSIZE)
      DOUBLE PRECISION   EHR(MSIZE),EHI(MSIZE)
*
*
      IF (EHI(ILIST).EQ.0.0D0) THEN
!$omp parallel do default(none)
!$omp+ shared(N,MNOW,LOCK,Q,VHR,ILIST,WV1)
!$omp+ private(I,J,S)
         DO I=1,N
            S=0.0D0
            DO J=1,MNOW-LOCK
               S=S+Q(I,J+LOCK)*VHR(J,ILIST)
            ENDDO
            WV1(I)=S
         ENDDO
!$omp end parallel do
         FLOP=FLOP+(MNOW-LOCK)*N*2.0D0
      ELSE
!$omp parallel do default(none)
!$omp+ shared(N,MNOW,LOCK,Q,VHR,ILIST,VHI,WV1)
!$omp+ private(I,J,S)
         DO I=1,N
            S=0.0D0
            DO J=1,MNOW-LOCK
               S=S+Q(I,J+LOCK)*(VHR(J,ILIST)+VHI(J,ILIST))
            ENDDO
            WV1(I)=S
         ENDDO
!$omp end parallel do
         FLOP=FLOP+(MNOW-LOCK)*N*3.0D0
      ENDIF
*
      RETURN
      END
