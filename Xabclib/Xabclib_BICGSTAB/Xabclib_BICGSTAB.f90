      SUBROUTINE Xabclib_BICGSTAB
     $                     (N, NNZ, IRP, ICOL, VAL, B, X,
     $                      PRECOND, NPRE, IATPARAM, RATPARAM,
     $                      WORK, LWORK, INFO)
*
*     Xabclib_BICGSTAB : 
*                     Compute the solution to a linear equation system
*                     by Itoh's right-preconditioned BICGSTAB algorithm.
*     October 2011
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
*  Xabclib_BICGSTAB computes the solution to a linear equation system 
*  by Itoh's right-preconditioned BiCGSTAB algorithm.
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
*           IATPARAM(31)    : Mat-Vec times.
*           IATPARAM(32)    : Krylov subspace iteration times.
*           IATPARAM(33)    : The flag of detectiong stagnation.
*           IATPARAM(34)    : Minimmum running iterations.
*           IATPARAM(50)    : Debug print control flag.
*
*  RATPARAM (input/output) DOUBLE PRECISION array, dimension (50)
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
*
*  WORK    (workspace) DOUBLE PRECISION array, dimension (LWORK)
*
*  LWORK   (input) INTEGER
*          The dimension of WORK.
*          LWORK >= (9)*N + (N-1)/2+1
*          
*  INFO    (output) INTEGER
*          Return code.
*          =   0:  Successful exit.
*          <   0:  if INFO = -i, the i-th argument had an illegal value.
*          = 100:  The preconditioner generation process failed.
*          = 200:  The BICG algorithm breaks down.
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
*  MAX_ITER   INTEGER
*             Maximum number of restart iterations.
*  STOP_TOL   DOUBLE PRECISION
*             Convergence criterion.
*  MAX_ETIME  DOUBLE PRECISION
*             Maximum elapsed time for the solver.
*  ICASE      INTEGER
*             Selected Algorithm NO. of matrix-vector product.
*  IDBG       INTEGER
*             Debug print control flag.
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
      DOUBLE PRECISION  UINF
      ALLOCATABLE ::    UINF(:)
      INTEGER           LUINF,NUM_SMP,JL,IERR_ALLOC
*
      INTEGER           IRT
      INTEGER           IP_IDIAGP,IP_Z,IP_ZT,IP_R,IP_RT
      INTEGER           IP_P,IP_PT,IP_Q,IP_QT,IP_R0
*
      INTEGER           MVCASE_USYM,KIND_PRECOND
      INTEGER           MV_AT_USYM,INFO_A
*
      DOUBLE PRECISION  ETIME0,RERR
*
      EXTERNAL OMP_GET_MAX_THREADS
      INTEGER  OMP_GET_MAX_THREADS
*
      IDBG = IATPARAM(50)
      JL   = IATPARAM(11)
*
*+++++ INITIALIZATION FOR RESTARTED GMRES ++++++++++++++++++++++++++++++
*
      INFO=0
      IERR_ALLOC=-1
      IRT=0
      IATPARAM(31)=0
      IATPARAM(32)=0
*--------------------------Floating Ope. counts
      RATPARAM(30)=0.0D0
*
*  CHECK INPUT PARAMETERS FOR GMRES(M)
*
      CALL Xabclib_BICGSTAB_PCHK(N,NNZ,IRP,ICOL,VAL,IATPARAM,RATPARAM,
     $                           NPRE,LWORK,INFO)
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
*
      IP_IDIAGP=1
      IP_Z     =IP_IDIAGP+ (N-1)/2+1
      IP_ZT    =IP_Z     +  N
      IP_R     =IP_ZT    +  N
      IP_RT    =IP_R     +  N
      IP_P     =IP_RT    +  N
      IP_PT    =IP_P     +  N
      IP_Q     =IP_PT    +  N
      IP_QT    =IP_Q     +  N
      IP_R0    =IP_QT    +  N
*
      NUM_SMP    = IATPARAM( 3)
*
      MV_AT_USYM = IATPARAM( 9)
*
      IF (IDBG .EQ. 1) THEN
         WRITE(6,*) 'NUM_SMP,MV_AT_USYM=',NUM_SMP,MV_AT_USYM
      ENDIF
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
         CALL OpenATI_DURMV(N,NNZ,IRP,ICOL,VAL,X,WORK(IP_Z),
     $           IATPARAM,RATPARAM,UINF,LUINF,INFO_A)
         IATPARAM( 9) = 0
      ELSE
* ----------------- In this case AT(MV) OFF.
         IATPARAM( 9) = 0
         MVCASE_USYM  = IATPARAM(10)
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
            INFO=600
            GO TO 9999
         END IF
         CALL OpenATI_DURMV_Setup(N,NNZ,IRP,IATPARAM,RATPARAM,
     &                            UINF,LUINF,INFO_A)
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
*+++ RESTARTED GMRES MAIN LOOP +++++++++++++++++++++++++++++++++++++++++
*
*
      CALL BICGSTAB  (N, NNZ, IRP, ICOL, VAL, B, X,
     $                 KIND_PRECOND,
     $                 PRECOND, NPRE, IATPARAM, RATPARAM,
     $                 WORK(IP_IDIAGP),
     $                 WORK(IP_Z),WORK(IP_ZT),WORK(IP_R),WORK(IP_RT),
     $                 WORK(IP_P),WORK(IP_PT),WORK(IP_Q),WORK(IP_QT),
     $                 WORK(IP_R0),
     $                 UINF,LUINF,NUM_SMP,
     $                 IDBG,INFO)
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
            WRITE(6,*) ' ## BICGSTAB PROCESS SUCCESSFULLY ENDED.'
         ELSE
            WRITE(6,*) ' ## BICGSTAB PROCESS FAILED.'
         END IF
         WRITE(6,6001) IATPARAM(23),RATPARAM(29),RATPARAM(32)
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
*=======================================================================
      SUBROUTINE BICGSTAB (N, NNZ, IRP, ICOL, VAL, B, X,
     $                      KIND_PRECOND,
     $                      PRECOND, NPRE, IATPARAM, RATPARAM,
     $                      IDIAGP,
     $                      R, RT, P, PT, S, ST, V, T, R0,
     $                      UINF,LUINF,NUM_SMP,
     $                      IDBG,INFO)
*
      IMPLICIT DOUBLE PRECISION ( A-H , O-Z )
*     ..
*     .. Scalar Auruments ..
      INTEGER           N, NNZ, KIND_PRECOND, NPRE, LUINF
      INTEGER           NUM_SMP, IDBG, INFO
*     ..
*     .. Array Arguments
      INTEGER           IRP(N+1), ICOL(NNZ), IDIAGP(N)
      DOUBLE PRECISION  VAL(NNZ), B(N), X(N), PRECOND(NPRE)
      DOUBLE PRECISION  R(N), RT(N), S(N), ST(N)
      DOUBLE PRECISION  P(N), PT(N), V(N), T (N), R0(N)
      DOUBLE PRECISION  UINF(LUINF)
*
      INTEGER           IATPARAM(50)
      DOUBLE PRECISION  RATPARAM(50)
*
*  =====================================================================
*  Purpose
*  =======
*
*  BICGSTAB
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
*  X       (output) DOUBLE PRECISION array, dimension (N)
*          The solution vector of Mx = y.
*
*  KIND_PRECOND (input) INTEGER
*
*
*  PRECOND (input) DOUBLE PRECISION array, dimension (NPRE)
*          Preconditioner data.
*
*  NPRE    (input) INTEGER
*          The dimension of PRECOND.
*          if KIND_PRECOND = 5,6, NPRE >=
*
*  IATPARAM (input/output) INTEGER array, dimension (50)
*
*  RATPARAM (input/output) DOUBLE PRECISION array, dimension(50)
*
*  IDIAGP  (input) INTEGER array, dimension (N)
*          The index of the diagonal elements in VAL.
*
*  R       (workspace) DOUBLE PRECISION array, dimension (N)
*
*  RT      (workspace) DOUBLE PRECISION array, dimension (N)
*
*  P       (workspace) DOUBLE PRECISION array, dimension (N)
*
*  PT      (workspace) DOUBLE PRECISION array, dimension (N)
*
*  S       (workspace) DOUBLE PRECISION array, dimension (N)
*
*  ST      (workspace) DOUBLE PRECISION array, dimension (N)
*
*  V       (workspace) DOUBLE PRECISION array, dimension (N)
*
*  T       (workspace) DOUBLE PRECISION array, dimension (N)
*
*  R0      (workspace) DOUBLE PRECISION array, dimension (N)
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
*  NUM_SMP (input) INTEGER
*          Number of threads.
*
*  IDBG    (input) INTEGER
*          Debug print control flag.
*
*  INFO    (output) INTEGER
*
*  =====================================================================
*     ..
*     .. Local variables
*
      INTEGER           MV_AT,ICASE
      INTEGER           I
      INTEGER           MAX_ITER,INFO_A,ITER
      DOUBLE PRECISION  OMP_GET_WTIME
      DOUBLE PRECISION  STOP_TOL,MAX_ETIME,ETIME1,ETIME2
      DOUBLE PRECISION  RNORM0,RHO0,BNORM,FLOP,TD,ALFA,BETA
      DOUBLE PRECISION  TMPS,TMPT,ETA,RNORM,RESID,RHO
*
      INTEGER           ISTG,ISTGCNT
      DOUBLE PRECISION  ETIME3,PRERESID,EMA,MIN_ETIME
*
      ISTGCNT   = 0
      PRERESID  = 1.0D0
      EMA       = 0.0D0
      ISTG      = IATPARAM(33)
      MPTH      = IATPARAM(6)
*
      MV_AT     = IATPARAM( 9)
      ICASE     = IATPARAM(10)
      IF (IDBG .EQ. 1) THEN
         WRITE(6,*) 'ICASE,MV_AT,NUM_SMP=',icase,MV_AT,num_smp
      ENDIF
*-----------------------------------------------------------------------
*        right preconditioned ILUT-BICG 
*-----------------------------------------------------------------------
      STOP_TOL  = RATPARAM(23)
      MAX_ETIME = RATPARAM(22)
      MAX_ITER  = IATPARAM(22)
      IF (MAX_ITER .LT. 0) MAX_ITER = N
*
      ETIME1= OMP_GET_WTIME()
      ETIME2= ETIME1
*
*  r = b - A x ,  r0 = M^(-1) r , rho0 = ||R0|| , Rt = R0 ;
*
      CALL OpenATI_DURMV(N,NNZ,IRP,ICOL,VAL,X,R,IATPARAM,RATPARAM,
     $        UINF,LUINF,INFO_A)
      IATPARAM(31)=IATPARAM(31)+1
      DO I=1,N
         R(I)=B(I)-R(I)
      END DO
*
      IF (KIND_PRECOND.GE.5) then
         CALL ILUTSOL('N',N,NNZ,IRP,ICOL,VAL,R0,R,IDIAGP,
     $                 KIND_PRECOND,PRECOND,IATPARAM,RATPARAM,INFO)
      ELSE
         CALL Xabclib_PCSLV(N,NNZ,IRP,ICOL,VAL,R0,R,IDIAGP,
     $                      PRECOND,IATPARAM,RATPARAM,INFO)
      END IF
*
      CALL  XNORM2(N,R ,RNORM0)
      CALL  XNORM2(N,R0,rho0)
      CALL  XNORM2(N,B,bnorm )
*
      IF (IDBG .EQ. 1) THEN
         WRITE(6,*) 'INITIAL RNORM0=',RNORM0
         WRITE(6,*) 'INITIAL BNORM,RESID0=',BNORM,RHO0
      ENDIF
*
      DO I=1,N
         RT(I)=R0(I)
      ENDDO
*
      FLOP=4.0D0*NNZ+6.0D0*N
*
      DO 1000 ITER = 1, MAX_ITER
*
         IF (ITER .EQ. 1) THEN
            DO I=1,N
               P(I) = RT(I)
            ENDDO
         ELSE
            DO I=1,N
               P(I) = RT(I) + beta*P(I)
            ENDDO
         END IF
*
         CALL OpenATI_DURMV(N,NNZ,IRP,ICOL,VAL,P,PT,IATPARAM,RATPARAM,
     $                      UINF,LUINF,INFO_A)
         IATPARAM(31)=IATPARAM(31)+1
         IATPARAM(32)=IATPARAM(32)+1
*
* solve M(V) =Pt
         IF (KIND_PRECOND.GE.5) then
            CALL ILUTSOL('N',N,NNZ,IRP,ICOL,VAL,V,PT,IDIAGP,
     $                    KIND_PRECOND,PRECOND,IATPARAM,RATPARAM,INFO)
         ELSE
            CALL Xabclib_PCSLV(N,NNZ,IRP,ICOL,VAL,V,PT,IDIAGP,
     $                         PRECOND,IATPARAM,RATPARAM,INFO)
         END IF
*
         CALL XDDOT(N,R0,V,td)
*
         IF (td .EQ. 0.0D0) THEN
            INFO = 250
            WRITE(6,*) 'BICGSTAB BREAKDOWN <r0,v>=',td
            GOTO 9000
         ENDIF
*
         alfa = rho0 / td
*
         DO I=1,N
            S (I) = R (I) - alfa*PT(I)
            ST(I) = RT(I) - alfa*V (I)
         ENDDO
*
*------------------check if ||ST|| is small enough?
         CALL  XNORM2(N,ST,tmps)
         IF (tmps .LT. 1.0D-15) THEN
            INFO = 200
            WRITE(6,*) 'BCGSTAB small S vector',tmps
            DO I=1,N
               X(I) = X(I) + alfa*P(I)
            ENDDO
            GOTO 9000
         ENDIF
*
         CALL OpenATI_DURMV(N,NNZ,IRP,ICOL,VAL,ST,T,IATPARAM,RATPARAM,
     $                      UINF,LUINF,INFO_A)
         IATPARAM(31)=IATPARAM(31)+1
*
         CALL XDDOT (N,T,S,td)
         CALL XNORMD(N,T,tmpt)
*
         IF (tmpt .EQ. 0.0D0) THEN
            INFO = 300
            WRITE(6,*) 'BICGSTAB BREAKDOWN in <T,T>=',tmpt
            GOTO 9000
         ENDIF
         eta = td / tmpt
*
         DO I=1,N
            X(I) = X(I) + alfa*P(I) + eta*ST(I)
            R(I) = S(I) - eta*T(I)
         ENDDO
*
         CALL  XNORM2(N,R,rnorm)
*
         RESID = RNORM / BNORM
*
         IF (IDBG .EQ. 1) THEN
            WRITE(6,*) ' ITER=',ITER,', RNORM=',RESID
         ENDIF
*
         ETIME3= ETIME2
         ETIME2= OMP_GET_WTIME()
         ETIME3= ETIME2 - ETIME3
         IF (MAX_ETIME .GT. 0.0D0) THEN
            IF (ETIME2-ETIME1+RATPARAM(31) .GT. MAX_ETIME) THEN
               INFO=400
               GO TO 9000
            END IF
         END IF
         IF (RESID .LE. STOP_TOL) THEN
            GOTO 9000
         ENDIF
*-- Detect Stagnation ---
*
        IF(ISTG.EQ.1) THEN
         CALL OpenATI_DAFSTG(ISTGCNT,EMA,RESID,PRERESID,STOP_TOL,
     $    ITER,MAX_ITER,ETIME2-ETIME1+RATPARAM(31),ETIME3,MAX_ETIME,
     $    IATPARAM,RATPARAM,INFO_A)
         IF ((ISTGCNT.GE.MPTH .AND. 
     $       ETIME2-ETIME1+RATPARAM(31) .GT. RATPARAM(34).AND.
     $       ITER .GT. IATPARAM(34)) .OR. 
     $        ISTGCNT.GT.MPTH*10) THEN
          INFO=1000
          IF (IDBG.EQ.1) WRITE(6,*)
     $    '    Stagnation of rerative residual occurred.'
          GO TO 9000
         ENDIF
         PRERESID=RESID
        ENDIF
*--------------------------
*
* solve M(Rt) = R
         IF (KIND_PRECOND.GE.5) then
            CALL ILUTSOL('N',N,NNZ,IRP,ICOL,VAL,RT,R,IDIAGP,
     $                    KIND_PRECOND,PRECOND,IATPARAM,RATPARAM,INFO)
         ELSE
            CALL Xabclib_PCSLV(N,NNZ,IRP,ICOL,VAL,RT,R,IDIAGP,
     $                         PRECOND,IATPARAM,RATPARAM,INFO)
         END IF
*
         CALL XDDOT(N,RT,R0,rho)
*
         DO I=1,N
            P(I) = P(I) - eta*V(I)
         ENDDO
*
         beta = (rho/rho0)*(alfa/eta)
         rho0 = rho
*
*
 1000 CONTINUE
*
 9000 CONTINUE
      IF (RESID.GT.STOP_TOL .AND. ITER.GE.MAX_ITER) THEN
         INFO=500
      END IF
      RATPARAM(28)=BNORM
      RATPARAM(29)=RNORM/BNORM
      IATPARAM(23)=ITER
      IATPARAM(31)=ITER*2
      IATPARAM(31)=ITER
      FLOP=FLOP+ITER*(24.0D0*N+8.0D0*NNZ)
      RATPARAM(30)=1.0D-9*(RATPARAM(30)+FLOP)
*
      RETURN
      END
