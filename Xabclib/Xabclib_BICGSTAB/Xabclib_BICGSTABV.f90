      SUBROUTINE Xabclib_BICGSTABV (N, NNZ, IRP, ICOL, VAL, B, X,
     $                      PRECOND, NPRE, IATPARAM, RATPARAM,
     $                      WORK, LWORK, INFO)
*
*     Xabclib_BICGSTABV :
*                  Compute the solution to a linear equation system
*                  by preconditioned BICGSTAB algorithm.
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
*  Xabclib_BICGSTABV computes the solution to a linear equation system 
*  by van der Vorst.
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
*          if IPCPARM(1) = 1, Workspace for the preconditioner used.
*          if IPCPARM(1) = 2, Already defined preconditioner data.
*          On exit,
*          Defined preconditioner data.
*
*  NPRE    (input) INTEGER
*          The dimension of PRECOND.
*          if KIND_PRECOND = 5,6, NPRE >=
*
*  IATPARAM (input/output) INTEGER array, dimension (50)
*
*  RATPARAM (input/output) DOUBLE PRECISION array, dimension (50)
*
*
*  WORK    (workspace) DOUBLE PRECISION array, dimension (LWORK)
*
*  LWORK   (input) INTEGER
*          The dimension of WORK.
*          LWORK >= (8)*N + (N-1)/2+1
*
*  INFO    (output) INTEGER
*          Return code.
*          =   0:  Successful exit.
*          <   0:  if INFO = -i, the i-th argument had an illegal value.
*          = 100:  The preconditioner generation process failed.
*          = 200:  The BICG algorithm breaks down.
*          = 400:  CPU time exceeds to the limit.
*          = 500:  Iteration count exceeds to the limit.
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
      DOUBLE PRECISION  ETIME0,ETIME1,ETIME2,OMP_GET_WTIME
*
      INTEGER           ICASE,IDBG
*
      INTEGER           KIND_PRECOND
      INTEGER           IRT,IP_IDIAGP,IP_Z,IP_ZT,IP_R
      INTEGER           IP_RT,IP_P,IP_PT,IP_Q,IP_QT
      INTEGER           MV_AT_USYM,MVCASE_USYM,INFO_A
      DOUBLE PRECISION  RERR
*
      DOUBLE PRECISION  UINF
      ALLOCATABLE ::    UINF(:)
      INTEGER           LUINF,NUM_SMP,JL,IERR_ALLOC
*
      EXTERNAL OMP_GET_MAX_THREADS
      INTEGER  OMP_GET_MAX_THREADS
*=======================================================================
*
      IDBG = IATPARAM(50)
      JL   = IATPARAM(11)
*
      INFO=0
      IERR_ALLOC=-1
      IRT=0
      IATPARAM(31)=0
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
*
      NUM_SMP    = IATPARAM( 3)
*
      MV_AT_USYM = IATPARAM( 9)
*
      IF (IDBG .EQ. 1) THEN
         WRITE(6,*) 'NUM_SMP,MV_AT_USYM=',NUM_SMP,MV_AT_USYM
      END IF
* ---------------------------------------------------------------
* -----       Matrix-vector product auto-tuning.   --------------
* ---------------------------------------------------------------
      IF (MV_AT_USYM .EQ. 1) THEN
* ----------------- In this case AT(MV) ON.
         IATPARAM( 9) = 3
         LUINF=INT(1.5*N)+INT(4.25*JL)+10
         ALLOCATE(UINF(LUINF),STAT=IERR_ALLOC)
         IF (IERR_ALLOC .NE. 0) THEN
            LUINF=INT(0.5*NUM_SMP)+1
            ALLOCATE(UINF(LUINF),STAT=IERR_ALLOC)
            IATPARAM( 9) = 2
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
            GOTO 9999
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
         IF (IDBG .EQ. 1) WRITE(6,*) ' # IDIAGPGEN, INFO=',INFO
         GO TO 9999
      END IF
*
*     IOPCTPC=0 : PRECONDITIONER NOT DEFINDED YET
*     IOPCTPC=1 : PRECONDITIONER ALREADY DEFINDED
*
      IOPTPC       = IATPARAM(24)
      KIND_PRECOND = IATPARAM(25)
*
*  GENERATE PRECONDITIONER FOR BICGSTABV(M)
*
      ETIME0= OMP_GET_WTIME()
*
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
      CALL BICGSTABV  (N, NNZ, IRP, ICOL, VAL, B, X,
     $                 KIND_PRECOND,
     $                 PRECOND, NPRE, IATPARAM,RATPARAM,
     $                 WORK(IP_IDIAGP),
     $                 WORK(IP_Z),WORK(IP_ZT),WORK(IP_R),WORK(IP_RT),
     $                 WORK(IP_P),WORK(IP_PT),WORK(IP_Q),WORK(IP_QT),
     $                 UINF,LUINF,NUM_SMP,
     $                 INFO)
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
      SUBROUTINE BICGSTABV (N, NNZ, IRP, ICOL, VAL, B, X,
     $                      KIND_PRECOND,
     $                      PRECOND, NPRE, IATPARAM,RATPARAM,
     $                      IDIAGP,
     $                      R, RT, P, PT, S, ST, V, T,
     $                      UINF,LUINF,NUM_SMP,
     $                      INFO)
*
      IMPLICIT DOUBLE PRECISION ( A-H , O-Z )
*     ..
*     .. Scalar Auruments ..
      INTEGER           N, NNZ, KIND_PRECOND, NPRE, LUINF
      INTEGER           NUM_SMP, INFO
*     ..
*     .. Array Arguments
      INTEGER           IRP(N+1), ICOL(NNZ)
      DOUBLE PRECISION  VAL(NNZ), B(N), X(N), PRECOND(NPRE)
*
      INTEGER           IDIAGP(N)
      DOUBLE PRECISION  R(N), RT(N), S(N), ST(N)
      DOUBLE PRECISION  P(N), PT(N), V(N), T (N)
      DOUBLE PRECISION  UINF(LUINF)
*
      INTEGER           IATPARAM(50)
      DOUBLE PRECISION  RATPARAM(50)
*
*  =====================================================================
*  Purpose
*  =======
*
*  BICGSTABV
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
*  INFO    (output) INTEGER
*
*  =====================================================================
*     ..
*     .. Local variables
*
      INTEGER          IDBG, MV_AT, ICASE, MAX_ITER
      INTEGER          I, INFO_A ,ITER
      DOUBLE PRECISION STOP_TOL
      DOUBLE PRECISION MAX_ETIME, ETIME1, ETIME2, OMP_GET_WTIME
      DOUBLE PRECISION RNORM0, BNORM, RHO, RHO0, TD, TMPRT, TMPZ
      DOUBLE PRECISION TMPS, TMPT, RNORM, RESID
      DOUBLE PRECISION ALFA, BETA, OMEGA
*
      INTEGER          IPCPARM(10)
      DOUBLE PRECISION RPCPARM(10)
*
      EXTERNAL OMP_GET_MAX_THREADS
      INTEGER  OMP_GET_MAX_THREADS
*
      IDBG      = IATPARAM(50)
      MV_AT     = IATPARAM( 9)
      ICASE     = IATPARAM(10)
*-----------------------------------------------------------------------
*        right preconditioned ILUT-BICG 
*-----------------------------------------------------------------------
      STOP_TOL  = RATPARAM(23)
      MAX_ETIME = RATPARAM(22)
      MAX_ITER  = IATPARAM(22)
      IF (MAX_ITER .LT. 0) MAX_ITER = N
*
      ETIME1= OMP_GET_WTIME()
*
*  r = b - A x
*
      CALL OpenATI_DURMV(N,NNZ,IRP,ICOL,VAL,X,R,IATPARAM,RATPARAM,
     $                   UINF,LUINF,INFO_A)
      IATPARAM(31)=IATPARAM(31)+1
      DO I=1,N
         R(I)=B(I)-R(I)
      END DO
*
      CALL  XNORM2(N,R,rnorm0)
      CALL  XNORM2(N,B,bnorm )
      IF (IDBG .EQ. 1) THEN
         WRITE(6,*) 'INITIAL BNORM,RESID0=',BNORM,RNORM0
      ENDIF
*
      DO I=1,N
         RT(I)=R(I)
      ENDDO
*
      rho0=1.0d0
*
      DO 1000 ITER = 1, MAX_ITER
*
         rho0 = rho
         CALL XDDOT(N,RT,R,rho)
*
         IF (ITER .EQ. 1) THEN
            DO I=1,N
               P(I) = R(I)
            ENDDO
         ELSE
            beta = (rho/rho0)*(alfa/omega)
            DO I=1,N
               P(I) = R(I) + beta*(P(I)-omega*V(I))
            ENDDO
         END IF
*
* solve M(Pt) = P
*
         if (kind_precond.ge.5) then
         CALL ILUTSOL('N',N,NNZ,IRP,ICOL,VAL,PT ,P ,IDIAGP,
     $                 KIND_PRECOND,PRECOND,IATPARAM,RATPARAM,INFO)
         else
         CALL Xabclib_PCSLV(N,NNZ,IRP,ICOL,VAL,PT,P,IDIAGP,
     $           PRECOND,IATPARAM,RATPARAM,INFO)
         endif
*
         CALL OpenATI_DURMV(N,NNZ,IRP,ICOL,VAL,PT,V,IATPARAM,RATPARAM,
     $                      UINF,LUINF,INFO_A)
         IATPARAM(31)=IATPARAM(31)+1
*
         CALL XDDOT(N,RT,V,td)
*
         IF (td .EQ. 0.0D0) THEN
            INFO = 250
            CALL  XNORM2(N,RT,tmprt)
            CALL  XNORM2(N,V,tmpz)
            write(6,*) '[Error] BCGSTAB breakdown rho',rho,tmpz,tmprt
            GOTO 9000
         ENDIF
         alfa = rho / td
*
         DO I=1,N
            S(I) = R(I) - alfa*V(I)
         ENDDO
*------------------check if ||S|| is small enough?
         CALL  XNORM2(N,S,tmps)
         IF (tmps .LT. 1.0d-15) THEN
            INFO = 200
            write(6,*) '[Error] BCGSTAB small S vector',tmps
            DO I=1,N
               X(I) = X(I) + alfa*PT(I)
            ENDDO
            GOTO 9000
         ENDIF
*
* solve M(St) = S
*
         if (kind_precond.ge.5) then
         CALL ILUTSOL('N',N,NNZ,IRP,ICOL,VAL,ST ,S ,IDIAGP,
     $                 KIND_PRECOND,PRECOND,IATPARAM,RATPARAM,INFO)
         else
         CALL Xabclib_PCSLV(N,NNZ,IRP,ICOL,VAL,ST,S,IDIAGP,
     $           PRECOND,IATPARAM,RATPARAM,INFO)
         end if
*
         CALL OpenATI_DURMV(N,NNZ,IRP,ICOL,VAL,ST,T,IATPARAM,RATPARAM,
     $                      UINF,LUINF,INFO_A)
         IATPARAM(31)=IATPARAM(31)+1
*
         CALL XDDOT (N,T,S,td)
         CALL XNORMD(N,T,tmpt)
*
         IF (TMPT .EQ. 0.0D0) THEN
            INFO = 300
            WRITE(6,*) '[Error] BCGSTAB BREAKDOWN IN OMEGA',tmpt
            GOTO 9000
         ENDIF
         omega = td / tmpt
*
         DO I=1,N
            X(I) = X(I) + alfa*PT(I) + omega*ST(I)
            R(I) = S(I) - omega*T(I)
         ENDDO
*
         CALL  XNORM2(N,R,rnorm)
*
         RESID = RNORM / BNORM
*
         IF (IDBG .EQ. 1) THEN
            WRITE(6,*) ' ITER=',ITER,', RNORM=',resid
         END IF
*
         ETIME2= OMP_GET_WTIME()
         IF (MAX_ETIME .GT. 0.0D0) THEN
            IF (ETIME2-ETIME1+RATPARAM(31) .GT. MAX_ETIME) THEN
               INFO=400
               GO TO 9000
            END IF
         END IF
         IF (RESID .LE. STOP_TOL ) THEN
            GOTO 9000
         ENDIF
*
*
 1000 CONTINUE
*
 9000 CONTINUE
      IF (RESID.GT.STOP_TOL .AND. ITER.GE.MAX_ITER) THEN
         INFO=500
      END IF
      RATPARAM(28)=BNORM
      RATPARAM(29)=RNORM/RNORM0
      IATPARAM(23)=ITER
*
      RETURN
      END
