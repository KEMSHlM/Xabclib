      SUBROUTINE OpenATI_DSRMV(N,NNZ,IRP,ICOL,VAL,X,Y,
     &                         IATPARAM,RATPARAM,WK,SINF,LSINF,INFO)
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*     ..
*     ..Scalar Arguments..
      INTEGER           N,NNZ,LSINF,INFO
*     ..
*     ..Array Arguments
      INTEGER           IRP(N+1),ICOL(NNZ)
      DOUBLE PRECISION  VAL(NNZ),X(N),Y(N),SINF(LSINF)
      DOUBLE PRECISION  WK(N,*)
*
      INTEGER           IATPARAM(50)
      DOUBLE PRECISION  RATPARAM(50)
*
*  =====================================================================
*  Purpose
*  =======
*
*  OpenATI_DSRMV computes matrix-vector products for real symmetric
*  sparse matrix with formatted by Harwell-Boeing format.
*
*        Y = A * X
*
*  Arguments
*  =========
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
*  X(N)    (input) DOUBLE PRECISION array
*          Input vector.
*
*  Y(N)    (output) DOUBLE PRECISION array
*          Output vector.
*
*  IATPARAM (input/output) INTEGER array, dimension (50)
*          Integer parameter set for OpenATI and Xabclib.
*          On entry,
*           IATPARAM(3)     : Number of threads.
*           IATPARAM(7)     : OpenATI_DSRMV auto-tuned on/off.
*           IATPARAM(8)     : OpenATI_DSRMV impl. method.
*           IATPARAM(50)    : Debug print control flag.
*          On exit,
*           IATPARAM(8)     : Fastest OpenATI_DSRMV impl. method.
*                             (if IATPARAM(7) = 1)
*
*  RATPARAM (input/output) DOUBLE PRECISION array, dimension (50)
*          Double precision parameter set for OpenATI and Xabclib.
*
*  WK(N,   (Workspase) DOUBLE PRECISION array
*  NUM_SMP)Workspace for parallel vector reduction array.
*          (NUM_SMP = IATPARAM(3))
*
*  SINF    (input/output) DOUBLE PRECISION array
*  (LSINF) Setup Information for Matrix-Vector product Impl.
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
*  INFO    (output) INTEGER
*          =   0:  Successful exit
*          = 100:  Invalid IATPARAM(8) value.
*          = 200:  Invalid IATPARAM(7) value.
*
*  =====================================================================
      INTEGER KIND_PRECOND
*
      NUM_SMP      = IATPARAM( 3)
      MV_AT_SYM    = IATPARAM( 7)
      ICASE        = IATPARAM( 8)
      KIND_PRECOND = IATPARAM(25)
*
      INFO = 0
*     
      IF      ( MV_AT_SYM.EQ.0 ) THEN
         IF ( ICASE.EQ.11) THEN
*
            IF (KIND_PRECOND .EQ. 2) THEN
               CALL OpenATI_DSRMV_21(N,NNZ,IRP,ICOL,VAL,X,Y)
            ELSE
               CALL OpenATI_DSRMV_11(N,NNZ,IRP,ICOL,VAL,X,Y)
            END IF   
            
*
         ELSE IF(ICASE .EQ. 12)THEN
*
            IP_K    =1
*     
            IF (KIND_PRECOND .EQ. 2) THEN
               CALL OpenATI_DSRMV_22(N,NNZ,IRP,ICOL,VAL,X,Y,NUM_SMP,
     &              SINF(IP_K))
            ELSE
               CALL OpenATI_DSRMV_12(N,NNZ,IRP,ICOL,VAL,X,Y,NUM_SMP,
     &              SINF(IP_K))
            END IF   
*     
         ELSE IF(ICASE .EQ. 13)THEN
*
            IP_JLS  =1
            IP_JLN  =IP_JLS+(N-1)/2+1
            IP_KMB  =IP_JLN+(N-1)/2+1
            IP_KWB  =IP_KMB+NUM_SMP/2+1
*
            IF (KIND_PRECOND .EQ. 2) THEN
               CALL OpenATI_DSRMV_23(N,NNZ,IRP,ICOL,VAL,X,Y,
     &              NUM_SMP,WK,SINF(IP_JLS),SINF(IP_JLN),
     &              SINF(IP_KWB),SINF(IP_KMB))
            ELSE
               CALL OpenATI_DSRMV_13(N,NNZ,IRP,ICOL,VAL,X,Y,
     &              NUM_SMP,WK,SINF(IP_JLS),SINF(IP_JLN),
     &              SINF(IP_KWB),SINF(IP_KMB))
            END IF   
*            
         ELSE
            INFO = 100
         ENDIF
*
      ELSE IF ( MV_AT_SYM.EQ.2  .OR. MV_AT_SYM.EQ.3 ) THEN
         IDBG=IATPARAM(50)
*-----------------------------------
         L=(1.0E+7)/(4*NNZ-2*N)
         IF (L .LE. 1) THEN
            L = 2
         ELSE
            L = MIN(5,L)
         END IF
         IFAC = L
*
*------------------------------------------------------
*     CASE=11 : Equal partitioning of N by NUM_SMP
*------------------------------------------------------
         JCASE=11
         CTIME0=OMP_GET_WTIME()
         DO I=1,IFAC
         CALL OpenATI_DSRMV_11(N,NNZ,IRP,ICOL,VAL,X,Y)
         END DO
         CTIME1=OMP_GET_WTIME()
         FAST=CTIME1-CTIME0
         FAST=FAST/DFLOAT(IFAC)         
         IF (IDBG .EQ. 1) THEN
            WRITE(6,1210) JCASE,2.0D-9*(2.0D0*NNZ-1.0D0*N)/FAST
         ENDIF
*------------------------------------------------------
*     CASE=12 :  Equal partitioning of NNZ by NUM_SMP
*------------------------------------------------------
*
         KCASE=12
         IATPARAM(8)=KCASE
         IP_K    =1
*
         CALL OpenATI_DSRMV_Setup(N,NNZ,IRP,ICOL,
     &        IATPARAM,RATPARAM,SINF,LSINF,INFO)
*
         CTIME0=OMP_GET_WTIME()
         DO I=1,IFAC
         CALL OpenATI_DSRMV_12(N,NNZ,IRP,ICOL,VAL,X,Y,NUM_SMP,
     &        SINF(IP_K))
         END DO
         CTIME1=OMP_GET_WTIME()
         TTMP=CTIME1-CTIME0
         TTMP=TTMP/DFLOAT(IFAC)
         IF (IDBG .EQ. 1) THEN
            WRITE(6,1210) KCASE,2.0D-9*(2.0D0*NNZ-1.0D0*N)/TTMP
         ENDIF
         IF (TTMP.LT.FAST) THEN
            JCASE=KCASE
            FAST=TTMP
         ENDIF

         IF (MV_AT_SYM.EQ.3) THEN
*------------------------------------------------------
*     CASE=13 : Equal partitioning of NNZ by NUM_SMP 
*                    and parallel vector reduction run
*------------------------------------------------------
*     
            KCASE=13
            IATPARAM(8)=KCASE
            IP_JLS  =1
            IP_JLN  =IP_JLS+(N-1)/2+1
            IP_KMB  =IP_JLN+(N-1)/2+1
            IP_KWB  =IP_KMB+NUM_SMP/2+1
*
            CALL OpenATI_DSRMV_Setup(N,NNZ,IRP,ICOL,
     &           IATPARAM,RATPARAM,SINF,LSINF,INFO)
*
            CTIME0=OMP_GET_WTIME()
            DO I=1,IFAC
            CALL OpenATI_DSRMV_13(N,NNZ,IRP,ICOL,VAL,X,Y,
     &           NUM_SMP,WK,SINF(IP_JLS),SINF(IP_JLN),
     &           SINF(IP_KWB),SINF(IP_KMB))
            END DO
            CTIME1=OMP_GET_WTIME()
            TTMP=CTIME1-CTIME0
            TTMP=TTMP/DFLOAT(IFAC)
            IF (IDBG .EQ. 1) THEN
               WRITE(6,1210) KCASE,2.0D-9*(2.0D0*NNZ-1.0D0*N)/TTMP
            ENDIF
            IF (TTMP.LT.FAST) THEN
               JCASE=KCASE
               FAST=TTMP
            ENDIF
         ENDIF
*
         ICASE=JCASE
         IATPARAM(8)=ICASE
         IF (IDBG .EQ. 1) THEN
            WRITE(6,*) ' >>> Select Automatic Tuning Mat*Vec <<< '
            WRITE(6,1200) ICASE,2.0D-9*(2.0D0*NNZ-1.0D0*N)/FAST
 1200       FORMAT(1H ,' >>> FASTEST Mat*Vec Impl.=',I2,2X,
     &           'PERFORMANCE=',F7.3,' [Gflops]')
 1210       FORMAT(1H ,' >>> TRYAL   Mat*Vec Impl.=',I2,2X,
     &           'PERFORMANCE=',F7.3,' [Gflops]')
         ENDIF
*
         IF(ICASE .EQ. 12) THEN
            CALL OpenATI_DSRMV_Setup(N,NNZ,IRP,ICOL,
     &           IATPARAM,RATPARAM,SINF,LSINF,INFO)
         END IF
*
      ELSE
         INFO = 200
      ENDIF
*
      RETURN
      END
*
*
      SUBROUTINE OpenATI_DSRMV_11(N,NNZ,IRP,ICOL,VAL,X,Y)
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*     ..
*     ..Scalar Arguments..
      INTEGER           N,NNZ
*     ..
*     ..Array Arguments 
      INTEGER           IRP(N+1),ICOL(NNZ)  
      DOUBLE PRECISION  VAL(NNZ),X(N),Y(N)
*
*  =====================================================================
*  Purpose
*  =======
*
*  OpenATI_DRSMV_11 computes matrix-vector products for real symmetric
*  sparse matrix with equal partitioning of N by NUM_SMP.
*
*        Y = A * X
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
*  X(N)    (input) DOUBLE PRECISION array
*          Input vector.
*
*  Y(N)    (output) DOUBLE PRECISION array
*          Output vector.
*
*================================================================
*
!$omp parallel do private(S,JJ,JC,I)
      DO I=1,N
         S=0.0D0
         DO JC=IRP(I),IRP(I+1)-1
            JJ=ICOL(JC)
            S=S+VAL(JC)*X(JJ)
         ENDDO
         Y(I)=S
      ENDDO
!$omp end parallel do
*
      DO I=1,N
         DO JC=IRP(I)+1,IRP(I+1)-1
            JJ=ICOL(JC)
            Y(JJ)=Y(JJ)+VAL(JC)*X(I)
         ENDDO
      ENDDO
*
      RETURN
      END
*
*
      SUBROUTINE OpenATI_DSRMV_12(N,NNZ,IRP,ICOL,VAL,X,Y,NUM_SMP,
     &     KMBORDER)
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*     ..
*     .. Scalar Arguments..
      INTEGER           N,NNZ,NUM_SMP
*     ..
*     ..Array Arguments
      INTEGER           IRP(N+1),ICOL(NNZ)
      INTEGER           KMBORDER(0:NUM_SMP)
      DOUBLE PRECISION  VAL(NNZ),X(N),Y(N)
*
*  =====================================================================
*  Purpose
*  =======
*
*  OpenATI_DRSMV_12 computes matrix-vector products for real symmetric
*  sparse matrix with equal partitioning of NNZ by NUM_SMP.
*
*        Y = A * X
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
*  X(N)    (input) DOUBLE PRECISION array
*          Input vector.
*
*  Y(N)    (output) DOUBLE PRECISION array
*          Output vector.
*
*  NUM_SMP (input) INTEGER
*          Number of threads 
*
*  KMBORDER(input) INTEGER array
*  (NUM    Partition infomation of matrix   
*  _SMP+1)
*  
*=====================================================================
*
!$omp parallel do private(S,K,I,JJ,JC)
         DO K=1,NUM_SMP
            DO I=KMBORDER(K-1)+1,KMBORDER(K)
               S=0.0D0
               DO JC=IRP(I),IRP(I+1)-1
                  JJ=ICOL(JC)
                  S=S+VAL(JC)*X(JJ)
               END DO
               Y(I)=S
            END DO
         END DO
!$omp end parallel do
*
         DO I=1,N
            DO JC=IRP(I)+1,IRP(I+1)-1
               JJ=ICOL(JC)
               Y(JJ)=Y(JJ)+VAL(JC)*X(I)
            END DO
         END DO
*
         RETURN
         END
*
*
      SUBROUTINE OpenATI_DSRMV_13(N,NNZ,IRP,ICOL,VAL,X,Y,NUM_SMP,
     &     WK,JLS,JLN,KWBORDER,KMBORDER)
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*     ..
*     ::Scalar Arguments..
      INTEGER           N,NNZ,NUM_SMP
*     ..
*     ..Array Arguments
      INTEGER           IRP(N+1),ICOL(NNZ)
      INTEGER           JLS(N),JLN(N)
      INTEGER           KWBORDER(0:NUM_SMP),KMBORDER(0:NUM_SMP)
      DOUBLE PRECISION  VAL(NNZ),X(N),Y(N)
*
      DOUBLE PRECISION  WK(N,NUM_SMP),AA
*
*  =====================================================================
*  Purpose
*  =======
*
*  OpenATI_DRSMV_13 computes matrix-vector products for real symmetric
*  sparse matrix with equal partitioning of NNZ by NUM_SMP and parallel
*  vector reduction run.
*
*        Y = A * X
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
*  X(N)    (input) DOUBLE PRECISION array
*          Input vector.
*
*  Y(N)    (output) DOUBLE PRECISION array
*          Output vector.
*
*  NUM_SMP (input) INTEGER
*          Number of threads 
*
*  WK(N,   (Workspase) DOUBLE PRECISION array
*  NUM_SMP)Workspace for parallel vector reduction array
*
*  JLS     (input) INTEGER array
*  (N)     Start point of vector reduction part
*
*  JLN     (input) INTEGER array
*  (N)     End point of vector reduction part
*
*  KMBORDER(input) INTEGER array
*  (NUM    Partition infomation of matrix
*   _SMP+1)
*
*  KWBORDER(input) INTEGER array
*  (NUM    Partition infomation in vector reduction array
*   _SMP+1)
*
*=====================================================================
*
!$omp parallel do private(S,XDIAG,AA,JJ,I,JC) 
         DO K=1,NUM_SMP
            DO I=KMBORDER(K-1)+1,N
              WK(I,K)=0.0D0
            ENDDO
*
            DO I=KMBORDER(K-1)+1,KMBORDER(K)
               XDIAG=X(I)
               S=VAL(IRP(I))*XDIAG
               DO JC=IRP(I)+1,IRP(I+1)-1
                  JJ=ICOL(JC)
                  AA=VAL(JC)
                  S=S+AA*X(JJ)
                  WK(JJ,K)=WK(JJ,K)+AA*XDIAG
               ENDDO
               Y(I)=S
            ENDDO
         ENDDO
!$omp end parallel do
*
!$omp parallel do private(S)
      DO K=1,NUM_SMP
         DO I=KWBORDER(K-1)+1,KWBORDER(K)
            S=0.0D0
            DO J=JLS(I),JLN(I)
               S=S+WK(I,J)
            ENDDO
            Y(I)=Y(I)+S
         ENDDO
      END DO
!$omp end parallel do
*
      RETURN
      END
*
      SUBROUTINE OpenATI_DSRMV_21(N,NNZ,IRP,ICOL,VAL,X,Y)
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*     ..
*     ..Scalar Arguments..
      INTEGER           N,NNZ
*     ..
*     ..Array Arguments
      INTEGER           IRP(N+1),ICOL(NNZ)
      DOUBLE PRECISION  VAL(NNZ),X(N),Y(N)
*
*  =====================================================================
*  Purpose
*  =======
*
*  OpenATI_DRSMV_21 computes matrix-vector products for real symmetric
*  sparse matrix with equal partitioning of N by NUM_SMP.
*
*        Y = A * X
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
*  X(N)    (input) DOUBLE PRECISION array
*          Input vector.
*
*  Y(N)    (output) DOUBLE PRECISION array
*          Output vector.
*
*================================================================
*
!$omp parallel do private(S,JJ,JC,I)
      DO I=1,N
         S=X(I)
         DO JC=IRP(I)+1,IRP(I+1)-1
            JJ=ICOL(JC)
            S=S+VAL(JC)*X(JJ)
         ENDDO
         Y(I)=S
      ENDDO
!$omp end parallel do
*
      DO I=1,N
         DO JC=IRP(I)+1,IRP(I+1)-1
            JJ=ICOL(JC)
            Y(JJ)=Y(JJ)+VAL(JC)*X(I)
         ENDDO
      ENDDO
*
      RETURN
      END
*
*
      SUBROUTINE OpenATI_DSRMV_22(N,NNZ,IRP,ICOL,VAL,X,Y,NUM_SMP,
     &     KMBORDER)
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*     ..
*     .. Scalar Arguments..
      INTEGER           N,NNZ,NUM_SMP
*     ..
*     ..Array Arguments
      INTEGER           IRP(N+1),ICOL(NNZ)
      INTEGER           KMBORDER(0:NUM_SMP)
      DOUBLE PRECISION  VAL(NNZ),X(N),Y(N)
*
*  =====================================================================
*  Purpose
*  =======
*
*  OpenATI_DRSMV_22 computes matrix-vector products for real symmetric
*  sparse matrix with equal partitioning of NNZ by NUM_SMP.
*
*        Y = A * X
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
*  X(N)    (input) DOUBLE PRECISION array
*          Input vector.
*
*  Y(N)    (output) DOUBLE PRECISION array
*          Output vector.
*
*  NUM_SMP (input) INTEGER
*          Number of threads
*
*  KMBORDER(input) INTEGER array
*  (NUM    Partition infomation of matrix
*  _SMP+1)
*
*=====================================================================
*
!$omp parallel do private(S,K,I,JJ,JC)
         DO K=1,NUM_SMP
            DO I=KMBORDER(K-1)+1,KMBORDER(K)
               S=X(I)
               DO JC=IRP(I)+1,IRP(I+1)-1
                  JJ=ICOL(JC)
                  S=S+VAL(JC)*X(JJ)
               END DO
               Y(I)=S
            END DO
         END DO
!$omp end parallel do
*
         DO I=1,N
            DO JC=IRP(I)+1,IRP(I+1)-1
               JJ=ICOL(JC)
               Y(JJ)=Y(JJ)+VAL(JC)*X(I)
            END DO
         END DO
*
         RETURN
         END
*
*
      SUBROUTINE OpenATI_DSRMV_23(N,NNZ,IRP,ICOL,VAL,X,Y,NUM_SMP,
     &     WK,JLS,JLN,KWBORDER,KMBORDER)
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*     ..
*     ::Scalar Arguments..
      INTEGER           N,NNZ,NUM_SMP
*     ..
*     ..Array Arguments
      INTEGER           IRP(N+1),ICOL(NNZ)
      INTEGER           JLS(N),JLN(N)
      INTEGER           KWBORDER(0:NUM_SMP),KMBORDER(0:NUM_SMP)
      DOUBLE PRECISION  VAL(NNZ),X(N),Y(N)
*
      DOUBLE PRECISION  WK(N,NUM_SMP),AA
*
*  =====================================================================
*  Purpose
*  =======
*
*  OpenATI_DRSMV_23 computes matrix-vector products for real symmetric
*  sparse matrix with equal partitioning of NNZ by NUM_SMP and parallel
*  vector reduction run.
*
*        Y = A * X
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
*  X(N)    (input) DOUBLE PRECISION array
*          Input vector.
*
*  Y(N)    (output) DOUBLE PRECISION array
*          Output vector.
*
*  NUM_SMP (input) INTEGER
*          Number of threads
*
*  WK(N,   (Workspase) DOUBLE PRECISION array
*  NUM_SMP)Workspace for parallel vector reduction array
*
*  JLS     (input) INTEGER array
*  (N)     Start point of vector reduction part
*
*  JLN     (input) INTEGER array
*  (N)     End point of vector reduction part
*
*  KMBORDER(input) INTEGER array
*  (NUM    Partition infomation of matrix
*   _SMP+1)
*
*  KWBORDER(input) INTEGER array
*  (NUM    Partition infomation in vector reduction array
*   _SMP+1)
*
*=====================================================================
*
!$omp parallel do private(S,XDIAG,AA,JJ,I,JC)
         DO K=1,NUM_SMP
            DO I=KMBORDER(K-1)+1,N
              WK(I,K)=0.0D0
            ENDDO
*
            DO I=KMBORDER(K-1)+1,KMBORDER(K)
               XDIAG=X(I)
               S=XDIAG
               DO JC=IRP(I)+1,IRP(I+1)-1
                  JJ=ICOL(JC)
                  AA=VAL(JC)
                  S=S+AA*X(JJ)
                  WK(JJ,K)=WK(JJ,K)+AA*XDIAG
               ENDDO
               Y(I)=S
            ENDDO
         ENDDO
!$omp end parallel do
*
!$omp parallel do private(S)
      DO K=1,NUM_SMP
         DO I=KWBORDER(K-1)+1,KWBORDER(K)
            S=0.0D0
            DO J=JLS(I),JLN(I)
               S=S+WK(I,J)
            ENDDO
            Y(I)=Y(I)+S
         ENDDO
      END DO
!$omp end parallel do
*
      RETURN
      END
*
      SUBROUTINE OpenATI_DSRMV_31(N,NNZ,IRP,ICOL,VAL,X,Y)
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*     ..
*     ..Scalar Arguments..
      INTEGER           N,NNZ
*     ..
*     ..Array Arguments
      INTEGER           IRP(N+1),ICOL(NNZ)
      DOUBLE PRECISION  VAL(NNZ),X(N),Y(N)
*
*  =====================================================================
*  Purpose
*  =======
*
*  OpenATI_DRSMV_31 computes matrix-vector products for real symmetric
*  sparse matrix with equal partitioning of N by NUM_SMP.
*
*        Y = A * X
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
*  X(N)    (input) DOUBLE PRECISION array
*          Input vector.
*
*  Y(N)    (output) DOUBLE PRECISION array
*          Output vector.
*
*================================================================
*
!$omp parallel do private(S,JJ,JC,I)
      DO I=1,N
         S=0.0D0
         DO JC=IRP(I)+1,IRP(I+1)-1
            JJ=ICOL(JC)
            S=S+VAL(JC)*X(JJ)
         ENDDO
         Y(I)=S
      ENDDO
!$omp end parallel do
*
      DO I=1,N
         DO JC=IRP(I)+1,IRP(I+1)-1
            JJ=ICOL(JC)
            Y(JJ)=Y(JJ)+VAL(JC)*X(I)
         ENDDO
      ENDDO
*
      RETURN
      END
*
*
      SUBROUTINE OpenATI_DSRMV_32(N,NNZ,IRP,ICOL,VAL,X,Y,NUM_SMP,
     &     KMBORDER)
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*     ..
*     .. Scalar Arguments..
      INTEGER           N,NNZ,NUM_SMP
*     ..
*     ..Array Arguments
      INTEGER           IRP(N+1),ICOL(NNZ)
      INTEGER           KMBORDER(0:NUM_SMP)
      DOUBLE PRECISION  VAL(NNZ),X(N),Y(N)
*
*  =====================================================================
*  Purpose
*  =======
*
*  OpenATI_DRSMV_32 computes matrix-vector products for real symmetric
*  sparse matrix with equal partitioning of NNZ by NUM_SMP.
*
*        Y = A * X
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
*  X(N)    (input) DOUBLE PRECISION array
*          Input vector.
*
*  Y(N)    (output) DOUBLE PRECISION array
*          Output vector.
*
*  NUM_SMP (input) INTEGER
*          Number of threads
*
*  KMBORDER(input) INTEGER array
*  (NUM    Partition infomation of matrix
*  _SMP+1)
*
*=====================================================================
*
!$omp parallel do private(S,K,I,JJ,JC)
         DO K=1,NUM_SMP
            DO I=KMBORDER(K-1)+1,KMBORDER(K)
               S=0.0D0
               DO JC=IRP(I)+1,IRP(I+1)-1
                  JJ=ICOL(JC)
                  S=S+VAL(JC)*X(JJ)
               END DO
               Y(I)=S
            END DO
         END DO
!$omp end parallel do
*
         DO I=1,N
            DO JC=IRP(I)+1,IRP(I+1)-1
               JJ=ICOL(JC)
               Y(JJ)=Y(JJ)+VAL(JC)*X(I)
            END DO
         END DO
*
         RETURN
         END
*
*
      SUBROUTINE OpenATI_DSRMV_33(N,NNZ,IRP,ICOL,VAL,X,Y,NUM_SMP,
     &     WK,JLS,JLN,KWBORDER,KMBORDER)
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*     ..
*     ::Scalar Arguments..
      INTEGER           N,NNZ,NUM_SMP
*     ..
*     ..Array Arguments
      INTEGER           IRP(N+1),ICOL(NNZ)
      INTEGER           JLS(N),JLN(N)
      INTEGER           KWBORDER(0:NUM_SMP),KMBORDER(0:NUM_SMP)
      DOUBLE PRECISION  VAL(NNZ),X(N),Y(N)
*
      DOUBLE PRECISION  WK(N,NUM_SMP),AA
*
*  =====================================================================
*  Purpose
*  =======
*
*  OpenATI_DRSMV_33 computes matrix-vector products for real symmetric
*  sparse matrix with equal partitioning of NNZ by NUM_SMP and parallel
*  vector reduction run.
*
*        Y = A * X
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
*  X(N)    (input) DOUBLE PRECISION array
*          Input vector.
*
*  Y(N)    (output) DOUBLE PRECISION array
*          Output vector.
*
*  NUM_SMP (input) INTEGER
*          Number of threads
*
*  WK(N,   (Workspase) DOUBLE PRECISION array
*  NUM_SMP)Workspace for parallel vector reduction array
*
*  JLS     (input) INTEGER array
*  (N)     Start point of vector reduction part
*
*  JLN     (input) INTEGER array
*  (N)     End point of vector reduction part
*
*  KMBORDER(input) INTEGER array
*  (NUM    Partition infomation of matrix
*   _SMP+1)
*
*  KWBORDER(input) INTEGER array
*  (NUM    Partition infomation in vector reduction array
*   _SMP+1)
*
*=====================================================================
*
!$omp parallel do private(S,XDIAG,AA,JJ,I,JC)
         DO K=1,NUM_SMP
            DO I=KMBORDER(K-1)+1,N
              WK(I,K)=0.0D0
            ENDDO
*
            DO I=KMBORDER(K-1)+1,KMBORDER(K)
               XDIAG=X(I)
               S=0.0D0
               DO JC=IRP(I)+1,IRP(I+1)-1
                  JJ=ICOL(JC)
                  AA=VAL(JC)
                  S=S+AA*X(JJ)
                  WK(JJ,K)=WK(JJ,K)+AA*XDIAG
               ENDDO
               Y(I)=S
            ENDDO
         ENDDO
!$omp end parallel do
*
!$omp parallel do private(S)
      DO K=1,NUM_SMP
         DO I=KWBORDER(K-1)+1,KWBORDER(K)
            S=0.0D0
            DO J=JLS(I),JLN(I)
               S=S+WK(I,J)
            ENDDO
            Y(I)=Y(I)+S
         ENDDO
      END DO
!$omp end parallel do
*
      RETURN
      END
