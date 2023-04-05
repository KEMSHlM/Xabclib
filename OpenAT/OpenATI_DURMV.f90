      SUBROUTINE OpenATI_DURMV(N,NNZ,IRP,ICOL,VAL,X,Y,IATPARAM,RATPARAM,
     &                         UINF,LUINF,INFO)
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
      INTEGER           JL

*     ..
*     ..Scalar Arguments..
      INTEGER           N,NNZ,ICASE,LUINF,NUM_SMP
*     ..
*     ..Array Arguments
      INTEGER           IRP(N+1),ICOL(NNZ)
      DOUBLE PRECISION  VAL(NNZ),X(N),Y(N),UINF(LUINF)
*
      INTEGER           IATPARAM(50)
      DOUBLE PRECISION  RATPARAM(50)
*
      EXTERNAL         OMP_GET_WTIME
      DOUBLE PRECISION OMP_GET_WTIME
*  =====================================================================
*  Purpose
*  =======
*
*  OpenATI_DURMV computes matrix-vector products for real unsymmetric
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
*  (N+1)   Diagonal pointer of the matrix in Harwell-Boeing format.
*
*  ICOL    (input) INTEGER array
*  (NNZ)   Column indecies of the matrix in Harwell-Boeing  format.
*
*  VAL     (input) DOUBLE PRECISION array
*  (NNZ)   Value of the matrix in Harwell-Boeing format.
*
*  X(N)    (input) DOUBLE PRECISION array
*          Input vector.
*
*  Y(N)    (output) DOUBLE PRECISION array
*          Output vector.
*
*  IATPARAM (input) INTEGER array, dimension (50)
*          Integer parameter set for OpenATI and Xabclib.
*           IATPARAM(3)     : Number of threads.
*           IATPARAM(9)     : OpenATI_DURMV auto-tuned on/off.
*           IATPARAM(10)    : OpenATI_DURMV impl. method.
*           IATPARAM(11)    : Number of Segment vector.
*
*  RATPARAM (input) DOUBLE PRECISION array, dimension (50)
*          Double precision parameter set for OpenATI and Xabclib.
*
*  UINF    (input/output) DOUBLE PRECISION array
*          Setup Information for Matrix-Vector product Impl.
*
*  LUINF   (input) INTEGER
*          The size of UINF.
*          case IATPARAM(9)=0 or 1 :
*             if IATPARAM(10) = 11 :
*                LUINF >= 0
*             else if IATPARAM(10) = 12 :
*                LUINF >= INT(0.5*NUM_SMP)+1
*             else if IATPARAM(10) = 13 :
*                LUINF >= INT(1.5*N)+INT(4.25*JL)+10
*             else if IATPARAM(10) = 21 :
*                LUINF >= INT(1.125*NNZ)+INT(2.125*JL)+10
*          case IATPARAM(9)=2 :
*             LUINF >= INT(0.5*NUM_SMP)+1
*          case IATPARAM(9)=3 or 4 :
*             LUINF >= INT(1.5*N)+INT(4.25*JL)+10
*          (NUM_SMP=IATPARAM(3), JL=IATPARAM(11))
*
*  INFO    (output) INTEGER
*          =   0:  Successful exit
*          = 100:  Invalid IATPARAM(10) value.
*          = 200:  Invalid IATPARAM(9) value.
*
*  =====================================================================
      NUM_SMP    = IATPARAM( 3)
      MV_AT_USYM = IATPARAM( 9)
      ICASE      = IATPARAM(10)
*------------------------------------------------------
*     IATPARAM(11)(JL) Automatic Set-up
*------------------------------------------------------
      IF (IATPARAM(9) .EQ. 1 .OR. IATPARAM(9) .EQ. 4) THEN
        IF (IATPARAM(11) .GT. NUM_SMP) THEN
           IATPARAM(11) = IATPARAM(11) - MOD(IATPARAM(11),NUM_SMP)
        END IF
      END IF
      JL         = IATPARAM(11)
*
      INFO = 0
*     
      IF      ( (MV_AT_USYM .EQ. 0) .OR. (MV_AT_USYM .EQ. 1) ) THEN
         IF ( ICASE .EQ.11 ) THEN
            CALL OpenATI_DURMV_11(N,NNZ,IRP,ICOL,VAL,X,Y)
         ELSE IF( ICASE .EQ. 12)THEN
*
            IP_K=1
*
            CALL OpenATI_DURMV_12(N,NNZ,IRP,ICOL,VAL,X,Y,
     &           NUM_SMP,UINF(IP_K))
*
         ELSE IF( ICASE .EQ. 13)THEN
*     
            IP_VSS   =1
            IP_S     =IP_VSS+N+JL
            IP_MF    =IP_S+JL
            IP_JFS   =IP_MF+(N+JL-1)/2+1
            IP_JY    =IP_JFS+JL/2+1
            IP_JS    =IP_JY+(JL-1)/2+1
            IP_JSF   =IP_JS+(JL-1)/2+1
            IP_PRE   =IP_JSF+(JL-1)/8+1
*
            IS=(NNZ-1)/JL+1
            CALL OpenATI_DURMV_13(N,NNZ,ICOL,VAL,X,Y,
     &           IS,UINF(IP_VSS),UINF(IP_S),UINF(IP_MF),UINF(IP_JFS),
     &           UINF(IP_JY),UINF(IP_JS),UINF(IP_JSF),UINF(IP_PRE),JL)      
*
         ELSE IF( ICASE .EQ. 21) THEN
*
            IS=(NNZ-1)/JL+1
*
            IP_VS    =1
            IP_S     =IP_VS+NNZ
            IP_JLA   =IP_S+JL
            IP_JST   =IP_JLA+(JL-1)/2+1
            IP_F     =IP_JST+(JL-1)/2+1
            IP_PRE   =IP_F+(NNZ-1)/8+1
*
            CALL OpenATI_DURMV_21(N,NNZ,ICOL,IRP,IS,VAL,X,Y,
     &           UINF(IP_VS),UINF(IP_S),UINF(IP_JLA),UINF(IP_JST),
     &           UINF(IP_F),UINF(IP_PRE),JL)
*
         ELSE
            INFO = 100
         ENDIF
*     
      ELSE IF ( MV_AT_USYM.EQ.2  .OR. MV_AT_USYM.EQ.3 
     &          .OR. MV_AT_USYM.EQ.4) THEN
         IDBG = IATPARAM(50)
*--------------------
         L=(1.0E+7)/(2*NNZ)
         IF (L .LE. 1) THEN
            L = 2
         ELSE
            L = MIN(5,L)
         END IF
         IFAC = L
*------------------------------------------------------
*     CASE=11 : Equal partitioning of N by NUM_SMP
*------------------------------------------------------
         JCASE=11
         CALL OpenATI_DURMV_11(N,NNZ,IRP,ICOL,VAL,X,Y)
         CTIME0=OMP_GET_WTIME()
         DO I=1,IFAC
            CALL OpenATI_DURMV_11(N,NNZ,IRP,ICOL,VAL,X,Y)
         END DO
         CTIME1=OMP_GET_WTIME()
         FAST=CTIME1-CTIME0
         FAST=FAST/DFLOAT(IFAC)
         IF (IDBG .EQ. 1) THEN
            WRITE(6,1210) JCASE,1.0D-9*(2.0D0*NNZ)/FAST
         ENDIF
*------------------------------------------------------
*     CASE=12 : Equal partitioning of NNZ by NUM_SMP
*------------------------------------------------------
*
         KCASE=12
         IATPARAM(10) = KCASE
         IP_K=1
*
         CALL OpenATI_DURMV_Setup(N,NNZ,IRP,IATPARAM,RATPARAM,
     &        UINF,LUINF,INFO)
*
         CTIME0=OMP_GET_WTIME()
         DO I=1,IFAC
            CALL OpenATI_DURMV_12(N,NNZ,IRP,ICOL,VAL,X,Y,
     &        NUM_SMP,UINF(IP_K))
         END DO
         CTIME1=OMP_GET_WTIME()
         TTMP=CTIME1-CTIME0
         TTMP=TTMP/DFLOAT(IFAC)
         IF (IDBG .EQ. 1) THEN
            WRITE(6,1210) KCASE,1.0D-9*(2.0D0*NNZ)/TTMP
         ENDIF
         IF (TTMP.LT.FAST) THEN
            JCASE=KCASE
            FAST=TTMP
         ENDIF
*------------------------------------------------------
*     CASE=13 : Improved segment scan
*------------------------------------------------------
*
         IF (MV_AT_USYM.EQ.3 .OR. MV_AT_USYM.EQ.4) THEN
            KCASE=13
            IATPARAM(10) = KCASE
            IP_VSS   =1
            IP_S     =IP_VSS+N+JL
            IP_MF    =IP_S+JL
            IP_JFS   =IP_MF+(N+JL-1)/2+1
            IP_JY    =IP_JFS+JL/2+1
            IP_JS    =IP_JY+(JL-1)/2+1
            IP_JSF   =IP_JS+(JL-1)/2+1
            IP_PRE   =IP_JSF+(JL-1)/8+1
*
            CALL OpenATI_DURMV_Setup(N,NNZ,IRP,IATPARAM,RATPARAM,
     &           UINF,LUINF,INFO)
*
            IS=(NNZ-1)/JL+1
            CTIME0=OMP_GET_WTIME()
            DO I=1,IFAC
               CALL OpenATI_DURMV_13(N,NNZ,ICOL,VAL,X,Y,
     &        IS,UINF(IP_VSS),UINF(IP_S),UINF(IP_MF),UINF(IP_JFS),
     &        UINF(IP_JY),UINF(IP_JS),UINF(IP_JSF),UINF(IP_PRE),JL)
            END DO
            CTIME1=OMP_GET_WTIME()
            TTMP=CTIME1-CTIME0
            TTMP=TTMP/DFLOAT(IFAC)
            IF (IDBG .EQ. 1) THEN
               WRITE(6,1210) KCASE,1.0D-9*(2.0D0*NNZ)/TTMP
            ENDIF
            IF (TTMP.LT.FAST) THEN
               JCASE=KCASE
               FAST=TTMP
            ENDIF
         ENDIF
*     
         ICASE=JCASE
         IATPARAM(10) = ICASE
         IF (IDBG .EQ. 1) THEN
            WRITE(6,*) ' >>> Select Automatic Tuning Mat*Vec <<< '
            WRITE(6,1200) ICASE,1.0D-9*(2.0D0*NNZ)/FAST
 1200       FORMAT(1H ,' >>> FASTEST Mat*Vec Impl.=',I2,2X,
     &           'PERFORMANCE=',F7.3,' [Gflops]')
 1210       FORMAT(1H ,' >>> TRYAL   Mat*Vec Impl.=',I2,2X,
     &           'PERFORMANCE=',F7.3,' [Gflops]')
         ENDIF
*
         IF (ICASE .EQ. 12) THEN
            CALL OpenATI_DURMV_Setup(N,NNZ,IRP,IATPARAM,RATPARAM,
     &           UINF,LUINF,INFO)
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
      SUBROUTINE OpenATI_DURMV_11(N,NNZ,IRP,ICOL,VAL,X,Y)
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*     ..
*     ..Scalar Argument..
      INTEGER           N,NNZ
*     ..
*     ..Array Argument
      INTEGER           IRP(N+1),ICOL(NNZ)
      DOUBLE PRECISION  X(N),Y(N),VAL(NNZ)
*      
      INTEGER           I,J_PTR
      DOUBLE PRECISION  S
*
*  =====================================================================
*  Purpose
*  =======
*
*  OpenATI_DURMV_11 computes matrix-vector products for real unsymmetric
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
*  (N+1)   Diagonal pointer of the matrix in Harwell-Boeing format.
*
*  ICOL    (input) INTEGER array
*  (NNZ)   Column indecies of the matrix in Harwell-Boeing format.
*
*  VAL     (input) DOUBLE PRECISION array
*  (NNZ)   Value of the matrix in Harwell-Boeing format.
*
*  X(N)    (input) DOUBLE PRECISION array
*          Input vector.
*
*  Y(N)    (output) DOUBLE PRECISION array
*          Output vector.
*
*  =====================================================================
*
!$omp parallel do private(S,J_PTR,I)
      DO I=1,N
         S=0.0D0
         DO J_PTR=IRP(I),IRP(I+1)-1
            S=S+VAL(J_PTR)*X(ICOL(J_PTR))
         END DO
         Y(I)=S
      END DO
!$omp end parallel do
*
      RETURN
      END
*
*
      SUBROUTINE OpenATI_DURMV_12(N,NNZ,IRP,ICOL,VAL,X,Y,
     &     NUM_SMP,KBORDER)
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*     ..
*     ..Scalar Argument
      INTEGER           N,NNZ,NUM_SMP
*     ..
*     ..Array Argument
      INTEGER          IRP(N+1),ICOL(NNZ),KBORDER(0:NUM_SMP)
      DOUBLE PRECISION X(N),Y(N),VAL(NNZ)
*
      INTEGER           I,K,J_PTR
      DOUBLE PRECISION  S
*
*=====================================================================
*  Purpose
*  =======
*
*  OpenATI_DURMV_12 computes matrix-vector products for real unsymmetric
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
*  IRP     (input) INTEGER 
*  (N+1)   Diagonal pointer of the matrix in Harwell-Boeing format.    
*
*  ICOL    (input) INTEGER array
*  (NNZ)   Column indecies of the matrix in Harwell-Boeing format.
*
*  VAL     (input) DOUBLE PRECISION array
*  (NNZ)   Value of the matrix in Harwell-Boeing format.
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
*  KBORDER (input) INTEGER array
*  (N+1)   partition infomation of matrix
*
*====================================================================
*
!$omp parallel do private(S,J_PTR,I)
      DO K=1,NUM_SMP
         DO I=KBORDER(K-1)+1,KBORDER(K)
            S=0.0D0
            DO J_PTR=IRP(I),IRP(I+1)-1
               S=S+VAL(J_PTR)*X(ICOL(J_PTR))
            END DO
            Y(I)=S
         END DO
      END DO
!$omp end parallel do
*
      RETURN
      END
*
*
      SUBROUTINE OpenATI_DURMV_13(N,NNZ,ICOL,VAL,X,Y,
     &     IS,VALSS,SUM,MFLAG,JFSTART,JYN,JSY,JSFLAG,PRESENT,JL)
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
      INTEGER           JL

*     ..
*     ..Scalar Argument
      INTEGER           IS,N,NNZ
*     ..
*     ..Array Argument
      INTEGER           ICOL(NNZ)
      DOUBLE PRECISION  X(N),Y(N),VAL(NNZ)
      DOUBLE PRECISION  VALSS(N+JL),SUM(JL)
      INTEGER           JFSTART(0:JL),JYN(JL),JSY(JL),MFLAG(N+JL)
      LOGICAL(1)        JSFLAG(JL),PRESENT(JL)
*
      DOUBLE PRECISION  S 
      INTEGER           I,J,K,I2,I4,II,II1
*
*=====================================================================
*  Purpose
*  =======
*
*  OpenATI_DURMV_13 computes matrix-vector products for real unsymmetric
*  sparse matrix with improved segment scan.
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
*  ICOL    (input) INTEGER array
*  (NNZ)   Column indecies of the matrix in Harwell-Boeing format.
*
*  VAL     (input) DOUBLE PRECISION array
*  (NNZ)   Value of the matrix in Harwell-Boeing format.
*
*  X(N)    (input) DOUBLE PRECISION array
*          Input vector.
*
*  Y(N)    (output) DOUBLE PRECISION array
*          Output vector.
*
*  IS      The order of SEGMENTED MATRIX 
*
*  VALSS   DOUBLE PRECISION array (NNZ)
*          Work area for segmented mat*vec
*
*  SUM     DOUBLE PRECISION array
*          Work area for spanning sum
*  
*  MFLAG   INTEGER array (N+JL)
*          Head element of each row in matrix and each column 
*          in segment array 
*
*  JFSTART INTEGER array (JL+1)
*          First location of last segment of each column
*
*  JYN     INTEGER array (JL)
*          Number of first location of segment of each column 
*
*  JSY     INTEGER array  (JL)
*          Sum of JYN   (JSY(j)=JYN(j-1)+JSY(j) , j=2,JL)
*         
*  JSFLAG  LOGICAL array (JL)
*          Flag that first element of each column
*
*  PRESENT LOGICAL array (JL)
*          Flag that first location of segment exists
*          in each column
*
*=====================================================================
*
*  MAT*VEC product of each column
*
!$OMP PARALLEL DO PRIVATE(S,K,I)
      DO J=1,JL
         DO K=JFSTART(J-1)+1,JFSTART(J)
            S=0.0D0
            DO I=MFLAG(K),MFLAG(K+1)-1
               S=S+VAL(I)*X(ICOL(I))
            END DO
            VALSS(K)=S
         END DO
      END DO
!$OMP END PARALLEL DO


*
*  Make spanning array   
*
!$OMP PARALLEL DO
      DO J=2,JL
         IF(JSFLAG(J) .EQV. .FALSE.)THEN
            SUM(J-1)=VALSS(JFSTART(J-1)+1)
         END IF
      END DO
!$OMP END PARALLEL DO
*
*  Sum up SUM array
*
      DO J=JL-1,1,-1
         IF(PRESENT(J+1) .EQV. .FALSE.)THEN
            SUM(J)=SUM(J)+SUM(J+1)
         END IF
      END DO
*
*  Spanning sum
*
!$OMP PARALLEL DO
      DO J=1,JL
         VALSS(JFSTART(J))=VALSS(JFSTART(J))+SUM(J)
      END DO
!$OMP END PARALLEL DO
*
*  Make output 
*
!$OMP PARALLEL DO PRIVATE(IN,JN,I)
      DO J=1,JL
         IN=JFSTART(J)
         JN=JSY(J)
         DO I=JYN(J),1,-1
            Y(JN)=VALSS(IN)
            IN=IN-1
            JN=JN-1
         END DO
      END DO
!$OMP END PARALLEL DO
*
      RETURN
      END
*
*
      SUBROUTINE OpenATI_DURMV_21(N,NNZ,ICOL,IRP,IS,VAL,X,Y,
     &     VALS,SUM,JLAST,JSTART,FLAG,PRESENT,JL)
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
      INTEGER           JL
*
*     ..Scalar Argumrnts..
      INTEGER           IS,N,NNZ
*     ..
*     ..Array Arguments
      LOGICAL(1)        FLAG(IS,JL),PRESENT(JL)
      INTEGER           JLAST(JL),JSTART(JL),IRP(N+1)
      INTEGER           ICOL(NNZ)
      DOUBLE PRECISION  VALS(IS,JL),SUM(JL),VAL(NNZ),X(N),Y(N)
*
      INTEGER           I,J,K
*
*=====================================================================
*  Purpose
*  =======
*
*  OpenATI_DURMV_21 computes matrix-vector products for real
*  unsymmetric sparse matrix with original segment scan.
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
*  ICOL    (input) INTEGER array
*  (NNZ)   Column indecies of the matrix in Harwell-Boeing format.
*
*  IRP     (input) INTEGER array
*  (N+1)   Diagonal pointer of the matrix in Harwell-Boeing format.
*
*  IS      The order of SEGMENTED MATRIX
*
*  VAL     (input) DOUBLE PRECISION array
*  (NNZ)   Value of the matrix in Harwell-Boeing format.
*
*  X(N)    (input) DOUBLE PRECISION array
*          Input vector.
*
*  Y(N)    (output) DOUBLE PRECISION array
*          Output vector.
*
*  VALS    DOUBLE PRECISION array (IS,JL)
*          Work area for segmented mat*vec
*
*  SUM     DOUBLE PRECISION array
*          Work area for SPANNING SUM
*
*  JLAST   INTEGER array (JL)
*          Number of last segmentation top of segmented column
*
*  JSTART  INTEGER array (JL)
*          MAT*VEC start point of each column
*
*  FLAG    LOGICAL array (IS,JL)
*          Head element of each row in matrix and each column
*          in segment array
*
*  PRESENT LOGICAL array (JL)
*          Flag that first location of segment exists
*          in each column
*
*=====================================================================
*
*  MAT*VEC product of each column
*
*
*  JLDASH is effective row number.
*
      JLDASH=(NNZ-1)/IS+1
!$omp parallel 
!$omp do PRIVATE(K)
      DO J=1,JLDASH
         K=ICOL(IS*(J-1)+JSTART(J))
         VALS(JSTART(J),J)=VAL(IS*(J-1)+JSTART(J))*X(K)
      END DO
!$omp end do
*
!$omp do PRIVATE(K,I)
      DO J=1,JLDASH
         DO I=JSTART(J)-1,1,-1
            K=ICOL(IS*(J-1)+I)
            IF(FLAG(I+1,J) .EQV. .true.)THEN
               VALS(I,J)=VAL(IS*(J-1)+I)*X(K)
            ELSE
               VALS(I,J)=VAL(IS*(J-1)+I)*X(K)+VALS(I+1,J)
            END IF
         END DO
      END DO
!$omp end do
!$omp end parallel
*
*  Make spanning array   
*
!$omp parallel do
      DO J=JLDASH,2,-1
         IF(FLAG(1,J) .EQV. .FALSE.)THEN
            SUM(J-1)=VALS(1,J)
         END IF
      END DO
!$omp end parallel do
*
*  Sum up sum array
*
      DO J=JLDASH-1,1,-1
         IF(PRESENT(J+1) .EQV. .FALSE.)THEN
            SUM(J)=SUM(J)+SUM(J+1)
         END IF
      END DO
*
*  Spanning sum
*
!$omp parallel do PRIVATE(I)
      DO J=1,JLDASH
         I=JLAST(J)
         IF(I .NE. 0 )THEN
            VALS(I,J)=VALS(I,J)+SUM(J)
         END IF
      end do
!$omp end parallel do
*
*  Make output
*
!$omp parallel do PRIVATE(J,I)
      DO K=1,N
         J=((IRP(K)-1)/IS)+1
         I=MOD(IRP(K)-1,IS)+1
         Y(K)=VALS(I,J)
      END DO
!$omp end parallel do
*
      RETURN
      END 

