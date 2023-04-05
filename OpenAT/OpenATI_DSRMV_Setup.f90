      SUBROUTINE OpenATI_DSRMV_Setup(N,NNZ,IRP,ICOL,
     &     IATPARAM,RATPARAM,SINF,LSINF,INFO)
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*     ..
*     ..Scalar Argument..
      INTEGER           N,NNZ,LSINF,INFO
*     ..
*     ..Array Argument..
      INTEGER           IRP(N+1),ICOL(NNZ)
      DOUBLE PRECISION  VAL(NNZ),X(N),Y(N)
*
      DOUBLE PRECISION  SINF(LSINF)
*
      INTEGER           IATPARAM(50)
      DOUBLE PRECISION  RATPARAM(50)
*
*  =====================================================================
*  Purpose
*  =======
*
*  Setup for OpenATI_DSRMV calculation
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
*  IATPARAM (input) INTEGER array, dimension (50)
*          Integer parameter set for OpenATI and Xabclib.
*           IATPARAM(3)     : Number of threads.
*           IATPARAM(8)     : OpenATI_DSRMV impl. method.
*
*  RATPARAM (input) DOUBLE PRECISION array, dimension (50)
*          Double precision parameter set for OpenATI and Xabclib.
*
*  SINF    (output) DOUBLE PRECISION array
*  (LSINF) Setup Information for OpenATI_DSRMV.
*
*  LSINF   (input) INTEGER
*          The size of SINF.
*          if IATPARAM(8) = 11 :
*             LSINF >= 0
*          else if IATPARAM(8) = 12 :
*             LSINF >= INT(0.5*NUM_SMP)+1
*          else if IATPARAM(8) = 13 :
*             LSINF >= N+NUM_SMP+3
*          (NUM_SMP=IATPARAM(3))
*
*  INFO    (output) INTEGER
*          =   0:  Successful exit
*          = 100:  Invalid IATPARAM(8) value.
*          = 200:  Invalid LSINF value.
*
*  =====================================================================
*---------------mat*vec
      INFO=0
*
      NUM_SMP = IATPARAM( 3)
      ICASE   = IATPARAM( 8)
*
      SELECT CASE (ICASE)
      CASE (11)
*   
      CASE (12)
*------------------------------------------------------
*     CASE=12 : Equal partitioning of NNZ by NUM_SMP
*------------------------------------------------------
*
         IP_K    =1
*
         IF (LSINF .LT. INT(0.5D0*NUM_SMP)+1) THEN
            INFO=200
            GO TO 9999
         END IF
*
         CALL OpenATI_DSRMV_Setup_12(N,NNZ,IRP,ICOL,
     &        NUM_SMP,SINF(IP_K))
*
*
      CASE (13)
*------------------------------------------------------
*     CASE=13 : Equal partitioning of NNZ by NUM_SMP
*                     and parallel vector reduction run
*------------------------------------------------------
*
         IF (LSINF .LT. N+NUM_SMP+3) THEN
            INFO=200
            GO TO 9999
         END IF
*
         IP_JLS  =1
         IP_JLN  =IP_JLS+(N-1)/2+1
         IP_KMB  =IP_JLN+(N-1)/2+1
         IP_KWB  =IP_KMB+NUM_SMP/2+1
*
         CALL OpenATI_DSRMV_Setup_13(N,NNZ,IRP,ICOL,
     &        NUM_SMP,SINF(IP_JLS),SINF(IP_JLN),
     &        SINF(IP_KMB),SINF(IP_KWB))
*
      CASE DEFAULT
         INFO=100
      END SELECT
*
 9999 CONTINUE
      RETURN
      END


      SUBROUTINE OpenATI_DSRMV_Setup_12(N,NNZ,IRP,ICOL,
     &     NUM_SMP,KMBORDER)
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*     ..
*     ..Scalar Arguments..
      INTEGER  N,NUM_SMP,NNZ
*     ..
*     ..Array Arguments
      INTEGER  ICOL(NNZ),IRP(N+1),KMBORDER(0:NUM_SMP)
*
      INTEGER  JNM(N)
*
*  =====================================================================
*  Purpose
*  =======
*
*  Setup for OpenATI_DSRMV calculation with equal partitioning
*  of N by NUM_SMP
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
*  NUM_SMP (input) INTEGER
*          Number of threads
*
*  KMBORDER(output) INTEGER array
*  (NUM    Partition infomation of matrix
*   _SMP+1)
*
*=======================================================================
*
*  Initialize array
*
!$OMP PARALLEL
!$OMP DO
      DO I=0,NUM_SMP
         KMBORDER(I)=0
      END DO
!$OMP END DO
!$OMP DO
      DO I=1,N
         JNM(I)=0
      END DO
!$OMP END DO
!$OMP END PARALLEL
*
*  NON-ZERO element count up each row 
*
!$OMP PARALLEL DO
      DO I=1,N
         JNM(I)=IRP(I+1)-IRP(I)
      END DO
!$OMP END PARALLEL DO
      IF(NUM_SMP .NE. 1)THEN
         CALL OpenATI_NZDIV(N,NNZ,NUM_SMP,JNM,KMBORDER)
      END IF
      KMBORDER(NUM_SMP)=N
      KMBORDER(0)=0
*
*
      RETURN
      END

      SUBROUTINE OpenATI_DSRMV_Setup_13(N,NNZ,IRP,ICOL,
     &     NUM_SMP,JLS,JLN,KMBORDER,KWBORDER)
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*     ..
*     ..Scalar Arguments..
      integer  N,NUM_SMP,NNZ
*     ..
*     ..Array Arguments
      INTEGER  ICOL(NNZ),IRP(N+1)
      INTEGER  JLS(N),JLN(N)
      INTEGER  KWBORDER(0:NUM_SMP)
      INTEGER  KMBORDER(0:NUM_SMP)
*
      INTEGER  LAST(NUM_SMP),ILVS(N)
      INTEGER  ILF(N),JNW(N),JNM(N)
*
*  =====================================================================
*  Purpose
*  =======
*
*  Setup for OpenATI_DSRMV calculation with Equal partitioning
*  of NNZ by NUM_SMP and parallel vector reduction run
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
*  NUM_SMP (input) INTEGER
*          Number of threads
*
*  JLS     (output) INTEGER array
*  (N)     Start point of vector reduction part
*
*  JLN     (output) INTEGER array
*  (N)     End point of vector reduction part
*
*  KMBORDER(output) INTEGER array
*  (NUM    Partition infomation of matrix
*   _SMP+1)
*
*  KWBORDER(output) INTEGER array
*  (NUM    Partition infomation in vector reduction array
*   _SMP+1) 
*
*=======================================================================
*
*  Initialize array
*
!$OMP PARALLEL
!$OMP DO
      DO I=0,NUM_SMP
         KWBORDER(I)=0
         KMBORDER(I)=0
      END DO
!$OMP END DO
!$OMP DO
      DO I=1,NUM_SMP
         LAST(I)=0
      END DO
!$OMP END DO
!$OMP DO
      DO I=1,N
         JNW(I)=0
         JNM(I)=0
         JLS(I)=0
         JLN(I)=0
         ILF(I)=0
         JLS(I)=NUM_SMP
      END DO
!$OMP END DO
!$OMP END PARALLEL 
*
*  NON-ZERO element count up each row
*
!$OMP PARALLEL DO
      DO I=1,N
         JNM(I)=IRP(I+1)-IRP(I)
      END DO
!$OMP END PARALLEL DO
      IF(NUM_SMP .NE. 1)THEN
         CALL OpenATI_NZDIV(N,NNZ,NUM_SMP,JNM,KMBORDER)
      END IF
      KMBORDER(NUM_SMP)=N
      KMBORDER(0)=0
*
!$OMP PARALLEL DO
      DO I=1,N
         ILF(I)=ICOL(IRP(I+1)-1)
      END DO
!$OMP END PARALLEL DO
*
*  Making LAST array and JLN array
*
!$OMP PARALLEL DO PRIVATE(K)
      DO K=1,NUM_SMP
         DO I=KMBORDER(K-1)+1,KMBORDER(K)
            IF(ILF(I) > LAST(K))THEN
               LAST(K)=ILF(I)
            END IF
            JLN(I)=K
         END DO
      END DO
!$OMP END PARALLEL DO
*
*  Make JLS aray
*
      DO K=NUM_SMP-1,1,-1
         DO I=1,N
            IF(LAST(K) .GE. I)THEN
               JLS(I)=K
            END IF
         END DO
      END DO
*
*  Vector reduction array partitioning of NNZ by NUM_SMP
*
      NWZ=0
!$OMP PARALLEL DO
      DO I=1,N
         JNW(I)=JLN(I)-JLS(I)+1
      END DO
!$OMP END PARALLEL DO
      DO I=1,N
         NWZ=NWZ+JNW(I)
      END DO
      IF(NUM_SMP .NE. 1)THEN
         CALL OpenATI_NZDIV(N,NWZ,NUM_SMP,JNW,KWBORDER)
      END IF
      KWBORDER(NUM_SMP)=N
      KWBORDER(0)=0
*
*
      RETURN
      END
