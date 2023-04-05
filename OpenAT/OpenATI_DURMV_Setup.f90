      SUBROUTINE OpenATI_DURMV_Setup(N,NNZ,IRP,IATPARAM,RATPARAM,
     &                               UINF,LUINF,INFO)
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
      INTEGER           JL
      INTEGER           NNZLMT
      PARAMETER         (NNZLMT=1908874111)
*     ..
*     .. Scalar Arguments ..
      INTEGER           N,NNZ,ICASE,LUINF,INFO
*     ..
*     ..Array Arguments 
      INTEGER           IRP(N+1)
      DOUBLE PRECISION  UINF(LUINF)
*
      INTEGER           IATPARAM(50)
      DOUBLE PRECISION  RATPARAM(50)
*
*  =====================================================================
*  Purpose
*  =======
*
*  Setup for OpenATI_DURMV calculation imformation
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
*  IATPARAM (input) INTEGER array, dimension (50)
*          Integer parameter set for OpenATI and Xabclib.
*           IATPARAM(3)     : Number of threads.
*           IATPARAM(10)    : OpenATI_DURMV impl. method.
*           IATPARAM(11)    : Number of Segment vector.
*
*  RATPARAM (input) DOUBLE PRECISION array, dimension (50)
*          Double precision parameter set for OpenATI and Xabclib.
*
*  UINF    (output) DOUBLE PRECISION
*          Setup Information for OpenATI_DURMV
*
*  LUINF   (output) INTEGER
*          The size of UINF.
*          if IATPARAM(10) = 11 :
*             LUINF >= 0
*          else if IATPARAM(10) = 12 :
*             LUINF >= INT(0.5*NUM_SMP)+1
*          else if IATPARAM(10) = 13 :
*             LUINF >= INT(1.5*N)+INT(4.25*JL)+10
*          else if IATPARAM(10) = 21 :
*             LUINF >= INT(1.125*NNZ)+INT(2.125*JL)+10
*          (NUM_SMP=IATPARAM(3), JL=IATPARAM(11))
*
*  INFO    (output) INTEGER
*          =   0:  Successful exit
*          = 100:  Invalid IATPARAM(10) value.
*          = 200:  Upper limit of NNZ.
*          = 300:  Invalid LUINF value.
*
*  =====================================================================
*
      INFO=0
*
      NUM_SMP    = IATPARAM( 3)
      ICASE      = IATPARAM(10)
*------------------------------------------------------
*     IATPARAM(11)(JL) Automatic Set-up
*------------------------------------------------------
      IF ((IATPARAM(9) .EQ. 1) .OR. (IATPARAM(9) .EQ. 4)) THEN
        IF (IATPARAM(11). GT. NUM_SMP) THEN
           IATPARAM(11) = IATPARAM(11) - MOD(IATPARAM(11),NUM_SMP)
        END IF
      END IF
      JL         = IATPARAM(11)
*
      SELECT CASE (ICASE)
      CASE (11)
*
      CASE (12)
*------------------------------------------------------
*     CASE=12 : Equal partitioning of NNZ by NUM_SMP
*------------------------------------------------------
*
         IF (LUINF .LT. INT(0.5D0*NUM_SMP)+1) THEN
            INFO=300
            GO TO 9999 
         END IF
*
         IP_K     =1
         CALL OpenATI_DURMV_Setup_12(N,NNZ,IRP,NUM_SMP,UINF(IP_K))
*     
      CASE (13)
*------------------------------------------------------
*     CASE=13 : Improved segment scan
*------------------------------------------------------
*
         IF (LUINF .LT. INT(1.5D0*N)+INT(4.25D0*JL)+10) THEN
            INFO=300
            GO TO 9999
         END IF
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
         CALL OpenATI_DURMV_Setup_13(N,NNZ,IRP,IS,
     &        UINF(IP_S),UINF(IP_MF),UINF(IP_JFS),UINF(IP_JY),
     &        UINF(IP_JS),UINF(IP_JSF),UINF(IP_PRE),JL)
*
      CASE (21)
*------------------------------------------------------
*       CASE=21 : Original segment scan
*------------------------------------------------------
*
         IF (NNZ .GT. NNZLMT) THEN
            INFO=200
         ELSE
*
            IF (LUINF .LT. INT(1.125D0*NNZ)+INT(2.125D0*JL)+10) THEN
               INFO=300
               GO TO 9999
            END IF
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
            CALL OpenATI_DURMV_Setup_21(N,NNZ,IRP,IS,UINF(IP_S),
     &           UINF(IP_JLA),UINF(IP_JST),UINF(IP_F),UINF(IP_PRE),JL)
         END IF
*
      CASE DEFAULT
         INFO=100
      END SELECT
*
 9999 CONTINUE
      RETURN
      END


      SUBROUTINE OpenATI_DURMV_Setup_12(N,NNZ,IRP,NUM_SMP,KBORDER)
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*     ..
*     .. Scalar Arguments ..
      INTEGER           N,NNZ,NUM_SMP
*     ..
*     ..Array Arguments
      INTEGER           IRP(N+1),KBORDER(0:NUM_SMP)
*
      INTEGER           DN(N)
      DOUBLE PRECISION  VALUE(NUM_SMP)
*
*  =====================================================================
*  Purpose
*  =======
*
*  Setup for OpenATI_DURMV calculation with equal partitioning
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
*  (N+1)   Diagonal pointer of the matrix in Harwell-Boeing format.
*
*  NUM_SMP (input) INTEGER
*          Number of threads
*
*  KBORDER (output) INTEGER array
*  (NUM    partition infomation of matrix
*   _SMP+1) 
*=======================================================================
*
*  Initialize array
*
      DO I=1,NUM_SMP
         VALUE(I)=0.0D0
      END DO
      DO I=0,NUM_SMP
         KBORDER(I)=0
      END DO
!$OMP PARALLEL DO
      DO I=1,N
         DN(I)=0
      END DO
!$OMP END PARALLEL DO
*
*  NON-ZERO element count up each row
*
!$OMP PARALLEL DO
      DO I=1,N
         DN(I)=DN(I)+IRP(I+1)-IRP(I)
      END DO
!$OMP END PARALLEL DO
*
      IF(NUM_SMP .NE. 1)THEN
         CALL OpenATI_NZDIV(N,NNZ,NUM_SMP,DN,KBORDER)
      END IF
      KBORDER(NUM_SMP)=N
      KBORDER(0)=0
*
      RETURN
      END
*
*
*
      SUBROUTINE OpenATI_DURMV_Setup_13(N,NNZ,IRP,IS,
     &     SUM,MFLAG,JFSTART,JYN,JSY,JSFLAG,PRESENT,JL)
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
      INTEGER          JL  
*     ..
*     ..Scalar Argument
      INTEGER          IS,N,NNZ
*     ..
*     ..Array Argument
      DOUBLE PRECISION SUM(JL)
      INTEGER          MFLAG(N+JL),IRP(N+1),JFSTART(0:JL)
      INTEGER          JYN(JL),JSY(JL)
      LOGICAL(1)       JSFLAG(JL),PRESENT(JL)
*
      INTEGER          I,J,K,ILL,JP,IY,DI,H
      INTEGER          JFS(JL)
*
*=====================================================================
*  Purpose
*  =======
*
*  Setup for OpenATI_DURMV calculation with improved segment scan.
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
*  IS      The order of SEGMENTED MATRIX
*
*  SUM     DOUBLE PRECISION array
*  (JL)    Work area for spanning sum
*
*  MFLAG   INTEGER array
*  (N+JL)  Head element of each row in matrix and each column
*          in segment array
*
*  JFSTART INTEGER array
*  (JL+1)  First location of last segment of each column
*
*  JYN     INTEGER array
*  (JL)    Number of first location of segment of each column
*
*  JSY     INTEGER array
*  (JL)    Sum of JYN   (JSY(j)=JYN(j-1)+JSY(j) , j=2,JL)
*
*  JSFLAG  LOGICAL array
*  (JL)    Flag that first element of each column
*
*  PRESENT LOGICAL array
*  (JL)    Flag that first location of segment exists
*          in each column
*
*=====================================================================
*
*  Initialize array
*
      DO I=1,JL
         JSFLAG(I)=.FALSE.
         JFS(I)=1
         PRESENT(I)=.FALSE.
         SUM(I)=0.0D0
         JFSTART(I)=0
      END DO
*
*  Making JFS,SMFLAG,PRESENT
*
      ILL=1
      JP=(IRP(1)/IS)+1
      DO K=1,N
*
*  Seach first location of segment using IRP
*
         J=((IRP(K)-1)/IS)+1
         I=MOD((IRP(K)-1),IS)+1
*
*  Making PRESENT  
*
         PRESENT(J)=.TRUE.
*
*  Making JFS  
*
         IF(J .NE. JP)THEN
            ILL=1
         END IF
         IF((ILL .EQ. 1) .AND. (I .NE. 1))THEN
            ILL=ILL+1
         END IF
         JFS(J)=ILL
         ILL=ILL+1
*
*  Making JSFLAG 
*
         IF(I .EQ. 1)THEN
            JSFLAG(J)=.TRUE.
         END IF
         JP=J
      ENDDO
*
*  Making FLAG
*
      H=1
      JP=1
      DO K=1,N
*
*  Seach first location of segment(='T') using IRP
*
         J=((IRP(K)-1)/IS)+1
*
         IF(JP .NE. J)THEN
            DO I=JP+1,J
               IF(JSFLAG(I) .EQV. .FALSE.)THEN
                  MFLAG(H)=(I-1)*IS+1
                  H=H+1
               END IF
            END DO
         END IF
*
         MFLAG(H)=IRP(K)
         H=H+1
         JP=J
      END DO
      MFLAG(H)=NNZ+1
*
*  Makinig JYN and JSY
*
      IY=0
      DO J=1,JL
         DI=JFS(J)
         IF(JSFLAG(J) .EQV. .FALSE.)THEN
            DI=DI-1
         END IF
         IY=IY+DI
         JYN(J)=DI
         JSY(J)=IY
      END DO
*
*  Making JFSTART
*
      DO I=1,JL
         JFSTART(I)=JFSTART(I)+JFS(I)
      END DO
      DO I=2,JL
         JFSTART(I)=JFSTART(I)+JFSTART(I-1)
      END DO
*
      JFSTART(0)=0
*
      RETURN
      END
*
*
      SUBROUTINE OpenATI_DURMV_Setup_21(N,NNZ,IRP,IS,SUM,JLAST,
     &     JSTART,FLAG,PRESENT,JL)
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
      INTEGER           JL
*     ..
*     ..Scalar Arguments..
      INTEGER           IS,N,NNZ 
*     ..
*     ..Array Arguments
      LOGICAL(1)        FLAG(IS,JL),PRESENT(JL)
      INTEGER           JLAST(JL),JSTART(JL),IRP(N+1)
      DOUBLE PRECISION  SUM(JL)
*
*
*=====================================================================
*  Purpose
*  =======
*
*  Setup for OpenATI_DURMV calculation with original segment scan.
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
*  IS      The order of SEGMENTED MATRIX
*
*  JLAST   INTEGER array
*  (JL)    Number of last segmentation top of segmented column
*
*  JSTART  INTEGER array
*  (JL)    MAT*VEC start point of each column
*
*  FLAG    LOGICAL array
*  (NZ)    Head element of each row in matrix and each column
*          in segment array
*
*  PRESENT LOGICAL array
*  (JL)    Flag that first location of segment exists
*          in each column
*
*=====================================================================
*
*  Array initialize
*
      DO J=1,JL
         JSTART(J)=IS
         PRESENT(J)=.FALSE.
         JLAST(J)=0
         SUM(J)=0.0D0
      END DO
*
*  Making JSTART
*
      I=MOD(NNZ-1,IS)+1
*
*  JLDASH is effective row number.
*
      JLDASH=(NNZ-1)/IS+1
      JSTART(JLDASH)=I
      DO J=JLDASH+1,JL
         JSTART(J)=0
      END DO
      DO J=1,JLDASH
         DO I=1,JSTART(J)
            FLAG(I,J)=.FALSE.
         END DO
      END DO
*
*  Making FLAG,PRESENT,JLAST 
*
      DO K=1,N
         J=((IRP(K)-1)/IS)+1
         I=MOD((IRP(K)-1),IS)+1
         JLAST(J)=I
         PRESENT(J)=.TRUE.
         FLAG(I,J)=.TRUE.
      END DO
*
*
      RETURN
      END
