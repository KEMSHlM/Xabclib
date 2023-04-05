      SUBROUTINE OpenATI_NZDIV(N,NNZ,NUM_SMP,DN,KBORDER)
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*     ..
*     ..Scalar Arguments..
      INTEGER N,NNZ,NUM_SMP
*     ..
*     ..Arry Arguments
      INTEGER  DN(N),KBORDER(0:NUM_SMP)
*
      INTEGER  OPT,IOMP,INUM
      INTEGER  NNNZ(NUM_SMP),JLL(0:NUM_SMP),VALUE(NUM_SMP)
*
*  =====================================================================
*  Purpose
*  =======
*   function for partition
*
*  Arguments
*  =========
*  N       (input) INTEGER
*          The order of the target block.  N >= 0.
*
*  NNZ     (input) INTEGER
*          Non-Zeros of the target block.  NNZ >= N.
*
*  NUM_SMP (input) INTEGER
*          Number of threads 
*
*  DN      (input) INTEGER
*  (N)     NON-ZERO number of each row in block
*
*  KBORDER (output) INTEGER array
*  (NUM    Partition infomation
*   _SMP+1)
*
*=======================================================================
*
      INUM=0
      IOMP=NUM_SMP
*
*  Calculate if NUM_SMP is multiplier of two
*
 1000 CONTINUE
      IF (IOMP .GE. 2) THEN
        IOMP=IOMP/2
        INUM=INUM+1
        GO TO 1000
      END IF
*
*  Divided of N to equation work balance of each thread 
*  case NUM_SMP=2**INUM 
*
      IF ((2**INUM) .EQ. NUM_SMP) THEN
        JLL(0)=0
        JLL(1)=N
        NNNZ(1)=NNZ
        DO II=1,INUM
          DO I=1,2**(II-1)
             CALL OpenATI_NZDIV2(N,NNNZ(I),2*I-1,JLL(I-1),
     &            JLL(I),DN,KBORDER,VALUE,NUM_SMP)
             KBORDER(I*2)=JLL(I)
          END DO
          DO J=1,NUM_SMP
             NNNZ(J)=VALUE(J)
          END DO
          DO J=0,NUM_SMP
             JLL(J)=KBORDER(J)
             END DO
        END DO
      ELSE
*
*  case NUM_SMP != 2**INUM
*
         CALL OpenATI_NZDIV3(N,NNZ,NUM_SMP,KBORDER,DN)
      END IF
*
      RETURN
      END
*
*
      SUBROUTINE OpenATI_NZDIV2(N,NNZ,K,JLS,JLF,DN,KBORDER,
     &     VALUE,NUM_SMP)
*
*     ..
*     ..Scalar Arguments..
      INTEGER :: NNZ,N,K,JLS,JLF,NUM_SMP
*     ..
*     ..Array Arguments
      INTEGER :: DN(N),KBORDER(0:N),VALUE(NUM_SMP)
*
      INTEGER :: OPT,DSUM1,DSUM2,I
*
*  =====================================================================
*  Purpose
*  =======
*   Function for division  case : OMP can be divided by 2  
*   Repeatition of half division
*
*
*  Arguments
*  =========
*  N       (input) INTEGER
*          The order of the target block.  N >= 0.
*
*  NNZ     (input) INTEGER
*          Non-Zeros of the target block.  NNZ >= N.
*
*  JLS     (input) INTEGER
*          Start line
*
*  JLF     (input) INTEGER
*          End line
*
*  DN      (input) INTEGER
*  (N)     NON-ZERO number of each row
*
*  KBORDER (output) INTEGER array
*  (NUM    Partition infomation
*   _SMP+1)
*   
*  VALUE   (output) INTEGER array
*  (NUM    Non-zero number of each block
*   _SMP)
*
*  NUM_SMP (input) INTEGER
*          Number of threads
*
*=======================================================================
*
      OPT=NNZ/2
      DSUM1=0
      DSUM2=0
*
*  Equating work balance   half to half
*
      DO I=JLS+1,JLF
         DSUM1=DSUM2
         DSUM2=DSUM2+DN(I)
         IF(DSUM2 > OPT)THEN
*
            IF(ABS(DSUM1-OPT) < ABS(DSUM2-OPT))THEN
               KBORDER(K)=I-1
               VALUE(K)=DSUM1
               GO TO 200
            ELSE
               KBORDER(K)=I
               VALUE(K)=DSUM2
               GO TO 200
            END IF
         END IF
      END DO
 200  CONTINUE
      VALUE(K+1)=NNZ-VALUE(K)
*
      RETURN
      END


      SUBROUTINE OpenATI_NZDIV3(N,NNZ,NUM_SMP,KBORDER,DN)
*
*     ..
*     ..Scalar Arguments..
      INTEGER  N,NNZ,NUM_SMP
*     ..
*     ..Array Arguments..
      INTEGER  KBORDER(0:NUM_SMP),DN(N)
*
      INTEGER  I,J,K
      INTEGER  DSUM1,DSUM2,OPT
*
*  =====================================================================
*  Purpose
*  =======
*   Function for division     case : OMP can not be divided by 2
*   divide one after another from top line
*
*
*  Arguments
*  =========
*  N       (input) INTEGER
*          The order of the target block.  N >= 0.
*
*  NNZ     (input) INTEGER
*          Non-Zeros of the target block.  NNZ >= N.
*
*  NUM_SMP (input) INTEGER
*          Number of threads
*
*  KBORDER (output) INTEGER array
*  (NUM    Partition infomation of matrix
*   _SMP)
*
*  DN      (input) INTEGER
*  (N)     NON-ZERO number of each row in block
*
*=======================================================================
*
      OPT=NNZ/NUM_SMP
      K=1
      DSUM1=0
      DSUM2=0
*
*  Equating work balance   one after another
*
      DO I=1,N
         DSUM1=DSUM2
         DSUM2=DSUM2+DN(I)
         IF(DSUM2 > OPT)THEN
            IF(K .EQ. NUM_SMP)THEN
               GO TO 300
            END IF
            IF(ABS(DSUM1-OPT) < ABS(DSUM2-OPT))THEN
               KBORDER(K)=I-1
               DSUM2=DN(I)
            ELSE
               KBORDER(K)=I
               DSUM2=0
            END IF
            K=K+1
            DSUM1=0
         END IF
      END DO
 300  CONTINUE
      KBORDER(NUM_SMP)=N
      KBORDER(0)=0
*
*
      RETURN
      END
