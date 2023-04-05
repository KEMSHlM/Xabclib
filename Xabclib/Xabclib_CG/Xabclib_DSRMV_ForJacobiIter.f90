      SUBROUTINE Xabclib_DSRMV_ForJacobiIter(N,NNZ,IRP,ICOL,VAL,X,Y,
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
*  Xabclib_DSRMV_ForJacobiIter computes matrix-vector products for real
*  symmetric sparse matrix with formatted by Harwell-Boeing format.
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
*           IATPARAM(8)     : OpenATI_DSRMV impl. method.
*           IATPARAM(50)    : Debug print control flag.
*          On exit,
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
*
      NUM_SMP      = IATPARAM( 3)
      ICASE        = IATPARAM( 8)
*
      INFO = 0
*
      IF ( ICASE.EQ.11) THEN
*
         CALL OpenATI_DSRMV_31(N,NNZ,IRP,ICOL,VAL,X,Y)
*
      ELSE IF(ICASE .EQ. 12)THEN
*
         IP_K    =1
*
         CALL OpenATI_DSRMV_32(N,NNZ,IRP,ICOL,VAL,X,Y,NUM_SMP,
     &                         SINF(IP_K))
*
      ELSE IF(ICASE .EQ. 13)THEN
*
         IP_JLS  =1
         IP_JLN  =IP_JLS+(N-1)/2+1
         IP_KMB  =IP_JLN+(N-1)/2+1
         IP_KWB  =IP_KMB+NUM_SMP/2+1
*
         CALL OpenATI_DSRMV_33(N,NNZ,IRP,ICOL,VAL,X,Y,
     &                         NUM_SMP,WK,SINF(IP_JLS),SINF(IP_JLN),
     &                         SINF(IP_KWB),SINF(IP_KMB))
*
      ELSE
         INFO = 100
      ENDIF
*
*
      RETURN
      END
*
