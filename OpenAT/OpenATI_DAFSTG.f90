      SUBROUTINE OpenATI_DAFSTG (ISTGCNT,EMA,
     $                           RERR,PERR,STOP_TOL,
     $                           ITER,MAX_ITER,
     $                           ETIME,EITRTIME,MAX_ETIME,
     $                           IATPARAM,RATPARAM,INFO)
      IMPLICIT DOUBLE PRECISION ( A-H , O-Z )
*     ..
*     .. Scalar Auruments ..
      INTEGER           ISTGCNT, ITER, MAX_ITER, INFO
      DOUBLE PRECISION  RERR, PERR, STOP_TOL, EMA
      DOUBLE PRECISION  ETIME, EITRTIME, MAX_ETIME
*     ..
*     .. Array Arguments
      INTEGER           IATPARAM(50)
      DOUBLE PRECISION  RATPARAM(50)
*     ..
*
*  =====================================================================
*  Purpose
*  =======
*
*  OpenATI_DAFSTG detects the stagnation of relative residual.
*
*  Arguments
*  =========
*
*  ISTGCNT (input/output) INTEGER
*          The counter for detectiong stagnation of rerlative residual.
*
*  EMA      (input/output) DOUBLE PRECISION
*           The exponential moving average of relative residual.
*
*  RERR    (input) DOUBLE PRECISION
*          The error of the approximate solution vector.
*          ( RERR = || b - A x_i ||_2 / ||  b ||_2 )
*
*  PERR    (input) DOUBLE PRECISION
*          The previous error of the approximate solution vector.
*
*  STOP_TOL(input) DOUBLE PRECISION
*          The convergence criterion.
*
*  ITER    (input) INTEGER
*          The number of iterations.
*
*  MAX_ITER(input) INTEGER
*          The maximum iterations.
*
*  ETIME   (input) DOUBLE PRECISION
*          The solver elapsed time. 
*
*  EITRTIME(input) DOUBLE PRECISION
*          The elapsed time per 1 step iteration.
*
*  MAX_TIME(input) DOUBLE PRECISION
*          The time tolerant.
*
*  IATPARAM (input/output) INTEGER array, dimension (50)
*
*  RATPARAM (input/output) DOUBLE PRECISION array, dimension (50)
*
*  INFO     (output) INTEGER
*           Error code.
*
*  Local valirables
*  ================
*
*  ALPHA   DOUBLE PRECISION
*          The coefficient of exponential moving average.
*
*  =====================================================================
*     ..
*     .. Local variables
*
      INTEGER          IPRETOL, IDBG
      DOUBLE PRECISION ALPHA, PRETOL, RES, PRERES, STL
*
*
*
*
      ALPHA    = RATPARAM(6)
      IDBG     = IATPARAM(50)
      ICNTADD  = 1
      IF (EMA.EQ.0.0D0) THEN
       EMA=ALPHA
      ENDIF
*
         RES=LOG10(RERR)
         PRERES=LOG10(PERR)
         STL=LOG10(STOP_TOL)
         EMA = ALPHA * ( RES - PRERES ) + ( 1 - ALPHA ) * EMA
*
         IF(MAX_ETIME.GT.0.0D0) THEN
*
          IPRETOL=MIN((MAX_ITER-ITER) * 1.0D0,
     $     (MAX_ETIME-ETIME)/(EITRTIME)+1.0D0)
          ICNTADD= EITRTIME * 100 / MAX_ETIME + 1
          ICNTADD= MIN(9,ICNTADD)
*
         ELSE
*
          IPRETOL=(MAX_ITER-ITER)
*
         ENDIF
*
         PRETOL=RES+(EMA*IPRETOL)
*
         IF (IDBG .EQ. 1) THEN
          WRITE(6,*)
     $    '     EMA=',EMA,'IPRE,PRETOL=',IPRETOL,PRETOL
         ENDIF
*
         IF (RERR.LT.STOP_TOL*1.0D20) THEN
*
          IF (PRETOL.LE.STL .AND. RES.LE.0.0D0) THEN
*
           ISTGCNT=0
* 
          ELSE
*
           ISTGCNT=ISTGCNT+ICNTADD
*
           IF (IDBG .EQ. 1) THEN
            WRITE(6,*) '   ISTGCNT=',ISTGCNT
           ENDIF
*
          ENDIF
*
         ELSE
*
          ISTGCNT=9999
*
         ENDIF
*
       INFO=0
       RETURN
      END
