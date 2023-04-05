      SUBROUTINE OpenATI_DAFRT(NSAMP,SAMP,IRT,IATPARAM,RATPARAM,INFO)
*
*     OpenATI_DAFRT : Determine whether subspace expand is required.
*
*     December 2008
*     December 2011
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER  IPARM(10)
      DOUBLE PRECISION  SAMP(NSAMP),RPARM(10)
      CHARACTER*4 PDEBUG
      INTEGER           IATPARAM(50)
      DOUBLE PRECISION  RATPARAM(50)
*
*
*  =====================================================================
*  Purpose
*  =======
*
*  OpenATI_DAFRT  : Determine whether subspace expand is required.
*
*
*  Arguments
*  =========
*
*  NSAMP   (input) INTEGER
*            The number of sampling data. 
*  SAMP    (input) DOUBLE PRECISION array
*  (NSAMP)   Sampling Data.
*
*  IRT     (output) INTEGER
*            If IATPARAM(4) = 1 : 
*               IRT=0  :  Don't need Subspace expand.
*               IRT=1  :  Need Subspace expand.
*
*  IATPARAM (input/output) Integer array, dimension (50)
*          Integer parameter set for OpenATI and Xabclib.
*           IATPARAM(4)     : Flag of Krylov subspace expand by
*                             MM-ratio.
*           IATPARAM(5)     : Value of Krylov subspace expand.
*
*  RATPARAM (input/output) DOUBLE PRECISION array, dimension (50)
*          Double precision parameter set for OpenATI and Xabclib.
*           RATPARAM(4)     : Threshhold of MM-ratio.
*           RATPARAM(5)     : Value of MM-ratio.
*
*  INFO    (output) INTEGER
*          =   0:  Successful exit
*
*  =====================================================================
*>>>>>>>>>>>>>>>>Get Debug print Env.
      IDBG = IATPARAM(50)
*
*---------- MSIZE_AT FLAG & MM-RATIO
      IPARM(1) = IATPARAM(4)
      RPARM(1) = RATPARAM(4)
*
      INFO=0
*
      IRT=0
      VMM=0.0D0
*
      IF (IPARM(1).NE.1) THEN
         RETURN
      ENDIF
*
      IF ( RPARM(1) .EQ. 0.0D0 ) THEN
         X=100.0D0
      ELSE
         X=RPARM(1)
      END IF
      THRESH=1.0D0/X
*
      SMAX=SAMP(1)
      SMIN=SAMP(1)
      DO J=2,NSAMP
         IF (SAMP(J).NE.0.0D0) THEN
            IF      ( SAMP(J) .GT. SMAX ) THEN
               SMAX=SAMP(J)
            ELSE IF ( SAMP(J) .LT. SMIN ) THEN
               SMIN=SAMP(J)
            ENDIF
         ENDIF
      ENDDO
*
      RMM=SMIN/SMAX
      IF ( RMM .GT.THRESH ) THEN
         IRT=1
      ENDIF
      IF (RMM.NE.0.0D0) THEN
         VMM=1.0D0 / RMM
      ENDIF
         RATPARAM(5)=VMM
      IF (IDBG.EQ.1) THEN
         WRITE(6,*) ' OpenATI_DAFRT: MM ratio=',VMM,' IRT=',IRT
      ENDIF
*
      RETURN
      END
