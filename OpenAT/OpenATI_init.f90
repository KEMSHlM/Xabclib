      SUBROUTINE OpenATI_init
     &                 (IATPARAM,RATPARAM,INFO)
*-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
      INTEGER          IATPARAM(50)
      DOUBLE PRECISION RATPARAM(50)
*
*  =====================================================================
*  Purpose
*  =======
*
*  OpenATI_init set default value of IATPARAM & RATPARAM
*
*
*  Arguments
*  =========
*
*  N         (input) INTEGER
*             vector length.
*
*  X(N)      (input/output) Double precision
*             input vector.
*
*
*  =====================================================================
      EXTERNAL OMP_GET_MAX_THREADS
      INTEGER  OMP_GET_MAX_THREADS
*
      DO I=1,50
         IATPARAM(I)=0
         RATPARAM(I)=0.0D0
      ENDDO
*-----------------------------------------------------------------------
*---------- OpenAT area [3:20] -----------------------------------------
*-----------------------------------------------------------------------
* int
      NUM_SMP     =  OMP_GET_MAX_THREADS()
      MMFLAG      =  0
      MDEL        =  5
      MPTH        = 10
      MV_AT_SYM   =  3
      MV_AT_USYM  =  3
      MVCASE_SYM  = 12
      MVCASE_USYM = 12
      JL          =128
      MGS         =  2
      MEMINFO     =  0
* double
      RATIO_MM    = 1.0D+2
      ALPHA       = 1.0D-2
*
*-----------------------------------------------------------------------
*---------- Xabclib area [21:50] ---------------------------------------
*-----------------------------------------------------------------------
* int
      MAX_ITER     = -1
      IPC_FLAG     =  1
      KIND_PRECOND =  4
      MSIZE        = 20
      MSIZE_INI    =  2
      IOPT         =  1
      ILUTFILL     =  5
      ISTG         =  0
      IMRI         =  0
* double
      STOP_TOL     =  1.0D-8
      TIME_MAX     = -1.0D0
      BD_ILU0      =  1.0D-8
      BD_ILUT      =  1.0D-8
      RMRT         =  0.0D0
*
*-----------------------------------------------------------------------
      IATPARAM(3)  = NUM_SMP
      IATPARAM(4)  = MMFLAG
      IATPARAM(5)  = MDEL
      IATPARAM(6)  = MPTH
      IATPARAM(7)  = MV_AT_SYM
      IATPARAM(8)  = MVCASE_SYM
      IATPARAM(9)  = MV_AT_USYM
      IATPARAM(10) = MVCASE_USYM
      IATPARAM(11) = JL
      IATPARAM(12) = MGS
      IATPARAM(14) = MEMINFO
*
      RATPARAM(4)  = RATIO_MM
      RATPARAM(6)  = ALPHA
*
*
      IATPARAM(22) = MAX_ITER
      IATPARAM(24) = IPC_FLAG
      IATPARAM(25) = KIND_PRECOND
      IATPARAM(26) = ILUTFILL
      IATPARAM(27) = MSIZE
      IATPARAM(28) = MSIZE_INI
      IATPARAM(30) = IOPT
      IATPARAM(33) = ISTG
      IATPARAM(34) = IMRI
*
      RATPARAM(23) = STOP_TOL
      RATPARAM(22) = TIME_MAX
      RATPARAM(25) = BD_ILU0
      RATPARAM(34) = RMRT
*
*
      INFO=0
      RETURN
      END
