      SUBROUTINE OpenATI_GET_VER(CVER)
*-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
      CHARACTER*20 CVER
*
*  =====================================================================
*  Purpose
*  =======
*
*  OpenATI_GET_VER get OpenATI Version
*
*
*  Arguments
*  =========
*
*  CVER      (output) Character*20
*              version of this OpenATI
*
*  =====================================================================
*
      include 'OpenATI_Ver.h'
      CVER = OpenATI_ver
*
      RETURN
      END
