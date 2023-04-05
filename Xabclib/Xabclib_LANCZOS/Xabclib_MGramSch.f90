      SUBROUTINE Xabclib_MGramSch(N,X,Q,MM)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION  X(N),Q(N,MM)
*
*  =====================================================================
*  Purpose
*  =======
*
*  Xabclib_MGramSch computes orthogonal vector by Modified Gram Schmidt
*  algorithm.
*
*      X(N)    |     Q(N,1:MM)
*             ---
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          The order of the matrix.  N >= 0.
*
*  X(N)    (input/output) DOUBLE PRECISION array
*          input vector. On exit, X is orthogonalize to Q.
*
*  Q(N,MM) (input) DOUBLE PRECISION array
*          Orthonormal vectors.
*
*  =====================================================================
*
      DO  J=1,MM
         S=0.0D0
!$omp parallel
!$omp do reduction(+:S)
         DO I=1,N
            S=S+Q(I,J)*X(I)
         ENDDO
!$omp end do
!$omp do
         DO I=1,N
            X(I)=X(I)-S*Q(I,J)
         ENDDO
!$omp end do
!$omp end parallel
      ENDDO
*
*
      RETURN
      END
