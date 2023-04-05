      PROGRAM MAIN
      IMPLICIT NONE
C
      INTEGER NMAX, NZMAX
      parameter (NMAX=268100,NZMAX=9400000)

      INTEGER NTMP,NZTMP,NEVTMP
      INTEGER IRPTMP(NMAX+1), ICOLTMP(NZMAX)
      DOUBLE PRECISION ATMP(NZMAX)
      INTEGER N, NZ, NEV, INFO
      INTEGER IRP, ICOL, IATPARAM
      ALLOCATABLE :: IRP(:), ICOL(:), IATPARAM(:)
      DOUBLE PRECISION A, E, V, RATPARAM
      ALLOCATABLE :: A(:), E(:), V(:), RATPARAM(:)
      DOUBLE PRECISION WK, O
      ALLOCATABLE :: WK(:), O(:)
      INTEGER I,ITEST,MAXP,IP,OMP_GET_THREAD_NUM
C
      EXTERNAL OMP_GET_NUM_THREADS,OMP_GET_MAX_THREADS
      INTEGER  OMP_GET_NUM_THREADS,OMP_GET_MAX_THREADS
C
      open(31,file='Input.param',status='OLD')
      read(31,*) itest
      close(31)
      CALL MATGEN(ITEST,NTMP,NZTMP,IRPTMP,ICOLTMP,ATMP)
      MAXP=OMP_GET_MAX_THREADS()
C
      NEVTMP=10
      WRITE(6,*) '++++++++++++++ Input Parameter List +++++++++++++++++'
      write(6,*) '+ itest =',itest
      WRITE(6,*) '+ Matrix Info. N=',NTMP,' NZ=',NZTMP
      WRITE(6,*) '+ OpenMP Number of MAX. Threads=',MAXP
      WRITE(6,*) '+++++++++++++++++++++++++++++++++++++++++++++++++++++'
C
!$omp parallel default(none)
!$omp+ private(N, NZ, IRP, ICOL, A, NEV, E, V, INFO)
!$omp+ private(IP, IATPARAM, RATPARAM)
!$omp+ private(WK, O)
!$omp+ shared(NTMP, NZTMP, IRPTMP, ICOLTMP, ATMP, NEVTMP, ITEST)
      IP=OMP_GET_THREAD_NUM()
*
      N=NTMP
      NZ=NZTMP
      NEV=NEVTMP
      ALLOCATE(IRP(N+1),ICOL(NZ), A(NZ))
      ALLOCATE(E(2*NEV), V(2*N*NEV))
      ALLOCATE(IATPARAM(50), RATPARAM(50))
      ALLOCATE(WK(2*N),O(NEV*N))
      DO I=1,N+1
         IRP(I)=IRPTMP(I)
      ENDDO
      DO I=1,NZ
         ICOL(I)=ICOLTMP(I)
         A(I)=ATMP(I)
      END DO

      CALL OpenATI_INIT(IATPARAM, RATPARAM, INFO)
      write(6,*)'*** OpenATI_EIGENSOLVE THREAD-SAFE TEST ***',IP

C      IATPARAM(50)=1
      IATPARAM(30)=2
      CALL OpenATI_EIGENSOLVE(N,NZ,IRP,ICOL,A,NEV,E,V,
     $                        IATPARAM,RATPARAM,INFO)
!$omp barrier
      write(6,*) 'OpenATI_EIGENSOLVE INFO=',INFO
*
      if (info.lt.0) THEN
         write(6,*) '  !!!! Parameter Error !!! Info=',INFO
         GOTO 9000
      else if (info .ne.0) then
         write(6,*) '  !!!! Breakdown Error !!! Info=',INFO
         GOTO 9000
      end if
      IF(ITEST.GT.300 .AND. ITEST.LE.321) THEN
         call resid(n,irp,icol,nz,a,nev,e,v,n,wk)
         call ORTHO(N,nev,V,N,O)
      ELSE IF (ITEST.GT.200 .AND. ITEST.LE.222)then
         call residz(n,irp,icol,nz,a,nev,e,v,n,wk)
      END IF
 9000 CONTINUE
      DEALLOCATE(IRP,ICOL, A)
      DEALLOCATE(E, V)
      DEALLOCATE(IATPARAM, RATPARAM)
      DEALLOCATE(WK,O)
*
!$omp barrier
!$omp end parallel
      STOP
      END
*
      subroutine resid(n,irp,icol,nz,a,nev,e,v,nv1,r)
      implicit real*8 (a-h,o-z)
      integer*4 irp(n+1),icol(nz)
      real*8    a(nz),e(nev),v(nv1,nev),r(n)
C>>>>>>>>>>
      resmax=0.0D0
      do 100 ic=1,nev
C---------------mat*vec
        do 200 i=1,n
          r(i)=0.0d0
  200   continue
        do 210 i=1,n
          s=a(irp(i))*v(i,ic)
          do 220 jc=irp(i)+1,irp(i+1)-1
            jj=icol(jc)
            s=s+a(jc)*v(jj,ic)
            r(jj)=r(jj)+a(jc)*v(i,ic)
  220     continue
          r(i)=r(i)+s
  210   continue
C
        do 230 i=1,n
          r(i)=r(i)-e(ic)*v(i,ic)
  230   continue
C
        zansa=0.0d0
        do 240 i=1,n
          zansa=zansa+r(i)*r(i)
  240   continue
        write(6,*) 'IC=',IC,'E=',e(ic),'RES=',sqrt(zansa)/abs(e(ic))
        resmax=max(resmax,sqrt(zansa)/abs(e(ic)))
  100 continue
C
      WRITE(6,*) '================================================'
      WRITE(6,*) '=== MAX RESID    =',resmax
      WRITE(6,*) '================================================'
C
      return
      end
C********************
      SUBROUTINE ORTHO(N,NV,V,NV1,O)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 V(NV1,NV),O(NV1,NV)
C
      ICHK=0
      DO 400 J=1,NV
        DO 500 I=1,J
          S=0.0D0
          DO 600 K=1,N
            S=S+V(K,I)*V(K,J)
  600     CONTINUE
          IF (I.EQ.J) THEN
            IF (DABS(DSQRT(S)-1.0D0) .GT.1.0D-12) THEN
              ICHK=1
              WRITE(6,*) '!!!NG!!! EIGENVECTOR=',J,' IS NOT NORMALIZED'
     &                   ,SQRT(S)
C             RETURN
            END IF
          END IF
          O(I,J)=S
  500   CONTINUE
  400 CONTINUE
      ERR=0.0D0
      DO 700 J=1,NV
        DO 800 I=1,J-1
          IF (I.NE.J) ERR=MAX(ERR,O(I,J))
  800   CONTINUE
  700 CONTINUE
      IF (ICHK.EQ.0) THEN
        WRITE(6,*) '!!! OK !!! EIGENVECTOR NORMALIZED'
      END IF
      WRITE(6,*) '================================================'
      WRITE(6,*) '=== ORTHOGONALITY=',ERR
      WRITE(6,*) '================================================'
C
      RETURN
      END
      subroutine matgen(itest,n,nz,irp,icol,a)
      implicit real*8 (a-h,o-z)
      integer*4 irp(*),icol(*)
      real*8    a(*)
      character filename*60
C
      IF      (itest .EQ. 207) THEN
         filename="../MatrixData/ecl32.dat"
      else if (itest.eq.301) then
         filename='../MatrixData/vibrobox.rb'
      end if
C
      if(itest.gt.300 .and. itest.le.321) then
         call matread(itest,filename,n,irp,icol,nz,a)
      else if(itest.gt.200 .and. itest.le.222) then
         OPEN(5,FILE=filename)
         READ(5,*) N,NZ
         READ(5,*) (IRP (I),I=1,N+1)
         READ(5,*) (ICOL(I),I=1,NZ)
         READ(5,*) (A (I),I=1,NZ)
         CLOSE(5)
      endif
C
      return
      end
      subroutine matread(itest,filename,ncol,colptr,rowind,nnzero,
     *                   values)
      implicit real*8 (a-h,o-z)
c
C     ================================================================
C     ... SAMPLE CODE FOR READING A SPARSE MATRIX IN STANDARD FORMAT
C     ================================================================

      CHARACTER      TITLE*72 , KEY*8    , MXTYPE*3 ,
     1               PTRFMT*16, INDFMT*16, VALFMT*20, RHSFMT*20
      INTEGER        TOTCRD, PTRCRD, INDCRD, VALCRD, RHSCRD,
     1               NROW  , NCOL  , NNZERO, NELTVL
      INTEGER        COLPTR (*), ROWIND (*)
      REAL*8         VALUES (*)
      character     filename*60
C
      lunit=23
      open(lunit,file=filename)
      if (itest.eq.308) then
      READ ( LUNIT, 1100 ) TITLE , KEY   ,
     1                     TOTCRD, PTRCRD, INDCRD, VALCRD, RHSCRD,
     2                     MXTYPE, NROW  , NCOL  , NNZERO,
     3                     PTRFMT, INDFMT, VALFMT, RHSFMT
 1100   FORMAT ( A72, A8 / 5I14 / A3, 11X, 3I14 / 2A16, 2A20 )
      READ ( LUNIT, * )
      else
      READ ( LUNIT, 1000 ) TITLE , KEY   ,
     1                     TOTCRD, PTRCRD, INDCRD, VALCRD, RHSCRD,
     2                     MXTYPE, NROW  , NCOL  , NNZERO, NELTVL,
     3                     PTRFMT, INDFMT, VALFMT, RHSFMT
 1000   FORMAT ( A72, A8 / 5I14 / A3, 11X, 4I14 / 2A16, 2A20 )
      endif
      write(6,*) '===> INPUT FILE NAME IS ',filename
      write(6,*) TITLE
      write(6,*) KEY
C     -------------------------
C     ... READ MATRIX STRUCTURE
C     -------------------------
      READ ( LUNIT, PTRFMT ) ( COLPTR (I), I = 1, NCOL+1 )
      READ ( LUNIT, INDFMT ) ( ROWIND (I), I = 1, NNZERO )
      IF  ( VALCRD .GT. 0 )  THEN
C         ----------------------
C         ... READ MATRIX VALUES
C         ----------------------
          READ ( LUNIT, VALFMT ) ( VALUES (I), I = 1, NNZERO )
      ENDIF
      return
      end
*
      subroutine residz(n,irp,icol,nz,a,nev,e,v,nv1,r)
      implicit real*8 (a-h,o-z)
      integer*4 irp(n+1),icol(nz)
      real*8    a(nz)
      complex*16   e(nev),v(nv1,nev),r(n)
      complex*16   s
C
      resmax=0.0D0
      do 100 ic=1,nev
C---------------mat*vec
        do 210 i=1,n
          s=dcmplx(0.0d0,0.0d0)
          do 220 jc=irp(i),irp(i+1)-1
            jj=icol(jc)
            s=s+a(jc)*v(jj,ic)
  220     continue
          r(i)=s
  210   continue
        do 230 i=1,n
          r(i)=r(i)-e(ic)*v(i,ic)
  230   continue
        zansa=0.0d0
        do 240 i=1,n
          zansa=zansa+dreal(conjg(r(i))*r(i))
  240   continue
        write(6,*) 'IC=',IC,'E=',e(ic),'RES=',sqrt(zansa)/abs(e(ic))
        resmax=max(resmax,sqrt(zansa)/abs(e(ic)))
  100 continue
      WRITE(6,*) '================================================'
      WRITE(6,*) '=== MAX RESID    =',resmax
      WRITE(6,*) '================================================'
      return
      end
