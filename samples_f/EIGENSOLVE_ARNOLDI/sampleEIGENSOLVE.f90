      program MAIN
      implicit double precision (a-h,o-z)
C
      parameter (nmax=300000,nzmax=9400000,nj=200)
C
      integer*4 irp(nmax+1),icol(nzmax),iopt(5)
      real*8    a(nzmax),e(2*nmax),v(2*nmax*nj),wk(2*nmax)
      real*8    o(nmax*nj)
C
      integer*4 iparm(10),iat(10)
      real*8    rparm(10)
C
      integer*4 iatparam(50)
      real*8    ratparam(50)
C
      external omp_get_num_threads,omp_get_max_threads
      integer*4 omp_get_num_threads,omp_get_max_threads
C
      open(31,file='Input.param',status='OLD')
      read(31,*) itest
      close(31)
C
      nev=10
C
      call matgen(itest,n,nz,irp,icol,a)
C
      MAXP=omp_get_max_threads()
C
      write(6,*) '+++++++++++++++++++++++++++++++++++++++++++++++++++++'
      write(6,*) '+ itest =',itest
      WRITE(6,*) '+ Matrix Info. N=',n,' NZ=',NZ
      write(6,*) '+ OpenMP Number of MAX. Threads=',MAXP
      write(6,*) '+++++++++++++++++++++++++++++++++++++++++++++++++++++'
C
      CALL OpenATI_INIT(IATPARAM, RATPARAM, INFO)
      IATPARAM(30)=2
C     IATPARAM(50)=1
      CALL OpenATI_EIGENSOLVE(N,NZ,IRP,ICOL,A,NEV,E,V,
     $                        IATPARAM,RATPARAM,INFO)
      write(6,*) 'OpenATI_EIGENSOLVE INFO=',INFO
*
      if (info.lt.0) THEN
         write(6,*) '  !!!! Parameter Error !!! Info=',INFO
         stop
      else if (info .ne.0) then
         write(6,*) '  !!!! Breakdown Error !!! Info=',INFO
         stop
      end if
      if(itest.gt.300 .and. itest.le.321) then
        call resid(n,irp,icol,nz,a,nev,e,v,n,wk)
        call ORTHO(N,nev,V,N,O)
      else if(itest.gt.200 .and. itest.le.222)then
        call residz(n,irp,icol,nz,a,nev,e,v,n,wk)
      end if
*
      stop
      end
      subroutine resid(n,irp,icol,nz,a,nev,e,v,nv1,r)
      implicit real*8 (a-h,o-z)
      integer*4 irp(n+1),icol(nz)
      real*8    a(nz),e(nev),v(nv1,nev),r(n)

C>>>>>>>>>>
      resmax=0.0D0
C<<<<<<<<<<

C
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
C       write(6,*) 'IC=',IC,'E=',e(ic),'RES=',sqrt(zansa)
        write(6,*) 'IC=',IC,'E=',e(ic),'RES=',sqrt(zansa)/abs(e(ic))

C>>>>>>>>>>
        resmax=max(resmax,sqrt(zansa)/abs(e(ic)))
C<<<<<<<<<<


C
  100 continue
C
C>>>>>>>>>>
      WRITE(6,*) '================================================'
      WRITE(6,*) '=== MAX RESID    =',resmax
      WRITE(6,*) '================================================'
C<<<<<<<<<<
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
      integer*4 irp(n+1),icol(*)
      real*8    a(*)
C
      character filename*60
C
      if (itest.eq.216) THEN
         filename='../MatrixData/ecl32.dat'
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

c
      character     filename*60
C
      lunit=23
      open(lunit,file=filename)
C    ------------------------
C     ... READ IN HEADER BLOCK
C     ------------------------

C      if (itest.le.304) then
C      READ ( LUNIT, 1000 ) TITLE , KEY   ,
C     1                     TOTCRD, PTRCRD, INDCRD, VALCRD, RHSCRD,
C     2                     MXTYPE, NROW  , NCOL  , NNZERO, NELTVL,
C     3                     PTRFMT, INDFMT, VALFMT, RHSFMT
C 1000 FORMAT ( A72, A8 / 5I14 / A3, 11X, 4I14 / 2A16, 2A20 )
C      else
C      READ ( LUNIT, 1100 ) TITLE , KEY   ,
C     1                     TOTCRD, PTRCRD, INDCRD, VALCRD, RHSCRD,
C     2                     MXTYPE, NROW  , NCOL  , NNZERO,
C     3                     PTRFMT, INDFMT, VALFMT, RHSFMT
C 1100 FORMAT ( A72, A8 / 5I14 / A3, 11X, 3I14 / 2A16, 2A20 )
C
C      READ ( LUNIT, * )
C      endif
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
C
C
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
C
        do 230 i=1,n
          r(i)=r(i)-e(ic)*v(i,ic)
  230   continue
C
        zansa=0.0d0
        do 240 i=1,n
          zansa=zansa+dreal(conjg(r(i))*r(i))
  240   continue
C       write(6,*) 'IC=',IC,'E=',e(ic),'RES=',sqrt(zansa)
        write(6,*) 'IC=',IC,'E=',e(ic),'RES=',sqrt(zansa)/abs(e(ic))
        resmax=max(resmax,sqrt(zansa)/abs(e(ic)))
C
  100 continue

      WRITE(6,*) '================================================'
      WRITE(6,*) '=== MAX RESID    =',resmax
      WRITE(6,*) '================================================'
C
      return
      end
