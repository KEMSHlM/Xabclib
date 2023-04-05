      program MAIN
      implicit double precision (a-h,o-z)
C
      parameter (nmax=268100,nzmax=9400000)
      parameter (nwkmax=4*nmax) 
      parameter (npremax=nmax) 
C
      integer*4 irp(nmax+1),icol(nzmax)
      real*8    a(nzmax),wk(nwkmax),precond(npremax)
      real*8    b(nmax),x(nmax)
      integer*4 nwk, npre
C
      integer*4 iatparam(50)
      real*8    ratparam(50)
C
      external omp_get_num_threads,omp_get_max_threads
      integer*4 omp_get_num_threads,omp_get_max_threads
C
      call OpenATI_INIT(IATPARAM,RATPARAM,INFO)
C
      maxiter=10000
      idbg=1
      stoptol=1.0d-8
      iatparam(25)=2
      iatparam(22)=maxiter
      iatparam(50)=idbg
      maxp=omp_get_max_threads()
      itest=301
      call matgen(itest,n,nz,irp,icol,a)
C
      write(6,*)'====================================================='
      write(6,*)'================= Xabclib_CG START =================='
      write(6,*)'====================================================='
      write(6,*) '++++++++++++++ Input Parameter List +++++++++++++++++'
      write(6,*) '+ Matrix Info. N=',n,' NZ=',NZ
      write(6,*) '+ Preconditioner        =',iatparam(25)
      write(6,*) '+ Convergence criterion =',ratparam(23)
      write(6,*) '+ MatVec AT             =',iatparam(7)
      write(6,*) '+++++++++++++++++++++++++++++++++++++++++++++++++++++'
C
      nwk=4*n
      npre=n
C
      do i=1,n
       x(i)=1.0d0
      end do
C
      call OpenATI_DSRMV_11(n,nz,irp,icol,a,x,b)
C
      do i=1,n
       x(i)=0.0d0
      end do
C
      call Xabclib_CG( n, nz, irp, icol, a, b, x,
     $                 precond, npre,
     $                 IATPARAM, RATPARAM,
     $                 WK, nwk, INFO )
C
      if (info.lt.0) THEN
         write(6,*) '  !!!! Parameter Error !!! Info=',INFO
         stop
      else if (info .ne.0) then
         write(6,*) '  !!!! Breakdown Error !!! Info=',INFO
         stop
      end if
      write(6,*) '<<<   Xabclib CG SUCCESSFUL EXIT    >>>'
      write(6,*) '<<<            RESULT               >>>'
      write(6,*) '  ----Mat*Vec Impl.         =',IATPARAM(10)
      write(6,*) '  ----Num. of Iteration     =',IATPARAM(23)
      write(6,*) '  ----Final residual        =',RATPARAM(29)
      write(6,*) '  ----Floating Operations   =',RATPARAM(30),'[Gflops]'
      write(6,*) '  ----Total solve time[sec] =',RATPARAM(32)
      write(6,*) '  -----------------------------------------'
C
 9999 CONTINUE
      stop
      end
C
C
C
      subroutine matgen(itest,n,nz,irp,icol,a)
      implicit real*8 (a-h,o-z)
      integer*4 irp(n+1),icol(*)
      real*8    a(*)
C
      character filename*60
C
      if (itest.eq.301) then
         filename='../MatrixData/vibrobox.rb'
      end if
C
      call matread(itest,filename,n,irp,icol,nz,a)
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

      if (itest.eq.308) then
      READ ( LUNIT, 1100 ) TITLE , KEY   ,
     1                     TOTCRD, PTRCRD, INDCRD, VALCRD, RHSCRD,
     2                     MXTYPE, NROW  , NCOL  , NNZERO,
     3                     PTRFMT, INDFMT, VALFMT, RHSFMT
 1100 FORMAT ( A72, A8 / 5I14 / A3, 11X, 3I14 / 2A16, 2A20 )
      READ ( LUNIT, * )
      else
      READ ( LUNIT, 1000 ) TITLE , KEY   ,
     1                     TOTCRD, PTRCRD, INDCRD, VALCRD, RHSCRD,
     2                     MXTYPE, NROW  , NCOL  , NNZERO, NELTVL,
     3                     PTRFMT, INDFMT, VALFMT, RHSFMT
 1000 FORMAT ( A72, A8 / 5I14 / A3, 11X, 4I14 / 2A16, 2A20 )
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
      close(lunit)
      return
      end
