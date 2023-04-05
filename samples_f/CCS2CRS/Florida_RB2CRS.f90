      PROGRAM MAIN
*
*-----------------------------------------------------------------------
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
      integer nzmax, nmax
      parameter (nzmax = 10000000, nmax = 700000)
      integer ptr (nmax), index (nzmax), n, nz, totcrd, ptrcrd,
     $        indcrd, valcrd, rhscrd, ncol, nrow, nrhs, row, col, p
      character title*72, key*30, type*3, ptrfmt*16,
     $        indfmt*16, valfmt*20, rhsfmt*20
      logical sym
      double precision Value (nzmax), skew, myrand
      character rhstyp*3
      integer nzrhs, nel
*
      integer IRP (nmax), ICOL (nzmax), IVAL(NZMAX)
      double precision VAL (nzmax)
*
      DATA ZERO,DCM8,DCM16,PT1 / 0.0D0,1.0D-8,1.0D-16,0.1D0 /
*
      PARAMETER (MGMAX=100,NIPC=1,NRPC=NMAX,NIWK=1)
      DIMENSION B(NMAX),X(NMAX),EPS(2),IOPT(3)
      DIMENSION IDP(NMAX),
     $          IWK(NMAX),RWK((MGMAX+3)*NMAX+2*(MGMAX+1)*(MGMAX+1))
*
      PARAMETER (NZMAX2=NMAX*(20*2+1),IWPC=NZMAX2*3/2+NMAX*3)
      DOUBLE PRECISION RGRPARM(10),RPARM(10),PRECOND(IWPC)
      INTEGER IGRPARM(10),JPARM(10),IPARM(10),INFO,IDIAGP(NMAX)
      INTEGER IAT(10)
*
      INTEGER IPCPARM(10)
      INTEGER IATPARAM(50)
      DOUBLE PRECISION RPCPARM(10)
      DOUBLE PRECISION Y(NMAX),Z(NMAX)
*
      REAL ETIME,TARRAY(2),TIME0,TIME1
      CHARACTER FNAME32 *100
*
      LDWORK=(MGMAX+3)*NMAX+2*(MGMAX+1)*(MGMAX+1)
*
      OPEN(15,FILE="param.dat")
      READ(15,5001) FNAME32
 5001 FORMAT(A100)
      OPEN(5,FILE=FNAME32)
      WRITE(6,*) '---------'
      WRITE(6,*) ' MATRIX DATA FILE = ',FNAME32
      WRITE(6,*) '---------'
*
*     OPEN(5,FILE='./memplus/memplus.rb')
*     OPEN(5,FILE='./chipcool0/chipcool0.rb')
*     OPEN(5,FILE='./wang3/wang3.rb')
*     OPEN(5,FILE='./chem_master1/chem_master1.rb')
*     OPEN(5,FILE='./epb3/epb3.rb')
*     OPEN(5,FILE='./matrix_9/matrix_9.rb')
*
*-----------------------------------------------------------------------
*       read header information from Harwell/Boeing matrix

        read (5, 10, err = 998)
     $          title, key,
     $          totcrd, ptrcrd, indcrd, valcrd,
     $          type, nrow, ncol, nz, nel,
     $          ptrfmt, indfmt, valfmt
        if (rhscrd .gt. 0) then
*          new Harwell/Boeing format:
           read (5, 20, err = 998) rhstyp,nrhs,nzrhs
           endif
10      format (a72, a8 / 4i14 / a3, 11x, 4i14 / 2a16, 2a20)
20      format (a3, 11x, 2i14)

        skew = 0.0
        if (type (2:2) .eq. 'Z' .or. type (2:2) .eq. 'z') skew = -1.0
        if (type (2:2) .eq. 'S' .or. type (2:2) .eq. 's') skew =  1.0
        sym = skew .ne. 0.0

        write (6, 30)
     $          title, key,
     $          totcrd, ptrcrd, indcrd, valcrd, rhscrd,
     $          type, nrow, ncol, nz, nel,
     $          ptrfmt, indfmt, valfmt, rhsfmt
        if (rhscrd .gt. 0) then
*          new Harwell/Boeing format:
           write (6, 40) rhstyp,nrhs,nzrhs
        endif
30      format (
     $          ' title: ', a72 /
     $          ' key: ', a8 /
     $          ' Lines: tot: ', i14,' ptr: ',i14,' ind: ',i14 /
     $          '        val: ', i14,' rhs: ',i14 /
     $          ' type: ', a3, ' nrow: ', i14, ' ncol: ', i14 /
     $          ' nz: ', i14, ' elements: ', i14, /
     $          ' ptrfmt: ', a20, ' rowfmt: ', a20, /
     $          ' valfmt: ', a20, ' rhsfmt: ', a20)
40      format (' rhstyp: ', a3, ' nrhs: ', i14, ' nzrhs: ', i14)
        write (0, *) ' sym: ', sym, ' skew: ', skew

        n = max (nrow, ncol)

        if (n .ge. nmax .or. nz .gt. nzmax) then
           write (0, *) ' Matrix too big!'
           stop
        endif

        read (5, ptrfmt, err = 998) (Ptr (p), p = 1, ncol+1)
        read (5, indfmt, err = 998) (Index (p), p = 1, nz)

        do 55 col = ncol+2, n+1
           Ptr (col) = Ptr (ncol+1)
55      continue

*       read the values, or create random-valued matrix
        write(6,*) 'valfmt=',valfmt
        if (valcrd .gt. 0) then
           read (5, valfmt, err = 998) (Value (p), p = 1, nz)
C          read (5, valfmt, err = 998) (IVAL (p), p = 1, nz)
C          do p = 1, nz
C             Value (p) = 1.0d0*IVAL(P)
C          enddo
        else
           Value (1) = myrand (0)
           do 50 p = 1, nz
              Value (p) = myrand (-1)
50         continue
        endif
*----------------------------------------------------------------------
      iatparam(50)=1
      t0=OMP_GET_WTIME()
      call OpenATI_DAFMC_CCS2CRS(iatparam,n,nz,ptr,index,value,
     & irp,icol,val)
      t1=OMP_GET_WTIME()
      WRITE(6,*) ' DATA CONVERSION CCS->CRS TIME=',SNGL(t1-t0)
C     INZ=1
C     DO I=1,N
C        IRP(I)=INZ
C        DO IC=1,N
C           DO IP=PTR(IC),PTR(IC+1)-1
C              IROW=INDEX(IP)
C              IF (I .EQ. IROW) THEN
C                 ICOL(INZ)=IC
C                 VAL(INZ)=VALUE(IP)
C                 INZ=INZ+1
C                 GO TO 100
C              END IF
C           END DO
C 100       CONTINUE
C        END DO
C     END DO
C     IRP(N+1)=INZ
*----------------------------------------------------------------------
*     OPEN(5,FILE='./memplus/memplus.rb')
*     OPEN(5,FILE='./chipcool0/chipcool0.rb')
*     OPEN(5,FILE='./wang3/wang3.rb')
*     OPEN(5,FILE='./chem_master1/chem_master1.rb')
*     OPEN(5,FILE='./epb3/epb3.rb')
*     OPEN(5,FILE='./matrix_9/matrix_9.rb')
*
*     OPEN(8,FILE='memplus.dat')
*     OPEN(8,FILE='chipcool0.dat')
*     OPEN(8,FILE='wang3.dat')
*     OPEN(8,FILE='chem_master1.dat')
*     OPEN(8,FILE='epb3.dat')
*     OPEN(8,FILE='matrix_9.dat')
      WRITE(6,*) 'N=',N,'NZ=',NZ   

      OPEN(8,FILE='AAA.dat')
*
      WRITE(8,*) N,NZ   
      WRITE(8,*) (IRP (I),I=1,N+1)  
      WRITE(8,*) (ICOL(I),I=1,NZ)  
      WRITE(8,*) (VAL (I),I=1,NZ)  
*
C     WRITE(8,*) (PTR (I),I=1,N+1)  
C     WRITE(8,*) (INDEX(I),I=1,NZ)  
C     WRITE(8,*) (VALUE (I),I=1,NZ)  
*
      CLOSE(8)
*----------------------------------------------------------------------
      STOP

  998   write (0,*) 'Read error: Harwell/Boeing matrix'
        stop

      END
c=== Myrand ============================================================
c
c  Derived from the FA01 routines in the MUPS package (CERFACS and/or
c  Harwell).  CERFACS and/or Harwell copyrights may apply.  Permission
c  granted to use this routine in the DEMO PROGRAM only.
c
c  DEMO PROGRAM.
c
c  random number generator
c  i = 0:  reinitialize the sequence
c  i >=0:  return 0 < x < 1
c  i < 0:  return -1 < x < 1

        double precision function myrand (i)
        integer i
        double precision seed, start, mfac, d2to32
        common /mrand/ seed
        parameter (start = 1431655765.d0,
     $             d2to32 = 4294967296.d0, mfac = 9228907.d0)

        if (i .eq. 0) then
c          reinitialize to known sequence
           seed = start
           endif
        seed = dmod (seed * mfac, d2to32)

        if (i .ge. 0) then
           myrand = (seed/d2to32)
        else
           myrand = 2 * (seed/d2to32) - 1
           endif
        return
        end
