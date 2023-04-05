      SUBROUTINE OpenATI_DAFMC_CCS2CRS
     &        (IATPARAM,N,NNZ,IPTR,INDEX,VALUE,IRP,ICOL,VAL)
*
*     OpenATI_DAFMC_CCS2CRS : Convert sparse matrix format. CCS to CRS
*
*     December 2011
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER          IPTR(N+1),INDEX(NNZ),IRP(N+1),ICOL(NNZ)
      DOUBLE PRECISION VALUE(NNZ),VAL(NNZ)
      INTEGER          IATPARAM(50)
*
      integer,allocatable::iord(:),mark(:),iwk(:)
      double precision,allocatable::rwk(:)
*
*
*  =====================================================================
*  Purpose
*  =======
*
*   OpenATI_DAFMC_CCS2CRS : Convert sparse matrix format. CCS to CRS
*
*
*  Arguments
*  =========
*  N       (input) INTEGER
*          The order of the matrix.  N >= 0.
*
*  NNZ     (input) INTEGER
*          Non-Zeros of the matrix.  NNZ >= N.
*
*  IPTR    (input) INTEGER array
*  (N+1)   Diagonal pointer of the matrix in CCS format.
*
*  INDEX   (input) INTEGER array
*  (NNZ)   Column indecies of the matrix in CCS format.
*
*  VALUE   (input) DOUBLE PRECISION array
*  (NNZ)   Value of the matrix in CCS format.
*
*  IRP     (output) INTEGER array
*  (N+1)   Diagonal pointer of the matrix in CRS format.
*
*  ICOL    (output) INTEGER array
*  (NNZ)   Column indecies of the matrix in CRS format.
*
*  VAL     (output) DOUBLE PRECISION array
*  (NNZ)   Value of the matrix in CRS format.
*
*  Local variables
*  ===============
*
*  =====================================================================
*>>>>>>>>>>>>>>>>Get Debug print Env.
      IDBG = IATPARAM(50)
*
      IF (IDBG.EQ.1) THEN
         WRITE(6,*) ' [OpenATI] DATA CONVERSION CCS->CRS START'
      END IF
C     INZ=1
C     DO I=1,N
C        IRP(I)=INZ
C        IF (IDBG.EQ.1) THEN
C          IF (MOD(I,1000).EQ.0) THEN
C             WRITE(6,*) '  [OpenATI] TRANSLATE ROW=',I,' LAST=',N
C          END IF
C        END IF
C        DO IC=1,N
C           DO IP=IPTR(IC),IPTR(IC+1)-1
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
C     return
*+++++++++++++++++++++++++++
      allocate(mark(n))
      mark(1:n)=0
      do i=1,n
         do k=iptr(i),iptr(i+1)-1
            j=index(k)
            mark(j)=mark(j)+1
         end do
      end do
      irp(1)=1
      do i=1,n
         irp(i+1)=irp(i)+mark(i)
      end do
      mark(1:n)=0
      do i=1,n
         do k=iptr(i),iptr(i+1)-1
            j=index(k)
            mark(j)=mark(j)+1
            icol(irp(j)+mark(j)-1)=i
            val (irp(j)+mark(j)-1)=value(k)
         end do
      end do
*
!$omp parallel default(none)
!$omp+ shared(n,irp,icol,val)
!$omp+ private(i,nl,il,iord,iwk,rwk)
      allocate(iord(n),iwk(n),rwk(n))
!$omp do
      do i=1,n
         nl=irp(i+1)-irp(i)
         call openati_qsorti(iord,nl,icol(irp(i)))
         do il=1,nl
            iwk(il)=icol(irp(i)-1+iord(il))
            rwk(il)=val (irp(i)-1+iord(il))
         enddo
         do il=1,nl
            icol(irp(i)-1+il)=iwk(il)
            val (irp(i)-1+il)=rwk(il)
         enddo
      enddo
      if (allocated(iord)) deallocate(iord)
      if (allocated(iwk)) deallocate(iwk)
      if (allocated(rwk)) deallocate(rwk)
!$omp end parallel
*+++++++++++++++++++++++++++
*+++++++++++++++++++++++++++
*+++++++++++++++++++++++++++
*----------------------------------------------------------------------
      IF (IDBG.EQ.1) THEN
         WRITE(6,*) ' [OpenATI] DATA CONVERSION CCS->CRS FINISHED'
      END IF
*----------------------------------------------------------------------
      if (allocated(mark)) deallocate(mark)
      RETURN
      END
