      PROGRAM MAIN

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      INTEGER NMAX,NNZMAX,ISMAX,JLMAX
      PARAMETER (NMAX=300000,NNZMAX=10000000)

      INTEGER             IRP(NMAX+1),ICOL(NNZMAX)
      DOUBLE PRECISION    VAL(NNZMAX),X(NMAX),Y(NMAX)
      DOUBLE PRECISION    WK(NNZMAX)

      DOUBLE PRECISION    SINF
      ALLOCATABLE         SINF(:)

      INTEGER             INFO
      INTEGER             NUM_SMP,omp_get_max_threads

      character           FNAME*60

      INTEGER             IATPARAM(50)
      DOUBLE PRECISION    RATPARAM(50)

      CALL OPENATI_INIT(IATPARAM,RATPARAM,INFO)

C>>>>>>>>>>>MV_AT OFF
C     IATPARAM( 7) = 0
C>>>>>>>>>>>MV_AT ON
C     IATPARAM( 7) = 2
      IATPARAM( 7) = 3
      IATPARAM(50) = 1

      FNAME='../MatrixData/vibrobox.rb'

      call matread(FNAME,n,irp,icol,nnz,VAL)

      DO I=1,N
        X(I)=1.0D-3*I/N
      ENDDO
      NUM_SMP=IATPARAM(3)

      WRITE(6,2000)FNAME
 2000 FORMAT(1H ,'****  INPUT MATRIX  = ',A30)
      WRITE(6,*)'****  N ,NZ         =',N,',',NNZ
      WRITE(6,*)'****  NUM_SMP       =',NUM_SMP
      WRITE(6,*)'****  AUTO-TUNED    =',IATPARAM(7)
      WRITE(6,*)'****  DEBUG WRITE   =',IATPARAM(50)

      call BANDCHK_SYM(N,IRP,ICOL,NNZ)

      LSINF=N+NUM_SMP+3
      ALLOCATE(SINF(LSINF))

      LOOP=100
      FLOP=(4.0D0*NNZ-2.0D0*N)*(1.0D-9)

      WRITE(6,*)''
      WRITE(6,*)'***** OpenATI DSRMV SAMPLE TEST START *****'

C***** auto-tuned off
      IF(IATPARAM(7) .EQ. 0) THEN
         IATPARAM(8)=11
         CALL OPENATI_DSRMV_SETUP(N,NNZ,IRP,ICOL,IATPARAM,RATPARAM,
     &                            SINF,LSINF,INFO)
         TIME1=OMP_GET_WTIME()
         DO I=1,LOOP
            CALL OPENATI_DSRMV(N,NNZ,IRP,ICOL,VAL,X,Y,
     &                         IATPARAM,RATPARAM,WK,SINF,LSINF,INFO)
         END DO
         TIME2=OMP_GET_WTIME()
         WRITE(6,1000)IATPARAM(8),FLOP/(TIME2-TIME1)*LOOP ,TIME2-TIME1

         IATPARAM(8)=12
         CALL OPENATI_DSRMV_SETUP(N,NNZ,IRP,ICOL,IATPARAM,RATPARAM,
     &                            SINF,LSINF,INFO)
         TIME1=OMP_GET_WTIME()
         DO I=1,LOOP
            CALL OPENATI_DSRMV(N,NNZ,IRP,ICOL,VAL,X,Y,
     &                         IATPARAM,RATPARAM,WK,SINF,LSINF,INFO)
         END DO
         TIME2=OMP_GET_WTIME()
         WRITE(6,1000)IATPARAM(8),FLOP/(TIME2-TIME1)*LOOP,TIME2-TIME1

         IATPARAM(8)=13
         CALL OPENATI_DSRMV_SETUP(N,NNZ,IRP,ICOL,IATPARAM,RATPARAM,
     &                            SINF,LSINF,INFO)
         TIME1=OMP_GET_WTIME()
         DO I=1,LOOP
         CALL OPENATI_DSRMV(N,NNZ,IRP,ICOL,VAL,X,Y,
     &                      IATPARAM,RATPARAM,WK,SINF,LSINF,INFO)
         END DO
         TIME2=OMP_GET_WTIME()
         WRITE(6,1000)IATPARAM(8),FLOP/(TIME2-TIME1)*LOOP,TIME2-TIME1
      ELSE IF(IATPARAM(7) .EQ. 2 .OR. IATPARAM(7) .EQ. 3) THEN
         CALL OPENATI_DSRMV(N,NNZ,IRP,ICOL,VAL,X,Y,
     &                      IATPARAM,RATPARAM,WK,SINF,LSINF,INFO)
         WRITE(6,*)
         WRITE(6,*)'****  OpenATI DSRMV END  ****'
         WRITE(6,*)'****  FASTEST MATVEC Impl. = ',IATPARAM(8)
         WRITE(6,*)
      END IF

 1000 FORMAT(1H ,' >>>    MAT*VEC IMPL.=',I2,2X,
     &     'PERFORMANCE=',F7.3,' [GFLOPS]',F7.3,'[SEC]')

      DEALLOCATE(SINF)

      WRITE(6,*)'***** OpenATI DSRMV SAMPLE TEST END   *****'
      WRITE(6,*)''

 10   STOP
      END


      subroutine matread(filename,ncol,colptr,rowind,nnzero,
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

      READ ( LUNIT, 1000 ) TITLE , KEY   ,
     1                     TOTCRD, PTRCRD, INDCRD, VALCRD, RHSCRD,
     2                     MXTYPE, NROW  , NCOL  , NNZERO, NELTVL,
     3                     PTRFMT, INDFMT, VALFMT, RHSFMT
 1000  FORMAT ( A72, A8 / 5I14 / A3, 11X, 4I14 / 2A16, 2A20 )

C      write(6,*) '===> INPUT FILE NAME IS ',filename
C      write(6,*) TITLE
C      write(6,*) KEY
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


      SUBROUTINE BANDCHK_SYM(N,LC,LL,NZ)
      IMPLICIT REAL*8 (A-H,O-Z)
C
      integer*4 LC(N+1),LL(NZ)
C
      MAXBAND=0
      MINBAND=N+1
      XMEAN=0.0D0
C
      DO 10 I=1,N
        LOCAL=0
        DO 20 J=LC(I)+1,LC(I+1)-1
          NEWJ=LL(J)
          IF (NEWJ .GT. I) then
            MAXBAND=MAX(MAXBAND,NEWJ-I)
            MINBAND=MIN(MINBAND,NEWJ-I)
            LOCAL  =MAX(LOCAL  ,NEWJ-I)
          END IF
 20          CONTINUE
        XMEAN=XMEAN+LOCAL
 10      CONTINUE
      XMEAN=XMEAN/n
      MEAN=NINT(XMEAN)
C
      IC=0
      MAXIC=0
      MINIC=N+1
      DO 30 I=1,N
         INUM=LC(I+1)-LC(I)
         IC=IC+INUM
         MAXIC=MAX(MAXIC,INUM)
         MINIC=MIN(MINIC,INUM)
 30       CONTINUE
      RC=1.0D0*IC/N
C
      WRITE(6,1000)MEAN,MAXBAND,MINBAND
 1000  FORMAT(1h ,'****  BAND INFO.        (MEAN,MAX,MIN)=',
     $            I12,',',I12,',',I12)
      WRITE(6,1100)RC,MAXIC,MINIC
 1100  FORMAT(1h ,'****  NUM. OF NON. ZERO (MEAN,MAX,MIN)=',
     $            F12.5,',',I12,',',I12)

      return
      end
