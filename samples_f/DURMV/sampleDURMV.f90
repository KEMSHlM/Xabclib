      PROGRAM MAIN

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)     

      INTEGER           NMAX,NNZMAX,JL
      PARAMETER         (NMAX=300000,NNZMAX=10000000)

      INTEGER           IRP(NMAX),ICOL(NNZMAX)
      DOUBLE PRECISION  VAL(NNZMAX),X(NMAX),Y(NMAX)

      DOUBLE PRECISION  UINF
      ALLOCATABLE       UINF(:)

      INTEGER           ICASE,INFO,LUINF
      INTEGER           NUM_SMP,OMP_GET_MAX_THREADS

      CHARACTER*60      FNAME

      INTEGER           IATPARAM(50)
      DOUBLE PRECISION  RATPARAM(50)

      CALL OpenATI_INIT(IATPARAM,RATPARAM,INFO)

C>>>>>>>>>>>MV_AT OFF
C     IATPARAM( 9) = 0
C>>>>>>>>>>>MV_AT ON
C     IATPARAM( 9) = 2
      IATPARAM( 9) = 3
      IATPARAM(50) = 1

      JL=128
      IATPARAM(11) = JL

      FNAME="../MatrixData/ecl32.dat"

      OPEN(32,FILE=FNAME)
      READ(32,*) N,NNZ
      READ(32,*) (IRP (I),I=1,N+1)
      READ(32,*) (ICOL(I),I=1,NNZ)
      READ(32,*) (VAL (I),I=1,NNZ)
      CLOSE(32)

      DO I=1,N
        X(I)=1.0D-3*I/N
      ENDDO
      NUM_SMP=OMP_GET_MAX_THREADS()

      WRITE(6,2000)FNAME
 2000 FORMAT(1H ,'****  INPUT MATRIX  = ',A30)
      WRITE(6,*)'****  N ,NZ         =',N,',',NNZ
      WRITE(6,*)'****  NUM_SMP       =',NUM_SMP
      WRITE(6,*)'****  AUTO-TUNED    =',IATPARAM(9)
      WRITE(6,*)'****  DEBUG WRITE   =',IATPARAM(50)

      CALL BANDCHK(N,IRP,ICOL,NNZ)

      LUINF=INT(1.125D0*NNZ)+INT(2.125D0*JL)+10
      ALLOCATE(UINF(LUINF))

      LOOP=100
      FLOP=2.0D0*NNZ*1.0D-9

      WRITE(6,*)''
      WRITE(6,*)'***** OpenATI DURMV SAMPLE TEST START *****'


C***** auto-tuned off
      IF(IATPARAM(9) .EQ. 0) THEN
         IATPARAM(10) = 11
         CALL OpenATI_DURMV(N,NNZ,IRP,ICOL,VAL,X,Y,
     $                      IATPARAM,RATPARAM,UINF,LUINF,INFO)
         TIME1=OMP_GET_WTIME()
         DO I=1,LOOP
            CALL OpenATI_DURMV(N,NNZ,IRP,ICOL,VAL,X,Y,
     $                         IATPARAM,RATPARAM,UINF,LUINF,INFO)
         END DO
         TIME2=OMP_GET_WTIME()
         WRITE(6,1000)ICASE,FLOP*LOOP/(TIME2-TIME1),TIME2-TIME1

         IATPARAM(10) = 12
         CALL OpenATI_DURMV_Setup(N,NNZ,IRP,IATPARAM,RATPARAM,
     &                            UINF,LUINF,INFO)
         TIME1=OMP_GET_WTIME()
         DO I=1,LOOP
            CALL OpenATI_DURMV(N,NNZ,IRP,ICOL,VAL,X,Y,
     $                         IATPARAM,RATPARAM,UINF,LUINF,INFO)
         END DO
         TIME2=OMP_GET_WTIME()
         WRITE(6,1000)ICASE,FLOP*LOOP/(TIME2-TIME1),TIME2-TIME1

         IATPARAM(10) = 13
         CALL OpenATI_DURMV_Setup(N,NNZ,IRP,IATPARAM,RATPARAM,
     &                            UINF,LUINF,INFO)
         TIME1=OMP_GET_WTIME()
         DO I=1,LOOP
            CALL OpenATI_DURMV(N,NNZ,IRP,ICOL,VAL,X,Y,
     $                         IATPARAM,RATPARAM,UINF,LUINF,INFO)
         END DO
         TIME2=OMP_GET_WTIME()
         WRITE(6,1010)ICASE,FLOP*LOOP/(TIME2-TIME1),TIME2-TIME1,JL
      ELSE IF(IATPARAM(9) .EQ. 2 .OR. IATPARAM(9) .EQ. 3) THEN
         CALL OpenATI_DURMV(N,NNZ,IRP,ICOL,VAL,X,Y,
     $                      IATPARAM,RATPARAM,UINF,LUINF,INFO)
         WRITE(6,*)
         WRITE(6,*)'****  OpenATI DURMV END  ****'
         WRITE(6,*)'****  FASTEST MATVEC Impl. = ',IATPARAM(10)
         WRITE(6,*)
      ENDIF

 1000 FORMAT(1H ,' >>>    Mat*Vec Impl.=',I2,2X,
     $     'PERFORMANCE=',F7.3,' [Gflops]',F7.3,'[sec]')
 1010 FORMAT(1H ,' >>>    Mat*Vec Impl.=',I2,2X,
     $     'PERFORMANCE=',F7.3,' [Gflops]',F7.3,'[sec] in JL=',I4)

      DEALLOCATE(UINF)

      WRITE(6,*)'***** OpenATI DURMV SAMPLE TEST END   *****'
      WRITE(6,*)''

 9999 STOP
      END


      SUBROUTINE BANDCHK(N,LC,LL,NZ)
      IMPLICIT REAL*8 (A-H,O-Z)
C
      integer*4 LC(N+1),LL(NZ)
C
      MAXBAND=0
      MINBAND=N+1
      XMEAN=0.0D0
C
      DO 10 I=1,N
        istart=LL(LC(I))
        iend  =LL(LC(I+1)-1)
        iband =iend-istart
        maxband=max(maxband,iband)
        minband=min(minband,iband)
        xmean=xmean+iband
   10 CONTINUE
      xmean=xmean/n
      mean=nint(xmean)
C
      IC=0
      MAXIC=0
      MINIC=N+1
      DO 30 I=1,N
         INUM=LC(I+1)-LC(I)
         IC=IC+INUM
         MAXIC=MAX(MAXIC,INUM)
         MINIC=MIN(MINIC,INUM)
   30 CONTINUE
      RC=1.0D0*IC/N

      WRITE(6,1000)MEAN,MAXBAND,MINBAND
 1000 FORMAT(1h ,'****  BAND INFO.        (MEAN,MAX,MIN)=',
     $            I12,',',I12,',',I12)
      WRITE(6,1100)RC,MAXIC,MINIC
 1100 FORMAT(1h ,'****  NUM. OF NON. ZERO (MEAN,MAX,MIN)=',
     $            F12.5,',',I12,',',I12)
      return
      end
