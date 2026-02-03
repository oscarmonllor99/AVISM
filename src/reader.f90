!************************************************************************
SUBROUTINE NOMFILE_MASCLET(ITER,FILNOM1,FILNOM2,FILNOM3, FILNOM4)
!************************************************************************

       IMPLICIT NONE
       INTEGER ITER
       CHARACTER*9 FILNOM1,FILNOM2,FILNOM4
       CHARACTER*10 FILNOM3
       CHARACTER*5 NOM
       INTEGER CONTA,I,N10,IT

       CONTA=0

       DO I=4,0,-1
         N10=10**I
         IT=ITER/N10 - CONTA
         CONTA=(CONTA+IT)*10
         NOM(5-I:5-I)=CHAR(48+IT)
       END DO

       FILNOM1='clus'//NOM
       FILNOM2='cldm'//NOM
       FILNOM3='grids'//NOM
       FILNOM4='clst'//NOM

!*************************************************************************
END SUBROUTINE NOMFILE_MASCLET
!*************************************************************************


!************************************************************************
SUBROUTINE NOMFILE_BIN_PART(ITER,FILNOM)
!************************************************************************

       IMPLICIT NONE
       INTEGER ITER
       CHARACTER*50 FILNOM
       CHARACTER*5 NOM
       INTEGER CONTA,I,N10,IT

       CONTA=0

       DO I=4,0,-1
         N10=10**I
         IT=ITER/N10 - CONTA
         CONTA=(CONTA+IT)*10
         NOM(5-I:5-I)=CHAR(48+IT)
       END DO

       FILNOM='bin_file_part'//NOM
       FILNOM=TRIM(FILNOM)

!*************************************************************************
END SUBROUTINE  NOMFILE_BIN_PART
!*************************************************************************

!************************************************************************
SUBROUTINE NOMFILE_BIN_GRID(ITER,FILNOM)
!************************************************************************

       IMPLICIT NONE
       INTEGER ITER
       CHARACTER*50 FILNOM
       CHARACTER*5 NOM
       INTEGER CONTA,I,N10,IT

       CONTA=0

       DO I=4,0,-1
         N10=10**I
         IT=ITER/N10 - CONTA
         CONTA=(CONTA+IT)*10
         NOM(5-I:5-I)=CHAR(48+IT)
       END DO

       FILNOM='bin_file_grid'//NOM
       FILNOM=TRIM(FILNOM)

!*************************************************************************
END SUBROUTINE  NOMFILE_BIN_GRID
!*************************************************************************



!*************************************************************************
SUBROUTINE READ_MASCLET(ITER, NDXYZ0, NL, T, ZETA) !output
!*************************************************************************

      USE COMMONDATA
      IMPLICIT NONE

! input variables
       INTEGER ITER

! local variables
       INTEGER IR,  IRR, LOW1,LOW2, N1, N2, N3
       INTEGER I, IX, J, K
       INTEGER CONTA, PARTIBAS
       INTEGER:: DIM1, DIM2, DIM3, DIM4, NHYXP, NHYYP, NHYZP, FLAG_DENS_P, PARTIREDP
       REAL*4 AAA,BBB,CCC, MAP
       REAL*4, ALLOCATABLE::UBAS(:)
       INTEGER,ALLOCATABLE::UBAS2(:)
       CHARACTER*9 FILNOM1,FILNOM2, FILNOM4
       CHARACTER*10 FILNOM3
       CHARACTER*25 FIL1,FIL2, FIL4
       CHARACTER*26 FIL3

! output variables
       INTEGER NL, NDXYZ0, NST0
       REAL*4 T,ZETA


!* ALLOCATE AMR VARIABLES

       ALLOCATE(NPATCH(0:NLEVELS),PARE(NPALEV))
       ALLOCATE(NPART(0:NLEVELS),NPARTST(0:NLEVELS))
       ALLOCATE(PATCHNX(NPALEV),PATCHNY(NPALEV),PATCHNZ(NPALEV))
       ALLOCATE(PATCHX(NPALEV),PATCHY(NPALEV),PATCHZ(NPALEV))
       ALLOCATE(PATCHRX(NPALEV),PATCHRY(NPALEV),PATCHRZ(NPALEV))

       PATCHNX=0
       PATCHNY=0
       PATCHNZ=0
       PATCHX=0
       PATCHY=0
       PATCHZ=0
       PATCHRX=0.0
       PATCHRY=0.0
       PATCHRZ=0.0
       NPATCH=0
       PARE=0

!*      READING DATA
       CALL NOMFILE_MASCLET(ITER,FILNOM1,FILNOM2,FILNOM3,FILNOM4)

       !WRITE(*,*) 'Reading iter',ITER,' ',FILNOM1,FILNOM2,FILNOM3

       FIL1='simu_masclet/'//FILNOM1 ! gas
       FIL2='simu_masclet/'//FILNOM2 ! dm
       FIL3='simu_masclet/'//FILNOM3 ! grid
       FIL4='simu_masclet/'//FILNOM4 ! stars

       OPEN (33,FILE=FIL3,STATUS='UNKNOWN',ACTION='READ')
       OPEN (31,FILE=FIL1, &
            STATUS='UNKNOWN',ACTION='READ',FORM='UNFORMATTED')
       OPEN (32,FILE=FIL2, &
            STATUS='UNKNOWN',ACTION='READ',FORM='UNFORMATTED')

       ! STARS NOT CONSIDERED
       !
       ! OPEN (34,FILE=FIL4, &
       !      STATUS='UNKNOWN',ACTION='READ',FORM='UNFORMATTED')

!*      GRID DATA
       READ(33,*) IRR,T,NL,MAP
       READ(33,*) ZETA
       READ(33,*) IR,NDXYZ0,NST0 
       WRITE(*,*) '     ZETA,NL,NDXYZ0,NST0,MAP', ZETA,NL,NDXYZ0,NST0,MAP

       DO IR=1,NL
       READ(33,*) IRR,NPATCH(IR),NPART(IR),NPARTST(IR)
       READ(33,*)

       IF (IR.NE.IRR) WRITE(*,*)'     Warning: fail in restart'
       LOW1=SUM(NPATCH(0:IR-1))+1
       LOW2=SUM(NPATCH(0:IR))
       DO I=LOW1,LOW2
        READ(33,*) PATCHNX(I),PATCHNY(I),PATCHNZ(I)
        READ(33,*) PATCHX(I),PATCHY(I),PATCHZ(I)
        READ(33,*) AAA,BBB,CCC
        PATCHRX(I)=AAA
        PATCHRY(I)=BBB
        PATCHRZ(I)=CCC
        READ(33,*) PARE(I)
       END DO
       END DO
       CLOSE(33)
       NPART(0)=NDXYZ0
       NPARTST(0)=NST0

       WRITE(*,*) '     NPART (LEV0, LEV>0, ALL):',NDXYZ0, SUM(NPART(1:NL)), SUM(NPART)

       NPARTT=INT(SUM(NPART(0:NL)), KIND=8)

       IF(NPARTT.GT.PARTIRED) THEN
         WRITE(*,*) '     ERROR: NPARTT > PARTIRED'
         STOP
       ENDIF
       
!***************** ALLOCATE AND INITIALIZE ARRAYS *******************************************

   !l=0 VARIABLES
       PARTIREDP=PARTIRED

       ALLOCATE(U1G(NHYX,NHYY,NHYZ), U1S(NHYX,NHYY,NHYZ),U1DM(NHYX,NHYY,NHYZ))
       ALLOCATE(U2G(NHYX,NHYY,NHYZ), U3G(NHYX,NHYY,NHYZ), U4G(NHYX,NHYY,NHYZ))
       ALLOCATE(U2PA(PARTIREDP), U3PA(PARTIREDP), U4PA(PARTIREDP))
       ALLOCATE(MASAP(PARTIREDP), RXPA(PARTIREDP), RYPA(PARTIREDP), RZPA(PARTIREDP))
       
       NHYXP=NHYX
       NHYYP=NHYY
       NHYZP=NHYZ
!$OMP PARALLEL DO SHARED(NHYXP,NHYYP,NHYZP,U1DM,U1G,U1S), &
!$OMP             PRIVATE(I,J,K)
      DO K=1,NHYZP
      DO J=1,NHYYP
      DO I=1,NHYXP
       U1DM(I,J,K)=-1.0 !min value for the density contrast
       U1G(I,J,K)=-1.0
       U1S(I,J,K)=-1.0
       U2G(I,J,K)=0.0
       U3G(I,J,K)=0.0
       U4G(I,J,K)=0.0
      END DO
      END DO
      END DO

      
!$OMP PARALLEL DO SHARED(PARTIREDP,U2PA,U3PA,U4PA,RXPA,RYPA,RZPA, &
!$OMP            MASAP), &
!$OMP            PRIVATE(I)
      DO I=1,PARTIREDP
       U2PA(I)=0.0 
       U3PA(I)=0.0
       U4PA(I)=0.0
       RXPA(I)=0.0  !DM vars
       RYPA(I)=0.0
       RZPA(I)=0.0
       MASAP(I)=0.0
      END DO

      !l>0 (AMR) VARIABLES

       DIM1=MAXVAL(PATCHNX)
       DIM2=MAXVAL(PATCHNY)
       DIM3=MAXVAL(PATCHNZ)
       DIM4=SUM(NPATCH(0:NL2))

       ALLOCATE(U11G(DIM1, DIM2, DIM3, DIM4))
       ALLOCATE(U12G(DIM1, DIM2, DIM3, DIM4))
       ALLOCATE(U13G(DIM1, DIM2, DIM3, DIM4))
       ALLOCATE(U14G(DIM1, DIM2, DIM3, DIM4))
       ALLOCATE(U11DM(DIM1, DIM2, DIM3, DIM4))
       ALLOCATE(U11S(DIM1, DIM2, DIM3, DIM4))

!$OMP PARALLEL DO SHARED(DIM1, DIM2, DIM3, DIM4, U11DM, U11G, U11S), PRIVATE(IX, I,J,K)
       DO I=1,DIM4
        DO K=1,DIM3
        DO J=1,DIM2
        DO IX=1,DIM1
         !matter components density contrast
         U11DM(IX,J,K,I)=-1.0
         U11G(IX,J,K,I)=-1.0
         U11S(IX,J,K,I)=-1.0
        END DO
        END DO
        END DO
       END DO

!************************************************************************

!*      GAS
       READ(31)
       IR=0
        N1=NHYX
        N2=NHYY
        N3=NHYZ
        READ(31) (((U1G(I,J,K),I=1,N1),J=1,N2),K=1,N3)
        READ(31) (((U2G(I,J,K),I=1,N1),J=1,N2),K=1,N3)
        READ(31) (((U3G(I,J,K),I=1,N1),J=1,N2),K=1,N3)
        READ(31) (((U4G(I,J,K),I=1,N1),J=1,N2),K=1,N3)
        READ(31) 
        READ(31) 
        READ(31) 
        !IF COOLING ACTIVATED
        READ(31) 
        READ(31) 
        !!!!!!!!
        READ(31) 
        !IF MAGNETIC FIELD CONSIDERED
        !READ(31)
        !READ(31)
        !READ(31)
      
       !gas is read only to the LEVMAX level
       DO IR=1,NL2
       LOW1=SUM(NPATCH(0:IR-1))+1
       LOW2=SUM(NPATCH(0:IR))
       DO I=LOW1,LOW2
        N1=PATCHNX(I)
        N2=PATCHNY(I)
        N3=PATCHNZ(I)
        IF(N1 .GT. DIM1 .OR. N2 .GT. DIM2 .OR. N3 .GT. DIM3 .OR. I .GT. DIM4) THEN
           WRITE(*,*) '     Bad dimension for U11G', N1,N2,N3,I
           STOP
        ENDIF
        READ(31) (((U11G(IX,J,K,I),IX=1,N1),J=1,N2),K=1,N3)
        READ(31) (((U12G(IX,J,K,I),IX=1,N1),J=1,N2),K=1,N3)
        READ(31) (((U13G(IX,J,K,I),IX=1,N1),J=1,N2),K=1,N3)
        READ(31) (((U14G(IX,J,K,I),IX=1,N1),J=1,N2),K=1,N3)
        READ(31) 
        READ(31)
        READ(31) 
        !IF COOLING ACTIVATED
        READ(31) 
        READ(31) 
        !!!!!!!!
        READ(31) 
        READ(31) 
        !IF MAGNETIC FIELD CONSIDERED
        !READ(31)
        !READ(31)
        !READ(31)
       END DO
       END DO

      CLOSE(31)
      write(*,*) '     CLUS read!'

!**     DARK MATTER
       READ(32)
       IR=0
        N1=NHYX
        N2=NHYY
        N3=NHYZ

        READ(32) (((U1DM(I,J,K),I=1,N1),J=1,N2),K=1,N3)
        READ(32) (RXPA(I),I=1,NDXYZ0)
        READ(32) (RYPA(I),I=1,NDXYZ0)
        READ(32) (RZPA(I),I=1,NDXYZ0)
        READ(32) (U2PA(I),I=1,NDXYZ0)
        READ(32) (U3PA(I),I=1,NDXYZ0)
        READ(32) (U4PA(I),I=1,NDXYZ0)
        READ(32) 
        CONTA=NDXYZ0
        MASAP(1:NDXYZ0)=MAP

        !WRITE(*,*) 'NPART(IR)=',IR, NPART(IR),CONTA

       PARTIBAS=MAXVAL(NPART)
       ALLOCATE(UBAS(PARTIBAS))
       ALLOCATE(UBAS2(PARTIBAS))

       !dark matter is read for all levels (since particles lie in all grid levels)
       !except for DM delta, which is like gas, read only to the LEVMAX level
       DO IR=1,NL
       LOW1=SUM(NPATCH(0:IR-1))+1
       LOW2=SUM(NPATCH(0:IR))
       DO I=LOW1,LOW2
        N1=PATCHNX(I)
        N2=PATCHNY(I)
        N3=PATCHNZ(I)
        IF (IR .LE. NL2) THEN
           READ(32) (((U11DM(IX,J,K,I),IX=1,N1),J=1,N2),K=1,N3)
        ELSE 
           READ(32) 
        ENDIF
       END DO

        UBAS=0.0
        UBAS2=0
        READ(32) (UBAS(IX),IX=1,NPART(IR))
        RXPA(CONTA+1:CONTA+NPART(IR))=UBAS(1:NPART(IR))
        READ(32) (UBAS(IX),IX=1,NPART(IR))
        RYPA(CONTA+1:CONTA+NPART(IR))=UBAS(1:NPART(IR))
        READ(32) (UBAS(IX),IX=1,NPART(IR))
        RZPA(CONTA+1:CONTA+NPART(IR))=UBAS(1:NPART(IR))
        READ(32) (UBAS(IX),IX=1,NPART(IR))
        U2PA(CONTA+1:CONTA+NPART(IR))=UBAS(1:NPART(IR))
        READ(32) (UBAS(IX),IX=1,NPART(IR))
        U3PA(CONTA+1:CONTA+NPART(IR))=UBAS(1:NPART(IR))
        READ(32) (UBAS(IX),IX=1,NPART(IR))
        U4PA(CONTA+1:CONTA+NPART(IR))=UBAS(1:NPART(IR))
        READ(32) (UBAS(IX),IX=1,NPART(IR))
        MASAP(CONTA+1:CONTA+NPART(IR))=UBAS(1:NPART(IR))
        READ(32) 

        IF (NPART(IR).GT.0) THEN
        !WRITE(*,*) 'ORIPA=',MAXVAL(ORIPA(CONTA+1:CONTA+NPART(IR))), &
        !                    MINVAL(ORIPA(CONTA+1:CONTA+NPART(IR)))
        END IF

        CONTA=CONTA+NPART(IR)
        WRITE(*,*) '     NPART(IR)=',IR,NPART(IR),CONTA

        END DO
       CLOSE(32)
       WRITE(*,*) '     CLDM read!'

        !*     delta to density
        U1DM=U1DM+1.
        U1G=U1G+1.
        U1S=U1S+1.
        U11DM=U11DM+1.
        U11G=U11G+1.
        U11S=U11S+1.
       !------------------------------------------

        DEALLOCATE(UBAS, UBAS2)

        CONTA=SUM(NPART(0:NL))
        WRITE(*,*) '     TOTAL PARTICLES IN ITER=', NPARTT

!*************************************************************************
END SUBROUTINE READ_MASCLET
!*************************************************************************


!*************************************************************************
SUBROUTINE READ_BINARY_PART_1(ITER, ZETA)
!READS BINARY FILE WITH ALL AVAILABLE PARTICLE DATA
!Input is assumed to be in: cMpc, km/s, Msun
!*************************************************************************
      USE COMMONDATA
      IMPLICIT NONE
      CHARACTER*50 FILENAME
      INTEGER ITER
      INTEGER(KIND=8) I
      REAL*4 ZETA
      
      !ALLOCATE PARTICLE ARRAYS
      ALLOCATE(U2PA(PARTIRED), U3PA(PARTIRED), U4PA(PARTIRED))
      ALLOCATE(MASAP(PARTIRED), RXPA(PARTIRED), RYPA(PARTIRED), RZPA(PARTIRED))

      !INITIALIZE PARTICLE ARRAYS
      !$OMP PARALLEL DO SHARED(PARTIRED,U2PA,U3PA,U4PA,RXPA,RYPA,RZPA, &
      !$OMP            MASAP), &
      !$OMP            PRIVATE(I)
      DO I=1,PARTIRED
       U2PA(I)=0.0 
       U3PA(I)=0.0
       U4PA(I)=0.0
       RXPA(I)=0.0  !DM vars
       RYPA(I)=0.0
       RZPA(I)=0.0
       MASAP(I)=0.0
      END DO

      CALL NOMFILE_BIN_PART(ITER, FILENAME)

      OPEN(UNIT=41, FILE='input_data/' // FILENAME, FORM = 'UNFORMATTED')
       READ(41) NPARTT, ZETA
       IF(NPARTT.GT.PARTIRED) THEN
         WRITE(*,*) '     ERROR: NPARTT > PARTIRED'
         STOP 
       ENDIF
       READ(41) (RXPA(I), I=1,NPARTT)
       READ(41) (RYPA(I), I=1,NPARTT)
       READ(41) (RZPA(I), I=1,NPARTT)
       READ(41) (U2PA(I), I=1,NPARTT)
       READ(41) (U3PA(I), I=1,NPARTT)
       READ(41) (U4PA(I), I=1,NPARTT)
       READ(41) (MASAP(I), I=1,NPARTT)
      CLOSE(41)

      !to masclet units
      MASAP = MASAP/UM
      U2PA = U2PA/UV
      U3PA = U3PA/UV
      U4PA = U4PA/UV

      WRITE(*,*) '     TOTAL PARTICLES IN ITER=', NPARTT

!************************************************************************
END SUBROUTINE READ_BINARY_PART_1
!************************************************************************

!*************************************************************************
SUBROUTINE READ_BINARY_PART_2(ITER, ZETA)
!SAME AS READ_BINARY_PART_1 BUT WITHOUT VELOCITIES
!Vel. Divergence will be reconstructed from the density field
!Input is assumed to be in: cMpc, Msun
!*************************************************************************
      USE COMMONDATA
      IMPLICIT NONE
      CHARACTER*50 FILENAME
      INTEGER ITER
      INTEGER(KIND=8) I
      REAL*4 ZETA
      
      !ALLOCATE PARTICLE ARRAYS
      ALLOCATE(U2PA(PARTIRED), U3PA(PARTIRED), U4PA(PARTIRED))
      ALLOCATE(MASAP(PARTIRED), RXPA(PARTIRED), RYPA(PARTIRED), RZPA(PARTIRED))

      !INITIALIZE PARTICLE ARRAYS
      !$OMP PARALLEL DO SHARED(PARTIRED,U2PA,U3PA,U4PA,RXPA,RYPA,RZPA, &
      !$OMP            MASAP), &
      !$OMP            PRIVATE(I)
      DO I=1,PARTIRED
       U2PA(I)=0.0 
       U3PA(I)=0.0
       U4PA(I)=0.0
       RXPA(I)=0.0  !DM vars
       RYPA(I)=0.0
       RZPA(I)=0.0
       MASAP(I)=0.0
      END DO

      CALL NOMFILE_BIN_PART(ITER, FILENAME)

      OPEN(UNIT=41, FILE='input_data/' // FILENAME, FORM = 'UNFORMATTED')
       READ(41) NPARTT, ZETA
       IF(NPARTT.GT.PARTIRED) THEN
         WRITE(*,*) '     ERROR: NPARTT > PARTIRED'
         STOP 
       ENDIF
       READ(41) (RXPA(I), I=1,NPARTT)
       READ(41) (RYPA(I), I=1,NPARTT)
       READ(41) (RZPA(I), I=1,NPARTT)
       READ(41) (MASAP(I), I=1,NPARTT)
      CLOSE(41)

      !to masclet units
      MASAP = MASAP/UM

      WRITE(*,*) '     TOTAL PARTICLES IN ITER=', NPARTT
      
!************************************************************************
END SUBROUTINE READ_BINARY_PART_2
!************************************************************************

!*************************************************************************
SUBROUTINE READ_BINARY_GRID(ITER, ZETA)
!*************************************************************************
      USE COMMONDATA
      IMPLICIT NONE
      CHARACTER*50 FILENAME
      INTEGER I,J,K, ITER
      REAL*4 ZETA

      !ALLOCATE PARTICLE ARRAYS
      ALLOCATE(U1G(NHYX,NHYY,NHYZ), U2G(NHYX,NHYY,NHYZ), U3G(NHYX,NHYY,NHYZ), U4G(NHYX,NHYY,NHYZ))

      CALL NOMFILE_BIN_GRID(ITER, FILENAME)

      OPEN(UNIT=41, FILE='input_data/' // FILENAME, FORM = 'UNFORMATTED')
       READ(41) ZETA
       READ(41) (((U1G(I,J,K),I=1,NHYX),J=1,NHYY),K=1,NHYZ)
       READ(41) (((U2G(I,J,K),I=1,NHYX),J=1,NHYY),K=1,NHYZ)
       READ(41) (((U3G(I,J,K),I=1,NHYX),J=1,NHYY),K=1,NHYZ)
       READ(41) (((U4G(I,J,K),I=1,NHYX),J=1,NHYY),K=1,NHYZ)
      CLOSE(41)

      !delta to density
       U1G = U1G + 1.0
      !v in c units, from km/s to c
       U2G = U2G/UV
       U3G = U3G/UV
       U4G = U4G/UV

!************************************************************************
END SUBROUTINE READ_BINARY_GRID
!************************************************************************

