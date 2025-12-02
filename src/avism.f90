!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! ABOUT:
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Author: Óscar Monllor-Berbegal 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Computational Cosmology Group, 
! Departament d'Astronomia i Astrofísica,
! Universitat de València
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!-------------------------------------------------------
!-------------------------------------------------------
!-------------------------------------------------------
!-------------------------------------------------------
       PROGRAM AVISM
!-------------------------------------------------------
!-------------------------------------------------------
!-------------------------------------------------------
!------------------------------------------------------
       USE COMMONDATA
       USE COSMOKDTREE
       USE PARTICLES
       USE VOIDFINDING
       
       IMPLICIT NONE

       ! ITER 
       INTEGER :: IFI, NITER, NDXYZ0
       INTEGER :: ITER, CONTA
       INTEGER :: FILES_PER_SNAP, PARTTYPEX
       REAL*4 :: MASSDM
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       
       ! COSMO
       REAL*4 :: UNTERCIO
       real*4 :: RODO,T0,RE0,ACHE,ACHE0,OMEGA0,ROCRIT
       REAL*4 :: T,ZETA
       REAL*4 :: RETE,ROTE,HTE
       REAL*4 :: MEANDENS, TOTALMASS
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       
       ! GRID
       REAL :: VCELL
       REAL*4 :: LADO
       INTEGER :: NXX,NYY,NZZ,NX2,NY2,NZ2
       REAL*4 :: DXX, DYY, DZZ
       INTEGER :: NL, IR
       INTEGER :: IX,JY,KZ,IXCO,JYCO,KZCO
       REAL*4 :: RX1, RY1, RZ1
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       
       ! VOIDS
       INTEGER :: NVOID, NVOIDP, NVOIDT, NVOID2
       INTEGER :: I, IV, IND, INDP
       REAL*8 :: VOLT_CLEAN, VOLM
       REAL*4 :: REQPP

       INTEGER, ALLOCATABLE, DIMENSION(:) :: INDICE2, INDICE
       REAL*4,  ALLOCATABLE, DIMENSION(:) :: VOL2
       INTEGER, ALLOCATABLE, DIMENSION(:) :: UVOID
       REAL*4,  ALLOCATABLE, DIMENSION(:) :: VOLNEW, UMEAN, MTOT, REQ, REQP
       REAL*4,  ALLOCATABLE, DIMENSION(:) :: XC, YC, ZC, GXC, GYC, GZC
       REAL*4,  ALLOCATABLE, DIMENSION(:) :: EPS, IP
       INTEGER, ALLOCATABLE, DIMENSION(:) :: NCELLV
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       ! TIMING
       integer*4 t1,t2,trate,tmax,clock_t0,clock_tf
       INTEGER DATE(3), TIME(3)
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       ! FILE
       CHARACTER*50 DIR, FILEO, FILEO1, FILEO2, FILEO3
       INTEGER :: exit_stat, cmd_stat
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       
       ! OMP
       INTEGER:: NUM, OMP_GET_NUM_THREADS
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       !MISC
       INTEGER :: FLAG_STOP
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !  !!!!!!!!!!!!!!
      !  !Void deconstruction variables
      !  INTEGER :: dNVOID
      !  INTEGER, ALLOCATABLE, DIMENSION(:) :: dUVOID
      !  INTEGER, ALLOCATABLE, DIMENSION(:) :: dINICIOX,dINICIOY,dINICIOZ
      !  INTEGER, ALLOCATABLE, DIMENSION(:) :: dFINALX,dFINALY,dFINALZ
      !  REAL, ALLOCATABLE, DIMENSION(:) :: dRINIXCO,dRINIYCO,dRINIZCO
      !  REAL, ALLOCATABLE, DIMENSION(:) :: dRFINXCO,dRFINYCO,dRFINZCO
      !  !!!!!!!!!!!!!!

       !kd-tree and SPH Related
       INTEGER :: FLAG_KD
       REAL, ALLOCATABLE :: HPART(:)
       REAL, ALLOCATABLE :: PART_DENS(:)
       type(KDTreeNode), pointer :: TREE
       REAL*4, ALLOCATABLE :: TREEPOINTS(:,:)
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       
       !Periodic boundary conditions
       INTEGER :: NPLUS
       INTEGER :: LOW1, LOW2
       INTEGER :: NPBC
       REAL :: SIZING
       REAL, ALLOCATABLE :: UBAS(:,:,:)
       REAL :: LPERIODIC(3)
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   


!*********************************************************************
       WRITE(*,*) 
       WRITE(*,*)
       WRITE(*,*) '           __      __  _____    _____   __  __    ' 
       WRITE(*,*) '     /\    \ \    / / |_   _|  / ____| |  \/  |   ' 
       WRITE(*,*) '    /  \    \ \  / /    | |   | (___   | \  / |   ' 
       WRITE(*,*) '   / /\ \    \ \/ /     | |    \___ \  | |\/| |   ' 
       WRITE(*,*) '  / ____ \    \  /     _| |_   ____) | | |  | |   ' 
       WRITE(*,*) ' /_/    \_\    \/     |_____| |_____/  |_|  |_|   ' 
       WRITE(*,*) '                                                  ' 
       WRITE(*,*) '                                                  ' 
       WRITE(*,*) 
       WRITE(*,*)
       CALL SLEEP(2)

       CALL IDATE(DATE)
       CALL ITIME(TIME)
       WRITE(*,*) 'DATE=',DATE(1),'/',DATE(2),'/',DATE(3)
       WRITE(*,*) 'TIME=',TIME(1),':',TIME(2),':',TIME(3) 

       !* TIME
       call system_clock(clock_t0,trate,tmax)

!*     NUMBER OF THREADS
      NUM=1
!$OMP PARALLEL SHARED(NUM)
!$OMP SINGLE
!$      NUM=OMP_GET_NUM_THREADS()
!$OMP END SINGLE NOWAIT
!$OMP END PARALLEL
      WRITE(*,*)'PROCESSOR =',NUM
      WRITE(*,*) 
      WRITE(*,*)


!*      OPENING FILES
       OPEN(1,FILE='config/voids.dat',STATUS='UNKNOWN',ACTION='READ')

!****************************************************
!*      READING INITIAL DATA
!****************************************************

       READ(1,*) 
       READ(1,*) !****************************************************
       READ(1,*) !*       General parameters block                      *
       READ(1,*) !****************************************************
       READ(1,*) 
       READ(1,*) FIRST,LAST,EVERY
       READ(1,*) 
       READ(1,*) NCOX,NCOY,NCOZ 
       READ(1,*) 
       READ(1,*) LEVMIN, LEVMAX 
       READ(1,*) 
       READ(1,*) ACHE,OMEGA0
       READ(1,*) 
       READ(1,*) LADO0
       READ(1,*) 
       READ(1,*) FLAG_PERIOD
       READ(1,*) !****************************************************
       READ(1,*) !*       Void finder input parameters                *  
       READ(1,*) !****************************************************
       READ(1,*) 
       READ(1,*) NVOID_MAX
       READ(1,*) 
       READ(1,*) DENS_THRE !used for the centers
       READ(1,*) 
       READ(1,*) DENS_THRE2 !used for the edges
       READ(1,*) 
       READ(1,*) GRAD_THRE
       READ(1,*) 
       READ(1,*) DIV_THRE
       READ(1,*) 
       READ(1,*) RMIN
       READ(1,*) 
       READ(1,*) RMIN_SUB
       READ(1,*) !****************************************************
       READ(1,*) !*       Type of data to process                      *
       READ(1,*) !****************************************************
       READ(1,*) 
       READ(1,*) FLAG_DATA
       READ(1,*) 
       READ(1,*) NHYX,NHYY,NHYZ
       READ(1,*) 
       READ(1,*) FILES_PER_SNAP,PARTTYPEX,MASSDM
       READ(1,*) !****************************************************
       READ(1,*) !*       Particle data handling parameters            *
       READ(1,*) !****************************************************
       READ(1,*)  
       READ(1,*) PARTIRED
       READ(1,*)  
       READ(1,*) FLAG_VEL_AVAILABLE
       READ(1,*)  
       READ(1,*) FLAG_VEL_INTERP
       READ(1,*)  
       READ(1,*) FLAG_DENS_INTERP
       READ(1,*)  
       READ(1,*) KNEIGHBOURS
       READ(1,*) !****************************************************
       READ(1,*) !*       Output management parameters                 *
       READ(1,*) !****************************************************
       READ(1,*)  
       READ(1,*) FLAG_WRITE_CUBES
      !  READ(1,*)  
      !  READ(1,*) FLAG_WRITE_PIECES
       READ(1,*) !****************************************************
       READ(1,*) !*        MASCLET ONLY (skip if another type of input is used)                *
       READ(1,*) !****************************************************
       READ(1,*)  
       READ(1,*) NPALEV, NLEVELS, NAMRX, NAMRY, NAMRZ
       READ(1,*)  
       READ(1,*) FLAG_DIV
       READ(1,*)  
       READ(1,*) FLAG_DENS
       CLOSE(1)

!****************************************************

       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! Check periodic boundary conditions consistency with k-d tree
#if periodic == 1
         IF (FLAG_PERIOD .EQ. 0) THEN
            WRITE(*,*) 'WARNING: PBCs not considered in voids.dat, but k-d tree has them activated!'
            WRITE(*,*) '         Compile with periodic=0!'
            STOP
         ENDIF
#else
         IF (FLAG_PERIOD .NE. 0) THEN
            WRITE(*,*) 'WARNING: PBCs considered in voids.dat, but k-d tree has them deactivated!'
            WRITE(*,*) '         Compile with periodic=1!'
            STOP
         ENDIF
#endif

       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! Allocate void finding variables
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ALLOCATE(INICIOX(NVOID_MAX), INICIOY(NVOID_MAX), INICIOZ(NVOID_MAX))
       ALLOCATE(PERIODICAL(NVOID_MAX))
       ALLOCATE(FINALX(NVOID_MAX), FINALY(NVOID_MAX), FINALZ(NVOID_MAX))
       ALLOCATE(ICX(NVOID_MAX),ICY(NVOID_MAX),ICZ(NVOID_MAX))
       ALLOCATE(VOL(NVOID_MAX))
       ALLOCATE(FATHER(NVOID_MAX))
       ALLOCATE(RINIXCO(NVOID_MAX),RINIYCO(NVOID_MAX),RINIZCO(NVOID_MAX))
       ALLOCATE(RFINXCO(NVOID_MAX),RFINYCO(NVOID_MAX),RFINZCO(NVOID_MAX))
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       !INPUT CHECKINGS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       !$!$ NEW, READ UNTIL LEVMAX IF FLAG_READ = 1
       IF (FLAG_READ .LT.0 .OR. FLAG_READ .GT. 1) STOP 'FLAG_READ must be 0 or 1'
       IF (FLAG_READ .EQ. 1) NL2 = LEVMAX
       IF (FLAG_READ .EQ. 0) NL2 = NLEVELS
       !$!$ FLAG_DATA 0,1,2,3
       IF (FLAG_DATA .LT.0 .OR. FLAG_DATA .GT. 3) STOP 'FLAG_DATA must be 0,1,2,3'
       !$!$ KNEIGHBOURS must be positive
       IF (KNEIGHBOURS .LT.0) STOP 'KNEIGHBOURS must be positive'
       !$!$ FLAG_VEL_AVAILABLE must be 0 or 1
       IF (FLAG_VEL_AVAILABLE .LT.0 .OR. FLAG_VEL_AVAILABLE .GT. 1) THEN
          STOP 'FLAG_VEL_AVAILABLE must be 0 or 1'
       ENDIF
       !for MASCLET and AREPO
       IF (FLAG_DATA .EQ.0 .OR. FLAG_DATA.EQ.3) THEN
         !$!$ FLAG_DIV must be 0, 1 or 2 
         IF (FLAG_DIV .LT.0 .OR. FLAG_DIV .GT. 2) STOP 'FLAG_DIV must be 0, 1 or 2'
         !$!$ FLAG_DENS must be 0, 1 or 2 
         IF (FLAG_DENS .LT.0 .OR. FLAG_DENS .GT. 2) STOP 'FLAG_DENS must be 0, 1 or 2'
       ENDIF
       !$!$ FLAG_VEL_INTERP must be 0 or 1
       IF (FLAG_VEL_INTERP .LT.0 .OR. FLAG_VEL_INTERP .GT. 1) STOP 'FLAG_VEL_INTERP must be 0 or 1'
       !$!$ FLAG_DENS_INTERP must be 0 or 1
       IF (FLAG_DENS_INTERP .LT.0 .OR. FLAG_DENS_INTERP .GT. 1) STOP 'FLAG_DENS_INTERP must be 0 or 1'
       !$!$ FLAG_WRITE_CUBES must be 0 or 1
       IF (FLAG_WRITE_CUBES .LT.0 .OR. FLAG_WRITE_CUBES .GT. 1) STOP 'FLAG_WRITE_CUBES must be 0 or 1'
      !  !$!$ FLAG_WRITE_PIECES must be 0 or 1
      !  IF (FLAG_WRITE_PIECES .LT.0 .OR. FLAG_WRITE_PIECES .GT. 1) STOP 'FLAG_WRITE_PIECES must be 0 or 1'
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       !$!$ Assume output written in dir "output_files"
       !$!$ Finder details will be written on the catalogue file
       DIR = 'output_files'

       !$!$ Create DIR if it does not exist
       call execute_command_line("mkdir -p " // trim(DIR), exitstat=exit_stat, cmdstat=cmd_stat)

       ! Check if the command worked (optional but recommended)
       if (exit_stat /= 0) then
         print *, "Error: Could not create output directory: " // trim(DIR)
         stop
       end if

       !$!$ initial prints
       IF (FLAG_DATA .EQ. 0) THEN
          WRITE(*,*) '-------> INPUT SIMULATION: MASCLET <-------'
       ELSE IF (FLAG_DATA .EQ. 1) THEN
          WRITE(*,*) '-------> INPUT PARTICLE/TRACER BINARY FILE <-------'
       ELSE IF (FLAG_DATA .EQ. 2) THEN
          WRITE(*,*) '-------> INPUT GRID BINARY FILE <-------'
       ELSE IF (FLAG_DATA .EQ. 3) THEN
          WRITE(*,*) '-------> INPUT SIMULATION: AREPO <-------'
          IF (PARTTYPEX .EQ. 1) THEN
             WRITE(*,*) 'Particle type: ** GAS **'
          ELSE
             WRITE(*,*) 'Particle type: ** DARK MATTER **'
             WRITE(*,*) 'Mass of DM particles:', MASSDM
          ENDIF
       ENDIF

       WRITE(*,*)

       !$!$ If GRID INPUT, allocate coarse grid
       IF (FLAG_DATA .EQ. 0 .OR. FLAG_DATA .EQ. 2) THEN
         NXX=NHYX
         NYY=NHYY
         NZZ=NHYZ
         LADO=LADO0-(LADO0/NXX)
         ALLOCATE(RADX(0:NXX+1))
         ALLOCATE(RADY(0:NYY+1))
         ALLOCATE(RADZ(0:NZZ+1))

         !*     GRID BUILDER  --> l=0 grid
         CALL MALLA(NXX,NYY,NZZ,LADO) !RADX, RADY, RADZ, DX, DY, DZ
         RADX0=RADX
         RADY0=RADY
         RADZ0=RADZ
         DX0=DX
         DY0=DY
         DZ0=DZ
       ENDIF
       !$!$!$!$!$!$!$!$!$!$!$!$!$!$!$!$!$!$

       WRITE(*,*)
       WRITE(*,*) '************ COSMOLOGY ******************** '    
       WRITE(*,*) 'h, omega_m', ACHE,OMEGA0
       WRITE(*,*)
       WRITE(*,*)
       WRITE(*,*) '************   INPUT  *********************' 
       WRITE(*,*) 'Iterations:         ', FIRST,LAST,EVERY
       IF (FLAG_DATA .EQ. 0 .OR. FLAG_DATA .EQ. 2) WRITE(*,*) 'Base mesh:  ', NHYX,NHYY,NHYZ,LADO0
       WRITE(*,*) 'Coarser grid        ', NCOX,NCOY,NCOZ
       WRITE(*,*) 'Threshold values    ', DENS_THRE, DENS_THRE2, GRAD_THRE
       WRITE(*,*) 'LEVMIN, LEVMAX, NL2:', LEVMIN, LEVMAX, NL2
       WRITE(*,*) 'Output directory:   ', TRIM(ADJUSTL(DIR))


!***********************************************************************
!***********************************************************************
!***********************************************************************
!      COSMOLOGICAL BACKGROUND (MASCLET UNITS):
!
!----- Summary:
!
!   a(z=0) = 1/10.98
!   8*pi*G = 1
!   c = 1
!
!   With these values, we get the following units for mass, lenght and time:
!   
!   l.u = 10.98 Mpc
!   l.t = 35.8 Myr
!   l.m = 9.1717E18 Msun
!
!***********************************************************************

       !numbers
       PI=DACOS(-1.D0)
       UNTERCIO=1.D0/3.D0

       !z = 0 constants
       ACHE0=ACHE*3.66D-3 !H0
       RODO=OMEGA0*3.D0*ACHE0**2 !rho_B(z = 0)
       RE0=1.0/10.98 !a(z=0)
       T0=2.D0*UNTERCIO/ACHE0 !t(z=0)

       !critical density in g/cm^3
       ROCRIT=1.879E-29*ACHE**2

!***********************************************************************
!***********************************************************************
!***********************************************************************



!***********************************************************************
      
       NITER=INT((LAST-FIRST)/EVERY) + 1

       WRITE(*,*)
       WRITE(*,*) 
       WRITE(*,*) '*******************************************************************'
       WRITE(*,*) '*******************************************************************'
       WRITE(*,*) '*******************************************************************'
       WRITE(*,*) '********************     Starting void finder     *****************'
       WRITE(*,*) '*******************************************************************'
       WRITE(*,*) '*******************************************************************'
       WRITE(*,*) '********************     Number of iterations:', NITER, '********'
       WRITE(*,*) '*******************************************************************'
       WRITE(*,*) '*******************************************************************'
       WRITE(*,*)
       
      
!*     LOOP OVER ITERATIONS
!*////////////////////////////////////
!*////////////////////////////////////
       DO IFI=1,NITER
!*////////////////////////////////////
!*////////////////////////////////////

       ITER=FIRST+EVERY*(IFI-1)
       
       WRITE(*,*)
       WRITE(*,*)
       WRITE(*,*)
       WRITE(*,*)
       WRITE(*,*) '////////////////////////////////////////////////////////////////'
       WRITE(*,*) '////////////////////////////////////////////////////////////////'
       WRITE(*,*) '////// NEW ITERATION:', ITER
       WRITE(*,*) '////////////////////////////////////////////////////////////////'
       WRITE(*,*) '////////////////////////////////////////////////////////////////'
       WRITE(*,*)
       WRITE(*,*) 'Reading input...' 

!*-------------------------------------------------------------------------------*
!*     READING AND SETTING DATA 
!*-------------------------------------------------------------------------------*
       
       IF (FLAG_DATA .EQ. 0) THEN
         WRITE(*,*) 'MASCLET data...'
         call system_clock(t1,trate,tmax)
         CALL READ_MASCLET(ITER, NDXYZ0, NL, T, ZETA)
         call system_clock(t2,trate,tmax)

       ELSE IF (FLAG_DATA .EQ. 1) THEN
         WRITE(*,*) 'Binary particle data...'
         IF(FLAG_VEL_AVAILABLE .EQ. 1) THEN
            WRITE(*,*) 'Velocities available: YES'
            call system_clock(t1,trate,tmax)
            CALL READ_BINARY_PART_1(ITER, ZETA)
            call system_clock(t2,trate,tmax)
            WRITE(*,*) 'Mass range:', MINVAL(MASAP(1:NPARTT)*UM), MAXVAL(MASAP(1:NPARTT)*UM)
         ELSE
            WRITE(*,*) 'Velocities available: NO'
            call system_clock(t1,trate,tmax)
            CALL READ_BINARY_PART_2(ITER, ZETA)
            call system_clock(t2,trate,tmax)
            WRITE(*,*) 'Mass range:', MINVAL(MASAP(1:NPARTT))*UM, MAXVAL(MASAP(1:NPARTT)*UM)
         ENDIF

       ELSE IF (FLAG_DATA .EQ. 2) THEN
         WRITE(*,*) 'Binary grid data...'
         call system_clock(t1,trate,tmax)
         CALL READ_BINARY_GRID(ITER, ZETA)
         call system_clock(t2,trate,tmax)

#if use_hdf5 == 1
       ELSE IF (FLAG_DATA .EQ. 3) THEN
         WRITE(*,*) 'AREPO data...'
         call system_clock(t1,trate,tmax)
         CALL READ_AREPO_HDF5(ITER,FILES_PER_SNAP,PARTTYPEX,MASSDM,ACHE)
         call system_clock(t2,trate,tmax)
         WRITE(*,*) 'Mass range:', MINVAL(MASAP(1:NPARTT)*UM), MAXVAL(MASAP(1:NPARTT)*UM)
#endif

       END IF
       
       !CHECK PARTICLES INSIDE BOUNDING BOX !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       IF (FLAG_DATA .EQ. 1 .OR. FLAG_DATA .EQ. 3) THEN
      
         FLAG_STOP = 0
         
         !$OMP PARALLEL SHARED(NPARTT, RXPA, RYPA, RZPA, LADO0), &
         !$OMP PRIVATE(I)
         !$OMP DO REDUCTION(+:FLAG_STOP)
         DO I=1,NPARTT
            IF ( (RXPA(I) .LT. -LADO0/2 .OR. RXPA(I) .GT. LADO0/2) .OR. &
                 (RYPA(I) .LT. -LADO0/2 .OR. RYPA(I) .GT. LADO0/2) .OR. &
                 (RZPA(I) .LT. -LADO0/2 .OR. RZPA(I) .GT. LADO0/2) ) THEN
                  FLAG_STOP = 1
            ENDIF
         ENDDO
         !$OMP END DO
         !$OMP END PARALLEL

         IF (FLAG_STOP .GE. 1) WRITE(*,*)
         IF (FLAG_STOP .GE. 1) WRITE(*,*)
         IF (FLAG_STOP .GE. 1) WRITE(*,*) 'Particles outside the bounding box! Check L0!'
         IF (FLAG_STOP .GE. 1) WRITE(*,*) 'Particles must be inside [-L0/2, L0/2] in each direction.'
         IF (FLAG_STOP .GE. 1) WRITE(*,*) 'STOPPING...'
         IF (FLAG_STOP .GE. 1) WRITE(*,*)
         IF (FLAG_STOP .GE. 1) STOP

       ENDIF
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


       WRITE(*,*) '///////////// Time (sec) spent during reading:', float(t2-t1)/1.e3
       WRITE(*,*)

       IF (ZETA.LT.0.0) ZETA=0.0

       !COSMOLOGY !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       !BACKGROUND DENSITY, SCALE FACTOR AND HUBBLE PARAMETER AT ZETA
       ROTE=RODO*(1.0+ZETA)**3 ! background density at z
       RETE=RE0/(1.0+ZETA) ! scale factor at z
       HTE=OMEGA0*(((RE0/RETE)**3) - 1.0) + 1.0
       HTE=ACHE*SQRT(HTE)
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       !WARNINGS FOR SPH KNEIGHBOURS
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       IF (FLAG_DATA .NE. 2) THEN

         WRITE(*,*) 'Mean number of particles per cell (COARSEST GRID):', REAL(NPARTT)/REAL(NCOX**3)
         WRITE(*,*) 'Number of kneighbours for SPH interpolation:', KNEIGHBOURS
         WRITE(*,*)
         IF ( (REAL(NPARTT)/REAL(NCOX*NCOY*NCOZ)) / REAL(KNEIGHBOURS) .LT. 0.01) THEN
            WRITE(*,*) 'WARNING!!!!!'
            WRITE(*,*) 'SPH kernel too big for the number of particles and grid: too smooth!'
            CALL SLEEP(3)
            WRITE(*,*) 'Decrease the number of nearest neighbours for SPH interpolation!'
            CALL SLEEP(3)
            WRITE(*,*)
         ENDIF

         IF ( REAL(NPARTT)/REAL(NCOX*NCOY*NCOZ) .LT. 0.1 ) THEN
            WRITE(*,*) 'WARNING!!!!!'
            WRITE(*,*) 'Poor number of tracer/particles compared to the grid...'
            CALL SLEEP(3)
            WRITE(*,*) 'Decrease the coarse grid size or provide a better sampling!'
            CALL SLEEP(3)
            WRITE(*,*)
         ENDIF

       ENDIF
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       WRITE(*,*)

       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! Building k-d tree
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       FLAG_KD = 0
       !MASCLET k-d tree only if dark-matter sph is required
       IF (FLAG_DATA .EQ. 0) THEN
         !DARK-MATTER used
         IF (FLAG_DENS .LE. 1 .OR. FLAG_DIV .LE. 1) THEN 
         !SPH used
         IF (FLAG_VEL_INTERP .EQ. 1 .OR. FLAG_DENS_INTERP .EQ. 1) THEN
            FLAG_KD = 1
         ENDIF
         ENDIF
       !PARTICLE DATA ALWAYS REQUIRES SPH SMOOTHING
       ELSE IF (FLAG_DATA .EQ. 1 .OR. FLAG_DATA .EQ. 3) THEN
         FLAG_KD = 1
       ENDIF

       IF (FLAG_KD .EQ. 1) THEN
         
         WRITE(*,*) 'Building k-d tree'

         ALLOCATE(HPART(NPARTT),TREEPOINTS(NPARTT,3), &
                  PART_DENS(NPARTT))

         HPART(:) = 0.
         PART_DENS(:) = 0.
         TREEPOINTS(:,1) = RXPA(1:NPARTT)
         TREEPOINTS(:,2) = RYPA(1:NPARTT)
         TREEPOINTS(:,3) = RZPA(1:NPARTT)

#if periodic == 1
            LPERIODIC = [LADO0, LADO0, LADO0]
            call system_clock(t1,trate,tmax)
            TREE => build_kdtree(TREEPOINTS,LPERIODIC)
            call system_clock(t2,trate,tmax)

#else
            call system_clock(t1,trate,tmax)
            TREE => build_kdtree(TREEPOINTS)
            call system_clock(t2,trate,tmax)
#endif

         WRITE(*,*) '///////////// Time (sec) spent building k-d tree ',float(t2-t1)/1.e3
         WRITE(*,*)

       END IF
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! OUTPUT FILE HEADERS
       CALL NOMFILE3(DIR,ITER,FILEO,FILEO1,FILEO2,FILEO3)

       !Void ASCII catalogue
       OPEN(UNIT=10, FILE=FILEO)
        
       !HEADER ON voids FILE
       NLEV = LEVMAX - LEVMIN + 1
       WRITE(10,*) NLEV, LEVMIN, LEVMAX, NCOX, NCOY, NCOZ, LADO0


       !Cubes ASCII catalogue
       IF (FLAG_WRITE_CUBES .EQ. 1) THEN
          OPEN(UNIT=12, FILE=FILEO2) 
          !HEADER ON cubes FILE: NLEV, LEVMIN, LEVMAX
          WRITE(12,*) NLEV, LEVMIN, LEVMAX
       END IF

      !  !Void pieces ASCII catalogue
      !  IF (FLAG_WRITE_PIECES .EQ. 1) THEN
      !     OPEN(UNIT=13, FILE=FILEO3) 
      !     !HEADER ON void pieces FILE: NLEV, LEVMIN, LEVMAX
      !     IF (FLAG_WRITE_PIECES .EQ. 1) THEN
      !        WRITE(13,*) NLEV, LEVMIN, LEVMAX
      !     END IF
      !  END IF
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!*-------------------------------------------------------------------------------*


       
!*     LOOP OVER DIFFERENT LEVELS 
!########################################
!########################################
       DO IR=LEVMIN, LEVMAX
!########################################
!########################################


          !CREATING GRID FOR THIS LEVEL
          !MASCLET or GRID CASE: levels come from base grid (level ordering matters)
          IF (FLAG_DATA .EQ. 0 .OR. FLAG_DATA .EQ. 2) THEN
           NXX=INT(REAL(NHYX)*(2.**(IR)))
           NYY=INT(REAL(NHYY)*(2.**(IR)))
           NZZ=INT(REAL(NHYZ)*(2.**(IR)))

          !PARTICLE CASE: level ordering does not matter, done with respect to levmin
          ELSE 
           NXX=INT(REAL(NCOX)*(2.**(IR-LEVMIN)))
           NYY=INT(REAL(NCOY)*(2.**(IR-LEVMIN)))
           NZZ=INT(REAL(NCOZ)*(2.**(IR-LEVMIN)))
          END IF   
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          WRITE(*,*)
          WRITE(*,*)
          WRITE(*,*)
          WRITE(*,*)
          WRITE(*,*)
          WRITE(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
          WRITE(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
          WRITE(*,*) '!!! NEW LEVEL:', IR                            
          WRITE(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
          WRITE(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
          WRITE(*,*)
          WRITE(*,*) 'Level:', IR, NXX, NYY, NZZ
          WRITE(*,*)


          IF(ALLOCATED(RADX) .EQV. .TRUE.) DEALLOCATE(RADX, RADY, RADZ)

          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !PERIODIC BOUNDARY CONDITIONS setup
          IF(FLAG_PERIOD .EQ. 0) THEN
            NPLUS = 0
            LOW1 = 1
            LOW2 = NXX
            LADO0PLUS = LADO0
            WRITE(*,*) 'PBCs are not applied...'
            WRITE(*,*) '    Grid', NXX, NYY, NZZ
            WRITE(*,*)
            ALLOCATE(RADX(0:NXX+1))
            ALLOCATE(RADY(0:NYY+1))
            ALLOCATE(RADZ(0:NZZ+1))
          ELSE
            SIZING = 50. !expected maximum size of voids outside the box
            DX = LADO0/NXX 
            !upper limit for PBC buffer
            SIZING = MIN(SIZING, LADO0/2.) 
            !lower limit for PBC buffer
            NPLUS = MAX( NXX/8, INT(SIZING/DX) ) 
            LOW1 = -NPLUS+1
            LOW2 = NXX+NPLUS
            !extended bounding box size
            LADO0PLUS = LADO0 + 2.*NPLUS*DX 
            ALLOCATE(RADX(LOW1:LOW2))
            ALLOCATE(RADY(LOW1:LOW2))
            ALLOCATE(RADZ(LOW1:LOW2))
            WRITE(*,*) 'Applying PBCs (extending grid)...'
            WRITE(*,*) '    Inner grid:', NXX, NYY, NZZ
            WRITE(*,*) '    New range:', LOW1, LOW2
            WRITE(*,*) '    Expanded L0:', LADO0, ' -> ', LADO0PLUS
            WRITE(*,*)
          ENDIF

          !number of cells in the extended grid
          NPBC = NXX + 2*NPLUS 

          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          LADO = LADO0-(LADO0/NXX)
          CALL MALLA(NXX,NYY,NZZ,LADO) ! --> RADX(0:NXX+1), DX

          !auxiliar array for inner grid calculation
          ALLOCATE(UBAS(NXX,NYY,NZZ))

          DXX=DX
          DYY=DY
          DZZ=DZ
          RX1=RADX(1)
          RY1=RADY(1)
          RZ1=RADZ(1)

          !cell volume
          VCELL=(LADO0/NXX)**3
          
          WRITE(*,*) 'GRID details:'
          WRITE(*,*) '     SIDE LENGTH, VCELL=',LADO0,VCELL
          WRITE(*,*) '     NX,DX,RADX(1),RADX(NX)=',NXX,DXX,RADX(1),RADX(NXX)
          WRITE(*,*)
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !CENTRES OF OUTER GRID (PBC)
          IF(FLAG_PERIOD .EQ. 1) THEN

            DO IX=-NPLUS+1,NXX+NPLUS
              IF (IX .LT. 1) THEN
                RADX(IX) = RADX(1) - (1-IX)*DXX
              ELSE IF (IX .GT. NXX) THEN
                RADX(IX) = RADX(NXX) + (IX-NXX)*DXX
              ENDIF
            ENDDO

            DO JY=-NPLUS+1,NYY+NPLUS
              IF (JY .LT. 1) THEN
                RADY(JY) = RADY(1) - (1-JY)*DYY
              ELSE IF (JY .GT. NYY) THEN
                RADY(JY) = RADY(NYY) + (JY-NYY)*DYY
              ENDIF
            ENDDO

            DO KZ=-NPLUS+1,NZZ+NPLUS
              IF (KZ .LT. 1) THEN
                RADZ(KZ) = RADZ(1) - (1-KZ)*DZZ
              ELSE IF (KZ .GT. NZZ) THEN
                RADZ(KZ) = RADZ(NZZ) + (KZ-NZZ)*DZZ
              ENDIF
            ENDDO

          ENDIF
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         !COMMON VARIABLES USED FOR VOID FINDING, UNIFORM GRID
          
          !DARK MATTER
          ALLOCATE(U1DMCO(NXX,NYY,NZZ)) 
          ALLOCATE(U2DMCO(NXX,NYY,NZZ))
          ALLOCATE(U3DMCO(NXX,NYY,NZZ))
          ALLOCATE(U4DMCO(NXX,NYY,NZZ))
          ALLOCATE(DIVERDMCO(NXX,NYY,NZZ))

          !GAS
          ALLOCATE(U1GCO(NXX,NYY,NZZ))
          ALLOCATE(U2GCO(NXX,NYY,NZZ))
          ALLOCATE(U3GCO(NXX,NYY,NZZ))
          ALLOCATE(U4GCO(NXX,NYY,NZZ))
          ALLOCATE(DIVERGCO(NXX,NYY,NZZ))

          !STARS (not used)
          ALLOCATE(U1SCO(NXX,NYY,NZZ))

          !$OMP PARALLEL DO SHARED(U1DMCO,U2DMCO,U3DMCO,U4DMCO,DIVERDMCO, &
          !$OMP                   U1GCO,U2GCO,U3GCO,U4GCO,DIVERGCO,U1SCO, &
          !$OMP                   NXX,NYY,NZZ), &
          !$OMP PRIVATE(IX,JY,KZ), DEFAULT(NONE)
          DO KZ=1,NZZ
          DO JY=1,NYY
          DO IX=1,NXX
            U1DMCO(IX,JY,KZ)=0.0 !DENSITY, NOT DENS.CONTRAST
            U2DMCO(IX,JY,KZ)=0.0
            U3DMCO(IX,JY,KZ)=0.0
            U4DMCO(IX,JY,KZ)=0.0
            DIVERDMCO(IX,JY,KZ)=0.0

            U1GCO(IX,JY,KZ)=0.0
            U2GCO(IX,JY,KZ)=0.0
            U3GCO(IX,JY,KZ)=0.0
            U4GCO(IX,JY,KZ)=0.0
            DIVERGCO(IX,JY,KZ)=0.0

            U1SCO(IX,JY,KZ)=0.0
          ENDDO
          ENDDO
          ENDDO
          
          !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          !voids-in-voids variables
          !lev=LEVMIN:  MARCAP = -1 in every cell, 
          !lev>LEVMIN:  MARCAP = ID(lev-1), tells to which lev-1 void does this cell belong

          IF(IR .EQ. LEVMIN) THEN
            ALLOCATE(MARCAP(NXX,NYY,NZZ))
            ALLOCATE(REQP(NVOID_MAX))
            MARCAP(:,:,:)=-1        

            !Father radius
            REQP(:)=0. 
          ENDIF

          !which 'l-1' void is the father of this 'l' void
          !it is updated inside VOIDFIND subroutine, using MARCAP
          FATHER(:)=0
          !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          

!*-------------------------------------------------------------------------------*
!*        FROM INPUT DATA TO UNIFORM GRID
!*-------------------------------------------------------------------------------*

          !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          !SPECIAL CASE: MASCLET
          !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          IF (FLAG_DATA .EQ. 0) THEN

            !---------------------------------------
            !  LEVEL < 0 ---------------------
            !---------------------------------------

            IF(IR .LT. 0) THEN
               !MASCLET, smooth l=0 grid 
                  WRITE(*,*) 'Starting smoothing...'
                  CALL SMOOTH(NXX,NYY,NZZ) 
                  WRITE(*,*) 'End of smoothing'
            ENDIF

            !---------------------------------------
            !  LEVEL = 0 -------------------
            !---------------------------------------

            IF(IR .EQ. 0) THEN
               U1DMCO=U1DM
               U1GCO=U1G
               U1SCO=U1S
               U2GCO=U2G
               U3GCO=U3G
               U4GCO=U4G
            ENDIF

            !INTERPOLATE PARTICLE VELOCITIES TO UNIFORM GRID
            IF(FLAG_DIV .LE. 1) THEN

               WRITE(*,*)
               WRITE(*,*) 'From lagrangian to eulerian...'
               IF (FLAG_VEL_INTERP .EQ. 0) WRITE(*,*) '     TSC KERNEL!'
               IF (FLAG_VEL_INTERP .EQ. 1) WRITE(*,*) '     SPH KERNEL!'

               !obtain velocities in the uniform grid comming from the particles
               call system_clock(t1,trate,tmax)
               IF (FLAG_VEL_INTERP .EQ. 0) CALL VEL_INTERP_TSC(NXX,NYY,NZZ,NPARTT, &
                                                RXPA,RYPA,RZPA,U2PA,U3PA,U4PA,U2DMCO,U3DMCO,U4DMCO)
               IF (FLAG_VEL_INTERP .EQ. 1) THEN

                  CALL PPART_DENS(NPARTT,TREE,MASAP,PART_DENS)

                  CALL VVEL_INTERP_SPH_VW(NXX,NYY,NZZ,NPARTT,TREE, &
                                 HPART,PART_DENS,MASAP,U2PA,U3PA,U4PA, &
                                 U2DMCO,U3DMCO,U4DMCO)

               ENDIF
               call system_clock(t2,trate,tmax)

               WRITE(*,*) '...done'
               WRITE(*,*) '///////////// Time (sec) spent during particle vel. interpolation (DM):', float(t2-t1)/1.e3
               WRITE(*,*)

            ENDIF

            WRITE(*,*)
            WRITE(*,*) 'Starting divergence calculation...'

            !---------------------------------------
            !  LEVEL <= 0 -------------------
            !---------------------------------------

            IF (IR .LE. 0) THEN
               !dm+gas 
               IF(FLAG_DIV .EQ. 0) THEN 
                  CALL DIVER_UNIFORM(NXX,NYY,NZZ,DXX,DYY,DZZ,U2DMCO,U3DMCO,U4DMCO,UBAS)
                  DIVERDMCO(1:NXX,1:NYY,1:NZZ) = UBAS
                  CALL DIVER_UNIFORM(NXX,NYY,NZZ,DXX,DYY,DZZ,U2GCO,U3GCO,U4GCO,UBAS)
                  DIVERGCO(1:NXX,1:NYY,1:NZZ) = UBAS
               !dm
               ELSE IF(FLAG_DIV .EQ. 1) THEN
                  CALL DIVER_UNIFORM(NXX,NYY,NZZ,DXX,DYY,DZZ,U2DMCO,U3DMCO,U4DMCO,UBAS)
                  DIVERDMCO(1:NXX,1:NYY,1:NZZ) = UBAS
               !gas
               ELSE IF(FLAG_DIV .EQ. 2) THEN
                  CALL DIVER_UNIFORM(NXX,NYY,NZZ,DXX,DYY,DZZ,U2GCO,U3GCO,U4GCO,UBAS)
                  DIVERGCO(1:NXX,1:NYY,1:NZZ) = UBAS
               ENDIF
            ENDIF

            !---------------------------------------
            !  LEVEL > 0 -------------------
            !---------------------------------------

            IF (IR .GE. 1) THEN

               !Computes Divergence using AMR hierarchy
               CALL DIVER_FINA_GAS(IR) !--> DIVR

               !Density and gas divergence calculation performed with AMR structure
               CALL VMESH(IR, NXX, NYY, NZZ, DXX, DYY, DZZ, RX1, RY1, RZ1, &
                     NHYX, NHYY, NHYZ) !--> UR, UGR, USR, DIVR, FLAGAMR (allocated within the subr)

               DO KZ=1, NZZ
                  DO JY=1,NYY
                     DO IX=1, NXX
                        U1DMCO(IX,JY,KZ)=UDMR(IX,JY,KZ)
                        U1GCO(IX,JY,KZ)=UGR(IX,JY,KZ)
                        U1SCO(IX,JY,KZ)=USR(IX,JY,KZ)
                        DIVERGCO(IX,JY,KZ)=DIVR(IX,JY,KZ)
                     ENDDO
                  ENDDO
               ENDDO

               DEALLOCATE(UDMR, UGR, USR, DIVR)

               !Particle divergence calculation performed in uniform grid
               IF (FLAG_DIV .EQ. 1 .OR. FLAG_DIV .EQ. 0) THEN
                  CALL DIVER_UNIFORM(NXX,NYY,NZZ,DXX,DYY,DZZ,U2DMCO,U3DMCO,U4DMCO,UBAS)
                  DIVERDMCO(1:NXX,1:NYY,1:NZZ) = UBAS
               ENDIF
               
            ENDIF

          !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          ! BINARY PART FILE or AREPO: hierarchy does not matter
          !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

          ELSE IF (FLAG_DATA .EQ. 1 .OR. FLAG_DATA .EQ. 3) THEN

          !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          IF(FLAG_VEL_AVAILABLE .EQ. 1) THEN
          !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            !Treated always internally as DM
            FLAG_DIV = 1 
            FLAG_DENS = 1

            WRITE(*,*)
            WRITE(*,*) 'From particles to grid ...'
            
            !---------------------------------------
            ! Particles velocities to grid
            !---------------------------------------

            WRITE(*,*)
            IF (FLAG_VEL_INTERP .EQ. 0) WRITE(*,*) 'Velocity TSC KERNEL!'
            IF (FLAG_VEL_INTERP .EQ. 1) WRITE(*,*) 'Velocity SPH KERNEL!'

            call system_clock(t1,trate,tmax)
            IF (FLAG_VEL_INTERP .EQ. 0) CALL VEL_INTERP_TSC(NXX,NYY,NZZ,NPARTT, &
                                       RXPA,RYPA,RZPA,U2PA,U3PA,U4PA,U2DMCO,U3DMCO,U4DMCO)
            IF (FLAG_VEL_INTERP .EQ. 1) THEN

               CALL PPART_DENS(NPARTT,TREE,MASAP,PART_DENS)

               CALL VVEL_INTERP_SPH_VW(NXX,NYY,NZZ,NPARTT,TREE, &
                              HPART,PART_DENS,MASAP,U2PA,U3PA,U4PA, &
                              U2DMCO,U3DMCO,U4DMCO)

            ENDIF              
            call system_clock(t2,trate,tmax)

            WRITE(*,*) '...done'
            WRITE(*,*) '///////////// Time (sec) spent during particle vel. interpolation (DM):', float(t2-t1)/1.e3   
            WRITE(*,*)
            !-----------------------------------------
            ! Particles to density
            !-----------------------------------------

            WRITE(*,*)
            IF (FLAG_DENS_INTERP .EQ. 0) WRITE(*,*) 'Density TSC KERNEL!'
            IF (FLAG_DENS_INTERP .EQ. 1) WRITE(*,*) 'Density SPH KERNEL!'

            call system_clock(t1,trate,tmax)
            IF(FLAG_DENS_INTERP.EQ.0) CALL DENS_INTERP_TSC(NXX,NYY,NZZ, &
                                       NPARTT,RXPA,RYPA,RZPA,MASAP,U1DMCO) !--> U1DMCO

            IF(FLAG_DENS_INTERP.EQ.1) THEN

               !IF SPH VELOCITY RECONSTRUCTION IS NOT CALLED PRIOR TO DDENS
               !WE NEED TO CALL H_DISTANCE TO COMPUTE HPART
               !(Although this situation should be avoided)

               IF(FLAG_VEL_INTERP .EQ. 0) THEN
                  WRITE(*,*) 'Calculating smoothing length...'
                  CALL H_DISTANCE(NXX,NYY,NZZ,NPARTT,TREE,HPART)
               ENDIF

               CALL DDENS_INTERP_SPH(NXX,NYY,NZZ,NPARTT, &
                                          RXPA,RYPA,RZPA,HPART,MASAP,U1DMCO) !--> U1DMCO
                                    
            ENDIF
            call system_clock(t2,trate,tmax)

            WRITE(*,*) '...done'
            WRITE(*,*) '///////////// Time (sec) spent during particle dens. interpolation (DM):', float(t2-t1)/1.e3
            WRITE(*,*)
            
            !U1DMCO is mass until this point, we need density
            U1DMCO = U1DMCO / (DXX*DYY*DZZ*RETE**3) !density in MASCLET units
            U1DMCO = U1DMCO / ROTE !overdensity with respect to background density at z
         
            !---------------------------------------
            ! Divergence calculation
            !---------------------------------------

            CALL DIVER_UNIFORM(NXX,NYY,NZZ,DXX,DYY,DZZ,U2DMCO,U3DMCO,U4DMCO,UBAS)
            DIVERDMCO(1:NXX,1:NYY,1:NZZ) = UBAS

            !---------------------------------------
            ! Deallocate PARTICLE input variables
            !---------------------------------------
            DEALLOCATE(U2PA, U3PA, U4PA, RXPA, RYPA, RZPA, MASAP)
            !---------------------------------------

          !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          ELSE !VELOCITIES NOT AVAILABLE
          !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            FLAG_DIV = -1 !minus indicates divergence will be
            !              reconstructed from density field 
            FLAG_DENS = 1

            WRITE(*,*)
            WRITE(*,*) 'From lagrangian to eulerian...'

            !-----------------------------------------
            ! Particles to density
            !-----------------------------------------

            WRITE(*,*)
            IF (FLAG_DENS_INTERP .EQ. 0) WRITE(*,*) 'Density TSC KERNEL!'
            call system_clock(t1,trate,tmax)
            IF(FLAG_DENS_INTERP.EQ.0) CALL DENS_INTERP_TSC(NXX,NYY,NZZ, &
                                       NPARTT,RXPA,RYPA,RZPA,MASAP,U1DMCO) !--> U1DMCO

            IF(FLAG_DENS_INTERP.EQ.1) THEN

               !SINCE THIS TIME VVEL IS NOT CALLED PRIOR TO DDENS
               !WE NEED TO CALL H_DISTANCE TO COMPUTE HPART
               WRITE(*,*) 'Calculating smoothing length...'
               CALL H_DISTANCE(NXX,NYY,NZZ,NPARTT,TREE,HPART)

               WRITE(*,*) 'Applying density SPH KERNEL!'

               CALL DDENS_INTERP_SPH(NXX,NYY,NZZ,NPARTT, &
                                     RXPA,RYPA,RZPA,HPART,MASAP,U1DMCO) !--> U1DMCO
                                    
            ENDIF
            call system_clock(t2,trate,tmax)

            WRITE(*,*) '...done'
            WRITE(*,*) '///////////// Time (sec) spent during particle dens. interpolation (DM):', float(t2-t1)/1.e3
            WRITE(*,*)
            !U1DMCO is mass until this point, we need density
            U1DMCO = U1DMCO / (DXX*DYY*DZZ*RETE**3) !density in MASCLET units
            U1DMCO = U1DMCO / ROTE !overdensity with respect to background density at z
            
          !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          END IF !if velocities are available or not
          !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

          !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          ! BINARY GRID FILE: hierarchy MATTERS
          !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           
          ELSE IF (FLAG_DATA .GE. 2) THEN

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !GRID DATA is lodaded onto U1G, U2G, U3G, U4G 
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            !Treated always internally as gas
            FLAG_DIV = 2
            FLAG_DENS = 2

            !---------------------------------------
            !  LEVEL < 0 ---------------------
            !---------------------------------------

            IF(IR .LT. 0) THEN
               !smooth l=0 grid 
                  WRITE(*,*) 
                  WRITE(*,*) 'Starting smoothing...'
                  CALL SMOOTH_2(NXX,NYY,NZZ) 
                  WRITE(*,*) 'End of smoothing'
                  WRITE(*,*)
            ENDIF

            !---------------------------------------
            !  LEVEL = 0 -------------------
            !---------------------------------------

            IF(IR .EQ. 0) THEN
               U1GCO=U1G
               U2GCO=U2G
               U3GCO=U3G
               U4GCO=U4G
            ENDIF

            !---------------------------------------
            ! LEVEL > 0 ----------------------------
            !---------------------------------------

            IF (IR .GT. 0) STOP 'LEVEL > 0 FOR BIN_GRID NOT IMPLEMENTED YET'

            !---------------------------------------
            ! Divergence calculation
            !---------------------------------------
            ! WRITE(*,*) MINVAL(U2GCO), MAXVAL(U2GCO)
            ! WRITE(*,*) MINVAL(U3GCO), MAXVAL(U3GCO)
            ! WRITE(*,*) MINVAL(U4GCO), MAXVAL(U4GCO)
            CALL DIVER_UNIFORM(NXX,NYY,NZZ,DXX,DYY,DZZ,U2GCO,U3GCO,U4GCO,UBAS)
            DIVERGCO(1:NXX,1:NYY,1:NZZ) = UBAS
            
          ENDIF
          !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

         !* Which density components to consider
          ALLOCATE(U1CO(LOW1:LOW2,LOW1:LOW2,LOW1:LOW2))
          IF(FLAG_DENS .EQ. 0) U1CO(1:NXX,1:NYY,1:NZZ) = U1DMCO + U1GCO
          IF(FLAG_DENS .EQ. 1) U1CO(1:NXX,1:NYY,1:NZZ) = U1DMCO
          IF(FLAG_DENS .EQ. 2) U1CO(1:NXX,1:NYY,1:NZZ) = U1GCO

         !* Until here, density is in rho_background units -> to units of mean density
          IF(FLAG_DATA .NE. 2) THEN
            U1CO = U1CO*ROTE
            MEANDENS = SUM(U1CO(1:NXX,1:NYY,1:NZZ)) / REAL(NXX*NYY*NZZ)
            U1CO = U1CO / MEANDENS 
            MEANDENS = MEANDENS * (UM / UL**3)
            WRITE(*,*) 'Mean density, background, fraction (Msun/Mpc^3):', MEANDENS, ROTE * (UM / UL**3), &
                        MEANDENS / (RODO * (UM / UL**3))

          ELSE !GRID input overdensity: assumed in rho_background units
            MEANDENS = ROTE * (UM / UL**3) !mean density in Msun/Mpc^3
          ENDIF
 
         !* Which divergence components to consider
          ALLOCATE(DIVERCO(LOW1:LOW2,LOW1:LOW2,LOW1:LOW2))
          IF (FLAG_DIV .EQ. 2) DIVERCO(1:NXX,1:NYY,1:NZZ) = DIVERGCO
          IF (FLAG_DIV .EQ. 1) DIVERCO(1:NXX,1:NYY,1:NZZ) = DIVERDMCO
          IF (FLAG_DIV .EQ. 0) THEN
            DO KZ=1,NZZ
               DO JY=1,NYY
                  DO IX=1,NXX
                     DIVERCO(IX,JY,KZ)=MIN(DIVERGCO(IX,JY,KZ), DIVERDMCO(IX,JY,KZ))
                  ENDDO
               ENDDO
            ENDDO
          ENDIF

          !SPECIAL CASE-> FLAG_DIV = -1, velocity reconstruction from density field
          IF (FLAG_DIV .EQ. -1) THEN
            DIVERCO(:,:,:) = - OMEGA0**(0.6)*RETE*HTE*(U1CO(:,:,:) - 1.0)
            WRITE(*,*) 'Velocity reconstructed from density field with linear theory'
          ENDIF

          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !PERIODIC BOUNDARY CONDITIONS
          !Fill ghost cells with periodic boundary conditions
            IF(FLAG_PERIOD .EQ. 1) THEN
               CALL REAL_MAKE_PERIODIC(NXX,NYY,NZZ,NPLUS,U1CO)
               CALL REAL_MAKE_PERIODIC(NXX,NYY,NZZ,NPLUS,DIVERCO)
            ENDIF
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!*-------------------------------------------------------------------------------*
!*-------------------------------------------------------------------------------*


!*-------------------------------------------------------------------------------*
!*      Marking cells suitable for being void centres
!*-------------------------------------------------------------------------------*

          !PBC variables
          !MARCAP2 extends MARCAP to the PBC grid
          IF (ALLOCATED(MARCAP2)) DEALLOCATE(MARCAP2)
          ALLOCATE(MARCAP2(LOW1:LOW2,LOW1:LOW2,LOW1:LOW2))
          MARCAP2(1:NXX,1:NYY,1:NZZ) = MARCAP
          CALL INT_MAKE_PERIODIC(NXX,NYY,NZZ,NPLUS,MARCAP2)

          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !--> FLAGV mark cells suitable for being void centres
          !--> FLAG_SUB tells which cells are allowed to grow a void 
          !             (all if l=levmin and only those with a void at lev-1 otherwise)
          CALL MARK_ALL(LOW1,LOW2)
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          CONTA=0
          DO KZ=LOW1,LOW2
            DO JY=LOW1,LOW2
               DO IX=LOW1,LOW2
                  IF(FLAGV(IX,JY,KZ).EQ.1) CONTA=CONTA+1
               ENDDO
            ENDDO
          ENDDO

          WRITE(*,*)
          WRITE(*,*) 'Volume fraction allowed to belong to a void: ', COUNT(FLAG_SUB .GT. 0)/(NPBC**3.)
          WRITE(*,*) 'Potential void centres, fraction ', CONTA, CONTA/(NPBC**3.)
          WRITE(*,*) 
          WRITE(*,*)


!*-------------------------------------------*
!*      VOIDFIND -> main subroutine
!*-------------------------------------------*

          WRITE(*,*)
          WRITE(*,*) '************************************************'
          WRITE(*,*) 'Finding voids...'
          call system_clock(t1,trate,tmax)
         
          CALL VOIDFIND(IR,LOW1,LOW2,DXX,DYY,DZZ,RX1,RY1,RZ1, &
                        NVOID_MAX,REQP,NVOID)

          call system_clock(t2,trate,tmax)

          !$!$ AFTER VOIDFIND:
          !common variables INICIOXYZ,FINALXYZ,RINIXYZCO,RFINXYZCO,ICXY
          !are updated with the cubes at this level

          DEALLOCATE(FLAGV, FLAG_SUB)

          WRITE(*,*) 'Found ',NVOID,' cubes'
          WRITE(*,*) '///////////// Time (sec) spent during voidfind:', float(t2-t1)/1.e3
          WRITE(*,*)
          
          !SORTING CUBES BY VOLUME

          ALLOCATE(INDICE(NVOID))
          ALLOCATE(INDICE2(NVOID))
          ALLOCATE(VOL2(NVOID))

          DO I=1, NVOID
             INDICE(I)=0
             INDICE2(I)=0
             VOL2(I)=0.D0
          ENDDO

          !$!$ lowest to highest
          CALL INDEXX(NVOID,VOL(1:NVOID),INDICE2)

          IV=NVOID
          DO I=1, NVOID
             VOL2(IV)=VOL(INDICE2(I))
             INDICE(IV)=INDICE2(I) !highest to lowest
             IV=IV-1
          ENDDO

          DEALLOCATE(INDICE2, VOL2)

          WRITE(*,*) 'Cubes sorted by volume:'
          WRITE(*,*) 'Largest (vol/Re):', VOL(INDICE(1)), &
                                          ((3.*VOL(INDICE(1)))/(4.*PI))**0.333
          WRITE(*,*) 'Smallest (vol/Re):', VOL(INDICE(NVOID)), &
                                          ((3.*VOL(INDICE(NVOID)))/(4.*PI))**0.333
          WRITE(*,*) '************************************************'
          WRITE(*,*)

          !allocate void properties
          ALLOCATE(XC(NVOID),YC(NVOID),ZC(NVOID)) !centre of void, maximum divergence
          ALLOCATE(GXC(NVOID),GYC(NVOID),GZC(NVOID)) !volume-weighted geometrical centre
          ALLOCATE(UVOID(NVOID))
          ALLOCATE(VOLNEW(NVOID))
          ALLOCATE(UMEAN(NVOID))
          ALLOCATE(MTOT(NVOID))
          ALLOCATE(NCELLV(NVOID))
          ALLOCATE(EPS(NVOID))
          ALLOCATE(IP(NVOID))
          ALLOCATE(REQ(NVOID))
          !to which void does a cell belong
          ALLOCATE(MARCA(LOW1:LOW2,LOW1:LOW2,LOW1:LOW2))

!*-------------------------------------------*
!*      From cubes to voids: MERGING
!*-------------------------------------------*


          WRITE(*,*) '************************************************'
          WRITE(*,*) 'Merging cubes ...'

          call system_clock(t1,trate,tmax)
          CALL MERGE_VOID(NVOID,INDICE,LOW1,LOW2,DXX,DYY,DZZ,VOLNEW,UVOID,GXC,GYC,GZC,NVOID2)
          call system_clock(t2,trate,tmax)

          WRITE(*,*) '///////////// Time (sec) spent during merging:', float(t2-t1)/1.e3
          WRITE(*,*)

          WRITE(*,*) 'Num. Cells in voids', COUNT(MARCA(1:NXX,1:NYY,1:NZZ) .GT. 0)
          WRITE(*,*) 'Number of voids after merging:', NVOID2
          WRITE(*,*) 'FF (Volume)', REAL(COUNT(MARCA(1:NXX,1:NYY,1:NZZ) .GT. 0))/(NXX*NYY*NZZ)
          WRITE(*,*) '************************************************'
          WRITE(*,*)

!*-------------------------------------------*
!*      POST-PROCESSING
!*-------------------------------------------*

          !Number of cells belonging to each void
          NCELLV(:)=0 
          DO KZ=LOW1,LOW2
            DO JY=LOW1,LOW2
               DO IX=LOW1,LOW2
                  IF(MARCA(IX,JY,KZ) .GT. 0.) THEN 
                     IND=MARCA(IX,JY,KZ)
                     IF(UVOID(IND) .EQ. -1 ) THEN 
                        NCELLV(IND)=NCELLV(IND)+1
                     ENDIF
                  ENDIF
               ENDDO
            ENDDO
          ENDDO
 

          !min num. of cells for keeping voids; for few cell voids inertia tensor not reliable
          DO I=1,NVOID
            IND = INDICE(I)
            IF(UVOID(IND) .GE. 0) CYCLE
            IF(NCELLV(IND) .LT. NCELL_MIN) THEN
               WHERE (MARCA  .EQ. IND)
                  MARCA = 0
               END WHERE
               UVOID(IND)=0
               WHERE (UVOID .EQ. IND)
                  UVOID = 0
               END WHERE
               VOLNEW(IND)=0.
               NCELLV(IND)=0
            ENDIF
          ENDDO

          !Compute 3D Ellipsoid fitting with GEOMETRICAL CENTRE
          !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          CALL SHAPE(LOW1,LOW2,NVOID,INDICE,NCELLV,UVOID,GXC,GYC,GZC,EPS,IP)
          !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

          !REDEFINE VOLNEW, WITH NCELLV
          VOLNEW(:)=0.
          DO IV=1, NVOID
            IND=INDICE(IV)
            IF(UVOID(IND) .GE. 0) CYCLE
            IND=INDICE(IV)
            VOLM=DBLE(NCELLV(IND))*DBLE(DXX*DYY*DZZ)
            VOLNEW(IND)=VOLM
          ENDDO  

          !min radius
          DO I=1,NVOID
            IND = INDICE(I)
            IF(UVOID(IND) .GE. 0) CYCLE
            IF(VOLNEW(IND).LT.4./3.*PI*RMIN**3) THEN
               WHERE (MARCA  .EQ. IND)
                  MARCA = 0
               END WHERE
               UVOID(IND)=0
               WHERE (UVOID .EQ. IND)
                  UVOID = 0
               END WHERE
               VOLNEW(IND)=0.
               NCELLV(IND)=0
            ENDIF
          ENDDO

          !Calculate centre for every void: maximum divergence
          !This corresponds to the geometrical centre of the cube to which
          !all the others are merged
          DO IV=1,NVOID
            IND=INDICE(IV)
            IF(UVOID(IND) .GE. 0) CYCLE
            XC(IND) = RADX(1)+(ICX(IND)-1)*DXX
            YC(IND) = RADY(1)+(ICY(IND)-1)*DYY
            ZC(IND) = RADZ(1)+(ICZ(IND)-1)*DZZ
          ENDDO
          

          !Voids are simply connected. Holes have to be filled.
          WRITE(*,*) '************************************************'
          WRITE(*,*) 'Filling holes within voids (Simple Connection)...'
          call system_clock(t1,trate,tmax)
         !  CALL SIMPLE_CONNECTION(NVOID,INDICE,UVOID,MARCA,LOW1,LOW2)
          CALL SIMPLE_CONNECTION_FAST(NVOID,MARCA,LOW1,LOW2)
          call system_clock(t2,trate,tmax)
          WRITE(*,*) '///////////// Time (sec) spent during hole filling:', float(t2-t1)/1.e3
          WRITE(*,*) '************************************************'
          WRITE(*,*)
          

          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !MORE POST-PROCESSING CAN BE DONE HERE

          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


          !PERIODIC BOUNDARY CONDITIONS: check which voids are the same
          !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
          IF(FLAG_PERIOD .EQ. 1) THEN
            WRITE(*,*) '************************************************'
            WRITE(*,*) 'Imposing periodic boundary conditions...'
               call system_clock(t1,trate,tmax)
               WRITE(*,*) NXX,NYY,NZZ
               CALL VOID_PERIODIC(NVOID,INDICE,UVOID,GXC,GYC,GZC, &
                                 VOLNEW,NXX,NYY,NZZ,DXX,DYY,DZZ, &
                                 LOW1,LOW2)
               call system_clock(t2,trate,tmax)
            WRITE(*,*) '///////////// Time (sec) spent during PBC:', float(t2-t1)/1.e3
            WRITE(*,*) '************************************************'
            WRITE(*,*)
            !AFTER MATCHING PERIODIC IMAGES AND MARCA AND UVOID UPDATED:
            ! - Only 1:NXX,1:NYY,1:NZZ cells are considered now
            ! - CENTRES AND SHAPE ARE NOT UPDATED, the inner voids remain
            !   unchanged in this matter
            ! - Recalculate now the number of cells 

          ENDIF !PBC
          !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

!*-------------------------------------------*
!*      GLOBAL (REAL) VOID PROPERTIES (NOW INSIDE 1:NXX,1:NYY,1:NZZ)
!*-------------------------------------------*

          !UPDATE NCELLV, UMEAN, MTOT, VOLNEW WITH FINAL MARCA
          NCELLV(:)=0 !number of cells in each void
          UMEAN(:)=0. !overdensity: dens/dens_B
          MTOT(:)=0. !total mass in each void
          VOLT_CLEAN=0 !total volume in real voids
          DO KZ=1,NZZ
            DO JY=1,NYY
               DO IX=1, NXX
                  IF(MARCA(IX,JY,KZ) .GT. 0.) THEN 
                     IND=MARCA(IX,JY,KZ)
                     IF(UVOID(IND) .EQ. -1 ) THEN 
                        !VOIDS
                        NCELLV(IND)=NCELLV(IND)+1
                        UMEAN(IND)=UMEAN(IND)+U1CO(IX,JY,KZ)
                        MTOT(IND)=MTOT(IND)+U1CO(IX,JY,KZ)*MEANDENS*DXX*DYY*DZZ*RETE**3 !in Msun
                        !TOTAL
                        VOLT_CLEAN=VOLT_CLEAN+DBLE(DXX*DYY*DZZ) !comoving
                     ENDIF
                  ENDIF
               ENDDO
            ENDDO
          ENDDO

          !mean density for every voids
          WHERE (NCELLV .GT. 0)
            UMEAN = UMEAN / REAL(NCELLV)
          END WHERE

          !REDEFINE VOLNEW, WITH NCELLV
          VOLNEW(:)=0.
          DO IV=1, NVOID
            IND=INDICE(IV)
            IF(UVOID(IND) .GE. 0) CYCLE
            IND=INDICE(IV)
            VOLM=DBLE(NCELLV(IND))*DBLE(DXX*DYY*DZZ)
            VOLNEW(IND)=VOLM
          ENDDO  

          !min radius
          DO I=1,NVOID
            IND = INDICE(I)
            IF(UVOID(IND) .GE. 0) CYCLE
            IF(VOLNEW(IND).LT.4./3.*PI*RMIN**3) THEN
               WHERE (MARCA  .EQ. IND)
                  MARCA = 0
               END WHERE
               UVOID(IND)=0
               WHERE (UVOID .EQ. IND)
                  UVOID = 0
               END WHERE
               VOLNEW(IND)=0.
               NCELLV(IND)=0
            ENDIF
          ENDDO

          !NUMBER OF VOIDS AFTER POST-PROCESSING
          NVOIDT=COUNT(UVOID .EQ. -1)

          WRITE(*,*)
          WRITE(*,*)
          WRITE(*,*) '\\\\\ AFTER POST-PROCESSING \\\\\'
          WRITE(*,*) '----------------------------------------------'
          WRITE(*,*) 'Number of voids (R>Rmin, Ncell>Nmin):', NVOIDT
         !  WRITE(*,*) 'Num. voids above R = 10 cMpc:', COUNT(((3.*VOLNEW)/(4.*PI))**(1./3.) .GT. 10.)
          WRITE(*,*) 'FF (Volume)', VOLT_CLEAN/LADO0**3
          WRITE(*,*) '----------------------------------------------'
          WRITE(*,*) 'MIN(RAD), MAX(RAD) [cMpc]:', MINVAL(((3.*VOLNEW)/(4.*PI))**(1./3.), MASK=(UVOID == -1)), &
                                            MAXVAL(((3.*VOLNEW)/(4.*PI))**(1./3.), MASK=(UVOID == -1))
          WRITE(*,*) 'MIN(DELTA), MAX(DELTA):', MINVAL(UMEAN-1., MASK=(UVOID == -1)), &
                                               MAXVAL(UMEAN-1., MASK=(UVOID == -1))
          WRITE(*,*) 'MIN(MASS), MAX(MASS) [Msun]:', MINVAL(MTOT, MASK=(UVOID == -1)), &
                                               MAXVAL(MTOT, MASK=(UVOID == -1))
          WRITE(*,*) '----------------------------------------------' 
          WRITE(*,*) '\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\'       
          WRITE(*,*)
          WRITE(*,*)    
          WRITE(*,*)  

         !  !Using MARCA, describe each VOID with a set of non-overlapping parallelepipeds
         !  !These parallelepipeds will be the "void pieces" with which we can search for
         !  !galaxies, etc. They describe the VOID structure.
         !  IF (FLAG_WRITE_PIECES .EQ. 1) THEN
         !    WRITE(*,*) '************************************************'
         !    WRITE(*,*) 'Deconstructing voids...'

         !    ALLOCATE(dUVOID(NVOID_MAX),dINICIOX(NVOID_MAX),dINICIOY(NVOID_MAX),dINICIOZ(NVOID_MAX), &
         !             dFINALX(NVOID_MAX),dFINALY(NVOID_MAX),dFINALZ(NVOID_MAX),dRINIXCO(NVOID_MAX), &
         !             dRINIYCO(NVOID_MAX),dRINIZCO(NVOID_MAX),dRFINXCO(NVOID_MAX),dRFINYCO(NVOID_MAX), &
         !             dRFINZCO(NVOID_MAX))

         !    call system_clock(t1,trate,tmax)
         !    CALL VOID_DECONSTRUCTION(NUM,NVOID,INDICE,UVOID,MARCA, &
         !                            NXX,NYY,NZZ,DXX,DYY,DZZ,RX1,RY1,RZ1, &
         !                            dNVOID,dUVOID, &
         !                            dINICIOX,dINICIOY,dINICIOZ, &
         !                            dFINALX,dFINALY,dFINALZ, &
         !                            dRINIXCO,dRINIYCO,dRINIZCO, &
         !                            dRFINXCO,dRFINYCO,dRFINZCO)
         !    call system_clock(t2,trate,tmax)
         !    WRITE(*,*) '///////////// Time (sec) spent during void deconstruction:', float(t2-t1)/1.e3
         !    WRITE(*,*) '************************************************'
         !    WRITE(*,*)
         !    WRITE(*,*)
         !  END IF

!*-------------------------------------------------------------*
!    OUTPUT
!*-------------------------------------------------------------*

          !Binary 3D maps
          OPEN(UNIT=11, FILE=FILEO1, FORM='UNFORMATTED') 

          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !3D arrays file
          WRITE(*,*) '************* OUTPUT ******************'
          WRITE(*,*) 'Writing 3D arrays:',NXX,NYY,NZZ
          WRITE(11) (((MARCA(IX,JY,KZ), IX=1,NXX), JY=1,NYY), KZ=1,NZZ)
          WRITE(11) (((-1. + U1CO(IX,JY,KZ), IX=1,NXX), JY=1,NYY), KZ=1,NZZ)
          WRITE(11) (((DIVERCO(IX,JY,KZ), IX=1,NXX), JY=1,NYY), KZ=1,NZZ)
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !!!! CATALOGUE file

          !^^^^^^^^^^^^^^HEADER^^^^^^^^^^^^^^^^
          WRITE(10,*) IR, NVOID, NVOIDT, NVOIDP, VOLT_CLEAN/(LADO0**3), MEANDENS


          !Optional: cubes file, written at the same time that catalogue file
          !^^^^^^^^^^^^^^HEADER(cubes)^^^^^^^^^^^^^^^^
          IF (FLAG_WRITE_CUBES .EQ. 1) THEN
             WRITE(12,*) IR, NVOID, NVOIDT
          END IF


          DO IV=1, NVOID
             IND=INDICE(IV)
         
             !cubes making up voids
             IF(FLAG_WRITE_CUBES .EQ. 1) THEN
                WRITE(12,*) IND, UVOID(IND), INICIOX(IND), FINALX(IND), INICIOY(IND), FINALY(IND), &
                                             INICIOZ(IND), FINALZ(IND), RINIXCO(IND), RFINXCO(IND), &
                                             RINIYCO(IND), RFINYCO(IND), RINIZCO(IND), RFINZCO(IND)
             END IF

             !catalogue: only main voids
             IF(UVOID(IND) .NE. -1) CYCLE

             VOLM = VOLNEW(IND)
             REQ(IND)=((3.*VOLM)/(4.*PI))**(1./3.) !cMpc

             !Find parent void:
             INDP=0
             REQPP=0
             IF(IR .GT. LEVMIN ) THEN
               INDP=FATHER(IND)
               REQPP=REQ(INDP)
             ENDIF

             WRITE(10,*) IND, XC(IND), YC(IND), ZC(IND), GXC(IND), GYC(IND), GZC(IND), &
                         VOLM, REQ(IND), UMEAN(IND), EPS(IND), &
                         IP(IND), INDP, REQPP, MTOT(IND)

          ENDDO
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


         !  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !  !!!! void pieces file
         !  !^^^^^^^^^^^^^^HEADER (Void pieces)^^^^^^^^^^^^^^^^
         !  !Optional: void pieces file
         !  IF (FLAG_WRITE_PIECES  .EQ. 1) THEN
         !    WRITE(13,*) IR, dNVOID, NVOIDT
         !    DO IV=1,dNVOID
         !       WRITE(13,*) IV, dUVOID(IV), dINICIOX(IV), dFINALX(IV), dINICIOY(IV), dFINALY(IV), &
         !                    dINICIOZ(IV), dFINALZ(IV), dRINIXCO(IV), dRFINXCO(IV), &
         !                    dRINIYCO(IV), dRFINYCO(IV), dRINIZCO(IV), dRFINZCO(IV)
         !    ENDDO
         !  END IF
         !  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


          WRITE(*,*)
          WRITE(*,*) 'Void list written in: ', FILEO
          WRITE(*,*) '3D arrays written in: ', FILEO1
          WRITE(*,*) 
          WRITE(*,*) 'Number of voids written:', NVOIDT
          WRITE(*,*) '************************************************'
          WRITE(*,*)
          WRITE(*,*)


!-------------------------------------------------------
! Store variables used in the next level (voids-in-voids)
!--------------------------------------------------------

      !~~~~~~~~~~~~~~~~~~~~~~~
      IF(IR .LT. LEVMAX) THEN
      !~~~~~~~~~~~~~~~~~~~~~~~
         
          DEALLOCATE(MARCAP, REQP)

          !Next grid size
          NX2=NXX*2
          NY2=NYY*2
          NZ2=NZZ*2

          WRITE(*,*) 'Storing variables for next level...'
          WRITE(*,*) 'Allocating new MARCAP', NX2, NY2, NZ2
          WRITE(*,*)

          ALLOCATE(MARCAP(NX2,NY2,NZ2))
          ALLOCATE(REQP(NVOID))

          DO IV=1, NVOID
             IND=INDICE(IV)
             IF(UVOID(IND) .NE. -1) CYCLE
             REQP(IND)=REQ(IND)
          ENDDO

          NVOIDP=NVOID

          DO KZ=1, NZ2
             DO JY=1, NY2
                DO IX=1, NX2
                   IXCO=INT((IX-1)/2.)+1
                   JYCO=INT((JY-1)/2.)+1
                   KZCO=INT((KZ-1)/2.)+1
                   MARCAP(IX,JY,KZ)=MARCA(IXCO,JYCO,KZCO)
                ENDDO
             ENDDO
          ENDDO

      !~~~~~~~~~~~~~~~~~~~~~~~
      ENDIF
      !~~~~~~~~~~~~~~~~~~~~~~~

       DEALLOCATE(UBAS)
       DEALLOCATE(U1CO, U1DMCO, U1GCO, U1SCO)
       DEALLOCATE(U2DMCO,U3DMCO,U4DMCO,U2GCO,U3GCO,U4GCO)
       DEALLOCATE(DIVERCO,DIVERGCO,DIVERDMCO)
       DEALLOCATE(MARCA)
       DEALLOCATE(XC,YC,ZC,GXC,GYC,GZC,UVOID,VOLNEW,UMEAN,MTOT,NCELLV,EPS,IP,REQ)
       DEALLOCATE(INDICE)
       DEALLOCATE(RADX,RADY,RADZ)

      !  IF (FLAG_WRITE_PIECES .EQ. 1) THEN
      !    DEALLOCATE(dUVOID,dRINIXCO,dRINIYCO,dRINIZCO,dRFINXCO,dRFINYCO,dRFINZCO)
      !    DEALLOCATE(dINICIOX,dINICIOY,dINICIOZ,dFINALX,dFINALY,dFINALZ)
      !  ENDIF


!########################################
!########################################
       ENDDO ! loop on all levels
!########################################
!########################################



       CLOSE(10)
       IF (FLAG_WRITE_CUBES .EQ. 1) CLOSE(12)
      !  IF (FLAG_WRITE_PIECES .EQ. 1) CLOSE(13)
       CLOSE(11)
       WRITE(*,*)
       WRITE(*,*)
       WRITE(*,*)
       WRITE(*,*) '////////////////////////////////////////////////////////////////'
       WRITE(*,*) '////////////////////////////////////////////////////////////////'
       WRITE(*,*) '//////////   ITERATION:', ITER, 'COMPLETED'
       WRITE(*,*) '////////////////////////////////////////////////////////////////'
       WRITE(*,*) '////////////////////////////////////////////////////////////////'
       WRITE(*,*)
       WRITE(*,*)

       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       !DEALLOCATE VOID FINDING VARIABLES
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       
      ! Void finding variables
       DEALLOCATE(MARCAP, REQP, MARCAP2)
       IF (ALLOCATED(DIVER)) DEALLOCATE(DIVER)  
       IF (ALLOCATED(DIVER0)) DEALLOCATE(DIVER0)
      
      ! Deallocate MASCLET input variables
       IF (FLAG_DATA .EQ. 0) THEN
         DEALLOCATE(NPATCH,PARE,NPART,NPARTST)
         DEALLOCATE(PATCHNX,PATCHNY,PATCHNZ)
         DEALLOCATE(PATCHX,PATCHY,PATCHZ)
         DEALLOCATE(PATCHRX,PATCHRY,PATCHRZ)
         DEALLOCATE(U1G,U1S,U1DM)
         DEALLOCATE(U2G,U3G,U4G)
         DEALLOCATE(U11G,U12G,U13G,U14G,U11S,U11DM)
         DEALLOCATE(U2PA,U3PA,U4PA,RXPA,RYPA,RZPA,MASAP)
       ENDIF

      !---------------------------------------
      ! Deallocate k-d tree variables (if needed)
      !--------------------------------------- 
       IF (FLAG_KD .EQ. 1) THEN
         DEALLOCATE(TREE)
         DEALLOCATE(TREEPOINTS)
         DEALLOCATE(HPART)
         DEALLOCATE(PART_DENS)
       ENDIF
      !---------------------------------------

!*////////////////////////////////////
!*////////////////////////////////////
      ENDDO !LOOP ON ITERATIONS: IFI=1, NITER
!*////////////////////////////////////
!*////////////////////////////////////

      DEALLOCATE(INICIOX,INICIOY,INICIOZ,FINALX,FINALY,FINALZ)
      DEALLOCATE(PERIODICAL)
      DEALLOCATE(ICX,ICY,ICZ)
      DEALLOCATE(VOL,FATHER)
      DEALLOCATE(RINIXCO,RINIYCO,RINIZCO)
      DEALLOCATE(RFINXCO,RFINYCO,RFINZCO)

      WRITE(*,*)
      WRITE(*,*) 
      WRITE(*,*) '************************************************'
      WRITE(*,*) '************************************************'
      WRITE(*,*) '      VOID-FINDING SUCCESSFULLY COMPLETED       '
      WRITE(*,*) '************************************************'
      WRITE(*,*) '************************************************'
      WRITE(*,*) '************************************************'



      !* END TIME
      call system_clock(clock_tf,trate,tmax)
      WRITE(*,*) 'TOTAL WALL TIME:',float(clock_tf-clock_t0)/1.e3,' seconds'

      CALL IDATE(DATE)
      CALL ITIME(TIME)
      WRITE(*,*) 'DATE=',DATE(1),'/',DATE(2),'/',DATE(3)
      WRITE(*,*) 'TIME=',TIME(1),':',TIME(2),':',TIME(3)

!-------------------------------------------------------
!-------------------------------------------------------
!-------------------------------------------------------
!-------------------------------------------------------
    END PROGRAM AVISM
!-------------------------------------------------------
!-------------------------------------------------------
!-------------------------------------------------------
!-------------------------------------------------------


!***********************************************************************
!*      SUBROUTINES IN EXTERNAL FILES                                  *
!***********************************************************************
!       Grid building
        INCLUDE 'grid.f90'
!       Reading data
        INCLUDE 'reader.f90'
#if use_hdf5 == 1
        INCLUDE 'reader_hdf5.f90'
#endif
!       Divergence calculation
        INCLUDE 'divergence.f90'
!       Passing from AMR to uniform grid
        INCLUDE 'amr2uniform.f90'
!       Void shapping
        INCLUDE 'shape.f90'
!       Basic numerical/array operations
        INCLUDE 'num.f90'
!       Name of output files
        INCLUDE 'nomfile.f90'