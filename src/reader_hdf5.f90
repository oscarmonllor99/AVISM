!***********************************************************************
SUBROUTINE READ_HDF5_SNAP(SIMTYPE,ITER,FILES_PER_SNAP,PARTTYPEX, &
                          MASSDM,ACHE,ZETA)
!***********************************************************************
!*  Reads particle data from Gadget-like HDF5 snapshots:
!*     SIMTYPE = 'AREPO'    (IllustrisTNG / Arepo)
!*     SIMTYPE = 'FLAMINGO' (SWIFT)
!*  Output (COMMONDATA): RXPA,RYPA,RZPA,U2PA,U3PA,U4PA,MASAP,NPARTT
!   in MASCLET internal units, 
!*  Positions in (-L/2, L/2).
!***********************************************************************
      use HDF5
      USE COMMONDATA
      implicit none

      !in / out
      CHARACTER(LEN=*), INTENT(IN) :: SIMTYPE
      INTEGER, INTENT(IN)  :: ITER            !snapshot number
      INTEGER, INTENT(IN)  :: FILES_PER_SNAP  !number of chunks
      INTEGER, INTENT(IN)  :: PARTTYPEX       !1 = gas, 2 = DM
      REAL*4,  INTENT(IN)  :: MASSDM          !DM particle mass [Msun]
      REAL*4,  INTENT(IN)  :: ACHE            !h = H0/100
      REAL*4,  INTENT(OUT) :: ZETA            !redshift

      !code-dependent configuration
      CHARACTER(LEN=100) :: PREFIX !snapshot naming prefix
      CHARACTER(LEN=8)   :: ITERFMT !iter formating specifier
      REAL*8  :: POSFAC     !file length unit -> cMpc
      REAL*8  :: MGASFAC    !file mass unit   -> Msun
      REAL*8  :: VELFAC     !file velocity    -> peculiar km/s
      LOGICAL :: VEL_SQRTA  !velocities stored as v/sqrt(a) !!! (true in Gadget/Arepo)

      !local
      INTEGER   :: IFILE, NTYPES, NFILES_HEADER
      INTEGER*8 :: I8, LOW1, LOW2, NPART_FILE
      CHARACTER(LEN=8)   :: ITER_STRING
      CHARACTER(LEN=32)  :: IFILE_STRING, PTYPE_STRING, GROUPNAME
      CHARACTER(LEN=200) :: FIL2

      !HID_T -> is the HDF5 type for file, group, dataset, attribute identifiers
      !HSIZE_T -> is the HDF5 type for dimensions and sizes
      INTEGER(HID_T) :: file_id, group_id, attr_id, dset_id, space_id
      INTEGER(HID_T) :: I8TYPE
      INTEGER        :: status
      INTEGER(HSIZE_T), DIMENSION(1) :: dims1d, maxdims1d
      INTEGER(HSIZE_T), DIMENSION(2) :: dims2d

      INTEGER*8, ALLOCATABLE :: NUMPART(:)
      REAL*8 :: ZETA8
      REAL*4, ALLOCATABLE :: SCR4(:), SCR42(:,:)

!-----------------------------------------------------------------------
!     Configuration per simulation code
!-----------------------------------------------------------------------

      !!! ADD NEW SIMULATION CODES HERE, if needed
      SELECT CASE (TRIM(SIMTYPE))
      CASE ('AREPO')
       PREFIX    = './simu_arepo/snap_'
       ITERFMT   = '(I3.3)'
       POSFAC    = 1.D-3/ACHE     !ckpc/h -> cMpc
       MGASFAC   = 1.D10/ACHE     !gas mass in 1e10 Msun/h -> Msun
       VEL_SQRTA = .TRUE.         !Velocities = v_pec / sqrt(a)
      CASE ('FLAMINGO')
       PREFIX    = './simu_flamingo/flamingo_'
       ITERFMT   = '(I4.4)'
       POSFAC    = 1.D0           !already cMpc (no h in SWIFT)!
       MGASFAC   = 1.D10          !gas mass in 1e10 Msun -> Msun (check /Units)
       VEL_SQRTA = .FALSE.        !peculiar km/s
      CASE DEFAULT
       WRITE(*,*) 'READ_HDF5_SNAP: unknown SIMTYPE = ', TRIM(SIMTYPE)
       STOP
      END SELECT

!-----------------------------------------------------------------------
!     Initialize the HDF5 library
!     HDF5 Fortran interface: required before using H5T_NATIVE_*,
!     H5F_ACC_RDONLY_F, h5kind_to_type...
!-----------------------------------------------------------------------
      CALL h5open_f(status)
      IF (status /= 0) THEN
       WRITE(*,*) 'READ_HDF5_SNAP: h5open_f failed, status =', status
       STOP
      END IF
      I8TYPE = h5kind_to_type(8, H5_INTEGER_KIND)

      IF (PARTTYPEX /= 1 .AND. PARTTYPEX /= 2) THEN
       WRITE(*,*) 'READ_HDF5_SNAP: unsupported PARTTYPEX =', PARTTYPEX
       WRITE(*,*) '  expected 1 (gas, /PartType0) or 2 (DM, /PartType1)'
       STOP
      END IF
      WRITE(PTYPE_STRING,'(I0)') PARTTYPEX-1
      GROUPNAME = '/PartType' // TRIM(PTYPE_STRING)

      !ALLOCATE + INITIALIZE PARTICLE ARRAYS
      ALLOCATE(U2PA(PARTIRED), U3PA(PARTIRED), U4PA(PARTIRED))
      ALLOCATE(MASAP(PARTIRED), RXPA(PARTIRED), RYPA(PARTIRED), &
               RZPA(PARTIRED))

      !$OMP PARALLEL DO SHARED(PARTIRED,U2PA,U3PA,U4PA, &
      !$OMP                    RXPA,RYPA,RZPA,MASAP), &
      !$OMP             PRIVATE(I8), DEFAULT(NONE)
      DO I8 = 1, PARTIRED
       U2PA(I8)  = 0.0
       U3PA(I8)  = 0.0
       U4PA(I8)  = 0.0
       RXPA(I8)  = 0.0
       RYPA(I8)  = 0.0
       RZPA(I8)  = 0.0
       MASAP(I8) = 0.0
      END DO

      WRITE(*,*) 'Simulation type:    ', TRIM(SIMTYPE)
      WRITE(*,*) 'Files per snapshot: ', FILES_PER_SNAP

      WRITE(ITER_STRING,ITERFMT) ITER

      LOW2 = 0
      !LOOP over snapshot chunks
      !###################################
      DO IFILE = 0, FILES_PER_SNAP-1
      !###################################

       !----------------------------------------------------------------
       !  File name
       !----------------------------------------------------------------
       FIL2 = TRIM(PREFIX) // TRIM(ITER_STRING)
       IF (FILES_PER_SNAP > 1) THEN
        WRITE(IFILE_STRING,'(I0)') IFILE
        FIL2 = TRIM(FIL2) // '.' // TRIM(IFILE_STRING)
       END IF
       FIL2 = TRIM(FIL2) // '.hdf5'

       WRITE(*,*) 'Reading iteration file: ', ITER, ' ', TRIM(FIL2)

       CALL h5fopen_f(FIL2, H5F_ACC_RDONLY_F, file_id, status)
       IF (status /= 0) THEN
        WRITE(*,*) 'READ_HDF5_SNAP: cannot open file ', TRIM(FIL2)
        STOP
       END IF

       !----------------------------------------------------------------
       !  Header
       !----------------------------------------------------------------
       CALL h5gopen_f(file_id, "/Header", group_id, status)
       IF (status /= 0) THEN
        WRITE(*,*) 'READ_HDF5_SNAP: cannot open /Header in ', TRIM(FIL2)
        STOP
       END IF

       ! NumPart_ThisFile: int32 x6 (Arepo) or int64 x6/x7 (SWIFT).
       CALL h5aopen_f(group_id, "NumPart_ThisFile", attr_id, status)
       IF (status /= 0) THEN
        WRITE(*,*) 'READ_HDF5_SNAP: no NumPart_ThisFile in ', TRIM(FIL2)
        STOP
       END IF
       !get attribute rank and dimensions
       CALL h5aget_space_f(attr_id, space_id, status)
       CALL h5sget_simple_extent_dims_f(space_id, dims1d, maxdims1d, &
                                        status)  ! status = rank or -1
       IF (status < 0) THEN
        WRITE(*,*) 'READ_HDF5_SNAP: bad NumPart_ThisFile in ', TRIM(FIL2)
        STOP
       END IF
       CALL h5sclose_f(space_id, status)

       !check particle types available
       NTYPES = INT(dims1d(1))
       IF (PARTTYPEX > NTYPES) THEN
        WRITE(*,*) 'READ_HDF5_SNAP: NumPart_ThisFile has only', NTYPES, &
                   ' entries'
        STOP
       END IF

       IF (ALLOCATED(NUMPART)) DEALLOCATE(NUMPART)
       ALLOCATE(NUMPART(NTYPES))
       CALL h5aread_f(attr_id, I8TYPE, NUMPART, dims1d, status)
       CALL h5aclose_f(attr_id, status)
       IF (status /= 0) THEN
        WRITE(*,*) 'READ_HDF5_SNAP: cannot read NumPart_ThisFile in ', &
                   TRIM(FIL2)
        STOP
       END IF
       NPART_FILE = NUMPART(PARTTYPEX)

       ! Redshift: float64 in both codes, read as double, then convert
       dims1d(1) = 1
       CALL h5aopen_f(group_id, "Redshift", attr_id, status)
       IF (status == 0) THEN
        CALL h5aread_f(attr_id, H5T_NATIVE_DOUBLE, ZETA8, dims1d, status)
        CALL h5aclose_f(attr_id, status)
       END IF
       IF (status /= 0) THEN
        WRITE(*,*) 'READ_HDF5_SNAP: cannot read Redshift in ', TRIM(FIL2)
        STOP
       END IF
       ZETA = REAL(ZETA8, KIND=4)

       !check consitency of FILES_PER_SNAP with header
       IF (IFILE == 0) THEN
        CALL h5aopen_f(group_id, "NumFilesPerSnapshot", attr_id, status)
        IF (status == 0) THEN
         dims1d(1) = 1
         CALL h5aread_f(attr_id, H5T_NATIVE_INTEGER, NFILES_HEADER, &
                        dims1d, status)
         IF (status == 0 .AND. NFILES_HEADER /= FILES_PER_SNAP) THEN
          WRITE(*,*) ' *** WARNING: FILES_PER_SNAP =', FILES_PER_SNAP, &
                     ' but header says NumFilesPerSnapshot =', &
                     NFILES_HEADER
         END IF
         CALL h5aclose_f(attr_id, status)
        END IF
       END IF

       CALL h5gclose_f(group_id, status)

       !----------------------------------------------------------------
       !  Particle block
       !----------------------------------------------------------------
       WRITE(*,*) NPART_FILE, 'particles'

       IF (NPART_FILE <= 0) THEN
        ! Chunks with no particles of this type have no group at all
        CALL h5fclose_f(file_id, status)
        CYCLE
       END IF

       LOW1 = LOW2 + 1
       LOW2 = LOW1 + NPART_FILE - 1

       IF (LOW2 > PARTIRED) THEN
        WRITE(*,*) 'READ_HDF5_SNAP: particle buffer too small.'
        WRITE(*,*) '  PARTIRED           =', PARTIRED
        WRITE(*,*) '  needed at least    =', LOW2
        WRITE(*,*) '  while reading file  ', TRIM(FIL2)
        CALL h5fclose_f(file_id, status)
        STOP
       END IF

       CALL h5gopen_f(file_id, TRIM(GROUPNAME), group_id, status)
       IF (status /= 0) THEN
        WRITE(*,*) 'READ_HDF5_SNAP: cannot open group ', &
                   TRIM(GROUPNAME), ' in ', TRIM(FIL2)
        STOP
       END IF

       ! Memory buffers set to REAL*4 with H5T_NATIVE_REAL
       ! Converts from whatever is on disk (float32 or float64), hence, no type error.
       ! (3,N) memory layout for Fortran-friendlyness. 
       dims2d(1) = 3
       dims2d(2) = NPART_FILE
       ALLOCATE(SCR42(3,NPART_FILE))

       CALL h5dopen_f(group_id, "Coordinates", dset_id, status)
       IF (status == 0) THEN
        CALL h5dread_f(dset_id, H5T_NATIVE_REAL, SCR42, dims2d, status)
        CALL h5dclose_f(dset_id, status)
       END IF
       IF (status /= 0) THEN
        WRITE(*,*) 'READ_HDF5_SNAP: cannot read Coordinates in ', &
                   TRIM(FIL2)
        STOP
       END IF
       RXPA(LOW1:LOW2) = SCR42(1,1:NPART_FILE)
       RYPA(LOW1:LOW2) = SCR42(2,1:NPART_FILE)
       RZPA(LOW1:LOW2) = SCR42(3,1:NPART_FILE)

       CALL h5dopen_f(group_id, "Velocities", dset_id, status)
       IF (status == 0) THEN
        CALL h5dread_f(dset_id, H5T_NATIVE_REAL, SCR42, dims2d, status)
        CALL h5dclose_f(dset_id, status)
       END IF
       IF (status /= 0) THEN
        WRITE(*,*) 'READ_HDF5_SNAP: cannot read Velocities in ', &
                   TRIM(FIL2)
        STOP
       END IF
       U2PA(LOW1:LOW2) = SCR42(1,1:NPART_FILE)
       U3PA(LOW1:LOW2) = SCR42(2,1:NPART_FILE)
       U4PA(LOW1:LOW2) = SCR42(3,1:NPART_FILE)

       DEALLOCATE(SCR42)

       IF (PARTTYPEX == 1) THEN
        ! Gas: individual masses, in file units for now
        dims1d(1) = NPART_FILE
        ALLOCATE(SCR4(NPART_FILE))
        CALL h5dopen_f(group_id, "Masses", dset_id, status)
        IF (status == 0) THEN
         CALL h5dread_f(dset_id, H5T_NATIVE_REAL, SCR4, dims1d, status)
         CALL h5dclose_f(dset_id, status)
        END IF
        IF (status /= 0) THEN
         WRITE(*,*) 'READ_HDF5_SNAP: cannot read Masses in ', TRIM(FIL2)
         STOP
        END IF
        MASAP(LOW1:LOW2) = SCR4(1:NPART_FILE)
        DEALLOCATE(SCR4)
       ELSE
        ! DM: same mass for every particle, already in Msun
        MASAP(LOW1:LOW2) = MASSDM
       END IF

       CALL h5gclose_f(group_id, status)
       CALL h5fclose_f(file_id, status)

      END DO

      IF (ALLOCATED(NUMPART)) DEALLOCATE(NUMPART)

      NPARTT = LOW2
      WRITE(*,*) '     TOTAL PARTICLES IN ITER=', NPARTT

      IF (NPARTT <= 0) THEN
       WRITE(*,*) 'READ_HDF5_SNAP: no particles of type ', PARTTYPEX, &
                  ' found in snapshot ', ITER
       RETURN
      END IF

!-----------------------------------------------------------------------
!     Unit conversions
!-----------------------------------------------------------------------
      ! file length unit -> cMpc, and (0,L) -> (-L/2,L/2)
      RXPA(1:NPARTT) = RXPA(1:NPARTT)*REAL(POSFAC,4) - LADO0/2.0
      RYPA(1:NPARTT) = RYPA(1:NPARTT)*REAL(POSFAC,4) - LADO0/2.0
      RZPA(1:NPARTT) = RZPA(1:NPARTT)*REAL(POSFAC,4) - LADO0/2.0

      ! gas masses: file units -> Msun
      IF (PARTTYPEX == 1) THEN
       MASAP(1:NPARTT) = MASAP(1:NPARTT)*REAL(MGASFAC,4)
      END IF

      ! velocities: -> peculiar km/s
      VELFAC = 1.D0
      IF (VEL_SQRTA) VELFAC = SQRT(1.D0/(1.D0+ZETA8))

      ! Msun and km/s -> MASCLET internal units
      MASAP(1:NPARTT) = MASAP(1:NPARTT)/UM
      U2PA(1:NPARTT)  = U2PA(1:NPARTT)*REAL(VELFAC,4)/UV
      U3PA(1:NPARTT)  = U3PA(1:NPARTT)*REAL(VELFAC,4)/UV
      U4PA(1:NPARTT)  = U4PA(1:NPARTT)*REAL(VELFAC,4)/UV

      WRITE(*,*) 'X range:', MINVAL(RXPA(1:NPARTT)), MAXVAL(RXPA(1:NPARTT))
      WRITE(*,*) 'Y range:', MINVAL(RYPA(1:NPARTT)), MAXVAL(RYPA(1:NPARTT))
      WRITE(*,*) 'Z range:', MINVAL(RZPA(1:NPARTT)), MAXVAL(RZPA(1:NPARTT))

      RETURN

!***********************************************************************
      END SUBROUTINE READ_HDF5_SNAP
!***********************************************************************
