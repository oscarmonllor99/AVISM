
!***********************************************************************
SUBROUTINE READ_AREPO_HDF5(ITER, FILES_PER_SNAP,PARTTYPEX,MASSDM,ACHE,ZETA)
!***********************************************************************
!*       Reads particle data of the simulation
!***********************************************************************
      use HDF5
      USE COMMONDATA
      implicit none 

      !in
      INTEGER, INTENT(IN)  :: ITER            ! snapshot number
      INTEGER, INTENT(IN)  :: FILES_PER_SNAP  ! number of chunks
      INTEGER, INTENT(IN)  :: PARTTYPEX       ! 1 = gas, 2 = DM
      REAL*4,  INTENT(IN)  :: MASSDM          ! DM particle mass [Msun]
      REAL*4,  INTENT(IN)  :: ACHE            ! h = H0/100
      REAL*4,  INTENT(OUT) :: ZETA            ! redshift

      !local
      INTEGER :: IFILE
      INTEGER*8 :: I8
      INTEGER*8 :: LOW1, LOW2
      INTEGER   :: NPART_FILE
      INTEGER   :: NFILES_HEADER

      CHARACTER(LEN=3)   :: ITER_STRING
      CHARACTER(LEN=32)  :: IFILE_STRING
      CHARACTER(LEN=32)  :: PTYPE_STRING
      CHARACTER(LEN=200) :: FIL1, FIL2
      CHARACTER(LEN=32)  :: GROUPNAME

      INTEGER(HID_T) :: file_id, group_id, attr_id, dset_id
      INTEGER        :: status

      INTEGER, DIMENSION(6) :: NumPart_ThisFile
      INTEGER(HSIZE_T), DIMENSION(1) :: dims1d
      INTEGER(HSIZE_T), DIMENSION(2) :: dims2d

      REAL*8 :: ZETA8

      REAL*4, ALLOCATABLE :: SCR4(:)
      REAL*4, ALLOCATABLE :: SCR42(:,:)

 
!     Initialize the HDF5 Fortran interface.
!     This is NOT optional here.  H5T_NATIVE_REAL, H5T_NATIVE_DOUBLE and
!     H5T_NATIVE_INTEGER are module variables of the HDF5 Fortran
!     bindings, and h5open_f is what assigns them their real datatype
!     IDs.  Without this call they are zero, and every h5dread_f /
      CALL h5open_f(status)
      IF (status /= 0) THEN
       WRITE(*,*) 'READ_AREPO_HDF5: h5open_f failed, status =', status
       STOP
      END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !PARTTYPEX check!
      IF (PARTTYPEX /= 1 .AND. PARTTYPEX /= 2) THEN
       WRITE(*,*) 'READ_AREPO_HDF5: unsupported PARTTYPEX =', PARTTYPEX
       WRITE(*,*) '  expected 1 (gas, /PartType0) or 2 (DM, /PartType1)'
       STOP
      END IF

      ! /PartType0 for PARTTYPEX=1, /PartType1 for PARTTYPEX=2.
      ! (NumPart_ThisFile is indexed 1..6 in Fortran, hence the offset.)
      WRITE(PTYPE_STRING,'(I0)') PARTTYPEX-1
      GROUPNAME = '/PartType' // TRIM(PTYPE_STRING)

      !ALLOCATE PARTICLE ARRAYS
      ALLOCATE(U2PA(PARTIRED), U3PA(PARTIRED), U4PA(PARTIRED))
      ALLOCATE(MASAP(PARTIRED), RXPA(PARTIRED), RYPA(PARTIRED), RZPA(PARTIRED))

      !INITIALIZE PARTICLE ARRAYS
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

      WRITE(*,*) 'Files per snapshot: ', FILES_PER_SNAP


      LOW2 = 0

      !LOOP over snapshot chunks (IFILE = 0..FILES_PER_SNAP-1)
      DO IFILE = 0, FILES_PER_SNAP-1

       !----------------------------------------------------------------
       !  File name
       !----------------------------------------------------------------
       WRITE(ITER_STRING,'(I3.3)') ITER
       FIL1 = './simu_arepo/snap_' // ITER_STRING

       IF (FILES_PER_SNAP == 1) THEN
        FIL2 = TRIM(FIL1)
       ELSE
        WRITE(IFILE_STRING,'(I0)') IFILE
        FIL2 = TRIM(FIL1) // '.' // TRIM(IFILE_STRING)
       END IF
       FIL2 = TRIM(FIL2) // '.hdf5'

       WRITE(*,*) 'Reading iteration file: ', ITER, ' ', TRIM(FIL2)

       CALL h5fopen_f(FIL2, H5F_ACC_RDONLY_F, file_id, status)
       IF (status /= 0) THEN
        WRITE(*,*) 'READ_AREPO_HDF5: cannot open file ', TRIM(FIL2)
        STOP
       END IF

       !----------------------------------------------------------------
       !  Header
       !----------------------------------------------------------------
       CALL h5gopen_f(file_id, "/Header", group_id, status)
       IF (status /= 0) THEN
        WRITE(*,*) 'READ_AREPO_HDF5: cannot open /Header in ', TRIM(FIL2)
        STOP
       END IF

       ! NumPart_ThisFile: 6 x 32-bit integer
       dims1d(1) = 6
       CALL h5aopen_f(group_id, "NumPart_ThisFile", attr_id, status)
       IF (status == 0) THEN
        CALL h5aread_f(attr_id, H5T_NATIVE_INTEGER, NumPart_ThisFile, &
                       dims1d, status)
        CALL h5aclose_f(attr_id, status)
       END IF
       IF (status /= 0) THEN
        WRITE(*,*) 'READ_AREPO_HDF5: cannot read NumPart_ThisFile in ', &
                   TRIM(FIL2)
        STOP
       END IF

       ! Redshift: stored as float64 in TNG.  Read it as a double and
       ! then narrow it, never straight into the REAL*4 ZETA.
       dims1d(1) = 1
       CALL h5aopen_f(group_id, "Redshift", attr_id, status)
       IF (status == 0) THEN
        CALL h5aread_f(attr_id, H5T_NATIVE_DOUBLE, ZETA8, dims1d, status)
        CALL h5aclose_f(attr_id, status)
       END IF
       IF (status /= 0) THEN
        WRITE(*,*) 'READ_AREPO_HDF5: cannot read Redshift in ', TRIM(FIL2)
        STOP
       END IF
       ZETA = REAL(ZETA8, KIND=4)

       ! Soft cross-check of the chunk count against the header, so a
       ! wrong FILES_PER_SNAP shows up as a clear warning rather than as
       ! a silently truncated particle list.
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
       NPART_FILE = NumPart_ThisFile(PARTTYPEX)
       WRITE(*,*) NPART_FILE, 'particles'

       IF (NPART_FILE <= 0) THEN
        ! Chunks with no particles of this type have no group at all.
        CALL h5fclose_f(file_id, status)
        CYCLE
       END IF

       LOW1 = LOW2 + 1
       LOW2 = LOW1 + NPART_FILE - 1

       IF (LOW2 > PARTIRED) THEN
        WRITE(*,*) 'READ_AREPO_HDF5: particle buffer too small.'
        WRITE(*,*) '  PARTIRED           =', PARTIRED
        WRITE(*,*) '  needed at least    =', LOW2
        WRITE(*,*) '  while reading file  ', TRIM(FIL2)
        CALL h5fclose_f(file_id, status)
        STOP
       END IF

       CALL h5gopen_f(file_id, TRIM(GROUPNAME), group_id, status)
       IF (status /= 0) THEN
        WRITE(*,*) 'READ_AREPO_HDF5: cannot open group ', &
                   TRIM(GROUPNAME), ' in ', TRIM(FIL2)
        STOP
       END IF

       ! Coordinates / Velocities memory layout is (3,N)
       dims2d(1) = 3
       dims2d(2) = NPART_FILE
       ALLOCATE(SCR42(3,NPART_FILE))

       CALL h5dopen_f(group_id, "Coordinates", dset_id, status)
       IF (status == 0) THEN
        CALL h5dread_f(dset_id, H5T_NATIVE_REAL, SCR42, dims2d, status)
        CALL h5dclose_f(dset_id, status)
       END IF
       IF (status /= 0) THEN
        WRITE(*,*) 'READ_AREPO_HDF5: cannot read Coordinates in ', &
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
        WRITE(*,*) 'READ_AREPO_HDF5: cannot read Velocities in ', &
                   TRIM(FIL2)
        STOP
       END IF
       U2PA(LOW1:LOW2) = SCR42(1,1:NPART_FILE)
       U3PA(LOW1:LOW2) = SCR42(2,1:NPART_FILE)
       U4PA(LOW1:LOW2) = SCR42(3,1:NPART_FILE)

       DEALLOCATE(SCR42)

       IF (PARTTYPEX == 1) THEN
        ! Gas: individual masses, in code units for now
        dims1d(1) = NPART_FILE
        ALLOCATE(SCR4(NPART_FILE))
        CALL h5dopen_f(group_id, "Masses", dset_id, status)
        IF (status == 0) THEN
         CALL h5dread_f(dset_id, H5T_NATIVE_REAL, SCR4, dims1d, status)
         CALL h5dclose_f(dset_id, status)
        END IF
        IF (status /= 0) THEN
         WRITE(*,*) 'READ_AREPO_HDF5: cannot read Masses in ', TRIM(FIL2)
         STOP
        END IF
        MASAP(LOW1:LOW2) = SCR4(1:NPART_FILE)
        DEALLOCATE(SCR4)

       ELSE
        ! DM: single mass for every particle, already in Msun.
        MASAP(LOW1:LOW2) = MASSDM
       END IF

       CALL h5gclose_f(group_id, status)
       CALL h5fclose_f(file_id, status)

      END DO

      NPARTT = LOW2
      WRITE(*,*) '     TOTAL PARTICLES IN ITER=', NPARTT

      IF (NPARTT <= 0) THEN
       WRITE(*,*) 'READ_AREPO_HDF5: no particles of type ', PARTTYPEX, &
                  ' found in snapshot ', ITER
       RETURN
      END IF


      !Unit conversions
      ! ckpc/h -> cMpc, and (0,L) -> (-L/2,L/2)
      RXPA(1:NPARTT) = RXPA(1:NPARTT)*1E-3/ACHE - LADO0/2.0
      RYPA(1:NPARTT) = RYPA(1:NPARTT)*1E-3/ACHE - LADO0/2.0
      RZPA(1:NPARTT) = RZPA(1:NPARTT)*1E-3/ACHE - LADO0/2.0

      ! Gas masses: code units (1e10 Msun/h) -> Msun.
      IF (PARTTYPEX == 1) THEN
       MASAP(1:NPARTT) = MASAP(1:NPARTT)*1E10/ACHE
      END IF

      ! Msun and km/s -> MASCLET internal units
      MASAP(1:NPARTT) = MASAP(1:NPARTT)/UM
      U2PA(1:NPARTT)  = U2PA(1:NPARTT)/UV
      U3PA(1:NPARTT)  = U3PA(1:NPARTT)/UV
      U4PA(1:NPARTT)  = U4PA(1:NPARTT)/UV

      RETURN

!***********************************************************************
END SUBROUTINE READ_AREPO_HDF5
!***********************************************************************


!***********************************************************************
SUBROUTINE READ_FLAMINGO_HDF5(ITER,FILES_PER_SNAP,PARTTYPEX,MASSDM,ACHE,ZETA)
!***********************************************************************
!*       Reads particle data of the simulation
!***********************************************************************
      use HDF5
      USE COMMONDATA
      implicit none 

      integer iter, files_per_snap, i
      integer*8 :: I8
      CHARACTER*4 ITER_STRING
      CHARACTER*200 IFILE_STRING
      INTEGER IFILE
      CHARACTER*200 FIL1,FIL2
      INTEGER PARTTYPEX
      REAL*4 ACHE,MASSDM !mass of DM particles in Msun
      REAL*4 :: ZETA

      integer(hid_t) :: file_id, group_id, attr_id, mem_space_id, file_space_id
      INTEGER(HID_T) :: memtype_id
      integer :: status
      integer(KIND=8), dimension(6) :: NumPart_ThisFile
      integer(hsize_t), dimension(1) :: dims1d
      integer(hsize_t), dimension(2) :: dims2d
      
      integer(KIND=8) :: LOW1, LOW2
      REAL*8,ALLOCATABLE::SCR4(:)
      REAL*8,ALLOCATABLE::SCR82(:,:)
      REAL*4,ALLOCATABLE::SCR42(:,:)


      !ALLOCATE PARTICLE ARRAYS
      ALLOCATE(U2PA(PARTIRED), U3PA(PARTIRED), U4PA(PARTIRED))
      ALLOCATE(MASAP(PARTIRED), RXPA(PARTIRED), RYPA(PARTIRED), RZPA(PARTIRED))

      !INITIALIZE PARTICLE ARRAYS
      !$OMP PARALLEL DO SHARED(PARTIRED,U2PA,U3PA,U4PA,RXPA,RYPA,RZPA, &
      !$OMP            MASAP), &
      !$OMP            PRIVATE(I8)
      DO I8=1,PARTIRED
       U2PA(I8)=0.0 
       U3PA(I8)=0.0
       U4PA(I8)=0.0
       RXPA(I8)=0.0  !DM vars
       RYPA(I8)=0.0
       RZPA(I8)=0.0
       MASAP(I8)=0.0
      END DO

      WRITE(*,*) 'Files per snapshot: ', FILES_PER_SNAP

      LOW2=0
      !###################################
      DO IFILE=0,FILES_PER_SNAP-1 
      !###################################

       !*     READING DATA
       WRITE(ITER_STRING, '(I4.4)') ITER
       FIL1 = './simu_flamingo/flamingo_' // ITER_STRING
       IF (FILES_PER_SNAP .EQ. 1) THEN
              FIL2 = FIL1
       ELSE
              WRITE(IFILE_STRING, '(I3)') IFILE
              FIL2 = TRIM(ADJUSTL(FIL1)) // '.' // TRIM(ADJUSTL(IFILE_STRING))
       END IF
       FIL2 = TRIM(ADJUSTL(FIL2)) // '.hdf5'

       ! Open the HDF5 file in read-only mode
       WRITE(*,*) 'Reading iteration file: ',ITER,' ', &             
                            TRIM(ADJUSTL(FIL2))

       CALL h5fopen_f(FIL2, H5F_ACC_RDONLY_F, file_id, status)
       IF (status /= 0) THEN
              PRINT *, "Error opening file: ", FIL2
              STOP
       END IF

       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       !header
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       CALL h5gopen_f(file_id, "/Header", group_id, status)
       !READ numpart
       CALL h5aopen_f(group_id, "NumPart_ThisFile", &                   
                            attr_id, status)
       dims1d(1) = 6
       CALL h5aget_type_f(attr_id, memtype_id, status)
       CALL h5aread_f(attr_id, memtype_id, NumPart_ThisFile,&                   
                            dims1d, status)

       WRITE(*,*) NumPart_ThisFile(PARTTYPEX), 'particles'
       LOW1=LOW2+1
       LOW2=LOW1+NumPart_ThisFile(PARTTYPEX)-1

       CALL h5aclose_f(attr_id, status)

       !READ ZETA
       CALL h5aopen_f(group_id, "Redshift", attr_id, status)
       CALL h5aget_type_f(attr_id, memtype_id, status)
       dims1d(1) = 1
       CALL h5aread_f(attr_id, memtype_id, ZETA, dims1d, status)

       CALL h5aclose_f(attr_id, status)
       CALL h5gclose_f(group_id, status)
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       dims1d(1) = NumPart_ThisFile(PARTTYPEX)
       dims2d(1) = NumPart_ThisFile(PARTTYPEX)
       dims2d(2) = 3

       CALL h5gopen_f(file_id, "/PartType1", group_id, status)
       if (status /= 0) then
              PRINT *, "Error opening group: /PartType1"
              CALL h5fclose_f(file_id, status)
              STOP
       end if
       
       ALLOCATE(SCR82(3,NumPart_ThisFile(PARTTYPEX)))

       ! WRITE(*,*) 'Reading positions ...'
       CALL h5dopen_f(group_id, "Coordinates", attr_id, status)
       CALL h5dget_type_f(attr_id, memtype_id, status)
       CALL h5dread_f(attr_id, memtype_id, SCR82, dims2d, status)
       RXPA(LOW1:LOW2)=SCR82(1,1:NumPart_ThisFile(PARTTYPEX))
       RYPA(LOW1:LOW2)=SCR82(2,1:NumPart_ThisFile(PARTTYPEX))
       RZPA(LOW1:LOW2)=SCR82(3,1:NumPart_ThisFile(PARTTYPEX))
       CALL h5dclose_f(attr_id, status)

       ALLOCATE(SCR42(3,NumPart_ThisFile(PARTTYPEX)))

       ! WRITE(*,*) 'Reading velocities ...'
       CALL h5dopen_f(group_id, "Velocities", attr_id, status)
       CALL h5dget_type_f(attr_id, memtype_id, status)
       CALL h5dread_f(attr_id, memtype_id, SCR42, dims2d, status)
       U2PA(LOW1:LOW2)=SCR42(1,1:NumPart_ThisFile(PARTTYPEX))
       U3PA(LOW1:LOW2)=SCR42(2,1:NumPart_ThisFile(PARTTYPEX))
       U4PA(LOW1:LOW2)=SCR42(3,1:NumPart_ThisFile(PARTTYPEX))
       CALL h5dclose_f(attr_id, status)

       DEALLOCATE(SCR82)
       DEALLOCATE(SCR42)

       IF (PARTTYPEX .EQ. 1) THEN
          ALLOCATE(SCR4(NumPart_ThisFile(PARTTYPEX)))
       
          CALL h5dopen_f(group_id, "Masses", attr_id, status)
          CALL h5dget_type_f(attr_id, memtype_id, status)
          CALL h5dread_f(attr_id, memtype_id, SCR4, dims1d, status)
          MASAP(LOW1:LOW2)=SCR4(1:NumPart_ThisFile(PARTTYPEX))
          CALL h5dclose_f(attr_id, status)

          DEALLOCATE(SCR4)

       ELSE IF (PARTTYPEX .EQ. 2) THEN
          MASAP = MASSDM ! mass of DM particles in Msun
       ENDIF
              
       CALL h5gclose_f(group_id, status)
       CALL h5fclose_f(file_id, status)

      !###################################
      END DO 
      !###################################

      NPARTT = LOW2
      WRITE(*,*) '     TOTAL PARTICLES IN ITER=', NPARTT

      ! From (0,L) to (-L/2,L/2) and from Mpc/h to Mpc
      RXPA = RXPA - LADO0/2.0
      RYPA = RYPA - LADO0/2.0
      RZPA = RZPA - LADO0/2.0

      WRITE(*,*) minval(RXPA), maxval(RXPA)
      WRITE(*,*) minval(RYPA), maxval(RYPA)
      WRITE(*,*) minval(RZPA), maxval(RZPA)

      !Now to masclet units
      MASAP = MASAP/UM
      U2PA = U2PA/UV
      U3PA = U3PA/UV
      U4PA = U4PA/UV

      WRITE(*,*) minval(MASAP), maxval(MASAP)
      WRITE(*,*) minval(U2PA), maxval(U2PA)
      WRITE(*,*) minval(U3PA), maxval(U3PA)
      WRITE(*,*) minval(U4PA), maxval(U4PA)

      RETURN

!***********************************************************************
END SUBROUTINE READ_FLAMINGO_HDF5
!***********************************************************************
