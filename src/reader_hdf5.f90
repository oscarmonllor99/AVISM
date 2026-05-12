
!***********************************************************************
SUBROUTINE READ_AREPO_HDF5(ITER, FILES_PER_SNAP,PARTTYPEX,MASSDM,ACHE,ZETA)
!***********************************************************************
!*       Reads particle data of the simulation
!***********************************************************************
      use HDF5
      USE COMMONDATA
      implicit none 

      integer iter, files_per_snap, i
      integer*8 :: I8
      CHARACTER*3 ITER_STRING
      CHARACTER*200 IFILE_STRING
      INTEGER IFILE
      CHARACTER*200 FIL1,FIL2
      INTEGER PARTTYPEX !parttype to read, 0 for gas, 1 for DM
      REAL*4 ACHE,MASSDM !mass of DM particles in Msun
      REAL*4 :: ZETA
      REAL*8 :: ZETA8

      integer(hid_t) :: file_id, group_id, attr_id, mem_space_id, file_space_id
      INTEGER(HID_T) :: memtype_id
      integer :: status
      integer, dimension(6) :: NumPart_ThisFile
      integer(hsize_t), dimension(1) :: dims1d
      integer(hsize_t), dimension(2) :: dims2d
      
      integer(KIND=8) :: LOW1, LOW2
      REAL*4,ALLOCATABLE::SCR4(:)
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

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      IF (PARTTYPEX .EQ. 1) THEN
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
       LOW2=0
       !###################################
       DO IFILE=0,FILES_PER_SNAP-1 
       !###################################

       !*     READING DATA
       WRITE(ITER_STRING, '(I3.3)') ITER
       FIL1 = './simu_arepo/snap_' // ITER_STRING
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
       CALL h5aread_f(attr_id, memtype_id, ZETA8, dims1d, status)
       ZETA = ZETA8
       CALL h5aclose_f(attr_id, status)

       CALL h5gclose_f(group_id, status)
       
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       dims1d(1) = NumPart_ThisFile(PARTTYPEX)
       dims2d(1) = NumPart_ThisFile(PARTTYPEX)
       dims2d(2) = 3

       CALL h5gopen_f(file_id, "/PartType0", group_id, status)
       if (status /= 0) then
              PRINT *, "Error opening group: /PartType0"
              CALL h5fclose_f(file_id, status)
              STOP
       end if
       
       ALLOCATE(SCR42(3,NumPart_ThisFile(PARTTYPEX)))

       ! WRITE(*,*) 'Reading positions ...'
       CALL h5dopen_f(group_id, "Coordinates", attr_id, status)
       CALL h5dget_type_f(attr_id, memtype_id, status)
       CALL h5dread_f(attr_id, memtype_id, SCR42, dims2d, status)
       RXPA(LOW1:LOW2)=SCR42(1,1:NumPart_ThisFile(PARTTYPEX))
       RYPA(LOW1:LOW2)=SCR42(2,1:NumPart_ThisFile(PARTTYPEX))
       RZPA(LOW1:LOW2)=SCR42(3,1:NumPart_ThisFile(PARTTYPEX))
       CALL h5dclose_f(attr_id, status)

       ! WRITE(*,*) 'Reading velocities ...'
       CALL h5dopen_f(group_id, "Velocities", attr_id, status)
       CALL h5dget_type_f(attr_id, memtype_id, status)
       CALL h5dread_f(attr_id, memtype_id, SCR42, dims2d, status)
       U2PA(LOW1:LOW2)=SCR42(1,1:NumPart_ThisFile(PARTTYPEX))
       U3PA(LOW1:LOW2)=SCR42(2,1:NumPart_ThisFile(PARTTYPEX))
       U4PA(LOW1:LOW2)=SCR42(3,1:NumPart_ThisFile(PARTTYPEX))
       CALL h5dclose_f(attr_id, status)

       DEALLOCATE(SCR42)

       ALLOCATE(SCR4(NumPart_ThisFile(PARTTYPEX)))
       ! WRITE(*,*) 'Reading masses ...'
       CALL h5dopen_f(group_id, "Masses", attr_id, status)
       CALL h5dget_type_f(attr_id, memtype_id, status)
       CALL h5dread_f(attr_id, memtype_id, SCR4, dims1d, status)
       MASAP(LOW1:LOW2)=SCR4(1:NumPart_ThisFile(PARTTYPEX))
       CALL h5dclose_f(attr_id, status)
       DEALLOCATE(SCR4)
              
       CALL h5gclose_f(group_id, status)
       CALL h5fclose_f(file_id, status)

       !###################################
       END DO 
       !###################################

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ELSE IF (PARTTYPEX .EQ. 2) THEN
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

       LOW2=0
       !###################################
       DO IFILE=0,FILES_PER_SNAP-1 
       !###################################

       !*     READING DATA
       WRITE(ITER_STRING, '(I3.3)') ITER
       FIL1 = './simu_arepo/snap_' // ITER_STRING
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
       
       ALLOCATE(SCR42(3,NumPart_ThisFile(PARTTYPEX)))

       ! WRITE(*,*) 'Reading positions ...'
       CALL h5dopen_f(group_id, "Coordinates", attr_id, status)
       CALL h5dget_type_f(attr_id, memtype_id, status)
       CALL h5dread_f(attr_id, memtype_id, SCR42, dims2d, status)
       RXPA(LOW1:LOW2)=SCR42(1,1:NumPart_ThisFile(PARTTYPEX))
       RYPA(LOW1:LOW2)=SCR42(2,1:NumPart_ThisFile(PARTTYPEX))
       RZPA(LOW1:LOW2)=SCR42(3,1:NumPart_ThisFile(PARTTYPEX))
       CALL h5dclose_f(attr_id, status)

       ! WRITE(*,*) 'Reading velocities ...'
       CALL h5dopen_f(group_id, "Velocities", attr_id, status)
       CALL h5dget_type_f(attr_id, memtype_id, status)
       CALL h5dread_f(attr_id, memtype_id, SCR42, dims2d, status)
       U2PA(LOW1:LOW2)=SCR42(1,1:NumPart_ThisFile(PARTTYPEX))
       U3PA(LOW1:LOW2)=SCR42(2,1:NumPart_ThisFile(PARTTYPEX))
       U4PA(LOW1:LOW2)=SCR42(3,1:NumPart_ThisFile(PARTTYPEX))
       CALL h5dclose_f(attr_id, status)

       DEALLOCATE(SCR42)

       ! WRITE(*,*) 'Assigning masses ...'
       MASAP = MASSDM ! mass of DM particles in Msun
              
       CALL h5gclose_f(group_id, status)
       CALL h5fclose_f(file_id, status)

       !###################################
       END DO 
       !###################################

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~       
      ENDIF 
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      NPARTT = LOW2
      WRITE(*,*) '     TOTAL PARTICLES IN ITER=', NPARTT

      ! From (0,L) to (-L/2,L/2) and from Mpc/h to Mpc
      RXPA = RXPA*1E-3/ACHE - LADO0/2.0
      RYPA = RYPA*1E-3/ACHE - LADO0/2.0
      RZPA = RZPA*1E-3/ACHE - LADO0/2.0

      !gas mass convert to M_sun
      IF (PARTTYPEX .EQ. 1) MASAP = MASAP*1E10/ACHE 

      !Now to masclet units
      MASAP = MASAP/UM
      U2PA = U2PA/UV
      U3PA = U3PA/UV
      U4PA = U4PA/UV

      RETURN

!***********************************************************************
END SUBROUTINE READ_AREPO_HDF5
!***********************************************************************



!***********************************************************************
SUBROUTINE READ_FLAMINGO_DMO_HDF5(ITER,FILES_PER_SNAP,MASSDM,ACHE,ZETA)
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

      !DARK MATTER
      PARTTYPEX = 2

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

       ! WRITE(*,*) 'Assigning masses ...'
       MASAP = MASSDM ! mass of DM particles in Msun
              
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
END SUBROUTINE READ_FLAMINGO_DMO_HDF5
!***********************************************************************
