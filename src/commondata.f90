
MODULE COMMONDATA
      IMPLICIT NONE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!HARD-CODDED PARAMETERS: 
      INTEGER, PARAMETER:: FLAG_READ=1 !1=read clus,cldm,clst files until IR=LEVMAX
      INTEGER, PARAMETER:: FLAG_NEXT_GRAD=1 !1= when ddens is above the threshold check the gradient also in the next cell
      INTEGER, PARAMETER:: FLAG_DIV_EDGE=1 !1=use condition on diverV to define the edges of voids
      INTEGER, PARAMETER:: NCELL_MIN=10 !min # of cells for keeping voids
!                                       ; for few cell voids inertia tensor is not reliable
      INTEGER, PARAMETER:: SIDE_MIN = 0.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!MASCLET units
      REAL*8 :: PI
      REAL*4, PARAMETER :: UL = 10.98          !to Mpc
      REAL*4, PARAMETER :: UV = 299792.458     !to km/s
      REAL*4, PARAMETER :: UM = 9.1717E18      !to in Msun

!variables (parameters) defined in the main unit
      !bounding box
       REAL*4 LADO0,LADO0PLUS

      !void finder parameters
       INTEGER FIRST,LAST,EVERY
       INTEGER NLEV
       INTEGER LEVMIN, LEVMAX
       REAL*4 DENS_THRE, DENS_THRE2, GRAD_THRE, DIV_THRE, RMIN, RMIN_SUB
       INTEGER NL2 !$!$ level of AMR hierarchy until which the files are read and arrays are allocated
       INTEGER FLAG_DIV !0=TOTAL, 1=DM, 2=GAS, -1=LINEAL RECONSTRUCTION
       INTEGER FLAG_DENS !0=TOTAL, 1=DM, 2=GAS
       INTEGER NCOX, NCOY, NCOZ
       INTEGER FLAG_VEL_INTERP !0=TSC, 1=SPH
       INTEGER FLAG_DENS_INTERP !0=TSC, 1=SPH
       INTEGER FLAG_WRITE_CUBES !0=NO, 1=YES
       INTEGER FLAG_WRITE_PIECES !0=NO, 1=YES
       INTEGER FLAG_DATA !0=MASCLET, 1=TEST1 ....
       INTEGER FLAG_PERIOD !0=NO, 1=YES
       INTEGER FLAG_VEL_AVAILABLE !0=NO, 1=YES
       INTEGER NHYX, NHYY, NHYZ, NPALEV, NLEVELS, NAMRX, NAMRY, NAMRZ
       INTEGER(KIND=8) PARTIRED
       INTEGER NVOID_MAX
       
!variables defined (and some of them allocated) in READ_MASCLET
       REAL*4, ALLOCATABLE:: U1G(:,:,:) !gas density contrast at ir=0
       REAL*4, ALLOCATABLE:: U2G(:,:,:) !gas velocity in x
       REAL*4, ALLOCATABLE:: U3G(:,:,:) !gas velocity in x
       REAL*4, ALLOCATABLE:: U4G(:,:,:) !gas velocity in x
       REAL*4, ALLOCATABLE:: U11G(:,:,:,:) !gas density contrast at ir>0
       REAL*4, ALLOCATABLE:: U12G(:,:,:,:) ! gas  x-velocity (eulerian) for refined levels
       REAL*4, ALLOCATABLE:: U13G(:,:,:,:) ! gas  y-velocity (eulerian) for refined levels
       REAL*4, ALLOCATABLE:: U14G(:,:,:,:) ! gas  z-velocity (eulerian) for refined levels

       REAL*4, ALLOCATABLE:: U1S(:,:,:) !stellar density contrast at ir=0
       REAL*4, ALLOCATABLE:: U11S(:,:,:,:)  !stellar density contrast at ir>0

       REAL*4, ALLOCATABLE:: U1DM(:,:,:) !DM density contrast at ir=0
       REAL*4, ALLOCATABLE:: U11DM(:,:,:,:)  !DM density contrast at ir>0

       REAL*4, ALLOCATABLE:: U2PA(:)  !DM particle velocity in x
       REAL*4, ALLOCATABLE:: U3PA(:)  !DM particle velocity in y
       REAL*4, ALLOCATABLE:: U4PA(:)  !DM particle velocity in z
       REAL*4, ALLOCATABLE:: MASAP(:) !DM particle mass
       REAL*4, ALLOCATABLE:: RXPA(:)  !DM particle position in x
       REAL*4, ALLOCATABLE:: RYPA(:)  !DM particle position in x
       REAL*4, ALLOCATABLE:: RZPA(:)  !DM particle position in x
       INTEGER(KIND=8) NPARTT   !$!$ Total number of DM particles

!MASCLET AMR GRID variables
       INTEGER, ALLOCATABLE, DIMENSION(:) :: NPATCH(:), PARE(:), NPART(:), NPARTST(:)
       INTEGER, ALLOCATABLE, DIMENSION(:) :: PATCHNX(:), PATCHNY(:), PATCHNZ(:)
       INTEGER, ALLOCATABLE, DIMENSION(:) :: PATCHX(:),  PATCHY(:),  PATCHZ(:)
       REAL*4, ALLOCATABLE, DIMENSION(:) :: PATCHRX(:), PATCHRY(:), PATCHRZ(:)

!defined/allocated in smooth
       REAL*4, ALLOCATABLE:: U1CO(:,:,:) !total density contrast, common grid
       REAL*4, ALLOCATABLE:: U1DMCO(:,:,:) !DM density contrast, common grid
       REAL*4, ALLOCATABLE:: U1GCO(:,:,:) !gas density contrast, common grid
       REAL*4, ALLOCATABLE:: U1SCO(:,:,:) !gas density contrast, common grid
       REAL*4, ALLOCATABLE:: U2DMCO(:,:,:) !DM x-velocity (eulerian), common grid
       REAL*4, ALLOCATABLE:: U3DMCO(:,:,:) !DM y-velocity (eulerian), common grid
       REAL*4, ALLOCATABLE:: U4DMCO(:,:,:) !DM z-velocity (eulerian), common grid
       REAL*4, ALLOCATABLE:: U2GCO(:,:,:) !gas x-velocity (eulerian), common grid
       REAL*4, ALLOCATABLE:: U3GCO(:,:,:) !gas y-velocity (eulerian), common grid
       REAL*4, ALLOCATABLE:: U4GCO(:,:,:) !gas z-velocity (eulerian), common grid
       REAL*4, ALLOCATABLE:: DIVERCO(:,:,:) !total divergence, considing gas and/or DM
       REAL*4, ALLOCATABLE:: DIVERDMCO(:,:,:) !DM divergence
       REAL*4, ALLOCATABLE:: DIVERGCO(:,:,:) !gas divergence

!variables defined in MALLA
       REAL*4 DX,DY,DZ
       REAL*4, ALLOCATABLE :: RADX(:), RADY(:),RADZ(:)
       REAL*4 DX0,DY0,DZ0
       REAL*4, ALLOCATABLE :: RADX0(:), RADY0(:),RADZ0(:)

!variables defined for subvoid management
       INTEGER, ALLOCATABLE:: FLAGV(:,:,:)
       INTEGER, ALLOCATABLE:: FLAG_SUB(:,:,:)
       INTEGER, ALLOCATABLE:: MARCAP(:,:,:) !inner (1:NX)
       INTEGER, ALLOCATABLE :: MARCAP2(:,:,:) !outer (low1:low2)

!variables defined in void_find
       INTEGER, ALLOCATABLE :: INICIOX(:),FINALX(:),INICIOY(:),FINALY(:), &
                                          INICIOZ(:), FINALZ(:), ICX(:), ICY(:), ICZ(:)
       INTEGER, ALLOCATABLE :: PERIODICAL(:)
       REAL*4, ALLOCATABLE :: VOL(:)
       INTEGER, ALLOCATABLE :: FATHER(:)
       REAL*4, ALLOCATABLE :: RINIXCO(:), RFINXCO(:), RINIYCO(:), RFINYCO(:), RINIZCO(:), RFINZCO(:)

!variables defined in overlapping
       INTEGER, ALLOCATABLE :: MARCA(:,:,:)

!variables defined in DIVER_FINA
       REAL*4, ALLOCATABLE :: DIVER0(:,:,:)
       REAL*4, ALLOCATABLE :: DIVER(:,:,:,:)

!variables defined in VMESH
       REAL*4, ALLOCATABLE:: UDMR(:,:,:),UGR(:,:,:), USR(:,:,:), DIVR(:,:,:)

END MODULE COMMONDATA