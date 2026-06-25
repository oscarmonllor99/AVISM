
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

MODULE PARTICLES
   IMPLICIT NONE
   PUBLIC

   !SPH kernel mode (1:M4, 2:Gaussian)
   INTEGER :: KERNEL_MODE = 1

! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

CONTAINS 

!************************************************************************
   SUBROUTINE DENS_INTERP_TSC(NX,NY,NZ,NPARTT,RXPA,RYPA,RZPA,MASAP,U)
!************************************************************************
      USE COMMONDATA, ONLY: DX,DY,DZ,RADX,RADY,RADZ,FLAG_PERIOD
      IMPLICIT NONE
      !input variables
      INTEGER NX,NY,NZ
      INTEGER(KIND=8) :: NPARTT, IP
      REAL*4, INTENT(IN) :: RXPA(:),RYPA(:),RZPA(:),MASAP(:)

      !local variables
      INTEGER I, J, K
      INTEGER IX, JY, KZ, I3, J3, K3
      REAL*4 RADX1, RADXNX, RADY1, RADYNY, RADZ1, RADZNZ, LENG
      REAL*4 VX(-1:1), VY(-1:1), VZ(-1:1)
      REAL*4 BAS, XTMP, YTMP, ZTMP
      REAL*4 INC_VAL

      !REAL*8 precision auxiliary U
      REAL*8, ALLOCATABLE :: U8(:,:,:)

      !output
      REAL*4 U(NX,NY,NZ)

      !smallest/largest coordinates of the coarse grid
      RADX1=RADX(1)-0.5*DX
      RADXNX=RADX(NX)+0.5*DX
      RADY1=RADY(1)-0.5*DY
      RADYNY=RADY(NY)+0.5*DY
      RADZ1=RADZ(1)-0.5*DZ
      RADZNZ=RADZ(NZ)+0.5*DZ
      LENG=RADXNX-RADX1

      ALLOCATE(U8(1:NX,1:NY,1:NZ))
      U8 = 0.0_8

   !$OMP PARALLEL DO SHARED(NPARTT,RXPA,RYPA,RZPA,RADX1,RADY1,RADZ1, &
   !$OMP     RADXNX,RADYNY,RADZNZ,LENG,DX,DY,DZ,NX,NY,NZ,MASAP,U8, &
   !$OMP     RADX,RADY,RADZ,FLAG_PERIOD), &
   !$OMP PRIVATE(IP,BAS,I,J,K,VX,VY,VZ,IX,JY,KZ,I3,J3,K3,XTMP,YTMP,ZTMP,INC_VAL), &
   !$OMP DEFAULT(NONE)
      DO IP=1,NPARTT
         ! Use local variables to avoid modifying intent(in) arrays
         XTMP = RXPA(IP)
         YTMP = RYPA(IP)
         ZTMP = RZPA(IP)

         !move by LENG particles outside the grid locally
         IF (FLAG_PERIOD .EQ. 1) THEN
            IF(XTMP.LT.RADX1)  XTMP=XTMP+LENG
            IF(XTMP.GE.RADXNX) XTMP=XTMP-LENG
            IF(YTMP.LT.RADY1)  YTMP=YTMP+LENG
            IF(YTMP.GE.RADYNY) YTMP=YTMP-LENG
            IF(ZTMP.LT.RADZ1)  ZTMP=ZTMP+LENG
            IF(ZTMP.GE.RADZNZ) ZTMP=ZTMP-LENG
         ENDIF

         !I,J,K: cell of the particle in the grid
         I=INT(((XTMP-RADX(1))/DX)+0.49999) + 1
         J=INT(((YTMP-RADY(1))/DY)+0.49999) + 1
         K=INT(((ZTMP-RADZ(1))/DZ)+0.49999) + 1

         BAS=ABS(RADX(I-1)-XTMP)
         VX(-1)=0.5*(1.5-BAS/DX)**2
         BAS=ABS(RADX(I)-XTMP)
         VX(0)=0.75-(BAS/DX)**2
         BAS=ABS(RADX(I+1)-XTMP)
         VX(1)=0.5*(1.5-BAS/DX)**2

         BAS=ABS(RADY(J-1)-YTMP)
         VY(-1)=0.5*(1.5-BAS/DY)**2
         BAS=ABS(RADY(J)-YTMP)
         VY(0)=0.75-(BAS/DY)**2
         BAS=ABS(RADY(J+1)-YTMP)
         VY(1)=0.5*(1.5-BAS/DY)**2

         BAS=ABS(RADZ(K-1)-ZTMP)
         VZ(-1)=0.5*(1.5-BAS/DZ)**2
         BAS=ABS(RADZ(K)-ZTMP)
         VZ(0)=0.75-(BAS/DZ)**2
         BAS=ABS(RADZ(K+1)-ZTMP)
         VZ(1)=0.5*(1.5-BAS/DZ)**2

         DO KZ=-1,1
         DO JY=-1,1
         DO IX=-1,1
            I3=I+IX
            J3=J+JY
            K3=K+KZ
            IF (FLAG_PERIOD .EQ. 1) THEN
               IF (I3.LT.1) I3=I3+NX
               IF (I3.GT.NX) I3=I3-NX
               IF (J3.LT.1) J3=J3+NY
               IF (J3.GT.NY) J3=J3-NY
               IF (K3.LT.1) K3=K3+NZ
               IF (K3.GT.NZ) K3=K3-NZ
            ENDIF

            INC_VAL = MASAP(IP)*VX(IX)*VY(JY)*VZ(KZ)
            !$OMP ATOMIC
            U8(I3,J3,K3) = U8(I3,J3,K3) + INC_VAL
         END DO
         END DO
         END DO
      ENDDO 

      U(1:NX,1:NY,1:NZ) = U8(1:NX,1:NY,1:NZ)
      DEALLOCATE(U8)

      !FOR THE MOMENT U IS MASS

!************************************************************************
   END SUBROUTINE DENS_INTERP_TSC
!************************************************************************


!********************************************************************************
   SUBROUTINE VEL_INTERP_TSC(NX,NY,NZ,NPARTT,RXPA,RYPA,RZPA, &
                     U2PA,U3PA,U4PA,U2,U3,U4)
!********************************************************************************
      USE COMMONDATA, ONLY: DX,DY,DZ,RADX,RADY,RADZ,FLAG_PERIOD

      IMPLICIT NONE
      !input variables
      INTEGER NX,NY,NZ
      INTEGER(KIND=8) :: NPARTT, LOW1, IP
      REAL*4, INTENT(IN) :: RXPA(:),RYPA(:),RZPA(:)
      REAL*4, INTENT(IN) :: U2PA(:),U3PA(:),U4PA(:)

      !local variables
      INTEGER I, J, K
      INTEGER IX, JY, KZ, I3, J3, K3
      REAL*4 RADX1, RADXNX, RADY1, RADYNY, RADZ1, RADZNZ, LENG
      REAL*4 VX(-1:1), VY(-1:1), VZ(-1:1)
      REAL*4 BAS, WEI, XTMP, YTMP, ZTMP
      REAL*4 INC_U2, INC_U3, INC_U4, INC_WT
      REAL*4, ALLOCATABLE :: WT(:,:,:)

      !output
      REAL*4 U2(NX,NY,NZ), U3(NX,NY,NZ), U4(NX,NY,NZ)

      ALLOCATE(WT(NX,NY,NZ))
      !$OMP PARALLEL DO SHARED(WT,NX,NY,NZ),PRIVATE(I,J,K),DEFAULT(NONE)
      DO K=1,NZ
      DO J=1,NY
      DO I=1,NX
         WT(I,J,K)=0.
      END DO
      END DO
      END DO

      RADX1=RADX(1)-0.5*DX
      RADXNX=RADX(NX)+0.5*DX
      RADY1=RADY(1)-0.5*DY
      RADYNY=RADY(NY)+0.5*DY
      RADZ1=RADZ(1)-0.5*DZ
      RADZNZ=RADZ(NZ)+0.5*DZ
      LENG=RADXNX-RADX1

      LOW1=NPARTT

   !$OMP PARALLEL DO SHARED(LOW1,RXPA,RYPA,RZPA,RADX1,RADY1,RADZ1, &
   !$OMP     RADXNX,RADYNY,RADZNZ,LENG,DX,DY,DZ,NX,NY,NZ, &
   !$OMP     U2PA,U3PA,U4PA,U2,U3,U4,WT,RADX,RADY,RADZ,FLAG_PERIOD), &
   !$OMP PRIVATE(IP,BAS,I,J,K,VX,VY,VZ,IX,JY,KZ,I3,J3,K3,XTMP,YTMP,ZTMP, &
   !$OMP INC_U2, INC_U3, INC_U4, INC_WT), DEFAULT(NONE)
      DO IP=1, LOW1
         XTMP = RXPA(IP)
         YTMP = RYPA(IP)
         ZTMP = RZPA(IP)

         IF(FLAG_PERIOD .EQ. 1) THEN
            IF(XTMP.LT.RADX1) XTMP=XTMP+LENG
            IF(XTMP.GE.RADXNX) XTMP=XTMP-LENG
            IF(YTMP.LT.RADY1) YTMP=YTMP+LENG
            IF(YTMP.GE.RADYNY) YTMP=YTMP-LENG
            IF(ZTMP.LT.RADZ1) ZTMP=ZTMP+LENG
            IF(ZTMP.GE.RADZNZ) ZTMP=ZTMP-LENG
         ENDIF

         I=INT(((XTMP-RADX(1))/DX)+0.49999) + 1
         J=INT(((YTMP-RADY(1))/DY)+0.49999) + 1
         K=INT(((ZTMP-RADZ(1))/DZ)+0.49999) + 1

         BAS=ABS(RADX(I-1)-XTMP)
         VX(-1)=0.5*(1.5-BAS/DX)**2
         BAS=ABS(RADX(I)-XTMP)
         VX(0)=0.75-(BAS/DX)**2
         BAS=ABS(RADX(I+1)-XTMP)
         VX(1)=0.5*(1.5-BAS/DX)**2

         BAS=ABS(RADY(J-1)-YTMP)
         VY(-1)=0.5*(1.5-BAS/DY)**2
         BAS=ABS(RADY(J)-YTMP)
         VY(0)=0.75-(BAS/DY)**2
         BAS=ABS(RADY(J+1)-YTMP)
         VY(1)=0.5*(1.5-BAS/DY)**2

         BAS=ABS(RADZ(K-1)-ZTMP)
         VZ(-1)=0.5*(1.5-BAS/DZ)**2
         BAS=ABS(RADZ(K)-ZTMP)
         VZ(0)=0.75-(BAS/DZ)**2
         BAS=ABS(RADZ(K+1)-ZTMP)
         VZ(1)=0.5*(1.5-BAS/DZ)**2

         DO KZ=-1,1
         DO JY=-1,1
         DO IX=-1,1
            I3=I+IX
            J3=J+JY
            K3=K+KZ
            IF(FLAG_PERIOD .EQ. 1) THEN
               IF (I3.LT.1) I3=I3+NX
               IF (I3.GT.NX) I3=I3-NX
               IF (J3.LT.1) J3=J3+NY
               IF (J3.GT.NY) J3=J3-NY
               IF (K3.LT.1) K3=K3+NZ
               IF (K3.GT.NZ) K3=K3-NZ
            ENDIF
            
            INC_WT = VX(IX)*VY(JY)*VZ(KZ)
            INC_U2 = U2PA(IP)*INC_WT
            INC_U3 = U3PA(IP)*INC_WT
            INC_U4 = U4PA(IP)*INC_WT

            !atomic for avoiding data racing
            !$OMP ATOMIC
            U2(I3,J3,K3) = U2(I3,J3,K3) + INC_U2
            !$OMP ATOMIC
            U3(I3,J3,K3) = U3(I3,J3,K3) + INC_U3
            !$OMP ATOMIC
            U4(I3,J3,K3) = U4(I3,J3,K3) + INC_U4
            !$OMP ATOMIC
            WT(I3,J3,K3) = WT(I3,J3,K3) + INC_WT
         END DO
         END DO
         END DO
      ENDDO 

      WRITE(*,*) '     Fraction of cells without particles (TSC):', REAL(COUNT(WT .EQ. 0.))/REAL(NX*NY*NZ)

   !$OMP PARALLEL DO SHARED(U2,U3,U4,WT,NX,NY,NZ), PRIVATE(I,J,K,WEI), DEFAULT(NONE)
      DO K=1,NZ
         DO J=1,NY
            DO I=1,NX
               WEI=WT(I,J,K)
               IF(WEI==0.0) WEI=1.0
               U2(I,J,K)=U2(I,J,K)/WEI
               U3(I,J,K)=U3(I,J,K)/WEI
               U4(I,J,K)=U4(I,J,K)/WEI
            END DO
         END DO
      END DO

      DEALLOCATE(WT)
!************************************************************************
   END SUBROUTINE VEL_INTERP_TSC
!************************************************************************


!********************************************************************************  
   SUBROUTINE H_DISTANCE(KNEIGHBOURS,NX,NY,NZ,NPARTT,TREE,SMASK,HPART)
   !this subroutine calculates the kernel distance h(x) for every cell x
   !required by DDENS_INTERP_SPH if VVEL_INTERP_SPH is not used
!********************************************************************************  
         USE COSMOKDTREE
         USE COMMONDATA, ONLY:DX,DY,DZ,RADX,RADY,RADZ

         IMPLICIT NONE
         !input variables
         INTEGER KNEIGHBOURS
         INTEGER NX,NY,NZ
         INTEGER(KIND=8) :: NPARTT, CONTA, I
         INTEGER*1 :: SMASK(:,:,:)

         !LOCAL
         INTEGER :: IX,JY,KZ

         !h smoothing length for each particle
         REAL,ALLOCATABLE::DIST(:)
         INTEGER(KIND=8),ALLOCATABLE::NEIGH(:)
         
         !query
         type(KDTreeNode), pointer :: TREE
         type(KDTreeResult) :: QUERY
         REAL*4 :: TAR(3)
         !SPH related
         REAL*4 :: HPART(NPARTT)

         !$OMP PARALLEL DO SHARED(KNEIGHBOURS,NX,NY,NZ,SMASK, &
         !$OMP                   TREE,DX,DY,DZ,RADX,RADY,RADZ,HPART,KERNEL_MODE) &
         !$OMP PRIVATE(I,IX,JY,KZ,QUERY,TAR,DIST,NEIGH,CONTA), &
         !$OMP SCHEDULE(DYNAMIC), DEFAULT(NONE)
         DO KZ=1,NZ
            DO JY=1,NY
               DO IX=1,NX

                  IF (SMASK(IX,JY,KZ) .EQ. 0) CYCLE
                  TAR(1) = RADX(IX)
                  TAR(2) = RADY(JY)
                  TAR(3) = RADZ(KZ)

                  !Search cell's kneighbours
                  ALLOCATE(DIST(KNEIGHBOURS), NEIGH(KNEIGHBOURS))
                  QUERY = knn_search(TREE, TAR, KNEIGHBOURS)
                  DIST = QUERY%dist
                  NEIGH = QUERY%idx
         
                  IF (DIST(KNEIGHBOURS).GT.DX) THEN
                     CONTA=KNEIGHBOURS 
                  ELSE 
                     DEALLOCATE(DIST,NEIGH)
                     !.true. stands for sorting
                     QUERY = ball_search(TREE, TAR, DX, .true.)
                     CONTA = SIZE(QUERY%dist)
                     ALLOCATE(DIST(CONTA), NEIGH(CONTA))
                     DIST = QUERY%dist
                     NEIGH = QUERY%idx
                  END IF

                  !!!! For SPH density interpolation !!!!!!!!!!!!
                  !CHANGE VALUE OF HPART FOR ALL CELL NEIGHBOURS
                  DO I=1,CONTA
                     !$OMP ATOMIC
                     HPART(NEIGH(I)) = MAX(HPART(NEIGH(I)), DIST(I))
                  END DO
                  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                  
                  DEALLOCATE(DIST,NEIGH)
               ENDDO
            ENDDO
         ENDDO
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!********************************************************************************
   END SUBROUTINE H_DISTANCE
!********************************************************************************


!********************************************************************************  
   SUBROUTINE H_RANDOMS(NRAND,NPARTT,RANDX,RANDY,RANDZ,TREE,HPART,HRAND,KMEAN)
   !ESTIMATES THE SMOOTHING LENGTH OF RANDOMS ACCORDIN TO NEIGHBOUR GALAXIES/PARTICLES
!********************************************************************************  
      USE COSMOKDTREE

      IMPLICIT NONE
      !input variables
      INTEGER(KIND=8) :: NRAND,NPARTT, I
      REAL*4, INTENT(IN) :: RANDX(:),RANDY(:),RANDZ(:)
      INTEGER :: KMEAN, J

      !query
      type(KDTreeNode), pointer :: TREE
      type(KDTreeResult) :: QUERY
      REAL*4 :: TAR(3)

      !SPH related
      REAL*4 :: HPART(NPARTT)
      REAL*4 :: HRAND(NRAND)
      REAL*4 :: HKERN

      !$OMP PARALLEL DO SHARED(NRAND,RANDX,RANDY,RANDZ,TREE,HRAND,HPART,KMEAN),&
      !$OMP PRIVATE(I,J,QUERY,TAR,HKERN), DEFAULT(NONE)
      DO I=1,NRAND
         TAR(1) = RANDX(I)
         TAR(2) = RANDY(I)
         TAR(3) = RANDZ(I)
         QUERY = knn_search(TREE, TAR, KMEAN)
         HKERN = 0.
         DO J=1,KMEAN
            HKERN = HKERN + HPART(QUERY%idx(J))
         ENDDO
         HRAND(I) = HKERN / REAL(KMEAN)
      ENDDO
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!********************************************************************************
   END SUBROUTINE H_RANDOMS
!********************************************************************************



!********************************************************************************  
   SUBROUTINE FILL_H_ZEROS(KNEIGHBOURS,NPARTT,RXPA,RYPA,RZPA,TREE,HPART)
   !fills zeros of hpart (particles that didnt contribute to the mean SPH kernel on velocity)
   !with distance to k-th neighbour
!********************************************************************************  
      USE COSMOKDTREE

      IMPLICIT NONE
      !input variables
      INTEGER(KIND=8) :: NPARTT, I
      INTEGER :: KNEIGHBOURS
      REAL*4, INTENT(IN) :: RXPA(:),RYPA(:),RZPA(:)

      !query
      type(KDTreeNode), pointer :: TREE
      type(KDTreeResult) :: QUERY
      REAL*4 :: TAR(3)

      !SPH related
      REAL*4 :: HPART(NPARTT)

      !$OMP PARALLEL DO SHARED(RXPA,RYPA,RZPA,NPARTT,TREE,HPART,KNEIGHBOURS),&
      !$OMP PRIVATE(I,QUERY,TAR), DEFAULT(NONE)
      DO I=1,NPARTT
         IF (HPART(I) > 0.) CYCLE
         TAR(1) = RXPA(I)
         TAR(2) = RYPA(I)
         TAR(3) = RZPA(I)
         QUERY = knn_search(TREE, TAR, KNEIGHBOURS)
         HPART(I) = QUERY%dist(KNEIGHBOURS)
      ENDDO
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!********************************************************************************
   END SUBROUTINE FILL_H_ZEROS
!********************************************************************************

!********************************************************************************
   SUBROUTINE DDENS_INTERP_SPH(NX,NY,NZ,NPARTT,RXPA,RYPA,RZPA, &
                  HPART,MASAP,SMASK,U)
!using HPART, smoothing lenght given
!by the furthest cell to which the particle contributes in the
!velocity interpolation
!********************************************************************************
         USE COMMONDATA, ONLY: DX,DY,DZ,RADX,RADY,RADZ

         IMPLICIT NONE
         !input variables
         INTEGER NX,NY,NZ
         INTEGER(KIND=8) :: NPARTT, I
         REAL*4, INTENT(IN) :: RXPA(:),RYPA(:),RZPA(:)
         REAL*4, INTENT(IN) :: MASAP(:)
         INTEGER*1 :: SMASK(:,:,:)

         !LOCAL
         INTEGER :: IX,JY,KZ,IX2,JY2,KZ2
         REAL*4 :: DX2, DY2, DZ2

         !grid position for each particle
         REAL :: RX1,RY1,RZ1

         !INTERPOLATION VARIABLES
         REAL*4 :: H,DDIST,D2,NORM
         INTEGER :: INCREX,INCREY,INCREZ
         
         !SPH related
         REAL*4 :: Q,Q2, INV_H, INV_H_SQ, H2_MAX, WEIGHT
         REAL*4 :: HPART(NPARTT)

         !REAL*8 precision auxiliary U
         REAL*8, ALLOCATABLE :: U8(:,:,:)

         !output
         REAL*4 U(NX,NY,NZ)

         !grid first cell center
         RX1 = RADX(1)
         RY1 = RADY(1)
         RZ1 = RADZ(1)   

         ALLOCATE(U8(1:NX,1:NY,1:NZ))
         U8 = 0.0_8

         !$OMP PARALLEL DO SHARED(NPARTT,HPART,RX1,RY1,RZ1,SMASK, &
         !$OMP                   RXPA,RYPA,RZPA,MASAP,RADX,RADY,RADZ, &
         !$OMP                   DX,DY,DZ,NX,NY,NZ,KERNEL_MODE,U8), &
         !$OMP PRIVATE(I,IX,JY,KZ,IX2,JY2,KZ2,H,DX2,DY2,DZ2, &
         !$OMP       INCREX,INCREY,INCREZ,DDIST,NORM,Q,Q2, &
         !$OMP       D2,H2_MAX,INV_H,INV_H_SQ,WEIGHT), &
         !$OMP SCHEDULE(DYNAMIC), DEFAULT(NONE)
         !$$$$$$$$$$$$$$$$$$$$
         DO I=1,NPARTT
         !$$$$$$$$$$$$$$$$$$$$
            IX = INT((RXPA(I)-RX1)/DX)+1
            JY = INT((RYPA(I)-RY1)/DY)+1
            KZ = INT((RZPA(I)-RZ1)/DZ)+1

            H = HPART(I)
            INV_H = 2./H
            INV_H_SQ = INV_H * INV_H
            H2_MAX = 4.*H*H 

            INCREX = INT(2*H/DX + 0.5) !0.5 accounts for the fact that
            INCREY = INT(2*H/DY + 0.5) !we are using the center of the cell
            INCREZ = INT(2*H/DZ + 0.5) !as the reference point

            !-------------------------------------
            !First, find normalisation factor
            NORM = 0.
            DO KZ2=KZ-INCREZ,KZ+INCREZ
               IF (KZ2.LT.1 .OR. KZ2.GT.NZ) CYCLE
               DZ2 = (RZPA(I)-RADZ(KZ2))**2
               DO JY2=JY-INCREY,JY+INCREY
                  IF (JY2.LT.1 .OR. JY2.GT.NY) CYCLE
                  DY2 = (RYPA(I)-RADY(JY2))**2
                  DO IX2=IX-INCREX,IX+INCREX
                     IF (IX2.LT.1 .OR. IX2.GT.NX) CYCLE
                     ! IF (SMASK(IX2,JY2,KZ2) .EQ. 0) CYCLE
                     DX2 = (RXPA(I)-RADX(IX2))**2

                     !DISTANCE BETWEEN PARTICLE AND GRID POINT
                     D2 = DX2+DY2+DZ2

                     !SKIP CUBE CORNERS
                     IF (D2 .GT. H2_MAX) CYCLE

                     !---------------------------------------------
                     ! Cubic spline kernel (M4)
                     IF (KERNEL_MODE .EQ. 1) THEN
                        DDIST = SQRT(D2)
                        Q = DDIST * INV_H
                        IF (Q.LE.1.0) THEN
                           WEIGHT=1.0-1.5*Q**2*(1.-0.5*Q)
                        ELSE IF (Q.LE.2.0) THEN
                           WEIGHT=0.25*(2.-Q)**3
                        ELSE
                           WEIGHT=0.0
                        END IF
                     ! Gaussian kernel
                     ELSE IF(KERNEL_MODE .EQ. 2) THEN
                        Q2 = D2 * INV_H_SQ
                        WEIGHT = EXP(-Q2)
                     !....
                     ENDIF
                     !---------------------------------------------
                     
                     NORM = NORM + WEIGHT

                  ENDDO
               ENDDO
            ENDDO

            ! Particle only contributes to its own cell
            IF (NORM .EQ. 0.) THEN
               IF(IX.GE.1 .AND. IX.LE.NX .AND. &
                  JY.GE.1 .AND. JY.LE.NY .AND. &
                  KZ.GE.1 .AND. KZ.LE.NZ) THEN
                  !$OMP ATOMIC
                  U8(IX,JY,KZ) =  U8(IX,JY,KZ) + REAL(MASAP(I),8)
               ENDIF
               CYCLE
            ENDIF

            !-------------------------------------
            !Now, interpolate
            DO KZ2=KZ-INCREZ,KZ+INCREZ
               IF (KZ2.LT.1 .OR. KZ2.GT.NZ) CYCLE
               DZ2 = (RZPA(I)-RADZ(KZ2))**2
               DO JY2=JY-INCREY,JY+INCREY
                  IF (JY2.LT.1 .OR. JY2.GT.NY) CYCLE
                  DY2 = (RYPA(I)-RADY(JY2))**2
                  DO IX2=IX-INCREX,IX+INCREX
                     IF (IX2.LT.1 .OR. IX2.GT.NX) CYCLE
                     ! IF (SMASK(IX2,JY2,KZ2) .EQ. 0) CYCLE
                     DX2 = (RXPA(I)-RADX(IX2))**2

                     !DISTANCE BETWEEN PARTICLE AND GRID POINT
                     D2 = DX2+DY2+DZ2

                     !SKIP CUBE CORNERS
                     IF (D2 .GT. H2_MAX) CYCLE

                     !---------------------------------------------
                     ! Cubic spline kernel (M4)
                     IF (KERNEL_MODE .EQ. 1) THEN
                        DDIST = SQRT(D2)
                        Q = DDIST * INV_H
                        IF (Q.LE.1.0) THEN
                           WEIGHT=1.0-1.5*Q**2*(1.-0.5*Q)
                        ELSE IF (Q.LE.2.0) THEN
                           WEIGHT=0.25*(2.-Q)**3
                        ELSE
                           WEIGHT=0.0
                        END IF
                     ! Gaussian kernel
                     ELSE IF(KERNEL_MODE .EQ. 2) THEN
                        Q2 = D2 * INV_H_SQ
                        WEIGHT = EXP(-Q2)
                     !....
                     ENDIF
                     !---------------------------------------------
                           
                     !Update
                     !$OMP ATOMIC
                     U8(IX2,JY2,KZ2) = U8(IX2,JY2,KZ2) + REAL(WEIGHT*MASAP(I)/NORM,8)

                  ENDDO
               ENDDO
            ENDDO
            !-------------------------------------

         !$$$$$$$$$$$$$$$$$$$$
         END DO
         !$$$$$$$$$$$$$$$$$$$$

         U(1:NX,1:NY,1:NZ) = U8(1:NX,1:NY,1:NZ)
         DEALLOCATE(U8)

         !FOR THE MOMENT U IS MASS

!************************************************************************
   END SUBROUTINE DDENS_INTERP_SPH
!************************************************************************


!********************************************************************************
   SUBROUTINE VVEL_INTERP_SPH(KNEIGHBOURS,NX,NY,NZ,NPARTT,TREE, &
                        HPART,MASAP,U2PA,U3PA,U4PA,SMASK,U2,U3,U4)
!********************************************************************************
         USE COSMOKDTREE 
         USE COMMONDATA, ONLY: DX,DY,DZ,RADX,RADY,RADZ

         IMPLICIT NONE
         !input variables
         INTEGER KNEIGHBOURS, NX, NY, NZ
         INTEGER(KIND=8) :: NPARTT, CONTA, I
         REAL*4, INTENT(IN) :: U2PA(:),U3PA(:),U4PA(:)
         REAL*4, INTENT(IN) :: MASAP(:)
         INTEGER*1 :: SMASK(:,:,:)

         !LOCAL
         INTEGER :: IX,JY,KZ
         REAL*8 BAS8,BAS8X,BAS8Y,BAS8Z

         !h smoothing length for each particle
         REAL*4 HKERN
         REAL,ALLOCATABLE::DIST(:)
         REAL*8,ALLOCATABLE::DIST8(:)
         INTEGER(KIND=8),ALLOCATABLE::NEIGH(:)
         
         !query
         type(KDTreeNode), pointer, intent(in) :: TREE
         type(KDTreeResult) :: QUERY
         REAL*4 :: TAR(3)
         !SPH related
         REAL*4 :: HPART(NPARTT)
         
         !output
         REAL*4 U2(NX,NY,NZ), U3(NX,NY,NZ), U4(NX,NY,NZ)

         !$OMP PARALLEL DO SHARED(KNEIGHBOURS,NX,NY,NZ,MASAP,U2PA,U3PA,U4PA,U2,U3,U4, &
         !$OMP                   TREE,DX,DY,DZ,RADX,RADY,RADZ,SMASK,HPART,KERNEL_MODE) &
         !$OMP PRIVATE(I,IX,JY,KZ,QUERY,TAR,DIST,DIST8,NEIGH, & 
         !$OMP               HKERN,BAS8,BAS8X,BAS8Y,BAS8Z,CONTA), &
         !$OMP SCHEDULE(DYNAMIC), DEFAULT(NONE)
         DO KZ=1,NZ
            DO JY=1,NY
               DO IX=1,NX

                  !avoid out-of-mask cells
                  IF (SMASK(IX,JY,KZ) .EQ. 0) CYCLE

                  TAR(1) = RADX(IX)
                  TAR(2) = RADY(JY)
                  TAR(3) = RADZ(KZ)

                  !Search cell's kneighbours
                  ALLOCATE(DIST(KNEIGHBOURS), NEIGH(KNEIGHBOURS))
                  QUERY = knn_search(TREE, TAR, KNEIGHBOURS)
                  DIST = QUERY%dist
                  NEIGH = QUERY%idx

                  IF (DIST(KNEIGHBOURS).GT.DX) THEN
                     CONTA=KNEIGHBOURS 
                  ELSE 
                     DEALLOCATE(DIST,NEIGH)
                     !.true. stands for sorting
                     QUERY = ball_search(TREE, TAR, DX, .true.)
                     CONTA = SIZE(QUERY%dist)
                     ALLOCATE(DIST(CONTA), NEIGH(CONTA))
                     DIST = QUERY%dist
                     NEIGH = QUERY%idx
                  END IF

                  !!!! For SPH density interpolation !!!!!!!!!!!!
                  !CHANGE VALUE OF HPART FOR ALL CELL NEIGHBOURS
                  DO I=1,CONTA
                     !$OMP ATOMIC
                     HPART(NEIGH(I)) = MAX(HPART(NEIGH(I)), DIST(I))
                  END DO
                  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                  !cell's smoothing length
                  HKERN=DIST(CONTA)
                  !KERNEL FUNCTION
                  DO I=1,CONTA
                     CALL KERNEL_FUNC(HKERN,DIST(I))
                  END DO

                  !mass weighting
                  DO I=1,CONTA
                     DIST(I)=DIST(I)*MASAP(NEIGH(I))
                  END DO

                  ALLOCATE(DIST8(CONTA))
                  DIST8 = REAL(DIST, KIND=8)

                  !averaging
                  BAS8=0.D0
                  BAS8X=0.D0
                  BAS8Y=0.D0
                  BAS8Z=0.D0
                  DO I=1,CONTA 
                     BAS8=BAS8+DIST8(I)
                     BAS8X=BAS8X+DIST8(I)*REAL(U2PA(NEIGH(I)),8)
                     BAS8Y=BAS8Y+DIST8(I)*REAL(U3PA(NEIGH(I)),8)
                     BAS8Z=BAS8Z+DIST8(I)*REAL(U4PA(NEIGH(I)),8)
                  END DO
                  U2(IX,JY,KZ)=BAS8X/BAS8
                  U3(IX,JY,KZ)=BAS8Y/BAS8
                  U4(IX,JY,KZ)=BAS8Z/BAS8
                  
                  DEALLOCATE(DIST,NEIGH)
                  DEALLOCATE(DIST8)

               ENDDO
            ENDDO
         ENDDO
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!************************************************************************
   END SUBROUTINE VVEL_INTERP_SPH
!************************************************************************


!********************************************************************************
   SUBROUTINE PPART_DENS(KNEIGHBOURS,NPARTT,TREE,MASAP,PART_DENS)
! Obtains the density at the position of each particle using Nearest Neighbours
!********************************************************************************
         USE COSMOKDTREE 
         USE COMMONDATA, ONLY: RXPA,RYPA,RZPA,PI

         IMPLICIT NONE
         !input variables
         INTEGER :: KNEIGHBOURS
         INTEGER(KIND=8) :: NPARTT, CONTA, I, J
         REAL*4, INTENT(IN) :: MASAP(:)
         !LOCAL
         REAL*8 BAS8

         !h smoothing length for each particle
         REAL*4 HKERN
         REAL*4,ALLOCATABLE::DIST(:)
         INTEGER(KIND=8),ALLOCATABLE::NEIGH(:)

         !query
         type(KDTreeNode), pointer, intent(in) :: TREE
         type(KDTreeResult) :: QUERY
         
         REAL*4 :: TAR(3)

         !output
         REAL PART_DENS(NPARTT)


         !$OMP PARALLEL DO SHARED(KNEIGHBOURS,MASAP,TREE,RXPA,RYPA,RZPA,PI, &
         !$OMP                     PART_DENS) &
         !$OMP PRIVATE(I,J,QUERY,TAR,DIST,NEIGH,BAS8,HKERN,CONTA), &
         !$OMP SCHEDULE(DYNAMIC), DEFAULT(NONE)
         DO I=1,NPARTT
            TAR(1) = RXPA(I)
            TAR(2) = RYPA(I)
            TAR(3) = RZPA(I)

            !Search cell's kneighbours
            ALLOCATE(DIST(KNEIGHBOURS), NEIGH(KNEIGHBOURS))
            QUERY = knn_search(TREE, TAR, KNEIGHBOURS)
            NEIGH = QUERY%idx
            DIST = QUERY%dist

            ! Kernel support: distance to the furthest neighbour
            HKERN = DIST(KNEIGHBOURS)
            CONTA = KNEIGHBOURS

            BAS8=0.D0
            DO J=1,CONTA
               BAS8=BAS8+MASAP(NEIGH(J))
            END DO
            PART_DENS(I) = BAS8 / ((4.D0/3.D0)*PI*HKERN**3)

            DEALLOCATE(DIST,NEIGH)
         ENDDO
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!************************************************************************
   END SUBROUTINE PPART_DENS
!************************************************************************


!********************************************************************************
   SUBROUTINE VVEL_INTERP_SPH_VW(KNEIGHBOURS,NX,NY,NZ,NPARTT,TREE, &
                        HPART,PART_DENS,MASAP,U2PA,U3PA,U4PA,SMASK,U2,U3,U4)
!
! SAME AS VVEL_INTERP_SPH BUT USING VOLUME WEIGHTING instead of MASS
!********************************************************************************
         USE COSMOKDTREE 
         USE COMMONDATA, ONLY: DX,DY,DZ,RADX,RADY,RADZ, &   
                              RXPA,RYPA,RZPA,PI

         IMPLICIT NONE
         !input variables
         INTEGER KNEIGHBOURS
         INTEGER NX,NY,NZ
         INTEGER(KIND=8) :: NPARTT, CONTA, I, J
         REAL*4, INTENT(IN) :: U2PA(:),U3PA(:),U4PA(:)
         REAL*4, INTENT(IN) :: MASAP(:)
         INTEGER*1 :: SMASK(:,:,:)

         !LOCAL
         INTEGER :: IX,JY,KZ
         REAL*8 BAS8,BAS88,BAS8X,BAS8Y,BAS8Z

         !h smoothing length for each particle
         REAL*4 HKERN
         REAL*4,ALLOCATABLE::DIST(:)
         REAL*8,ALLOCATABLE::DIST8(:)
         INTEGER(KIND=8),ALLOCATABLE::NEIGH(:)
         
         !query
         type(KDTreeNode), pointer, intent(in) :: TREE
         type(KDTreeResult) :: QUERY
         
         REAL*4 :: TAR(3)
         !SPH related
         REAL*4 :: HPART(NPARTT)
         REAL*4 :: PART_DENS(NPARTT)
         
         !output
         REAL*4 U2(NX,NY,NZ), U3(NX,NY,NZ), U4(NX,NY,NZ)

         !$OMP PARALLEL DO SHARED(KNEIGHBOURS,NX,NY,NZ,MASAP,U2PA,U3PA,U4PA,U2,U3,U4, &
         !$OMP                   TREE,DX,DY,DZ,RADX,RADY,RADZ,PART_DENS,SMASK,HPART,KERNEL_MODE) &
         !$OMP PRIVATE(I,IX,JY,KZ,QUERY,TAR,DIST,DIST8,NEIGH, &
         !$OMP               HKERN,BAS8,BAS8X,BAS8Y,BAS8Z,CONTA), &
         !$OMP SCHEDULE(DYNAMIC), DEFAULT(NONE)
         DO KZ=1,NZ
            DO JY=1,NY
               DO IX=1,NX

                  !avoid out-of-mask cells
                  IF (SMASK(IX,JY,KZ) .EQ. 0) CYCLE

                  TAR(1) = RADX(IX)
                  TAR(2) = RADY(JY)
                  TAR(3) = RADZ(KZ)

                  !Search cell's kneighbours
                  ALLOCATE(DIST(KNEIGHBOURS), NEIGH(KNEIGHBOURS))
                  QUERY = knn_search(TREE, TAR, KNEIGHBOURS)
                  DIST = QUERY%dist
                  NEIGH = QUERY%idx

                  IF (DIST(KNEIGHBOURS).GT.DX) THEN
                     CONTA=KNEIGHBOURS 
                  ELSE 
                     DEALLOCATE(DIST,NEIGH)
                     !.true. stands for sorting
                     QUERY = ball_search(TREE, TAR, DX, .true.)
                     CONTA = SIZE(QUERY%dist)
                     ALLOCATE(DIST(CONTA), NEIGH(CONTA))
                     DIST = QUERY%dist
                     NEIGH = QUERY%idx
                  END IF

                  !!!! For SPH density interpolation !!!!!!!!!!!!
                  !CHANGE VALUE OF HPART FOR ALL CELL NEIGHBOURS
                  DO I=1,CONTA
                     !$OMP ATOMIC
                     HPART(NEIGH(I)) = MAX(HPART(NEIGH(I)), DIST(I))
                  END DO
                  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                  !cell's smoothing length
                  HKERN=DIST(CONTA)

                  !KERNEL FUNCTION
                  DO I=1,CONTA
                     CALL KERNEL_FUNC(HKERN,DIST(I))
                  END DO

                  !VOLUME-weighting
                  DO I=1,CONTA
                     DIST(I)=DIST(I) * (MASAP(NEIGH(I)) / &
                                        PART_DENS(NEIGH(I)) )
                  END DO

                  ALLOCATE(DIST8(CONTA))
                  DIST8 = REAL(DIST, KIND=8)

                  !averaging
                  BAS8=0.D0
                  BAS8X=0.D0
                  BAS8Y=0.D0
                  BAS8Z=0.D0
                  DO I=1,CONTA 
                     BAS8=BAS8+DIST8(I)
                     BAS8X=BAS8X+DIST8(I)*U2PA(NEIGH(I))
                     BAS8Y=BAS8Y+DIST8(I)*U3PA(NEIGH(I))
                     BAS8Z=BAS8Z+DIST8(I)*U4PA(NEIGH(I))
                  END DO

                  U2(IX,JY,KZ)=BAS8X/BAS8
                  U3(IX,JY,KZ)=BAS8Y/BAS8
                  U4(IX,JY,KZ)=BAS8Z/BAS8
                  
                  DEALLOCATE(DIST,NEIGH)
                  DEALLOCATE(DIST8)
               ENDDO
            ENDDO
         ENDDO
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!************************************************************************
   END SUBROUTINE VVEL_INTERP_SPH_VW
!************************************************************************



!************************************************************************
   SUBROUTINE KERNEL_FUNC(H,DIST)
!************************************************************************
!* DIST contains initially the distance (particle to cell)
!* and it is updated with the (unnormalised) value of the kernel
   IMPLICIT NONE
   REAL H,DIST,Q
      
   Q=DIST/(H/2) !q

   ! Cubic spline kernel (M4)
   IF (KERNEL_MODE .EQ. 1) THEN
      IF (Q.LE.1.0) THEN
         DIST=1.0-1.5*Q**2*(1.-0.5*Q)
      ELSE IF (Q.LE.2.0) THEN
         DIST=0.25*(2.-Q)**3
      ELSE
         DIST=0.0
      END IF
      
   ! Gaussian kernel
   ELSE IF(KERNEL_MODE .EQ. 2) THEN
      DIST = EXP(-Q**2)

   !....
   ENDIF
      
   RETURN
!************************************************************************
   END SUBROUTINE KERNEL_FUNC
!************************************************************************


! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

END MODULE PARTICLES

! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~