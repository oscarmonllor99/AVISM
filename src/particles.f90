
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

MODULE PARTICLES
   IMPLICIT NONE
   PUBLIC

! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

CONTAINS 

!************************************************************************
   SUBROUTINE DENS_INTERP_TSC(NX,NY,NZ,NPARTT,RXPA,RYPA,RZPA,MASAP,U)
!************************************************************************
      USE COMMONDATA, ONLY: DX,DY,DZ,RADX,RADY,RADZ,PARTIRED
      IMPLICIT NONE
      !input variables
      INTEGER NX,NY,NZ
      INTEGER(KIND=8) :: NPARTT, LOW1, IP
      REAL*4 RXPA(PARTIRED),RYPA(PARTIRED),RZPA(PARTIRED),MASAP(PARTIRED)

      !local variables
      INTEGER I, J, K
      INTEGER IX, JY, KZ, I3, J3, K3
      REAL*4 RADX1, RADXNX, RADY1, RADYNY, RADZ1, RADZNZ, LENG
      REAL*4 VX(-1:1), VY(-1:1), VZ(-1:1)
      REAL*4 BAS

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

      DO I=-1,1
         VX(I)=0.0
         VY(I)=0.0
         VZ(I)=0.0
      END DO

      LOW1=NPARTT

   !$OMP PARALLEL DO SHARED(LOW1,RXPA,RYPA,RZPA,RADX1,RADY1,RADZ1, &
   !$OMP     RADXNX,RADYNY,RADZNZ,LENG,DX,DY,DZ,NX,NY,NZ,MASAP,U, &
   !$OMP     RADX,RADY,RADZ), &
   !$OMP PRIVATE(IP,BAS,I,J,K,VX,VY,VZ,IX,JY,KZ,I3,J3,K3), &
   !$OMP DEFAULT(NONE)
      DO IP=1, LOW1

         !move by LENG particles outside the grid
         IF(RXPA(IP).LT.RADX1) RXPA(IP)=RXPA(IP)+LENG
         IF(RXPA(IP).GE.RADXNX) RXPA(IP)=RXPA(IP)-LENG
         IF(RYPA(IP).LT.RADY1) RYPA(IP)=RYPA(IP)+LENG
         IF(RYPA(IP).GE.RADYNY) RYPA(IP)=RYPA(IP)-LENG
         IF(RZPA(IP).LT.RADZ1) RZPA(IP)=RZPA(IP)+LENG
         IF(RZPA(IP).GE.RADZNZ) RZPA(IP)=RZPA(IP)-LENG

         !initial position
         !I,J,K: cell of the particle in the grid
         BAS=RXPA(IP)
         I=INT(((BAS-RADX(1))/DX)+0.49999) + 1
         BAS=RYPA(IP)
         J=INT(((BAS-RADY(1))/DY)+0.49999) + 1
         BAS=RZPA(IP)
         K=INT(((BAS-RADZ(1))/DZ)+0.49999) + 1

   !***  TSC

         BAS=ABS(RADX(I-1)-RXPA(IP))
         VX(-1)=0.5*(1.5-BAS/DX)**2
         BAS=ABS(RADX(I)-RXPA(IP))
         VX(0)=0.75-(BAS/DX)**2
         BAS=ABS(RADX(I+1)-RXPA(IP))
         VX(1)=0.5*(1.5-BAS/DX)**2

         BAS=ABS(RADY(J-1)-RYPA(IP))
         VY(-1)=0.5*(1.5-BAS/DY)**2
         BAS=ABS(RADY(J)-RYPA(IP))
         VY(0)=0.75-(BAS/DY)**2
         BAS=ABS(RADY(J+1)-RYPA(IP))
         VY(1)=0.5*(1.5-BAS/DY)**2

         BAS=ABS(RADZ(K-1)-RZPA(IP))
         VZ(-1)=0.5*(1.5-BAS/DZ)**2
         BAS=ABS(RADZ(K)-RZPA(IP))
         VZ(0)=0.75-(BAS/DZ)**2
         BAS=ABS(RADZ(K+1)-RZPA(IP))
         VZ(1)=0.5*(1.5-BAS/DZ)**2


         DO KZ=-1,1
         DO JY=-1,1
         DO IX=-1,1
         I3=I+IX
         J3=J+JY
         K3=K+KZ
         IF (I3.LT.1) I3=I3+NX
         IF (I3.GT.NX) I3=I3-NX
         IF (J3.LT.1) J3=J3+NY
         IF (J3.GT.NY) J3=J3-NY
         IF (K3.LT.1) K3=K3+NZ
         IF (K3.GT.NZ) K3=K3-NZ

         U(I3,J3,K3)=MASAP(IP)*VX(IX)*VY(JY)*VZ(KZ)+U(I3,J3,K3)

         END DO
         END DO
         END DO

      ENDDO !loop on particles

!************************************************************************
   END SUBROUTINE DENS_INTERP_TSC
!************************************************************************


!********************************************************************************
   SUBROUTINE VEL_INTERP_TSC(NX,NY,NZ,NPARTT,RXPA,RYPA,RZPA, &
                     U2PA,U3PA,U4PA,U2,U3,U4)
!********************************************************************************
      USE COMMONDATA, ONLY: DX,DY,DZ,RADX,RADY,RADZ,PARTIRED

      IMPLICIT NONE
      !input variables
      INTEGER NX,NY,NZ
      INTEGER(KIND=8) :: NPARTT, LOW1, IP
      REAL*4 RXPA(PARTIRED),RYPA(PARTIRED),RZPA(PARTIRED)
      REAL*4 U2PA(PARTIRED),U3PA(PARTIRED),U4PA(PARTIRED)

      !local variables
      INTEGER I, J, K
      INTEGER IX, JY, KZ, I3, J3, K3
      REAL*4 RADX1, RADXNX, RADY1, RADYNY, RADZ1, RADZNZ, LENG
      REAL*4 VX(-1:1), VY(-1:1), VZ(-1:1)
      REAL*4 BAS, WEI
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

      !smallest/largest coordinates of the coarse grid
      RADX1=RADX(1)-0.5*DX
      RADXNX=RADX(NX)+0.5*DX
      RADY1=RADY(1)-0.5*DY
      RADYNY=RADY(NY)+0.5*DY
      RADZ1=RADZ(1)-0.5*DZ
      RADZNZ=RADZ(NZ)+0.5*DZ
      LENG=RADXNX-RADX1

   !****   from lagrangian velocities (U2PA, U3PA, U4PA)  to eulerian using TSC
      DO I=-1,1
         VX(I)=0.0
         VY(I)=0.0
         VZ(I)=0.0
      END DO

      !LOW1=SUM(NPART(0:NLEVELS))
      LOW1=NPARTT

   !$OMP PARALLEL DO SHARED(LOW1,RXPA,RYPA,RZPA,RADX1,RADY1,RADZ1, &
   !$OMP     RADXNX,RADYNY,RADZNZ,LENG,DX,DY,DZ,NX,NY,NZ, &
   !$OMP     U2PA,U3PA,U4PA,U2,U3,U4,WT,RADX,RADY,RADZ), &
   !$OMP PRIVATE(IP,BAS,I,J,K,VX,VY,VZ,IX,JY,KZ,I3,J3,K3), &
   !$OMP DEFAULT(NONE)
      DO IP=1, LOW1

         !move by LENG particles outside the grid
         IF(RXPA(IP).LT.RADX1) RXPA(IP)=RXPA(IP)+LENG
         IF(RXPA(IP).GE.RADXNX) RXPA(IP)=RXPA(IP)-LENG
         IF(RYPA(IP).LT.RADY1) RYPA(IP)=RYPA(IP)+LENG
         IF(RYPA(IP).GE.RADYNY) RYPA(IP)=RYPA(IP)-LENG
         IF(RZPA(IP).LT.RADZ1) RZPA(IP)=RZPA(IP)+LENG
         IF(RZPA(IP).GE.RADZNZ) RZPA(IP)=RZPA(IP)-LENG

         !initial position
         !I,J,K: cell of the particle in the grid
         BAS=RXPA(IP)
         I=INT(((BAS-RADX(1))/DX)+0.49999) + 1
         BAS=RYPA(IP)
         J=INT(((BAS-RADY(1))/DY)+0.49999) + 1
         BAS=RZPA(IP)
         K=INT(((BAS-RADZ(1))/DZ)+0.49999) + 1

   !***  using a cell in cloud scheme (TSC) to transform lagrangian velocities (DM part) into eulerian
         BAS=ABS(RADX(I-1)-RXPA(IP))
         VX(-1)=0.5*(1.5-BAS/DX)**2
         BAS=ABS(RADX(I)-RXPA(IP))
         VX(0)=0.75-(BAS/DX)**2
         BAS=ABS(RADX(I+1)-RXPA(IP))
         VX(1)=0.5*(1.5-BAS/DX)**2

         BAS=ABS(RADY(J-1)-RYPA(IP))
         VY(-1)=0.5*(1.5-BAS/DY)**2
         BAS=ABS(RADY(J)-RYPA(IP))
         VY(0)=0.75-(BAS/DY)**2
         BAS=ABS(RADY(J+1)-RYPA(IP))
         VY(1)=0.5*(1.5-BAS/DY)**2

         BAS=ABS(RADZ(K-1)-RZPA(IP))
         VZ(-1)=0.5*(1.5-BAS/DZ)**2
         BAS=ABS(RADZ(K)-RZPA(IP))
         VZ(0)=0.75-(BAS/DZ)**2
         BAS=ABS(RADZ(K+1)-RZPA(IP))
         VZ(1)=0.5*(1.5-BAS/DZ)**2


         DO KZ=-1,1
         DO JY=-1,1
         DO IX=-1,1
         I3=I+IX
         J3=J+JY
         K3=K+KZ
         IF (I3.LT.1) I3=I3+NX
         IF (I3.GT.NX) I3=I3-NX
         IF (J3.LT.1) J3=J3+NY
         IF (J3.GT.NY) J3=J3-NY
         IF (K3.LT.1) K3=K3+NZ
         IF (K3.GT.NZ) K3=K3-NZ

         U2(I3,J3,K3)=U2PA(IP)*VX(IX)*VY(JY)*VZ(KZ)+U2(I3,J3,K3)
         U3(I3,J3,K3)=U3PA(IP)*VX(IX)*VY(JY)*VZ(KZ)+U3(I3,J3,K3)
         U4(I3,J3,K3)=U4PA(IP)*VX(IX)*VY(JY)*VZ(KZ)+U4(I3,J3,K3)
         WT(I3,J3,K3)=VX(IX)*VY(JY)*VZ(KZ)+WT(I3,J3,K3)

         END DO
         END DO
         END DO

      ENDDO !loop on particles

      WRITE(*,*) '     Fraction of cells without particles (TSC):', REAL(COUNT(WT .EQ. 0.))/REAL(NX*NY*NZ)

      DO K=1,NX
      DO J=1,NY
      DO I=1,NZ
      WEI=WT(I,J,K)
      IF(WT(I,J,K)==0.) WEI=1.D0
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
   SUBROUTINE H_DISTANCE(NX,NY,NZ,NPARTT,TREE,SMASK,HPART)
   !this subroutine calculates the kernel distance h(x) for every cell x
   !required by DDENS_INTERP_SPH if VVEL_INTERP_SPH is not used
!********************************************************************************  
         USE COSMOKDTREE
         USE COMMONDATA, ONLY: KNEIGHBOURS,DX,DY,DZ,RADX,RADY,RADZ,PARTIRED

         IMPLICIT NONE
         !input variables
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
         !$OMP                   TREE,DX,DY,DZ,RADX,RADY,RADZ) &
         !$OMP PRIVATE(I,IX,JY,KZ,QUERY,TAR,DIST,NEIGH, & 
         !$OMP               CONTA), REDUCTION(MAX:HPART) &
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
                     IF (DIST(I).GT.HPART(NEIGH(I))) THEN
                        HPART(NEIGH(I)) = DIST(I)
                     END IF
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
   SUBROUTINE DDENS_INTERP_SPH(NX,NY,NZ,NPARTT,RXPA,RYPA,RZPA, &
                  HPART,MASAP,U)
!using HPART, smoothing lenght given
!by the furthest cell to which the particle contributes in the
!velocity interpolation
!********************************************************************************
         USE COMMONDATA, ONLY: DX,DY,DZ,RADX,RADY,RADZ,PARTIRED

         IMPLICIT NONE
         !input variables
         INTEGER NX,NY,NZ
         INTEGER(KIND=8) :: NPARTT, I
         REAL*4 RXPA(PARTIRED),RYPA(PARTIRED),RZPA(PARTIRED)
         REAL*4 MASAP(PARTIRED)

         !LOCAL
         INTEGER :: IX,JY,KZ,IX2,JY2,KZ2

         !grid position for each particle
         REAL :: RX1,RY1,RZ1
         INTEGER :: XGRIDPOS(NPARTT),YGRIDPOS(NPARTT),ZGRIDPOS(NPARTT)

         !INTERPOLATION VARIABLES
         REAL*4 :: H,DDIST,NORM
         INTEGER :: INCREX,INCREY,INCREZ
         
         !SPH related
         REAL*4 :: HPART(NPARTT)

         ! !REAL*8 precision auxiliary U
         ! REAL*8, ALLOCATABLE :: U8(:,:,:)

         !output
         REAL*4 U(NX,NY,NZ)

         !! Locate particles in grid !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !grid first cell center
         RX1 = RADX(1)
         RY1 = RADY(1)
         RZ1 = RADZ(1)   
         !$OMP PARALLEL DO SHARED(NPARTT,RXPA,RYPA,RZPA,DX,DY,DZ, &
         !$OMP                    RX1,RY1,RZ1,XGRIDPOS,YGRIDPOS,ZGRIDPOS), &
         !$OMP            PRIVATE(I,IX,JY,KZ), DEFAULT(NONE)
         DO I=1,NPARTT 
            IX = INT((RXPA(I)-RX1)/DX)+1
            JY = INT((RYPA(I)-RY1)/DY)+1
            KZ = INT((RZPA(I)-RZ1)/DZ)+1
            XGRIDPOS(I) = IX
            YGRIDPOS(I) = JY
            ZGRIDPOS(I) = KZ
         END DO
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         ! !Alloc U8
         ! ALLOCATE(U8(NX,NY,NZ))

         !$OMP PARALLEL DO SHARED(NPARTT,XGRIDPOS,YGRIDPOS,ZGRIDPOS,HPART, &
         !$OMP                   RXPA,RYPA,RZPA,MASAP,RADX,RADY,RADZ, &
         !$OMP                   DX,DY,DZ,NX,NY,NZ,U), &
         !$OMP PRIVATE(I,IX,JY,KZ,IX2,JY2,KZ2,H,INCREX,INCREY,INCREZ,DDIST,NORM), &
         !$OMP SCHEDULE(DYNAMIC), DEFAULT(NONE)
         DO I=1,NPARTT
            IX = XGRIDPOS(I)
            JY = YGRIDPOS(I)
            KZ = ZGRIDPOS(I)
            H = HPART(I)
            INCREX = INT(2*H/DX + 0.5) !0.5 accounts for the fact that
            INCREY = INT(2*H/DY + 0.5) !we are using the center of the cell
            INCREZ = INT(2*H/DZ + 0.5) !as the reference point

            !First, find normalisation factor
            NORM = 0.
            DO KZ2=KZ-INCREZ,KZ+INCREZ
               IF (KZ2.LT.1 .OR. KZ2.GT.NZ) CYCLE
               DO JY2=JY-INCREY,JY+INCREY
                  IF (JY2.LT.1 .OR. JY2.GT.NY) CYCLE
                  DO IX2=IX-INCREX,IX+INCREX
                     IF (IX2.LT.1 .OR. IX2.GT.NX) CYCLE
                        !DISTANCE BETWEEN PARTICLE AND GRID POINT
                        DDIST = ( (RXPA(I)-RADX(IX2))**2 + &
                                  (RYPA(I)-RADY(JY2))**2 + &
                                  (RZPA(I)-RADZ(KZ2))**2 )**0.5

                        !KERNEL FUNCTION
                        CALL KERNEL_FUNC_2(H,DDIST)
                        
                        NORM = NORM + DDIST
                  ENDDO
               ENDDO
            ENDDO

            ! Particle only contributes to its own cell
            IF (NORM .EQ. 0.) THEN
               !$OMP ATOMIC
               U(IX,JY,KZ) =  U(IX,JY,KZ) + MASAP(I)
               CYCLE
            ENDIF

            !Now, interpolate
            DO KZ2=KZ-INCREZ,KZ+INCREZ
               IF (KZ2.LT.1 .OR. KZ2.GT.NZ) CYCLE
               DO JY2=JY-INCREY,JY+INCREY
                  IF (JY2.LT.1 .OR. JY2.GT.NY) CYCLE
                  DO IX2=IX-INCREX,IX+INCREX
                     IF (IX2.LT.1 .OR. IX2.GT.NX) CYCLE
                        !DISTANCE BETWEEN PARTICLE AND GRID POINT
                        DDIST = ( (RXPA(I)-RADX(IX2))**2 + &
                                 (RYPA(I)-RADY(JY2))**2 + &
                                 (RZPA(I)-RADZ(KZ2))**2 )**0.5

                        !KERNEL FUNCTION
                        CALL KERNEL_FUNC_2(H,DDIST)
                        
                        !Update
                        !$OMP ATOMIC
                        U(IX2,JY2,KZ2) = U(IX2,JY2,KZ2) + DDIST*MASAP(I)/NORM
                  ENDDO
               ENDDO
            ENDDO
         END DO

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         ! !CONVERT REAL*8 TO REAL*4 AFTER INTERPOLATION
         ! U = REAL(U8, KIND=4)
         ! DEALLOCATE(U8)

         !FOR THE MOMENT U IS MASS

!************************************************************************
   END SUBROUTINE DDENS_INTERP_SPH
!************************************************************************


!********************************************************************************
   SUBROUTINE VVEL_INTERP_SPH(NX,NY,NZ,NPARTT,TREE, &
                        HPART,MASAP,U2PA,U3PA,U4PA,SMASK,U2,U3,U4)
!********************************************************************************
         USE COSMOKDTREE 
         USE COMMONDATA, ONLY: KNEIGHBOURS,DX,DY,DZ,RADX,RADY,RADZ,PARTIRED

         IMPLICIT NONE
         !input variables
         INTEGER NX,NY,NZ
         INTEGER(KIND=8) :: NPARTT, CONTA, I
         REAL*4 U2PA(PARTIRED),U3PA(PARTIRED),U4PA(PARTIRED)
         REAL*4 MASAP(PARTIRED)
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
         !$OMP                   TREE,DX,DY,DZ,RADX,RADY,RADZ,SMASK) &
         !$OMP PRIVATE(I,IX,JY,KZ,QUERY,TAR,DIST,DIST8,NEIGH, & 
         !$OMP               HKERN,BAS8,BAS8X,BAS8Y,BAS8Z,CONTA), REDUCTION(MAX:HPART) &
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
                     IF (DIST(I).GT.HPART(NEIGH(I))) THEN
                        HPART(NEIGH(I)) = DIST(I)
                     END IF
                  END DO
                  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                  !cell's smoothing length
                  HKERN=DIST(CONTA)
                  !KERNEL FUNCTION
                  CALL KERNEL_FUNC(CONTA,CONTA,HKERN,DIST)

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
   END SUBROUTINE VVEL_INTERP_SPH
!************************************************************************


!********************************************************************************
   SUBROUTINE PPART_DENS(NPARTT,TREE,MASAP,PART_DENS)
! Obtains the density at the position of each particle using Nearest Neighbours
!********************************************************************************
         USE COSMOKDTREE 
         USE COMMONDATA, ONLY: KNEIGHBOURS,PARTIRED, RXPA,RYPA,RZPA,PI

         IMPLICIT NONE
         !input variables
         INTEGER(KIND=8) :: NPARTT, CONTA, I, J
         REAL*4 MASAP(PARTIRED)
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
   SUBROUTINE VVEL_INTERP_SPH_VW(NX,NY,NZ,NPARTT,TREE, &
                        HPART,PART_DENS,MASAP,U2PA,U3PA,U4PA,SMASK,U2,U3,U4)
!
! SAME AS VVEL_INTERP_SPH BUT USING VOLUME WEIGHTING instead of MASS
!********************************************************************************
         USE COSMOKDTREE 
         USE COMMONDATA, ONLY: KNEIGHBOURS,DX,DY,DZ,RADX,RADY,RADZ,PARTIRED, &   
                              RXPA,RYPA,RZPA,PI

         IMPLICIT NONE
         !input variables
         INTEGER NX,NY,NZ
         INTEGER(KIND=8) :: NPARTT, CONTA, I, J
         REAL*4 U2PA(PARTIRED),U3PA(PARTIRED),U4PA(PARTIRED)
         REAL*4 MASAP(PARTIRED)
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
         !$OMP                   TREE,DX,DY,DZ,RADX,RADY,RADZ,PART_DENS,SMASK) &
         !$OMP PRIVATE(I,IX,JY,KZ,QUERY,TAR,DIST,DIST8,NEIGH, & 
         !$OMP               HKERN,BAS8,BAS8X,BAS8Y,BAS8Z,CONTA), REDUCTION(MAX:HPART) &
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
                     IF (DIST(I).GT.HPART(NEIGH(I))) THEN
                        HPART(NEIGH(I)) = DIST(I)
                     END IF
                  END DO
                  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                  !cell's smoothing length
                  HKERN=DIST(CONTA)
                  !KERNEL FUNCTION
                  CALL KERNEL_FUNC(CONTA,CONTA,HKERN,DIST)

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
   SUBROUTINE KERNEL_FUNC(N,N2,W,DIST)
!************************************************************************
   !*     DIST contains initially the distance (particle to cell), and it is
   !*     updated with the (unnormalised) value of the kernel
         IMPLICIT NONE
         INTEGER(KIND=8) I,N,N2 ! N is the dimension of the array dist;
                     ! N2, the actual number of particles filled in
         REAL W,DIST(N) ! W is the total width of the kernel 
         REAL H

         REAL DISTS

   !     Cubic spline kernel (M4)
         H=W/2.0
         DO I=1,N2
         DISTS=DIST(I)/H !q
         IF (DISTS.LE.1.0) THEN
         DIST(I)=1.0-1.5*DISTS**2*(1.-0.5*DISTS)
         ELSE IF (DISTS.LE.2.0) THEN
         DIST(I)=0.25*(2.-DISTS)**3
         ELSE
         DIST(I)=0.0
         END IF
         END DO

         RETURN
!************************************************************************
   END SUBROUTINE KERNEL_FUNC
!************************************************************************


!************************************************************************
   SUBROUTINE KERNEL_FUNC_2(W,DIST)
!************************************************************************
   !*     DIST contains initially the distance (particle to cell)
   !*    and it is updated with the (unnormalised) value of the kernel
            IMPLICIT NONE
            REAL W,DIST,H
      
   !     Cubic spline kernel (M4)
            H=W/2.0
            DIST=DIST/H !q
            IF (DIST.LE.1.0) THEN
            DIST=1.0-1.5*DIST**2*(1.-0.5*DIST)
            ELSE IF (DIST.LE.2.0) THEN
            DIST=0.25*(2.-DIST)**3
            ELSE
            DIST=0.0
            END IF
      
            RETURN
!************************************************************************
   END SUBROUTINE KERNEL_FUNC_2
!************************************************************************


! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

END MODULE PARTICLES

! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~