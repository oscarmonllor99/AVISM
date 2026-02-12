
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

MODULE VOIDFINDING
   IMPLICIT NONE
   PUBLIC

! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

CONTAINS 

!********************************************************************************
SUBROUTINE MARK_ALL(LOW1,LOW2) !valid for all hierarchies
!********************************************************************************
     USE COMMONDATA, ONLY: MARCAP2, FLAGV, FLAG_SUB, DENS_THRE, U1CO, DIVERCO
     IMPLICIT NONE
     INTEGER LOW1,LOW2
     INTEGER IX, JY, KZ, F1, F2

     !MARKS CELL AS CANDIDATE FOR CENTRE
     ALLOCATE(FLAGV(LOW1:LOW2,LOW1:LOW2,LOW1:LOW2))
     !TELLS WHICH CELLS ARE ALLOWED TO GROW A SUBVOID
     ALLOCATE(FLAG_SUB(LOW1:LOW2,LOW1:LOW2,LOW1:LOW2))
     !--> if lev=LEVMIN: FLAG_SUB=1 in every cell

     FLAGV(:,:,:)=0
     FLAG_SUB(:,:,:)=0

     DO KZ=LOW1, LOW2
        DO JY=LOW1, LOW2
           DO IX=LOW1, LOW2
              IF(MARCAP2(IX,JY,KZ) .NE. 0 ) THEN 

                   FLAG_SUB(IX,JY,KZ) = 1 !can be part of a subvoid

                   F1 = 0
                   F2 = 0

                   !check density
                   IF (U1CO(IX,JY,KZ) .LT. DENS_THRE+1.) F1 = 1   

                   !check divergence
                   IF (DIVERCO(IX,JY,KZ) .GT. 0.) F2 = 1

                   !mark as void cell
                   IF(F1 .EQ. 1 .AND. F2 .EQ. 1) FLAGV(IX,JY,KZ)=1 
                ENDIF
           ENDDO
        ENDDO
     ENDDO
!***************************************************************************
END SUBROUTINE MARK_ALL
!***************************************************************************


!********************************************************************************
SUBROUTINE VOIDFIND(LEV,LOW1,LOW2,DDX,DDY,DDZ, &
                    RX1,RY1,RZ1,NVOID1,REQP,NVOID)
!Serial routine 
!********************************************************************************
USE COMMONDATA
IMPLICIT NONE
!input variables
INTEGER:: LEV, NVOID1
REAL*4:: DDX, DDY, DDZ, RX1, RY1, RZ1
INTEGER, ALLOCATABLE :: MARCA_AUX(:,:,:)
REAL*4, DIMENSION(NVOID1):: REQP
!local variables
INTEGER:: NXYZ, II, I1, I, J, K
INTEGER:: L1, IX, JY, KZ, IXX, JYY, KZZ
REAL*4:: DIJK, VV, REQ0
REAL*4 DDENS
REAL*4 DDENS2
REAL*4, ALLOCATABLE:: DDD(:), DDD2(:)
REAL*4 A,B,C
REAL*4 :: DIVERCO2(LOW1:LOW2,LOW1:LOW2,LOW1:LOW2)
INTEGER, ALLOCATABLE:: DDDX(:), DDDY(:), DDDZ(:), &
    DDDX2(:), DDDY2(:), DDDZ2(:), INDICE(:), INDICE2(:)
INTEGER:: FLAG_DIV_P, FLAG_NEXT_GRAD_P, NVOID_MAX_P, NSEC_P
INTEGER:: INDV
INTEGER :: FLAG,FLAG1,FLAG2,FLAGX1,FLAGX2,FLAGY1,FLAGY2,FLAGZ1,FLAGZ2
INTEGER, ALLOCATABLE :: FLAG1_GRID(:,:,:),FLAGX1_GRID(:,:,:),FLAGX2_GRID(:,:,:)
INTEGER, ALLOCATABLE :: FLAGY1_GRID(:,:,:),FLAGY2_GRID(:,:,:),FLAGZ1_GRID(:,:,:),FLAGZ2_GRID(:,:,:)

!parameters
INTEGER NSEC

!PBC:
INTEGER :: LOW1, LOW2

!output variables
INTEGER:: NVOID

WRITE(*,*) 'Within void finder, level:', LEV, 'LOW1, LOW2:', LOW1, LOW2

!$$$$$$$$$$$$
!Set NSEC 
NSEC=1
!$$$$$$$$$$$$

!max divergence
NXYZ=0
DIVERCO2=DIVERCO*REAL(FLAGV)

!new DIVERCO2 is =DIVERCO if DIVERCO>0 and 0 elsewhere
NXYZ=COUNT(DIVERCO2.GT.0)
WRITE(*,*) '       NXYZ:', NXYZ

ALLOCATE(DDD(NXYZ))
ALLOCATE(DDD2(NXYZ))
ALLOCATE(DDDX2(NXYZ))
ALLOCATE(DDDY2(NXYZ))
ALLOCATE(DDDZ2(NXYZ))
ALLOCATE(DDDX(NXYZ))
ALLOCATE(DDDY(NXYZ))
ALLOCATE(DDDZ(NXYZ))
ALLOCATE(INDICE2(NXYZ))
ALLOCATE(INDICE(NXYZ))

!$OMP PARALLEL DO SHARED(NXYZ,DDD, DDD2,DDDX2,DDDZ2,DDDX,DDDY,DDDZ,INDICE, INDICE2), &
!$OMP     PRIVATE(I)
DO I=1, NXYZ
 DDD(I)=0.0 !DIVERCO collapsed to void candidate cells (1D array)
 DDD2(I)=0.0 !DDD ordered from largest to smallest, temp
 DDDX2(I)=0 !x index temp
 DDDY2(I)=0 !y index temp
 DDDZ2(I)=0 !z index temp
 DDDX(I)=0 !x index
 DDDY(I)=0 !y index
 DDDZ(I)=0 !z index
 INDICE2(I)=0 !indices of sorted diver (1=smallest, NXYZ=largest)
 INDICE(I)=0 !indices of sorted diver (1=largest, NXYZ=smallest)
ENDDO

!*parallel
!variables for parallelization:
FLAG_NEXT_GRAD_P=FLAG_NEXT_GRAD
FLAG_DIV_P=FLAG_DIV_EDGE
NVOID_MAX_P=NVOID_MAX
NSEC_P=NSEC

!$OMP PARALLEL DO SHARED(NVOID_MAX_P,INICIOX,FINALX,INICIOY,FINALY, &
!$OMP     INICIOZ, FINALZ, ICX, ICY, ICZ), &
!$OMP     PRIVATE(I)
DO I=1, NVOID_MAX_P
 INICIOX(I)=0
 FINALX(I)=0
 INICIOY(I)=0
 FINALY(I)=0
 INICIOZ(I)=0
 FINALZ(I)=0
 ICX(I)=0
 ICY(I)=0
 ICZ(I)=0
ENDDO

!collapse DIVERCO onto DDD
II=0
DO K=LOW1,LOW2
DO J=LOW1,LOW2
DO I=LOW1,LOW2
DIJK=0.0
DIJK=DIVERCO2(I,J,K) !max diver
IF(DIJK .GT. 0.0) THEN ! div
 II=II+1 
 DDD(II)=DIJK
 DDDX(II)=I
 DDDY(II)=J
 DDDZ(II)=K
ENDIF
ENDDO
ENDDO
ENDDO

IF (II.NE.NXYZ) THEN
 WRITE(*,*) '       WARNING,NXYZ', II, NXYZ
 STOP
END IF

!!!  ORDERING DIVERGENCE of candidate cells: 
!!!                from largest to smallest: We want centers in cells of max diver
CALL INDEXX(NXYZ,DDD(1:NXYZ),INDICE2(1:NXYZ)) !sort from smallest to largest
write(*,*)'       Min, max values for divergence:', minval(DDD), maxval(DDD)

!max diver: reverse order
I1=NXYZ

!INDICE2 --> DDD ordered from smallest to largest
!INDICE --> DDD ordered from largest to smallest

DO I=1,NXYZ
   DDD2(I1)=DDD(INDICE2(I)) !DIVERCO ordered from largest to smallest
   DDDX2(I1)=DDDX(INDICE2(I))
   DDDY2(I1)=DDDY(INDICE2(I))
   DDDZ2(I1)=DDDZ(INDICE2(I))
   INDICE(I1)=INDICE2(I)
   I1=I1-1
END DO

!$OMP PARALLEL DO SHARED(NXYZ,DDD,DDDX,DDDY,DDDZ), &
!$OMP     PRIVATE(I)
DO I=1,NXYZ
 DDD(I)=0.0
 DDDX(I)=0
 DDDY(I)=0
 DDDZ(I)=0
END DO

!$OMP PARALLEL DO SHARED(NXYZ,DDD,DDD2,DDDX2,DDDY2,DDDZ2,&
!$OMP     DDDX,DDDY,DDDZ,INDICE), &
!$OMP     PRIVATE(I)
DO I=1,NXYZ
 DDD(I)=DDD2(I)
 DDDX(I)=DDDX2(I)
 DDDY(I)=DDDY2(I)
 DDDZ(I)=DDDZ2(I)
END DO
!DDD ordered from largest to smallest and same with DDDX, DDDY, DDDZ
DEALLOCATE(DDD2, DDDX2, DDDY2, DDDZ2, INDICE2)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!THIS VARIABLES WILL TELL IF A CELL CAN BE EXPANDED AS A VOID IN EACH DIRECTION
ALLOCATE(FLAGX1_GRID(LOW1:LOW2,LOW1:LOW2,LOW1:LOW2))
ALLOCATE(FLAGX2_GRID(LOW1:LOW2,LOW1:LOW2,LOW1:LOW2))
ALLOCATE(FLAGY1_GRID(LOW1:LOW2,LOW1:LOW2,LOW1:LOW2))
ALLOCATE(FLAGY2_GRID(LOW1:LOW2,LOW1:LOW2,LOW1:LOW2))
ALLOCATE(FLAGZ1_GRID(LOW1:LOW2,LOW1:LOW2,LOW1:LOW2))
ALLOCATE(FLAGZ2_GRID(LOW1:LOW2,LOW1:LOW2,LOW1:LOW2))
ALLOCATE(FLAG1_GRID(LOW1:LOW2,LOW1:LOW2,LOW1:LOW2))

!$OMP PARALLEL SHARED(LOW1,LOW2,FLAGX1_GRID,FLAGX2_GRID, &
!$OMP FLAGY1_GRID,FLAGY2_GRID,FLAGZ1_GRID,FLAGZ2_GRID,FLAG1_GRID), &
!$OMP     PRIVATE(IX,JY,KZ), DEFAULT(NONE)
!$OMP DO ORDERED 
   DO KZ=LOW1,LOW2
      DO JY=LOW1,LOW2
         DO IX=LOW1,LOW2
            FLAG1_GRID(IX,JY,KZ)=0
            FLAGX1_GRID(IX,JY,KZ)=0
            FLAGX2_GRID(IX,JY,KZ)=0
            FLAGY1_GRID(IX,JY,KZ)=0
            FLAGY2_GRID(IX,JY,KZ)=0
            FLAGZ1_GRID(IX,JY,KZ)=0
            FLAGZ2_GRID(IX,JY,KZ)=0
         ENDDO
      ENDDO
   ENDDO
!$OMP END DO
!$OMP END PARALLEL


!$OMP PARALLEL SHARED(LOW1,LOW2,DIVERCO,U1CO,FLAG_DIV_P,FLAG_NEXT_GRAD_P,DIV_THRE,&
!$OMP FLAGX1_GRID,FLAGX2_GRID,FLAGY1_GRID,FLAGY2_GRID,FLAGZ1_GRID,FLAGZ2_GRID, &
!$OMP FLAG1_GRID,MARCAP2,REQP,RMIN_SUB,FLAG_SUB,DDX,DDY,DDZ,GRAD_THRE,DENS_THRE2,NSEC_P), &
!$OMP     PRIVATE(II,IX,JY,KZ,REQ0,INDV,DDENS,DDENS2), DEFAULT(NONE)
!$OMP DO REDUCTION(+:FLAGX1_GRID,FLAGX2_GRID,FLAGY1_GRID,FLAGY2_GRID, &
!$OMP                FLAGZ1_GRID,FLAGZ2_GRID,FLAG1_GRID)
   DO KZ=LOW1,LOW2
      DO JY=LOW1,LOW2
         DO IX=LOW1,LOW2
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !!!!!!!!!!!! Search for subvoids only in large voids
            INDV=MARCAP2(IX, JY, KZ) !ID of the parent void when looking for subvoids
                                    !   =-1 when looking for voids
            REQ0=0.
            IF(INDV .GT. 0) THEN 
               REQ0=REQP(INDV)
               IF(REQ0 .LT. RMIN_SUB) THEN
                  FLAG1_GRID(IX,JY,KZ)=1
               ENDIF
            ENDIF
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            !+X
            !!!!!!!!!!!!
            IF(IX .GE. LOW2-NSEC_P) FLAGX2_GRID(IX,JY,KZ)=1
            IF(FLAGX2_GRID(IX,JY,KZ) .EQ. 0) THEN
               II = IX+1
               !+X GRADIENT
               IF(II .EQ. LOW2) THEN
                  !only check divergence and density
                  IF(FLAG_DIV_P==1 .AND. DIVERCO(II,JY,KZ) .LT. DIV_THRE) FLAGX2_GRID(IX,JY,KZ)=1
                  IF(U1CO(II,JY,KZ) .GE. DENS_THRE2+1.) FLAGX2_GRID(IX,JY,KZ)=1
               ELSE
                  DDENS = (U1CO(II+1,JY,KZ)-U1CO(II-1,JY,KZ))/(2.D0*DDX)
                  IF(DDENS .GE. GRAD_THRE) FLAGX2_GRID(IX,JY,KZ)=1
                  !check gradient in the next cell
                  IF(FLAG_NEXT_GRAD_P==1 .AND. FLAGX2_GRID(IX,JY,KZ)==1 .AND. II .LT. LOW2-1) THEN
                     DDENS2=(U1CO(II+2,JY,KZ)-U1CO(II,JY,KZ))/(2.D0*DDX)
                     IF(DDENS2 .LT. GRAD_THRE) FLAGX2_GRID(IX,JY,KZ)=0
                  ENDIF
                  IF(FLAG_DIV_P==1 .AND. DIVERCO(II,JY,KZ) .LT. DIV_THRE) FLAGX2_GRID(IX,JY,KZ)=1
                  IF(U1CO(II,JY,KZ) .GE. DENS_THRE2+1.) FLAGX2_GRID(IX,JY,KZ)=1
               ENDIF
            ENDIF
            !!!!!!!!!!!!

            !-X
            !!!!!!!!!!!!
            IF(IX .LE. LOW1+NSEC_P) FLAGX1_GRID(IX,JY,KZ)=1
            IF(FLAGX1_GRID(IX,JY,KZ) .EQ. 0) THEN
               II = IX-1
               !+X GRADIENT
               IF(II .EQ. LOW1) THEN
                  !only check divergence and density
                  IF(FLAG_DIV_P==1 .AND. DIVERCO(II,JY,KZ) .LT. DIV_THRE) FLAGX1_GRID(IX,JY,KZ)=1
                  IF(U1CO(II,JY,KZ) .GE. DENS_THRE2+1.) FLAGX1_GRID(IX,JY,KZ)=1
               ELSE
                  DDENS = (U1CO(II-1,JY,KZ)-U1CO(II+1,JY,KZ))/(2.D0*DDX)
                  IF(DDENS .GE. GRAD_THRE) FLAGX1_GRID(IX,JY,KZ)=1
                  !check gradient in the next cell
                  IF(FLAG_NEXT_GRAD_P==1 .AND. FLAGX1_GRID(IX,JY,KZ)==1 .AND. II .GT. LOW1+1) THEN 
                     DDENS2=(U1CO(II-2,JY,KZ)-U1CO(II,JY,KZ))/(2.D0*DDX)
                     IF(DDENS2 .LT. GRAD_THRE) FLAGX1_GRID(IX,JY,KZ)=0
                  ENDIF
                  IF(FLAG_DIV_P==1 .AND. DIVERCO(II,JY,KZ) .LT. DIV_THRE) FLAGX1_GRID(IX,JY,KZ)=1
                  IF(U1CO(II,JY,KZ) .GE. DENS_THRE2+1.) FLAGX1_GRID(IX,JY,KZ)=1
               ENDIF
            ENDIF
            !!!!!!!!!!!!

            !+Y
            !!!!!!!!!!!!
            IF(JY .GE. LOW2-NSEC_P) FLAGY2_GRID(IX,JY,KZ)=1
            IF(FLAGY2_GRID(IX,JY,KZ) .EQ. 0) THEN
               II = JY+1
               !+Y GRADIENT
               IF(II .EQ. LOW2) THEN
                  !only check divergence and density
                  IF(FLAG_DIV_P==1 .AND. DIVERCO(IX,II,KZ) .LT. DIV_THRE) FLAGY2_GRID(IX,JY,KZ)=1
                  IF(U1CO(IX,II,KZ) .GE. DENS_THRE2+1.) FLAGY2_GRID(IX,JY,KZ)=1
               ELSE
                  DDENS = (U1CO(IX,II+1,KZ)-U1CO(IX,II-1,KZ))/(2.D0*DDY)
                  IF(DDENS .GE. GRAD_THRE) FLAGY2_GRID(IX,JY,KZ)=1
                  !check gradient in the next cell
                  IF(FLAG_NEXT_GRAD_P==1 .AND. FLAGY2_GRID(IX,JY,KZ)==1 .AND. II .LT. LOW2-1) THEN 
                     DDENS2=(U1CO(IX,II+2,KZ)-U1CO(IX,II,KZ))/(2.D0*DDY)
                     IF(DDENS2 .LT. GRAD_THRE) FLAGY2_GRID(IX,JY,KZ)=0
                  ENDIF
                  IF(FLAG_DIV_P==1 .AND. DIVERCO(IX,II,KZ) .LT. DIV_THRE) FLAGY2_GRID(IX,JY,KZ)=1
                  IF(U1CO(IX,II,KZ) .GE. DENS_THRE2+1.) FLAGY2_GRID(IX,JY,KZ)=1
               ENDIF
            ENDIF
            !!!!!!!!!!!!

            !-Y
            !!!!!!!!!!!!
            IF(JY .LE. LOW1+NSEC_P) FLAGY1_GRID(IX,JY,KZ)=1
            IF(FLAGY1_GRID(IX,JY,KZ) .EQ. 0) THEN
               II = JY-1
               !+Y GRADIENT
               IF(II .EQ. LOW1) THEN
                  !only check divergence and density
                  IF(FLAG_DIV_P==1 .AND. DIVERCO(IX,II,KZ) .LT. DIV_THRE) FLAGY1_GRID(IX,JY,KZ)=1
                  IF(U1CO(IX,II,KZ) .GE. DENS_THRE2+1.) FLAGY1_GRID(IX,JY,KZ)=1
               ELSE
                  DDENS = (U1CO(IX,II-1,KZ)-U1CO(IX,II+1,KZ))/(2.D0*DDY)
                  IF(DDENS .GE. GRAD_THRE) FLAGY1_GRID(IX,JY,KZ)=1
                  !check gradient in the next cell
                  IF(FLAG_NEXT_GRAD_P==1 .AND. FLAGY1_GRID(IX,JY,KZ)==1 .AND. II .GT. LOW1+1) THEN 
                     DDENS2=(U1CO(IX,II-2,KZ)-U1CO(IX,II,KZ))/(2.D0*DDY)
                     IF(DDENS2 .LT. GRAD_THRE) FLAGY1_GRID(IX,JY,KZ)=0
                  ENDIF
                  IF(FLAG_DIV_P==1 .AND. DIVERCO(IX,II,KZ) .LT. DIV_THRE) FLAGY1_GRID(IX,JY,KZ)=1
                  IF(U1CO(IX,II,KZ) .GE. DENS_THRE2+1.) FLAGY1_GRID(IX,JY,KZ)=1
               ENDIF
            ENDIF
            !!!!!!!!!!!!

            !+Z
            !!!!!!!!!!!!
            IF(KZ .GE. LOW2-NSEC_P) FLAGZ2_GRID(IX,JY,KZ)=1
            IF(FLAGZ2_GRID(IX,JY,KZ) .EQ. 0) THEN
               II = KZ+1
               !+Z GRADIENT
               IF(II .EQ. LOW2) THEN
                  !only check divergence and density
                  IF(FLAG_DIV_P==1 .AND. DIVERCO(IX,JY,II) .LT. DIV_THRE) FLAGZ2_GRID(IX,JY,KZ)=1
                  IF(U1CO(IX,JY,II) .GE. DENS_THRE2+1.) FLAGZ2_GRID(IX,JY,KZ)=1
               ELSE
                  DDENS = (U1CO(IX,JY,II+1)-U1CO(IX,JY,II-1))/(2.D0*DDZ)
                  IF(DDENS .GE. GRAD_THRE) FLAGZ2_GRID(IX,JY,KZ)=1
                  !check gradient in the next cell
                  IF(FLAG_NEXT_GRAD_P==1 .AND. FLAGZ2_GRID(IX,JY,KZ)==1 .AND. II .LT. LOW2-1) THEN 
                     DDENS2=(U1CO(IX,JY,II+2)-U1CO(IX,JY,II))/(2.D0*DDZ)
                     IF(DDENS2 .LT. GRAD_THRE) FLAGZ2_GRID(IX,JY,KZ)=0
                  ENDIF
                  IF(FLAG_DIV_P==1 .AND. DIVERCO(IX,JY,II) .LT. DIV_THRE) FLAGZ2_GRID(IX,JY,KZ)=1
                  IF(U1CO(IX,JY,II) .GE. DENS_THRE2+1.) FLAGZ2_GRID(IX,JY,KZ)=1
               ENDIF
            ENDIF
            !!!!!!!!!!!!

            !-Z
            !!!!!!!!!!!!
            IF(KZ .LE. LOW1+NSEC_P) FLAGZ1_GRID(IX,JY,KZ)=1
            IF(FLAGZ1_GRID(IX,JY,KZ) .EQ. 0) THEN
               II = KZ-1
               !+Z GRADIENT
               IF(II .EQ. LOW1) THEN
                  !only check divergence and density
                  IF(FLAG_DIV_P==1 .AND. DIVERCO(IX,JY,II) .LT. DIV_THRE) FLAGZ1_GRID(IX,JY,KZ)=1
                  IF(U1CO(IX,JY,II) .GE. DENS_THRE2+1.) FLAGZ1_GRID(IX,JY,KZ)=1
               ELSE
                  DDENS = (U1CO(IX,JY,II-1)-U1CO(IX,JY,II+1))/(2.D0*DDZ)
                  IF(DDENS .GE. GRAD_THRE) FLAGZ1_GRID(IX,JY,KZ)=1
                  !check gradient in the next cell
                  IF(FLAG_NEXT_GRAD_P==1 .AND. FLAGZ1_GRID(IX,JY,KZ)==1 .AND. II .GT. LOW1+1) THEN 
                     DDENS2=(U1CO(IX,JY,II-2)-U1CO(IX,JY,II))/(2.D0*DDZ)
                     IF(DDENS2 .LT. GRAD_THRE) FLAGZ1_GRID(IX,JY,KZ)=0
                  ENDIF
                  IF(FLAG_DIV_P==1 .AND. DIVERCO(IX,JY,II) .LT. DIV_THRE) FLAGZ1_GRID(IX,JY,KZ)=1
                  IF(U1CO(IX,JY,II) .GE. DENS_THRE2+1.) FLAGZ1_GRID(IX,JY,KZ)=1
               ENDIF
            ENDIF
            !!!!!!!!!!!!

         ENDDO
      ENDDO
   ENDDO
!$OMP END DO
!$OMP END PARALLEL

WRITE(*,*) '       Conditions for cube expansion set'
WRITE(*,*) '       Starting search of void zones'

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  
!$$$$$$$$   loop on potential centers, starting from cells with max diver    $$$$$$$$$$$$$$$

!CUBE COUNT
NVOID=0

!INDICATES IF A CELL IS ALREADY PART OF A CUBE
ALLOCATE(MARCA_AUX(LOW1:LOW2,LOW1:LOW2,LOW1:LOW2))
MARCA_AUX = 0

DO L1=1, NXYZ 

   IX=DDDX(L1)    !!!celda centro
   JY=DDDY(L1)
   KZ=DDDZ(L1)

   !If cell (IX, JY, KZ) falls within another void I exclude it
   FLAG1=0 !
   IF (MARCA_AUX(IX,JY,KZ) .GT. 0) FLAG1=1

   !NSEC -> security border for avoiding voids at the edges
   IF(IX .LE. LOW1+NSEC .OR. JY .LE. LOW1+NSEC .OR. KZ .LE. LOW1+NSEC .OR. &
      IX .GE. LOW2-NSEC .OR. JY .GE. LOW2-NSEC .OR. KZ .GE. LOW2-NSEC) FLAG1=1
   
   !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ 
   IF(FLAG1==0) THEN !cell of the center does not fall in other voids and not at the edge

   NVOID=NVOID+1
   IF(NVOID>NVOID_MAX) THEN
      WRITE(*,*) '       NVOID > NVOID_MAX!! Increase NVOID_MAX '
      STOP
   ENDIF

   FLAG=0 ! if FLAG=1 means that void can not be expanded,
   FLAGX1=0
   FLAGX2=0
   FLAGY1=0
   FLAGY2=0
   FLAGZ1=0
   FLAGZ2=0

   INICIOX(NVOID)=IX !void dimensions
   FINALX(NVOID)=IX
   INICIOY(NVOID)=JY
   FINALY(NVOID)=JY
   INICIOZ(NVOID)=KZ
   FINALZ(NVOID)=KZ

   INDV=MARCAP2(IX,JY,KZ) !ID of the parent void when looking for subvoids
   !                        =-1 when looking for voids
   FATHER(NVOID)=INDV !ID of the parent void 

   IF(FLAG1_GRID(IX,JY,KZ)==1) FLAG=1 !means parent void is too small
   
   !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   !%%%%%   start expanding the cubes in all the directions    %%%%%%%%%%%%%%%%%%%

   FLAG2=0 !for excluding voids after DO WHILE
   !FLAG -> for stopping expansion in DO WHILE

   DO WHILE(FLAG==0) ! if FLAG=1 means that void can not be expanded,
                     ! because does not satisfy the void conditions in at least one face

      !!!!!    expand along  X
      !X!X!X!X!X!X!X!X!X!X!X!X!X!X!X!X!X!X!X!X!X!X!X!X
      !X!X!X!X!X!X!X!X!X!X!X!X!X!X!X!X!X!X!X!X!X!X!X!X
         OUTX:DO KZZ=INICIOZ(NVOID), FINALZ(NVOID)
         DO JYY=INICIOY(NVOID), FINALY(NVOID)
            !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
               IF(FLAGX1==0 .OR. FLAGX2==0) THEN 
               !+X
                  IF(FINALX(NVOID) .GE. LOW2-NSEC)  THEN
                     FLAGX2=1
                     FLAG2 = 1
                  ELSE
                     II=FINALX(NVOID)
                     IF(FLAG_SUB(II+1, JYY, KZZ).EQ.0 .OR. & 
                        MARCAP2(II+1,JYY,KZZ).NE.INDV) FLAGX2=1 !If next cell is from other parent void
                     IF(FLAGX2 .EQ. 0) THEN
                        FLAGX2=FLAGX2_GRID(II,JYY,KZZ) !Check if it can be expanded in this direction
                     ENDIF
                  ENDIF
               !-X
                  IF(INICIOX(NVOID) .LE. LOW1+NSEC)  THEN
                     FLAGX1=1
                     FLAG2=1
                  ELSE
                     II=INICIOX(NVOID)
                     IF(FLAG_SUB(II-1, JYY, KZZ).EQ.0 .OR. &
                        MARCAP2(II-1,JYY,KZZ).NE.INDV) FLAGX1=1
                     IF(FLAGX1 .EQ. 0) THEN
                        FLAGX1=FLAGX1_GRID(II,JYY,KZZ)
                     ENDIF
                  ENDIF
            !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
               ENDIF
         ENDDO
         ENDDO OUTX

      !!!!!    expand along  Y
      !X!X!X!X!X!X!X!X!X!X!X!X!X!X!X!X!X!X!X!X!X!X!X!X
      !X!X!X!X!X!X!X!X!X!X!X!X!X!X!X!X!X!X!X!X!X!X!X!X
         OUTY:DO KZZ=INICIOZ(NVOID), FINALZ(NVOID)
         DO IXX=INICIOX(NVOID), FINALX(NVOID)
            !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
               IF(FLAGY1==0 .OR. FLAGY2==0) THEN
               !+Y
                  IF(FINALY(NVOID) .GE. LOW2-NSEC)  THEN
                     FLAGY2=1
                     FLAG2 = 1
                  ELSE
                     II=FINALY(NVOID)
                     IF(FLAG_SUB(IXX,II+1,KZZ).EQ.0 .OR. &
                        MARCAP2(IXX,II+1,KZZ).NE.INDV) FLAGY2=1
                     IF (FLAGY2 .EQ. 0) THEN
                        FLAGY2=FLAGY2_GRID(IXX,II,KZZ)
                     ENDIF
                  ENDIF
               !-Y
                  IF(INICIOY(NVOID) .LE. LOW1+NSEC)  THEN
                     FLAGY1=1
                     FLAG2 = 1
                  ELSE
                     II=INICIOY(NVOID)
                     IF(FLAG_SUB(IXX,II-1,KZZ).EQ.0 .OR. &
                        MARCAP2(IXX,II-1,KZZ).NE.INDV) FLAGY1=1
                     IF (FLAGY1 .EQ. 0) THEN
                        FLAGY1=FLAGY1_GRID(IXX,II,KZZ)
                     ENDIF
                  ENDIF
            !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
               ENDIF
         ENDDO
         ENDDO OUTY
         
      !!!!!    expand along  Z
      !X!X!X!X!X!X!X!X!X!X!X!X!X!X!X!X!X!X!X!X!X!X!X!X
      !X!X!X!X!X!X!X!X!X!X!X!X!X!X!X!X!X!X!X!X!X!X!X!X
         OUTZ:DO JYY=INICIOY(NVOID), FINALY(NVOID)
         DO IXX=INICIOX(NVOID), FINALX(NVOID)
            !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-  
               IF(FLAGZ1==0 .OR. FLAGZ2==0) THEN
               !+Z
                  IF(FINALZ(NVOID) .GE. LOW2-NSEC)  THEN
                     FLAGZ2=1
                     FLAG2 = 1
                  ELSE
                     II=FINALZ(NVOID)
                     IF(FLAG_SUB(IXX,JYY,II+1).EQ.0 .OR. &
                        MARCAP2(IXX,JYY,II+1).NE.INDV) FLAGZ2=1
                     IF(FLAGZ2 .EQ. 0) THEN
                        FLAGZ2=FLAGZ2_GRID(IXX,JYY,II)
                     ENDIF
                  ENDIF
               !-Z
                  IF(INICIOZ(NVOID) .LE. LOW1+NSEC)  THEN
                     FLAGZ1=1
                     FLAG2 = 1
                  ELSE
                     II=INICIOZ(NVOID)
                     IF(FLAG_SUB(IXX,JYY,II-1).EQ.0 .OR. &
                        MARCAP2(IXX,JYY,II-1).NE.INDV) FLAGZ1=1
                     IF(FLAGZ1 .EQ. 0) THEN
                        FLAGZ1=FLAGZ1_GRID(IXX,JYY,II)
                     ENDIF
                  ENDIF
            !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
               ENDIF
         ENDDO
         ENDDO OUTZ
         
      !X!X!X!X!X!X!X!X!X!X!X!X!X!X!X!X!X!X!X!X!X!X!X!X
      !X!X!X!X!X!X!X!X!X!X!X!X!X!X!X!X!X!X!X!X!X!X!X!X

      ! !%%%  FLAGXYZ=0: face of the parallelepiped can be expanded --> move by one row
      ! IF(FLAGX1 == 0) INICIOX(NVOID)=INICIOX(NVOID)-1
      ! IF(FLAGX2 == 0) FINALX(NVOID)=FINALX(NVOID)+1
      ! IF(FLAGY1 == 0) INICIOY(NVOID)=INICIOY(NVOID)-1
      ! IF(FLAGY2 == 0) FINALY(NVOID)=FINALY(NVOID)+1
      ! IF(FLAGZ1 == 0) INICIOZ(NVOID)=INICIOZ(NVOID)-1
      ! IF(FLAGZ2 == 0) FINALZ(NVOID)=FINALZ(NVOID)+1

      !%%% EXPAND SQUARE-LIKE
      IF(FLAGX1 == 0 .AND. FLAGX2 == 0 .AND. &
         FLAGY1 == 0 .AND. FLAGY2 == 0 .AND. &
         FLAGZ1 == 0 .AND. FLAGZ2 == 0) THEN
         INICIOX(NVOID)=INICIOX(NVOID)-1
         FINALX(NVOID)=FINALX(NVOID)+1
         INICIOY(NVOID)=INICIOY(NVOID)-1
         FINALY(NVOID)=FINALY(NVOID)+1
         INICIOZ(NVOID)=INICIOZ(NVOID)-1
         FINALZ(NVOID)=FINALZ(NVOID)+1
      ELSE
         FLAGX1=1
         FLAGX2=1
         FLAGY1=1
         FLAGY2=1
         FLAGZ1=1
         FLAGZ2=1
      ENDIF

      !%%% if edges close to the box boundaries stop expansion
      IF(INICIOX(NVOID) .LE. LOW1+NSEC) FLAGX1=1 
      IF(FINALX(NVOID) .GE.  LOW2-NSEC) FLAGX2=1
      IF(INICIOY(NVOID) .LE. LOW1+NSEC) FLAGY1=1
      IF(FINALY(NVOID) .GE.  LOW2-NSEC) FLAGY2=1
      IF(INICIOZ(NVOID) .LE. LOW1+NSEC) FLAGZ1=1
      IF(FINALZ(NVOID) .GE.  LOW2-NSEC) FLAGZ2=1

      !%%% flag=1 in all the directions: no more expansion is possible
      IF(FLAGX1==1 .AND. FLAGX2==1 .AND. &
         FLAGY1==1 .AND. FLAGY2==1 .AND. &
         FLAGZ1==1 .AND. FLAGZ2==1) FLAG=1 !no faces to be expanded --> exit loop

   ENDDO !void expansion: do while (FLAG==0) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   IF(INICIOX(NVOID) .GE. FINALX(NVOID) .AND. &
      INICIOY(NVOID) .GE. FINALY(NVOID) .AND. &
      INICIOZ(NVOID) .GE. FINALZ(NVOID)) FLAG2=1 !exclude one-cell side voids

   IF(INICIOX(NVOID) .GT. FINALX(NVOID) .OR. &
      INICIOY(NVOID) .GT. FINALY(NVOID) .OR. &
      INICIOZ(NVOID) .GT. FINALZ(NVOID)) FLAG2=1 !include one-cell size voids
   
   !%%% Delete paralleliped voids with small sides to avoid percollation
   A=(FINALX(NVOID)-INICIOX(NVOID)+1)*DDX
   B=(FINALY(NVOID)-INICIOY(NVOID)+1)*DDY
   C=(FINALZ(NVOID)-INICIOZ(NVOID)+1)*DDZ
   VV=A*B*C

   A = MIN(A,B)
   A = MIN(A,C)
   IF(A .LT. SIDE_MIN) FLAG2=1

   !check if void is at the edge of the box
   IF(INICIOX(NVOID) .EQ. LOW1 .OR. FINALX(NVOID) .EQ. LOW2 .OR. &
      INICIOY(NVOID) .EQ. LOW1 .OR. FINALY(NVOID) .EQ. LOW2 .OR. &
      INICIOZ(NVOID) .EQ. LOW1 .OR. FINALZ(NVOID) .EQ. LOW2) FLAG2=1

   IF(FLAG2==0) THEN
      !%%% save void quantities
      ICX(NVOID)=IX
      ICY(NVOID)=JY
      ICZ(NVOID)=KZ
      !%%% save cube borders (not boundary cells centres)
      RINIXCO(NVOID)=RX1+(INICIOX(NVOID)-1)*DDX - 0.5*DDX 
      RFINXCO(NVOID)=RX1+(FINALX(NVOID)-1 )*DDX + 0.5*DDX
      RINIYCO(NVOID)=RY1+(INICIOY(NVOID)-1)*DDY - 0.5*DDY
      RFINYCO(NVOID)=RY1+(FINALY(NVOID)-1 )*DDY + 0.5*DDY
      RINIZCO(NVOID)=RZ1+(INICIOZ(NVOID)-1)*DDZ - 0.5*DDZ
      RFINZCO(NVOID)=RZ1+(FINALZ(NVOID)-1 )*DDZ + 0.5*DDZ
      VOL(NVOID)=VV
      MARCA_AUX(INICIOX(NVOID):FINALX(NVOID), &
                INICIOY(NVOID):FINALY(NVOID), &
                INICIOZ(NVOID):FINALZ(NVOID))=NVOID   
   ELSE
      NVOID=NVOID-1
   ENDIF

   ENDIF !if flag1=0
   !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ 

ENDDO !loop on center cells: L =1, NXYZ !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ 
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ 
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

DEALLOCATE(FLAGX1_GRID, FLAGX2_GRID, FLAGY1_GRID, FLAGY2_GRID, FLAGZ1_GRID, FLAGZ2_GRID)
DEALLOCATE(FLAG1_GRID)
!***************************************************************************
END SUBROUTINE VOIDFIND
!***************************************************************************



!*******************************************************************************************
SUBROUTINE MERGE_VOID(NVOID,INDICE,LOW1,LOW2,DXX,DYY,DZZ,VOLNEW,UVOID,GXC,GYC,GZC,NVOID2)
!******************************************************************************************* 
     USE COMMONDATA
     USE COSMOKDTREE
     IMPLICIT NONE
!input variables
     INTEGER:: NVOID
     REAL*4:: DXX, DYY, DZZ
     INTEGER:: LOW1,LOW2
!local variables
     INTEGER I, II, IND0, IND1, IND2, J
     INTEGER NBUFF
     INTEGER NVOID2
     INTEGER:: IX, JY, KZ
     INTEGER:: FLAGX, FLAGY, FLAGZ
     INTEGER, DIMENSION(NVOID) ::  INDICE
     REAL*4:: RX1, RX2, RY1, RY2, RZ1, RZ2
     INTEGER :: INDP, INDP2
     INTEGER, DIMENSION(NVOID):: MAJOR
!search
     INTEGER:: IX1, IX2, JY1, JY2, KZ1, KZ2, IX3, IX4, JY3, JY4, KZ3, KZ4
     INTEGER:: DIST,DIST2
!output variables
     REAL*4, DIMENSION(NVOID)::GXC,GYC,GZC !geometrical center of the void
     REAL*4 :: TREEPOINTS(NVOID,3)
     REAL*4, DIMENSION(NVOID):: VOLNEW
     INTEGER, DIMENSION(NVOID):: UVOID
!k-d tree
     REAL :: LPERIODIC(3)
     type(KDTreeNode), pointer :: TREE
     type(KDTreeResult) :: QUERY
     INTEGER :: CONTA
     REAL*4 :: TAR(3), RTHIS, RR

     !DX--> DDX
     !NCOX --> NXX
!*----------------------------------------------------------*
!*      cross-match all voids to find all the overlappings
!*----------------------------------------------------------*
     GXC(:)=0 !center of the void
     GYC(:)=0
     GZC(:)=0
     VOLNEW(:)=0. !new volume of the void

     !Geometrical centre of void and volume
     DO I=1,NVOID
         !define center
         GXC(I)=RADX(1)+(ICX(I)-1)*DXX
         GYC(I)=RADY(1)+(ICY(I)-1)*DYY
         GZC(I)=RADZ(1)+(ICZ(I)-1)*DZZ   
         !Volume
         RX1=RINIXCO(I)
         RX2=RFINXCO(I)
         RY1=RINIYCO(I)
         RY2=RFINYCO(I)
         RZ1=RINIZCO(I)
         RZ2=RFINZCO(I)
         VOLNEW(I)=(RX2-RX1)*(RY2-RY1)*(RZ2-RZ1)
     ENDDO
     
     !Build k-d tree for fast neighbour search
     TREEPOINTS(:,1)=GXC(:)
     TREEPOINTS(:,2)=GYC(:)
     TREEPOINTS(:,3)=GZC(:)

#if periodic == 1
      LPERIODIC = [LADO0PLUS, LADO0PLUS, LADO0PLUS]
      TREE => build_kdtree(TREEPOINTS,LPERIODIC)
#else
      TREE => build_kdtree(TREEPOINTS)
#endif 

     UVOID(:)=-1 !if 0, void has been merged or is unable to be merged, 
                 !if -1 it needs to be merged
                 !if >0, indicates to which void it has been merged
     NVOID2=NVOID !new number of voids
     MARCA = 0

     !Buffer for overlapping. If 1, takes into account voids next to each other.
     NBUFF = 1

     !If Major = 1 this void has become a major and does not need to be seen again
     MAJOR(:) = 0

     !%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%
     !%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%
     DO I=1,NVOID

        IND1=INDICE(I) !INDICE(1) --> largest void

        !----------------------------------------------------------
        ! Belongs to a major void
        IF(UVOID(IND1) .GT. 0) THEN
        !----------------------------------------------------------
            IND0=UVOID(IND1) !if void has been already merged, consider the void it has been merged to
            INDP=FATHER(IND0) !father of the void

        !----------------------------------------------------------
        !Able to become a major
        ELSE
        !----------------------------------------------------------

            IND0=IND1
            INDP=FATHER(IND0) !father of the void

            !CHECK IF, INCREASING THE VOID BY X CELLS IN EACH DIRECTION, THERE IS ALREADY A VOID
            !IF MORE THAN ONE, THE CLOSEST ONE WILL BE CONSIDERED
            !THE NUMBER 1 COMES FRON THE MEDIAN OF THE MINIMUM DISTANCE BETWEEN A VOID
            !AND THE NON-OVERLAPPING ONES THAT A BIGGER THAN IT. THIS NUMBER HAS BEEN OBTAINED
            !PERFORMING AN STUDY AND IT IS THE SAME FOR 64^3, 128^3 AND 256^3 BOXES

            IX1 = INICIOX(IND1) - 1 
            IX2 = FINALX(IND1)  + 1 
            JY1 = INICIOY(IND1) - 1 
            JY2 = FINALY(IND1)  + 1 
            KZ1 = INICIOZ(IND1) - 1 
            KZ2 = FINALZ(IND1)  + 1
            IX1=MAX(LOW1,IX1)
            IX2=MIN(LOW2,IX2)
            JY1=MAX(LOW1,JY1)
            JY2=MIN(LOW2,JY2)
            KZ1=MAX(LOW1,KZ1)
            KZ2=MIN(LOW2,KZ2)

            DIST=HUGE(0)
            DIST2 = 0
            DO KZ=KZ1, KZ2
               DO JY=JY1, JY2
                  DO IX=IX1, IX2
                     IF(MARCA(IX,JY,KZ) .GT. 0) THEN
                        !ONLY CONSIDER VOIDS THAT BELONG TO THE SAME FATHER
                        INDP2 = FATHER(MARCA(IX,JY,KZ))
                        IF(INDP2 .NE. INDP) CYCLE  
                        !IF THERE IS A CLOSE BROTHER, YOU BELONG TO HIM
                        DIST2 =(ICX(IND1)-IX)**2 + (ICY(IND1)-JY)**2 + (ICZ(IND1)-KZ)**2
                        IF( DIST2 .LT. DIST) THEN
                           IND0=MARCA(IX,JY,KZ)
                           DIST2 = DIST
                        ENDIF 
                     ENDIF
                  ENDDO
               ENDDO
            ENDDO

            !IF ANOTHER VOID, DEFINE OWNERSHIP
            IF(IND0 .NE. IND1) THEN
               UVOID(IND1)=IND0
               VOLNEW(IND0)=VOLNEW(IND0)+VOLNEW(IND1)

               !cells belong IND0
               DO KZ=INICIOZ(IND1), FINALZ(IND1)
                  DO JY=INICIOY(IND1), FINALY(IND1)
                     DO IX=INICIOX(IND1), FINALX(IND1)
                        IF(MARCA(IX,JY,KZ) .EQ. 0) MARCA(IX,JY, KZ)=IND0
                     ENDDO
                  ENDDO
               ENDDO

               !redefine center of IND0 with Volume Weighted Center
               GXC(IND0) = (GXC(IND0)*VOLNEW(IND0) + GXC(IND1)*VOLNEW(IND1))/(VOLNEW(IND0)+VOLNEW(IND1))
               GYC(IND0) = (GYC(IND0)*VOLNEW(IND0) + GYC(IND1)*VOLNEW(IND1))/(VOLNEW(IND0)+VOLNEW(IND1))
               GZC(IND0) = (GZC(IND0)*VOLNEW(IND0) + GZC(IND1)*VOLNEW(IND1))/(VOLNEW(IND0)+VOLNEW(IND1))

               NVOID2=NVOID2-1

            !IT IS A MAJOR VOID
            ELSE

               !cells belong IND0
               DO KZ=INICIOZ(IND0), FINALZ(IND0)
                  DO JY=INICIOY(IND0), FINALY(IND0)
                     DO IX=INICIOX(IND0), FINALX(IND0)
                        IF(MARCA(IX,JY,KZ) .EQ. 0) MARCA(IX,JY, KZ)=IND0
                     ENDDO
                  ENDDO
               ENDDO

               MAJOR(IND0) = 1

            ENDIF

        !----------------------------------------------------------
        ENDIF
        !----------------------------------------------------------

        !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
        !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

        ! Cells belonging to this cube, plus NBUFF expanding to detect 
        ! voids touching it but not overlapping
        IX1 = INICIOX(IND1) - NBUFF
        IX2 = FINALX(IND1)  + NBUFF
        JY1 = INICIOY(IND1) - NBUFF
        JY2 = FINALY(IND1)  + NBUFF
        KZ1 = INICIOZ(IND1) - NBUFF
        KZ2 = FINALZ(IND1)  + NBUFF

        !Nearby voids that can overlap
        TAR(1) = GXC(IND1)
        TAR(2) = GYC(IND1)
        TAR(3) = GZC(IND1) 

        !Sphere containing this cube
        RTHIS = SQRT(3.)/2.*VOLNEW(IND1)**(1./3.)

        !query radius
        RR = (2*RTHIS) + 2.1*SQRT(3.)*DXX*NBUFF

        QUERY = ball_search(TREE, TAR, RR, .true.) 
        CONTA = SIZE(QUERY%idx)

         ! DO J=I+1,NVOID
         !   IND2 = INDICE(J)        
         DO II=1,CONTA
           IND2=QUERY%idx(II)

           IF(IND2.EQ.IND1) CYCLE !it is the same void
           IF(MAJOR(IND2) .EQ. 1) CYCLE !it is a major void
           IF(UVOID(IND2) .GE. 0) CYCLE !it belongs to another void
           IF(INDP .NE. FATHER(IND2)) CYCLE !only consider voids with the same father

           FLAGX=0
           FLAGY=0
           FLAGZ=0

           IX3 = INICIOX(IND2) - NBUFF
           IX4 = FINALX(IND2)  + NBUFF
           JY3 = INICIOY(IND2) - NBUFF
           JY4 = FINALY(IND2)  + NBUFF
           KZ3 = INICIOZ(IND2) - NBUFF
           KZ4 = FINALZ(IND2)  + NBUFF

           !Check if the voids overlap with integer coordinates
           IF (IX1.LE.IX4 .AND. IX3.LE.IX2) FLAGX=1
           IF (JY1.LE.JY4 .AND. JY3.LE.JY2) FLAGY=1
           IF (KZ1.LE.KZ4 .AND. KZ3.LE.KZ2) FLAGZ=1

          !***************************************************
          IF(FLAGX==1 .AND. FLAGY==1 .AND. FLAGZ==1) THEN !overlapping

            UVOID(IND2)=IND0

            !update volume
            VOLNEW(IND0)=VOLNEW(IND0)+VOLNEW(IND2)
            
            !cells belong now to IND0
            DO KZ=INICIOZ(IND2), FINALZ(IND2)
               DO JY=INICIOY(IND2), FINALY(IND2)
                  DO IX=INICIOX(IND2), FINALX(IND2)
                     !marca is only updated in marca = 0 cells
                     !thus, if a cell already belongs to a void, it will not be updated
                     IF(MARCA(IX,JY,KZ) .EQ. 0) MARCA(IX,JY, KZ)=IND0
                  ENDDO
               ENDDO
            ENDDO

            !redefine center of IND0 with Volume Weighted Center
            GXC(IND0) = (GXC(IND0)*VOLNEW(IND0) + GXC(IND2)*VOLNEW(IND2))/(VOLNEW(IND0)+VOLNEW(IND2))
            GYC(IND0) = (GYC(IND0)*VOLNEW(IND0) + GYC(IND2)*VOLNEW(IND2))/(VOLNEW(IND0)+VOLNEW(IND2))
            GZC(IND0) = (GZC(IND0)*VOLNEW(IND0) + GZC(IND2)*VOLNEW(IND2))/(VOLNEW(IND0)+VOLNEW(IND2))

            !do not consider IND2 anymore
            NVOID2=NVOID2-1

          ENDIF
          !***************************************************

        ENDDO
        !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
        !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

      ENDDO
      !%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%
      !%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%

      !Deallocate k-d tree
      call deallocate_kdtree(TREE)

!***************************************************************************
END SUBROUTINE MERGE_VOID
!***************************************************************************





!***************************************************************************
SUBROUTINE SIMPLE_CONNECTION_FAST(NVOID,MARCA,LOW1,LOW2)
!Breadth-First Search [O(Ncell)] algorithm to find the connected area outside
!a void. Everything that is not connected to this area and does
!not belong to the void, is a hole that must be filled.
!***************************************************************************

   IMPLICIT NONE
   !input variables
   INTEGER :: NVOID, LOW1, LOW2
   INTEGER, DIMENSION(NVOID) :: INDICE, UVOID
   !local
   INTEGER IX,JY,KZ,IND,I,II
   INTEGER :: COUNTER
   INTEGER, ALLOCATABLE, DIMENSION(:,:,:) :: AUX_GRID
   !BREADTH-FIRST SEARCH
   INTEGER :: MAX_QUEUE_SIZE
   INTEGER, ALLOCATABLE, DIMENSION(:,:) :: QUEUE
   INTEGER :: HEAD, TAIL
   !output variable
   INTEGER :: MARCA(LOW1:LOW2,LOW1:LOW2,LOW1:LOW2)
   
   WRITE(*,*) 'Breaddth-First Search algorithm starting...'

   !allocate queue
   MAX_QUEUE_SIZE = (LOW2-LOW1+1)**3
   ALLOCATE(QUEUE(3,MAX_QUEUE_SIZE))

   !GRID WHERE BFS IS PERFORMED
   ALLOCATE(AUX_GRID(LOW1:LOW2,LOW1:LOW2,LOW1:LOW2))
   AUX_GRID(:,:,:) = 0
   WHERE (MARCA .GT. 0) AUX_GRID = 1

   !Nothing in the borders can be a hole or belong to a void: outside

   !BREADTH-FIRST SEARCH algorithm
   !Starting from the box border, that we know does not belong to any void or hole,
   !find the outside connected area filling with 2s

   HEAD = 1
   TAIL = 1
   QUEUE(1,HEAD) = LOW1
   QUEUE(2,HEAD) = LOW1
   QUEUE(3,HEAD) = LOW1
   AUX_GRID(LOW1,LOW1,LOW1) = 2

   !+++++++++++++++++++++++++++
   DO WHILE (HEAD .LE. TAIL)
   !+++++++++++++++++++++++++++

      !CHECK NEXT CELL IN THE QUEUE
      IX = QUEUE(1,HEAD)
      JY = QUEUE(2,HEAD)
      KZ = QUEUE(3,HEAD)
      HEAD = HEAD + 1

      !PROCESS NEIGHBOURS IN 6 POSSIBLE DIRECTIONS
      !-X
      II = IX-1
      IF (II .GE. LOW1) THEN 
      IF (AUX_GRID(II,JY,KZ) .EQ. 0) THEN
         AUX_GRID(II,JY,KZ) = 2
         TAIL = TAIL + 1
         QUEUE(1,TAIL) = II
         QUEUE(2,TAIL) = JY
         QUEUE(3,TAIL) = KZ
      ENDIF
      ENDIF

      !+X
      II = IX+1
      IF (II .LE. LOW2) THEN
      IF (AUX_GRID(II,JY,KZ) .EQ. 0) THEN
         AUX_GRID(II,JY,KZ) = 2
         TAIL = TAIL + 1
         QUEUE(1,TAIL) = II
         QUEUE(2,TAIL) = JY
         QUEUE(3,TAIL) = KZ
      ENDIF
      ENDIF

      !-Y
      II = JY-1
      IF (II .GE. LOW1) THEN
      IF (AUX_GRID(IX,II,KZ) .EQ. 0) THEN
         AUX_GRID(IX,II,KZ) = 2
         TAIL = TAIL + 1
         QUEUE(1,TAIL) = IX
         QUEUE(2,TAIL) = II
         QUEUE(3,TAIL) = KZ
      ENDIF
      ENDIF

      !+Y
      II = JY+1
      IF (II .LE. LOW2) THEN
      IF (AUX_GRID(IX,II,KZ) .EQ. 0) THEN
         AUX_GRID(IX,II,KZ) = 2
         TAIL = TAIL + 1
         QUEUE(1,TAIL) = IX
         QUEUE(2,TAIL) = II
         QUEUE(3,TAIL) = KZ
      ENDIF
      ENDIF

      !-Z
      II = KZ-1
      IF (II .GE. LOW1) THEN
      IF (AUX_GRID(IX,JY,II) .EQ. 0) THEN
         AUX_GRID(IX,JY,II) = 2
         TAIL = TAIL + 1
         QUEUE(1,TAIL) = IX
         QUEUE(2,TAIL) = JY
         QUEUE(3,TAIL) = II
      ENDIF
      ENDIF

      !+Z
      II = KZ+1
      IF (II .LE. LOW2) THEN
      IF (AUX_GRID(IX,JY,II) .EQ. 0) THEN
         AUX_GRID(IX,JY,II) = 2
         TAIL = TAIL + 1
         QUEUE(1,TAIL) = IX
         QUEUE(2,TAIL) = JY
         QUEUE(3,TAIL) = II
      ENDIF
      ENDIF

   !+++++++++++++++++++++++++++
   ENDDO !while
   !+++++++++++++++++++++++++++

   !NOW, EVERYTHING AFTER THE BREADTH-FIRST SEARCH THAT HAS AUX_GRID=0
   !IS A HOLE AND HENCE MUST BE FILLED
   DO KZ=LOW1,LOW2
      DO JY=LOW1,LOW2
         DO IX=LOW1,LOW2
            !IF IT IS A HOLE
            IF(AUX_GRID(IX,JY,KZ) .GT. 0) CYCLE
            !FILL HOLE FINDING ENCLOSING VOID
            !+X
            II = IX
            DO WHILE(II .LE. LOW2 .AND. AUX_GRID(II,JY,KZ) .EQ. 0)
               II = II + 1
            ENDDO
            IND = MARCA(II,JY,KZ)
            MARCA(IX,JY,KZ) = IND
         
         ENDDO
      ENDDO
   ENDDO

   DEALLOCATE(AUX_GRID)

!***************************************************************************
END SUBROUTINE SIMPLE_CONNECTION_FAST
!***************************************************************************





! !***************************************************************************
! SUBROUTINE SIMPLE_CONNECTION(NVOID,INDICE,UVOID,MARCA,LOW1,LOW2)
! !Breadth-First Search [O(Ncell)] algorithm to find the connected area outside
! !a void. Everything that is not connected to this area and does
! !not belong to the void, is a hole that must be filled.
! !***************************************************************************

!    IMPLICIT NONE
!    !input variables
!    INTEGER :: NVOID, LOW1, LOW2
!    INTEGER, DIMENSION(NVOID) :: INDICE, UVOID
!    !local
!    INTEGER IX,JY,KZ,IND,I,II
!    INTEGER IX1,JY1,KZ1,IX2,JY2,KZ2
!    INTEGER NXVOID, NYVOID, NZVOID
!    INTEGER :: COUNTER
!    INTEGER, ALLOCATABLE, DIMENSION(:,:,:) :: AUX_GRID
!    !BREADTH-FIRST SEARCH
!    INTEGER :: MAX_QUEUE_SIZE
!    INTEGER, ALLOCATABLE, DIMENSION(:,:) :: QUEUE
!    INTEGER :: HEAD, TAIL
!    !output variable
!    INTEGER :: MARCA(LOW1:LOW2,LOW1:LOW2,LOW1:LOW2)
   
!    WRITE(*,*) 'Breaddth-First Search algorithm starting...'

!    !allocate queue
!    MAX_QUEUE_SIZE = (LOW2-LOW1+1)**3
!    ALLOCATE(QUEUE(3,MAX_QUEUE_SIZE))


!    DO I=1,NVOID
!    !-----------------------
!       IND = INDICE(I)
!       IF(UVOID(IND) .GE. 0) CYCLE !ONLY MAIN VOIDS

!       !GET MINIMUM AND MAXIMUM IX, JY, KZ WITH MARCA(IX,JY,KZ)=IND
!       IX1=LOW2+1
!       IX2=LOW1-1
!       JY1=LOW2+1
!       JY2=LOW1-1
!       KZ1=LOW2+1
!       KZ2=LOW1-1
!       DO KZ=LOW1,LOW2
!          DO JY=LOW1,LOW2
!             DO IX=LOW1,LOW2
!                IF(MARCA(IX,JY,KZ) .EQ. IND) THEN
!                   IF(IX .LT. IX1) IX1=IX
!                   IF(IX .GT. IX2) IX2=IX
!                   IF(JY .LT. JY1) JY1=JY
!                   IF(JY .GT. JY2) JY2=JY
!                   IF(KZ .LT. KZ1) KZ1=KZ
!                   IF(KZ .GT. KZ2) KZ2=KZ
!                ENDIF
!             ENDDO
!          ENDDO
!       ENDDO
      
!       IF(IX2.LT.IX1 .OR. JY2.LT.JY1 .OR. KZ2.LT.KZ1) THEN
!          !VOID HAS NO CELLS LEFT, IT HAS BEEN FILLED BY OTHER VOID
!          UVOID(IND) = 0
!          CYCLE
!       ENDIF

!       NXVOID = IX2-IX1+1
!       NYVOID = JY2-JY1+1
!       NZVOID = KZ2-KZ1+1

!       !Consider a buffer of 1 cell in each direction
!       ALLOCATE(AUX_GRID(0:NXVOID+1,0:NYVOID+1,0:NZVOID+1))
!       AUX_GRID(:,:,:) = 0
!       !WHERE MARCA .EQ. IND, AUX_GRID = 1
!       DO KZ=KZ1,KZ2
!          DO JY=JY1,JY2
!             DO IX=IX1,IX2
!                IF(MARCA(IX,JY,KZ) .EQ. IND) AUX_GRID(IX-IX1+1,JY-JY1+1,KZ-KZ1+1) = 1
!             ENDDO
!          ENDDO
!       ENDDO

!       !BREADTH-FIRST SEARCH algorithm
!       !Starting from the void border, that we know does not belong to the void
!       !find the outside connected area filling with 2s
!       HEAD = 1
!       TAIL = 1
!       QUEUE(1,HEAD) = 0
!       QUEUE(2,HEAD) = 0
!       QUEUE(3,HEAD) = 0
!       AUX_GRID(0,0,0) = 2

!       !+++++++++++++++++++++++++++
!       DO WHILE (HEAD .LE. TAIL)
!       !+++++++++++++++++++++++++++

!          !CHECK NEXT CELL IN THE QUEUE
!          IX = QUEUE(1,HEAD)
!          JY = QUEUE(2,HEAD)
!          KZ = QUEUE(3,HEAD)
!          HEAD = HEAD + 1

!          !PROCESS NEIGHBOURS IN 6 POSSIBLE DIRECTIONS
!          !-X
!          II = IX-1
!          IF (II .GE. 0) THEN 
!          IF (AUX_GRID(II,JY,KZ) .EQ. 0) THEN
!             AUX_GRID(II,JY,KZ) = 2
!             TAIL = TAIL + 1
!             QUEUE(1,TAIL) = II
!             QUEUE(2,TAIL) = JY
!             QUEUE(3,TAIL) = KZ
!          ENDIF
!          ENDIF

!          !+X
!          II = IX+1
!          IF (II .LE. NXVOID+1) THEN
!          IF (AUX_GRID(II,JY,KZ) .EQ. 0) THEN
!             AUX_GRID(II,JY,KZ) = 2
!             TAIL = TAIL + 1
!             QUEUE(1,TAIL) = II
!             QUEUE(2,TAIL) = JY
!             QUEUE(3,TAIL) = KZ
!          ENDIF
!          ENDIF

!          !-Y
!          II = JY-1
!          IF (II .GE. 0) THEN
!          IF (AUX_GRID(IX,II,KZ) .EQ. 0) THEN
!             AUX_GRID(IX,II,KZ) = 2
!             TAIL = TAIL + 1
!             QUEUE(1,TAIL) = IX
!             QUEUE(2,TAIL) = II
!             QUEUE(3,TAIL) = KZ
!          ENDIF
!          ENDIF

!          !+Y
!          II = JY+1
!          IF (II .LE. NYVOID+1) THEN
!          IF (AUX_GRID(IX,II,KZ) .EQ. 0) THEN
!             AUX_GRID(IX,II,KZ) = 2
!             TAIL = TAIL + 1
!             QUEUE(1,TAIL) = IX
!             QUEUE(2,TAIL) = II
!             QUEUE(3,TAIL) = KZ
!          ENDIF
!          ENDIF

!          !-Z
!          II = KZ-1
!          IF (II .GE. 0) THEN
!          IF (AUX_GRID(IX,JY,II) .EQ. 0) THEN
!             AUX_GRID(IX,JY,II) = 2
!             TAIL = TAIL + 1
!             QUEUE(1,TAIL) = IX
!             QUEUE(2,TAIL) = JY
!             QUEUE(3,TAIL) = II
!          ENDIF
!          ENDIF

!          !+Z
!          II = KZ+1
!          IF (II .LE. NZVOID+1) THEN
!          IF (AUX_GRID(IX,JY,II) .EQ. 0) THEN
!             AUX_GRID(IX,JY,II) = 2
!             TAIL = TAIL + 1
!             QUEUE(1,TAIL) = IX
!             QUEUE(2,TAIL) = JY
!             QUEUE(3,TAIL) = II
!          ENDIF
!          ENDIF

!       !+++++++++++++++++++++++++++
!       ENDDO !while
!       !+++++++++++++++++++++++++++

!       !NOW, EVERYTHING AFTER THE BREADTH-FIRST SEARCH THAT HAS AUX_GRID=0
!       !IS A HOLE AND HENCE MUST BE FILLED
!       DO KZ=KZ1,KZ2
!          DO JY=JY1,JY2
!             DO IX=IX1,IX2
!                IF(AUX_GRID(IX-IX1+1,JY-JY1+1,KZ-KZ1+1) .EQ. 0) MARCA(IX,JY,KZ) = IND
!             ENDDO
!          ENDDO
!       ENDDO

!       DEALLOCATE(AUX_GRID)

!    !-----------------------
!    ENDDO
!    !-----------------------

! !***************************************************************************
! END SUBROUTINE SIMPLE_CONNECTION
! !***************************************************************************


!***************************************************************************
SUBROUTINE VOID_PERIODIC(NVOID,INDICE,UVOID,GXC,GYC,GZC,VOLNEW, &
                        NX,NY,NZ,DX,DY,DZ,LOW1,LOW2)
!Checks voids having cells outside the box and finds the periodic images
!***************************************************************************
   USE COMMONDATA, ONLY: LADO0, LADO0PLUS, MARCA, PI, RADX, RADY, RADZ
   IMPLICIT NONE
   !input variables
   INTEGER :: NVOID, NX, NY, NZ, LOW1, LOW2
   INTEGER, DIMENSION(NVOID) :: INDICE, UVOID, INDICE2
   REAL, DIMENSION(NVOID) :: GXC, GYC, GZC, VOLNEW
   REAL :: DX, DY, DZ
   !local
   INTEGER :: COUNTER, COUNTER2
   INTEGER :: IND,I,J,IND2,IX,JY,KZ,FLAG,IXX,JYY,KZZ
   INTEGER :: IGCX,IGCY,IGCZ
   REAL :: XC1,YC1,ZC1
   REAL :: INTERSEC
   ! REAL :: TOL, RR
   REAL, ALLOCATABLE :: LCX(:),LCY(:),LCZ(:)
   INTEGER, ALLOCATABLE :: LINDEXING(:)
   INTEGER :: NLOCAL


   !FIND MAJOR VOIDS INSIDE OR INTERSECTING THE BOX
   INDICE2 = 0
   NLOCAL = 0
   DO I=1,NVOID
      IND = INDICE(I)
      IF (UVOID(IND) .GE. 0) CYCLE !ONLY MAIN VOIDS

      NLOCAL = NLOCAL + 1
      INDICE2(NLOCAL) = IND

      !CHECK VOID HAS COMPONENTS INSIDE THE BOX
      FLAG = 0
      DO KZ=1,NZ
         DO JY=1,NY
            DO IX=1,NX
               IF(MARCA(IX,JY,KZ) .EQ. IND) THEN
                  FLAG = 1
                  EXIT
               ENDIF
            ENDDO
         ENDDO
      ENDDO
      
      !DELETE (COMPLETELY) OUTSIDE VOIDS
      IF(FLAG.EQ.0) THEN
         VOLNEW(IND) = 0.
         UVOID(IND) = 0
         WHERE(UVOID .EQ. IND)
            UVOID = 0
         END WHERE
         WHERE (MARCA .EQ. IND)
            MARCA = 0
         END WHERE

      !SAVE VOIDS INSIDE OR INTERSECTING THE BOX
      ELSE
         NLOCAL = NLOCAL + 1
         INDICE2(NLOCAL) = IND
      ENDIF
   ENDDO



   ALLOCATE(LCX(NLOCAL),LCY(NLOCAL),LCZ(NLOCAL),LINDEXING(NLOCAL))
   DO I=1,NLOCAL
      IND = INDICE2(I)
      LINDEXING(I) = IND
      LCX(I) = GXC(IND)
      LCY(I) = GYC(IND)
      LCZ(I) = GZC(IND)
   ENDDO



   DO J=1,NLOCAL
      IND = LINDEXING(J)

      XC1 = LCX(J)
      YC1 = LCY(J)
      ZC1 = LCZ(J)

      !CHECK WHAT IS THE PERIODICITY OF THE VOID
      IX = 0
      IF(XC1 .LT. -LADO0/2) IX = +1
      IF(XC1 .GT. LADO0/2) IX = -1
      JY = 0
      IF(YC1 .LT. -LADO0/2) JY = +1
      IF(YC1 .GT. LADO0/2) JY = -1
      KZ = 0
      IF(ZC1 .LT. -LADO0/2) KZ = +1
      IF(ZC1 .GT. LADO0/2) KZ = -1

      !Cycle inner voids
      IF(IX.EQ.0 .AND. JY.EQ.0 .AND. KZ.EQ.0) CYCLE

      !For outer voids, calculate the intersection of their cells with the
      !corresponding (bigger) inner void
      
      !1-Identify where the center lies
      IGCX = INT((GXC(IND)-RADX(1))/DX) + 1
      IGCY = INT((GYC(IND)-RADY(1))/DY) + 1
      IGCZ = INT((GZC(IND)-RADZ(1))/DZ) + 1

      !2-Inner void containing the periodically shifted center
      IND2 = MARCA(IGCX+IX*NX,IGCY+JY*NY,IGCZ+KZ*NZ)
      !This should not happen, but just in case
      IF (IND2 .EQ. 0) CYCLE

      !Loop over all cells belonging to IND, shifting them periodically
      !and check the fraction of IND belonging to IND2
      COUNTER = 0
      COUNTER2 = 0
      DO KZZ=LOW1,LOW2
         DO JYY=LOW1,LOW2
            DO IXX=LOW1,LOW2
               IF(MARCA(IXX,JYY,KZZ) .EQ. IND) THEN
                  !Bounds check
                  IF(IXX+IX*NX.GE.LOW1 .AND. IXX+IX*NX.LE.LOW2 .AND. &
                     JYY+JY*NY.GE.LOW1 .AND. JYY+JY*NY.LE.LOW2 .AND. &
                     KZZ+KZ*NZ.GE.LOW1 .AND. KZZ+KZ*NZ.LE.LOW2) THEN
                     !Check if cell belongs to IND2
                     COUNTER2 = COUNTER2 + 1
                     IF(MARCA(IXX+IX*NX,JYY+JY*NY,KZZ+KZ*NZ) .EQ. IND2) THEN
                        COUNTER = COUNTER + 1
                     ENDIF
                  ENDIF
               ENDIF
            ENDDO
         ENDDO
      ENDDO

      !Fraction of IND belonging to IND2
      IF (COUNTER2 .EQ. 0) THEN
         INTERSEC = 0
      ELSE
         INTERSEC = REAL(COUNTER)/REAL(COUNTER2)
      ENDIF

      !If the fraction is larger than X, then is the same void
      IF(INTERSEC .GT. 0.4) CALL CHANGE_BELONG(IND2,IND)

   ENDDO


   ! !CLEAN VOIDS OUTSIDE THAT HAVEN'T BEEN MERGED
   ! !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! DO J=1,NLOCAL
   !    IND = LINDEXING(J)

   !    IF (UVOID(IND) .GE. 0) CYCLE

   !    XC1 = LCX(J)
   !    YC1 = LCY(J)
   !    ZC1 = LCZ(J)

   !    FLAG = 0
   !    !Check if void lies outside the box
   !    IF(XC1 .LT. -LADO0/2. .OR. XC1 .GT. LADO0/2.) FLAG = 1
   !    IF(YC1 .LT. -LADO0/2. .OR. YC1 .GT. LADO0/2.) FLAG = 1
   !    IF(ZC1 .LT. -LADO0/2. .OR. ZC1 .GT. LADO0/2.) FLAG = 1
   !    IF(FLAG.EQ.0) CYCLE

   !    !CLEAN
   !    UVOID(IND) = 0
   !    VOLNEW(IND) = 0.
   !    WHERE (MARCA .EQ. IND)
   !       MARCA = 0
   !    END WHERE

   !    WHERE(UVOID .EQ. IND)
   !       UVOID = 0
   !    END WHERE

   ! ENDDO
   ! !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   DEALLOCATE(LCX,LCY,LCZ,LINDEXING)

   CONTAINS

      SUBROUTINE CHANGE_BELONG(I,II)
         IMPLICIT NONE
         INTEGER :: I,II

         UVOID(II) = I
         WHERE (UVOID .EQ. II)
            UVOID = I
         END WHERE

         WHERE (MARCA .EQ. II)
            MARCA = I
         END WHERE

         VOLNEW(II) = 0.

      END SUBROUTINE CHANGE_BELONG

!***************************************************************************
END SUBROUTINE VOID_PERIODIC
!***************************************************************************


!***************************************************************************
SUBROUTINE VOID_DECONSTRUCTION(NTH,NVOID,INDICE,UVOID,MARCA, &
                                  NX,NY,NZ,DX,DY,DZ,RX1,RY1,RZ1, &
                                  dNVOID,dUVOID, &
                                  dINICIOX,dINICIOY,dINICIOZ, &
                                  dFINALX,dFINALY,dFINALZ, &
                                  dRINIXCO,dRINIYCO,dRINIZCO, &
                                  dRFINXCO,dRFINYCO,dRFINZCO)
!**************************************************************************
!IT IS A SIMPLER VERSION OF VOIDFIND TO FIND A SET OF NON-OVERLAPPING
!PARALELEPIPEDS THAT DECONSTRUCT ALL VOIDS
!**************************************************************************
   USE COMMONDATA, ONLY: NVOID_MAX
   IMPLICIT NONE
   !input variables
   INTEGER :: NVOID, NX, NY, NZ
   INTEGER, DIMENSION(NVOID) :: INDICE, UVOID
   INTEGER, DIMENSION(NX,NY,NZ) :: MARCA
   REAL :: RX1, RY1, RZ1, DX, DY, DZ
   !local
   INTEGER, ALLOCATABLE, DIMENSION(:,:,:) :: MARCA_AUX
   INTEGER :: NXVOID, NYVOID, NZVOID
   INTEGER :: IXX,JYY,KZZ,IND,I,II,IX1,IX2,JY1,JY2,KZ1,KZ2,IX,JY,KZ
   INTEGER :: INIX,INIY,INIZ,IFIX,IFIY,IFIZ
   INTEGER :: FLAG,FLAGX1,FLAGX2,FLAGY1,FLAGY2,FLAGZ1,FLAGZ2
   REAL :: DIST, DIST2
   !output
   INTEGER :: dNVOID
   INTEGER, DIMENSION(NVOID_MAX) :: dUVOID
   INTEGER, DIMENSION(NVOID_MAX) :: dINICIOX,dINICIOY,dINICIOZ
   INTEGER, DIMENSION(NVOID_MAX) :: dFINALX,dFINALY,dFINALZ
   REAL, DIMENSION(NVOID_MAX) :: dRINIXCO,dRINIYCO,dRINIZCO
   REAL, DIMENSION(NVOID_MAX) :: dRFINXCO,dRFINYCO,dRFINZCO
   !explicit parallelization
   INTEGER :: ID, NTH, OMP_GET_THREAD_NUM
   INTEGER, ALLOCATABLE, DIMENSION(:) :: N2
   INTEGER, ALLOCATABLE, DIMENSION(:,:) :: UVOIDPAR
   INTEGER, ALLOCATABLE, DIMENSION(:,:) :: INICIOX,INICIOY,INICIOZ
   INTEGER, ALLOCATABLE, DIMENSION(:,:) :: FINALX,FINALY,FINALZ
   REAL, ALLOCATABLE, DIMENSION(:,:) :: RINIXCO,RINIYCO,RINIZCO
   REAL, ALLOCATABLE, DIMENSION(:,:) :: RFINXCO,RFINYCO,RFINZCO

   !$$$$$$$$$$$$$$$$
   !$PARALLELIZE$
   ALLOCATE(N2(NTH))
   ALLOCATE(UVOIDPAR(NTH,NVOID_MAX))
   ALLOCATE(INICIOX(NTH,NVOID_MAX),INICIOY(NTH,NVOID_MAX),INICIOZ(NTH,NVOID_MAX))
   ALLOCATE(FINALX(NTH,NVOID_MAX),FINALY(NTH,NVOID_MAX),FINALZ(NTH,NVOID_MAX))
   ALLOCATE(RINIXCO(NTH,NVOID_MAX),RINIYCO(NTH,NVOID_MAX),RINIZCO(NTH,NVOID_MAX))
   ALLOCATE(RFINXCO(NTH,NVOID_MAX),RFINYCO(NTH,NVOID_MAX),RFINZCO(NTH,NVOID_MAX))
   !$$$$$$$$$$$$$$$$

   !TOTAL NUMBER OF PARALELEPIPEDS DECONSTRUCTING ALL VOIDS
   N2(:) = 0 !dNVOID for each thread
  
   !INITIALIZE ARRAYS
   !Before reduction
   N2(:) = 0
   UVOIDPAR(:,:) = 0
   INICIOX(:,:) = 0
   INICIOY(:,:) = 0
   INICIOZ(:,:) = 0
   FINALX(:,:) = 0
   FINALY(:,:) = 0
   FINALZ(:,:) = 0
   RINIXCO(:,:) = 0
   RINIYCO(:,:) = 0
   RINIZCO(:,:) = 0
   RFINXCO(:,:) = 0
   RFINYCO(:,:) = 0
   RFINZCO(:,:) = 0

   !Final
   dNVOID = 0
   !INITIALIZE ARRAYS
   dUVOID(:) = 0
   dINICIOX(:) = 0
   dINICIOY(:) = 0
   dINICIOZ(:) = 0
   dFINALX(:) = 0
   dFINALY(:) = 0
   dFINALZ(:) = 0
   dRINIXCO(:) = 0
   dRINIYCO(:) = 0
   dRINIZCO(:) = 0
   dRFINXCO(:) = 0
   dRFINYCO(:) = 0
   dRFINZCO(:) = 0

   !$OMP PARALLEL DO SHARED(NTH,INDICE,UVOID,NX,NY,NZ,MARCA, &
   !$OMP DX,DY,DZ,RX1,RY1,RZ1,NVOID, &
   !$OMP N2,UVOIDPAR,INICIOX,INICIOY,INICIOZ,FINALX,FINALY,FINALZ, &
   !$OMP RINIXCO,RINIYCO,RINIZCO,RFINXCO,RFINYCO,RFINZCO), &
   !$OMP PRIVATE(I,ID,IND,IX1,IX2,JY1,JY2,KZ1,KZ2,IX,JY,KZ), &
   !$OMP PRIVATE(NXVOID,NYVOID,NZVOID,MARCA_AUX), &
   !$OMP PRIVATE(DIST,DIST2,II,IXX,JYY,KZZ,FLAG), &
   !$OMP PRIVATE(INIX,INIY,INIZ,IFIX,IFIY,IFIZ), &
   !$OMP PRIVATE(FLAGX1,FLAGX2,FLAGY1,FLAGY2,FLAGZ1,FLAGZ2), &
   !$OMP DEFAULT(NONE), SCHEDULE(DYNAMIC)
   !-----------------------
   DO I=1,NVOID
   !-----------------------
      !$$$$$$$$$$$$
      ID = OMP_GET_THREAD_NUM() + 1
      !$$$$$$$$$$$$

      IND = INDICE(I)
      IF(UVOID(IND) .GE. 0) CYCLE !ONLY MAIN VOIDS

      !GET MINIMUM AND MAXIMUM IX, JY, KZ WITH MARCA(IX,JY,KZ)=IND
      IX1=NX
      IX2=1
      JY1=NY
      JY2=1
      KZ1=NZ
      KZ2=1
      DO KZ=1,NZ
         DO JY=1,NY
            DO IX=1,NX
               IF(MARCA(IX,JY,KZ) .EQ. IND) THEN
                  IF(IX .LT. IX1) IX1=IX
                  IF(IX .GT. IX2) IX2=IX
                  IF(JY .LT. JY1) JY1=JY
                  IF(JY .GT. JY2) JY2=JY
                  IF(KZ .LT. KZ1) KZ1=KZ
                  IF(KZ .GT. KZ2) KZ2=KZ
               ENDIF
            ENDDO
         ENDDO
      ENDDO

      !BOX CONTAINING THE VOID
      NXVOID = IX2-IX1+1
      NYVOID = JY2-JY1+1
      NZVOID = KZ2-KZ1+1
      ALLOCATE(MARCA_AUX(NXVOID,NYVOID,NZVOID))
      MARCA_AUX(1:NXVOID,1:NYVOID,1:NZVOID) = MARCA(IX1:IX2,JY1:JY2,KZ1:KZ2)


      !NOW, ITERATE OVER ALL CELLS INSIDE THE BOX CONTAINING THE VOID
      !START BUILDING A PARALELEPIPED FROM THE CELL FURTHEST FROM THE BORDER

      !&&&&&&&&&&&&&&&&&&
      DO WHILE(COUNT(MARCA_AUX .EQ. IND) .GT. 0)
      !&&&&&&&&&&&&&&&&&&

         DIST2 = 0
         DIST = 0
         INIX = 0
         INIY = 0
         INIZ = 0
         IFIX = 0
         IFIY = 0
         IFIZ = 0

         !FIRST, FIND THE FURTHEST CELL FROM THE BORDER
         DO KZ=1,NZVOID
            DO JY=1,NYVOID
               DO IX=1,NXVOID

                  IF(MARCA_AUX(IX,JY,KZ) .NE. IND) CYCLE

                  DIST = MAX(NXVOID,NYVOID,NZVOID)
                  !!!!!!!!!!!!!!!!!!!!!!!!!!!
                  !DISTANCE TO THE BORDER IN EACH
                  !!!!!!!!!!!!!!!!!!!!!!!!!!!
                  !-X
                  II=IX
                  DO WHILE(MARCA_AUX(II,JY,KZ) .EQ. IND .AND. II .GT. 1)
                     II=II-1
                  ENDDO
                  DIST = MIN(DIST,ABS(REAL(IX-II+1)))

                  !+X
                  II=IX
                  DO WHILE(MARCA_AUX(II,JY,KZ) .EQ. IND .AND. II .LT. NXVOID)
                     II=II+1
                  ENDDO
                  DIST = MIN(DIST,ABS(REAL(IX-II+1)))

                  !-Y
                  II=JY
                  DO WHILE(MARCA_AUX(IX,II,KZ) .EQ. IND .AND. II .GT. 1)
                     II=II-1
                  ENDDO
                  DIST = MIN(DIST,ABS(REAL(JY-II+1)))

                  !+Y
                  II=JY
                  DO WHILE(MARCA_AUX(IX,II,KZ) .EQ. IND .AND. II .LT. NYVOID)
                     II=II+1
                  ENDDO
                  DIST = MIN(DIST,ABS(REAL(JY-II+1)))

                  !-Z
                  II=KZ
                  DO WHILE(MARCA_AUX(IX,JY,II) .EQ. IND .AND. II .GT. 1)
                     II=II-1
                  ENDDO
                  DIST = MIN(DIST,ABS(REAL(KZ-II+1)))

                  !+Z
                  II=KZ
                  DO WHILE(MARCA_AUX(IX,JY,II) .EQ. IND .AND. II .LT. NZVOID)
                     II=II+1
                  ENDDO
                  DIST = MIN(DIST,ABS(REAL(KZ-II+1)))
                  !!!!!!!!!!!!!!!!!!!!!!!!!!!

                  IF(DIST .GT. DIST2) THEN
                     DIST2 = DIST
                     INIX = IX
                     INIY = JY
                     INIZ = KZ
                     IFIX = IX
                     IFIY = JY
                     IFIZ = KZ
                  ENDIF

               ENDDO
            ENDDO
         ENDDO

         IF(INIX .EQ. 0 .AND. INIY .EQ. 0 .AND. INIZ .EQ. 0 .AND. &
            IFIX .EQ. 0 .AND. IFIY .EQ. 0 .AND. IFIZ .EQ. 0) EXIT

         !New parallelepiped
         N2(ID) = N2(ID) + 1

         !HAVING THE FURTHEST CELL, BUILD THE PARALELEPIPED
         UVOIDPAR(ID,N2(ID)) = IND
         
         !Void can be expanded at least in one direction
         FLAG = 0
         FLAGX1 = 0
         FLAGX2 = 0
         FLAGY1 = 0
         FLAGY2 = 0
         FLAGZ1 = 0
         FLAGZ2 = 0
         
         !&&&&&&
         DO WHILE(FLAG==0)
         !&&&&&&
            !!!!!    expand along  X
            !X!X!X!X!X!X!X!X!X!X!X!X!X!X!X!X!X!X!X!X!X!X!X!X
            OUTX:DO KZZ=INIZ,IFIZ
            DO JYY=INIY,IFIY
                  !+X
                  II = IFIX
                  IF (II .EQ. NXVOID) THEN
                     FLAGX2 = 1
                  ELSE
                     IF(MARCA_AUX(II+1,JYY,KZZ) .NE. IND) FLAGX2 = 1
                  ENDIF
                  !-X
                  II = INIX
                  IF (II .EQ. 1) THEN
                     FLAGX1 = 1
                  ELSE
                     IF(MARCA_AUX(II-1,JYY,KZZ) .NE. IND) FLAGX1 = 1
                  ENDIF
            ENDDO
            ENDDO OUTX

            !!!!!    expand along  Y
            !X!X!X!X!X!X!X!X!X!X!X!X!X!X!X!X!X!X!X!X!X!X!X!X
            OUTY:DO KZZ=INIZ,IFIZ
            DO IXX=INIX,IFIX
                  !+Y
                  II = IFIY
                  IF(II .EQ. NYVOID) THEN
                     FLAGY2 = 1
                  ELSE
                     IF(MARCA_AUX(IXX,II+1,KZZ) .NE. IND) FLAGY2 = 1
                  ENDIF
                  !-Y
                  II = INIY
                  IF(II .EQ. 1) THEN
                     FLAGY1 = 1
                  ELSE
                     IF(MARCA_AUX(IXX,II-1,KZZ) .NE. IND) FLAGY1 = 1
                  ENDIF
            ENDDO
            ENDDO OUTY

            !!!!!!    expand along  Z
            !X!X!X!X!X!X!X!X!X!X!X!X!X!X!X!X!X!X!X!X!X!X!X!X
            OUTZ:DO JYY=INIY,IFIY
            DO IXX=INIX,IFIX
                  !+Z
                  II = IFIZ
                  IF(II .EQ. NZVOID) THEN
                     FLAGZ2 = 1
                  ELSE
                     IF(MARCA_AUX(IXX,JYY,II+1) .NE. IND) FLAGZ2 = 1
                  ENDIF
                  !-Z
                  II = INIZ
                  IF(II .EQ. 1) THEN
                     FLAGZ1 = 1
                  ELSE
                     IF(MARCA_AUX(IXX,JYY,II-1) .NE. IND) FLAGZ1 = 1
                  ENDIF
            ENDDO
            ENDDO OUTZ

            !%%%  FLAGXYZ=0: face of the parallelepiped can be expanded --> move by one row
            IF(FLAGX1 == 0) INIX = INIX - 1
            IF(FLAGX2 == 0) IFIX = IFIX + 1
            IF(FLAGY1 == 0) INIY = INIY - 1
            IF(FLAGY2 == 0) IFIY = IFIY + 1
            IF(FLAGZ1 == 0) INIZ = INIZ - 1
            IF(FLAGZ2 == 0) IFIZ = IFIZ + 1

            !%%% flag=1 in all the directions: no more expansion is possible
            IF(FLAGX1==1 .AND. FLAGX2==1 .AND. &
               FLAGY1==1 .AND. FLAGY2==1 .AND. &
               FLAGZ1==1 .AND. FLAGZ2==1) FLAG=1 !no faces to be expanded --> exit loop
         !&&&&&&
         ENDDO
         !&&&&&&

         !Mark the cells of the parallelepiped as non-void
         MARCA_AUX(INIX:IFIX,INIY:IFIY,INIZ:IFIZ) = 0

         !Bring the parallelepiped to the original coordinates
         INIX = INIX + IX1 - 1
         INIY = INIY + JY1 - 1
         INIZ = INIZ + KZ1 - 1
         IFIX = IFIX + IX1 - 1
         IFIY = IFIY + JY1 - 1
         IFIZ = IFIZ + KZ1 - 1
         
         !Store the parallelepiped
         INICIOX(ID,N2(ID)) = INIX
         INICIOY(ID,N2(ID)) = INIY
         INICIOZ(ID,N2(ID)) = INIZ
         FINALX(ID,N2(ID)) = IFIX
         FINALY(ID,N2(ID)) = IFIY
         FINALZ(ID,N2(ID)) = IFIZ

         !Find cell center coordinates
         RINIXCO(ID,N2(ID)) = RX1 + DX*(INIX-1) - DX/2
         RINIYCO(ID,N2(ID)) = RY1 + DY*(INIY-1) - DY/2
         RINIZCO(ID,N2(ID)) = RZ1 + DZ*(INIZ-1) - DZ/2
         RFINXCO(ID,N2(ID)) = RX1 + DX*(IFIX-1) + DX/2
         RFINYCO(ID,N2(ID)) = RY1 + DY*(IFIY-1) + DY/2
         RFINZCO(ID,N2(ID)) = RZ1 + DZ*(IFIZ-1) + DZ/2

      !&&&&&&&&&&&&&&&&&&
      ENDDO
      !&&&&&&&&&&&&&&&&&&

      DEALLOCATE(MARCA_AUX)

   !-----------------------
   ENDDO
   !-----------------------

   !########################
   !REDUCTION
   !########################
   DO ID=1,NTH
      DO I=1,N2(ID)
         dNVOID = dNVOID + 1
         dUVOID(dNVOID) = UVOIDPAR(ID,I)
         dINICIOX(dNVOID) = INICIOX(ID,I)
         dINICIOY(dNVOID) = INICIOY(ID,I)
         dINICIOZ(dNVOID) = INICIOZ(ID,I)
         dFINALX(dNVOID) = FINALX(ID,I)
         dFINALY(dNVOID) = FINALY(ID,I)
         dFINALZ(dNVOID) = FINALZ(ID,I)
         dRINIXCO(dNVOID) = RINIXCO(ID,I)
         dRINIYCO(dNVOID) = RINIYCO(ID,I)
         dRINIZCO(dNVOID) = RINIZCO(ID,I)
         dRFINXCO(dNVOID) = RFINXCO(ID,I)
         dRFINYCO(dNVOID) = RFINYCO(ID,I)
         dRFINZCO(dNVOID) = RFINZCO(ID,I)
      ENDDO
   ENDDO
   !########################
   !########################

   DEALLOCATE(N2)
   DEALLOCATE(UVOIDPAR)
   DEALLOCATE(INICIOX,INICIOY,INICIOZ)
   DEALLOCATE(FINALX,FINALY,FINALZ)
   DEALLOCATE(RINIXCO,RINIYCO,RINIZCO)
   DEALLOCATE(RFINXCO,RFINYCO,RFINZCO)

!***************************************************************************
END SUBROUTINE VOID_DECONSTRUCTION
!***************************************************************************


! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

END MODULE VOIDFINDING

! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~