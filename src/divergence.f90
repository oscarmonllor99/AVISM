!***************************************************************************
SUBROUTINE DIVER_UNIFORM(NX,NY,NZ,DX,DY,DZ,V2,V3,V4,DIVER,FLAG_PERIODIC)
!***************************************************************************
IMPLICIT NONE
!in 
INTEGER :: FLAG_PERIODIC
INTEGER NX, NY, NZ
REAL*4 DX, DY, DZ
REAL*4 :: V2(NX,NY,NZ), V3(NX,NY,NZ), V4(NX,NY,NZ)
!local
REAL*4, ALLOCATABLE :: DIFFX(:,:,:), DIFFY(:,:,:), DIFFZ(:,:,:)
!out
REAL*4 :: DIVER(NX,NY,NZ)

DIVER(:,:,:) = 0.
ALLOCATE(DIFFX(NX,NY,NZ), DIFFY(NX,NY,NZ), DIFFZ(NX,NY,NZ))
DIFFX = 0.
DIFFY = 0.
DIFFZ = 0.

DIFFX(2:NX-1,:,:) = (V2(3:NX,:,:) - V2(1:NX-2,:,:))/2.
DIFFY(:,2:NY-1,:) = (V3(:,3:NY,:) - V3(:,1:NY-2,:))/2.
DIFFZ(:,:,2:NZ-1) = (V4(:,:,3:NZ) - V4(:,:,1:NZ-2))/2.

!BOUNDARY EXTRAPOLATION
IF (FLAG_PERIODIC .EQ. 0) THEN
   DIFFX(1,:,:)  = 2.*(V2(2,:,:)-V2(1,:,:)) - DIFFX(2,:,:)
   DIFFX(NX,:,:) = 2.*(V2(NX,:,:)-V2(NX-1,:,:)) - DIFFX(NX-1,:,:)

   DIFFY(:,1,:)  = 2.*(V3(:,2,:)-V3(:,1,:)) - DIFFY(:,2,:)
   DIFFY(:,NY,:) = 2.*(V3(:,NY,:)-V3(:,NY-1,:)) - DIFFY(:,NY-1,:)

   DIFFZ(:,:,1)  = 2.*(V4(:,:,2)-V4(:,:,1)) - DIFFZ(:,:,2)
   DIFFZ(:,:,NZ) = 2.*(V4(:,:,NZ)-V4(:,:,NZ-1)) - DIFFZ(:,:,NZ-1)

!PERIODIC BOUNDARY CONDITIONS
ELSE IF (FLAG_PERIODIC .EQ. 1) THEN
   DIFFX(1,:,:)  =  (V2(2,:,:) - V2(NX,:,:))/2.
   DIFFX(NX,:,:) =  (V2(1,:,:) - V2(NX-1,:,:))/2.
   DIFFY(:,1,:)  =  (V3(:,2,:) - V3(:,NY,:))/2.
   DIFFY(:,NY,:) =  (V3(:,1,:) - V3(:,NY-1,:))/2.
   DIFFZ(:,:,1)  =  (V4(:,:,2) - V4(:,:,NZ))/2.
   DIFFZ(:,:,NZ) =  (V4(:,:,1) - V4(:,:,NZ-1))/2.
ENDIF

DIVER = DIFFX/DX + DIFFY/DY + DIFFZ/DZ

DEALLOCATE(DIFFX, DIFFY, DIFFZ)
!***************************************************************************
END SUBROUTINE DIVER_UNIFORM
!**************************************************************************


!***************************************************************************
SUBROUTINE DIVER_FINA_GAS(IRR)
!***************************************************************************
       USE COMMONDATA
       IMPLICIT NONE
!input variables
       INTEGER IRR

!local variables
       INTEGER IR, I, IX, JY, KZ, LOW1, LOW2, N1, N2, N3, NX, NY, NZ
       INTEGER L1, L2, L3, CR1, CR2, CR3, DIM1, DIM2, DIM3, DIM4
       REAL*4 DXR, DYR, DZR
       REAL*4 BAS21, BAS32, BAS43
       REAL*4, ALLOCATABLE:: DIVER_AUX(:,:,:,:)


       NX=INT(REAL(NHYX)*(2.**IRR))
       NY=INT(REAL(NHYX)*(2.**IRR))
       NZ=INT(REAL(NHYX)*(2.**IRR))

      !  WRITE(*,*) NX, NY, NZ

!*-------------------------------------
!*      LEVEL -1 
!*-------------------------------------

       IF(IRR .LT. 0) THEN

          DO KZ=1,NZ
             DO JY=1,NY
                DO IX=1,NX
                   DIVERDMCO(IX,JY,KZ)=0.0
                   DIVERGCO(IX,JY,KZ)=0.0
                ENDDO
             ENDDO
          ENDDO

          !divergence only for cells 2,N-1
          DO KZ=2,NZ-1
             DO JY=2,NY-1
                DO IX=2,NX-1

                   !velocity divergence total
                   BAS21=U2DMCO(IX+1,JY,KZ)-U2DMCO(IX-1,JY,KZ)
                   BAS21=BAS21/(2.0*DX)

                   BAS32=U3DMCO(IX,JY+1,KZ)-U3DMCO(IX,JY-1,KZ)
                   BAS32=BAS32/(2.0*DY)

                   BAS43=U4DMCO(IX,JY,KZ+1)-U4DMCO(IX,JY,KZ-1)
                   BAS43=BAS43/(2.0*DZ)

                   DIVERDMCO(IX,JY,KZ)=BAS21+BAS32+BAS43

                   !velocity divergence for gas
                   BAS21=U2GCO(IX+1,JY,KZ)-U2GCO(IX-1,JY,KZ)
                   BAS21=BAS21/(2.0*DX)

                   BAS32=U3GCO(IX,JY+1,KZ)-U3GCO(IX,JY-1,KZ)
                   BAS32=BAS32/(2.0*DY)

                   BAS43=U4GCO(IX,JY,KZ+1)-U4GCO(IX,JY,KZ-1)
                   BAS43=BAS43/(2.0*DZ)

                   DIVERGCO(IX,JY,KZ)=BAS21+BAS32+BAS43

                ENDDO
             ENDDO


          ENDDO

          DIVERDMCO(1,:,:)=DIVERDMCO(2,:,:)
          DIVERDMCO(NX,:,:)=DIVERDMCO(NX-1,:,:)
          DIVERDMCO(:,1,:)=DIVERDMCO(:,2,:)
          DIVERDMCO(:,NY,:)=DIVERDMCO(:,NY-1,:)
          DIVERDMCO(:,:,1)=DIVERDMCO(:,:,2)
          DIVERDMCO(:,:,NZ)=DIVERDMCO(:,:,NZ-1)
          DIVERGCO(1,:,:)=DIVERGCO(2,:,:)
          DIVERGCO(NX,:,:)=DIVERGCO(NX-1,:,:)
          DIVERGCO(:,1,:)=DIVERGCO(:,2,:)
          DIVERGCO(:,NY,:)=DIVERGCO(:,NY-1,:)
          DIVERGCO(:,:,1)=DIVERGCO(:,:,2)
          DIVERGCO(:,:,NZ)=DIVERGCO(:,:,NZ-1)

       ELSE IF (IRR .GE. 0) THEN 

!*-------------------------------*
!*      LEVEL 0 (DIVER0)
!*-------------------------------*
          IF(ALLOCATED(DIVER0) .EQV. .TRUE.) DEALLOCATE(DIVER0)
          ALLOCATE(DIVER0(NHYX, NHYY, NHYZ)) !COMMON

          DO KZ=2,NHYZ-1
             DO JY=2,NHYY-1
                DO IX=2,NHYX-1

                   BAS21=U2G(IX+1,JY,KZ)-U2G(IX-1,JY,KZ)
                   BAS21=BAS21/(2.0*DX0)

                   BAS32=U3G(IX,JY+1,KZ)-U3G(IX,JY-1,KZ)
                   BAS32=BAS32/(2.0*DY0)

                   BAS43=U4G(IX,JY,KZ+1)-U4G(IX,JY,KZ-1)
                   BAS43=BAS43/(2.0*DZ0)

                   DIVER0(IX,JY,KZ)=BAS21+BAS32+BAS43

                END DO
             END DO
          END DO


          DIVER0(1,:,:)=DIVER0(2,:,:)
          DIVER0(NHYX,:,:)=DIVER0(NHYX-1,:,:)
          DIVER0(:,1,:)=DIVER0(:,2,:)
          DIVER0(:,NHYY,:)=DIVER0(:,NHYY-1,:)
          DIVER0(:,:,1)=DIVER0(:,:,2)
          DIVER0(:,:,NHYZ)=DIVER0(:,:,NHYZ-1)

          IF (IRR .GT. 0) THEN
!*-------------------------------------
!*      LEVEL > 0  (DIVER)
!*-------------------------------------
            
             DIM1=MAXVAL(PATCHNX)
             DIM2=MAXVAL(PATCHNY)
             DIM3=MAXVAL(PATCHNZ)
             DIM4=SUM(NPATCH(0:NL2))
             IF(ALLOCATED(DIVER) .EQV. .TRUE.) DEALLOCATE(DIVER)
             ALLOCATE(DIVER(DIM1, DIM2, DIM3, DIM4)) !COMMON

             DIVER=0.0

             DO IR=1,NL2

                DXR=0.0
                DYR=0.0
                DZR=0.0

                DXR=DX0/(2.0**IR)
                DYR=DY0/(2.0**IR)
                DZR=DZ0/(2.0**IR)

                LOW1=SUM(NPATCH(0:IR-1))+1
                LOW2=SUM(NPATCH(0:IR))
                DO I=LOW1, LOW2

                   N1=PATCHNX(I)
                   N2=PATCHNY(I)
                   N3=PATCHNZ(I)

                   DO KZ=2,N3-1 !exclude the borders of the patches
                      DO JY=2,N2-1
                         DO IX=2,N1-1

                            BAS21=0.0
                            BAS32=0.0
                            BAS43=0.0

                            BAS21=U12G(IX+1,JY,KZ,I)-U12G(IX-1,JY,KZ,I)
                            BAS21=BAS21/(2.0*DXR)

                            BAS32=U13G(IX,JY+1,KZ,I)-U13G(IX,JY-1,KZ,I)
                            BAS32=BAS32/(2.0*DYR)

                            BAS43=U14G(IX,JY,KZ+1,I)-U14G(IX,JY,KZ-1,I)
                            BAS43=BAS43/(2.0*DZR)

                            DIVER(IX,JY,KZ,I)=BAS21+BAS32+BAS43

                         END DO
                      END DO
                   END DO

                END DO
             END DO

!*-------------------------------------
!*      LEVEL == 1, with LEVEL = 0 values
!*-------------------------------------

             ALLOCATE(DIVER_AUX(DIM1, DIM2, DIM3, DIM4)) !LOCAL

             DIVER_AUX=0.0  

             IR=1

             LOW1=SUM(NPATCH(0:IR-1))+1
             LOW2=SUM(NPATCH(0:IR))
             DO I=LOW1, LOW2

                N1=PATCHNX(I)
                N2=PATCHNY(I)
                N3=PATCHNZ(I)
                L1=PATCHX(I)
                L2=PATCHY(I)
                L3=PATCHZ(I)

                DO KZ=1,N3
                   DO JY=1,N2
                      DO IX=1,N1

                         BAS21=0.0
                         BAS32=0.0
                         BAS43=0.0

                         IF (IX.LT.3.OR.IX.GT.N1-2.OR. &
                              JY.LT.3.OR.JY.GT.N2-2.OR. &
                              KZ.LT.3.OR.KZ.GT.N3-2) THEN

                            CR1=INT((IX+1)/2)+L1-1
                            CR2=INT((JY+1)/2)+L2-1
                            CR3=INT((KZ+1)/2)+L3-1

                            BAS21=U2G(CR1+1,CR2,CR3)-U2G(CR1-1,CR2,CR3)
                            BAS21=BAS21/(2.0*DX0)

                            BAS32=U3G(CR1,CR2+1,CR3)-U3G(CR1,CR2-1,CR3)
                            BAS32=BAS32/(2.0*DY0)

                            BAS43=U4G(CR1,CR2,CR3+1)-U4G(CR1,CR2,CR3-1)
                            BAS43=BAS43/(2.0*DZ0)

                            DIVER_AUX(IX,JY,KZ,I)=BAS21+BAS32+BAS43

                         END IF

                      END DO
                   END DO
                END DO
             END DO


!*--------------------------------------------------
!*      LEVEL == 1, put LEVEL = 0 values in borders
!*--------------------------------------------------
             DO IR=1, 1 !NL

                DXR=DX0/(2.0**IR)
                DYR=DY0/(2.0**IR)
                DZR=DZ0/(2.0**IR)

                LOW1=SUM(NPATCH(0:IR-1))+1
                LOW2=SUM(NPATCH(0:IR))
                DO I=LOW1, LOW2

                   N1=PATCHNX(I)
                   N2=PATCHNY(I)
                   N3=PATCHNZ(I)

                   DO KZ=1, N3
                      DO JY=1, N2
                         DO IX=1, N1


                            IF(IX.EQ.2.OR.IX.EQ.N1-1.OR. &
                               JY.EQ.2.OR.JY.EQ.N2-1.OR. &
                               KZ.EQ.2.OR.KZ.EQ.N3-1) THEN
                               DIVER(IX,JY,KZ,I)=0.5*(DIVER_AUX(IX,JY,KZ,I)+ &
                                                      DIVER(IX,JY,KZ,I))
                            END IF


                            IF(IX.EQ.1.OR.IX.EQ.N1.OR. &
                               JY.EQ.1.OR.JY.EQ.N2.OR. &
                               KZ.EQ.1.OR.KZ.EQ.N3) THEN
                               DIVER(IX,JY,KZ,I)=DIVER_AUX(IX,JY,KZ,I)
                            END IF

                         END DO
                      END DO
                   END DO
                END DO
             END DO !LOOP ON LEVELS

             DEALLOCATE(DIVER_AUX)

          ENDIF !IRR > 0

          ENDIF !IR test

          RETURN
          !

!***************************************************************************
END SUBROUTINE DIVER_FINA_GAS
!***************************************************************************