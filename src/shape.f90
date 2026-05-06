!*************************************************************
SUBROUTINE SORT(D,N,NP)
!*************************************************************

       IMPLICIT NONE

       integer N,NP,I,J,K
       real*4  D(NP)
       real*4 P

       DO I=1,N-1
       K=I
       P=D(I)
       DO J=I+1, N
               IF(D(J).GE.P) THEN
                  K=J
                  P=D(J)
               END IF
       END DO
       IF(K.NE.I) THEN
       D(K)=D(I)
       D(I)=P

       END IF

       END DO
       RETURN
!***************************************************************
END SUBROUTINE SORT
!***************************************************************


!*************************************************************
!* This subroutine computes all eigenvalues and eigenvectors *
!* of a real symmetric square matrix A(N,N). On output, ele- *
!* ments of A above the diagonal are destroyed. D(N) returns *
!* the eigenvalues of matrix A. V(N,N) contains, on output,  *
!* the eigenvectors of A by columns. THe normalization to    *
!* unity is made by main program before printing results.    *
!* NROT returns the number of Jacobi matrix rotations which  *
!* were required.                                            *
!* --------------------------------------------------------- *
!* Ref.:"NUMERICAL RECIPES, Cambridge University Press, 1986,*
!*       chap. 11, pages 346-348" [BIBLI 08].                *
!*************************************************************
Subroutine JACOBI(A,N,D,NROT)

       implicit none
       integer N,NROT,ip,iq,ialloc,i,j
       real*4  A(1:N,1:N),D(1:N)
       real*4, pointer :: B(:), Z(:)
       real*4  c,g,h,s,sm,t,tau,theta,tresh

       allocate(B(1:100))   !,stat=ialloc)
       allocate(Z(1:100))   !,stat=ialloc)
       do ip=1, N
         B(ip)=A(ip,ip)
         D(ip)=B(ip)
         Z(ip)=0.d0
       end do
       NROT=0
       do i=1, 50
         sm=0.d0
         do ip=1, N-1           !sum off-diagonal elements
            do iq=ip+1, N
               sm=sm+ABS(A(ip,iq))
            end do
         end do
         if(sm==0.d0) return    !normal return
         if(i.lt.4) then
            tresh=0.2d0*sm**2
         else
            tresh=0.d0
         end if
         do ip=1, N-1
            do iq=ip+1, N
               g=100.d0*ABS(A(ip,iq))
!      after 4 sweeps,skip the rotation if the off-diag element is small

       if((i.gt.4).and.(ABS(D(ip))+g.eq.ABS(D(ip))) &
          .and.(ABS(D(iq))+g.eq.ABS(D(iq)))) then

       A(ip,iq)=0.d0
       else if(ABS(A(ip,iq)).gt.tresh) then
                  h=D(iq)-D(ip)
       if(ABS(h)+g.eq.ABS(h)) then
          t=A(ip,iq)/h
       else
          theta=0.5d0*h/A(ip,iq)
          t=1.d0/(ABS(theta)+SQRT(1.d0+theta**2))
          if(theta.lt.0.d0) t=-t
       end if
       c=1.d0/SQRT(1.d0+t**2)
       s=t*c
          tau=s/(1.d0+c)
       h=t*A(ip,iq)
       Z(ip)=Z(ip)-h
       Z(iq)=Z(iq)+h
       D(ip)=D(ip)-h
       D(iq)=D(iq)+h
       A(ip,iq)=0.d0
       do j=1, ip-1
          g=A(j,ip)
          h=A(j,iq)
          A(j,ip)=g-s*(h+g*tau)
          A(j,iq)=h+s*(g-h*tau)
       end do
       do j=ip+1, iq-1
          g=A(ip,j)
          h=A(j,iq)
          A(ip,j)=g-s*(h+g*tau)
          A(j,iq)=h+s*(g-h*tau)
       end do
       do j=iq+1, N
          g=A(ip,j)
          h=A(iq,j)
          A(ip,j)=g-s*(h+g*tau)
          A(iq,j)=h+s*(g-h*tau)
       end do

       NROT=NROT+1
    end if                   !if ((i.gt.4)...
 end do                    !main iq loop
end do                    !main ip loop
do ip=1, N
   B(ip)=B(ip)+Z(ip)
       D(ip)=B(ip)
       Z(ip)=0.d0
    end do
 end do                    !main i loop
!c       pause ' 50 iterations !'
       return
!***************************************************************************
end subroutine JACOBI
!***************************************************************************



!***************************************************************************
SUBROUTINE VOID_SHAPE(LOW1,LOW2,NVOID,INDICE,NCELLV,UVOID,GXC,GYC,GZC,EPS,IP)
!***************************************************************************
USE COMMONDATA
IMPLICIT NONE
!input variables
INTEGER:: LOW1,LOW2,NVOID
!INTEGER, DIMENSION(NX,NY,NZ):: MARCA
INTEGER, DIMENSION(NVOID):: INDICE, UVOID, NCELLV
REAL*4,DIMENSION(NVOID):: GXC, GYC, GZC
!local variables
INTEGER:: IX,JY,KZ, IND0, I1,I2, NROT, IV, IND
REAL*4:: RR, DELTA, AA, BB, CC, VOLM, VELL
REAL*4, DIMENSION(3):: RK, AXIS, BASEIGENVAL
REAL, DIMENSION(3,3,NVOID):: INERTIA
!output variables
REAL*4, DIMENSION(NVOID) :: EPS, IP
INTEGER, PARAMETER:: NCELLV_MIN=0

EPS(:)=0.
IP(:)=0.
INERTIA(:,:,:)=0.

!inertia tensor for each void
DO KZ=LOW1, LOW2
   DO JY=LOW1, LOW2
      DO IX=LOW1, LOW2
         IND0=MARCA(IX,JY,KZ)
         IF(IND0 .GT. 0) THEN !filter void cell
            RK(1)=RADX(IX)-GXC(IND0)
            RK(2)=RADY(JY)-GYC(IND0)
            RK(3)=RADZ(KZ)-GZC(IND0)

            RR=SQRT(RK(1)**2.+RK(2)**2.+RK(3)**2.)
            DO I2=1,3
               DO I1=1, 3
                  IF(I1 .EQ. I2) DELTA=1 !KHRONAKER DELTA
                  IF(I1 .NE. I2) DELTA=0
                  INERTIA(I1,I2,IND0)=INERTIA(I1,I2,IND0)+DELTA*RR*RR-RK(I1)*RK(I2)
               ENDDO
            ENDDO

         ENDIF
      ENDDO
   ENDDO
ENDDO

DO IV=1,NVOID
   IND=INDICE(IV)
   IF(UVOID(IND) .NE. -1) CYCLE

   !normalize inertia
   DO I2=1,3
      DO I1=1,3
         INERTIA(I1,I2,IND)=INERTIA(I1,I2,IND)/REAL(NCELLV(IND))
      ENDDO
   ENDDO

!diagonalize tensor
   BASEIGENVAL(:)=0.
   CALL JACOBI(INERTIA(:,:,IND),3,BASEIGENVAL,NROT)

!principal axes
   AXIS(1)=SQRT(2.5*(ABS(BASEIGENVAL(2)+BASEIGENVAL(3)-BASEIGENVAL(1))))
   AXIS(2)=SQRT(2.5*(ABS(BASEIGENVAL(3)+BASEIGENVAL(1)-BASEIGENVAL(2))))
   AXIS(3)=SQRT(2.5*(ABS(BASEIGENVAL(1)+BASEIGENVAL(2)-BASEIGENVAL(3))))

   CALL SORT(AXIS,3,3) !from largest to smallest
   AA=AXIS(1)
   BB=AXIS(2)
   CC=AXIS(3)

   EPS(IND)=1.-CC/AA

   VOLM=NCELLV(IND)*DX*DY*DZ ! ACTUAL VOLUME OF THE VOID

   VELL=(4./3)*PI*AA*BB*CC ! VOLUME OF THE ELLIPSOIDAL FIT

   IP(IND)=VOLM/VELL
   IF(IP(IND) .GT. 1.) THEN
      IP(IND)=-99
      EPS(IND)=-99.
   ENDIF

ENDDO
!***************************************************************************
END SUBROUTINE VOID_SHAPE
!***************************************************************************



