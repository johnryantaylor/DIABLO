       SUBROUTINE USER_RHS_CHAN_PHYSICAL
       include 'header'
! Here, you can add terms to the right hand side
! of the momentum and scalar equations.  
! The right hand side forcing arrays, CF2, CF2, CF3, CFTH
! are in Fourier space.  The velocity and scalars are available 
! in physical space.
! S1 is available as a working variable

       integer i,j,k,n

       real*8 alpha

!       CALL SLIP_VEL

      RETURN 
      END


       SUBROUTINE USER_RHS_CHAN_FOURIER
       include 'header'
! Here, you can add terms to the right hand side
! of the momentum and scalar equations.  
! The right hand side forcing arrays, CF2, CF2, CF3, CFTH
! are in Fourier space.  The velocity and scalars are available 
! in Fourier space.
! S1 is available as a working variable

       integer i,j,k,n

! Advection owing to thermal wind
      IF ((FLAVOR.eq.'Front').and.(I_RO.ne.0.d0)) THEN
      DO N=1,N_TH
! Loop over all scalars

! Add thermal wind advection to the momentum equations
      do j=JSTART,JEND
      do k=0,TNKZ
      do i=0,NXP-1
        CF1(I,K,J)=CF1(I,K,J)
     &           -(DRHODX(N)*RI(N)*GYF(J)/I_RO)
     &                            *CIKZ(K)*CU1(I,K,J)
     &           -(-1.d0*DRHODZ(N)*RI(N)*GYF(J)/I_RO)
     &                            *CIKX(I)*CU1(I,K,J)
     &                   -(-1.d0*DRHODZ(N)*RI(N)/I_RO)
     &                      *0.5d0*(CU2(I,K,J)+CU2(I,K,J+1))
        CF3(I,K,J)=CF3(I,K,J)
     &           -(DRHODX(N)*RI(N)*GYF(J)/I_RO)
     &                           *CIKZ(K)*CU3(I,K,J)
     &           -(-1.d0*DRHODZ(N)*RI(N)*GYF(J)/I_RO)
     &                           *CIKX(I)*CU3(I,K,J)
     &                   -(DRHODX(N)*RI(N)/I_RO)
     &                      *0.5d0*(CU2(I,K,J)+CU2(I,K,J+1))
      end do
      end do
      end do

      do j=2,NY
      do k=0,TNKZ
      do i=0,NXP-1
        CF2(I,K,J)=CF2(I,K,J)
     &            -(DRHODX(N)*RI(N)*GY(J)/I_RO)
     &                     *CIKZ(K)*CU2(I,K,J)
     &            -(-1.d0*DRHODZ(N)*RI(N)*GY(J)/I_RO)
     &                     *CIKX(I)*CU2(I,K,J)
      end do
      end do
      end do

! Add advection by thermal wind to the scalar equations
      DO J=JSTART_TH(N),JEND_TH(N)
        DO K=0,TNKZ
          DO I=0,NXP-1
            CFTH(I,K,J,N)=CFTH(I,K,J,N)
     &     -(RI(N)/I_RO)*DRHODX(N)*GYF(J)
     &                 *CIKZ(K)*CTH(I,K,J,N)
     &     -(RI(N)/I_RO)*-1.d0*DRHODZ(N)*GYF(J)
     &                 *CIKX(I)*CTH(I,K,J,N)
          END DO
        END DO
      END DO

! End do N_TH
      END DO

      do j=2,NY
      do k=0,TNKZ
      do i=0,NXP-1
        CF2(I,K,J)=CF2(I,K,J)
     &            -(DRHODX(N)*RI(N)*GY(J)/I_RO)
     &                     *CIKZ(K)*CU2(I,K,J)
     &            -(-1.d0*DRHODZ(N)*RI(N)*GY(J)/I_RO)
     &                     *CIKX(I)*CU2(I,K,J)
      end do
      end do
      end do

! Add advection by thermal wind to the scalar equations
      DO J=JSTART_TH(N),JEND_TH(N)
        DO K=0,TNKZ
          DO I=0,NXP-1
            CFTH(I,K,J,N)=CFTH(I,K,J,N)
     &     -(RI(N)/I_RO)*DRHODX(N)*GYF(J)
     &                 *CIKZ(K)*CTH(I,K,J,N)
     &     -(RI(N)/I_RO)*-1.d0*DRHODZ(N)*GYF(J)
     &                 *CIKX(I)*CTH(I,K,J,N)
          END DO
        END DO
      END DO

! End do N_TH
      END DO

      END IF

! Add sponge layer forcing
      DO N=1,N_TH
        CALL SPONGE_TH(N)
      END DO
      CALL SPONGE_VEL

      RETURN 
      END


      SUBROUTINE USER_RHS_PER_PHYSICAL
      include 'header'
C Optionally, add forcing terms to the right hand side
C of the momentum and scalar evolution equations
C Here, CF1, CF2, CF3, and CFTH are the forcing arrays in Fourier space
C U1, U2, U3, and TH are the velocity and scalar variables in Physical space

      integer i,j,k,n

C For forced homogeneous turbulence:
C Add some forcing to the system to keep the Batchelor scale fixed
      EK=0.d0
      DO J=0,NYM
        DO K=0,NZM
          DO I=0,NXM
              EK=EK+U1(I,K,J)**2.d0+U2(I,K,J)**2.d0+U3(I,K,J)**2.d0
          END DO
        END DO
      END DO
C Note, that each cell has the same volume, so we can just average over all points
      EK=EK/dble(NX*NY*NZ)

! Scale EK by an amount to compensate for dissipation from 2/3 de-aliasing:
      EK=0.8d0*EK
      DO J=0,NYM
        DO K=0,NZM
          DO I=0,NXM
            S1(I,K,J)=(EPSILON_TARGET/EK)*U1(I,K,J)
          END DO
        END DO
      END DO
      CALL FFT_XZY_TO_FOURIER(S1,CS1)
      DO J=0,TNKY
        DO K=0,TNKZ
          DO I=0,NKX
            CF1(I,K,J)=CF1(I,K,J)+CS1(I,K,J)
          END DO
        END DO
      END DO
      DO J=0,NYM
        DO K=0,NZM
          DO I=0,NXM
            S1(I,K,J)=(EPSILON_TARGET/EK)*U2(I,K,J)
          END DO
        END DO
      END DO
      CALL FFT_XZY_TO_FOURIER(S1,CS1)
      DO J=0,TNKY
        DO K=0,TNKZ
          DO I=0,NKX
            CF2(I,K,J)=CF2(I,K,J)+CS1(I,K,J)
          END DO
        END DO
      END DO
      DO J=0,NYM
        DO K=0,NZM
          DO I=0,NXM
            S1(I,K,J)=(EPSILON_TARGET/EK)*U3(I,K,J)
          END DO
        END DO
      END DO
      CALL FFT_XZY_TO_FOURIER(S1,CS1)
      DO J=0,TNKY
        DO K=0,TNKZ
          DO I=0,NKX
            CF3(I,K,J)=CF3(I,K,J)+CS1(I,K,J)
          END DO
        END DO
      END DO

      RETURN
      END

      SUBROUTINE USER_RHS_DUCT_PHYSICAL
      include 'header'
      RETURN
      END

      SUBROUTINE USER_RHS_CAVITY_PHYSICAL
      include 'header'
      RETURN
      END 


      SUBROUTINE SPONGE_TH(N)
! This subroutine applies a sponge relaxation (Rayleigh damping) towards a
! specified background state for the temperature field
! The intention is to allow an open boundary
      include 'header'
      integer i,j,k,n
      real*8 L_sponge,L_bottom
      real*8 SPONGE_AMP

! The following variables will store the background state
      real*8 TH_0(-1:NY+1)

      real*8 RI_B(0:NY+1)

! This variable will hold the forcing rate
      real*8 SPONGE_SIGMA(0:NY+1)

! Set the amplitude of the sponge
      SPONGE_AMP=0.005d0
! Set the top of the sponge layer in physical units
      L_sponge=-40.d0
        DO J=0,NY+1
! Quadratic damping at lower wall
         if (GYF(J).lt.L_sponge) then
           SPONGE_SIGMA(j)=SPONGE_AMP*((L_sponge-GYF(J))
     &       /(L_sponge-L_bottom))**2.d0
         else
           SPONGE_SIGMA(j)=0.d0
         end if
        END DO

! Set the profile for relaxing the mean TH
      DO J=0,NY+1
        TH_0(J)=TH_BC_YMIN_C1(N)*GYF(J)
      END DO       

! For MLI latmix
      if (n.eq.1) then
      TH_0(0)=0.d0
      do j=1,NY+1
        RI_B(J)=20.d0
        TH_0(J)=TH_0(J-1)
     &          +DY(J)*RI_B(J)*(RI(N)*DRHODX(N))**2.d0
     &           /I_RO**2.d0/RI(N)
      end do
      else
        do j=0,NY+1
          TH_0(j)=0.d0
        end do
      end if

! Add damping to R-K terms
! Damp the perturbations towards 0
      do k=0,TNKZ
         do i=0,NXP-1
           if ((RANKZ.ne.0).or.(i.ne.0).or.(k.ne.0)) then
           do j=JSTART_TH(N),JEND_TH(N)
              CFTH(i,k,j,n)=CFTH(i,k,j,n)
     &                 -SPONGE_SIGMA(j)*(CTH(i,k,j,n)-0.)
           end do
           end if
         end do
      end do
! Damp the mean gradient towards TH_0
        do j=JSTART_TH(N),JEND_TH(N)
          CFTH(0,0,j,n)=CFTH(0,0,j,n)-SPONGE_SIGMA(j)
     &          *(CTH(0,0,j,n)-TH_0(J))
        end do
      end if

      return
      end
  
      SUBROUTINE SPONGE_VEL
! This subroutine applies a sponge relaxation (Rayleigh damping) towards a
! specified background state
! The intention is to allow an open boundary
      include 'header'
      integer i,j,k

      real*8 L_sponge,L_bottom
      real*8 SPONGE_AMP

! The following variables will store the background state
      real*8 U1_0(-1:NY+1), U2_0(0:NY+1), U3_0(-1:NY+1)

! This variable will hold the forcing rate
      real*8 SPONGE_SIGMA(0:NY+1)

! Set the amplitude of the sponge
      SPONGE_AMP=0.005d0
! Set the top of the sponge layer in physical units
      L_sponge=-40.d0
! Set the bottom of the computational domain in physical units
      L_bottom=-50.d0
        DO J=0,NY+1
! Quadratic damping at lower wall
         if (GYF(J).lt.L_sponge) then
           SPONGE_SIGMA(j)=SPONGE_AMP*((L_sponge-GYF(J))
     &       /(L_sponge-L_bottom))**2.d0
         else
           SPONGE_SIGMA(j)=0.d0
         end if
        END DO

! Set the background state
! Here, set the background to be geostrophic, with a linear temperature profile
      do j=0,NY+1
        U1_0(j)=0.d0
        U3_0(j)=0.d0
      end do
      do j=0,NY+1
        U2_0(j)=0.d0
      end do

! Add damping function to explicit R-K 
       do k=0,TNKZ
         do i=0,NXP-1
           if ((RANKZ.ne.0).or.(i.ne.0).or.(k.ne.0)) then
           do j=jstart,jend
             CF1(I,K,J)=CF1(I,K,J)-SPONGE_SIGMA(j)*(CU1(i,k,j)-0.d0)
             CF3(I,K,J)=CF3(I,K,J)-SPONGE_SIGMA(j)*(CU3(i,k,j)-0.d0)
           end do
           do j=2,NY
             CF2(I,K,J)=CF2(I,K,J)-
     &        0.5*(SPONGE_SIGMA(j)+SPONGE_SIGMA(j+1))*(CU2(i,k,j)-0.d0)
           end do
           end if   
         end do
      end do
! Damp mean flow
      if (RANKZ.eq.0) then
      do j=jstart,jend
        CF1(0,0,j)=CF1(0,0,j)-SPONGE_SIGMA(j)*(CU1(0,0,j)-U1_0(j))
        CF3(0,0,j)=CF3(0,0,j)-SPONGE_SIGMA(j)*(CU3(0,0,j)-U3_0(j))
      end do
      do j=2,NY
        CF2(0,0,j)=CF2(0,0,j)-SPONGE_SIGMA(j)*(CU2(0,0,j)-U2_0(j))
      end do
      end if

      return
      end

      SUBROUTINE SLIP_VEL
! This subroutine adds advection by a slip velocity to some scalars
      include 'header'
      integer i,j,k,n,J1_TH(1:N_TH),J2_TH(1:N_TH)
      real*8 W_S(0:NY+1,1:N_TH)
   
! Set indices corresponding to start and end of GYF grid
      do n=1,N_TH 
      IF (RANKY.eq.NPROCY-1) then
! We are at the upper wall
            J1_TH(n)=JSTART_TH(n)
            J2_TH(n)=NY-1
      ELSE IF (RANKY.eq.0) then
! We are at the lower wall
            J1_TH(n)=2
            J2_TH(n)=JEND_TH(n)
      ELSE
! We are on a middle process
            J1_TH(n)=JSTART_TH(n)
            J2_TH(n)=JEND_TH(n)
      END IF
      end do

! First, set the slip velocity
      do j=0,NY+1
        W_S(j,1)=0.d0
        W_S(j,2)=0.d0
        W_S(j,3)=0.00005d0
        W_S(j,4)=0.0005d0
        W_S(j,5)=0.005d0
      end do

      IF (RANKY.eq.NPROCY-1) THEN
! We are on a process at the top boundary
! Set the slip velocity to zero at GY(NY) (and ghost cells)
      do n=1,N_TH
        W_S(NY,n)=0.d0
        W_S(NY+1,n)=0.d0
      end do
      ELSE IF (RANKY.eq.0) THEN
! We are on a process at the bottom boundary
! Set the slip velocit to zero at GY(2) (and ghost cells)
      do n=1,N_TH
        W_S(0,n)=0.d0
        W_S(1,n)=0.d0
        W_S(2,n)=0.d0
      end do
      END IF

      do n=1,N_TH
        DO J=J1_TH(N),J2_TH(N)
          DO K=0,NZP-1
            DO I=0,NXM
! Central differencing
!              S1(I,K,J)=
!     &     ((TH(I,K,J+1,N)*W_S(J+1,N) + TH(I,K,J,N)*W_S(J+1,N)
!     &    -TH(I,K,J,N)*W_S(J,N)-TH(I,K,J-1,N)*W_S(J,n))/(2.d0*DYF(J)))
! Second order Upwinding
!              S1(I,K,J)=(W_S(J+1,N)*TH(I,K,J,N)
!     &               -W_S(J,N)*(TH(I,K,J,N)+TH(I,K,J-1,N))/2.d0)
!     &                 /(GYF(j)-GY(j))
! First order upwinding 
              S1(I,K,J)=(W_S(J+1,N)*TH(I,K,J,N)
     &               -W_S(J,N)*TH(I,K,J-1,N))
     &                 /(GYF(j)-GYF(j-1))

!              S1(I,K,J)=0.5d0*(W_S(J+1,N)+W_S(J,N))
!     &              *(TH(I,K,J,N)-TH(I,K,J-1,N))/(GYF(J)-GYF(J-1))
            END DO
          END DO
        END DO
        CALL FFT_XZ_TO_FOURIER(S1,CS1,0,NY+1)
        DO J=J1_TH(N),J2_TH(N)
          DO K=0,TNKZ
            DO I=0,NXP-1
              CFTH(I,K,J,N)=CFTH(I,K,J,N) - CS1(I,K,J)
            END DO
          END DO
        END DO
      end do

      RETURN 
      END
 
=======
        do j=0,NY+1
          CFTH(0,0,j,n)=CFTH(0,0,j,n)-SPONGE_SIGMA(j)
     &          *(CTH(0,0,j,n)-TH_0(J))
        end do

      return
      end
  
      SUBROUTINE SPONGE_VEL
! This subroutine applies a sponge relaxation (Rayleigh damping) towards a
! specified background state
! The intention is to allow an open boundary
      include 'header'
      integer i,j,k

      real*8 L_sponge,L_bottom
      real*8 SPONGE_AMP

! The following variables will store the background state
      real*8 U1_0(-1:NY+1), U2_0(0:NY+1), U3_0(-1:NY+1)

! This variable will hold the forcing rate
      real*8 SPONGE_SIGMA(0:NY+1)

! Set the amplitude of the sponge
      SPONGE_AMP=0.0001d0
! Set the top of the sponge layer in physical units
      L_sponge=-120.d0
! Set the bottom of the computational domain in physical units
      L_bottom=-140.d0
        DO J=0,NY+1
! Quadratic damping at lower wall
         if (GYF(J).lt.L_sponge) then
           SPONGE_SIGMA(j)=SPONGE_AMP*((L_sponge-GYF(J))
     &       /(L_sponge-L_bottom))**2.d0
         else
           SPONGE_SIGMA(j)=0.d0
         end if
        END DO

! Set the background state
! Here, set the background to be geostrophic, with a linear temperature profile
      do j=0,NY+1
        U1_0(j)=0.d0
        U3_0(j)=0.d0
      end do
      do j=0,NY+1
        U2_0(j)=0.d0
      end do

! Add damping function to explicit R-K 
       do k=0,TNKZ
         do i=0,NKX
           if ((i.ne.0).or.(k.ne.0)) then
           do j=jstart,jend
             CF1(I,K,J)=CF1(I,K,J)-SPONGE_SIGMA(j)*(CU1(i,k,j)-0.d0)
             CF3(I,K,J)=CF3(I,K,J)-SPONGE_SIGMA(j)*(CU3(i,k,j)-0.d0)
           end do
           do j=1,NY
             CF2(I,K,J)=CF2(I,K,J)-
     &        0.5*(SPONGE_SIGMA(j)+SPONGE_SIGMA(j+1))*(CU2(i,k,j)-0.d0)
           end do
           end if   
         end do
      end do
! Damp mean flow
      do j=jstart,jend
        CF1(0,0,j)=CF1(0,0,j)-SPONGE_SIGMA(j)*(CU1(0,0,j)-U1_0(j))
        CF3(0,0,j)=CF3(0,0,j)-SPONGE_SIGMA(j)*(CU3(0,0,j)-U3_0(j))
      end do
      do j=1,NY
        CF2(0,0,j)=CF2(0,0,j)-SPONGE_SIGMA(j)*(CU2(0,0,j)-U2_0(j))
      end do

      return
      end



