      subroutine les_chan
C This subroutine models the terms owing to the subgrid scale stress
C if the computation is to be treated as an LES not a DNS
C This subroutine should be called when the velocity is in fourier space 
C in the periodic directions, on output, the velocity will be 
C in physical space.
C It is assumed that the test filter and the LES filter are performed
C by the same operation
C On output S1 should contain |S| which may be used again in les_chan_th
C if for the subgrid scalar dissipation

      include 'header'
      include 'header_les'

      REAL*8 TH_real(0:NX+1,0:NZP+1,0:NY+1)
      COMPLEX*16 CTMP(0:NXP,0:NZ+1,0:NY+1)

      integer N
      REAL*8 DELTA_TMP_Y, DELTA_TMP_YF
      REAL*8 C_lb
      REAL*8 EPS_DELTA
      parameter (EPS_DELTA=0.000000001d0)
      parameter (C_lb=0.0625d0) !1/16
      integer i,j,k,l,m,ij

      real*8 S1_mean(0:NY+1)
      real*8 NU_T_mean(0:NY+1)
      real*8 EPS_SGS1_MEAN(0:NY+1)
      real*8 U3_bar(0:NY+1)
      real*8 U1_bar(0:NY+1)
 
      real*8 C_SMAG
      parameter (C_SMAG=0.17d0)
      real*8 C_AMD
      parameter (C_AMD=0.2887d0) 
      real*8 DELTA_Y(0:NY+1),DELTA_YF(0:NY+1) 
      real*8 deltax,deltay,deltaz
      real*8 alpha_sgs
      real*8 denominator_sum

! Array for writing HDF5 files
      real*8 Diag(1:NY)
      character*20 gname

      character*35 FNAME

! Array to store the velocity index for each component of the strain rate tensor
      integer U_index1(6)
      integer U_index2(6)

! Here, alpha is the test/LES filter width ratio
      parameter (alpha_sgs=2.449d0)
! beta is the LES/grid filter width ratio
      beta_sgs=3.d0

! Set the indices that are used when adding the off-diagnoal SGS stress terms
      IF (RANKY.eq.NPROCY-1) then
! We are at the upper wall
            J1=JSTART
            J2=NY-1
      ELSE IF (RANKY.eq.0) then
! We are at the lower wall
            J1=2
            J2=JEND
      ELSE
! We are on a middle process
            J1=JSTART
            J2=JEND
      END IF

! First, for all models, apply boundary conditions to the velocity field
! (fill ghost cells) to ensure accurate calculation of gradients
C Apply Boundary conditions to velocity field
      IF (USE_MPI) THEN
        CALL APPLY_BC_VEL_MPI
        CALL GHOST_CHAN_MPI
      ELSE
        CALL APPLY_BC_VEL_LOWER
        CALL APPLY_BC_VEL_UPPER
      END IF

! If we are using Neuman boundary conditions, over-write the values of the
! velocity at the ghost cells so that the LES model doesn't use the large
! velocity gradient
      CALL APPLY_BC_LES

      if (LES_MODEL_TYPE.EQ.1) then
!     Constant Smagorinsky model

! First, compute the rate of strain tensor S_ij

      call compute_strain_chan


! Convert the velocity to physical space
      call FFT_XZ_TO_PHYSICAL(CU1,U1,0,NY+1)
      call FFT_XZ_TO_PHYSICAL(CU2,U2,0,NY+1)
      call FFT_XZ_TO_PHYSICAL(CU3,U3,0,NY+1)

! Compute |S| at GYF points, store in S1
! Interpolation to GYF points is easy since by definition
! GYF points are exactly midway between neighboring GY points
       DO J=JSTART,JEND
        DO K=0,NZP-1
          DO I=0,NXM
            S1(I,K,J)=SQRT(
     &                2.d0*Sij1(I,K,J)**2.d0
     &               +4.d0*(0.5d0*(Sij4(I,K,J+1)+Sij4(I,K,J)))**2.d0
     &               +4.d0*Sij5(I,K,J)**2.d0
     &               +2.d0*Sij2(I,K,J)**2.d0
     &               +4.d0*(0.5d0*(Sij6(I,K,J+1)+Sij6(I,K,J)))**2.d0 
     &               +2.d0*Sij3(I,K,J)**2.d0 )
          END DO
        END DO
       END DO

! Compute |S| at GY points and store in TEMP
       DO J=2,NY+1
        DO K=0,NZP-1
          DO I=0,NXM
            TEMP(I,K,J)=SQRT(
     &                2.d0*((Sij1(I,K,J)*DYF(j-1)+Sij1(I,K,J-1)*DYF(j))
     &                               /(2.d0*DY(j)))**2.d0
     &               +4.d0*Sij4(I,K,J)**2.d0
     &               +4.d0*((Sij5(I,K,J)*DYF(j-1)+Sij5(I,K,J-1)*DYF(j))
     &                               /(2.d0*DY(j)))**2.d0
     &               +2.d0*((Sij2(I,K,J)*DYF(j-1)+Sij2(I,K,J-1)*DYF(j))
     &                               /(2.d0*DY(j)))**2.d0
     &               +4.d0*Sij6(I,K,J)**2.d0
     &               +2.d0*((Sij3(I,K,J)*DYF(j-1)+Sij3(I,K,J-1)*DYF(j))
     &                               /(2.d0*DY(j)))**2.d0
     &                )
          END DO
        END DO
       END DO


! Now, compute |S|*S_ij, storing in Sij
! First compute at GYF points 
       DO J=JSTART,JEND
        DO K=0,NZP-1
          DO I=0,NXM
            Sij1(I,K,J)=S1(I,K,J)*Sij1(I,K,J)
            Sij5(I,K,J)=S1(I,K,J)*Sij5(I,K,J)
! Sij2 is added through an implicit eddy viscosity
            Sij2(I,K,J)=0.d0
            Sij3(I,K,J)=S1(I,K,J)*Sij3(I,K,J)
          END DO
        END DO
       END DO


! Now, compute at |S|*S_ij at GY points
       DO J=2,NY+1
        DO K=0,NZP-1
          DO I=0,NXM
! The terms dU1/dy and dU3/dy in CSij4(:,:,:) and CSij6(:,:,:) respectively
! are subtracted out from Sij here since they are treated implicitly
! in eddy viscosity terms
            Sij4(I,K,J)=TEMP(I,K,J)
     &        *(Sij4(I,K,J)-0.5*(U1(I,K,J)-U1(I,K,J-1))/DY(j))
            Sij6(I,K,J)=TEMP(I,K,J)
     &        *(Sij6(I,K,J)-0.5*(U3(I,K,J)-U3(I,K,J-1))/DY(j))
          END DO
        END DO
       END DO

! We now have |S|*S_ij stored in Sij in Physical space

! Compute the filter lengthscale
! Absorb -2.d0*C_SMAG**2.d0 here for effienciency
       DO J=1,NY+1
! At GYF points:
! Constant Smagorinsky
        DELTA_YF(J)=-2.d0*C_SMAG**2.d0
     &     *(DX(1)*beta_sgs*DYF(J)*2.d0*DZ(1)*beta_sgs)**(2.d0/3.d0)
! Wall Damping
!        DELTA_YF(J)=
!     &    -2.d0*(0.1d0*(1.d0-exp((-GYF(J)/(NU*25.d0))**3.d0)))**2.d0
!     &            *(DX(1)*beta_sgs*DYF(J)*2.d0*DZ(1)*beta_sgs)**(2.d0/3.d0)

       END DO

       DO J=1,NY+1
! At GY points:
! Constant Smagorinsky
        DELTA_Y(J)=-2.d0*C_SMAG**2.d0
     &        *(DX(1)*beta_sgs*DY(J)*2.d0*DZ(1)*beta_sgs)**(2.d0/3.d0)
! Wall Damping
!        DELTA_Y(J)=
!     &    -2.d0*(0.1d0*(1.d0-exp((-GY(J)/(NU*25.d0))**3.d0)))**2.d0
!     &            *(DX(1)*beta_sgs*DY(J)*2.d0*DZ(1)*beta_sgs)**(2.d0/3.d0)
       END DO

! Get the eddy viscosity at GY points
! NU_T = (C_S^2 * DELTA^2)*|S|
       DO J=2,NY
        DO K=0,NZP-1
          DO I=0,NXM
           NU_T(I,K,J)=-0.5d0*DELTA_Y(J)*TEMP(I,K,J)
          END DO
        END DO
       END DO

! Now that we have calculated NU_T, set the value at the ghost cells
! by sharing with neighboring processes.  This subroutine also sets
! the value of NU_T to zero at both walls
      CALL GHOST_LES_MPI

! Convert the stress tensor to Fourier space


      CALL FFT_XZ_TO_FOURIER(Sij1,CSij1,0,NY+1)
! Sij2 is added through an implicit eddy viscosity
!      CALL FFT_XZ_TO_FOURIER(Sij2,CSij2,0,NY+1)
      CALL FFT_XZ_TO_FOURIER(Sij3,CSij3,0,NY+1)
      CALL FFT_XZ_TO_FOURIER(Sij4,CSij4,0,NY+1)
      CALL FFT_XZ_TO_FOURIER(Sij5,CSij5,0,NY+1)
      CALL FFT_XZ_TO_FOURIER(Sij6,CSij6,0,NY+1)

! Now, compute TAU, store in the corresponging Sij
       DO J=1,NY+1
        DO K=0,TNKZ
          DO I=0,NXP-1
            CSij1(I,K,J)=DELTA_YF(J)*CSij1(I,K,J)
            CSij5(I,K,J)=DELTA_YF(J)*CSij5(I,K,J)
! CSij2(:,:,:) is added through an implicit eddy viscosity
!            CSij2(I,K,J)=DELTA_YF(J)*CSij2(I,K,J)
            CSij3(I,K,J)=DELTA_YF(J)*CSij3(I,K,J)
          END DO
        END DO
       END DO
       DO J=2,NY+1
        DO K=0,TNKZ
          DO I=0,NXP-1
            CSij4(I,K,J)=DELTA_Y(J)*CSij4(I,K,J)
            CSij6(I,K,J)=DELTA_Y(J)*CSij6(I,K,J)
          END DO
        END DO
       END DO


! tau_ij is now contained in CSij in Fourier space

      else if ((LES_MODEL_TYPE.EQ.2).or.(LES_MODEL_TYPE.eq.3)) then
! Here, use a dynamic smagorinsky model with or without scale similar part

        stop 'ERROR: LES_MODEL_TYPE=2, 3 not supported yet in MPI'  

      else if (LES_MODEL_TYPE.EQ.4) then
!     Anisotrophic minimum-dissipation model Rozema

! First, compute the rate of strain tensor S_ij

      call compute_strain_chan

! Second, compute the rate of rotation tensor omega_ij 

      call compute_rotation_chan

! Convert the velocity to physical space
      call FFT_XZ_TO_PHYSICAL(CU1,U1,0,NY+1)
      call FFT_XZ_TO_PHYSICAL(CU2,U2,0,NY+1)
      call FFT_XZ_TO_PHYSICAL(CU3,U3,0,NY+1)

!Set filter length (based on grid size) in x,z directions
!Based on constant Smag code
         deltax = (DX(1)*beta_sgs)**2.d0
         deltaz = (DZ(1)*beta_sgs)**2.d0

! Compute max{0,-delta_k*(I3-I4)}/(I1-I2) at GYF points, store in S1
! Interpolation to GYF points is easy since by definition
! GYF points are exactly midway between neighboring GY points
       DO J=JSTART,JEND
!Set filter length (based on grid size) in y direction
!Based on constant Smag code        
         deltay=(DYF(J)*2.d0)**2.d0
        DO K=0,NZP-1
          DO I=0,NXM

!First calculate delta_k*I3 and store it in S1.  
            S1(I,K,J)=
     &               deltax*Sij1(I,K,J)**3.d0
     &              +deltay*Sij2(I,K,J)**3.d0
     &              +deltaz*Sij3(I,K,J)**3.d0
     &       +(2.d0*deltax+deltay)*Sij1(I,K,J)
     &              *(0.5d0*(Sij4(I,K,J+1)+Sij4(I,K,J)))**2.d0
     &       +(2.d0*deltax+deltaz)*Sij1(I,K,J)*Sij5(I,K,J)**2.d0
     &       +(2.d0*deltay+deltax)*Sij2(I,K,J)
     &              *(0.5d0*(Sij4(I,K,J+1)+Sij4(I,K,J)))**2.d0
     &       +(2.d0*deltay+deltaz)*Sij2(I,K,J)
     &              *(0.5d0*(Sij6(I,K,J+1)+Sij6(I,K,J)))**2.d0
     &       +(2.d0*deltaz+deltax)*Sij3(I,K,J)*Sij5(I,K,J)**2.d0
     &       +(2.d0*deltaz+deltay)*Sij3(I,K,J)
     &              *(0.5d0*(Sij6(I,K,J+1)+Sij6(I,K,J)))**2.d0
     &       +2.d0*(deltax+deltay+deltaz)
     &              *(0.5d0*(Sij4(I,K,J+1)+Sij4(I,K,J)))
     &              *Sij5(I,K,J)*(0.5d0*(Sij6(I,K,J+1)+Sij6(I,K,J)))

!Then calculate -delta_k*(I3-I4) = -S1+delta_k*I4 and store in S1.

            S1(I,K,J)=(-S1(I,K,J)
     &        -deltay*Sij1(I,K,J)
     &                    *(0.5d0*(omgij4(I,K,J+1)+omgij4(I,K,J)))**2.d0
     &        -deltaz*Sij1(I,K,J)*omgij5(I,K,J)**2.d0
     &        -deltax*Sij2(I,K,J)
     &                    *(0.5d0*(omgij4(I,K,J+1)+omgij4(I,K,J)))**2.d0
     &        -deltaz*Sij2(I,K,J)
     &                    *(0.5d0*(omgij6(I,K,J+1)+omgij6(I,K,J)))**2.d0
     &        -deltax*Sij3(I,K,J)*omgij5(I,K,J)**2.d0
     &        -deltay*Sij3(I,K,J)
     &                    *(0.5d0*(omgij6(I,K,J+1)+omgij6(I,K,J)))**2.d0
     &      -2.d0*deltaz*(0.5d0*(Sij4(I,K,J+1)+Sij4(I,K,J)))
     &            *omgij5(I,K,J)*(0.5d0*(omgij6(I,K,J+1)+omgij6(I,K,J)))
     &      +2.d0*deltay*Sij5(I,K,J)
     &           *(0.5d0*(omgij4(I,K,J+1)+omgij4(I,K,J)))
     &           *(0.5d0*(omgij6(I,K,J+1)+omgij6(I,K,J)))
     &      -2.d0*deltax*(0.5d0*(Sij6(I,K,J+1)+Sij6(I,K,J)))
     &           *(0.5d0*(omgij4(I,K,J+1)+omgij4(I,K,J)))*omgij5(I,K,J)
     & )

           IF (S1(I,K,J) <= 0.0d0) THEN
            S1(I,K,J)=0.0d0
           ELSE
            S1(I,K,J)=S1(I,K,J)/(Sij1(I,K,J)**2.d0
     &                          +Sij2(I,K,J)**2.d0
     &                          +Sij3(I,K,J)**2.d0
     &            +2.d0*(0.5d0*(Sij4(I,K,J+1)+Sij4(I,K,J)))**2.d0
     &            +2.d0*(0.5d0*(omgij4(I,K,J+1)+omgij4(I,K,J)))**2.d0
     &            +2.d0*Sij5(I,K,J)**2.d0
     &            +2.d0*omgij5(I,K,J)**2.d0
     &            +2.d0*(0.5d0*(Sij6(I,K,J+1)+Sij6(I,K,J)))**2.d0
     &            +2.d0*(0.5d0*(omgij6(I,K,J+1)+omgij6(I,K,J)))**2.d0)
           END IF

          END DO
        END DO
       END DO


! Compute max{0,-delta_k*(I3-I4)}/(I1-I2) at GY points and store in TEMP
       DO J=2,NY+1
!Set filter length (based on grid size) in y direction
!Based on constant Smag code        
         deltay=(DY(J)*2.d0)**2.d0
        DO K=0,NZP-1
          DO I=0,NXM

!First calculate delta_k*I3 and store it in TEMP.
            TEMP(I,K,J)=
     &      deltax*((Sij1(I,K,J)*DYF(j-1)+Sij1(I,K,J-1)*DYF(j))
     &                               /(2.d0*DY(j)))**3.d0
     &     +deltay*((Sij2(I,K,J)*DYF(j-1)+Sij2(I,K,J-1)*DYF(j))
     &                               /(2.d0*DY(j)))**3.d0
     &     +deltaz*((Sij3(I,K,J)*DYF(j-1)+Sij3(I,K,J-1)*DYF(j))
     &                               /(2.d0*DY(j)))**3.d0
     &     +(2.d0*deltax+deltay)
     &       *((Sij1(I,K,J)*DYF(j-1)+Sij1(I,K,J-1)*DYF(j))/(2.d0*DY(j)))
     &        *Sij4(I,K,J)**2.d0
     &     +(2.d0*deltax+deltaz)
     &       *((Sij1(I,K,J)*DYF(j-1)+Sij1(I,K,J-1)*DYF(j))/(2.d0*DY(j)))
     & *((Sij5(I,K,J)*DYF(j-1)+Sij5(I,K,J-1)*DYF(j))/(2.d0*DY(j)))**2.d0
     &     +(2.d0*deltay+deltax)
     &       *((Sij2(I,K,J)*DYF(j-1)+Sij2(I,K,J-1)*DYF(j))/(2.d0*DY(j)))
     &        *Sij4(I,K,J)**2.d0
     &     +(2.d0*deltay+deltaz)
     &       *((Sij2(I,K,J)*DYF(j-1)+Sij2(I,K,J-1)*DYF(j))/(2.d0*DY(j)))
     &        *Sij6(I,K,J)**2.d0
     &     +(2.d0*deltaz+deltax)*
     &       *((Sij3(I,K,J)*DYF(j-1)+Sij3(I,K,J-1)*DYF(j))/(2.d0*DY(j)))
     & *((Sij5(I,K,J)*DYF(j-1)+Sij5(I,K,J-1)*DYF(j))/(2.d0*DY(j)))**2.d0
     &     +(2.d0*deltaz+deltay)
     &       *((Sij3(I,K,J)*DYF(j-1)+Sij3(I,K,J-1)*DYF(j))/(2.d0*DY(j)))
     &        *Sij6(I,K,J)**2.d0
     &     +2.d0*(deltax+deltay+deltaz)*Sij4(I,K,J)*Sij6(I,K,J)
     &        *((Sij5(I,K,J)*DYF(j-1)+Sij5(I,K,J-1)*DYF(j))/(2.d0*DY(j))
     & )

!Then calculate -(I3-I4) = -TEMP+delta_k*I4 and store in TEMP.
            TEMP(I,K,J)=-TEMP(I,K,J)
     &     -deltay*((Sij1(I,K,J)*DYF(j-1)+Sij1(I,K,J-1)*DYF(j))
     &                               /(2.d0*DY(j)))
     &        *omgij4(I,K,J)**2.d0
     &     -deltaz*((Sij1(I,K,J)*DYF(j-1)+Sij1(I,K,J-1)*DYF(j))
     &                               /(2.d0*DY(j)))
     &        *((omgij5(I,K,J)*DYF(j-1)+omgij5(I,K,J-1)*DYF(j))
     &                               /(2.d0*DY(j)))**2.d0
     &     -deltax*((Sij2(I,K,J)*DYF(j-1)+Sij2(I,K,J-1)*DYF(j))
     &                               /(2.d0*DY(j)))
     &        *omgij4(I,K,J)**2.d0
     &     -deltaz*((Sij2(I,K,J)*DYF(j-1)+Sij2(I,K,J-1)*DYF(j))
     &                               /(2.d0*DY(j)))
     &        *omgij6(I,K,J)**2.d0
     &     -deltax*((Sij3(I,K,J)*DYF(j-1)+Sij3(I,K,J-1)*DYF(j))
     &                               /(2.d0*DY(j)))
     &        *((omgij5(I,K,J)*DYF(j-1)+omgij5(I,K,J-1)*DYF(j))
     &                               /(2.d0*DY(j)))**2.d0
     &     -deltay*((Sij3(I,K,J)*DYF(j-1)+Sij3(I,K,J-1)*DYF(j))
     &                               /(2.d0*DY(j)))
     &        *omgij6(I,K,J)**2.d0
     &     -2.d0*deltaz*Sij4(I,K,J)
     &        *((omgij5(I,K,J)*DYF(j-1)+omgij5(I,K,J-1)*DYF(j))
     &                               /(2.d0*DY(j)))
     &        *omgij6(I,K,J)
     &     +2.d0*deltay*((Sij5(I,K,J)*DYF(j-1)+Sij5(I,K,J-1)*DYF(j))
     &                               /(2.d0*DY(j)))
     &        *omgij4(I,K,J)*omgij6(I,K,J)
     &     -2.d0*deltax*Sij6(I,K,J)*omgij4(I,K,J)
     &         *((omgij5(I,K,J)*DYF(j-1)+omgij5(I,K,J-1)*DYF(j))
     &                               /(2.d0*DY(j)))


            IF (TEMP(I,K,J) <= 0.0d0) THEN
              TEMP(I,K,J)=0.0d0
            ELSE
              TEMP(I,K,J)=TEMP(I,K,J)/
     &              (((Sij1(I,K,J)*DYF(j-1)+Sij1(I,K,J-1)*DYF(j))
     &                               /(2.d0*DY(j)))**2.d0
     &              +((Sij2(I,K,J)*DYF(j-1)+Sij2(I,K,J-1)*DYF(j))
     &                               /(2.d0*DY(j)))**2.d0
     &              +((Sij3(I,K,J)*DYF(j-1)+Sij3(I,K,J-1)*DYF(j))
     &                               /(2.d0*DY(j)))**2.d0
     &            +2.d0*Sij4(I,K,J)**2.d0
     &            +2.d0*omgij4(I,K,J)**2.d0
     &            +2.d0*((Sij5(I,K,J)*DYF(j-1)+Sij5(I,K,J-1)*DYF(j))
     &                               /(2.d0*DY(j)))**2.d0
     &            +2.d0*((omgij5(I,K,J)*DYF(j-1)+omgij5(I,K,J-1)*DYF(j))
     &                               /(2.d0*DY(j)))**2.d0
     &            +2.d0*Sij6(I,K,J)**2.d0
     &            +2.d0*omgij6(I,K,J)**2.d0 )

            END IF

          END DO
        END DO
       END DO

! Now, compute max{0,-delta_k*(I3-I4)}/(I1-I2)*S_ij, storing in Sij
! First compute at GYF points 
       DO J=JSTART,JEND
        DO K=0,NZP-1
          DO I=0,NXM
            Sij1(I,K,J)=S1(I,K,J)*Sij1(I,K,J)
            Sij5(I,K,J)=S1(I,K,J)*Sij5(I,K,J)
! Sij2 is added through an implicit eddy viscosity
            Sij2(I,K,J)=0.d0
            Sij3(I,K,J)=S1(I,K,J)*Sij3(I,K,J)
          END DO
        END DO
       END DO


! Now, compute at max{0,-delta_k*(I3-I4)}/(I1-I2)*S_ij at GY points
       DO J=2,NY+1
        DO K=0,NZP-1
          DO I=0,NXM
! The terms dU1/dy and dU3/dy in CSij4(:,:,:) and CSij6(:,:,:) respectively
! are subtracted out from Sij here since they are treated implicitly
! in eddy viscosity terms
            Sij4(I,K,J)=TEMP(I,K,J)
     &        *(Sij4(I,K,J)-0.5*(U1(I,K,J)-U1(I,K,J-1))/DY(j))
            Sij6(I,K,J)=TEMP(I,K,J)
     &        *(Sij6(I,K,J)-0.5*(U3(I,K,J)-U3(I,K,J-1))/DY(j))
          END DO
        END DO
       END DO

! We now have max{0,-delta_k*(I3-I4)}/(I1-I2)*S_ij stored in Sij in Physical space

! Compute  -2.d0*C_AMD**2.d0 here for efficiency
       DO J=1,NY+1
! At GYF points:
! AMD (based off constant Smag)
        DELTA_YF(J)=-2.d0*C_AMD**2.d0
!     &     *(DX(1)*beta_sgs*DYF(J)*2.d0*DZ(1)*beta_sgs)**(2.d0/3.d0)
! Wall Damping
!        DELTA_YF(J)=
!     &    -2.d0*(0.1d0*(1.d0-exp((-GYF(J)/(NU*25.d0))**3.d0)))**2.d0
!     &            *(DX(1)*beta_sgs*DYF(J)*2.d0*DZ(1)*beta_sgs)**(2.d0/3.d0)

       END DO

       DO J=1,NY+1
! At GY points:
! AMD (based off Constant Smagorinsky)
        DELTA_Y(J)=-2.d0*C_AMD**2.d0
!     &        *(DX(1)*beta_sgs*DY(J)*2.d0*DZ(1)*beta_sgs)**(2.d0/3.d0)
! Wall Damping
!        DELTA_Y(J)=
!     &    -2.d0*(0.1d0*(1.d0-exp((-GY(J)/(NU*25.d0))**3.d0)))**2.d0
!     &            *(DX(1)*beta_sgs*DY(J)*2.d0*DZ(1)*beta_sgs)**(2.d0/3.d0)
       END DO


! Get the eddy viscosity at GY points
! NU_T = (C_S^2 * DELTA^2)*|S|
       DO J=2,NY
        DO K=0,NZP-1
          DO I=0,NXM
           NU_T(I,K,J)=-0.5d0*DELTA_Y(J)*TEMP(I,K,J)
          END DO
        END DO
       END DO

! Now that we have calculated NU_T, set the value at the ghost cells
! by sharing with neighboring processes.  This subroutine also sets
! the value of NU_T to zero at both walls
      CALL GHOST_LES_MPI

! Convert the stress tensor to Fourier space


      CALL FFT_XZ_TO_FOURIER(Sij1,CSij1,0,NY+1)
! Sij2 is added through an implicit eddy viscosity
!      CALL FFT_XZ_TO_FOURIER(Sij2,CSij2,0,NY+1)
      CALL FFT_XZ_TO_FOURIER(Sij3,CSij3,0,NY+1)
      CALL FFT_XZ_TO_FOURIER(Sij4,CSij4,0,NY+1)
      CALL FFT_XZ_TO_FOURIER(Sij5,CSij5,0,NY+1)
      CALL FFT_XZ_TO_FOURIER(Sij6,CSij6,0,NY+1)

! Now, compute TAU, store in the corresponging Sij
       DO J=1,NY+1
        DO K=0,TNKZ
          DO I=0,NXP-1
            CSij1(I,K,J)=DELTA_YF(J)*CSij1(I,K,J)
            CSij5(I,K,J)=DELTA_YF(J)*CSij5(I,K,J)
! CSij2(:,:,:) is added through an implicit eddy viscosity
!            CSij2(I,K,J)=DELTA_YF(J)*CSij2(I,K,J)
            CSij3(I,K,J)=DELTA_YF(J)*CSij3(I,K,J)
          END DO
        END DO
       END DO
       DO J=2,NY+1
        DO K=0,TNKZ
          DO I=0,NXP-1
            CSij4(I,K,J)=DELTA_Y(J)*CSij4(I,K,J)
            CSij6(I,K,J)=DELTA_Y(J)*CSij6(I,K,J)
          END DO
        END DO
       END DO


! tau_ij is now contained in CSij in Fourier space


      else if (LES_MODEL_TYPE.EQ.5) then
!     Anisotrophic minimum-dissipation model Verstappen

! First, compute the rate of strain tensor S_ij

      call compute_strain_chan

! Then compute the rate of strain tensor invariant SI_ij

      call compute_strain_chan_invariant

! Then compute the rate of rotation tensor invariant omega_ij 

      call compute_rotation_chan_invariant

!---------------------------------------------------------!
! Compute the mean velocity and TH
! Copied from Sr. SAVE_STATS_CHAN 
! Get the mean value of the velocities
       IF (RANKZ.EQ.0) THEN
         ume=dble(CU1(0,0,:))
         vme=dble(CU2(0,0,:))
         wme=dble(CU3(0,0,:))
c         DO n=1,N_TH
c            thme(:,n)=dble(CTH(0,0,:,n)) 
c         END DO
       END IF
      CALL MPI_BCAST(ume,NY+2,MPI_DOUBLE_PRECISION,0,
     &     MPI_COMM_Z,ierror)
      CALL MPI_BCAST(vme,NY+2,MPI_DOUBLE_PRECISION,0,
     &     MPI_COMM_Z,ierror)
      CALL MPI_BCAST(wme,NY+2,MPI_DOUBLE_PRECISION,0,
     &     MPI_COMM_Z,ierror)
c      IF (N_TH.GT.0) CALL MPI_BCAST(thme,(NY+2)*N_TH,
c     &     MPI_DOUBLE_PRECISION,0,MPI_COMM_Z,ierror)
!!!   We do not need thme here since local theta is needed


       N=1 ! Only buoyancy/temperature field
       CTMP(:,:,:)=CTH(:,:,:,N)
       CALL FFT_XZ_TO_PHYSICAL(CTMP,TH_real,0,NY+1)

! Compute dthetady in fourier space
c       DO J=1,NY+1
c        DO K=0,TNKZ
c          DO I=0,NXP-1
! Cdthetady at GY points
c            Cdthetady(I,K,J,N)=(CTH(I,K,J,N)-CTH(I,K,J-1,N))/DY(J)
! Cdthetady at GYF points
c            Cdthetady_YF(I,K,J,N)=(CTH(I,K,J+1,N)-CTH(I,K,J-1,N))
c     &                             /(GYF(J+1)-GYF(J-1)) 
c          END DO
c        END DO
c       END DO

! dthetady at GY points in physical space
c       CS1(:,:,:)=Cdthetady(:,:,:,N)
c       CALL FFT_XZ_TO_PHYSICAL(CS1,S1,0,NY+1)
c       dthetady(:,:,:,N)=S1(:,:,:)
! dthetady at GYF points in physical space
c       CS1(:,:,:)=Cdthetady_YF(:,:,:,N)
c       CALL FFT_XZ_TO_PHYSICAL(CS1,S1,0,NY+1)
c       dthetady_YF(:,:,:,N)=S1(:,:,:)

!=========================================================!

! Convert the velocity to physical space
      call FFT_XZ_TO_PHYSICAL(CU1,U1,0,NY+1)
      call FFT_XZ_TO_PHYSICAL(CU2,U2,0,NY+1)
      call FFT_XZ_TO_PHYSICAL(CU3,U3,0,NY+1)

! Compute max{0,-(I3-I4)}/(I1-I2) at GYF points, store in S1
! Interpolation to GYF points is easy since by definition
! GYF points are exactly midway between neighboring GY points
       DO J=JSTART,JEND
        DO K=0,NZP-1
          DO I=0,NXM

!First calculate I3 and store it in S1.  
            S1(I,K,J)=
     &               SIij1(I,K,J)**3.d0
     &              +SIij2(I,K,J)**3.d0
     &              +SIij3(I,K,J)**3.d0
     &    +3.d0*SIij1(I,K,J)*(0.5d0*(SIij4(I,K,J+1)+SIij4(I,K,J)))**2.d0
     &    +3.d0*SIij1(I,K,J)*SIij5(I,K,J)**2.d0
     &    +3.d0*SIij2(I,K,J)*(0.5d0*(SIij4(I,K,J+1)+SIij4(I,K,J)))**2.d0
     &    +3.d0*SIij2(I,K,J)*(0.5d0*(SIij6(I,K,J+1)+SIij6(I,K,J)))**2.d0
     &    +3.d0*SIij3(I,K,J)*SIij5(I,K,J)**2.d0
     &    +3.d0*SIij3(I,K,J)*(0.5d0*(SIij6(I,K,J+1)+SIij6(I,K,J)))**2.d0
     &    +6.d0*(0.5d0*(SIij4(I,K,J+1)+SIij4(I,K,J)))*SIij5(I,K,J)
     &             *(0.5d0*(SIij6(I,K,J+1)+SIij6(I,K,J)))

!Then calculate -(I3-I4) = -S1+I4 and store in S1.

            S1(I,K,J)=(-S1(I,K,J)
     &      -SIij1(I,K,J)*(0.5d0*(omgij4(I,K,J+1)+omgij4(I,K,J)))**2.d0
     &      -SIij1(I,K,J)*omgij5(I,K,J)**2.d0
     &      -SIij2(I,K,J)*(0.5d0*(omgij4(I,K,J+1)+omgij4(I,K,J)))**2.d0
     &      -SIij2(I,K,J)*(0.5d0*(omgij6(I,K,J+1)+omgij6(I,K,J)))**2.d0
     &      -SIij3(I,K,J)*omgij5(I,K,J)**2.d0
     &      -SIij3(I,K,J)*(0.5d0*(omgij6(I,K,J+1)+omgij6(I,K,J)))**2.d0
     &      -2.d0*(0.5d0*(SIij4(I,K,J+1)+SIij4(I,K,J)))*omgij5(I,K,J)
     &           *(0.5d0*(omgij6(I,K,J+1)+omgij6(I,K,J)))!
     &      +2.d0*SIij5(I,K,J)*(0.5d0*(omgij4(I,K,J+1)+omgij4(I,K,J)))
     &           *(0.5d0*(omgij6(I,K,J+1)+omgij6(I,K,J)))
     &      -2.d0*(0.5d0*(SIij6(I,K,J+1)+SIij6(I,K,J)))
     &           *(0.5d0*(omgij4(I,K,J+1)+omgij4(I,K,J)))*omgij5(I,K,J)
     & )

        IF (S1(I,K,J) <= 0.0d0) THEN
          S1(I,K,J)=0.0d0   
        ELSE
          S1(I,K,J)=S1(I,K,J)/(SIij1(I,K,J)**2.d0
     &                          +SIij2(I,K,J)**2.d0 
     &                          +SIij3(I,K,J)**2.d0 
     &            +2.d0*(0.5d0*(SIij4(I,K,J+1)+SIij4(I,K,J)))**2.d0 
     &            +2.d0*(0.5d0*(omgij4(I,K,J+1)+omgij4(I,K,J)))**2.d0 
     &            +2.d0*SIij5(I,K,J)**2.d0
     &            +2.d0*omgij5(I,K,J)**2.d0
     &            +2.d0*(0.5d0*(SIij6(I,K,J+1)+SIij6(I,K,J)))**2.d0     
     &            +2.d0*(0.5d0*(omgij6(I,K,J+1)+omgij6(I,K,J)))**2.d0)
        END IF

          END DO
        END DO
       END DO



! Compute max{0,-(I3-I4)}/(I1-I2) at GY points and store in TEMP
       DO J=2,NY+1
        DO K=0,NZP-1
          DO I=0,NXM

!First calculate I3 and store it in TEMP.
            TEMP(I,K,J)=
     &                ((SIij1(I,K,J)*DYF(j-1)+SIij1(I,K,J-1)*DYF(j))
     &                               /(2.d0*DY(j)))**3.d0
     &               +((SIij2(I,K,J)*DYF(j-1)+SIij2(I,K,J-1)*DYF(j))
     &                               /(2.d0*DY(j)))**3.d0
     &               +((SIij3(I,K,J)*DYF(j-1)+SIij3(I,K,J-1)*DYF(j))
     &                               /(2.d0*DY(j)))**3.d0
     &          +3.d0*((SIij1(I,K,J)*DYF(j-1)+SIij1(I,K,J-1)*DYF(j))
     &                               /(2.d0*DY(j)))
     &               *SIij4(I,K,J)**2.d0
     &          +3.d0*((SIij1(I,K,J)*DYF(j-1)+SIij1(I,K,J-1)*DYF(j))
     &                               /(2.d0*DY(j)))
     &               *((SIij5(I,K,J)*DYF(j-1)+SIij5(I,K,J-1)*DYF(j))
     &                               /(2.d0*DY(j)))**2.d0           
     &          +3.d0*((SIij2(I,K,J)*DYF(j-1)+SIij2(I,K,J-1)*DYF(j))
     &                               /(2.d0*DY(j)))
     &               *SIij4(I,K,J)**2.d0
     &          +3.d0*((SIij2(I,K,J)*DYF(j-1)+SIij2(I,K,J-1)*DYF(j))
     &                               /(2.d0*DY(j)))
     &               *SIij6(I,K,J)**2.d0
     &          +3.d0*((SIij3(I,K,J)*DYF(j-1)+SIij3(I,K,J-1)*DYF(j))
     &                               /(2.d0*DY(j)))
     &               *((SIij5(I,K,J)*DYF(j-1)+SIij5(I,K,J-1)*DYF(j))
     &                               /(2.d0*DY(j)))**2.d0 
     &          +3.d0*((SIij3(I,K,J)*DYF(j-1)+SIij3(I,K,J-1)*DYF(j))
     &                               /(2.d0*DY(j)))
     &               *SIij6(I,K,J)**2.d0
     &          +6.d0*SIij4(I,K,J)*SIij6(I,K,J)
     &               *((SIij5(I,K,J)*DYF(j-1)+SIij5(I,K,J-1)*DYF(j))
     &                               /(2.d0*DY(j)))

!Then calculate -(I3-I4) = -TEMP+I4 and store in TEMP.
            TEMP(I,K,J)=-TEMP(I,K,J)
     &         -((SIij1(I,K,J)*DYF(j-1)+SIij1(I,K,J-1)*DYF(j))
     &                               /(2.d0*DY(j)))
     &         *omgij4(I,K,J)**2.d0
     &        -((SIij1(I,K,J)*DYF(j-1)+SIij1(I,K,J-1)*DYF(j))
     &                               /(2.d0*DY(j)))
     &         *((omgij5(I,K,J)*DYF(j-1)+omgij5(I,K,J-1)*DYF(j))
     &                               /(2.d0*DY(j)))**2.d0
     &        -((SIij2(I,K,J)*DYF(j-1)+SIij2(I,K,J-1)*DYF(j))
     &                               /(2.d0*DY(j))) 
     &         *omgij4(I,K,J)**2.d0
     &        -((SIij2(I,K,J)*DYF(j-1)+SIij2(I,K,J-1)*DYF(j))
     &                               /(2.d0*DY(j))) 
     &         *omgij6(I,K,J)**2.d0
     &        -((SIij3(I,K,J)*DYF(j-1)+SIij3(I,K,J-1)*DYF(j))
     &                               /(2.d0*DY(j)))
     &         *((omgij5(I,K,J)*DYF(j-1)+omgij5(I,K,J-1)*DYF(j))
     &                               /(2.d0*DY(j)))**2.d0
     &        -((SIij3(I,K,J)*DYF(j-1)+SIij3(I,K,J-1)*DYF(j))
     &                               /(2.d0*DY(j)))
     &         *omgij6(I,K,J)**2.d0
     &        -2.d0*SIij4(I,K,J)
     &         *((omgij5(I,K,J)*DYF(j-1)+omgij5(I,K,J-1)*DYF(j))
     &                               /(2.d0*DY(j)))
     &         *omgij6(I,K,J)
     &        +2.d0*((SIij5(I,K,J)*DYF(j-1)+SIij5(I,K,J-1)*DYF(j))
     &                               /(2.d0*DY(j)))
     &         *omgij4(I,K,J)*omgij6(I,K,J)
     &        -2.d0*SIij6(I,K,J)*omgij4(I,K,J)
     &         *((omgij5(I,K,J)*DYF(j-1)+omgij5(I,K,J-1)*DYF(j))
     &                               /(2.d0*DY(j)))

        IF (TEMP(I,K,J) <= 0.0d0) THEN
          TEMP(I,K,J)=0.0d0
        ELSE
        TEMP(I,K,J)=TEMP(I,K,J)/
     &              (((SIij1(I,K,J)*DYF(j-1)+SIij1(I,K,J-1)*DYF(j))
     &                               /(2.d0*DY(j)))**2.d0
     &              +((SIij2(I,K,J)*DYF(j-1)+SIij2(I,K,J-1)*DYF(j))
     &                               /(2.d0*DY(j)))**2.d0
     &              +((SIij3(I,K,J)*DYF(j-1)+SIij3(I,K,J-1)*DYF(j))
     &                               /(2.d0*DY(j)))**2.d0
     &            +2.d0*SIij4(I,K,J)**2.d0
     &            +2.d0*omgij4(I,K,J)**2.d0
     &            +2.d0*((SIij5(I,K,J)*DYF(j-1)+SIij5(I,K,J-1)*DYF(j))
     &                               /(2.d0*DY(j)))**2.d0
     &            +2.d0*((omgij5(I,K,J)*DYF(j-1)+omgij5(I,K,J-1)*DYF(j))
     &                               /(2.d0*DY(j)))**2.d0
     &            +2.d0*SIij6(I,K,J)**2.d0
     &            +2.d0*omgij6(I,K,J)**2.d0 )

        END IF

          END DO
        END DO
       END DO


! Now, compute max{0,-(I3-I4)}/(I1-I2)*S_ij, storing in Sij
! First compute at GYF points 
       DO J=JSTART,JEND
        DO K=0,NZP-1
          DO I=0,NXM
            Sij1(I,K,J)=S1(I,K,J)*Sij1(I,K,J)
            Sij5(I,K,J)=S1(I,K,J)*Sij5(I,K,J)
! Sij2 is added through an implicit eddy viscosity
            Sij2(I,K,J)=0.d0
            Sij3(I,K,J)=S1(I,K,J)*Sij3(I,K,J)
          END DO
        END DO
       END DO


! Now, compute at max{0,-(I3-I4)}/(I1-I2)*S_ij at GY points
       DO J=2,NY+1
        DO K=0,NZP-1
          DO I=0,NXM
! The terms dU1/dy and dU3/dy in CSij4(:,:,:) and CSij6(:,:,:) respectively
! are subtracted out from Sij here since they are treated implicitly
! in eddy viscosity terms
            Sij4(I,K,J)=TEMP(I,K,J)
     &        *(Sij4(I,K,J)-0.5*(U1(I,K,J)-U1(I,K,J-1))/DY(j))
            Sij6(I,K,J)=TEMP(I,K,J)
     &        *(Sij6(I,K,J)-0.5*(U3(I,K,J)-U3(I,K,J-1))/DY(j))
          END DO
        END DO
       END DO

! We now have max{0,-(I3-I4)}/(I1-I2)*S_ij stored in Sij in Physical space




! Compute the filter lengthscale
! Absorb -2.d0*C_AMD**2.d0 here for efficiency
       DO J=1,NY+1
! At GYF points:
! AMD (based off constant Smag)
        DELTA_YF(J)=-2.d0*C_AMD**2.d0
     &     *3/(1/(DX(1)*beta_sgs)**2.d0+1/(DYF(J)*2.d0)**2.d0
     &        +1/(DZ(1)*beta_sgs)**2.d0)
!     &     *(DX(1)*beta_sgs*DYF(J)*2.d0*DZ(1)*beta_sgs)**(2.d0/3.d0)
! Wall Damping
!        DELTA_YF(J)=
!     &    -2.d0*(0.1d0*(1.d0-exp((-GYF(J)/(NU*25.d0))**3.d0)))**2.d0
!     &            *(DX(1)*beta_sgs*DYF(J)*2.d0*DZ(1)*beta_sgs)**(2.d0/3.d0)

       END DO

       DO J=1,NY+1
! At GY points:
! AMD (based off Constant Smagorinsky)
        DELTA_Y(J)=-2.d0*C_AMD**2.d0
     &     *3/(1/(DX(1)*beta_sgs)**2.d0+1/(DY(J)*2.d0)**2.d0
     &        +1/(DZ(1)*beta_sgs)**2.d0)
!     &        *(DX(1)*beta_sgs*DY(J)*2.d0*DZ(1)*beta_sgs)**(2.d0/3.d0)
! Wall Damping
!        DELTA_Y(J)=
!     &    -2.d0*(0.1d0*(1.d0-exp((-GY(J)/(NU*25.d0))**3.d0)))**2.d0
!     &            *(DX(1)*beta_sgs*DY(J)*2.d0*DZ(1)*beta_sgs)**(2.d0/3.d0)
       END DO

! Get the eddy viscosity at GY points
! NU_T = (C_S^2 * DELTA^2)*|S|
       N=1
       DO J=2,NY
        DO K=0,NZP-1
          DO I=0,NXM

           DELTA_TMP_Y = -2.d0*C_AMD**2.d0
     &     *3/(1/(DX(1)*beta_sgs)**2.d0+1/(DY(J)*2.d0)**2.d0
     &        +1/(DZ(1)*beta_sgs)**2.d0 
     &        + C_lb*MAX(0.0,(TH_real(I,K,J)-TH_real(I,K,J-1))/DY(J))/(
     &          ((U1(I,K,J)-ume(j))**2.d0*DYF(J-1) 
     &            +(U1(I,K,J-1)-ume(J-1))**2.d0*DYF(J))/(2.d0*DY(J))
     &         +((U3(I,K,J)-wme(j))**2.d0*DYF(J-1)
     &            +(U3(I,K,J-1)-wme(J-1))**2.d0*DYF(J))/(2.d0*DY(J))
     &         +((U2(I,K,J)-vme(J))**2.d0)+EPS_DELTA )
     &         )

           NU_T(I,K,J)=-0.5d0*DELTA_TMP_Y*TEMP(I,K,J)
c           NU_T(I,K,J)=-0.5d0*DELTA_Y(J)*TEMP(I,K,J)

          END DO
        END DO
       END DO

! Now that we have calculated NU_T, set the value at the ghost cells
! by sharing with neighboring processes.  This subroutine also sets
! the value of NU_T to zero at both walls
      CALL GHOST_LES_MPI

! Convert the stress tensor to Fourier space


      CALL FFT_XZ_TO_FOURIER(Sij1,CSij1,0,NY+1)
! Sij2 is added through an implicit eddy viscosity
!      CALL FFT_XZ_TO_FOURIER(Sij2,CSij2,0,NY+1)
      CALL FFT_XZ_TO_FOURIER(Sij3,CSij3,0,NY+1)
      CALL FFT_XZ_TO_FOURIER(Sij4,CSij4,0,NY+1)
      CALL FFT_XZ_TO_FOURIER(Sij5,CSij5,0,NY+1)
      CALL FFT_XZ_TO_FOURIER(Sij6,CSij6,0,NY+1)

! Now, compute TAU, store in the corresponging Sij
       DO J=1,NY+1
        DO K=0,TNKZ
          DO I=0,NXP-1
           IF (J.EQ.NY+1) THEN
            DELTA_TMP_YF = -2.d0*C_AMD**2.d0
     &      *3/(1/(DX(1)*beta_sgs)**2.d0+1/(DYF(J)*2.d0)**2.d0
     &        +1/(DZ(1)*beta_sgs)**2.d0
     &        +C_lb*MAX(0.0,(TH_real(I,K,J)-TH_real(I,K,J-1))/
     &              (GYF(J)-GYF(J-1)))/(
     &          (U1(I,K,J)-ume(J))**2.d0+(U3(I,K,J)-wme(J))**2.d0
     &              + 0.5*((U2(I,K,J)-vme(J))**2.d0 +
     &                     (U2(I,K,J)-vme(J))**2.d0 )+EPS_DELTA)
     &         )
           ELSE
            DELTA_TMP_YF = -2.d0*C_AMD**2.d0
     &      *3/(1/(DX(1)*beta_sgs)**2.d0+1/(DYF(J)*2.d0)**2.d0
     &        +1/(DZ(1)*beta_sgs)**2.d0
     &        +C_lb*MAX(0.0,(TH_real(I,K,J+1)-TH_real(I,K,J-1))/
     &              (GYF(J+1)-GYF(J-1)))/(
     &          (U1(I,K,J)-ume(J))**2.d0+(U3(I,K,J)-wme(J))**2.d0
     &              + 0.5*((U2(I,K,J)-vme(J))**2.d0 +
     &                     (U2(I,K,J+1)-vme(J+1))**2.d0 )+EPS_DELTA)
     &         )
           ENDIF
c            CSij1(I,K,J)=DELTA_YF(J)*CSij1(I,K,J)
            CSij1(I,K,J)=DELTA_TMP_YF*CSij1(I,K,J)
c            CSij5(I,K,J)=DELTA_YF(J)*CSij5(I,K,J)
            CSij5(I,K,J)=DELTA_TMP_YF*CSij5(I,K,J)
! CSij2(:,:,:) is added through an implicit eddy viscosity
!            CSij2(I,K,J)=DELTA_YF(J)*CSij2(I,K,J)
c            CSij3(I,K,J)=DELTA_YF(J)*CSij3(I,K,J)
            CSij3(I,K,J)=DELTA_TMP_YF*CSij3(I,K,J)

          END DO
        END DO
       END DO
       DO J=2,NY+1
        DO K=0,TNKZ
          DO I=0,NXP-1

           DELTA_TMP_Y = -2.d0*C_AMD**2.d0
     &     *3/(1/(DX(1)*beta_sgs)**2.d0+1/(DY(J)*2.d0)**2.d0
     &        +1/(DZ(1)*beta_sgs)**2.d0
     &        +C_lb*MAX(0.0,(TH_real(I,K,J)-TH_real(I,K,J-1))/DY(J))/(
     &          ((U1(I,K,J)-ume(j))**2.d0*DYF(J-1)
     &            +(U1(I,K,J-1)-ume(J-1))**2.d0*DYF(J))/(2.d0*DY(J))
     &         +((U3(I,K,J)-wme(j))**2.d0*DYF(J-1)
     &            +(U3(I,K,J-1)-wme(J-1))**2.d0*DYF(J))/(2.d0*DY(J))
     &         +((U2(I,K,J)-vme(J))**2.d0)+EPS_DELTA )
     &         )



c            CSij4(I,K,J)=DELTA_Y(J)*CSij4(I,K,J)
            CSij4(I,K,J)=DELTA_TMP_Y*CSij4(I,K,J)
c            CSij6(I,K,J)=DELTA_Y(J)*CSij6(I,K,J)
            CSij6(I,K,J)=DELTA_TMP_Y*CSij6(I,K,J)
          END DO
        END DO
       END DO


! tau_ij is now contained in CSij in Fourier space

      end if

! Now, add the subgrid scale forcing to CFi
! (This includes the subgrid scale stress as an explicit R-K term


      DO J=J1,J2
        DO K=0,TNKZ
          DO I=0,NXP-1
            CF1(I,K,J)=CF1(I,K,J)
     &                -CIKX(I)*CSij1(I,K,J)
     &                -(CSij4(I,K,J+1)-CSij4(I,K,J))/DYF(j)
     &                -CIKZ(K)*CSij5(I,K,J)
            CF3(I,K,J)=CF3(I,K,J)
     &                -CIKX(I)*CSij5(I,K,J)
     &                -(CSij6(I,K,J+1)-CSij6(I,K,J))/DYF(J)
     &                -CIKZ(K)*CSij3(I,K,J)   
          END DO
        END DO
      END DO
      DO J=2,NY
       DO K=0,TNKZ
         DO I=0,NXP-1
           CF2(I,K,J)=CF2(I,K,J)
     &                -CIKX(I)*CSij4(I,K,J)
! Sij2 is added through an implict eddy viscosity
!     &                -(CSij2(I,K,J)-CSij2(I,K,J-1))/DY(j)
     &                -CIKZ(K)*CSij6(I,K,J)
          END DO
        END DO
      END DO

! Add back the Cross-terms using an explicit treatment
      DO K=0,NZP-1
        DO I=0,NXM
          DO J=2,NY
            TEMP(I,K,J)=NU_T(I,K,J)*(U1(I,K,J)-U1(I,K,J-1))/DY(J)
          END DO
        END DO
      END DO
      CALL FFT_XZ_TO_FOURIER(TEMP,CTEMP,2,NY)
      DO K=0,TNKZ
        DO I=0,NXP-1
          DO J=2,NY
            CF2(I,K,J)=CF2(I,K,J)+CIKX(I)*CTEMP(I,K,J)
          END DO
        END DO
      END DO
      DO K=0,NZP-1
        DO I=0,NXM
          DO J=2,NY
            TEMP(I,K,J)=NU_T(I,K,J)*(U3(I,K,J)-U3(I,K,J-1))/DY(J)
          END DO
        END DO
      END DO
      CALL FFT_XZ_TO_FOURIER(TEMP,CTEMP,2,NY)
      DO K=0,TNKZ
        DO I=0,NXP-1
          DO J=2,NY
            CF2(I,K,J)=CF2(I,K,J)+CIKZ(K)*CTEMP(I,K,J)
          END DO
        END DO
      END DO

! Periodically, output mean quantities
      IF ((MOD(TIME_STEP,SAVE_STATS_INT).EQ.0).AND.(RK_STEP.EQ.1)) THEN
! Get plane averages
        do J=0,NY+1
          S1_mean(J)=0.d0
          NU_T_mean(J)=0.d0
          do I=0,NXM
          do K=0,NZP-1
            S1_mean(J)=S1_mean(J)+S1(I,K,J)
            NU_T_mean(J)=NU_T_mean(J)+NU_T(I,K,J)
          end do
          end do
        end do
! compute part of sgs term
! U1*F1+U2*F2+U3*F3
! compute part of CF1
      DO J=J1,J2
        DO K=0,TNKZ
          DO I=0,NXP-1
            CTEMP(I,K,J)=
     &                -CIKX(I)*CSij1(I,K,J)
     &                -(CSij4(I,K,J+1)-CSij4(I,K,J))/DYF(J)
     &                -CIKZ(K)*CSij5(I,K,J)
          END DO
        END DO
      END DO
        CALL FFT_XZ_TO_PHYSICAL(CTEMP,TEMP,0,NY+1)
        do J=0,NY+1
          EPS_SGS1_MEAN(J)=0.d0
          do I=0,NXM
          do K=0,NZP-1
            EPS_SGS1_MEAN(J)=EPS_SGS1_MEAN(J)+TEMP(I,K,J)*U1(I,K,J)
          end do
          end do
        end do
! compute part of CF3
      DO J=J1,J2
        DO K=0,TNKZ
          DO I=0,NXP-1
            CTEMP(I,K,J)=
     &                -CIKX(I)*CSij5(I,K,J)
     &                -(CSij6(I,K,J+1)-CSij6(I,K,J))/DYF(J)
     &                -CIKZ(K)*CSij3(I,K,J)
          END DO
        END DO
      END DO
        CALL FFT_XZ_TO_PHYSICAL(CTEMP,TEMP,0,NY+1)
        do J=0,NY+1
          do I=0,NXM
          do K=0,NZP-1
            EPS_SGS1_MEAN(J)=EPS_SGS1_MEAN(J)+TEMP(I,K,J)*U3(I,K,J)
          end do
          end do
        end do
! compute part of CF2
      DO J=2,NY
       DO K=0,TNKZ
         DO I=0,NXP-1
           CTEMP(I,K,J)=
     &                -CIKX(I)*CSij4(I,K,J)
! Sij2 is added through an implict eddy viscosity
!     &                -(CSij2(I,K,J)-CSij2(I,K,J-1))/DY(j)
     &                -CIKZ(K)*CSij6(I,K,J)
          END DO
        END DO
      END DO
        CALL FFT_XZ_TO_PHYSICAL(CTEMP,TEMP,0,NY+1)
        do J=0,NY+1
          do I=0,NXM
          do K=0,NZP-1
            EPS_SGS1_MEAN(J)=EPS_SGS1_MEAN(J)+TEMP(I,K,J)*U2(I,K,J)
          end do
          end do
        end do
! Add back the Cross-terms using an explicit treatment
! compute second part of CF2
      DO K=0,NZP-1
        DO I=0,NXM
          DO J=2,NY
            TEMP(I,K,J)=NU_T(I,K,J)*(U1(I,K,J)-U1(I,K,J-1))/DY(J)
          END DO
        END DO
      END DO
      CALL FFT_XZ_TO_FOURIER(TEMP,CTEMP,2,NY)
      DO K=0,TNKZ
        DO I=0,NXP-1
          DO J=2,NY
            CTEMP(I,K,J)=CIKX(I)*CTEMP(I,K,J)
          END DO
        END DO
      END DO
        CALL FFT_XZ_TO_PHYSICAL(CTEMP,TEMP,0,NY+1)
        do J=0,NY+1
          do I=0,NXM
          do K=0,NZP-1
            EPS_SGS1_MEAN(J)=EPS_SGS1_MEAN(J)+TEMP(I,K,J)*U2(I,K,J)
          end do
          end do
        end do
      DO K=0,NZP-1
        DO I=0,NXM
          DO J=2,NY
            TEMP(I,K,J)=NU_T(I,K,J)*(U3(I,K,J)-U3(I,K,J-1))/DY(J)
          END DO
        END DO
      END DO
      CALL FFT_XZ_TO_FOURIER(TEMP,CTEMP,2,NY)
      DO K=0,TNKZ
        DO I=0,NXP-1
          DO J=2,NY
            CTEMP(I,K,J)=CIKZ(K)*CTEMP(I,K,J)
          END DO
        END DO
      END DO
        CALL FFT_XZ_TO_PHYSICAL(CTEMP,TEMP,0,NY+1)
        do J=0,NY+1
          do I=0,NXM
          do K=0,NZP-1
            EPS_SGS1_MEAN(J)=EPS_SGS1_MEAN(J)+TEMP(I,K,J)*U2(I,K,J)
          end do
          end do
        end do
      call mpi_allreduce(mpi_in_place,S1_mean,NY+2,MPI_DOUBLE_PRECISION,
     &     MPI_SUM,MPI_COMM_Z,ierror)
      call mpi_allreduce(mpi_in_place,NU_T_mean,NY+2
     &    ,MPI_DOUBLE_PRECISION,
     &     MPI_SUM,MPI_COMM_Z,ierror)
      call mpi_allreduce(mpi_in_place,EPS_SGS1_MEAN,NY+2
     &    ,MPI_DOUBLE_PRECISION,
     &     MPI_SUM,MPI_COMM_Z,ierror)

        do j=0,NY+1
          S1_mean(j)=S1_mean(j)/dble(NX*NZ)
          NU_T_mean(j)=NU_T_mean(j)/dble(NX*NZ)
          EPS_SGS1_MEAN(j)=EPS_SGS1_MEAN(j)/dble(NX*NZ)
        end do

#ifdef HDF5
      FNAME='mean_les.h5'

      gname='time'
      call WriteHDF5_real(FNAME,gname,TIME)

      IF (RANKZ.EQ.0) THEN

      gname='gyf'
      Diag=gyf(1:NY)
      call WriteStatH5(FNAME,gname,Diag)

      gname='nu_sgs'
      Diag=NU_T_mean(1:NY)
      call WriteStatH5(FNAME,gname,Diag)

      gname='eps_sgs1'
      Diag=EPS_SGS1_MEAN(1:NY)
      call WriteStatH5(FNAME,gname,Diag)

      END IF

#else
! Here we aren't using HDF5, so write to text file
      IF (RANKZ.EQ.0) THEN
      IF (USE_MPI) THEN
        FNAME='mean_les'//trim(MPI_IO_NUM)//'.txt'
      ELSE
        FNAME='mean_les.txt'
      END IF
      open(42,file=FNAME,form='formatted',status='unknown')
        write(42,*) TIME_STEP,TIME,DELTA_T
        do j=1,NY
          write(42,420) j,GYF(J),
     &      NU_T_mean(J),EPS_SGS1_MEAN(J)
        end do
420     format(I3,' ',2(F30.20,' '))
      END IF
#endif

      END IF

999   continue

      RETURN
      END

      subroutine compute_strain_chan
C This subroutine computes S_ij for the filtered velocity field
C The input velocity field should be in fourier space in the periodic
C directions.
C For use in the LES model in channel flow (2 periodic directions)
      include 'header'
      include 'header_les'

      integer I,J,K,ij

      DO J=1,NY
        DO K=0,TNKZ
          DO I=0,NXP-1
            CSij1(I,K,J)=CIKX(I)*CU1(I,K,J)
            CSij2(I,K,J)=(CU2(I,K,J+1)-CU2(I,K,J))/DYF(j)
            CSij3(I,K,J)=CIKZ(K)*CU3(I,K,J)
            CSij5(I,K,J)=0.5d0*(CIKZ(K)*CU1(I,K,J)
     &                  +CIKX(I)*CU3(I,K,J))
          END DO
        END DO
      END DO
      DO J=1,NY+1
        DO K=0,TNKZ
          DO I=0,NXP-1
            CSij4(I,K,J)=0.5d0*((CU1(I,K,J)-CU1(I,K,J-1))/DY(j)
     &                          +CIKX(I)*CU2(I,K,J) )
            CSij6(I,K,J)=0.5d0*(CIKZ(K)*CU2(I,K,J)
     &                         +(CU3(I,K,J)-CU3(I,K,J-1))/DY(j) )
          END DO
        END DO
      END DO


       CALL FFT_XZ_TO_PHYSICAL(CSij1,Sij1,0,NY+1)
       CALL FFT_XZ_TO_PHYSICAL(CSij2,Sij2,0,NY+1)
       CALL FFT_XZ_TO_PHYSICAL(CSij3,Sij3,0,NY+1)
       CALL FFT_XZ_TO_PHYSICAL(CSij4,Sij4,0,NY+1)
       CALL FFT_XZ_TO_PHYSICAL(CSij5,Sij5,0,NY+1)
       CALL FFT_XZ_TO_PHYSICAL(CSij6,Sij6,0,NY+1)

! We now have S_ij in Physical space

      RETURN
      END


      subroutine compute_rotation_chan
C This subroutine computes omgij (omega_ij) for the filtered velocity field
C The input velocity field should be in fourier space in the periodic
C directions.
C For use in the LES model in channel flow (2 periodic directions)
      include 'header'
      include 'header_les'

      integer I,J,K,ij

      

      DO J=1,NY
        DO K=0,TNKZ
          DO I=0,NXP-1
            CSij1(I,K,J)=CIKX(I)*CU1(I,K,J)
            CSij2(I,K,J)=(CU2(I,K,J+1)-CU2(I,K,J))/DYF(j)
            CSij3(I,K,J)=CIKZ(K)*CU3(I,K,J)
            Comgij5(I,K,J)=0.5d0*(CIKZ(K)*CU1(I,K,J)
     &                  -CIKX(I)*CU3(I,K,J))
C ^ this is omega_ik, this is anti-cyclic permutation
          END DO
        END DO
      END DO
      DO J=1,NY+1
        DO K=0,TNKZ
          DO I=0,NXP-1
            IF (DY(j).ne.0.d0) then
            Comgij4(I,K,J)=0.5d0*((CU1(I,K,J)-CU1(I,K,J-1))/DY(j)
     &                          -CIKX(I)*CU2(I,K,J) )
C ^ this is omega_ij, cyclic permutation
            Comgij6(I,K,J)=0.5d0*(CIKZ(K)*CU2(I,K,J)
     &                         -(CU3(I,K,J)-CU3(I,K,J-1))/DY(j) )
C ^ this is omega_jk, cyclic permutation
            END IF
          END DO
        END DO
      END DO


C       CALL FFT_XZ_TO_PHYSICAL(CSij1,Sij1,0,NY+1)
C       CALL FFT_XZ_TO_PHYSICAL(CSij2,Sij2,0,NY+1)
C       CALL FFT_XZ_TO_PHYSICAL(CSij3,Sij3,0,NY+1)

       CALL FFT_XZ_TO_PHYSICAL(Comgij4,omgij4,0,NY+1)
       CALL FFT_XZ_TO_PHYSICAL(Comgij5,omgij5,0,NY+1)
       CALL FFT_XZ_TO_PHYSICAL(Comgij6,omgij6,0,NY+1)

! We now have S_ij in Physical space

      RETURN
      END

      subroutine compute_strain_chan_invariant
C This subroutine computes S_ij for the filtered velocity field
C The input velocity field should be in fourier space in the periodic
C directions.
C For use in the LES model in channel flow (2 periodic directions)
      include 'header'
      include 'header_les'

      integer I,J,K,ij
      real*8 deltax,deltay,deltaz

!Set filter length (based on grid size) in x,z directions
!Based on constant Smag code
         deltax = (DX(1)*beta_sgs)
         deltaz = (DZ(1)*beta_sgs)
     
      DO J=1,NY
!Set filter length (based on grid size) in y direction
!Based on constant Smag code        
         deltay=(DYF(J)*2.d0)
        DO K=0,TNKZ
          DO I=0,NXP-1
            CSIij1(I,K,J)=CIKX(I)*CU1(I,K,J)
            CSIij2(I,K,J)=(CU2(I,K,J+1)-CU2(I,K,J))/DYF(j) 
            CSIij3(I,K,J)=CIKZ(K)*CU3(I,K,J)  
            CSIij5(I,K,J)=0.5d0*(CIKZ(K)*CU1(I,K,J)
     &                                           *(deltaz/deltax)
     &                  +CIKX(I)*CU3(I,K,J)
     &                                           *(deltax/deltaz) )
          END DO
        END DO
      END DO
      DO J=1,NY+1
         deltay=(DY(J)*2.d0)
        DO K=0,TNKZ
          DO I=0,NXP-1
            CSIij4(I,K,J)=0.5d0*((CU1(I,K,J)-CU1(I,K,J-1))/DY(j)
     &                                           *(deltay/deltax)
     &                          +CIKX(I)*CU2(I,K,J) 
     &                                           *(deltax/deltay) ) 
            CSIij6(I,K,J)=0.5d0*(CIKZ(K)*CU2(I,K,J)
     &                                           *(deltaz/deltay)
     &                         +(CU3(I,K,J)-CU3(I,K,J-1))/DY(j) 
     &                                           *(deltay/deltaz) )
          END DO
        END DO
      END DO  


       CALL FFT_XZ_TO_PHYSICAL(CSIij1,SIij1,0,NY+1)
       CALL FFT_XZ_TO_PHYSICAL(CSIij2,SIij2,0,NY+1)
       CALL FFT_XZ_TO_PHYSICAL(CSIij3,SIij3,0,NY+1)
       CALL FFT_XZ_TO_PHYSICAL(CSIij4,SIij4,0,NY+1)
       CALL FFT_XZ_TO_PHYSICAL(CSIij5,SIij5,0,NY+1)
       CALL FFT_XZ_TO_PHYSICAL(CSIij6,SIij6,0,NY+1)

! We now have S_ij in Physical space

      RETURN
      END



      subroutine compute_rotation_chan_invariant
C This subroutine computes omgij (omega_ij) for the filtered velocity field
C The input velocity field should be in fourier space in the periodic
C directions.
C For use in the LES model in channel flow (2 periodic directions)
      include 'header'
      include 'header_les'

      integer I,J,K,ij
      real*8 deltax,deltay,deltaz

!Set filter length (based on grid size) in x,z directions
!Based on constant Smag code
         deltax = (DX(1)*beta_sgs)
         deltaz = (DZ(1)*beta_sgs)

      DO J=1,NY
         deltay=(DYF(J)*2.d0)
        DO K=0,TNKZ
          DO I=0,NXP-1
C            CSij1(I,K,J)=CIKX(I)*CU1(I,K,J)
C            CSij2(I,K,J)=(CU2(I,K,J+1)-CU2(I,K,J))/DYF(j)
C            CSij3(I,K,J)=CIKZ(K)*CU3(I,K,J)
            Comgij5(I,K,J)=0.5d0*(CIKZ(K)*CU1(I,K,J)
     &                                           *(deltaz/deltax)
     &                  -CIKX(I)*CU3(I,K,J)
     &                                           *(deltax/deltaz) ) 
C ^ this is omega_ik, this is anti-cyclic permutation
          END DO
        END DO
      END DO
      DO J=1,NY+1
        deltay=(DY(J)*2.d0)
        if (DY(j).ne.0.d0) then
        DO K=0,TNKZ
          DO I=0,NXP-1
            Comgij4(I,K,J)=0.5d0*((CU1(I,K,J)-CU1(I,K,J-1))/DY(j)
     &                                           *(deltay/deltax)
     &                          -CIKX(I)*CU2(I,K,J) 
     &                                           *(deltax/deltay) )
C ^ this is omega_ij, cyclic permutation
            Comgij6(I,K,J)=0.5d0*(CIKZ(K)*CU2(I,K,J)
     &                                           *(deltaz/deltay)
     &                         -(CU3(I,K,J)-CU3(I,K,J-1))/DY(j) 
     &                                           *(deltay/deltaz) )
C ^ this is omega_jk, cyclic permutation
          END DO
        END DO
        end if
      END DO 


C       CALL FFT_XZ_TO_PHYSICAL(CSij1,Sij1,0,NY+1)
C       CALL FFT_XZ_TO_PHYSICAL(CSij2,Sij2,0,NY+1)
C       CALL FFT_XZ_TO_PHYSICAL(CSij3,Sij3,0,NY+1)
       CALL FFT_XZ_TO_PHYSICAL(Comgij4,omgij4,0,NY+1)
       CALL FFT_XZ_TO_PHYSICAL(Comgij5,omgij5,0,NY+1)
       CALL FFT_XZ_TO_PHYSICAL(Comgij6,omgij6,0,NY+1)

! We now have S_ij in Physical space

      RETURN
      END


 
      subroutine les_filter_chan(A,jstart,jend)
! This subroutine applies the les filter to the input field
! The indices to the start and end of the array in the y-direction
! are also inputted to make the routine cablable of filtering fields
! at either GYF or GY points.
! The array that is passed should be in physical space
      integer i,j,k,n
      integer jstart,jend,NX,NY,NZ,NXM,NZM,N_TH
      INCLUDE 'grid_def'

      real*8 A(0:NX+1,0:NZ+1,0:NY+1)
      real*8 B(0:NX-1,0:NZ-1,0:NY+1)

      integer im2(0:NX-1),im1(0:NX-1),ip1(0:NX+1),ip2(0:NX+2)
      integer km2(0:NZ-1),km1(0:NZ-1),kp1(0:NZ+1),kp2(0:NZ+2)

! These are the weights for the filtering operation used
      real*8 W0,W1,W2,Wm1,Wm2,Wm1_j,W0_j,W1_j


! The following is for the 3-point trapezoidal rule, alpha*beta=sqrt(6)
!      Wm2=0.d0
!      Wm1=1.d0/4.d0
!      W0=1.d0/2.d0
!      W1=1.d0/4.d0
!      W2=0.d0
      Wm1_j=1.d0/4.d0  
      W0_j=1.d0/2.d0
      W1_j=1.d0/4.d0
! The following is for the 5-point trapezoidal rule, alpha*beta=9
      Wm2=1.d0/8.d0
      Wm1=1.d0/4.d0
      W0=1.d0/4.d0
      W1=1.d0/4.d0
      W2=1.d0/8.d0

      NXM=NX-1
      NZM=NZ-1

!      do j=0,NY+1
!        do k=0,NZM
!          do i=0,NXM
!            B(i,k,j)=A(i,k,j)
!          end do
!        end do
!      end do

! Filter in the periodic directions using cshift
! Note, cshift if not used since it appears to be slower
! Apply filter in the x-direction
!      B=Wm2*CSHIFT(B,-2,1)+Wm1*CSHIFT(B,-1,1)+W0*B+W1*CSHIFT(B,1,1)
!     &       +W2*CSHIFT(B,2,1)

! Filter using more efficient F77 syntax:
! Set up array to loop around periodic directions
      do i=2,NXM
        im2(i)=i-2
      end do
      im2(1)=NXM
      im2(0)=NX-2
      do i=1,NXM
        im1(i)=i-1
      end do
      im1(0)=NXM
      do i=0,NX-2
        ip1(i)=i+1
      end do
      ip1(NXM)=0
      do i=0,NX-3
        ip2(i)=i+2    
      end do
      ip2(NX-2)=0
      ip2(NXM)=1

      do j=jstart,jend
        do k=0,NZM
          do i=0,NXM
            B(i,k,j)=Wm2*A(im2(i),k,j)+Wm1*A(im1(i),k,j)+W0*A(i,k,j)
     &         +W1*A(ip1(i),k,j)+W2*A(ip2(i),k,j)
          end do
        end do  
      end do
 
! Apply filter in the z-diretion
!      B=Wm2*CSHIFT(B,-2,2)+Wm1*CSHIFT(B,-1,2)+W0*B+W1*CSHIFT(B,1,2)
!     &       +W2*CSHIFT(B,2,2)
! Filter using more efficient F77 syntax:
! Set up array to loop around periodic directions
      do k=2,NZM
        km2(k)=k-2
      end do
      km2(1)=NZM
      km2(0)=NZ-2
      do k=1,NZM
        km1(k)=k-1
      end do
      km1(0)=NZM
      do k=0,NZ-2
        kp1(k)=k+1
      end do
      kp1(NZM)=0
      do k=0,NZ-3
        kp2(k)=k+2    
      end do
      kp2(NZ-2)=0
      kp2(NZM)=1

      do j=jstart,jend
        do k=0,NZM
          do i=0,NXM
            A(i,k,j)=Wm2*B(i,km2(k),j)+Wm1*B(i,km1(k),j)+W0*B(i,k,j)
     &         +W1*B(i,kp1(k),j)+W2*B(i,kp2(k),j)
          end do
        end do  
      end do

! Apply filter in the vertical direction at all physical cells
! (filter is not applied to ghost cells, but the values of the ghost cells
! is used for averaging)
!      B(:,:,jstart+1:jend-1)=W0*B(:,:,jstart:jend-2)
!     &                      +W1*B(:,:,jstart+1:jend-1)
!     &                      +W2*B(:,:,jstart+2:jend)
! Use more efficient F77 syntax:
!       do j=jstart+1,jend-1
!         do k=0,NZM
!           do i=0,NXM
!             B(i,k,j)=Wm1_j*B(i,k,j-1)+W0_j*B(i,k,j)+W1_j*B(i,k,j+1)  
!           end do
!         end do
!       end do

!      do j=jstart,jend
!        do k=0,NZM
!          do i=0,NXM
!            A(i,k,j)=B(i,k,j)
!          end do
!        end do
!      end do

      return
      end


      subroutine les_filter_chan_fourier(A,jstart,jend)
! This subroutine applies the les filter to the input field
! The filter is a spectral cutoff filter
! The indices to the start and end of the array in the y-direction
! are also inputted to make the routine cablable of filtering fields
! at either GYF or GY points.
! The array that is passed should be in physical space
      integer jstart,jend,NX,NZ,NY,NXM,NZM,i,j,k,N_TH

      INCLUDE 'grid_def'

      real*8 PI, LX, LZ
      integer NKX,NKZ,TNKZ

      real*8 KX(0:NX/3),KZ(0:2*(NZ/3))  

      real*8 A(0:NX+1,0:NZ+1,0:NY+1)

      real alpha

      real*8 B(0:NX+1,0:NZ+1,0:NY+1)

      complex*16 CB(0:NX/2,0:NZ+1,0:NY+1)

ccc      equivalence (B,CB) !by J.Liu, 18/Jan/2020

! Set the ratio of filter scales
      parameter (alpha=2.d0)

      NXM=NX-1
      NZM=NZ-1

      PI = 4. * ATAN(1.0)

      LX=PI
      LZ=2.d0*PI
   
! Get the wavenumber vectors:
        NKX=NX/3
        DO I=0,NKX
          KX(I)=I*(2.*PI)/LX
        END DO

        NKZ=NZ/3
        TNKZ=NKZ*2
        DO K=0,NKZ
          KZ(K)=K*(2.*PI)/LZ
        END DO
        DO K=1,NKZ
          KZ(TNKZ+1-K)=-K*(2.*PI)/LZ
        END DO

       do j=0,NY+1
         do k=0,NZM
           do i=0,NXM
             B(i,k,j)=A(i,k,j)
           end do
         enddo
       end do 


! Convert to fourier space
      call fft_xz_to_fourier(B,CB,jstart,jend)

! Perform the filtering
      do j=jstart,jend
        do k=0,TNKZ
          do i=0,NKX
            if (sqrt(KX(I)**2.d0+KZ(K)**2.d0)
     &       .gt.sqrt(KX(NKX)**2.d0+KZ(NKZ)**2.d0)/alpha) then
                CB(i,k,j)=0.d0
            end if
          end do
        end do
      end do

! Now, convert back to physical space
      call fft_xz_to_physical(CB,B,jstart,jend)

       do j=jstart,jend 
         do k=0,NZM
           do i=0,NXM
             A(i,k,j)=B(i,k,j)
           end do
         end do
       end do 

      return
      end

      SUBROUTINE APPLY_BC_LES
      include 'header'
      integer i,j,k

! If we are using Neuman boundary conditions, over-write the values of the
! velocity at the ghost cells so that the LES model doesn't use the large
! velocity gradient
      IF (U_BC_YMAX.eq.1) THEN
        IF ((RANKY.eq.NPROCY-1).or.(.NOT.USE_MPI)) THEN
! We are on process at the upper wall
          DO K=0,TNKZ
           DO I=0,NXP-1
             CU1(I,K,NY)=CU1(I,K,NY-1)
             CU1(I,K,NY+1)=CU1(I,K,NY)
           END DO
          END DO
        END IF
      END IF

      IF (W_BC_YMAX.eq.1) THEN
        IF ((RANKY.eq.NPROCY-1).or.(.NOT.USE_MPI)) THEN
! We are on process at the upper wall
          DO K=0,TNKZ
           DO I=0,NXP-1
             CU3(I,K,NY)=CU3(I,K,NY-1)
             CU3(I,K,NY+1)=CU3(I,K,NY)
           END DO
          END DO
        END IF
      END IF

      IF (U_BC_YMIN.eq.1) THEN
        IF ((RANKY.eq.0).or.(.NOT.USE_MPI)) THEN
! We are on a process at the bottom wall          
         DO K=0,TNKZ
           DO I=0,NXP-1
             CU1(I,K,1)=CU1(I,K,2)
             CU1(I,K,0)=CU1(I,K,1)
           END DO
         END DO
       END IF
      END IF

      IF (W_BC_YMIN.eq.1) THEN
        IF ((RANKY.eq.0).or.(.NOT.USE_MPI)) THEN
! We are on a process at the bottom wall          
         DO K=0,TNKZ
           DO I=0,NXP-1
             CU3(I,K,1)=CU3(I,K,2)
             CU3(I,K,0)=CU3(I,K,1)
           END DO
         END DO
       END IF
      END IF

      RETURN
      END






