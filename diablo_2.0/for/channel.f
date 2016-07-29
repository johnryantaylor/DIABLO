C******************************************************************************|
C channel.f, the channel-flow solvers for diablo.                  VERSION 1.0
C This solver was written by John R. Taylor.
C******************************************************************************|
!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE INIT_CHAN
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
C This subroutine initializes various variables for the serial channel
C flow subroutine. This is not called if MPI is turned on
C Initialize any constants here

      INCLUDE 'header'
      INTEGER I,J,K,N

      PI=4.D0*ATAN(1.D0)

! Defined starting and ending indeces in the wall-bounded direction
! These will be used for loops in the GYF variables
! The start and end of these loops depends on the bounary conditions
! For example, Direchlet BCs do not need to update the wall values

        IF (U_BC_YMIN.EQ.0) THEN
          JSTART=2
        ELSE IF (U_BC_YMIN.EQ.1) THEN
          JSTART=1
        ELSE
          JSTART=2
        END IF
! Now, set the indexing for the scalar equations
        DO N=1,N_TH
          IF (TH_BC_YMIN(N).EQ.0) THEN
            JSTART_TH(N)=2
          ELSE IF (TH_BC_YMIN(N).EQ.1) THEN
            JSTART_TH(N)=1
          ELSE
            JSTART_TH(N)=2
          END IF
        END DO
        IF (U_BC_YMAX.EQ.0) THEN
          JEND=NY-1
        ELSE IF (U_BC_YMAX.EQ.1) THEN
          JEND=NY
        ELSE
          JEND=NY-1
        END IF

! Set the upper and lower limits of timestepping of the scalar equations
        DO N=1,N_TH
        IF (TH_BC_YMAX(N).EQ.0) THEN
          JEND_TH(N)=NY-1
        ELSE IF (TH_BC_YMAX(N).EQ.1) THEN
          JEND_TH(N)=NY
        ELSE
          JEND_TH(N)=NY-1
        END IF
        END DO

      RETURN
      END

C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE RK_CHAN_1
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
C Main time-stepping algorithm for the channel-flow case.
C This algorithm uses Crank-Nicolson for viscous/diffusive terms and
C and 3rd order Runge-Kutta for the rest of the terms
C INPUTS  (in Fourier space):  CUi, CP, and (if k>1) CFi at (k-1)  (for i=1,2,3)
C OUTPUTS (in Fourier space):  CUi, CP, and (if k<3) CFi at (k)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      INCLUDE 'header'

      INTEGER I,J,K,N,ISTART      
      REAL*8 TEMP1, TEMP2, TEMP3, TEMP4, TEMP5, UBULK

C Define the constants that are used in the time-stepping
C For reference, see Numerical Renaissance
      TEMP1=NU * H_BAR(RK_STEP) / 2.0
      TEMP2=H_BAR(RK_STEP) / 2.0
      TEMP3=ZETA_BAR(RK_STEP) * H_BAR(RK_STEP)
      TEMP4=H_BAR(RK_STEP)
      TEMP5=BETA_BAR(RK_STEP) * H_BAR(RK_STEP)

C First, we will compute the explicit RHS terms and store in Ri
C Note, Momentum equation and hence the RHS is evaluated at the
C corresponding velocity points.

C Store the old velocity in the RHS vector
      DO J=JSTART,JEND
        DO K=0,TNKZ
          DO I=0,NXP-1
            CR1(I,K,J)=CU1(I,K,J)
            CR3(I,K,J)=CU3(I,K,J)
          END DO
        END DO
      END DO
      DO J=2,NY 
        DO K=0,TNKZ
          DO I=0,NXP-1
            CR2(I,K,J)=CU2(I,K,J)
          END DO
        END DO
      END DO

C Add the R-K term from the rk-1 step 
      IF (RK_STEP .GT. 1) THEN
        DO J=JSTART,JEND
          DO K=0,TNKZ
            DO I=0,NXP-1
              CR1(I,K,J)=CR1(I,K,J)+TEMP3*CF1(I,K,J)
              CR3(I,K,J)=CR3(I,K,J)+TEMP3*CF3(I,K,J)
            END DO
          END DO
        END DO
        DO J=2,NY
          DO K=0,TNKZ
            DO I=0,NXP-1
              CR2(I,K,J)=CR2(I,K,J)+TEMP3*CF2(I,K,J)
            END DO
          END DO
        END DO
      END IF
          
C Take the y-derivative of the pressure at GY points in Fourier space
      DO J=2,NY
        DO K=0,TNKZ
          DO I=0,NXP-1
            CS1(I,K,J)=(CP(I,K,J) - CP(I,K,J-1)) / DY(J)
          END DO
        END DO
      END DO

C Add the pressure gradient to the RHS as explicit Euler
      DO J=JSTART,JEND
        DO K=0,TNKZ
          DO I=0,NXP-1
            CR1(I,K,J)=CR1(I,K,J)-TEMP4*(CIKX(I)*CP(I,K,J))
            CR3(I,K,J)=CR3(I,K,J)-TEMP4*(CIKZ(K)*CP(I,K,J))
          END DO
        END DO
      END DO

      DO J=2,NY
        DO K=0,TNKZ
          DO I=0,NXP-1
            CR2(I,K,J)=CR2(I,K,J)-TEMP4*CS1(I,K,J)
          END DO
        END DO
      END DO

C Here, add the constant, forcing pressure gradient
C There are two built-in ways of doing this
C F_TYPE=1 -> Constant pressure gradient in the x-direction
C F_TYPE=2 -> Oscillatory pressure gradient in the x-direction
C OTHER VALUES -> No forcing added
      IF (F_TYPE.EQ.1) THEN 
C Add forcing for a constant pressure gradient
        DO J=JSTART,JEND
          IF (RANKZ.eq.0) CR1(0,0,J)=CR1(0,0,J)-TEMP4*PX0
        END DO
      ELSE IF (F_TYPE.EQ.2) THEN
C Oscillatory pressure gradient
        DO J=JSTART,JEND
           IF (RANKZ.eq.0) CR1(0,0,J)=CR1(0,0,J)-
     &          TEMP4*(PX0+AMP_OMEGA0*cos(OMEGA0*TIME))
        END DO
C End if forcing type
      END IF

C Now compute the term R-K term Ai
C Compile terms of Ai in CFi which will be saved for next time step
C First, store the horizontal viscous terms in CFi
C Note that BETA is the order of the Laplacian operator
C e.g. BETA=1 for second order viscosity/diffusivity
C or   BETA=2, etc. for hyper viscosity/diffusivity
      DO J=JSTART,JEND
        DO K=0,TNKZ
          DO I=0,NXP-1
            CF1(I,K,J)=-NU * (KX2(I)+KZ2(K))**BETA * CU1(I,K,J) 
            CF3(I,K,J)=-NU * (KX2(I)+KZ2(K))**BETA * CU3(I,K,J) 
          END DO
        END DO
      END DO
      DO J=2,NY
        DO K=0,TNKZ
          DO I=0,NXP-1
            CF2(I,K,J)=-NU * (KX2(I)+KZ2(K))**BETA * CU2(I,K,J) 
          END DO 
        END DO
      END DO

! Add the terms owing to the system rotation (Coriolis terms) 
! Assume that the flow is on an f-plane
      DO K=0,TNKZ
        DO I=0,NXP-1
          DO J=JSTART,JEND
            CF1(I,K,J)=CF1(I,K,J)+I_RO*CU3(I,K,J)
            CF3(I,K,J)=CF3(I,K,J)-I_RO*CU1(I,K,J)
          END DO
        END DO
      END DO

! Do for each scalar
      DO N=1,N_TH

! If a scalar contributes to the denisty, RI is not equal to zero and
! add the buoyancy term as explicit R-K.  Don't add the 0,0 mode in the 
! y-direction, which corresponds to the plane-average.
! The plane averaged density balances the hydrostatic pressure 
      DO J=2,NY
        DO K=1,TNKZ
          DO I=0,NXP-1
! Use second order interpolation
             CF2(I,K,J)=CF2(I,K,J)+RI(N)*GRAV_Y*
     &      (CTH(I,K,J,N)*DYF(J-1)+CTH(I,K,J-1,N)*DYF(J))/(2.d0*DY(J))
          END DO
        END DO
        K=0
        IF (RANKZ.eq.0) THEN
           ISTART=1
        ELSE
           ISTART=0
        END IF
        DO I=ISTART,NXP-1
             CF2(I,K,J)=CF2(I,K,J)+RI(N)*GRAV_Y*
     &      (CTH(I,K,J,N)*DYF(J-1)+CTH(I,K,J-1,N)*DYF(J))/(2.d0*DY(J))
        END DO
      END DO

      DO J=JSTART,JEND
        DO K=0,TNKZ
          DO I=0,NXP-1
            CF1(I,K,J)=CF1(I,K,J)+RI(N)*GRAV_X*CTH(I,K,J,N)
          END DO
        END DO
      END DO

      DO J=JSTART,JEND
        DO K=0,TNKZ
          DO I=0,NXP-1
            CF3(I,K,J)=CF3(I,K,J)+RI(N)*GRAV_Z*CTH(I,K,J,N)
          END DO
        END DO
      END DO

! Now, compute the RHS vector for the scalar equations
! Since TH is defined at horizontal velocity points, the
! scalar update equation will be very similar to the horizontal
! velocity update.

! We will store the RHS scalar terms in CRTH, RTH
! The k-1 term for the R-K stepping is saved in FTH, CFTH

! First, build the RHS vector, use CRTH
      DO J=JSTART_TH(N),JEND_TH(N)
        DO K=0,TNKZ
          DO I=0,NXP-1
          CRTH(I,K,J,N)=CTH(I,K,J,N)
         ENDDO
       END DO
      END DO
! Add term from k-2 step to free up CFTH variable
      IF (RK_STEP .GT. 1) THEN
        DO J=JSTART_TH(N),JEND_TH(N)
          DO K=0,TNKZ
            DO I=0,NXP-1
              CRTH(I,K,J,N)=CRTH(I,K,J,N)+TEMP3*CFTH(I,K,J,N)
            END DO
          END DO
        END DO
       END IF

! Now compute the explicit R-K term Ai
! Compile terms of Ai in CFi which will be saved for next time step
      DO J=JSTART_TH(N),JEND_TH(N)
        DO K=0,TNKZ
          DO I=0,NXP-1
            CFTH(I,K,J,N)=-(NU/PR(N))*(KX2(I)+KZ2(K))**BETA*CTH(I,K,J,N)
          END DO
        END DO
      END DO

! Add advection acting on the background scalar gradient if present
      DO J=JSTART_TH(N),JEND_TH(N)
        DO K=0,TNKZ
          DO I=0,NXP-1
            CFTH(I,K,J,N)=CFTH(I,K,J,N)-CU1(I,K,J)*DRHODX(N)
     &                                 -CU3(I,K,J)*DRHODZ(N)
          END DO
        END DO
      END DO

C End do number of passive scalars (N_TH)
      END DO

C Optionally, add user forcing to the right hand side
C Here, we have U1, U2, U3, and TH in Fourier space
      CALL USER_RHS_CHAN_FOURIER

C If we are considering an LES, then add the subgrid scale stress:
C Here, velocity and CFi should be in Fourier space
C The subgrid scale stress is added to CFi:   CFi=CFi - d/dx_i tau_ij

      IF (LES.AND.((.NOT.CREATE_NEW_FLOW).OR.(TIME_STEP.GT.100))) THEN
C If we have created new flow with random perturbations, wait for a
C spinup before applying the subgrid model for stability purposes
C In the process, Ui is converted to physical space
          call les_chan

C APPLY constant SGS Prandlt number
         DO N=1,N_TH
         do j=1,NY+1
           do k=0,NZP-1
             do i=0,NXM
               KAPPA_T(I,K,J,N)=1.d0*NU_T(I,K,J)
             end do
           end do
         end do
         end do
        DO N=1,N_TH
          CS1(:,:,:)=CTH(:,:,:,N)
          CALL FFT_XZ_TO_PHYSICAL(CS1,S1,0,NY+1)
          TH(:,:,:,N)=S1(:,:,:)
        END DO

      ELSE 
C If the subgrid model hasn't been called, then it is necessary to 
C convert to physical space.
        CALL FFT_XZ_TO_PHYSICAL(CU1,U1,0,NY+1)
        CALL FFT_XZ_TO_PHYSICAL(CU2,U2,0,NY+1)
        CALL FFT_XZ_TO_PHYSICAL(CU3,U3,0,NY+1)
! Transform THETA to physical space for computation of nonlinear terms
! Here pass the first location in memory of the array for scalar n

        DO N=1,N_TH
          CS1(:,:,:)=CTH(:,:,:,N)
          CALL FFT_XZ_TO_PHYSICAL(CS1,S1,0,NY+1)
          TH(:,:,:,N)=S1(:,:,:)
        END DO      

      END IF


C Compute the nonlinear products in physical space, then transform
C back to Fourier space to compute the derivative.
C Here, we compute the horizontal derivatives of the nonlinear terms
C which will be treated with RKW3.  
C Do terms one at a time to save on memory
C U1*U3
      DO J=JSTART,JEND
        DO K=0,NZP-1
          DO I=0,NXM
            S1(I,K,J)=U3(I,K,J)*U1(I,K,J)
          END DO
        END DO
      END DO
      
      CALL FFT_XZ_TO_FOURIER(S1,CS1,0,NY+1)
      
      DO J=JSTART,JEND
        DO K=0,TNKZ
          DO I=0,NXP-1
            CF1(I,K,J)=CF1(I,K,J) - CIKZ(K) * CS1(I,K,J) 
            CF3(I,K,J)=CF3(I,K,J) - CIKX(I) * CS1(I,K,J) 
          END DO
        END DO
      END DO

C U1*U1
      DO J=JSTART,JEND
        DO K=0,NZP-1
          DO I=0,NXM
            S1(I,K,J)=U1(I,K,J)*U1(I,K,J)
          END DO
        END DO
      END DO
      
      CALL FFT_XZ_TO_FOURIER(S1,CS1,0,NY+1)
      
      DO J=JSTART,JEND
        DO K=0,TNKZ
          DO I=0,NXP-1
            CF1(I,K,J)=CF1(I,K,J) - CIKX(I) * CS1(I,K,J) 
          END DO
        END DO
      END DO

C U3*U3
      DO J=JSTART,JEND
        DO K=0,NZP-1
          DO I=0,NXM
            S1(I,K,J)=U3(I,K,J)*U3(I,K,J)
          END DO
        END DO
      END DO
      
      CALL FFT_XZ_TO_FOURIER(S1,CS1,0,NY+1)
      
      DO J=JSTART,JEND
        DO K=0,TNKZ
          DO I=0,NXP-1
            CF3(I,K,J)=CF3(I,K,J) - CIKZ(K) * CS1(I,K,J) 
          END DO
        END DO
      END DO


C U1*U2
      DO J=2,NY
        DO K=0,NZP-1
          DO I=0,NXM
            S1(I,K,J)=((DYF(J)*U1(I,K,J)
     &                +DYF(J-1)*U1(I,K,J-1))/(2.*DY(J))) 
     &                *U2(I,K,J)
          END DO
        END DO
      END DO
      
      CALL FFT_XZ_TO_FOURIER(S1,CS1,0,NY+1)
      
      DO J=2,NY
        DO K=0,TNKZ
          DO I=0,NXP-1
            CF2(I,K,J)=CF2(I,K,J) - CIKX(I) * CS1(I,K,J) 
          END DO
        END DO
      END DO

C U3*U2
      DO J=2,NY
        DO K=0,NZP-1
          DO I=0,NXM
            S1(I,K,J)=((DYF(J)*U3(I,K,J)
     &                +DYF(J-1)*U3(I,K,J-1))/(2.*DY(J))) 
     &                *U2(I,K,J)
          END DO
        END DO
      END DO
      
      CALL FFT_XZ_TO_FOURIER(S1,CS1,0,NY+1)
      
      DO J=2,NY
        DO K=0,TNKZ
          DO I=0,NXP-1
            CF2(I,K,J)=CF2(I,K,J) - CIKZ(K) * CS1(I,K,J)
          END DO
        END DO
      END DO

! Add the vertical derivative term
      DO J=JSTART,JEND
        DO K=0,NZP-1
          DO I=0,NXM
            S1(I,K,J)=
     &     (U1(I,K,J+1)*U2(I,K,J+1) + U1(I,K,J)*U2(I,K,J+1)
     &     - U1(I,K,J)*U2(I,K,J) - U1(I,K,J-1)*U2(I,K,J))/(2.d0*DYF(J))
          END DO
        END DO
      END DO

      CALL FFT_XZ_TO_FOURIER(S1,CS1,0,NY+1)

      DO J=JSTART,JEND
        DO K=0,TNKZ
          DO I=0,NXP-1
            CF1(I,K,J)=CF1(I,K,J) - CS1(I,K,J) 
         END DO
        END DO
      END DO

! Add the vertical derivative term explicitly
      DO J=JSTART,JEND
        DO K=0,NZP-1
          DO I=0,NXM
            S1(I,K,J)=
     &     (U3(I,K,J+1)*U2(I,K,J+1) + U3(I,K,J)*U2(I,K,J+1)
     &     - U3(I,K,J)*U2(I,K,J) - U3(I,K,J-1)*U2(I,K,J))/(2.d0*DYF(J))
          END DO
        END DO
      END DO

      CALL FFT_XZ_TO_FOURIER(S1,CS1,0,NY+1)

      DO J=JSTART,JEND
        DO K=0,TNKZ
          DO I=0,NXP-1
            CF3(I,K,J)=CF3(I,K,J) - CS1(I,K,J) 
          END DO
        END DO
      END DO

! Add the vertical derivative term explicitly
      DO J=2,NY
        DO K=0,NZP-1
          DO I=0,NXM
            S1(I,K,J)=
     &    (0.25d0*(U2(I,K,J)+U2(I,K,J+1))**2.d0
     &    -0.25d0*(U2(I,K,J)+U2(I,K,J-1))**2.d0)/DY(J)
          END DO
        END DO
      END DO

      CALL FFT_XZ_TO_FOURIER(S1,CS1,0,NY+1)

      DO J=2,NY
        DO K=0,TNKZ
          DO I=0,NXP-1
            CF2(I,K,J)=CF2(I,K,J) - CS1(I,K,J)
          END DO
        END DO
      END DO


C -- At this point, we are done computing the nonlinear terms --

C Optionally, add user forcing to the right hand side
C Here, we have U1, U2, U3, and TH in physical space
      CALL USER_RHS_CHAN_PHYSICAL

C Finally, Add CFi to CRi
      DO J=JSTART,JEND
        DO K=0,TNKZ
          DO I=0,NXP-1
            CR1(I,K,J)=CR1(I,K,J) + TEMP5 * CF1(I,K,J)
            CR3(I,K,J)=CR3(I,K,J) + TEMP5 * CF3(I,K,J)
          END DO
        END DO
      END DO
      DO J=2,NY
        DO K=0,TNKZ
          DO I=0,NXP-1
            CR2(I,K,J)=CR2(I,K,J) + TEMP5 * CF2(I,K,J)
          END DO
        END DO
      END DO

C Convert RHS terms to physical space
      CALL FFT_XZ_TO_PHYSICAL(CR1,R1,0,NY+1)                 
      CALL FFT_XZ_TO_PHYSICAL(CR2,R2,2,NY)                 
      CALL FFT_XZ_TO_PHYSICAL(CR3,R3,0,NY+1)                 

C Compute the vertical viscous term in physical space and add to RHS
C This is the explicit part of the Crank-Nicolson term
C NU_V_SCALE is a coefficient defined as the ratio of the vertical
C to horizontal viscosity and can be used to add anisotropic viscosity
      DO J=JSTART,JEND
        DO K=0,NZP-1
          DO I=0,NXM
            R1(I,K,J)=R1(I,K,J)+TEMP1*NU_V_SCALE*
     &        (  ((U1(I,K,J+1) - U1(I,K,J)) / DY(J+1)  
     &           -(U1(I,K,J)   - U1(I,K,J-1)) / DY(J)) /DYF(J)  )
            R3(I,K,J)=R3(I,K,J)+TEMP1*NU_V_SCALE*
     &        (  ((U3(I,K,J+1) - U3(I,K,J)) / DY(J+1) 
     &           -(U3(I,K,J)   - U3(I,K,J-1)) / DY(J)) /DYF(J)  )
          END DO
        END DO
      END DO
      DO J=2,NY 
        DO K=0,NZP-1
          DO I=0,NXM
            R2(I,K,J)=R2(I,K,J)+TEMP1*NU_V_SCALE*
     &        (  ((U2(I,K,J+1) - U2(I,K,J))  / DYF(J) 
     &           -(U2(I,K,J)   - U2(I,K,J-1))/ DYF(J-1))/DY(J)  )
          END DO
        END DO
      END DO

C If we are using a subgrid model, add the eddy viscosity term
C This is an added viscosity that will be treated just like the 
C molecular viscosity with Crank-Nicolson for the vertical derivatives
      IF (LES) then
C Note, NU_T is defined at GY points
      DO J=JSTART,JEND
        DO K=0,NZP-1
          DO I=0,NXM
            R1(I,K,J)=R1(I,K,J)+TEMP2*
     &        (  (NU_T(I,K,J+1) * (U1(I,K,J+1) - U1(I,K,J)) / DY(J+1)  
     &         -  NU_T(I,K,J) * (U1(I,K,J)   - U1(I,K,J-1)) / DY(J))
     &               /DYF(J)  )
            R3(I,K,J)=R3(I,K,J)+TEMP2*
     &        (  (NU_T(I,K,J+1) * (U3(I,K,J+1) - U3(I,K,J)) / DY(J+1) 
     &        - NU_T(I,K,J) * (U3(I,K,J)   - U3(I,K,J-1)) / DY(J)) 
     &              /DYF(J)  )
          END DO
        END DO
      END DO
! Here, interpolate NU_T to GYF points
      DO J=2,NY 
        DO K=0,NZP-1
          DO I=0,NXM
            R2(I,K,J)=R2(I,K,J)+TEMP2*
     &     ((0.5d0*(NU_T(I,K,J)+NU_T(I,K,J+1))*(U2(I,K,J+1)-U2(I,K,J))
     &                                              / DYF(J) 
     &    -0.5d0*(NU_T(I,K,J)+NU_T(I,K,J-1))*(U2(I,K,J)-U2(I,K,J-1))
     &                                          / DYF(J-1))   /DY(J)  )
          END DO
        END DO
      END DO
      END IF

C -- Here, we are done with computation of Velocity RHS, explicit terms --

C Now, build the explicit RHS terms for the passive scalar(s)

      DO N=1,N_TH
! Do for each scalar:

! Compute the nonlinear terms that are present in the explicit term A
! U1*TH
      DO J=JSTART_TH(N),JEND_TH(N)
        DO K=0,NZP-1
          DO I=0,NXM
            S1(I,K,J)=TH(I,K,J,N)*U1(I,K,J)
          END DO
        END DO
      END DO
      CALL FFT_XZ_TO_FOURIER(S1,CS1,0,NY+1)
      DO J=JSTART_TH(N),JEND_TH(N)
        DO K=0,TNKZ
          DO I=0,NXP-1
            CFTH(I,K,J,N)=CFTH(I,K,J,N) - CIKX(I) * CS1(I,K,J)
          END DO
        END DO
      END DO
! U3*TH 
      DO J=JSTART_TH(N),JEND_TH(N)
        DO K=0,NZP-1
          DO I=0,NXM
            S1(I,K,J)=TH(I,K,J,N)*U3(I,K,J)
          END DO
        END DO
      END DO
      CALL FFT_XZ_TO_FOURIER(S1,CS1,0,NY+1)
      DO J=JSTART_TH(N),JEND_TH(N)
        DO K=0,TNKZ
          DO I=0,NXP-1
            CFTH(I,K,J,N)=CFTH(I,K,J,N) - CIKZ(K) * CS1(I,K,J)
          END DO
        END DO
      END DO

! We are done with the horizontal derivatives of the nonlinear terms
! Add the vertical derivative term explicitly
      DO J=JSTART_TH(N),JEND_TH(N)
        DO K=0,NZP-1
          DO I=0,NXM
            S1(I,K,J)=
     &     (TH(I,K,J+1,N)*U2(I,K,J+1) + TH(I,K,J,N)*U2(I,K,J+1)
     &    -TH(I,K,J,N)*U2(I,K,J)-TH(I,K,J-1,N)*U2(I,K,J))/(2.d0*DYF(J))
          END DO
        END DO
      END DO
      CALL FFT_XZ_TO_FOURIER(S1,CS1,0,NY+1)
      DO J=JSTART_TH(N),JEND_TH(N)
        DO K=0,TNKZ
          DO I=0,NXP-1
            CFTH(I,K,J,N)=CFTH(I,K,J,N) - CS1(I,K,J)
          END DO
        END DO
      END DO

! Add CFTH to the RHS vector CRTH
      DO J=JSTART_TH(N),JEND_TH(N)
        DO K=0,TNKZ
          DO I=0,NXP-1
            CRTH(I,K,J,N)=CRTH(I,K,J,N) + TEMP5 * CFTH(I,K,J,N)
          END DO
        END DO
      END DO
! Done with computation of RHS, explicit terms for the THETA equation
! Transform back to physical space

      CS1(:,:,:)=CRTH(:,:,:,N)
      CALL FFT_XZ_TO_PHYSICAL(CS1,S1,0,NY+1)
      RTH(:,:,:,N)=S1(:,:,:)

! Compute the Explicit part of the Crank-Nicolson terms for the TH equation
! First, the vertical derivative viscous term
      DO J=JSTART_TH(N),JEND_TH(N)
        DO K=0,NZP-1
          DO I=0,NXM
            RTH(I,K,J,N)=RTH(I,K,J,N)+(TEMP1/PR(N))*NU_V_SCALE*(
     &            ((TH(I,K,J+1,N) - TH(I,K,J,N)) / DY(J+1)
     &            -(TH(I,K,J,N) - TH(I,K,J-1,N)) / DY(J)) / DYF(J) )
          END DO
        END DO
      END DO
! If we are using a subgrid model (LES) then add the eddy diffusivity here
! Note, KAPPA_T is defined at GY points
      IF (LES) THEN
      DO J=JSTART_TH(N),JEND_TH(N)
        DO K=0,NZP-1
          DO I=0,NXM
            RTH(I,K,J,N)=RTH(I,K,J,N)+TEMP2*(
     &     (KAPPA_T(I,K,J+1,N)*(TH(I,K,J+1,N)-TH(I,K,J,N))/DY(J+1)
     &     -KAPPA_T(I,K,J,N)*(TH(I,K,J,N)-TH(I,K,J-1,N))/DY(J))/DYF(J))
          END DO
        END DO
      END DO  
      END IF

C -- Now, timestep the passive scalar equation --
C      We solve the the passive scalar before the velocity so that
C      it is advected with the velocity from the previous R-K step
C      which we have already made divergence free 
 
! Solve the implicit equation for THETA
! Note that the system size is NY+1, but only 1..NY are used

! Initialize the matrix used to store implicit coefficients
      DO J=0,NY+1
        DO I=0,NXM
          MATL(I,J)=0.
          MATD(I,J)=1.
          MATU(I,J)=0.
          VEC(I,J)=0.
        END DO
      END DO 
    
! Build implicit matrix
! Use quasi-second order interpolation for TH on GY points
      DO K=0,NZP-1
        DO J=JSTART_TH(N),JEND_TH(N)
          DO I=0,NXM
            MATL(I,J) = -(TEMP1/PR(N)*NU_V_SCALE) / (DY(J)*DYF(J))
            MATD(I,J) = 1. + (TEMP1/PR(N)*NU_V_SCALE) / (DY(J+1)*DYF(J))
     &           +(TEMP1/PR(N)*NU_V_SCALE) / (DY(J)*DYF(J))
            MATU(I,J)=-(TEMP1/PR(N)*NU_V_SCALE) / (DY(J+1)*DYF(J))
! Define RHS vector
            VEC(I,J)=RTH(I,K,J,N)
          END DO
        END DO
! IF using a subgrid model (LES) then add the eddy diffusivity part implicitly
        IF (LES) THEN
        DO J=JSTART_TH(N),JEND_TH(N)
          DO I=0,NXM   
            MATL(I,J) = MATL(I,J) - TEMP2 * KAPPA_T(I,K,J,N) 
     &                                 / (DY(J)*DYF(J))
            MATD(I,J) = MATD(I,J)+ TEMP2 * KAPPA_T(I,K,J+1,N)
     &                                 / (DY(J+1)*DYF(J))
     &                           + TEMP2 * KAPPA_T(I,K,J,N)
     &                                 / (DY(J)*DYF(J))
            MATU(I,J) = MATU(I,J)- TEMP2 * KAPPA_T(I,K,J+1,N)
     &                                / (DY(J+1)*DYF(J))
          END DO
        END DO
        END IF

! If we are using MPI, then solve the implicit system in separate forward
! and backward sweeps for efficiency
          IF (USE_MPI) THEN

             CALL APPLY_BC_TH_MPI(MATL,MATD,MATU,VEC,N)
! If we are using MPI, split the implicit solve into foward and
! backward sweeps for efficiency
             CALL THOMAS_FORWARD_REAL_MPI(MATL,MATD,MATU,VEC,NY,NX)
             CALL THOMAS_BACKWARD_REAL_MPI(MATL,MATD,MATU,VEC,NY,NX)
          ELSE
! Else we are running in serial mode
             CALL APPLY_BC_TH_LOWER(MATL,MATD,MATU,VEC,N)
             CALL APPLY_BC_TH_UPPER(MATL,MATD,MATU,VEC,N)
             CALL THOMAS_REAL(MATL,MATD,MATU,VEC,NY+1,NXM)
! Apply the boundary conditions to our linear system
          END IF

        DO J=JSTART_TH(N),JEND_TH(N)
          DO I=0,NXM
            TH(I,K,J,N)=VEC(I,J)
          END DO
        END DO

! END do k
      END DO 

! End do number of passive scalars
        END DO
        
C Initialize the matrix to zeros to be used for implicit solves
C Note that the system size is NY+1, but only 1..NY are used

! Initialize the matrix used to store implicit coefficients
      DO J=0,NY+1
        DO I=0,NXM
          MATL(I,J)=0.
          MATD(I,J)=1.
          MATU(I,J)=0.
          VEC(I,J)=0.
        END DO
      END DO 

C Build implicit matrix for U2
      DO K=0,NZP-1
        DO J=2,NY
          DO I=0,NXM
            MATL(I,J)= -TEMP1*NU_V_SCALE/(DYF(J-1)*DY(J))
            MATD(I,J)=1.+TEMP1*NU_V_SCALE/(DYF(J)*DY(J)) 
     &                 + TEMP1*NU_V_SCALE/(DYF(J-1)*DY(J)) 
            MATU(I,J)= -TEMP1*NU_V_SCALE/(DYF(J)*DY(J))
            VEC(I,J)=R2(I,K,J)
          END DO 
        END DO
        IF (LES) THEN
! IF using a subgrid model (LES) then add the eddy viscosity part implicitly
        DO J=2,NY
          DO I=0,NXM
            MATL(I,J) = MATL(I,J) 
     &      - TEMP2 * 0.5d0*(NU_T(I,K,J)+NU_T(I,K,J-1))/(DYF(J-1)*DY(J))
            MATD(I,J) = MATD(I,J) 
     &      + TEMP2 * 0.5d0*(NU_T(I,K,J)+NU_T(I,K,J+1))/(DYF(J)*DY(J))
     &      + TEMP2 * 0.5d0*(NU_T(I,K,J)+NU_T(I,K,J-1))/(DYF(J-1)*DY(J))
            MATU(I,J) = MATU(I,J) 
     &      - TEMP2 * 0.5d0*(NU_T(I,K,J)+NU_T(I,K,J+1))/(DYF(J)*DY(J))
          END DO
        END DO
        END IF

        IF (USE_MPI) THEN

! First, apply the boundary conditions
          CALL APPLY_BC_U2_MPI(MATL,MATD,MATU,VEC)
! If we are using MPI, split the implicit solve into forward and
! backward sweeps for efficiency
          CALL THOMAS_FORWARD_REAL_MPI(MATL,MATD,MATU,VEC,NY,NX)
          CALL THOMAS_BACKWARD_REAL_MPI(MATL,MATD,MATU,VEC,NY,NX)
        ELSE
C Else, we are running in serial mode
C Set the boundary conditions for U2
          CALL APPLY_BC_2_LOWER(MATL,MATD,MATU,VEC)
          CALL APPLY_BC_2_UPPER(MATL,MATD,MATU,VEC)

C Now, solve the tridiagonal system for U2(i,:,k)
          CALL THOMAS_REAL(MATL,MATD,MATU,VEC,NY+1,NXM)
        END IF

        DO J=1,NY+1
          DO I=0,NXM
            U2(I,K,J)=VEC(I,J)
          END DO
        END DO
! End do k
      END DO 

C Solve for U1
C Note, here the matrix will be indexed from 1...NY+1 corresponding to U1(0:NY)

! Initialize the matrix used to store implicit coefficients
      DO J=0,NY+1
        DO I=0,NXM
          MATL(I,J)=0.
          MATD(I,J)=1.
          MATU(I,J)=0.
          VEC(I,J)=0.
        END DO
      END DO 

C Build the implicit system of equations for U1 
      DO K=0,NZP-1
        DO J=JSTART,JEND
          DO I=0,NXM
            MATL(I,J)=-TEMP1*NU_V_SCALE/(DY(J)*DYF(J))
            MATD(I,J)=1.-TEMP1*NU_V_SCALE*(-1./(DY(J+1)*DYF(J))
     &         -1./(DY(J)*DYF(J))) 
            MATU(I,J)=-TEMP1*NU_V_SCALE/(DY(J+1)*DYF(J))
            VEC(I,J)=R1(I,K,J)
          END DO
        END DO
! IF using a subgrid model (LES) then add the eddy viscosity part implicitly
        IF (LES) THEN
        DO J=JSTART,JEND
          DO I=0,NXM
            MATL(I,J) = MATL(I,J) - TEMP2 * NU_T(I,K,J) 
     &                               / (DY(J)*DYF(J))
            MATD(I,J) = MATD(I,J) + TEMP2 * NU_T(I,K,J+1)
     &                              / (DY(J+1)*DYF(J))
     &                            + TEMP2 * NU_T(I,K,J)
     &                              / (DY(J)*DYF(J))
            MATU(I,J) = MATU(I,J) - TEMP2 * NU_T(I,K,J+1)
     &                             / (DY(J+1)*DYF(J))
          END DO
        END DO
        END IF

        IF (USE_MPI) THEN
! First, apply the boundary conditions
          CALL APPLY_BC_U1_MPI(MATL,MATD,MATU,VEC)
! If we are using MPI, split the implicit solve into forward and
! backward sweeps for efficiency
          CALL THOMAS_FORWARD_REAL_MPI(MATL,MATD,MATU,VEC,NY,NX)
          CALL THOMAS_BACKWARD_REAL_MPI(MATL,MATD,MATU,VEC,NY,NX)
        ELSE
C Else, we are running in serial mode
C Set the boundary conditions for U1
          CALL APPLY_BC_1_LOWER(MATL,MATD,MATU,VEC)
          CALL APPLY_BC_1_UPPER(MATL,MATD,MATU,VEC)
C Now, solve the tridiagonal system for U1(:,k,:)
          CALL THOMAS_REAL(MATL,MATD,MATU,VEC,NY+1,NXM)
        END IF

        DO J=JSTART-1,JEND+1
          DO I=0,NXM
            U1(I,K,J)=VEC(I,J)
          END DO
        END DO

! End do k
      END DO

! Initialize the matrix used to store implicit coefficients
      DO J=0,NY+1
        DO I=0,NXM
          MATL(I,J)=0.
          MATD(I,J)=1.
          MATU(I,J)=0.
          VEC(I,J)=0.
        END DO
      END DO 

C Solve for U3
C Note, here the matrix will be indexed from 1...NY+1 corresponding to U1(0:NY)
C Build the implicit system of equations for U3
      DO K=0,NZP-1
        DO J=JSTART,JEND
          DO I=0,NXM
            MATL(I,J)=-TEMP1*NU_V_SCALE/(DY(J)*DYF(J))
            MATD(I,J)=1.-TEMP1*NU_V_SCALE*(-1./(DY(J+1)*DYF(J))
     &         -1./(DY(J)*DYF(J)))
            MATU(I,J)=-TEMP1*NU_V_SCALE/(DY(J+1)*DYF(J))
            VEC(I,J)=R3(I,K,J)
          END DO
        END DO
! IF using a subgrid model (LES) then add the eddy viscosity part implicitly
        IF (LES) THEN
        DO J=JSTART,JEND
          DO I=0,NXM
            MATL(I,J) = MATL(I,J) - TEMP2 * NU_T(I,K,J)
     &                               / (DY(J)*DYF(J))
            MATD(I,J) = MATD(I,J) + TEMP2 * NU_T(I,K,J+1)
     &                              / (DY(J+1)*DYF(J))
     &                            + TEMP2 * NU_T(I,K,J)
     &                              / (DY(J)*DYF(J))
            MATU(I,J) = MATU(I,J) - TEMP2 * NU_T(I,K,J+1)
     &                             / (DY(J+1)*DYF(J))
          END DO
        END DO
        END IF

        IF (USE_MPI) THEN
! First, apply the boundary conditions
          CALL APPLY_BC_U3_MPI(MATL,MATD,MATU,VEC)
! If we are using MPI, split the implicit solve into forward and
! backward sweeps for efficiency
          CALL THOMAS_FORWARD_REAL_MPI(MATL,MATD,MATU,VEC,NY,NX)
          CALL THOMAS_BACKWARD_REAL_MPI(MATL,MATD,MATU,VEC,NY,NX)
        ELSE
C Else, we are running in serial mode
C Set the boundary conditions for U3
          CALL APPLY_BC_3_LOWER(MATL,MATD,MATU,VEC)
          CALL APPLY_BC_3_UPPER(MATL,MATD,MATU,VEC)
C Now, solve the tridiagonal system for U3(i,:,k)
          CALL THOMAS_REAL(MATL,MATD,MATU,VEC,NY+1,NXM)
        END IF

        DO J=JSTART-1,JEND+1
          DO I=0,NXM
            U3(I,K,J)=VEC(I,J)
          END DO
        END DO
! End do k
      END DO

C -- Done getting U1hat, U2hat, U3hat at new RK Step --

C If we are on the final RK step, optionally update the timestep
C based on the CFL criteria. This won't affect the current timestep
C since the TEMP1, etc. variables have already been set using
C the current timestep
      IF (VARIABLE_DT.and.(RK_STEP.eq.3)
     &      .and.(MOD(TIME_STEP,UPDATE_DT).EQ.0)) THEN
        CALL COURANT
      END IF

! Transform TH and U to Fourier Space 
      CALL FFT_XZ_TO_FOURIER(U1,CU1,0,NY+1)
      CALL FFT_XZ_TO_FOURIER(U2,CU2,0,NY+1)
      CALL FFT_XZ_TO_FOURIER(U3,CU3,0,NY+1)
      DO N=1,N_TH
        S1(:,:,:)=TH(:,:,:,N)
        CALL FFT_XZ_TO_FOURIER(S1,CS1,0,NY+1)
        CTH(:,:,:,N)=CS1(:,:,:)
      END DO

C Begin second step of the Fractional Step algorithm, making u divergence free
C The following subroutine projects Uhat onto divergence free space

      CALL REM_DIV_CHAN

C Now, phi is stored in CR1, use this to update the pressure field
C Note, here we divide by H_BAR since it was absorbed into PHI in REM_DIV
      DO J=JSTART,JEND
        DO K=0,TNKZ
          DO I=0,NXP-1
            CP(I,K,J)=CP(I,K,J)+CR1(I,K,J)/TEMP4
          END DO
        END DO
      END DO

      ! Fix disparities at the boundary due to the thoms algorithm in parallel
      IF (USE_MPI) THEN
         CALL GHOST_CHAN_MPI
      END IF

      RETURN
      END

C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE RK_CHAN_2
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
C Alternative time-stepping algorithm for the channel-flow case.
C This algorithm uses Crank-Nicolson for all viscous terms and 
C third order Runge-Kutta for all nonlinear terms
C INPUTS  (in Fourier space):  CUi, P, and (if k>1) CFi at (k-1)  (for i=1,2,3)
C OUTPUTS (in Fourier space):  CUi, P, and (if k<3) CFi at (k)
C Each RK step, there are 11 FFT calls. 11 storage variables are used.     
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|

      INCLUDE 'header'

      INTEGER I,J,K,N      
      REAL*8 TEMP1, TEMP2, TEMP3, TEMP4, TEMP5, UBULK

      ! STOP -----------------------------------
      IF (RANK.EQ.0) 
     &     write(*,*) ' RK_CHAN_2 not supported yet '
      call mpi_finalize(ierror)
      stop 
      ! ----------------------------------------


      RETURN
      END

C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE REM_DIV_CHAN
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      
C Compute varphi, store in variable CR1.
C Solves for phi in computational space
C H_BAR has been absorbed into PHI, so we are solving for H_BAR*PHI

      INCLUDE 'header'
      INTEGER I,J,K
 
C First, Initialize the matrix components
      DO J=0,NY+1
        DO I=0,NXP-1
          MATL_C(I,J)=0.
          MATD_C(I,J)=1.
          MATU_C(I,J)=0.
          VEC_C(I,J)=(0.,0.)
        END DO
      END DO

C The 2d FFT of Ui should have been taken and stored in CUi
C Solving for phi amounts to solving a tridiagonal system
C First, construct the system to be solved
      DO K=0,TNKZ
        DO J=1,NY
          DO I=0,NXP-1
            MATL_C(I,J)=1./(DY(J)*DYF(J))
            MATD_C(I,J)=-KX2(I)-KZ2(K)
     &         -1./(DY(J+1)*DYF(J))-1./(DY(J)*DYF(J))
            MATU_C(I,J)=1./(DY(J+1)*DYF(J))
          END DO
        END DO

C Now, create the RHS vector
        DO J=1,NY         
          DO I=0,NXP-1
            VEC_C(I,J)=(CIKX(I)*CU1(I,K,J) 
     &            + (CU2(I,K,J+1)-CU2(I,K,J))/DYF(J) 
     &            + CIKZ(K)*CU3(I,K,J))
          END DO
        END DO

        IF (USE_MPI) THEN
C If we are using the MPI package...
          CALL APPLY_BC_REM_DIV_MPI(MATL_C,MATD_C,MATU_C,VEC_C,K)
C First, do all forward sweeps
          CALL THOMAS_FORWARD_COMPLEX_MPI(MATL_C,MATD_C,MATU_C,VEC_C
     &                                 ,NY,NXP)
C Now, do the backward sweeps
          CALL THOMAS_BACKWARD_COMPLEX_MPI(MATL_C,MATD_C,MATU_C,VEC_C
     &                                 ,NY,NXP)
        ELSE
C Else we are running in serial mode
        DO I=0,NKX
          IF ((K.EQ.0).AND.(I.EQ.0)) THEN
C Use homogeneous dirichlet BCS for kx=kz=0 component at bottom wall
C Otherwise the matrix will be singular
            MATL_C(I,1)=0. 
            MATD_C(I,1)=1.
            MATU_C(I,1)=0.
            VEC_C(I,1)=(0.,0.)

            MATL_C(I,NY)=1.
            MATD_C(I,NY)=-1.
            MATU_C(I,NY)=0.
            VEC_C(I,NY)=(0.,0.)
          ELSE
C Use Dirichlet boundary conditions, dp/dz=0 at walls
            MATL_C(I,1)=0.
            MATD_C(I,1)=1.
            MATU_C(I,1)=-1.
            VEC_C(I,1)=(0.,0.)

            MATL_C(I,NY)=1.
            MATD_C(I,NY)=-1.
            MATU_C(I,NY)=0.
            VEC_C(I,NY)=(0.,0.)
          END IF
        END DO
C Now solve the tridiagonal system for phi, store in CR1
        CALL THOMAS_COMPLEX(MATL_C,MATD_C,MATU_C,VEC_C,NY,NKX)
        END IF


        DO J=1,NY
          DO I=0,NXP-1
            CR1(I,K,J)=VEC_C(I,J)
          END DO
        END DO

      END DO

C Now, Solve for CUi, the divergenceless velocity field
      DO J=1,NY
        DO K=0,TNKZ
          DO I=0,NXP-1
            CU1(I,K,J)=CU1(I,K,J)-CIKX(I)*CR1(I,K,J)
            CU3(I,K,J)=CU3(I,K,J)-CIKZ(K)*CR1(I,K,J)           
          END DO
        END DO
      END DO
      DO J=2,NY
        DO K=0,TNKZ
          DO I=0,NXP-1
            CU2(I,K,J)=CU2(I,K,J)-(CR1(I,K,J)
     &             -CR1(I,K,J-1))/DY(J)
          END DO
        END DO
      END DO

      RETURN
      END

C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE POISSON_P_CHAN
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
C We have CUi, need to compute CP.  Solve tridiagonal system exactly

      INCLUDE 'header'

      INTEGER I,J,K,N
      
	if (flavor.eq.'Basic') then
      IF (RANK.EQ.0) 
     &          WRITE(*,*) 'COMPUTING CP FROM CUI'
	end if


C First, construct the RHS vector, (dui/dxj)(duj/dxi) 
      DO J=2,NY
        DO K=0,TNKZ
          DO I=0,NXP-1 ! NKX
            CF1(I,K,J)=CIKX(I)*CU1(I,K,J)
            CF2(I,K,J)=(CU2(I,K,J+1)-CU2(I,K,J))/DYF(J)
            CF3(I,K,J)=CIKZ(K)*CU3(I,K,J)
          END DO
        END DO
      END DO

      CALL FFT_XZ_TO_PHYSICAL(CF1,F1,0,NY+1)
      CALL FFT_XZ_TO_PHYSICAL(CF2,F2,0,NY+1)
      CALL FFT_XZ_TO_PHYSICAL(CF3,F3,0,NY+1)
      
      DO J=2,NY
        DO K=0,NZP-1
          DO I=0,NXM
            F1(I,K,J)=F1(I,K,J)**2.
            F2(I,K,J)=F2(I,K,J)**2.
            F3(I,K,J)=F3(I,K,J)**2.
          END DO
        END DO
      END DO

      CALL FFT_XZ_TO_FOURIER(F1,CF1,0,NY+1)
      CALL FFT_XZ_TO_FOURIER(F2,CF2,0,NY+1)
      CALL FFT_XZ_TO_FOURIER(F3,CF3,0,NY+1)

C Now we have the diagonal terms, add to the rhs term
      DO J=2,NY
        DO K=0,TNKZ
          DO I=0,NXP-1 ! NKX
            CS1(I,K,J)=CF1(I,K,J)+CF2(I,K,J)+CF3(I,K,J)
          END DO
        END DO
      END DO

C Now get the first of the off-diagonal terms
      DO J=2,NY
        DO K=0,TNKZ
          DO I=0,NXP-1 ! NKX
            CF1(I,K,J)=(CU1(I,K,J+1)-CU1(I,K,J-1))/(2.*DYF(J))
            CF2(I,K,J)=CIKX(I)*0.5*(CU2(I,K,J)+CU2(I,K,J+1))
          END DO
        END DO
      END DO

      CALL FFT_XZ_TO_PHYSICAL(CF1,F1,0,NY+1)
      CALL FFT_XZ_TO_PHYSICAL(CF2,F2,0,NY+1)

C Compute product
      DO J=2,NY
        DO K=0,NZP-1
          DO I=0,NXM
            F1(I,K,J)=2.*F1(I,K,J)*F2(I,K,J)
          END DO
        END DO
      END DO
      
      CALL FFT_XZ_TO_FOURIER(F1,CF1,0,NY+1)

C Add to RHS term
      DO J=2,NY
        DO K=0,TNKZ
          DO I=0,NXP-1 ! NKX 
            CS1(I,K,J)=CS1(I,K,J)+CF1(I,K,J)
          END DO
        END DO
      END DO

C Now get the second of the off-diagonal terms
      DO J=2,NY
        DO K=0,TNKZ
          DO I=0,NXP-1 ! NKX
            CF1(I,K,J)=(CU3(I,K,J+1)-CU3(I,K,J-1))/(2.*DYF(J))
            CF2(I,K,J)=CIKZ(K)*0.5*(CU2(I,K,J)+CU2(I,K,J+1))
          END DO
        END DO
      END DO

C Convert to Physical space
      CALL FFT_XZ_TO_PHYSICAL(CF1,F1,0,NY+1)
      CALL FFT_XZ_TO_PHYSICAL(CF2,F2,0,NY+1)

C Compute product
      DO J=2,NY
        DO K=0,NZP-1
          DO I=0,NXM
            F1(I,K,J)=2.*F1(I,K,J)*F2(I,K,J)
          END DO
        END DO
      END DO

      CALL FFT_XZ_TO_FOURIER(F1,CF1,0,NY+1)

C Add to RHS term
      DO J=2,NY
        DO K=0,TNKZ
          DO I=0,NXP-1  ! NKX
            CS1(I,K,J)=CS1(I,K,J)+CF1(I,K,J)
          END DO
        END DO
      END DO

C Now get the third of the off-diagonal terms
      DO J=2,NY
        DO K=0,TNKZ
          DO I=0,NXP-1  ! NKX
            CF1(I,K,J)=CIKZ(K)*CU1(I,K,J)
            CF2(I,K,J)=CIKX(I)*CU3(I,K,J)
          END DO
        END DO
      END DO

      CALL FFT_XZ_TO_PHYSICAL(CF1,F1,0,NY+1)
      CALL FFT_XZ_TO_PHYSICAL(CF2,F2,0,NY+1)
      
C Compute product
      DO J=2,NY
        DO K=0,NZP-1
          DO I=0,NXM
            F1(I,K,J)=2.*F1(I,K,J)*F2(I,K,J)
          END DO
        END DO
      END DO

      CALL FFT_XZ_TO_FOURIER(F1,CF1,0,NY+1)

C Add to RHS term
      DO J=2,NY
        DO K=0,TNKZ
          DO I=0,NXP-1  ! NKX
            CS1(I,K,J)=CS1(I,K,J)+CF1(I,K,J)
          END DO
        END DO
      END DO     
    
C Finally, if the buoyancy force is active, then we need to add
C the contribution of the density to the pressure.  Note that the
C plane averaged density and the corresponding hydrostatic part of the
C pressure have been cancelled, so skip the 0,0 mode
      DO N=1,N_TH
      DO J=2,NY
        DO K=0,TNKZ
          DO I=0,NXP-1  
            IF ((RANKZ.NE.0).OR.(I.NE.0).or.(K.NE.0)) THEN
              CS1(I,K,J)=CS1(I,K,J)+RI(N)*
     &          (CTH(I,K,J+1,N)-CTH(I,K,J-1,N))/(GYF(J+1)-GYF(J-1))
            END IF
          END DO
        END DO
      END DO
      END DO

C Now, the RHS term should be stored in CS1     

C Construct the tridiagonal system in Fourier space to solve for CP
C First, zero the vectors
      DO J=0,NY+1
        DO I=0,NXP-1
          MATL_C(I,J)=0.d0
          MATD_C(I,J)=1.d0
          MATU_C(I,J)=0.d0
          VEC_C(I,J)=(0.,0.)
        END DO
      END DO

      DO K=0,TNKZ
        DO J=2,NY
          DO I=0,NXP-1
            MATL_C(I,J)=1./(DY(J)*DYF(J))
            MATD_C(I,J)=-KX2(I)-KZ2(K)-1./(DY(J+1)*DYF(J))
     &                    -1./(DY(J)*DYF(J))
            MATU_C(I,J)=1./(DY(J+1)*DYF(J))   
            VEC_C(I,J)=-1.*CS1(I,K,J)
          END DO
        END DO

        IF (USE_MPI) THEN
          CALL APPLY_BC_POISSON_MPI(MATL_C,MATD_C,MATU_C,VEC_C,K)
C First, do the forward sweeps
          CALL THOMAS_FORWARD_COMPLEX_MPI(MATL_C,MATD_C,MATU_C,VEC_C
     &                                 ,NY,NXP)
C Now, do the backwared sweeps to put the solution in VEC_C
          CALL THOMAS_BACKWARD_COMPLEX_MPI(MATL_C,MATD_C,MATU_C,VEC_C
     &                                  ,NY,NXP)
        ELSE
C Else we are running in serial mode
C Apply BCs
        DO I=0,NKX
C Use dirichlet boundary condition at the lower wall to
C prevent the tridiagonal matrix from becomming singular for i,k=0
          IF ((I.EQ.0).AND.(K.EQ.0)) THEN
            MATD_C(I,1)=1.
            MATU_C(I,1)=0.
            VEC_C(I,1)=(0.,0.)
            MATD_C(I,NY)=-1.
            MATL_C(I,NY)=1.
            VEC_C(I,NY)=(0.,0.)
          ELSE
! Here, apply Neumann boundary conditions (dp/dz=0) at the walls
            MATD_C(I,1)=1.
            MATU_C(I,1)=-1.
            VEC_C(I,1)=(0.,0.)
            MATD_C(I,NY)=-1.
            MATL_C(I,NY)=1.
            VEC_C(I,NY)=(0.,0.)
          END IF
        END DO
C Now, solve for CP
        CALL THOMAS_COMPLEX(MATL_C,MATD_C,MATU_C,VEC_C,NY,NKX)
        END IF

        DO J=1,NY
          DO I=0,NXP-1
            CP(I,K,J)=VEC_C(I,J)
          END DO
        END DO
      END DO

      RETURN
      END

C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE INPUT_CHAN
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      INCLUDE 'header'
      REAL    VERSION, CURRENT_VERSION
      INTEGER I,J,K,N

! Read in input parameters specific for channel flow case
      OPEN (11,file='input_chan.dat',form='formatted',status='old')      
C Read input file.

      CURRENT_VERSION=2.0
      READ(11,*)
      READ(11,*)
      READ(11,*)
      READ(11,*)
      READ(11,*) VERSION
      IF (VERSION .NE. CURRENT_VERSION) 
     &         STOP 'Wrong input data format input_chan'
      READ(11,*)
      READ(11,*) TIME_AD_METH
      READ(11,*) 
      READ(11,*) LES_MODEL_TYPE
      READ(11,*)
      READ(11,*) IC_TYPE, KICK
      READ(11,*)
      READ(11,*) I_RO
      READ(11,*) 
      READ(11,*) GRAV_X, GRAV_Y, GRAV_Z
      READ(11,*)
      READ(11,*) F_TYPE, UBULK0, PX0, OMEGA0, AMP_OMEGA0
      READ(11,*)
      READ(11,*) U_BC_YMIN, U_BC_YMIN_C1, U_BC_YMIN_C2, U_BC_YMIN_C3
      READ(11,*) 
      READ(11,*) V_BC_YMIN, V_BC_YMIN_C1, V_BC_YMIN_C2, V_BC_YMIN_C3
      READ(11,*)
      READ(11,*) W_BC_YMIN, W_BC_YMIN_C1, W_BC_YMIN_C2, W_BC_YMIN_C3
      READ(11,*)
      READ(11,*) U_BC_YMAX, U_BC_YMAX_C1, U_BC_YMAX_C2, U_BC_YMAX_C3
      READ(11,*)
      READ(11,*) V_BC_YMAX, V_BC_YMAX_C1, V_BC_YMAX_C2, V_BC_YMAX_C3
      READ(11,*)
      READ(11,*) W_BC_YMAX, W_BC_YMAX_C1, W_BC_YMAX_C2, W_BC_YMAX_C3
      READ(11,*)
! Read in boundary conditions and background gradients for the N_TH scalars
      DO N=1,N_TH
        READ(11,*)
        READ(11,*) TH_BC_YMIN(N),TH_BC_YMIN_C1(N),TH_BC_YMIN_C2(N)
     &             ,TH_BC_YMIN_C3(N)
        READ(11,*)
        READ(11,*) TH_BC_YMAX(N),TH_BC_YMAX_C1(N),TH_BC_YMAX_C2(N)
     &             ,TH_BC_YMAX_C3(N)
        READ(11,*)
        READ(11,*) DRHODX(N), DRHODZ(N)
      END DO

      RETURN
      END


C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE INIT_CHAN_MOVIE
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      INCLUDE 'header'
      CHARACTER*55 FNAME
      INTEGER I,J,K

! Set parameters for writing "movie" files with 2d slices

! First, read in the center of each 2D slice
            OPEN(unit=650,file='MOVIE.dat',status='old',
     &           form='formatted')
            READ(650, *) XcMovie, YcMovie, ZcMovie
            CLOSE(650)

!     Get the indices
            NxMovie=int(XcMovie*NX/LX)
            NyMovie=nint(YcMovie*(((NY-1)*NPROCS)/LY))
            NzMovie=int(ZcMovie*NZ/LZ)

            RankZMovie = int(NzMovie/NZP)
            NzMovie    = NzMovie-RankZMovie*NZP

            RankYMovie=-1
            IF (GYF(JSTART).LE.YcMovie.and.GYF(JEND+1).GT.YcMovie) THEN
               RankYMovie=RANKY
               I=1
               do while(.not.
     &              (GYF(I).LE.YcMovie .and. GYF(I+1).GT.YcMovie))
                  I=I+1
               end do
               NyMovie=I;
            END IF

            if (RANKY.eq.RankYMovie .and. RANKZ.eq.RankZMovie) then
               write(*,*) 'Movie Parameters, RANK:', RANK
               write(*,*) '    Xc: ', GX(NxMovie), ' (NxMovie: ',
     &              NxMovie, ')'
               write(*,*) '    Yc: ', GYF(NyMovie)
     &              , ' (NyMovie: ', NyMovie, ')'
               write(*,*) '    Zc: ', GZ(RankZMovie*NZP+NzMovie),
     &              ' (NzMovie: ', RankZMovie*NZP+NzMovie, ')'

            END IF
 
      RETURN
      END

C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE CREATE_GRID_CHAN
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      INCLUDE 'header'
      CHARACTER*55 FNAME
      INTEGER I,J,K

         IF (RANK.EQ.0) 
     &     WRITE (6,*) 'Fourier in X'
         DO I=0,NX
           GX(I)=(I*LX)/NX
           DX(I)=LX/NX
           IF (VERBOSITY .GT. 3 .AND. RANK.EQ.0) 
     &          WRITE(6,*) 'GX(',I,') = ',GX(I)
         END DO
         IF (RANK.EQ.0) 
     &        WRITE (6,*) 'Fourier in Z'
         DO K=0,NZ
           GZ(K)=(K*LZ)/NZ
           DZ(K)=LZ/NZ
           IF (RANK.EQ.0 .AND. VERBOSITY .GT. 3) 
     &          WRITE(6,*) 'GZ(',K,') = ',GZ(K)
         END DO
         IF (RANK.EQ.0) 
     &        WRITE (6,*) 'Finite-difference in Y'

         IF (RANK.EQ.0) 
     &        write(*,*) 'USE_MPI: ',USE_MPI

         FNAME='grid.h5'
         if (FNAME(len_trim(FNAME)-2:len_trim(FNAME)).eq.".h5") then
#ifdef HDF5
            if (USE_MPI) then
            call mpi_barrier(MPI_COMM_WORLD,ierror)
            end if
            call ReadGridHDF5(FNAME,2)
#else
            IF (RANK.EQ.0) THEN
            write(*,*) ' **** ERROR ******************************'
            write(*,*) ' Program not compiled with HDF5 libraries.'
            END IF
            stop 
#endif
         else 
         IF (USE_MPI) THEN
           FNAME='./ygrid'//trim(MPI_IO_NUM)//'.txt'
           IF (RANK.EQ.0) THEN
           write(*,*) 'FNAME: ',FNAME
           write(*,*) 'MPI_IO_NUM: ****',trim(MPI_IO_NUM),'*****'
           END IF
         END IF

         OPEN (30,file=FNAME,form='formatted',status='old')
         READ (30,*) NY_T
C Check to make sure that grid file is the correct dimensions
         IF (NY_T.ne.NY) THEN
           IF (RANK.EQ.0) 
     &           WRITE(6,*) 'NY, NY_T',NY,NY_T
           STOP 'Error: ygrid.txt wrong dimensions'
         END IF
         DO J=1,NY+1
           READ(30,*) GY(j)
           IF (VERBOSITY .GT. 3 .AND. RANK.EQ.0) 
     &          WRITE(6,*) 'GY(',J,') = ',GY(J)
         END DO
         DO J=1,NY
           READ(30,*) GYF(j)
           IF (RANK.EQ.0 .AND. VERBOSITY .GT. 3)
     &          WRITE(6,*) 'GYF(',J,') = ',GYF(J)
         END DO
         CLOSE(30)

         IF (USE_MPI) THEN
           CALL GHOST_GRID_MPI
         ELSE
C Define ghost cells
           GYF(0)=2.d0*GYF(1)-GYF(2)
           GYF(NY+1)=2.d0*GYF(NY)-GYF(NYM)
           GY(0)=2.d0*GY(1)-GY(2)
         END IF
         end if

C Define grid spacing
         DO J=1,NY+1
           DY(J)=(GYF(J)-GYF(J-1))
         END DO
         DO J=1,NY
           DYF(J)=(GY(J+1)-GY(J))
         END DO
         DYF(NY+1)=DYF(NY)

         RETURN 
         END


 
C----*|--.---------.---------.---------.---------.---------.---------.-|------
      SUBROUTINE APPLY_BC_1_LOWER(MATL,MATD,MATU,VEC)
C----*|--.---------.---------.---------.---------.---------.---------.-|-----
      INCLUDE 'header'
      INTEGER I

C Bottom Wall:
      IF (U_BC_YMIN.EQ.0) THEN
C Dirichlet
        DO I=0,NXM
          MATL(I,0)=0. 
          MATD(I,0)=1.
          MATU(I,0)=0.                   
          VEC(I,0)=0.

          MATL(I,1)=0. 
          MATD(I,1)=1.
          MATU(I,1)=0.                   
          VEC(I,1)=U_BC_YMIN_C1 
        END DO
      ELSE
C Neumann
        DO I=0,NXM
          MATL(I,0)=0.
          MATD(I,0)=1.
          MATU(I,0)=0.
          VEC(I,0)=0.
        END DO
        DO I=0,NXM
          MATL(I,1)=0.
          MATD(I,1)=-1.
          MATU(I,1)=1.
          VEC(I,1)=DY(2)*U_BC_YMIN_C1
        END DO

      END IF

      RETURN 
      END

C----*|--.---------.---------.---------.---------.---------.---------.-|------
      SUBROUTINE APPLY_BC_1_LOWER_C(MATL_C,MATD_C,MATU_C,VEC_C)
C----*|--.---------.---------.---------.---------.---------.---------.-|-----
      INCLUDE 'header'
      INTEGER I

C Bottom Wall:
      IF (U_BC_YMIN.EQ.0) THEN
C Dirichlet
        DO I=0,NKX
          MATL_C(I,0)=0. 
          MATD_C(I,0)=1.
          MATU_C(I,0)=0.                   
          VEC_C(I,0)=0.

          MATL_C(I,1)=0. 
          MATD_C(I,1)=1.
          MATU_C(I,1)=0.                   
          VEC_C(I,1)=U_BC_YMIN_C1 
        END DO
      ELSE
C Neumann
        DO I=0,NKX
          MATL_C(I,0)=0.
          MATD_C(I,0)=1.
          MATU_C(I,0)=0.
          VEC_C(I,0)=0.
        END DO
        DO I=0,NKX
          MATL_C(I,1)=0.
          MATD_C(I,1)=-1.
          MATU_C(I,1)=1.
          VEC_C(I,1)=DY(2)*U_BC_YMIN_C1
        END DO

      END IF

      RETURN 
      END


C----*|--.---------.---------.---------.---------.---------.---------.-|----
      SUBROUTINE APPLY_BC_1_UPPER(MATL,MATD,MATU,VEC)
C----*|--.---------.---------.---------.---------.---------.---------.-|--
      INCLUDE 'header'
      INTEGER I

C Top wall
      IF (U_BC_YMAX.EQ.0) THEN
C Dirichlet
        DO I=0,NXM
          MATL(I,NY+1)=0.
          MATD(I,NY+1)=1.
          MATU(I,NY+1)=0.
          VEC(I,NY+1)=0.

          MATL(I,NY)=0.
          MATD(I,NY)=1.
          MATU(I,NY)=0.
          VEC(I,NY)=U_BC_YMAX_C1
        END DO
      ELSE
C Neumann
        DO I=0,NXM
          MATL(I,NY)=-1.
          MATD(I,NY)=1.
          MATU(I,NY)=0.
          VEC(I,NY)=DY(NY)*U_BC_YMAX_C1
        END DO
        DO I=0,NXM
          MATL(I,NY+1)=0.
          MATD(I,NY+1)=1.
          MATU(I,NY+1)=0.
          VEC(I,NY+1)=0.
        END DO

      END IF

      RETURN
      END

C----*|--.---------.---------.---------.---------.---------.---------.-|----
      SUBROUTINE APPLY_BC_1_UPPER_C(MATL_C,MATD_C,MATU_C,VEC_C)
C----*|--.---------.---------.---------.---------.---------.---------.-|--
      INCLUDE 'header'
      INTEGER I

C Top wall
      IF (U_BC_YMAX.EQ.0) THEN
C Dirichlet
        DO I=0,NKX
          MATL_C(I,NY+1)=0.
          MATD_C(I,NY+1)=1.
          MATU_C(I,NY+1)=0.
          VEC_C(I,NY+1)=0.

          MATL_C(I,NY)=0.
          MATD_C(I,NY)=1.
          MATU_C(I,NY)=0.
          VEC_C(I,NY)=U_BC_YMAX_C1
        END DO
      ELSE
C Neumann
        DO I=0,NKX
          MATL_C(I,NY+1)=0.
          MATD_C(I,NY+1)=1.
          MATU_C(I,NY+1)=0.
          VEC_C(I,NY+1)=0.
        END DO
        DO I=0,NKX
          MATL_C(I,NY)=-1.
          MATD_C(I,NY)=1.
          MATU_C(I,NY)=0.
          VEC_C(I,NY)=DY(NY)*U_BC_YMAX_C1
        END DO

      END IF

      RETURN
      END


C----*|--.---------.---------.---------.---------.---------.---------.-|---
      SUBROUTINE APPLY_BC_2_LOWER(MATL,MATD,MATU,VEC)
C----*|--.---------.---------.---------.---------.---------.---------.-|--
      INCLUDE 'header'
      INTEGER I

C Bottom Wall:
      IF (V_BC_YMIN.EQ.0) THEN
C Dirichlet
        DO I=0,NXM
          MATL(I,1)=0.d0 
          MATD(I,1)=1.d0
          MATU(I,1)=0.d0                   
          VEC(I,1)=V_BC_YMIN_C1 

          MATL(I,2)=0.d0 
          MATD(I,2)=1.d0
          MATU(I,2)=0.d0                   
          VEC(I,2)=V_BC_YMIN_C1 
        END DO
      ELSE IF (V_BC_YMIN.EQ.1) THEN
C Neumann
        DO I=0,NXM
          MATD(I,1)=-1.d0
          MATU(I,1)=1.d0
          MATL(I,1)=0.d0
          VEC(I,1)=DYF(1)*V_BC_YMIN_C1
        END DO
      END IF

C The following is only a placeholder, this row is used for U1 and U3
      DO I=0,NXM
        MATL(I,0) = 0.
        MATD(I,0) = 1.
        MATU(I,0) = 0.
        VEC(I,0) = 0.
      END DO

      RETURN
      END

C----*|--.---------.---------.---------.---------.---------.---------.-|---
      SUBROUTINE APPLY_BC_2_LOWER_C(MATL_C,MATD_C,MATU_C,VEC_C)
C----*|--.---------.---------.---------.---------.---------.---------.-|--
      INCLUDE 'header'
      INTEGER I

C Bottom Wall:
      IF (V_BC_YMIN.EQ.0) THEN
C Dirichlet
        DO I=0,NKX
          MATL_C(I,1)=0.d0 
          MATD_C(I,1)=1.d0
          MATU_C(I,1)=0.d0                   
          VEC_C(I,1)=V_BC_YMIN_C1 

          MATL_C(I,2)=0.d0 
          MATD_C(I,2)=1.d0
          MATU_C(I,2)=0.d0                   
          VEC_C(I,2)=V_BC_YMIN_C1 
        END DO
      ELSE IF (V_BC_YMIN.EQ.1) THEN
C Neumann
        DO I=0,NKX
          MATD_C(I,1)=-1.d0
          MATU_C(I,1)=1.d0
          MATL_C(I,1)=0.d0
          VEC_C(I,1)=DYF(1)*V_BC_YMIN_C1
        END DO
      END IF

C The following is only a placeholder, this row is used for U1 and U3
      DO I=0,NKX
        MATL_C(I,0) = 0.
        MATD_C(I,0) = 1.
        MATU_C(I,0) = 0.
        VEC_C(I,0) = 0.
      END DO

      RETURN
      END

 
C----*|--.---------.---------.---------.---------.---------.---------.-|--
      SUBROUTINE APPLY_BC_2_UPPER(MATL,MATD,MATU,VEC)
C----*|--.---------.---------.---------.---------.---------.---------.-|--
      INCLUDE 'header'
      INTEGER I
C Top wall
      IF (V_BC_YMAX.EQ.0) THEN
C Dirichlet
        DO I=0,NXM
          MATL(I,NY+1)=0.
          MATD(I,NY+1)=1.
          MATU(I,NY+1)=0.
          VEC(I,NY+1)=V_BC_YMAX_C1
          
          MATL(I,NY)=0.
          MATD(I,NY)=1.
          MATU(I,NY)=0.
          VEC(I,NY)=V_BC_YMAX_C1
        END DO
      ELSE IF (V_BC_YMAX.EQ.1) THEN
C Neumann
        DO I=0,NXM
          MATL(I,NY+1)=-1.
          MATD(I,NY+1)=1.
          MATU(I,NY+1)=0.
          VEC(I,NY+1)=DYF(NY)*V_BC_YMAX_C1
        END DO
      END IF      
      RETURN
      END

C----*|--.---------.---------.---------.---------.---------.---------.-|--
      SUBROUTINE APPLY_BC_2_UPPER_C(MATL_C,MATD_C,MATU_C,VEC_C)
C----*|--.---------.---------.---------.---------.---------.---------.-|--
      INCLUDE 'header'
      INTEGER I
C Top wall
      IF (V_BC_YMAX.EQ.0) THEN
C Dirichlet
        DO I=0,NKX
          MATL_C(I,NY+1)=0.
          MATD_C(I,NY+1)=1.
          MATU_C(I,NY+1)=0.
          VEC_C(I,NY+1)=V_BC_YMAX_C1
          
          MATL_C(I,NY)=0.
          MATD_C(I,NY)=1.
          MATU_C(I,NY)=0.
          VEC_C(I,NY)=V_BC_YMAX_C1
        END DO
      ELSE IF (V_BC_YMAX.EQ.1) THEN
C Neumann
        DO I=0,NKX
          MATL_C(I,NY+1)=-1.
          MATD_C(I,NY+1)=1.
          MATU_C(I,NY+1)=0.
          VEC_C(I,NY+1)=DYF(NY)*V_BC_YMAX_C1
        END DO      
      END IF
      RETURN
      END

C----*|--.---------.---------.---------.---------.---------.---------.-|--
      SUBROUTINE APPLY_BC_3_LOWER(MATL,MATD,MATU,VEC)
C----*|--.---------.---------.---------.---------.---------.---------.-|--
      INCLUDE 'header'
      INTEGER I

C Bottom Wall:
      IF (W_BC_YMIN.EQ.0) THEN
C Dirichlet
        DO I=0,NXM
          MATL(I,0)=0. 
          MATD(I,0)=1.
          MATU(I,0)=0.                   
          VEC(I,0)=0.

          MATL(I,1)=0. 
          MATD(I,1)=1.
          MATU(I,1)=0.                   
          VEC(I,1)=W_BC_YMIN_C1
        END DO
      ELSE
C Neumann
        DO I=0,NXM
          MATL(I,0)=0.
          MATD(I,0)=1.
          MATU(I,0)=0.
          VEC(I,0)=0.
        END DO
        DO I=0,NXM
          MATL(I,1)=0.
          MATD(I,1)=-1.
          MATU(I,1)=1.
          VEC(I,1)=DY(2)*W_BC_YMIN_C1
        END DO

      END IF

      RETURN
      END

C----*|--.---------.---------.---------.---------.---------.---------.-|--
      SUBROUTINE APPLY_BC_3_LOWER_C(MATL_C,MATD_C,MATU_C,VEC_C)
C----*|--.---------.---------.---------.---------.---------.---------.-|--
      INCLUDE 'header'
      INTEGER I

C Bottom Wall:
      IF (W_BC_YMIN.EQ.0) THEN
C Dirichlet
        DO I=0,NKX
          MATL_C(I,0)=0. 
          MATD_C(I,0)=1.
          MATU_C(I,0)=0.                   
          VEC_C(I,0)=0.

          MATL_C(I,1)=0. 
          MATD_C(I,1)=1.
          MATU_C(I,1)=0.                   
          VEC_C(I,1)=W_BC_YMIN_C1
        END DO
      ELSE
C Neumann
        DO I=0,NKX
          MATL_C(I,0)=0.
          MATD_C(I,0)=1.
          MATU_C(I,0)=0.
          VEC_C(I,0)=0.
        END DO
        DO I=0,NKX
          MATL_C(I,1)=0.
          MATD_C(I,1)=-1.
          MATU_C(I,1)=1.
          VEC_C(I,1)=DY(2)*W_BC_YMIN_C1
        END DO

      END IF

      RETURN
      END

C----*|--.---------.---------.---------.---------.---------.---------.-|--
      SUBROUTINE APPLY_BC_3_UPPER(MATL,MATD,MATU,VEC)
C----*|--.---------.---------.---------.---------.---------.---------.-|--
      INCLUDE 'header'
      INTEGER I

C Top wall
      IF (W_BC_YMAX.EQ.0) THEN
C Dirichlet
        DO I=0,NXM
          MATL(I,NY+1)=0.
          MATD(I,NY+1)=1.
          MATU(I,NY+1)=0.
          VEC(I,NY+1)=0.

          MATL(I,NY)=0.
          MATD(I,NY)=1.
          MATU(I,NY)=0.
          VEC(I,NY)=W_BC_YMAX_C1
        END DO
      ELSE
C Neumann
        DO I=0,NXM
          MATL(I,NY)=-1.
          MATD(I,NY)=1.
          MATU(I,NY)=0.
          VEC(I,NY)=DY(NY)*W_BC_YMAX_C1
        END DO
        DO I=0,NXM
          MATL(I,NY+1)=0.
          MATD(I,NY+1)=1.
          MATU(I,NY+1)=0.
          VEC(I,NY+1)=0.
        END DO

      END IF

      RETURN 
      END


C----*|--.---------.---------.---------.---------.---------.---------.-|--
      SUBROUTINE APPLY_BC_3_UPPER_C(MATL_C,MATD_C,MATU_C,VEC_C)
C----*|--.---------.---------.---------.---------.---------.---------.-|--
      INCLUDE 'header'
      INTEGER I

C Top wall
      IF (W_BC_YMAX.EQ.0) THEN
C Dirichlet
        DO I=0,NKX
          MATL_C(I,NY+1)=0.
          MATD_C(I,NY+1)=1.
          MATU_C(I,NY+1)=0.
          VEC_C(I,NY+1)=0.

          MATL_C(I,NY)=0.
          MATD_C(I,NY)=1.
          MATU_C(I,NY)=0.
          VEC_C(I,NY)=W_BC_YMAX_C1
        END DO
      ELSE
C Neumann
        DO I=0,NKX
          MATL_C(I,NY)=-1.
          MATD_C(I,NY)=1.
          MATU_C(I,NY)=0.
          VEC_C(I,NY)=DY(NY)*W_BC_YMAX_C1
        END DO
        DO I=0,NKX
          MATL_C(I,NY+1)=0.
          MATD_C(I,NY+1)=1.
          MATU_C(I,NY+1)=0.
          VEC_C(I,NY+1)=0.
        END DO

      END IF

      RETURN
      END

C----*|--.---------.---------.---------.---------.---------.---------.-|--
      subroutine APPLY_BC_TH_LOWER(MATL,MATD,MATU,VEC,N)
C----*|--.---------.---------.---------.---------.---------.---------.-|--
      include 'header'
      integer i,N
! Bottom Wall:
      if (TH_BC_YMIN(N).eq.0) then
! Dirichlet
        do i=0,NXM
          MATL(i,0)=0. 
          MATD(i,0)=1.
          MATU(i,0)=0.                   
          VEC(i,0)=0.

          MATL(i,1)=0. 
          MATD(i,1)=1.
          MATU(i,1)=0.                   
          VEC(i,1)=TH_BC_YMIN_C1(N)
        end do
      else
! Neumann
! NOTE: BC enforced at GY(2)
        do i=0,NXM
          MATL(i,1)=0.
          MATD(i,1)=-1.
          MATU(i,1)=1.
          VEC(i,1)=DY(2)*TH_BC_YMIN_C1(N)
        end do
        do i=0,NXM
          MATL(i,0)=0.
          MATD(i,0)=1.
          MATU(i,0)=0.
          VEC(i,0)=0.
        end do
      end if
      RETURN
      END

C----*|--.---------.---------.---------.---------.---------.---------.-|--
      subroutine APPLY_BC_TH_LOWER_C(MATL_C,MATD_C,MATU_C,VEC_C,N)
C----*|--.---------.---------.---------.---------.---------.---------.-|--
      include 'header'
      integer i,N
! Bottom Wall:
      if (TH_BC_YMIN(N).eq.0) then
! Dirichlet
        do i=0,NKX
          MATL_C(i,0)=0. 
          MATD_C(i,0)=1.
          MATU_C(i,0)=0.                   
          VEC_C(i,0)=0.

          MATL_C(i,1)=0. 
          MATD_C(i,1)=1.
          MATU_C(i,1)=0.                   
          VEC_C(i,1)=TH_BC_YMIN_C1(N)
        end do
      else
! Neumann
! NOTE: BC enforced at GY(2)
        do i=0,NXM
          MATL(i,1)=0.
          MATD(i,1)=-1.
          MATU(i,1)=1.
          VEC(i,1)=DY(2)*TH_BC_YMIN_C1(N)
        end do
        do i=0,NXM
          MATL(i,0)=0.
          MATD(i,0)=1.
          MATU(i,0)=0.
          VEC(i,0)=0.
        end do

      end if
      RETURN
      END

C----*|--.---------.---------.---------.---------.---------.---------.-|--
      subroutine APPLY_BC_TH_UPPER(MATL,MATD,MATU,VEC,N)
C----*|--.---------.---------.---------.---------.---------.---------.-|--
      include 'header'
      integer i,N
! Top wall
      if (TH_BC_YMAX(N).eq.0) then
! Dirichlet
        do i=0,NXM
          MATL(i,NY+1)=0.
          MATD(i,NY+1)=1.
          MATU(i,NY+1)=0.
          VEC(i,NY+1)=0.

          MATL(i,NY)=0.
          MATD(i,NY)=1.
          MATU(i,NY)=0.
          VEC(i,NY)=TH_BC_YMAX_C1(N)
        end do
      else
! Neumann
! NOTE: BC enforced at GY(NY)
        do i=0,NXM
          MATL(i,NY)=-1.
          MATD(i,NY)=1.
          MATU(i,NY)=0.
          VEC(i,NY)=DY(NY)*TH_BC_YMAX_C1(N)
        end do
        do i=0,NXM
          MATL(i,NY+1)=0.
          MATD(i,NY+1)=1.
          MATU(i,NY+1)=0.
          VEC(i,NY+1)=0.
        end do

      end if
      return
      end

C----*|--.---------.---------.---------.---------.---------.---------.-|--
      subroutine APPLY_BC_TH_UPPER_C(MATL_C,MATD_C,MATU_C,VEC_C,N)
C----*|--.---------.---------.---------.---------.---------.---------.-|--
      include 'header'
      integer i,N
! Top wall
      if (TH_BC_YMAX(N).eq.0) then
! Dirichlet
        do i=0,NKX
          MATL_C(i,NY+1)=0.
          MATD_C(i,NY+1)=1.
          MATU_C(i,NY+1)=0.
          VEC_C(i,NY+1)=0.

          MATL_C(i,NY)=0.
          MATD_C(i,NY)=1.
          MATU_C(i,NY)=0.
          VEC_C(i,NY)=TH_BC_YMAX_C1(N)
        end do
      else
! Neumann
! NOTE: BC enforced at GY(NY)

        do i=0,NKX
          MATL_C(i,NY)=-1.
          MATD_C(i,NY)=1.
          MATU_C(i,NY)=0.
          VEC_C(i,NY)=DY(NY)*TH_BC_YMAX_C1(N)
        end do
        do i=0,NKX
          MATL_C(i,NY+1)=0.
          MATD_C(i,NY+1)=1.
          MATU_C(i,NY+1)=0.
          VEC_C(i,NY+1)=0.
        end do


      end if
      return
      end

C----*|--.---------.---------.---------.---------.---------.---------.-|--
      SUBROUTINE APPLY_BC_VEL_LOWER
C----*|--.---------.---------.---------.---------.---------.---------.-|--
C This subroutine is called after initializing the flow
C It sets the appropriate boundary conditions including ghost cell values
C  on the velocity field in Fourier space
      INCLUDE 'header'
      INTEGER I,K      

C Now, apply the boundary conditions depending on the type specified 
      IF (U_BC_YMIN.EQ.0) THEN
C Dirichlet 
C Start with zero
         DO K=0,TNKZ
           DO I=0,NXP-1
             CU1(I,K,0)=0.d0
             CU1(I,K,1)=0.d0
           END DO
         END DO
C Now, set only the mean
         IF (RANKZ.EQ.0) THEN
           CU1(0,0,1)=U_BC_YMIN_C1
           CU1(0,0,0)=U_BC_YMIN_C1
         END IF
      ELSE IF (U_BC_YMIN.EQ.1) THEN
C Neumann
         DO K=0,TNKZ
           DO I=0,NXP-1
             CU1(I,K,1)=CU1(I,K,2)
             CU1(I,K,0)=CU1(I,K,1)
           END DO
         END DO
C Now, Apply BC to mean
         IF (RANKZ.EQ.0) THEN
           CU1(0,0,1)=CU1(0,0,2)-DY(2)*U_BC_YMIN_C1
           CU1(0,0,0)=CU1(0,0,1)-DY(1)*U_BC_YMIN_C1
         END IF
      ELSE
         STOP 'Error: U_BC_YMIN must be 0, or 1'
      END IF

      IF (W_BC_YMIN.EQ.0) THEN
C Dirichlet
C Start with zero
         DO K=0,TNKZ
           DO I=0,NXP-1
             CU3(I,K,0)=0.d0
             CU3(I,K,1)=0.d0
           END DO
         END DO
C Now, set only the mean
         IF (RANKZ.EQ.0) THEN
           CU3(0,0,1)=W_BC_YMIN_C1
           CU3(0,0,0)=W_BC_YMIN_C1
         END IF
      ELSE IF (W_BC_YMIN.EQ.1) THEN
C Neumann
         DO K=0,TNKZ
           DO I=0,NXP-1
             CU3(I,K,1)=CU3(I,K,2)
             CU3(I,K,0)=CU3(I,K,1)
           END DO
         END DO
C Now, Apply BC to mean
         IF (RANKZ.EQ.0) THEN
           CU3(0,0,1)=CU3(0,0,2)-DY(2)*W_BC_YMIN_C1
           CU3(0,0,0)=CU3(0,0,1)-DY(1)*W_BC_YMIN_C1
         END IF
      ELSE
         STOP 'Error: W_BC_YMIN must be 0, or 1' 
      END IF

      IF (V_BC_YMIN.EQ.0) THEN
C Dirichlet
C Set the vertical velocity at GYF(1) (halfway between GY(2) and GY(1))
         DO K=0,TNKZ
           DO I=0,NXP-1
             CU2(I,K,1)=2.d0*V_BC_YMIN_C1-CU2(I,K,2)  
             CU2(I,K,0)=CU2(I,K,1)  
           END DO
         END DO
      ELSE IF (V_BC_YMIN.EQ.1) THEN
C Neumann
         DO K=0,TNKZ
           DO I=0,NXP-1
             CU2(I,K,1)=CU2(I,K,2)
             CU2(I,K,0)=CU2(I,K,1)
           END DO
         END DO
         IF (RANKZ.EQ.0) THEN
           CU2(0,0,1)=CU2(0,0,2)-DYF(1)*V_BC_YMIN_C1 
           CU2(0,0,0)=CU2(0,0,1)-DYF(1)*V_BC_YMIN_C1 
         END IF

      ELSE IF (V_BC_YMIN.EQ.2) THEN
C Upstream-travelling wave proposed by Speyer/Kim
C (initialize as zero)
         IF (RANKZ.EQ.0) CU2(0,0,1)=-CU2(0,0,2)
      ELSE
         STOP 'Error: V_BC_YMIN must be 0, 1, or 2'
      END IF

      RETURN
      END

C----*|--.---------.---------.---------.---------.---------.---------.-|--
      SUBROUTINE APPLY_BC_VEL_UPPER
C----*|--.---------.---------.---------.---------.---------.---------.-|--
C This subroutine is called after initializing the flow
C It sets the appropriate boundary conditions including ghost cell values
C  on the velocity field in Fourier space
      INCLUDE 'header'
      INTEGER I,K      

! Now, apply boundary conditions to the top of the domain
      IF (U_BC_YMAX.EQ.0) THEN
C Dirichlet 
C Start with zero
         DO K=0,TNKZ
           DO I=0,NXP-1
             CU1(I,K,NY)=0.d0
             CU1(I,K,NY+1)=0.d0
           END DO
         END DO
C Now, set only the mean
         IF (RANKZ.EQ.0) THEN
           CU1(0,0,NY)=U_BC_YMAX_C1
           CU1(0,0,NY+1)=U_BC_YMAX_C1
         END IF
      ELSE IF (U_BC_YMAX.EQ.1) THEN
C Neumann
         DO K=0,TNKZ
           DO I=0,NXP-1
             CU1(I,K,NY)=CU1(I,K,NY-1)
             CU1(I,K,NY+1)=CU1(I,K,NY)
           END DO
         END DO
C Now, Apply BC to mean
         IF (RANKZ.EQ.0) THEN
           CU1(0,0,NY)=CU1(0,0,NY-1)+DY(NY)*U_BC_YMAX_C1
           CU1(0,0,NY+1)=CU1(0,0,NY)+DY(NY)*U_BC_YMAX_C1
         END IF
      ELSE
         STOP 'Error: U_BC_YMAX must be 0, or 1'
      END IF

      IF (W_BC_YMAX.EQ.0) THEN
C Dirichlet
C Start with zero
         DO K=0,TNKZ
           DO I=0,NXP-1
             CU3(I,K,NY)=0.d0
             CU3(I,K,NY+1)=0.d0
           END DO
         END DO
C Now, set only the mean
         IF (RANKZ.EQ.0) THEN
           CU3(0,0,NY)=W_BC_YMAX_C1
           CU3(0,0,NY+1)=W_BC_YMAX_C1
         END IF
C Ghost cell not used
         CU3(0,0,NY+1)=0.d0
      ELSE IF (W_BC_YMAX.EQ.1) THEN
C Neumann
         DO K=0,TNKZ
           DO I=0,NXP-1
             CU3(I,K,NY)=CU3(I,K,NY-1)
             CU3(I,K,NY+1)=CU3(I,K,NY)
           END DO
         END DO
C Now, Apply BC to mean
         IF (RANKZ.EQ.0) THEN
           CU3(0,0,NY)=CU3(0,0,NY-1)+DY(NY)*W_BC_YMAX_C1
           CU3(0,0,NY+1)=CU3(0,0,NY)+DY(NY)*W_BC_YMAX_C1
         END IF
      ELSE
        STOP 'Error: W_BC_YMAX must be 0, or 1'
      END IF

      IF (V_BC_YMAX.EQ.0) THEN
C Dirichlet
C Set the vertical velocity at GYF(NY) (halfway between GY(NY) and GY(NY+1))
         DO K=0,TNKZ
           DO I=0,NXP-1
             CU2(0,0,NY+1)=2.d0*V_BC_YMAX_C1-CU2(0,0,NY)
           END DO
         END DO
      ELSE IF (V_BC_YMAX.EQ.1) THEN
C Neumann
         DO K=0,TNKZ
           DO I=0,NXP-1
             CU2(I,K,NY)=CU2(I,K,NY-1)
             CU2(I,K,NY+1)=CU2(I,K,NY)
           END DO
         END DO
C Now, Apply BC to mean
         IF (RANKZ.EQ.0) THEN
           CU2(0,0,NY+1)=CU2(0,0,NY)+DY(NY)*V_BC_YMAX_C1
         END IF
      ELSE IF (V_BC_YMAX.EQ.2) THEN
C Upstream-travelling wave proposed by Speyer/Kim
C (initialize as zero gradient)
         IF (RANKZ.EQ.0) CU2(0,0,NY+1)=-CU2(0,0,NY)
      ELSE
         STOP 'Error: V_BC_YMAX must be 0, 1, or 2'
      END IF

      RETURN
      END


C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE THOMAS_REAL(A,B,C,G,NY,NX)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
C Uses the Thomas algorithm to solve Ax=b for tridiagonal A
C The RHS vector and solution are real
C Input lower, main, and upper diagonals, ld, md, ud, and rhs x
C Returns solution in x
C The indexing should be done by ROW, ie.
C [ b1  c1   0   0   0 ...
C [ a2  b2  c2   0   0 ...
C [  0  a3  b3   c3  0 ...

      INTEGER I, J, NX, NY
      REAL*8 A(0:NX,0:NY), B(0:NX,0:NY), C(0:NX,0:NY), G(0:NX,0:NY)

      DO J=0,NY-1
        DO I=0,NX
          A(I,J+1)=-A(I,J+1)/B(I,J)
          B(I,J+1)=B(I,J+1)+A(I,J+1)*C(I,J)
          G(I,J+1)=G(I,J+1)+A(I,J+1)*G(I,J)
        END DO
      END DO
      DO I=0,NX
        G(I,NY)=G(I,NY)/B(I,NY)
      END DO
      DO J=NY-1,0,-1
        DO I=0,NX
          G(I,J)=(G(I,J)-C(I,J)*G(I,J+1))/B(I,J)
        END DO
      END DO

      RETURN
      END

C----*|--.---------.---------.---------.---------.---------.---------.-|-------|    
      SUBROUTINE THOMAS_COMPLEX(A,B,C,G,NY,NX)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|

C Uses the Thomas algorithm to solve Ax=b for tridiagonal A
C The RHS vector and solution is complex
C Input lower, main, and upper diagonals, ld, md, ud, and rhs x
C Returns solution in x
C The indexing should be done by ROW, ie.
C [ b1  c1   0   0   0 ...
C [ a2  b2  c2   0   0 ...
C [  0  a3  b3   c3  0 ...

      INTEGER I, J, K, NY, NX
      REAL*8 A(0:NX,0:NY), B(0:NX,0:NY), C(0:NX,0:NY)
      COMPLEX*16 G(0:NX,0:NY)

      DO J=0,NY-1
        DO I=0,NX
          A(I,J+1)=-A(I,J+1)/B(I,J)
          B(I,J+1)=B(I,J+1)+A(I,J+1)*C(I,J)
          G(I,J+1)=G(I,J+1)+A(I,J+1)*G(I,J)
        END DO
      END DO
      DO I=0,NX
        G(I,NY)=G(I,NY)/B(I,NY)
      END DO
      DO I=0,NX
        DO J=NY-1,0,-1
          G(I,J)=(G(I,J)-C(I,J)*G(I,J+1))/B(I,J)
        END DO
      END DO

      RETURN
      END






