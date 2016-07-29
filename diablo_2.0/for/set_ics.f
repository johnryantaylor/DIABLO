C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE CREATE_FLOW_CHAN
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      INCLUDE 'header'


      INTEGER I,J,K,N
      REAL*8 RNUM1,RNUM2,RNUM3
      INTEGER, DIMENSION(:), ALLOCATABLE :: seed

C Initialize the random number generator
      CALL RANDOM_SEED(SIZE = K)
      Allocate (seed(1:K))
      do k=1,K
        seed(k)=RANK+k+999
      end do
      CALL RANDOM_SEED(PUT = seed)

C UBULK0 and KICK should be set in input.dat

C IC_TYPE is set in input_chan.dat and can be used to easily
C control which initial condition is used.  A few examples
C are given here. These can be modified, or new types can be
C added 

       IF (IC_TYPE.eq.0) then
C Parabolic profile for laminar closed channel flow
       DO J=0,NY
         DO K=0,NZP-1
           DO I=0,NXM
             U1(I,K,J)=(3./2.)*UBULK0*(1.d0-GYF(J)**2.)
             U2(I,K,J)=0.
             U3(I,K,J)=0.
           END DO
         END DO
      END DO
      else if (IC_TYPE.eq.1) then
C Laminar profile for open channel flow :
       DO K=0,NZP-1
         DO I=0,NXM
           DO J=1,NY
             U1(I,K,J)=-(3./2.)*UBULK0*GYF(J)**2.+3.*UBULK0*GYF(J)
             U2(I,K,J)=0.
             U3(I,K,J)=0.
           END DO
           U1(I,K,0)=0.
           U3(I,K,0)=0.
           U1(I,K,NY+1)=0.
           U3(I,K,NY+1)=0.
         END DO
      END DO
      else if (IC_TYPE.eq.2) then
C Linear profile for laminar Couette flow:
       DO J=0,NY
         DO K=0,NZP-1
           DO I=0,NXM
             U1(I,K,J)=gyf(j)
             U2(I,K,J)=0.
             U3(I,K,J)=0.
           END DO
         END DO
      END DO
      else if (IC_TYPE.eq.3) then
C Tanh shear layer
       DO J=0,NY
         DO K=0,NZP-1
           DO I=0,NXM
             U1(I,K,J)=TANH(GYF(J))
             U2(I,K,J)=0.d0
             U3(I,K,J)=0.d0
           END DO
         END DO
       END DO
      else if (IC_TYPE.eq.4) then
C For Front
C Initialize in thermal wind balance: 
       DO J=0,NY
         DO K=0,NZP-1
           DO I=0,NXM
             U1(I,K,J)=0.d0
             U2(I,K,J)=0.d0
             U3(I,K,J)=0.d0
           END DO
         END DO
       END DO
      else
        WRITE(*,*) 'WARNING, unsupported IC_TYPE in CREATE_FLOW'
      end if
C Add random noise in physical space
      CALL RANDOM_NUMBER(RNUM1)
      CALL RANDOM_NUMBER(RNUM1)
      CALL RANDOM_NUMBER(RNUM1)

      DO J=0,NY+1
       DO K=0,NZP-1
         DO I=0,NXM
           CALL RANDOM_NUMBER(RNUM1)
           U1(I,K,J)=U1(I,K,J)+KICK*(RNUM1-0.5d0)
           CALL RANDOM_NUMBER(RNUM1)
           U2(I,K,J)=U2(I,K,J)+KICK*(RNUM1-0.5d0)
           CALL RANDOM_NUMBER(RNUM1)
           U3(I,K,J)=U3(I,K,J)+KICK*(RNUM1-0.5d0)
         END DO
        END DO
      END DO

C Zero the ghost cells
       IF (.NOT.USE_MPI) THEN
       DO K=0,NZM
         DO I=0,NXM
           U1(I,K,0)=0.
           U2(I,K,0)=0.
           U3(I,K,0)=0.
           U1(I,K,NY+1)=0.
           U2(I,K,NY+1)=0.
           U3(I,K,NY+1)=0.
         END DO
      END DO
      END IF

C Convert to Fourier space      
      CALL FFT_XZ_TO_FOURIER(U1,CU1,0,NY+1)
      CALL FFT_XZ_TO_FOURIER(U2,CU2,0,NY+1)
      CALL FFT_XZ_TO_FOURIER(U3,CU3,0,NY+1)
      CALL FFT_XZ_TO_FOURIER(P,CP,0,NY+1)

! Optionally, add random noise in Fourier space instead
!      DO I=0,NXP-1
!        DO J=1,NY
!          DO K=0,TNKZ
C Now, give the velocity field a random perturbation
!            CALL RANDOM_NUMBER(RNUM1)
!            CALL RANDOM_NUMBER(RNUM2)
!            CU1(I,K,J)=CU1(I,K,J)
!     &           +CMPLX((RNUM1-0.5d0),(RNUM2-0.5d0))*KICK
!            CALL RANDOM_NUMBER(RNUM1)
!            CALL RANDOM_NUMBER(RNUM2)
!            CU2(I,K,J)=CU2(I,K,J)
!     &           +CMPLX((RNUM1-0.5d0),(RNUM2-0.5d0))*KICK
!            CALL RANDOM_NUMBER(RNUM1)
!            CALL RANDOM_NUMBER(RNUM2)
!            CU3(I,K,J)=CU3(I,K,J)
!     &           +CMPLX((RNUM1-0.5d0),(RNUM2-0.5d0))*KICK
!          END DO
!          IF (TNKZ.EQ.0) THEN
! Here, In the 2d case we want to add a kick to the mean in z
!            K=0         
!            CALL RANDOM_NUMBER(RNUM1)
!            CALL RANDOM_NUMBER(RNUM2)
!            CALL RANDOM_NUMBER(RNUM3)

!            IF (IC_TYPE.eq.3) THEN
!              CU1(I,K,J)=CU1(I,K,J)
!     &             +(RNUM1-0.5)*KICK*EXP(-(GYF(J)*20.d0)**2.d0)
!              CU2(I,K,J)=CU2(I,K,J)
!     &             +(RNUM1-0.5)*KICK*EXP(-(GYF(J)*20.d0)**2.d0)
!              CU3(I,K,J)=CU3(I,K,J)
!     &             +(RNUM1-0.5)*KICK*EXP(-(GYF(J)*20.d0)**2.d0)
!            ELSE
!              CU1(I,K,J)=CU1(I,K,J)+(RNUM1-0.5)*KICK
!              CU2(I,K,J)=CU2(I,K,J)+(RNUM2-0.5)*KICK
!              CU3(I,K,J)=CU3(I,K,J)+(RNUM3-0.5)*KICK
!            END IF
!          END IF 
!        END DO
!      END DO

      IF (USE_MPI) THEN
        CALL GHOST_CHAN_MPI
      END IF
C Apply Boundary conditions to velocity field
      IF (USE_MPI) THEN
        CALL APPLY_BC_VEL_MPI
      ELSE
        CALL APPLY_BC_VEL_LOWER
        CALL APPLY_BC_VEL_UPPER
      END IF

C Remove the divergence of the velocity field
      CALL REM_DIV_CHAN

      IF (USE_MPI) THEN
        CALL GHOST_CHAN_MPI
      END IF

C Get the pressure from the poisson equation
!      CALL POISSON_P_CHAN
! Fix for the pressure
!      IF (USE_MPI) THEN
!        CALL GHOST_CHAN_MPI
!      END IF

C Save various statistics to keep track of the initial condition
      CALL SAVE_STATS_CHAN(.FALSE.)

      RETURN
      END


C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE CREATE_TH_CHAN
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
C Initialize the scalar fields
C In this subroutine, you should initialize each scalar field for the
C particular problem of interest

      INCLUDE 'header'
      INTEGER I,J,K,N
! A variable for Front case...
      REAL*8 RI_B(0:NY+1)

      DO N=1,N_TH
        IF (CREATE_NEW_TH(N)) THEN

      IF (IC_TYPE.eq.0) THEN
! As an example, initialize TH1 with a sine in x
       DO K=0,NZP-1
         DO I=0,NXM
           DO J=1,NY
             TH(I,K,J,N)=sin(2.d0*PI*GX(I)/LX)/(4.d0*PI**2.d0)
           END DO
         END DO
       END DO
       ELSE IF ((IC_TYPE.eq.1).or.(IC_TYPE.eq.2)) THEN
! Initialize with a linear profile using the bcs
       DO K=0,NZP-1
         DO I=0,NXM
           IF ((TH_BC_YMIN(N).EQ.0).AND.(TH_BC_YMAX(N).EQ.0)) THEN
               DO J=1,NY
               IF (GYF(J).LE.2.0) THEN
                 TH(I,K,J,N)=(TH_BC_YMAX_C1(N)-TH_BC_YMIN_C1(N))
     &                *(GYF(J)+1.)/2.0+TH_BC_YMIN_C1(N)
               ELSE
                 TH(I,K,J,N)=TH_BC_YMAX_C1(N)
               END IF
             END DO
           ELSE IF ((TH_BC_YMIN(N).EQ.1)
     &            .AND.(TH_BC_YMAX(N).EQ.1)) THEN
             DO J=1,NY
! Linear profile with slope corresponding to upper value
                TH(I,K,J,N)=TH_BC_YMAX_C1(N)*GYF(J)
              END DO
           ELSE
             IF (RANK.EQ.0) then
                WRITE(*,*) 'WARNING, THETA INITIALIZED TO ZERO ...'
                WRITE(*,*) 'CREATE AN INITIAL VALUE IN CREATE_FLOW_CHAN'
             end if
           END IF
         END DO
        END DO
        ELSE IF (IC_TYPE.eq.3) then
! Tanh profile 
       DO K=0,NZP-1
         DO I=0,NXM
           DO J=1,NY
             TH(I,K,J,N)=TANH(GYF(J))
           END DO
         END DO
       END DO
       ELSE IF (IC_TYPE.eq.4) THEN
! For Front case, specify given RI_B profile
       DO K=0,NZP-1
         DO I=0,NXM
          TH(I,K,0,N)=0.d0
           DO J=1,NY
             if (GYF(J).lt.-60.d0) then
               RI_B(J)=20.d0
               TH(I,K,J,N)=(GYF(J)-GYF(1))*
     &                    RI_B(J)*(RI(N)*DRHODX(N))**2.d0
     &                    /I_RO**2.d0/RI(N)
             else
                   RI_B(J)=1.0d0
                   TH(I,K,J,N)=(GYF(J)+60.d0)*
     &                    RI_B(J)*(RI(N)*DRHODX(N))**2.d0
     &                    /I_RO**2.d0/RI(N)
     &                   +(-60+140.d0)*20.d0*(RI(N)*DRHODX(N))**2.d0
     &                    /I_RO**2.d0/RI(N)
             end if
          END DO 
        END DO
      END DO

       ELSE
        WRITE(*,*) 'WARNING, unsupported IC_TYPE in CREATE_FLOW'
        END IF


      S1(:,:,:)=TH(:,:,:,N)
      CALL FFT_XZ_TO_FOURIER(S1,CS1,0,NY+1)
      CTH(:,:,:,N)=CS1(:,:,:)

      END IF
      END DO

      RETURN
      END


C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE CREATE_FLOW_PER
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      INCLUDE 'header'
      INTEGER I, J, K
      REAL*8 RNUM1,RNUM2,RNUM3
      REAL*8 K0

C For an initial vortex, define the location of the centerline
      REAL*8 XC(0:NY+1),ZC(0:NY+1)

      WRITE(6,*) 'Creating new flow from scratch.'

C Initialize random number generator
      CALL RANDOM_SEED

      IF (IC_TYPE.eq.0) THEN
C Initizlize the flow using a Taylor-Green vortex
C Nondimensionalize with U0 and 1/kappa
      DO J=0,NYM
        DO K=0,NZM
          DO I=0,NXM
            U1(I,K,J)=cos(2*pi*(GY(J))/LY)
     &               *cos(2*pi*(GX(I))/LX)
     &               *SIN(2*pi*(GZ(K))/LZ)
            U2(I,K,J)=0.d0
            U3(I,K,J)=-cos(2*pi*(GY(J))/LY)
     &               *sin(2*pi*(GX(I))/LX)
     &               *COS(2*pi*(GZ(K))/LZ)
          END DO
        END DO
      END DO
      ELSE IF (IC_TYPE.eq.1) THEN
C Start with an ideal vortex centered in the domain
      DO J=0,NYM
        XC(J)=LX/2.
        ZC(J)=LZ/2.
        DO K=0,NZM
          DO I=0,NXM
            IF ((GX(I)-XC(j))**2.+(GZ(K)-ZC(j))**2..gt.0.1) then
! If we aren't too close to the vortex center
              U1(I,K,J)=-1.d0*(GZ(K)-ZC(j))
     &                /((GX(I)-XC(j))**2.+(GZ(K)-ZC(j))**2.)
              U3(I,K,J)=1.d0*(GX(I)-XC(j))
     &                /((GX(I)-XC(j))**2.+(GZ(K)-ZC(j))**2.)
              U2(I,K,J)=0.d0
            ELSE
! Otherwise:
              U1(I,K,J)=-1.d0*(GZ(K)-ZC(j))
     &                /0.1
              U3(I,K,J)=1.d0*(GX(I)-XC(j))
     &                /0.1
              U2(I,K,J)=0.d0
            END IF
          END DO
        END DO
      END DO
      END IF
! Add random noise in Fourier space

      CALL FFT_XZY_TO_FOURIER(U1,CU1)
      CALL FFT_XZY_TO_FOURIER(U2,CU2)
      CALL FFT_XZY_TO_FOURIER(U3,CU3)
      DO J=1,TNKY
        DO K=1,TNKZ
          DO I=1,NKX
            CALL RANDOM_NUMBER(RNUM1)
            CALL RANDOM_NUMBER(RNUM2)
            CALL RANDOM_NUMBER(RNUM3)
            K0=sqrt(KX(I)**2.d0+KY(J)**2.d0+KZ(K)**2.d0)
     &        /sqrt(KX(1)**2.d0+KY(1)**2.d0+KZ(1)**2.d0)
            CU1(I,K,J)=CU1(I,K,J)+(RNUM1-0.5)*KICK/K0
            CU2(I,K,J)=CU2(I,K,J)+(RNUM1-0.5)*KICK/K0
            CU3(I,K,J)=CU3(I,K,J)+(RNUM1-0.5)*KICK/K0
          end do
        end do
      end do
! get the initial energy in low wavenumbers
      CALL FFT_XZY_TO_PHYSICAL(CU1,U1)
      CALL FFT_XZY_TO_PHYSICAL(CU2,U2)
      CALL FFT_XZY_TO_PHYSICAL(CU3,U3)
      EK0=0.d0
      DO J=0,NYM
        DO K=0,NZM
          DO I=0,NXM
              EK0=EK0+U1(I,K,J)**2.d0+U2(I,K,J)**2.d0+U3(I,K,J)**2.d0
          END DO
        END DO
      END DO
      write(*,*) 'EK0: ',EK0
      IF (N_TH.gt.0) THEN
!      EPSILON_TARGET=((1.d0/DX(1))**4.d0)*(NU**3.d0)*(PR(1))**(-2.d0)
      EPSILON_TARGET=((1.d0/DX(1))**4.d0)*(NU**3.d0)*(100.d0)**(-2.d0)
      write(*,*) 'EPSILON_TARGET: ',EPSILON_TARGET
      write(*,*) 'TARGET KOLMOGOROV SCALE: ',
     &         (NU**3.d0/epsilon_target)**(0.25d0)
      END IF
      CALL FFT_XZY_TO_FOURIER(U1,CU1)
      CALL FFT_XZY_TO_FOURIER(U2,CU2)
      CALL FFT_XZY_TO_FOURIER(U3,CU3)


      CALL REM_DIV_PER
      CALL POISSON_P_PER

      CALL SAVE_STATS_PER(.FALSE.)

      RETURN
      END
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE CREATE_TH_PER
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      INCLUDE 'header'
      INTEGER I,J,K,N

C Note, Since stratification is not permitted in the periodic flow field
C Any background stratification must be added to the governing equations

      DO N=1,N_TH
      IF (CREATE_NEW_TH(N)) THEN
        DO J=0,NYM
          DO K=0,NZM
            DO I=0,NXM
C Example: Gaussian patch centered in the domain
         TH(I,K,J,N)=EXP(-((GX(I)-LX/2)*10.d0)**2.d0
     &                   -((GY(J)-LY/2)*10.d0)**2.d0
     &                   -((GZ(K)-LZ/2)*10.d0)**2.d0)
            END DO
          END DO
        END DO
       CALL FFT_XZY_TO_FOURIER(TH(0,0,0,N),CTH(0,0,0,N))

       END IF
       END DO

       RETURN
       END

C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE CREATE_FLOW_DUCT
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|

      RETURN
      END


C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE CREATE_FLOW_CAV
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|

      RETURN
      END





