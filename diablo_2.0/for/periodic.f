C******************************************************************************|
C periodic.f, the fully-periodic-box solvers for diablo.           VERSION 0.9
C
C These solvers were written by Tom Bewley and John Taylor.
C******************************************************************************|

C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE INIT_PER
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|

      INCLUDE 'header'

C For forcing of homogeneous turbulence, set the target dissipation rate
      IF (N_TH.gt.0) THEN
      EPSILON_TARGET=((1.d0/DX(1))**4.d0)*(NU**3.d0)*(PR(1))**(-2.d0)
      write(*,*) 'EPSILON_TARGET: ',EPSILON_TARGET
      END IF
      
      RETURN
      END

C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE RK_PER_1
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
C Main time-stepping algorithm for the fully periodic case
C INPUTS  (in Fourier space):  CUi, P, and (if k>1) CFi at (k-1)  (for i=1,2,3)
C OUTPUTS (in Fourier space):  CUi, P, and (if k<3) CFi at (k)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      INCLUDE 'header'
      REAL*8 TEMP1, TEMP2, TEMP3, TEMP4, TEMP5, TEMP6
      INTEGER I,J,K,N

C Compute the RHS in Fourier Space, CRi.

C First, define the constants used for time-stepping
      TEMP1=NU*H_BAR(RK_STEP)/2.0
      TEMP2=BETA_BAR(RK_STEP)*H_BAR(RK_STEP)
      TEMP3=ZETA_BAR(RK_STEP)*H_BAR(RK_STEP)
      TEMP4=H_BAR(RK_STEP)
      TEMP5=0.0

      DO J=0,TNKY
        DO K=0,TNKZ
          DO I=0,NKX
C Start with the explicit part of the Crank-Nicolson viscous term and
C  the pressure gradient treated with Explicit Euler:
            TEMP5=1-TEMP1*(KX2(I)+KY2(J)+KZ2(K))
            CR1(I,K,J)=TEMP5*CU1(I,K,J)-TEMP4*(CIKX(I)*CP(I,K,J))
            CR2(I,K,J)=TEMP5*CU2(I,K,J)-TEMP4*(CIKY(J)*CP(I,K,J))
            CR3(I,K,J)=TEMP5*CU3(I,K,J)-TEMP4*(CIKZ(K)*CP(I,K,J))
C For each scalar, start with the explict part of the Crank-Nicolson
C diffusive term for each scalar
            DO N=1,N_TH
              TEMP6=1-(TEMP1/PR(N))*(KX2(I)+KY2(J)+KZ2(K))
              CRTH(I,K,J,N)=TEMP6*CTH(I,K,J,N)
            END DO
          END DO
        END DO
        IF (RK_STEP .GT. 1) THEN
          DO K=0,TNKZ
            DO I=0,NKX
C Add the term: ZETA_BAR(RK_STEP)*R(U(RK_STEP-1))
              CR1(I,K,J)=CR1(I,K,J)+TEMP3*CF1(I,K,J)
              CR2(I,K,J)=CR2(I,K,J)+TEMP3*CF2(I,K,J)
              CR3(I,K,J)=CR3(I,K,J)+TEMP3*CF3(I,K,J)
C Do the same for each scalar:
              DO N=1,N_TH
                CRTH(I,K,J,N)=CRTH(I,K,J,N)+TEMP3*CFTH(I,K,J,N)
              END DO
            END DO
          END DO
        END IF
      END DO
C If we are considering a linear background scalar gradient then add
C the term owing to advection of the background state.
C This allows us to consider a mean scalar gradient (ie stratification)
C even though the vertical boundary conditions are periodic.
C (In this case the passive scalar is a perturbation from a linear
C gradient. This gradient and the vertical domain size are used to
C make the passive scalar nondimensional, so here the nondimensional
C gradient is equal to one
      DO N=1,N_TH
      IF (BACKGROUND_GRAD(N)) THEN
C If there is a background scalar gradient add advection term:
      DO J=0,TNKY
        DO K=0,TNKZ
          DO I=0,NKX
            CFTH(I,K,J,N)=-CU2(I,K,J)
          END DO
        END DO
      END DO 
      ELSE
C Otherwise don't
      DO J=0,TNKY
        DO K=0,TNKZ
          DO I=0,NKX
            CFTH(I,K,J,N)=0.D0
          END DO
        END DO
      END DO 
      END IF
      END DO

C Inverse transform to physical space to compute the nonlinear terms
      CALL FFT_XZY_TO_PHYSICAL(CU1,U1)
      CALL FFT_XZY_TO_PHYSICAL(CU2,U2)
      CALL FFT_XZY_TO_PHYSICAL(CU3,U3)

      DO N=1,N_TH 
        CALL FFT_XZY_TO_PHYSICAL(CTH(0,0,0,N),TH(0,0,0,N))
      END DO

C Compute the nonlinear terms for the passive scalar equation
C Do this before the nonlinear momentum terms to use Fi as a working
C  array before using it for the momentum equation.

      DO N=1,N_TH
        DO J=0,NYM
          DO K=0,NZM
            DO I=0,NXM
              F1(I,K,J)=U1(I,K,J)*TH(I,K,J,N)
              F2(I,K,J)=U2(I,K,J)*TH(I,K,J,N)
              F3(I,K,J)=U3(I,K,J)*TH(I,K,J,N)
            END DO
          END DO
        END DO
        CALL FFT_XZY_TO_FOURIER(F1,CF1)
        CALL FFT_XZY_TO_FOURIER(F2,CF2)
        CALL FFT_XZY_TO_FOURIER(F3,CF3)
        DO J=0,TNKY
          DO K=0,TNKZ
            DO I=0,NKX
              CFTH(I,K,J,N)=CFTH(I,K,J,N)-CIKX(I)*CF1(I,K,J)
     &                      -CIKY(J)*CF2(I,K,J)
     &                      -CIKZ(K)*CF3(I,K,J)
     
            END DO
          END DO
        END DO
      END DO
C The RHS vector for the TH equation is now ready

C Compute the nonlinear terms for the momentum equations
      DO J=0,NYM
        DO K=0,NZM
          DO I=0,NXM
            F1(I,K,J)=U1(I,K,J)*U1(I,K,J)
            F2(I,K,J)=U1(I,K,J)*U2(I,K,J)
            F3(I,K,J)=U1(I,K,J)*U3(I,K,J)
            S1(I,K,J)=U2(I,K,J)*U2(I,K,J)
          END DO
        END DO
      END DO
      CALL FFT_XZY_TO_FOURIER(F1,CF1)
      CALL FFT_XZY_TO_FOURIER(F2,CF2)
      CALL FFT_XZY_TO_FOURIER(F3,CF3)
      CALL FFT_XZY_TO_FOURIER(S1,CS1)
C Here we start constructing the R-K terms in CFi
C Note, that the order of the following operations are important
      DO J=0,TNKY
        DO K=0,TNKZ
          DO I=0,NKX
            CF1(I,K,J)=-CIKX(I)*CF1(I,K,J)
     *                 -CIKY(J)*CF2(I,K,J)
     *                 -CIKZ(K)*CF3(I,K,J)
            CF2(I,K,J)=-CIKX(I)*CF2(I,K,J)
     *                 -CIKY(J)*CS1(I,K,J)
            CF3(I,K,J)=-CIKX(I)*CF3(I,K,J)
          END DO
        END DO
       END DO
C At this point, F1,F2,F3 contain some of the nonlinear terms

C Compute the remaining nonlinear terms
C We cannot use F1,F2,F3 for storage, so use S1 as the working variable
       DO J=0,NYM
        DO K=0,NZM
          DO I=0,NXM
            S1(I,K,J)=U2(I,K,J)*U3(I,K,J)
          END DO
        END DO 
       END DO
       CALL FFT_XZY_TO_FOURIER(S1,CS1)
       DO J=0,TNKY
         DO K=0,TNKZ
          DO I=0,NKX
            CF2(I,K,J)=CF2(I,K,J)-CIKZ(K)*CS1(I,K,J)
            CF3(I,K,J)=CF3(I,K,J)-CIKY(J)*CS1(I,K,J)
          END DO
         END DO
       END DO
       DO J=0,NYM
         DO K=0,NZM
          DO I=0,NXM
            S1(I,K,J)=U3(I,K,J)*U3(I,K,J)
          END DO
         END DO
       END DO
       CALL FFT_XZY_TO_FOURIER(S1,CS1)
       DO J=0,TNKY
         DO K=0,TNKZ
          DO I=0,NKX
            CF3(I,K,J)=CF3(I,K,J)-CIKZ(K)*CS1(I,K,J)
          END DO
         END DO
       END DO
C Done with the computation of nonlinear terms

C Optionally, add forcing terms to the RHS 
C Here we have U1, U2, U3, and TH in physical space
       CALL USER_RHS_PER_PHYSICAL

C If the scalar is active (RI NE 0), add the bouyancy forcing term
C  as explicit R-K
       DO N=1,N_TH
C Fisrt, convert back to Fourier space
       CALL FFT_XZY_TO_FOURIER(TH(0,0,0,N),CTH(0,0,0,N))
       DO J=0,TNKY
         DO K=0,TNKZ
           DO I=0,NKX
             CF2(I,K,J)=CF2(I,K,J)+RI(N)*CTH(I,K,J,N)
           END DO
         END DO
       END DO
       END DO

C Now, add the R-K terms to the RHS
      DO J=0,TNKY
        DO K=0,TNKZ
          DO I=0,NKX
            CR1(I,K,J)=CR1(I,K,J)+TEMP2*CF1(I,K,J)
            CR2(I,K,J)=CR2(I,K,J)+TEMP2*CF2(I,K,J)
            CR3(I,K,J)=CR3(I,K,J)+TEMP2*CF3(I,K,J)
          END DO
        END DO
      END DO

C Add R-K terms for the TH equation to the RHS
      DO N=1,N_TH
        DO J=0,TNKY
          DO K=0,TNKZ
            DO I=0,NKX
              CRTH(I,K,J,N)=CRTH(I,K,J,N)+TEMP2*CFTH(I,K,J,N)
            END DO
          END DO
        END DO
      END DO

C Computation of CRi complete.

C Now solve the implicit system for the intermediate field.
C (In the fully-periodic case, this is easy!)
      DO J=0,TNKY
        DO K=0,TNKZ
          DO I=0,NKX
            TEMP5=1+TEMP1*(KX2(I)+KY2(J)+KZ2(K))
            CU1(I,K,J)=CR1(I,K,J)/TEMP5
            CU2(I,K,J)=CR2(I,K,J)/TEMP5
            CU3(I,K,J)=CR3(I,K,J)/TEMP5
            DO N=1,N_TH
              TEMP6=1+(TEMP1/PR(N))*(KX2(I)+KY2(J)+KZ2(K)) 
              CTH(I,K,J,N)=CRTH(I,K,J,N)/TEMP6
            END DO
          END DO
        END DO
      END DO
C First step of the Fractional Step algorithm complete.


C Begin second step of the Fractional Step algorithm, making u divergence free.

C Compute varphi, store in the variable CR1, and project velocity field
      CALL REM_DIV_PER
C Then, update the pressure with phi.
      DO J=0,TNKY
        DO K=0,TNKZ
          DO I=0,NKX
            CP(I,K,J)=CP(I,K,J)+CR1(I,K,J)/TEMP4
          END DO
        END DO
      END DO
C Second step of the Fractional Step algorithm complete.

      RETURN
      END

C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE RK_PER_2
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
C Alternative time-stepping algorithm for the fully periodic case.
C INPUTS  (in Fourier space):  CUi, P, and (if k>1) CFi at (k-1)  (for i=1,2,3)
C OUTPUTS (in Fourier space):  CUi, P, and (if k<3) CFi at (k)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|

      RETURN
      END

C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE REM_DIV_PER
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      INCLUDE 'header'
      INTEGER I,J,K
      REAL*8  TEMP5

C Compute phi, store in the variable CR1.
C Note the coefficient H_BAR is absorbed into phi.
      DO J=0,TNKY
        DO K=0,TNKZ
          DO I=0,NKX
            TEMP5=-(KX2(I)+KY2(J)+KZ2(K)+EPS)
            CR1(I,K,J)=(CIKX(I)*CU1(I,K,J)+CIKY(J)*CU2(I,K,J)+
     *                  CIKZ(K)*CU3(I,K,J))/TEMP5
          END DO
        END DO
C Then update the CUi to make velocity field divergence-free.
        DO K=0,TNKZ
          DO I=0,NKX
            CU1(I,K,J)=CU1(I,K,J)-CIKX(I)*CR1(I,K,J)
            CU2(I,K,J)=CU2(I,K,J)-CIKY(J)*CR1(I,K,J)
            CU3(I,K,J)=CU3(I,K,J)-CIKZ(K)*CR1(I,K,J)
          END DO
        END DO
      END DO

      RETURN
      END

C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE POISSON_P_PER
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      INCLUDE 'header'
      INTEGER I,J,K

      DO J=0,TNKY
        DO K=0,TNKZ
          DO I=0,NKX
            CR1(I,K,J)=CIKX(I)*CU1(I,K,J)
            CR2(I,K,J)=CIKZ(K)*CU3(I,K,J)
          END DO
        END DO
      END DO

      CALL FFT_XZY_TO_PHYSICAL(CR1,R1)
      CALL FFT_XZY_TO_PHYSICAL(CR2,R2)

      DO J=0,NYM
        DO K=0,NZM
          DO I=0,NXM
            P(I,K,J)=R1(I,K,J)*R1(I,K,J)+R1(I,K,J)*R2(I,K,J)+
     *               R2(I,K,J)*R2(I,K,J)
          END DO
        END DO
      END DO
      DO J=0,TNKY
        DO K=0,TNKZ
          DO I=0,NKX
            CF1(I,K,J)=CIKY(J)*CU1(I,K,J)
            CF2(I,K,J)=CIKX(I)*CU2(I,K,J)
            CF3(I,K,J)=CIKZ(K)*CU1(I,K,J)
            CR1(I,K,J)=CIKX(I)*CU3(I,K,J)
            CR2(I,K,J)=CIKZ(K)*CU2(I,K,J)
            CR3(I,K,J)=CIKY(J)*CU3(I,K,J)
          END DO
        END DO
      END DO
      CALL FFT_XZY_TO_PHYSICAL(CF1,F1)
      CALL FFT_XZY_TO_PHYSICAL(CF2,F2)
      CALL FFT_XZY_TO_PHYSICAL(CF3,F3)
      CALL FFT_XZY_TO_PHYSICAL(CR1,R1)
      CALL FFT_XZY_TO_PHYSICAL(CR2,R2)
      CALL FFT_XZY_TO_PHYSICAL(CR3,R3)
      DO J=0,NYM
        DO K=0,NZM
          DO I=0,NXM
            P(I,K,J)=2*(P(I,K,J)+F1(I,K,J)*F2(I,K,J)
     *                          +F3(I,K,J)*R1(I,K,J)
     *                          +R2(I,K,J)*R3(I,K,J))
          END DO
        END DO
      END DO

      CALL FFT_XZY_TO_FOURIER(P,CP)
      DO J=0,TNKY
        DO K=0,TNKZ
          DO I=0,NKX
            CP(I,K,J)=CP(I,K,J)/(KX2(I)+KY2(J)+KZ2(K)+EPS)
          END DO
        END DO
      END DO
      CP(0,0,0)=0.0

      RETURN
      END


C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE CREATE_GRID_PER
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      INCLUDE 'header'
      INTEGER I,J,K

         IF (FLAVOR.NE.'Ensemble') WRITE (6,*) 'Fourier in X'
         DO I=0,NX
           GX(I)=(I*LX)/NX
           DX(I)=LX/NX
           IF (VERBOSITY .GT. 3) WRITE(6,*) 'GX(',I,') = ',GX(I)
         END DO
         IF (FLAVOR.NE.'Ensemble') WRITE (6,*) 'Fourier in Z'
         DO K=0,NZ
           GZ(K)=(K*LZ)/NZ
           DZ(K)=LZ/NZ
           IF (VERBOSITY .GT. 3) WRITE(6,*) 'GZ(',K,') = ',GZ(K)
         END DO
         IF (FLAVOR.NE.'Ensemble') WRITE (6,*) 'Fourier in Y'
         DO J=0,NY
           GY(J)=(J*LY)/NY
           DY(J)=LY/NY
           IF (VERBOSITY .GT. 3) WRITE(6,*) 'GY(',J,') = ',GY(J)
         END DO

         RETURN
         END


C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE INPUT_PER
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      INCLUDE 'header'

      INTEGER N
      REAL    VERSION, CURRENT_VERSION

! Read in input parameters specific for channel flow case
      OPEN (11,file='input_per.dat',form='formatted',status='old')
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
      DO N=1,N_TH
        READ(11,*)
        READ(11,*) BACKGROUND_GRAD(N)
      END DO

      RETURN
      END

C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE SAVE_STATS_PER(FINAL)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      INCLUDE 'header'
      CHARACTER*35 FNAME
      LOGICAL FINAL
          real*8 mean_u1(0:NYM)
          real*8 mean_u2(0:NYM)
          real*8 mean_u3(0:NYM)
          real*8 mean_th(0:NYM,1:N_TH)
          real*8 mean_p(0:NYM)


      integer i,j,k,n
      real*8 uc, ubulk
	
! Note that this routine uses CFi and CFTH for storage, so it should
! only be used between full R-K timesteps

      WRITE(6,*) 'Saving flow statistics.'

      if (FINAL) then
! We are done with the simulation
! Write out statistics to a file
        open(20,file='stats.txt',form='formatted',status='unknown')
        do j=0,NYM
          write(20,201) j,GYF(j),UBAR(j),VBAR(j),WBAR(j)
        end do
201     format(I3,',',F16.9,',',F16.9,',',F16.9,',',F16.9)
        do n=1,N_TH
        do j=0,NYM
          write(20,202) j,GYF(j),THBAR(j,n)
        end do
        end do
202     format(I3,',',F16.9,',',F16.9)

        else
! We are in the middle of a run, compile statistics

! Store the velocity in Fourier space in CRi
      do j=0,TNKY
        do k=0,TNKZ
          do i=0,NKX
            CR1(i,k,j)=CU1(i,k,j)
            CR2(i,k,j)=CU2(i,k,j)
            CR3(i,k,j)=CU3(i,k,j)
          end do
        end do
      end do

! Compute the vertical gradients in fourier space, store in CFi
      do j=0,TNKY
        do k=0,TNKZ
          do i=0,NKX
            CF1(i,k,j)=CIKY(j)*CU1(i,k,j)
            CF2(i,k,j)=CIKY(j)*CU2(i,k,j)
            CF3(i,k,j)=CIKY(j)*CU3(i,k,j)
          end do
        end do
      end do
! Save the  scalar in fourier space in CFTH
      do n=1,N_TH
        do j=0,TNKY
          do k=0,TNKZ
            do i=0,NKX
              CFTH(i,k,j,n)=CTH(i,k,j,n)
            end do
          end do
        end do
      end do

! Now, convert the velocity and vertical gradients to physical space
      CALL FFT_XZY_TO_PHYSICAL(CU1,U1)
      CALL FFT_XZY_TO_PHYSICAL(CU2,U2)
      CALL FFT_XZY_TO_PHYSICAL(CU3,U3)
      do n=1,N_TH
        CALL FFT_XZY_TO_PHYSICAL(CTH(0,0,0,n),TH(0,0,0,n))
      end do
      CALL FFT_XZY_TO_PHYSICAL(CF1,F1)
      CALL FFT_XZY_TO_PHYSICAL(CF2,F2)
      CALL FFT_XZY_TO_PHYSICAL(CF3,F3)

! First get the number of samples taken so far
      NSAMPLES=NSAMPLES+1
! Get the mean velocity

! First, get the mean in physical space
      do j=0,NYM
        mean_u1(j)=0.d0
        mean_u2(j)=0.d0
        mean_u3(j)=0.d0
        do n=1,N_TH
          mean_th(j,n)=0.d0
        end do
        mean_p(j)=0.d0
        do i=0,NXM
          do k=0,NZM
            mean_u1(j)=mean_u1(j)+U1(i,k,j)
            mean_u2(j)=mean_u2(j)+U2(i,k,j)
            mean_u3(j)=mean_u3(j)+U3(i,k,j)
            do n=1,N_TH
              mean_th(j,n)=mean_th(j,n)+TH(i,k,j,n)
            end do
            mean_p(j)=mean_p(j)+P(i,k,j)
          end do
        end do
        mean_u1(j)=mean_u1(j)/dble(NX*NZ)
        mean_u2(j)=mean_u2(j)/dble(NX*NZ)
        mean_u3(j)=mean_u3(j)/dble(NX*NZ)
         do n=1,N_TH
          mean_th(j,n)=mean_th(j,n)/dble(NX*NZ)
        end do
        mean_p(j)=mean_p(j)/dble(NX*NZ)
      end do
     
      do j=0,NYM
        UBAR(j)=(1./float(NSAMPLES))*mean_u1(j)
     &      +((float(NSAMPLES)-1.)/float(NSAMPLES))*UBAR(j)
        VBAR(j)=(1./float(NSAMPLES))*mean_u2(j)
     &      +((float(NSAMPLES)-1.)/float(NSAMPLES))*VBAR(j)
        WBAR(j)=(1./float(NSAMPLES))*mean_u3(J)
     &      +((float(NSAMPLES)-1.)/float(NSAMPLES))*WBAR(j)
        do n=1,N_TH
          THBAR(j,n)=(1./float(NSAMPLES))*mean_th(j,n)
     &      +((float(NSAMPLES)-1.)/float(NSAMPLES))*THBAR(j,n)
        end do
      end do

! Get the turbulent kinetic energy at each level 
      do j=0,NYM
        urms(j)=0.
        vrms(j)=0.
        wrms(j)=0.
      do k=0,NZM
      do i=0,NXM 
        urms(j)=urms(j)+(abs(U1(i,k,j)-mean_u1(j)))**2.
        vrms(j)=vrms(j)+(abs(U2(i,k,j)-mean_u2(j)))**2.
        wrms(j)=wrms(j)+(abs(U3(i,k,j)-mean_u3(j)))**2.
      end do
      end do
        urms(j)=sqrt(urms(j)/(float(NZ)*float(NX)))
        vrms(j)=sqrt(vrms(j)/(float(NZ)*float(NX)))
        wrms(j)=sqrt(wrms(j)/(float(NZ)*float(NX)))
      end do 
! Get the bulk rms value
      urms_b=0.
      do j=1,NYM
        urms_b=urms_b+0.5*(urms(j)+urms(j-1))*(GY(j)-GY(j-1))
      end do
      urms_b=urms_b/LY
! If we are in 2d:
      if (NY.eq.1) urms_b=urms(0)
       
! Compute the Reynolds stress and mean velocity gradient
      do j=0,NYM
        uv(j)=0. 
        uw(j)=0.
        wv(j)=0.
      do k=0,NZM
      do i=0,NXM
        uv(j)=uv(j)+(U1(i,k,j)-mean_u1(j))
     +    *(U2(i,k,j)-mean_u2(j))
        wv(j)=wv(j)+(U3(i,k,j)-mean_u3(j))
     +    *(U2(i,k,j)-mean_u2(j))
        uw(j)=uw(j)+(U1(i,k,j)-mean_u1(j))
     +    *(U3(i,k,j)-mean_u3(j))
      end do
      end do
        uv(j)=uv(j)/(float(NZ)*float(NX))
        uw(j)=uw(j)/(float(NZ)*float(NX))
        wv(j)=wv(j)/(float(NZ)*float(NX))
      end do
              
! Get the y-derivative of the mean velocity at GYF points
      do j=1,NY-2
        dudy(j)=(mean_u1(j+1)-mean_u1(j-1))/(GY(j+1)-GY(j-1))
        dwdy(j)=(mean_u3(j+1)-mean_u3(j-1))/(GY(j+1)-GY(j-1))
      end do
      j=0
        dudy(j)=(mean_u1(j+1)-mean_u1(NYM))/(2.d0*(GY(j+1)-GY(j)))
        dwdy(j)=(mean_u3(j+1)-mean_u3(NYM))/(2.d0*(GY(j+1)-GY(j)))
      j=NYM
        dudy(j)=(mean_u1(0)-mean_u1(j-1))/(2.d0*(GY(j)-GY(j-1)))
        dwdy(j)=(mean_u3(0)-mean_u3(j-1))/(2.d0*(GY(j)-GY(j-1)))

! Get the mean square shear
      do j=0,NYM
        shear(j)=0.d0
        do k=0,NZM
          do i=0,NXM
            shear(j)=shear(j)+F1(i,k,j)**2.d0+F3(i,k,j)**2.d0
          end do
        end do
        shear(j)=shear(j)/dble(NX*NZ)
      end do

! Write out the bulk rms velocity
      write(*,*) '<U_rms>: ',urms_b
      write(*,*) 'VERBOSITY: ',VERBOSITY

! Write out the mean statistics at each time
      open(40,file='mean.txt',form='formatted',status='unknown')
      write(40,*) TIME_STEP,TIME,DELTA_T,UBULK
      do j=0,NYM
        write(40,401) j,GY(J),mean_u1(j)
     +      ,mean_u2(j)
     +      ,mean_u3(j),urms(j),vrms(j),wrms(j)
     +      ,uv(j),uw(j),wv(j),dudy(j),dwdy(j),mean_p(j),shear(j)
      end do

401   format(I3,' ',14(F20.9,' '))



      if (MOVIE) then
! Output a 2d slice through the velocity field for animation in matlab
        open(79,file='movie_vel_xz.txt',status='unknown'
     &      ,form='formatted')
        do i=0,NXM
        do k=0,NZM
          write(79,*) U1(I,K,0)**2.+U3(I,K,0)**2.+U2(I,K,0)**2.
        end do
        end do

        open(80,file='movie_vel_xy.txt',status='unknown'
     &        ,form='formatted')
        do i=0,NXM
        do j=0,NYM
          write(80,*) U1(I,0,J)**2.+U3(I,0,J)**2.+U2(I,0,J)**2.
        end do
        end do

! This file will contain a single plane and is used in conjunction with
! the matlab script 'realtime_movie' to visualize data during
! simulation
        open (76,file='temp.txt',status='unknown',form='formatted')
        do K=0,NZM
          write(76,*) gz(k)
        end do
        do I=0,NXM
        do K=0,NZM
          write(76,*) U1(I,K,0)**2.+U3(I,K,0)**2.+U2(I,K,0)**2.
        end do
        end do
        close (76)
        CALL SYSTEM('mv temp.txt ../post_process/matlab/latest_slice.txt
     &')

      end if 

! Do over the number of passive scalars
      do n=1,N_TH

      do j=0,NYM
        thrms(j,n)=0.
      do k=0,NZM
      do i=0,NXM
        thrms(j,n)=thrms(j,n)+(abs(TH(i,k,j,n)-mean_th(j,n)))**2.
      end do
      end do
        thrms(j,n)=sqrt(thrms(j,n)/(float(NZ)*float(NX)))
      end do
! Compute the Reynolds stress and mean velocity gradient
      do j=0,NYM
        thv(j,n)=0.
      do k=0,NZM
      do i=0,NXM
       thv(j,n)=thv(j,n)+(TH(i,k,j,n)-mean_th(j,n))
     +    *(U2(i,k,n)-mean_u2(j))
      end do
      end do
      thv(j,n)=thv(j,n)/(float(NZ)*float(NX))
      end do

! Get the y-derivative of the mean velocity at GYF points
      do j=2,NY
        dthdy(j,n)=(CRTH(0,0,j+1,n)-CRTH(0,0,j-1,n))/(2.*DYF(j))
      end do
      do j=1,NY-2
        dthdy(j,n)=(mean_th(j+1,n)-mean_th(j-1,n))/(GY(j+1)-GY(j-1))
      end do
      j=0
       dthdy(j,n)=(mean_th(j+1,n)-mean_th(NYM,n))/(2.d0*(GY(j+1)-GY(j)))
      j=NYM
       dthdy(j,n)=(mean_th(0,n)-mean_th(j-1,n))/(2.d0*(GY(j)-GY(j-1)))


      if (MOVIE) then
! Output a 2d slice through the scalar field for animation in matlab
        if (n.eq.1) then
! Chose which scalar is to be outputted
        open(75,file='movie1_xy.txt',status='unknown',form='formatted')
        do i=0,NXM
        do j=0,NYM
          write(75,*) TH(I,0,J,n)
        end do
        end do
        open(77,file='movie1_xz.txt',status='unknown',form='formatted')
        do i=0,NXM
        do k=0,NZM
          write(77,*) TH(I,K,0,n)
        end do
        end do
        end if
        if (n.eq.2) then
! Chose which scalar is to be outputted
        open(85,file='movie2_xy.txt',status='unknown',form='formatted')
         
        do i=0,NXM
        do j=0,NYM
          write(85,*) TH(I,0,J,n)
        end do
        end do
        open(87,file='movie2_xz.txt',status='unknown',form='formatted')
        do i=0,NXM
        do k=0,NZM
          write(87,*) TH(I,K,0,n)
        end do
        end do
        end if
        if (n.eq.3) then
! Chose which scalar is to be outputted
        open(95,file='movie3_xy.txt',status='unknown',form='formatted')
        do i=0,NXM
        do j=0,NYM
          write(95,*) TH(I,0,J,n)
        end do
        end do
        open(97,file='movie3_xz.txt',status='unknown',form='formatted')
        do i=0,NXM
        do k=0,NZM
          write(97,*) TH(I,K,0,n)
        end do
        end do
        end if
      end if 


! End do over number of passive scalars, n
      end do

! Convert back to Fourier space
      do n=1,N_TH
        call fft_xzy_to_fourier(TH(0,0,0,n),CTH(0,0,0,n))
      end do

! Write out the mean statistics at each time
      open(41,file='mean_th.txt',form='formatted',status='unknown')
      write(41,*) TIME_STEP,TIME,DELTA_T,UBULK
      do n=1,N_TH 
      do j=0,NYM
        write(41,402) j,GYF(J),mean_th(j,n)
     &      ,dthdy(j,n),thrms(j,n),thv(j,n),pe_diss(j,n)
      end do
      end do

402   format(I3,' ',6(F30.25,' '))

      write(*,*) 'VERBOSITY: ',VERBOSITY
      if (VERBOSITY.gt.4) then 
      write(*,*) 'Outputting info for gnuplot...'
      open (unit=10, file="solution")
      do i=2,NXM
        do j=2,NYM
          write (10,*) i, j, U1(i,0,j)
        end do
        write (10,*) ""
      end do
      close (10)
      call system ('gnuplot <gnuplot.in') 
      end if

C Convert velocity back to Fourier space
      call fft_xzy_to_fourier(U1,CU1)
      call fft_xzy_to_fourier(U2,CU2)
      call fft_xzy_to_fourier(U3,CU3)

      end if

      call tkebudget_per
	
      RETURN
      END


      subroutine filter_per
C This subroutine applies a filter to the highest wavenumbers
C It should be applied to the scalars in Fourier space
C The filter used is a sharpened raised cosine filter

      include 'header'

      integer I,J,K,N

! Variables for horizontal filtering
      real*8 sigma0

C Set the filtering constants for the all directions
      DO N=1,N_TH
      DO i=0,NKX
       DO k=0,TNKZ
        DO j=0,TNKY
          sigma0=0.5d0*(1.d0+
     &       cos(sqrt((KX(i)*LX*1.d0/float(NX))**2.d0
     &            +(KZ(k)*LZ*1.d0/float(NZ))**2.d0
     &            +(KY(j)*LY*1.d0/float(NY))**2.d0)))
! Apply a sharpened raised cosine filter
          CTH(i,k,j,n)=CTH(i,k,j,n)*
     &          sigma0**4.d0*(35.d0-84.d0*sigma0
     &        +70.d0*sigma0**2.d0-20.d0*sigma0**3.d0)
        END DO
       END DO
      END DO
      END DO

       return
       end


      subroutine tkebudget_per
! NOte, it is important to only run this routine after complete R-K
!  time advancement since F1 is overwritten which is needed between R-K steps

      include 'header'

      integer i,j,k
      real epsilon_mean,k_eta

! Compute the turbulent dissipation rate, epsilon=nu*<du_i/dx_j du_i/dx_j>
      do j=0,NYM
        epsilon(j)=0.
      end do
! Store du/dx in CS1
      do j=0,TNKY
      do k=0,TNKZ
      do i=0,NKX
        CS1(i,k,j)=CIKX(i)*CU1(i,k,j)
      end do
      end do
      end do
! Convert to physical space
      call fft_xzy_to_physical(CS1,S1)
      do j=0,NYM
      do k=0,NZM
      do i=0,NXM
        epsilon(j)=epsilon(j)+(S1(i,k,j)**2.0)
      end do
      end do
      end do
! Store dv/dx in CS1
      do j=0,TNKY
      do k=0,TNKZ
      do i=0,NKX
        CS1(i,k,j)=CIKX(i)*CU2(i,k,j)
      end do
      end do
      end do
! Convert to physical space
      call fft_xzy_to_physical(CS1,S1)
      do j=0,NYM
      do k=0,NZM
      do i=0,NXM
        epsilon(j)=epsilon(j)+0.5*(S1(i,k,j)**2.0)
      end do
      end do
      end do
! Compute du/dy at GYF gridpoints, note remove mean
      do j=0,TNKY
      do k=0,TNKZ
      do i=0,NKX
        CF1(i,k,j)=CIKY(j)*CU1(i,k,j)
      end do
      end do
      end do
! Convert to physical space
      call fft_xzy_to_physical(CF1,F1)
      do j=0,NYM
      do k=0,NZM
      do i=0,NXM
        epsilon(j)=epsilon(j)+0.5*(F1(i,k,j)**2.0)
! Cross term dvdx*dudy
        epsilon(j)=epsilon(j)+(S1(i,k,j)*F1(i,k,j))
      end do
      end do
      end do
! Store dw/dx in CS1
      do j=0,TNKY
      do k=0,TNKZ
      do i=0,NKX
        CS1(i,k,j)=CIKX(i)*CU3(i,k,j)
      end do
      end do
      end do
! Convert to physical space
      call fft_xzy_to_physical(CS1,S1)
      do j=0,NYM
      do k=0,NZM
      do i=0,NXM
        epsilon(j)=epsilon(j)+0.5*(S1(i,k,j)**2.0)
      end do
      end do
      end do
! Compute du/dz at GYF gridpoints, note remove mean
! Store du/dz in CS1
      do j=0,TNKY
      do k=0,TNKZ
      do i=0,NKX
        CF1(i,k,j)=CIKZ(k)*CU1(i,k,j)
      end do
      end do
      end do
! Convert to physical space
      call fft_xzy_to_physical(CF1,F1)
      do j=0,NYM
      do k=0,NZM
      do i=0,NXM
        epsilon(j)=epsilon(j)+0.5*(F1(i,k,j)**2.0)
! Cross term dudz*dwdx
        epsilon(j)=epsilon(j)+S1(i,k,j)*F1(i,k,j)
      end do
      end do
      end do
! Compute dv/dy at GYF gridpoints, note remove mean
      do j=0,TNKY
      do k=0,TNKZ
      do i=0,NKX
        CS1(i,k,j)=CIKY(j)*CU2(i,k,j)
      end do
      end do
      end do
! Convert to physical space
      call fft_xzy_to_physical(CS1,S1)
      do j=0,NYM
      do k=0,NZM
      do i=0,NXM
        epsilon(j)=epsilon(j)+(S1(i,k,j)**2.0)
      end do
      end do
      end do
! Compute dw/dy at GYF gridpoints, note remove mean
      do j=0,TNKY
      do k=0,TNKZ
      do i=0,NKX
        CS1(i,k,j)=CIKY(j)*CU3(i,k,j)
      end do
      end do
      end do
! Convert to physical space 
      call fft_xzy_to_physical(CS1,S1)
      do j=0,NYM
      do k=0,NZM
      do i=0,NXM
        epsilon(j)=epsilon(j)+0.5*(S1(i,k,j)**2.0)
      end do
      end do
      end do
! Store dv/dz in CF1
      do j=0,TNKY
      do k=0,TNKZ
      do i=0,NKX
        CF1(i,k,j)=CIKZ(k)*CU2(i,k,j)
      end do
      end do
      end do
! Convert to physical space
      call fft_xzy_to_physical(CF1,F1)
      do j=0,NYM
      do k=0,NZM
      do i=0,NXM
        epsilon(j)=epsilon(j)+0.5*(F1(i,k,j)**2.0)
! Cross term dvdz*dwdy
        epsilon(j)=epsilon(j)+S1(i,k,j)*F1(i,k,j)
      end do
      end do
      end do
! Store dw/dz in CS1
      do j=0,TNKY
      do k=0,TNKZ
      do i=0,NKX
        CS1(i,k,j)=CIKZ(k)*CU3(i,k,j)
      end do
      end do
      end do
! Convert to physical space
      call fft_xzy_to_physical(CS1,S1)
      do j=0,NYM
      do k=0,NZM
      do i=0,NXM
        epsilon(j)=epsilon(j)+(S1(i,k,j)**2.0)
      end do
      end do
      end do
      do j=0,NYM
        epsilon(j)=epsilon(j)/float(NX*NZ)
      end do


! Write out the bulk rms velocity
      write(*,*) '<U_rms>: ',urms_b


! Write out the mean statistics at each time
      open(45,file='tke.txt',form='formatted',status='unknown')
      write(45,*) TIME_STEP,TIME,DELTA_T
      do j=0,NYM
        write(45,401) j,GY(J),epsilon(j)
      end do
401   format(I3,' ',2(F20.9,' '))


! Get Kolmogorov wavelength
      epsilon_mean=NU*SUM(epsilon(0:NYM))/dble(NY) 

      k_eta=2.d0*PI*(NU**3.d0/epsilon_mean)**(-0.25d0)

      write(*,*) 'Kolmogorov scale: ',2.d0*PI/k_eta

      return
      end






