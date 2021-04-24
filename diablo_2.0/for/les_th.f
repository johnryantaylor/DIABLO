      subroutine les_chan_th
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
      include 'header_les_th'

      integer i,j,k,l,m,ij,n

      real*8 S1_mean(0:NY+1)
      real*8 NU_T_mean(0:NY+1)
      real*8 KAPPA_T_mean(0:NY+1)
 
      real*8 C_SMAG
!      parameter (C_SMAG=0.13d0)
      parameter (C_SMAG=0.0433d0)
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

! First, for all models, apply boundary conditions to the velocity field
! (fill ghost cells) to ensure accurate calculation of gradients
C Apply Boundary conditions to velocity field
      IF (USE_MPI) THEN
        CALL APPLY_BC_VEL_MPI
        CALL APPLY_BC_SCALAR_MPI
        CALL GHOST_CHAN_MPI
      ELSE
        CALL APPLY_BC_VEL_LOWER
        CALL APPLY_BC_VEL_UPPER
        CALL APPLY_BC_SCALAR_LOWER
        CALL APPLY_BC_SCALAR_UPPER
      END IF

! If we are using Neuman boundary conditions, over-write the values of the
! velocity at the ghost cells so that the LES model doesn't use the large
! velocity gradient
      CALL APPLY_BC_LES
      CALL APPLY_BC_SCALAR_LES

      if (LES_MODEL_TYPE.EQ.1) then
!     Constant Smagorinsky model

! APPLY constant SGS Prandlt number
         DO N=1,N_TH
         do j=1,NY+1
           do k=0,NZP-1
             do i=0,NXM
               KAPPA_T(I,K,J,N)=1.d0*NU_T(I,K,J)
             end do
           end do
         end do
         end do

! Convert the velocity to physical space
      call FFT_XZ_TO_PHYSICAL(CU1,U1,0,NY+1)
      call FFT_XZ_TO_PHYSICAL(CU2,U2,0,NY+1)
      call FFT_XZ_TO_PHYSICAL(CU3,U3,0,NY+1)

! Convert the scalar to physical space

        DO N=1,N_TH
          CS1(:,:,:)=CTH(:,:,:,N)
          CALL FFT_XZ_TO_PHYSICAL(CS1,S1,0,NY+1)
          TH(:,:,:,N)=S1(:,:,:)
        END DO

      else if ((LES_MODEL_TYPE.EQ.2).or.(LES_MODEL_TYPE.eq.3)) then
! Here, use a dynamic smagorinsky model with or without scale similar part

      stop 'ERROR: LES_MODEL_TYPE=2, 3 not supported yet in MPI'  

      else if (LES_MODEL_TYPE.EQ.4) then
!     Anisotrophic minimum-dissipation model Rozema

! APPLY constant SGS Prandlt number
         DO N=1,N_TH
         do j=1,NY+1
           do k=0,NZP-1
             do i=0,NXM
               KAPPA_T(I,K,J,N)=1.d0*NU_T(I,K,J)
             end do
           end do
         end do
         end do


      else if (LES_MODEL_TYPE.EQ.5) then
!     Anisotrophic minimum-dissipation model Verstappen

! Compute all the velocity and scalar  gradients 

      call compute_all_gradients_chan

! Convert the velocity to physical space
      call FFT_XZ_TO_PHYSICAL(CU1,U1,0,NY+1)
      call FFT_XZ_TO_PHYSICAL(CU2,U2,0,NY+1)
      call FFT_XZ_TO_PHYSICAL(CU3,U3,0,NY+1)

! Convert the scalar to physical space

        DO N=1,N_TH
          CS1(:,:,:)=CTH(:,:,:,N)
          CALL FFT_XZ_TO_PHYSICAL(CS1,S1,0,NY+1)
          TH(:,:,:,N)=S1(:,:,:)
        END DO

         deltax = (DX(1)*beta_sgs)
         deltaz = (DZ(1)*beta_sgs)

         DO N=1,N_TH

! First compute at GYF points 
!      DO J=JSTART,JEND
      DO J=JSTART_TH(N),JEND_TH(N)
!Set filter length (based on grid size) in y direction
!Based on constant Smag code        
         deltay=(DYF(J)*2.d0)
       DO K=0,NZP-1
        DO I=0,NXM
!First calculate numerator and store it in S1_th.
            S1_th(I,K,J)=-(
     &        deltax**2*du1dx(I,K,J)*dthetadx(I,K,J,N)**2
     &       +deltay**2*du2dy(I,K,J)
     &         *(0.5d0*(dthetady(I,K,J+1,N)+dthetady(I,K,J,N)))**2
     &       +deltaz**2*du3dz(I,K,J)*dthetadz(I,K,J,N)**2
     &       +(deltay**2*0.5d0*(du1dy(I,K,J+1)+du1dy(I,K,J))
     &         +deltax**2*0.5d0*(du2dx(I,K,J+1)+du2dx(I,K,J)))
     &         *dthetadx(I,K,J,N)
     &         *0.5d0*(dthetady(I,K,J+1,N)+dthetady(I,K,J,N))
     &       +(deltaz**2*du1dz(I,K,J)+deltax**2*du3dx(I,K,J))
     &         *dthetadx(I,K,J,N)
     &         *dthetadz(I,K,J,N)
     &       +(deltaz**2*0.5d0*(du2dz(I,K,J+1)+du2dz(I,K,J))
     &         +deltay**2*0.5d0*(du3dy(I,K,J+1)+du3dy(I,K,J)))
     &         *0.5d0*(dthetady(I,K,J+1,N)+dthetady(I,K,J,N))
     &         *dthetadz(I,K,J,N)
     &        )

        IF (S1_th(I,K,J) <= 0.0d0) THEN
          S1_th(I,K,J)=0.0d0
        ELSE

        S1_th(I,K,J)=S1_th(I,K,J)/
     &  ((deltax*dthetadx(I,K,J,N))**2
     &  +(deltay*0.5d0*(dthetady(I,K,J+1,N)+dthetady(I,K,J,N)))**2
     &  +(deltaz*dthetadz(I,K,J,N))**2)

        END IF


          END DO
        END DO
      END DO

! Compute kappa_e at GY points and store in TEMP_th
      DO J=2,NY+1
!Set filter length (based on grid size) in y direction
!Based on constant Smag code        
         deltay=(DY(J)*2.d0)
       DO K=0,NZP-1
        DO I=0,NXM
!First calculate numerator and store it in TEMP_th.
            TEMP_th(I,K,J)=-(
     &        deltax**2*(du1dx(I,K,J)*DYF(j-1)+du1dx(I,K,J-1)*DYF(j))
     &                               /(2.d0*DY(j))
     &         *((dthetadx(I,K,J,N)*DYF(j-1)+dthetadx(I,K,J-1,N)*DYF(j))
     &                               /(2.d0*DY(j)))**2
     &       +deltay**2*(du2dy(I,K,J)*DYF(j-1)+du2dy(I,K,J-1)*DYF(j))
     &                               /(2.d0*DY(j))
     &         *(dthetady(I,K,J,N))**2
     &       +deltaz**2*(du3dz(I,K,J)*DYF(j-1)+du3dz(I,K,J-1)*DYF(j))
     &                               /(2.d0*DY(j))
     &         *((dthetadz(I,K,J,N)*DYF(j-1)+dthetadz(I,K,J-1,N)*DYF(j))
     &                               /(2.d0*DY(j)))**2
     &       +(deltay**2*du1dy(I,K,J)
     &         +deltax**2*du2dx(I,K,J))
     &         *(dthetadx(I,K,J,N)*DYF(j-1)+dthetadx(I,K,J-1,N)*DYF(j))
     &                               /(2.d0*DY(j))
     &         *(dthetady(I,K,J,N))
     &       +(deltaz**2*(du1dz(I,K,J)*DYF(j-1)+du1dz(I,K,J-1)*DYF(j))
     &                               /(2.d0*DY(j))
     &         +deltax**2*(du3dx(I,K,J)*DYF(j-1)+du3dx(I,K,J-1)*DYF(j))
     &                               /(2.d0*DY(j)))
     &         *(dthetadx(I,K,J,N)*DYF(j-1)+dthetadx(I,K,J-1,N)*DYF(j))
     &                               /(2.d0*DY(j))
     &         *(dthetadz(I,K,J,N)*DYF(j-1)+dthetadz(I,K,J-1,N)*DYF(j))
     &                               /(2.d0*DY(j))
     &       +(deltaz**2*du2dz(I,K,J)
     &         +deltay**2*du3dy(I,K,J))
     &         *(dthetady(I,K,J,N))
     &         *(dthetadz(I,K,J,N)*DYF(j-1)+dthetadz(I,K,J-1,N)*DYF(j))
     &                               /(2.d0*DY(j))   
     &         )

        IF (TEMP_th(I,K,J) <= 0.0d0) THEN
          TEMP_th(I,K,J)=0.0d0
        ELSE

        TEMP_th(I,K,J)=TEMP_th(I,K,J)/
     &  ((deltax*(dthetadx(I,K,J,N)*DYF(j-1)+dthetadx(I,K,J-1,N)*DYF(j))
     &                               /(2.d0*DY(j)))**2
     &  +(deltay*dthetady(I,K,J,N))**2
     &  +(deltaz*(dthetadz(I,K,J,N)*DYF(j-1)+dthetadz(I,K,J-1,N)*DYF(j))
     &                               /(2.d0*DY(j)))**2)

        END IF



          END DO
        END DO
      END DO


! Now, compute S1_th*dthetadx_i, storing in dthetadx_i
! Need only to compute at GYF points
      DO J=JSTART_TH(N),JEND_TH(N)
       DO K=0,NZP-1
        DO I=0,NXM
      dthetadx(I,K,J,N)=S1_th(I,K,J)*dthetadx(I,K,J,N)
! dthetady is added through an implicit eddy diffusivity
      dthetadz(I,K,J,N)=S1_th(I,K,J)*dthetadz(I,K,J,N)
          END DO
        END DO
      END DO

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


! Get the eddy diffusivity at GY points
      DO J=2,NY
        DO K=0,NZP-1
          DO I=0,NXM
           KAPPA_T(I,K,J,N)=-0.5d0*DELTA_Y(J)*TEMP_th(I,K,J)
          END DO
        END DO
      END DO

! Now that we have calculated KAPPA_T, set the value at the ghost cells
! by sharing with neighboring processes.  This subroutine also sets
! the value of KAPPA_T to zero at both walls
      CALL GHOST_LES_MPI_KAPPA_T

! Convert the scalar flux tensor to Fourier space

!        DO N=1,N_TH
          S1(:,:,:)=dthetadx(:,:,:,N)
          CALL FFT_XZ_TO_FOURIER(S1,CS1,0,NY+1)
          Cdthetadx(:,:,:,N)=CS1(:,:,:)
! dthetady is added through an implicit eddy diffusivity
          S1(:,:,:)=dthetadz(:,:,:,N)
          CALL FFT_XZ_TO_FOURIER(S1,CS1,0,NY+1)
          Cdthetadz(:,:,:,N)=CS1(:,:,:)
!        END DO

! Now, compute TAU, store in the corresponging dthetadx_i
      DO J=1,NY+1
        DO K=0,TNKZ
          DO I=0,NXP-1
            Cdthetadx(I,K,J,N)=0.5d0*DELTA_YF(J)*Cdthetadx(I,K,J,N)
! Cthetady(:,:,:,N) is added through an implicit eddy diffusivity
!            Cdthetady(I,K,J,N)=DELTA_YF(J)*Cdthetady(I,K,J,N)
            Cdthetadz(I,K,J,N)=0.5d0*DELTA_YF(J)*Cdthetadz(I,K,J,N)
          END DO
        END DO
      END DO

        END DO

      end if


         DO N=1,N_TH
      DO J=JSTART_TH(N),JEND_TH(N)
        DO K=0,TNKZ
          DO I=0,NXP-1
            CFTH(I,K,J,N)=CFTH(I,K,J,N)
     &                -CIKX(I)*Cdthetadx(I,K,J,N)
! dthetady is added through implicit eddy diffusivity
     &                -CIKZ(K)*Cdthetadz(I,K,J,N)
          END DO
        END DO
      END DO
         END DO

! Periodically, output mean quantities
      IF ((MOD(TIME_STEP,SAVE_STATS_INT).EQ.0).AND.(RK_STEP.EQ.1)) THEN
! Get plane averages
        do j=0,NY+1
          S1_mean(j)=0.d0
          NU_T_mean(j)=0.d0
          KAPPA_T_mean(J)=0.d0
          do i=0,NXM
          do k=0,NZP-1
            S1_mean(j)=S1_mean(j)+S1(I,K,J)
            NU_T_mean(j)=NU_T_mean(j)+NU_T(I,K,J)
! Only set up for single scalar (N = 1)
            KAPPA_T_mean(j)=KAPPA_T_mean(j)+KAPPA_T(I,K,J,1)
          end do
          end do
        end do

      call mpi_allreduce(mpi_in_place,S1_mean,NY+2,MPI_DOUBLE_PRECISION,
     &     MPI_SUM,MPI_COMM_Z,ierror)
      call mpi_allreduce(mpi_in_place,NU_T_mean,NY+2
     &    ,MPI_DOUBLE_PRECISION,
     &     MPI_SUM,MPI_COMM_Z,ierror)
      call mpi_allreduce(mpi_in_place,KAPPA_T_mean,NY+2
     &    ,MPI_DOUBLE_PRECISION,
     &     MPI_SUM,MPI_COMM_Z,ierror)

        do j=0,NY+1
          S1_mean(j)=S1_mean(j)/dble(NX*NZ)
          NU_T_mean(j)=NU_T_mean(j)/dble(NX*NZ)
          KAPPA_T_mean(j)=KAPPA_T_mean(j)/dble(NX*NZ)
        end do

#ifdef HDF5
      FNAME='mean_les.h5'
      IF (RANKZ.eq.0) then
      gname='kappa_sgs'
      Diag=KAPPA_T_mean(1:NY)
      call WriteStatH5(FNAME,gname,Diag)
      END IF

#else

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
     &      NU_T_mean(J),
     &      KAPPA_T_mean(J)
        end do
      END IF
420     format(I3,' ',2(F20.9,' '),2(F20.9,' '))
      END IF

#endif
      END IF

999   continue

      RETURN
      END


      subroutine compute_all_gradients_chan
C This subroutine computes all gradients for the filtered velocity field
C and scalar field (testing with one scalar field only).
C The input velocity + scalar field should be in fourier space in the periodic
C directions.
C For use in the LES model in channel flow (2 periodic directions)
      include 'header'
      include 'header_les_th'

      integer I,J,K,ij,N

      DO J=1,NY
        DO K=0,TNKZ
          DO I=0,NXP-1


            Cdu1dx(I,K,J)=CIKX(I)*CU1(I,K,J)
            Cdu1dz(I,K,J)=CIKZ(I)*CU1(I,K,J)

            Cdu2dx(I,K,J)=CIKX(I)*CU2(I,K,J)
            Cdu2dy(I,K,J)=(CU2(I,K,J+1)-CU2(I,K,J))/DYF(J)
            Cdu2dz(I,K,J)=CIKZ(I)*CU2(I,K,J)

            Cdu3dx(I,K,J)=CIKX(I)*CU3(I,K,J)
            Cdu3dz(I,K,J)=CIKZ(K)*CU3(I,K,J)
 
          END DO
        END DO
      END DO
      DO J=1,NY+1
        DO K=0,TNKZ
          DO I=0,NXP-1
          
            Cdu1dy(I,K,J)=(CU1(I,K,J)-CU1(I,K,J-1))/DY(J)

            Cdu3dy(I,K,J)=(CU3(I,K,J)-CU3(I,K,J-1))/DY(J)

          END DO
        END DO
      END DO



         DO N=1,N_TH

      DO J=1,NY
        DO K=0,TNKZ
          DO I=0,NXP-1

            Cdthetadx(I,K,J,N)=CIKX(I)*CTH(I,K,J,N)
            Cdthetadz(I,K,J,N)=CIKZ(I)*CTH(I,K,J,N)
          
         END DO
        END DO
      END DO
      DO J=1,NY+1
        DO K=0,TNKZ
          DO I=0,NXP-1
            
            Cdthetady(I,K,J,N)=(CTH(I,K,J,N)-CTH(I,K,J-1,N))/DY(J)

          END DO
        END DO
      END DO



         END DO


       CALL FFT_XZ_TO_PHYSICAL(Cdu1dx,du1dx,0,NY+1)
       CALL FFT_XZ_TO_PHYSICAL(Cdu2dy,du2dy,0,NY+1)
       CALL FFT_XZ_TO_PHYSICAL(Cdu3dz,du3dz,0,NY+1)

       CALL FFT_XZ_TO_PHYSICAL(Cdu2dx,du2dx,0,NY+1)
       CALL FFT_XZ_TO_PHYSICAL(Cdu3dx,du3dx,0,NY+1)

       CALL FFT_XZ_TO_PHYSICAL(Cdu1dz,du1dz,0,NY+1)
       CALL FFT_XZ_TO_PHYSICAL(Cdu2dz,du2dz,0,NY+1)

       CALL FFT_XZ_TO_PHYSICAL(Cdu1dy,du1dy,0,NY+1)
       CALL FFT_XZ_TO_PHYSICAL(Cdu3dy,du3dy,0,NY+1)

        DO N=1,N_TH
          CS1(:,:,:)=Cdthetadx(:,:,:,N)
          CALL FFT_XZ_TO_PHYSICAL(CS1,S1,0,NY+1)
          dthetadx(:,:,:,N)=S1(:,:,:)
       
          CS1(:,:,:)=Cdthetady(:,:,:,N)
          CALL FFT_XZ_TO_PHYSICAL(CS1,S1,0,NY+1)
          dthetady(:,:,:,N)=S1(:,:,:)
       
          CS1(:,:,:)=Cdthetadz(:,:,:,N)
          CALL FFT_XZ_TO_PHYSICAL(CS1,S1,0,NY+1)
          dthetadz(:,:,:,N)=S1(:,:,:)
       
        END DO

!We now have all the gradients in physical space


      RETURN
      END




      SUBROUTINE APPLY_BC_SCALAR_LES
      include 'header'
      integer i,j,k,N

      DO N=1,N_TH

! If we are using Neuman boundary conditions, over-write the values of the
! scalar  at the ghost cells so that the LES model doesn't use the large
! scalar gradient
      IF (TH_BC_YMAX(N).eq.1) THEN
        IF ((RANKY.eq.NPROCY-1).or.(.NOT.USE_MPI)) THEN
! We are on process at the upper wall
          DO K=0,TNKZ
           DO I=0,NXP-1
             CTH(I,K,NY,N)=CTH(I,K,NY-1,N)
             CTH(I,K,NY+1,N)=CTH(I,K,NY,N)
           END DO
          END DO
        END IF
      END IF

      IF (TH_BC_YMIN(N).eq.1) THEN
        IF ((RANKY.eq.0).or.(.NOT.USE_MPI)) THEN
! We are on a process at the bottom wall          
         DO K=0,TNKZ
           DO I=0,NXP-1
             CTH(I,K,1,N)=CTH(I,K,2,N)
             CTH(I,K,0,N)=CTH(I,K,1,N)
           END DO
         END DO
       END IF
      END IF

      END DO

      RETURN
      END




