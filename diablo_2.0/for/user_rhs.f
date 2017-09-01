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

! for plume test (sponge_vel)
      REAL*8 sigma, LYC, LYP, U0, TH0, tau_sponge, RNUM1, PIPE
      REAL*8 sigma_y(0:NY+1),U0_y(0:NY+1)

      sigma=(10.d-3)/2.d0
      LYC=0.05d0
      LYP=0.004d0
      U0=0.228d0
      TH0=0.1659d0

! Optionally, perturb radius of plume
! Also, scale U0 to keep the buoyancy flux constant
      DO J=2,NY
       sigma_y(J)=sigma*(1.d0+0.1d0*sin(2.d0*PI*(GY(J)-U0*TIME)/0.05d0))
       U0_y(J)=U0*(sigma**2.d0)/(sigma_y(J)**2.d0)
      END DO

! Ramp down the forcing timescale for the first second
      if (TIME.lt.1.d0) then
        tau_sponge=100.d0-TIME*(100.d0-0.1d0)/1.d0
      else
         tau_sponge = 0.1d0
      end if

! Create the damping function
      DO J=2,NY
        DO K=0,NZP-1
          DO I=0,NXM
           CALL RANDOM_NUMBER(RNUM1)
! Gaussian
             S1(I,K,J)=(U2(I,K,J)-U0*
     &     EXP(-((GX(I)-LX/2.d0)**2.d0+(GZ(RANKZ*NZP+K)-LZ/2.d0)**2.d0)
     &                  /(2.d0*(sigma)**2.d0)) 
     &                 *(1.d0+2.d0*(RNUM1-0.5d0)))
     &               *(1.d0-tanh((GY(J)-LYC)/LYP))/2.d0

          END DO
        END DO
      END DO

      CALL FFT_XZ_TO_FOURIER(S1,CS1,0,NY+1)

      DO J=2,NY
        DO K=0,TNKZ
          DO I=0,NXP-1
            CF2(I,K,J)=CF2(I,K,J) - CS1(I,K,J)/tau_sponge
          END DO
        END DO
      END DO

! Create the damping function
!      DO N=1,N_TH
      DO N=1,1
      DO J=JSTART_TH(N),JEND_TH(N)
        DO K=0,NZP-1
          DO I=0,NXM
           CALL RANDOM_NUMBER(RNUM1)
! Gaussian
             S1(I,K,J)=(TH(I,K,J,N)-TH0*
     &     EXP(-((GX(I)-LX/2.d0)**2.d0+(GZ(RANKZ*NZP+K)-LZ/2.d0)**2.d0)
     &                       /(2.d0*(sigma)**2.d0)) 
     &                 *(1.d0+2.d0*(RNUM1-0.5d0)))
     &               *(1.d0-tanh((GY(J)-LYC)/LYP))/2.d0

          END DO
        END DO
      END DO

      CALL FFT_XZ_TO_FOURIER(S1,CS1,0,NY+1)

      DO J=JSTART_TH(N),JEND_TH(N)
        DO K=0,TNKZ
          DO I=0,NXP-1
            CFTH(I,K,J,N)=CFTH(I,K,J,N) - CS1(I,K,J)/tau_sponge
          END DO
        END DO
      END DO
      END DO

! For example, to add a linear damping term (e.g. -alpha*U) to the RHS:
!       alpha=-0.1d0
!       DO J=JSTART,JEND
!         DO K=0,NZP-1
!           DO I=0,NXM
!             S1(I,K,J)=-alpha*U1(I,K,J)
!           END DO
!         END DO
!       END DO
! Convert to Fourier space
!       CALL FFT_XZ_TO_FOURIER(S1,CS1,0,NY+1)
!       DO J=JSTART,JEND
!         DO K=0,TNKY 
!           DO I=0,NXP-1
!             CF1(I,K,J)=CF1(I,K,J)+CS1(I,K,J)
!           END DO
!         END DO
!       END DO

! For U2 do this...
! Note that the only thing that changes are the bounds of the J index
!       DO J=2,NY
!         DO K=0,NZP-1
!           DO I=0,NXM
!             S1(I,K,J)=-alpha*U2(I,K,J)
!           END DO
!         END DO
!       END DO
! Convert to Fourier space
!       CALL FFT_XZ_TO_FOURIER(S1,CS1,0,NY+1)
!       DO J=2,NY
!         DO K=0,TNKY
!           DO I=0,NXP-1
!             CF2(I,K,J)=CF2(I,K,J)+CS1(I,K,J)
!           END DO
!         END DO
!       END DO

! For scalars, do this...
! Loop over all scalars
!       DO N=1,N_TH
!       DO J=JSTART,JEND
!         DO K=0,NZP-1
!           DO I=0,NXM
!             S1(I,K,J)=-alpha*TH(I,K,J,N)
!           END DO
!         END DO
!       END DO
! Convert to Fourier space
!       CALL FFT_XZ_TO_FOURIER(S1,CS1,0,NY+1)
!       DO J=JSTART,JEND
!         DO K=0,TNKY
!           DO I=0,NXP-1
!             CFTH(I,K,J,N)=CFTH(I,K,J,N)+CS1(I,K,J)
!           END DO
!         END DO
!       END DO
!      END DO 

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

       real*8 alpha

! For example, to add a linear damping term (e.g. -alpha*U) to the RHS:
!       alpha=-0.1d0
!       DO J=JSTART,JEND
!         DO K=0,TNKY 
!           DO I=0,NXP-1
!             CF1(I,K,J)=CF1(I,K,J)-alpha*CU1(I,K,J)
!           END DO
!         END DO
!       END DO

! For U2 do this...
! Note that the only thing that changes are the bounds of the J index
!       DO J=2,NY
!         DO K=0,TNKY
!           DO I=0,NXP-1
!             CF2(I,K,J)=CF2(I,K,J)-alpha*CU2(I,K,J)
!           END DO
!         END DO
!       END DO

! For scalars, do this...
!       DO J=JSTART,JEND
!         DO K=0,TNKY
!           DO I=0,NXP-1
!             CFTH(I,K,J,N)=CFTH(I,K,J,N)-alpha*CTH(I,K,J,N)
!           END DO
!         END DO
!       END DO
!      END DO 

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






