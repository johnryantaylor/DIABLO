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






