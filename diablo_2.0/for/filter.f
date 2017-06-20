      subroutine filter_chan(n)
C This subroutine applies a filter to the highest wavenumbers
C It should be applied to the scalar in Fourier space
C The filter used is a sharpened raised cosine filter in the horizontal
C and a fourth order implicit compact filter in the vertical, with the
C parameter alpha determining the width of the vertical filtering window

      include 'header'

      integer I,J,K,js,je,n

! Variables for horizontal filtering
      real*8 sigma(0:NXP-1,0:TNKZ),sigma0

! Variables for vertical filtering
      real*8 alpha
      parameter (alpha=0.0d0)
! Parameters for a larger stencil filter
      real*8 f_a,f_b,f_c

C Set the filtering constants for the horizontal direction
      DO i=0,NXP-1
       DO k=0,TNKZ
        sigma0=0.5d0*(1.d0+
     &       cos(sqrt((KX(i)*LX*1.d0/float(NX))**2.d0
     &            +(KZ(k)*LZ*1.d0/float(NZ))**2.d0)))
! Apply a sharpened raised cosine filter
        sigma(i,k)=sigma0**4.d0*(35.d0-84.d0*sigma0
     &        +70.d0*sigma0**2.d0-20.d0*sigma0**3.d0)
       END DO
      END DO

C Do the spectral filtering in the horizontal
        DO K=0,TNKZ
          DO I=0,NXP-1
            DO J=JSTART_TH(N),JEND_TH(N)
              CTH(I,K,J,n)=CTH(I,K,J,n)*sigma(i,k)
            END DO
          END DO 
        END DO

C Filter the passive scalar, TH in the vertical direction
C Set the filtering constants
!      f_a=(1.d0/8.d0)*(5.d0+6.d0*alpha)
!      f_b=0.5d0*(1.d0+2.d0*alpha)
!      f_c=(-1.d0/8.d0)*(1.d0-2.d0*alpha)      
C First, zero the tridiagonal matrix components
!      DO I=0,NKX
!        DO J=1,NY
!          MATD_C(I,J)=1.d0
!          MATL_C(I,J)=0.d0
!          MATU_C(I,J)=0.d0
!          VEC_C(I,J)=0.d0
!        END DO
!      END DO  
!      DO K=1,TNKZ
!        DO I=1,NKX
C Construct the centered difference terms
!          DO J=2,NY-1
!            MATL_C(I,J)=alpha
!            MATD_C(I,J)=1.d0
!            MATU_C(I,J)=alpha
!            VEC_C(I,J)=f_a*CTH(I,K,J,n)
!     &                +(f_b/2.d0)*(CTH(I,K,J+1,n)+CTH(I,K,J-1,n))
!     &                +(f_c/2.d0)*(CTH(I,K,J+2,n)+CTH(I,K,J-2,n))
!          END DO
C Now, construct the equations for the boundary nodes
!          J=1
!            MATL_C(I,J)=0.d0
!            MATD_C(I,J)=1.d0
!            MATU_C(I,J)=0.d0
!            VEC_C(I,J)=CTH(I,K,J,n)
!          J=NY
!            MATL_C(I,J)=0.d0
!            MATD_C(I,J)=1.d0
!            MATU_C(I,J)=0.d0
!            VEC_C(I,J)=CTH(I,K,J,n)
!         END DO
C Now, solve the tridiagonal system
!         CALL THOMAS_COMPLEX(MATL_C,MATD_C,MATU_C,VEC_C,NY,NKX)
!         DO I=1,NKX
!           DO J=JSTART_TH(N),JEND_TH(N)
!             CTH(I,K,J,n)=VEC_C(I,J)
!           END DO
!         END DO
C END DO K  
!       END DO


       return
       end 





