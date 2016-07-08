C******************************************************************************|
C cavity.f, the cavity-flow solvers for diablo.                    VERSION 0.9
C
C These solvers were written by ? and ? (spring 2001).
C******************************************************************************|

C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE INIT_CAV
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|

      RETURN
      END

C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE RK_CAV_1
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
C Main time-stepping algorithm for the channel-flow case.
C INPUTS  (in Fourier space):  CUi, P, and (if k>1) CFi at (k-1)  (for i=1,2,3)
C OUTPUTS (in Fourier space):  CUi, P, and (if k<3) CFi at (k)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|

      RETURN
      END

C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE RK_CAV_2
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
C Alternative time-stepping algorithm for the channel-flow case.
C INPUTS  (in Fourier space):  CUi, P, and (if k>1) CFi at (k-1)  (for i=1,2,3)
C OUTPUTS (in Fourier space):  CUi, P, and (if k<3) CFi at (k)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|

      RETURN
      END

C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE REM_DIV_CAV
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|

      RETURN
      END

C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE POISSON_P_CAV
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|

      RETURN
      END

C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE CREATE_GRID_CAV
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      include 'header'
      integer I,J,K

         WRITE (6,*) 'Finite-difference in X'
         OPEN (30,file='xgrid.txt',form='formatted',status='old')
         READ (30,*) NX_T
C Check to make sure that grid file is the correct dimensions
         IF (NX_T.ne.NX) THEN
           WRITE(6,*) 'NX, NX_T',NX,NX_T
           STOP 'Error: xgrid.txt wrong dimensions'
         END IF
         DO I=1,NX+1
           READ(30,*) GX(I)
           IF (VERBOSITY .GT. 3) WRITE(6,*) 'GX(',I,') = ',GX(I)
         END DO 
         DO I=1,NX
           READ(30,*) GXF(I)
           IF (VERBOSITY .GT. 3) WRITE(6,*) 'GXF(',I,') = ',GXF(I)
         END DO
C Define ghost cells, if needed for this grid...
         GXF(0)=2.d0*GXF(1)-GXF(2)
         GXF(NX+1)=2.d0*GXF(NX)-GXF(NXM)
         GX(0)=2.d0*GX(1)-GX(2)
C Define the grid spacings
         DO I=1,NX+1
           DX(I)=(GXF(I)-GXF(I-1))
         END DO
         DO I=1,NX
           DXF(I)=(GX(I+1)-GX(I))
         END DO
         CLOSE(30)
         WRITE (6,*) 'Finite-difference in Z'
         OPEN (30,file='zgrid.txt',form='formatted',status='old')
         READ (30,*) NZ_T
C Check to make sure that grid file is the correct dimensions
         IF (NZ_T.ne.NZ) THEN
           WRITE(6,*) 'NZ, NZ_T',NZ,NZ_T
           STOP 'Error: zgrid.txt wrong dimensions'
         END IF
         DO K=1,NZ+1
           READ(30,*) GZ(k)
           IF (VERBOSITY .GT. 3) WRITE(6,*) 'GZ(',K,') = ',GZ(K)
         END DO
         DO K=1,NZ
           READ(30,*) GZF(k)
           IF (VERBOSITY .GT. 3) WRITE(6,*) 'GZF(',K,') = ',GZF(K)
         END DO
C Define ghost cells, if needed for this grid...
         GZF(0)=2.d0*GZF(1)-GZF(2)
         GZF(NZ+1)=2.d0*GZF(NZ)-GZF(NZM)
         GZ(0)=2.d0*GZ(1)-GZ(2)
C Define grid spacing
         DO K=1,NZ+1
           DZ(K)=(GZF(K)-GZF(K-1))
         END DO
         DO K=1,NZ
           DZF(K)=(GZ(K+1)-GZ(K))
         END DO
         CLOSE(30)

         WRITE (6,*) 'Finite-difference in Y'
         OPEN (30,file='./ygrid.txt',form='formatted',status='old')
         READ (30,*) NY_T
C Check to make sure that grid file is the correct dimensions
         IF (NY_T.ne.NY) THEN
           WRITE(6,*) 'NY, NY_T',NY,NY_T
           STOP 'Error: ygrid.txt wrong dimensions'
         END IF
         DO J=1,NY+1
           READ(30,*) GY(j)
           IF (VERBOSITY .GT. 3) WRITE(6,*) 'GY(',J,') = ',GY(J)
         END DO
         DO J=1,NY
           READ(30,*) GYF(j)
           IF (VERBOSITY .GT. 3) WRITE(6,*) 'GYF(',J,') = ',GYF(J)
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

C Define grid spacing
         DO J=1,NY+1
           DY(J)=(GYF(J)-GYF(J-1))
         END DO
         DO J=1,NY
           DYF(J)=(GY(J+1)-GY(J))
         END DO
       RETURN
       END

C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE INPUT_CAV 
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      RETURN
      END

C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE SAVE_STATS_CAV(FINAL)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      LOGICAL FINAL

      RETURN
      END
      

