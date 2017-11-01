      SUBROUTINE GHOST_CHAN_MPI
! This subroutine is part of the MPI package for the channel flow
! Diablo package.
! Here, we define a set of ghost cells on each process
! the ghost cells contain information from the neighboring nodes
! and allow us to compute finite differences over the local gridpoints.
! We need to update the contents of the ghost cells at the start of
! each Runge-Kutta substep

      include 'header'

      integer i,j,k,N

! Define the arrays that will be used for data packing.  This makes the
! communication between processes more efficient by only requiring one
! send and recieve.
! The communication will be done in Fourier space, so these arrays should
! be complex arrays to match the velocity
! The size of the buffer array is 0:NKX,0:TNKZ,# of variables
      COMPLEX*16 OCPACK(0:NXP-1,0:TNKZ,4+N_TH)
      COMPLEX*16 ICPACK(0:NXP-1,0:TNKZ,4+N_TH)

! If we are using more than one processor, then we need to pass data

      IF (NPROCY.gt.1) THEN

! First, Pass data up the chain to higher ranked processes

      IF (RANKY.eq.0) THEN
! If we are the lowest ranked process, then we don't need to recieve
! data at the lower ghost cells, these will be filled with boundary
! condition information
        DO K=0,TNKZ
          DO I=0,NXP-1
            OCPACK(I,K,1)=CU1(I,K,NY)
            OCPACK(I,K,2)=CU2(I,K,NY)
            OCPACK(I,K,3)=CU3(I,K,NY)
            OCPACK(I,K,4)=CP(I,K,NY)
            DO N=1,N_TH
              OCPACK(I,K,4+N)=CTH(I,K,NY,N)
            END DO
          END DO
        END DO
! Now, we have packed the data into a compact array, pass the data up
        CALL MPI_SEND(OCPACK,(4+N_TH)*(NXP)*(TNKZ+1)
     &               ,MPI_DOUBLE_COMPLEX
     &               ,RANKY+1,1,MPI_COMM_Y,IERROR)

! End if RANK=0
      ELSE IF (RANKY.lt.NPROCY-1) THEN
! Here, we are one of the middle processes and we need to pass data
! up and recieve data from below
        DO K=0,TNKZ
          DO I=0,NXP-1
            OCPACK(I,K,1)=CU1(I,K,NY)
            OCPACK(I,K,2)=CU2(I,K,NY)
            OCPACK(I,K,3)=CU3(I,K,NY)
            OCPACK(I,K,4)=CP(I,K,NY)
            DO N=1,N_TH
              OCPACK(I,K,4+N)=CTH(I,K,NY,N)
            END DO
          END DO
        END DO
! Use MPI_SENDRECV since we need to recieve and send data
        CALL MPI_SEND(OCPACK,(4+N_TH)*(NXP)*(TNKZ+1)
     &               ,MPI_DOUBLE_COMPLEX
     &               ,RANKY+1,1,MPI_COMM_Y,IERROR)

        CALL MPI_RECV(ICPACK,(4+N_TH)*(NXP)*(TNKZ+1)
     &               ,MPI_DOUBLE_COMPLEX
     &               ,RANKY-1,1,MPI_COMM_Y,STATUS,IERROR)
! Now, unpack the data that we have recieved
        DO K=0,TNKZ
          DO I=0,NXP-1
            CU1(I,K,1)=ICPACK(I,K,1)
            CU2(I,K,1)=ICPACK(I,K,2)
            CU3(I,K,1)=ICPACK(I,K,3)
            CP(I,K,1)=ICPACK(I,K,4)
            DO N=1,N_TH
              CTH(I,K,1,N)=ICPACK(I,K,4+N)
            END DO
          END DO
        END DO

      ELSE
! Otherwise, we must be the uppermost process with RANK=NPROCS-1
! Here, we need to recieve data from below, but don't need to send data up
        CALL MPI_RECV(ICPACK,(4+N_TH)*(NXP)*(TNKZ+1)
     &               ,MPI_DOUBLE_COMPLEX
     &               ,RANKY-1,1,MPI_COMM_Y,STATUS,IERROR)
! Unpack the data that we have recieved
        DO K=0,TNKZ
          DO I=0,NXP-1
            CU1(I,K,1)=ICPACK(I,K,1)
            CU2(I,K,1)=ICPACK(I,K,2)
            CU3(I,K,1)=ICPACK(I,K,3)
            CP(I,K,1)=ICPACK(I,K,4)
            DO N=1,N_TH
              CTH(I,K,1,N)=ICPACK(I,K,4+N)
            END DO
          END DO
        END DO
      END IF
      
! AT this point we have passed data up the chain
      IF (RANKY.eq.NPROCY-1) THEN
! If we are the higest ranked process, then we don't need to recieve
! data at the upper ghost cells, these will be filled with boundary
! condition information
        DO K=0,TNKZ
          DO I=0,NXP-1
            OCPACK(I,K,1)=CU1(I,K,2)
            OCPACK(I,K,2)=CU2(I,K,2)
            OCPACK(I,K,3)=CU3(I,K,2)
            OCPACK(I,K,4)=CP(I,K,2)
            DO N=1,N_TH
              OCPACK(I,K,4+N)=CTH(I,K,2,N)
            END DO
          END DO
        END DO
! Now, we have packed the data into a compact array, pass the data up
        CALL MPI_SEND(OCPACK,(4+N_TH)*(NXP)*(TNKZ+1)
     &               ,MPI_DOUBLE_COMPLEX
     &               ,RANKY-1,3,MPI_COMM_Y,IERROR)
      ELSE IF (RANKY.GT.0) THEN
! Here, we are one of the middle processes and we need to pass data
! down and recieve data from above us
        DO K=0,TNKZ
          DO I=0,NXP-1
            OCPACK(I,K,1)=CU1(I,K,2)
            OCPACK(I,K,2)=CU2(I,K,2)
            OCPACK(I,K,3)=CU3(I,K,2)
            OCPACK(I,K,4)=CP(I,K,2)
            DO N=1,N_TH
              OCPACK(I,K,4+N)=CTH(I,K,2,N)
            END DO
          END DO
        END DO

        CALL MPI_SEND(OCPACK,(4+N_TH)*(NXP)*(TNKZ+1)
     &               ,MPI_DOUBLE_COMPLEX
     &               ,RANKY-1,3,MPI_COMM_Y,IERROR)

        CALL MPI_RECV(ICPACK,(4+N_TH)*(NXP)*(TNKZ+1)
     &               ,MPI_DOUBLE_COMPLEX
     &               ,RANKY+1,3,MPI_COMM_Y,STATUS,IERROR)
! Now, unpack the data that we have recieved
        DO K=0,TNKZ
          DO I=0,NXP-1
            CU1(I,K,NY+1)=ICPACK(I,K,1)
            CU2(I,K,NY+1)=ICPACK(I,K,2)
            CU3(I,K,NY+1)=ICPACK(I,K,3)
            CP(I,K,NY+1)=ICPACK(I,K,4)
            DO N=1,N_TH
              CTH(I,K,NY+1,N)=ICPACK(I,K,4+N)
            END DO
          END DO
        END DO
      ELSE
! Here, we must be the lowest process (RANK=0) and we need to recieve
! data from above
        CALL MPI_RECV(ICPACK,(4+N_TH)*(NXP)*(TNKZ+1)
     &               ,MPI_DOUBLE_COMPLEX
     &               ,RANKY+1,3,MPI_COMM_Y,STATUS,IERROR)
! Unpack the data that we have recieved
        DO K=0,TNKZ
          DO I=0,NXP-1
            CU1(I,K,NY+1)=ICPACK(I,K,1)
            CU2(I,K,NY+1)=ICPACK(I,K,2)
            CU3(I,K,NY+1)=ICPACK(I,K,3)
            CP(I,K,NY+1)=ICPACK(I,K,4)
            DO N=1,N_TH
              CTH(I,K,NY+1,N)=ICPACK(I,K,4+N)
            END DO
          END DO
        END DO
      END IF

      END IF

      RETURN
      END


      SUBROUTINE GHOST_LES_MPI
! This subroutine is part of the MPI package for the LES subroutine
! Here, after calculating the SGS viscosity, NU_T on each core,
! We need to share the ghost cells between neighboring processors

      include 'header'

      integer i,j,k,N

! Define the arrays that will be used for data packing.  This makes the
! communication between processes more efficient by only requiring one
! send and recieve.
! The communication will be done in Fourier space, so these arrays should
! be complex arrays to match the velocity
! The size of the buffer array is 0:NXM,0:NZP-1
      REAL*8 OCPACK(0:NXM,0:NZP-1)
      REAL*8 ICPACK(0:NXM,0:NZP-1)

! If we are using more than one processor, then we need to pass data

      IF (NPROCY.gt.1) THEN

! First, Pass data up the chain to higher ranked processes

      IF (RANKY.eq.0) THEN
! If we are the lowest ranked process, then we don't need to recieve
! data at the lower ghost cells. Instead, set NU_T=0 at the lower wall
        DO K=0,NZP-1
          DO I=0,NXM
            NU_T(I,K,1)=0.d0
            NU_T(I,K,2)=0.d0
          END DO
        END DO

! Pass data up to the next process from GY(NY)
        DO K=0,NZP-1
          DO I=0,NXM
            OCPACK(I,K)=NU_T(I,K,NY)
          END DO
        END DO
! Now, we have packed the data into a compact array, pass the data up
        CALL MPI_SEND(OCPACK,(NXM+1)*(NZP)
     &               ,MPI_DOUBLE_PRECISION
     &               ,RANKY+1,1,MPI_COMM_Y,IERROR)

! End if RANK=0
      ELSE IF (RANKY.lt.NPROCY-1) THEN
! Here, we are one of the middle processes and we need to pass data
! up and recieve data from below
        DO K=0,NZP-1
          DO I=0,NXM
            OCPACK(I,K)=NU_T(I,K,NY)
          END DO
        END DO
! Use MPI_SENDRECV since we need to recieve and send data
        CALL MPI_SEND(OCPACK,(NXM+1)*(NZP)
     &               ,MPI_DOUBLE_PRECISION
     &               ,RANKY+1,1,MPI_COMM_Y,IERROR)

        CALL MPI_RECV(ICPACK,(NXM+1)*(NZP)
     &               ,MPI_DOUBLE_PRECISION
     &               ,RANKY-1,1,MPI_COMM_Y,STATUS,IERROR)
! Now, unpack the data that we have recieved
        DO K=0,NZP-1
          DO I=0,NXM
            NU_T(I,K,1)=ICPACK(I,K)
          END DO
        END DO

      ELSE
! Otherwise, we must be the uppermost process with RANK=NPROCS-1
! Here, we need to recieve data from below, but don't need to send data up
        CALL MPI_RECV(ICPACK,(NXM+1)*(NZP)
     &               ,MPI_DOUBLE_PRECISION
     &               ,RANKY-1,1,MPI_COMM_Y,STATUS,IERROR)
! Unpack the data that we have recieved
        DO K=0,NZP-1
          DO I=0,NXM
            NU_T(I,K,1)=ICPACK(I,K)
          END DO
        END DO
      END IF
      
! Now, we have hit the top process.  Set the BCs and pass data down

      IF (RANKY.eq.NPROCY-1) THEN
! If we are the higest ranked process, then we don't need to recieve
! data at the upper ghost cells.
! Set NU_T=0 at the upper wall
        DO K=0,NZP-1
          DO I=0,NXM
            NU_T(I,K,NY)=0.d0
            NU_T(I,K,NY+1)=0.d0
          END DO
        END DO

! Now, send data down the chain
        DO K=0,NZP-1
          DO I=0,NXM
            OCPACK(I,K)=NU_T(I,K,2)
          END DO
        END DO
! Now, we have packed the data into a compact array, pass the data up
        CALL MPI_SEND(OCPACK,(NXM+1)*(NZP)
     &               ,MPI_DOUBLE_PRECISION
     &               ,RANKY-1,3,MPI_COMM_Y,IERROR)
      ELSE IF (RANKY.GT.0) THEN
! Here, we are one of the middle processes and we need to pass data
! down and recieve data from above us
        DO K=0,NZP-1
          DO I=0,NXM
            OCPACK(I,K)=NU_T(I,K,2)
          END DO
        END DO

        CALL MPI_SEND(OCPACK,(NXM+1)*(NZP)
     &               ,MPI_DOUBLE_PRECISION
     &               ,RANKY-1,3,MPI_COMM_Y,IERROR)

        CALL MPI_RECV(ICPACK,(NXM+1)*(NZP)
     &               ,MPI_DOUBLE_PRECISION
     &               ,RANKY+1,3,MPI_COMM_Y,STATUS,IERROR)
! Now, unpack the data that we have recieved
        DO K=0,NZP-1
          DO I=0,NXM
            NU_T(I,K,NY+1)=ICPACK(I,K)
          END DO
        END DO
      ELSE
! Here, we must be the lowest process (RANK=0) and we need to recieve
! data from above
        CALL MPI_RECV(ICPACK,(NXM+1)*(NZP)
     &               ,MPI_DOUBLE_PRECISION
     &               ,RANKY+1,3,MPI_COMM_Y,STATUS,IERROR)
! Unpack the data that we have recieved
        DO K=0,NZP-1
          DO I=0,NXM
            NU_T(I,K,NY+1)=ICPACK(I,K)
          END DO
        END DO
      END IF

      ELSE
! HERE, NPROCY=1, so just apply the boundary conditions to set NU=0 at the
! top and bottom walls
        DO K=0,NZP-1
          DO I=0,NXM
            NU_T(I,K,1)=0.d0
            NU_T(I,K,2)=0.d0
            NU_T(I,K,NY)=0.d0
            NU_T(I,K,NY+1)=0.d0
          END DO
        END DO

      END IF

      RETURN
      END


      SUBROUTINE GHOST_GRID_MPI
! This subroutine is part of the MPI package for the channel flow
! Diablo package.
! Here, we define a set of ghost cells on each process
! the ghost cells contain information from the neighboring nodes
! and allow us to compute finite differences over the local gridpoints.
! We need to update the contents of the ghost cells at the start of
! each Runge-Kutta substep

      include 'header'

      integer i,j,k,N

      real*8 OCPACK(3),ICPACK(3)

! First, Pass data up the chain to higher ranked processes

      IF (RANKY.eq.0) THEN

! Set the lower ghost cells
        GYF(0)=2.d0*GYF(1)-GYF(2)
        GY(0)=2.d0*GY(1)-GY(2)

        OCPACK(1)=GYF(NY-1)
        OCPACK(2)=GYF(NY)
        OCPACK(3)=GY(NY)

        CALL MPI_SEND(OCPACK,3
     &               ,MPI_DOUBLE_PRECISION
     &               ,RANKY+1,1,MPI_COMM_Y,IERROR)

! End if RANK=0
      ELSE IF (RANKY.lt.NPROCY-1) THEN
! Here, we are one of the middle processes and we need to pass data
! up and recieve data from below
        OCPACK(1)=GYF(NY-1)
        OCPACK(2)=GYF(NY)
        OCPACK(3)=GY(NY)

        CALL MPI_SEND(OCPACK,3
     &               ,MPI_DOUBLE_PRECISION
     &               ,RANKY+1,1,MPI_COMM_Y,IERROR)

        CALL MPI_RECV(ICPACK,3
     &               ,MPI_DOUBLE_PRECISION
     &               ,RANKY-1,1,MPI_COMM_Y,STATUS,IERROR)
! Now, unpack the data that we have recieved
        GYF(0)=ICPACK(1)
        GYF(1)=ICPACK(2)
        GY(1)=ICPACK(3)

      ELSE
! Otherwise, we must be the uppermost process with RANK=NPROCS-1
! Set the top ghost cell
        GYF(NY+1)=2.*GYF(NY)-GYF(NYM)

! Here, we need to recieve data from below, but don't need to send data up
        CALL MPI_RECV(ICPACK,3
     &               ,MPI_DOUBLE_PRECISION
     &               ,RANKY-1,1,MPI_COMM_Y,STATUS,IERROR)
! Now, unpack the data that we have recieved
        GYF(0)=ICPACK(1)
        GYF(1)=ICPACK(2)
        GY(1)=ICPACK(3)

      END IF
! AT this point we have passed data up the chain

      IF (RANKY.eq.NPROCY-1) THEN
! If we are the higest ranked process, then we don't need to recieve
! data at the upper ghost cells, these will be filled with boundary
! condition information

        OCPACK(1)=GYF(2)
        OCPACK(2)=GY(2)


! Now, we have packed the data into a compact array, pass the data up
        CALL MPI_SEND(OCPACK,3
     &               ,MPI_DOUBLE_PRECISION
     &               ,RANKY-1,3,MPI_COMM_Y,IERROR)
      ELSE IF (RANKY.GT.0) THEN
! Here, we are one of the middle processes and we need to pass data
! down and recieve data from above us

        OCPACK(1)=GYF(2)
        OCPACK(2)=GY(2)
        CALL MPI_SEND(OCPACK,3
     &               ,MPI_DOUBLE_PRECISION
     &               ,RANKY-1,3,MPI_COMM_Y,IERROR)
        CALL MPI_RECV(ICPACK,3
     &               ,MPI_DOUBLE_PRECISION
     &               ,RANKY+1,3,MPI_COMM_Y,STATUS,IERROR)
! Now, unpack the data that we have recieved
        GYF(NY+1)=ICPACK(1)
        GY(NY+1)=ICPACK(2)

      ELSE
! Here, we must be the lowest process (RANK=0) and we need to recieve
! data from above
        CALL MPI_RECV(ICPACK,3
     &               ,MPI_DOUBLE_PRECISION
     &               ,RANKY+1,3,MPI_COMM_Y,STATUS,IERROR)
! Unpack the data that we have recieved
        GYF(NY+1)=ICPACK(1)
        GY(NY+1)=ICPACK(2)

      END IF

      RETURN
      END

C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE THOMAS_FORWARD_REAL_MPI(A,B,C,G,INY,INX)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
C This version solves for one full row, then passes the data
C This subroutine performs the backward sweep of the Thomas algorithm
C Thomas algorithm solves Ax=b for tridiagonal A
C The RHS vector and solution are real
C Input lower, main, and upper diagonals, ld, md, ud, and rhs x
C Returns solution in x
C The indexing should be done by ROW, ie.
C [ b1  c1   0   0   0 ...
C [ a2  b2  c2   0   0 ...
C [  0  a3  b3   c3  0 ...

      INCLUDE 'header'

      INTEGER I,J,INX,INY
      REAL*8 A(0:INX-1,0:INY+1), B(0:INX-1,0:INY+1), C(0:INX-1,0:INY+1)
     &         , G(0:INX-1,0:INY+1)

      REAL*8 OCPACK(0:INX-1,4),ICPACK(0:INX-1,4)


      IF (RANKY.NE.0) THEN
C If we aren't the lowest process, then wait for data
        CALL MPI_RECV(OCPACK,4*INX,MPI_DOUBLE_PRECISION,RANKY-1,12
     &               ,MPI_COMM_Y,status,ierror)
C Unpack the data
        DO I=0,INX-1
        A(I,1)=OCPACK(I,1)
        B(I,1)=OCPACK(I,2)
        C(I,1)=OCPACK(I,3)
        G(I,1)=OCPACK(I,4)
        END DO
! If we aren't the lowest process, start at J=2
        DO J=2,INY
          DO I=0,INX-1
            A(I,J)=-A(I,J)/B(I,J-1)
            B(I,J)=B(I,J)+A(I,J)*C(I,J-1)
            G(I,J)=G(I,J)+A(I,J)*G(I,J-1)
          END DO
        END DO
      ELSE
! Here, we are the lowest process, start solving at J=1
      DO J=1,INY
        DO I=0,INX-1
        A(I,J)=-A(I,J)/B(I,J-1)
        B(I,J)=B(I,J)+A(I,J)*C(I,J-1)
        G(I,J)=G(I,J)+A(I,J)*G(I,J-1)
        END DO
      END DO
      END IF

      IF (RANKY.NE.NPROCY-1) THEN
        DO I=0,INX-1
        ICPACK(I,1)=A(I,INY)
        ICPACK(I,2)=B(I,INY)
        ICPACK(I,3)=C(I,INY)
        ICPACK(I,4)=G(I,INY)
        END DO
        CALL MPI_SEND(ICPACK,4*INX,MPI_DOUBLE_PRECISION,RANKY+1,12
     &               ,MPI_COMM_Y,ierror)
      ELSE
! Here, we are at the upper process, so solve one more row containing
! the boundary conditions
        J=INY+1
        DO I=0,INX-1
        A(I,J)=-A(I,J)/B(I,J-1)
        B(I,J)=B(I,J)+A(I,J)*C(I,J-1)
        G(I,J)=G(I,J)+A(I,J)*G(I,J-1)
        END DO
      END IF

      RETURN
      END



C----*|--.---------.---------.---------.---------.---------.---------.-|------
      SUBROUTINE THOMAS_FORWARD_COMPLEX_MPI(A,B,C,G,INY,INX)
C----*|--.---------.---------.---------.---------.---------.---------.-|-----|
C This subroutine performs the backward sweep of the Thomas algorithm
C Thomas algorithm solves Ax=b for tridiagonal A
C The RHS vector and solution are real
C Input lower, main, and upper diagonals, ld, md, ud, and rhs x
C Returns solution in x
C The indexing should be done by ROW, ie.
C [ b1  c1   0   0   0 ...
C [ a2  b2  c2   0   0 ...
C [  0  a3  b3   c3  0 ...

      INCLUDE 'header'

      INTEGER I,J,INX,INY
      REAL*8 A(0:INX,0:INY+1), B(0:INX,0:INY+1), C(0:INX,0:INY+1)
      COMPLEX*16 G(0:INX,0:INY+1)

      COMPLEX*16 OCPACK(4),ICPACK(4)

      DO I=0,INX

      IF (RANKY.NE.0) THEN
C If we aren't the lowest process, then wait for data
        CALL MPI_RECV(OCPACK,4,MPI_DOUBLE_COMPLEX,RANKY-1,13
     &               ,MPI_COMM_Y,status,ierror)
C Unpack the data
        A(I,1)=REAL(OCPACK(1))
        B(I,1)=REAL(OCPACK(2))
        C(I,1)=REAL(OCPACK(3))
        G(I,1)=OCPACK(4)
      END IF

      DO J=2,INY
        A(I,J)=-A(I,J)/B(I,J-1)
        B(I,J)=B(I,J)+A(I,J)*C(I,J-1)
        G(I,J)=G(I,J)+A(I,J)*G(I,J-1)
      END DO

      IF (RANKY.NE.NPROCY-1) THEN
        ICPACK(1)=A(I,INY)
        ICPACK(2)=B(I,INY)
        ICPACK(3)=C(I,INY)
        ICPACK(4)=G(I,INY)
        CALL MPI_SEND(ICPACK,4,MPI_DOUBLE_COMPLEX,RANKY+1,13
     &               ,MPI_COMM_Y,ierror)
      END IF


      END DO

      RETURN
      END

C----*|--.---------.---------.---------.---------.---------.---------.-|------
      SUBROUTINE THOMAS_BACKWARD_REAL_MPI(A,B,C,G,INY,INX)
C----*|--.---------.---------.---------.---------.---------.---------.-|------
C This subroutine performs the backward sweep of the Thomas algorithm
C Thomas algorithm solves Ax=b for tridiagonal A
C The RHS vector and solution are real
C Input lower, main, and upper diagonals, ld, md, ud, and rhs x
C Returns solution in x
C The indexing should be done by ROW, ie.
c ￼￼￼C[b1c1 0 C[a2 b2 c2 C[0a3b3
c 0 0... 0 0... c30...

      INCLUDE 'header'

      INTEGER I,J,INX,INY
      REAL*8 A(0:INX-1,0:INY+1), B(0:INX-1,0:INY+1), C(0:INX-1,0:INY+1)
     &     , G(0:INX-1,0:INY+1)
      REAL*8 ICPACK(1),OCPACK(1)
      DO I = 0,INX-1
       IF (RANKY.NE.NPROCY-1) THEN
! If we aren’t the highest process, then wait for data
        CALL MPI_RECV(OCPACK,1,MPI_DOUBLE_PRECISION,RANKY+1,10
     &               ,MPI_COMM_Y,status,ierror)
          G(I,INY+1)=OCPACK(1)
       ELSE
! Else, if we are the highest process, compute the solution at j=INY
          G(I,INY+1)=G(I,INY+1)/B(I,INY+1)
       END IF
! All processes solve from INY..1
      DO J=INY,0,-1
        G(I,J)=(G(I,J)-C(I,J)*G(I,J+1))/B(I,J)
      END DO
      IF (RANKY.NE.0) THEN
          ICPACK(1)=G(I,2)
        CALL MPI_SEND(ICPACK,1,MPI_DOUBLE_PRECISION,RANKY-1,10
     &               ,MPI_COMM_Y,ierror)
      END IF
      END DO
      RETURN
      END

C----*|--.---------.---------.---------.---------.---------.---------.-|------
      SUBROUTINE THOMAS_BACKWARD_COMPLEX_MPI(A,B,C,G,INY,INX)
C----*|--.---------.---------.---------.---------.---------.---------.-|-----|
C This subroutine performs the backward sweep of the Thomas algorithm
C Thomas algorithm solves Ax=b for tridiagonal A
C The RHS vector and solution are real
C Input lower, main, and upper diagonals, ld, md, ud, and rhs x
C Returns solution in x
C The indexing should be done by ROW, ie.
C [ b1  c1   0   0   0 ...
C [ a2  b2  c2   0   0 ...
C [  0  a3  b3   c3  0 ...

      INCLUDE 'header'

      INTEGER I,J,INX,INY
      REAL*8 A(0:INX,0:INY+1), B(0:INX,0:INY+1), C(0:INX,0:INY+1)
      COMPLEX*16 G(0:INX,0:INY+1)

      DO I=0,INX

      IF (RANKY.NE.NPROCY-1) THEN
C If we aren't the highest process, then wait for data
        CALL MPI_RECV(G(I,INY+1),1,MPI_DOUBLE_COMPLEX,RANKY+1,11
     &               ,MPI_COMM_Y,status,ierror)
        J=INY
        G(I,J)=(G(I,J)-C(I,J)*G(I,J+1))/B(I,J)
      ELSE
C Else, if we are the highest process, then compute the solution at j=INY
        G(I,INY)=G(I,INY)/B(I,INY)
      END IF

      DO J=INY-1,1,-1
        G(I,J)=(G(I,J)-C(I,J)*G(I,J+1))/B(I,J)
      END DO

      IF (RANKY.NE.0) THEN
        CALL MPI_SEND(G(I,2),1,MPI_DOUBLE_COMPLEX,RANKY-1,11
     &               ,MPI_COMM_Y,ierror)
      END IF

      END DO

      RETURN
      END

!----*|--.---------.---------.---------.---------.---------.---------.-|------
      SUBROUTINE INIT_MPI
C----*|--.---------.---------.---------.---------.---------.---------.-|-----|
      INCLUDE 'header'
 
      INTEGER bl(2),disp(2),types(2)

      INTEGER IPROCS,TYPE1,COMM_CART
      INTEGER DIMS(2),PERDIM(2)
      INTEGER MFLAG(2)
      
      INTEGER I,J,K,XI,ZI

      COMPLEX*16 TMP(0:NX/2,0:NZP+1,0:NY+1)
      CHARACTER*55 FNAME

      INTEGER         FFTW_FORWARD,      FFTW_BACKWARD,
     *                FFTW_ESTIMATE,     FFTW_MEASURE,
     *                FFTW_OUT_OF_PLACE, FFTW_IN_PLACE,
     *                FFTW_USE_WISDOM,   FFTW_THREADSAFE
      PARAMETER(      FFTW_FORWARD=-1,      FFTW_BACKWARD=1,
     *                FFTW_ESTIMATE=0,      FFTW_MEASURE=1,
     *                FFTW_OUT_OF_PLACE=0,  FFTW_IN_PLACE=8,
     *                FFTW_USE_WISDOM=16,   FFTW_THREADSAFE=128 )

! This subroutine initializes all mpi variables

c$$$      integer(C_INTPTR_T) :: L,M
c$$$      type(C_PTR) :: plan, cdata
c$$$      complex(C_DOUBLE_COMPLEX), pointer :: data(:,:)
c$$$      integer(C_INTPTR_T) :: i, j, alloc_local, local_M, local_j_offset
c$$$      INTEGER(C_INTPTR_T) NZP,iNZ,alloc_local

      CALL MPI_INIT(IERROR)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,IPROCS,IERROR)
      CALL MPI_COMM_RANK(MPI_COMM_WORLD,RANK,IERROR)

      IF (NPROCS.NE.IPROCS) THEN
         IF (RANK.EQ.0) WRITE(*,*) 'ERROR: compiled with ',NPROCS,
     &        'cores, running with ',IPROCS,' cores.'
         CALL MPI_FINALIZE(ierror)
         stop 
      END IF

      IF (MOD(NPROCS,NPROCY).NE.0) THEN
         IF (RANK.EQ.0) WRITE(*,*) ' Error. NPROCS is not a ',
     &        'multiple of NPROCY'
         CALL MPI_FINALIZE(ierror)
         stop 
      END IF

      DIMS(2)=NPROCY
      DIMS(1)=NPROCZ
      MFLAG(:)=.FALSE.

      NXPP=int(NX/(2*NPROCZ))

      call MPI_CART_CREATE(MPI_COMM_WORLD,2,DIMS,PERDIM,.FALSE.,
     &     COMM_CART,IERROR)
      ! In PERDIM I put the information for the remain_dims 
      PERDIM=(/0,1/)
      call MPI_CART_SUB(COMM_CART,PERDIM,MPI_COMM_Y,IERROR)
      PERDIM=(/1,0/)
      call MPI_CART_SUB(COMM_CART,PERDIM,MPI_COMM_Z,IERROR)

      call MPI_COMM_RANK(MPI_COMM_Y,RANKY,IERROR)
      call MPI_COMM_RANK(MPI_COMM_Z,RANKZ,IERROR)

      if (RANK.eq.0) write(*,'(1A,1I8)') 'MPI Initialised. NPROCS: ',
     &     NPROCS
      call MPI_BARRIER(MPI_COMM_WORLD,IERROR)
      write(*,'(1A,4I8)') 'RANK,RANKY,RANKZ: ',RANK,RANKY,RANKZ
      
      !call fftw_mpi_init()
      !alloc_local=fftw_mpi_local_size_2d(NZ,NX, MPI_COMM_Z,
      !&     NZP,iNZ)
      !write(*,*) NZP,iNZ,alloc_local

      !NZP=NZ/NPROCSZ
      !NXP=NX/(2*NPROCSZ)
      
C$$$      PI=4.*atan(1.0)
C$$$      do i=0,NX-1
C$$$         do j=1,1
C$$$            do k=0,NZP-1
C$$$               xi= i    
C$$$               zi= k+NZP*RANKZ
C$$$               V(i,k,j)=sin(2*PI/NX*xi+2*2*PI/NZ*zi)
C$$$            end do
C$$$         end do
C$$$      end do

c$$$c     ------------------------------
c$$$c     Define FFT
c$$$c     ------------------------------
c$$$      CALL RFFTWND_F77_CREATE_PLAN(FFTW_X_TO_F_PLAN, 1, NX,                                 
c$$$     *        FFTW_FORWARD,  FFTW_MEASURE  ) 
c$$$      CALL  FFTWND_F77_CREATE_PLAN(FFTW_Z_TO_F_PLAN, 1, NZ,
c$$$     *        FFTW_FORWARD,  FFTW_MEASURE + FFTW_IN_PLACE )
c$$$      CALL RFFTWND_F77_CREATE_PLAN(FFTW_X_TO_P_PLAN, 1, NX,                                 
c$$$     *        FFTW_BACKWARD,  FFTW_MEASURE  ) 
c$$$      CALL  FFTWND_F77_CREATE_PLAN(FFTW_Z_TO_P_PLAN, 1, NZ,
c$$$     *        FFTW_BACKWARD,  FFTW_MEASURE + FFTW_IN_PLACE )


c     ------------------------------
c     Define datatypes
c     ------------------------------

      ! Box full x to z (1)
      call MPI_TYPE_VECTOR(NZP,NXP,NX/2+1,
     &     MPI_DOUBLE_COMPLEX,TYPE1,ierror) 
      call MPI_TYPE_COMMIT(TYPE1,ierror)                                              

      bl(1:2)=(/1, 1/)
      disp(1:2)= (/0, NXP*16/)  
      types=(/TYPE1, MPI_UB/)

      call MPI_TYPE_STRUCT(2,bl,disp,types,XY2ZY_1,ierror)
      call MPI_TYPE_COMMIT(XY2ZY_1,ierror)
      call MPI_TYPE_FREE(TYPE1,ierror)

      ! Box full x to z (2)
      call MPI_TYPE_VECTOR(NZP,NXP,NXP+1,
     &     MPI_DOUBLE_COMPLEX,TYPE1,ierror) 
      call MPI_TYPE_COMMIT(TYPE1,ierror)                                              

      bl(1:2)=(/1, 1/)
      disp(1:2)= (/0, NZP*(NXP+1)*16/)  
      types=(/TYPE1, MPI_UB/)

      call MPI_TYPE_STRUCT(2,bl,disp,types,XY2ZY_2,ierror)
      call MPI_TYPE_COMMIT(XY2ZY_2,ierror)
      call MPI_TYPE_FREE(TYPE1,ierror)

c     ///////////////////////////////////
c     OTHER POSSIBLE TRANSPOSES!!!!
c     ///////////////////////////////////
c$$$      ! Box full x to z (2)
c$$$      call MPI_TYPE_VECTOR(NXP,1,NZ,
c$$$     &     MPI_DOUBLE_COMPLEX,TYPE1,ierror) 
c$$$      call MPI_TYPE_COMMIT(TYPE1,ierror)                                              
c$$$
c$$$      bl(1:2)=(/1, 1/)
c$$$      disp(1:2)= (/0, 16/)  
c$$$      types=(/TYPE1, MPI_UB/)
c$$$
c$$$      call MPI_TYPE_STRUCT(2,bl,disp,types,XY2ZY_2,ierror)
c$$$      call MPI_TYPE_COMMIT(XY2ZY_2,ierror)

c$$$      call MPI_TYPE_VECTOR(NZP,NXP,NXP+1,
c$$$     &     MPI_DOUBLE_COMPLEX,XY2ZY_2,ierror) 
c$$$      call MPI_TYPE_COMMIT(XY2ZY_2,ierror)                                     
c$$$      ! Box full x to z
c$$$      call MPI_TYPE_VECTOR(NZ,NXP,NX/2+1,
c$$$     &     MPI_DOUBLE_COMPLEX,XY2ZY_1,ierr) 
c$$$      call MPI_TYPE_COMMIT(XY2ZY_1,ierr)                                              
c$$$
c$$$      ! Box full z to x
c$$$      call MPI_TYPE_VECTOR(NX/2,NX/2,NX/2+1,
c$$$     &     MPI_DOUBLE_COMPLEX,XY2ZY,ierr) 
c$$$      call MPI_TYPE_COMMIT(XY2ZY,ierr)                                      
         


c$$$c
c$$$c     CHECK POINTS
c$$$c
c$$$      FNAME='start.h5'
c$$$      call readHDF5(FNAME)
c$$$
c$$$!      tvar=U1
c$$$!      FNAME='test.h5'
c$$$!      call WriteHDF5_var_real(FNAME)
c$$$
c$$$      FNAME='out.h5'
c$$$      call writeHDF5(FNAME)  
c$$$
c$$$      call mpi_finalize(ierror)
c$$$      stop 
c$$$     

C Set a string to determine which input/output files to use
C When MPI is used, each process will read/write to files with the
C number of their rank (+1) appended to the end of the file.
C The string MPI_IO_NUM will be used to define the RANK+1 for each process
        IF (NPROCY.lt.10) THEN
          MPI_IO_NUM=CHAR(MOD(RANKY+1,10)+48)
        ELSE IF (NPROCY.lt.100) THEN
          MPI_IO_NUM=CHAR(MOD(RANKY+1,100)/10+48)
     &             //CHAR(MOD(RANKY+1,10)+48)
        ELSE IF (NPROCY.lt.1000) THEN
          MPI_IO_NUM=CHAR(MOD(RANKY+1,1000)/100+48)
     &             //CHAR(MOD(RANKY+1,100)/10+48)
     &             //CHAR(MOD(RANKY+1,10)+48)
        ELSE IF (NPROCY.lt.10000) THEN
          MPI_IO_NUM=CHAR(MOD(RANKY+1,10000)/1000+48)
     &             //CHAR(MOD(RANKY+1,1000)/100+48)
     &             //CHAR(MOD(RANKY+1,100)/10+48)
     &             //CHAR(MOD(RANKY+1,10)+48)
        ELSE
           WRITE(6,*) 'ERROR, NPROCS>10,000, Unsupported problem size'
        END IF

      RETURN
      END


!----*|--.---------.---------.---------.---------.---------.---------.-|------
      SUBROUTINE INIT_CHAN_MPI
C----*|--.---------.---------.---------.---------.---------.---------.-|-----|
C Initialize any constants here
      INCLUDE 'header'

      INTEGER J, N

      IF (RANK.EQ.0) 
     &     write(*,*) '*******IN INIT_CHAN_MPI *********'

      PI=4.D0*ATAN(1.D0)

! Set the upper and lower bounds for timestepping
      IF (RANKY.eq.0) THEN
        write(*,*) 'U_BC_YMIN: ',U_BC_YMIN
        JEND=NY
        IF (U_BC_YMIN.EQ.0) THEN
          JSTART=2
        ELSE IF (U_BC_YMIN.EQ.1) THEN
          JSTART=1
        ELSE
          JSTART=2
        END IF
! Now, set the indexing for the scalar equations
        DO N=1,N_TH
          JEND_TH(N)=NY
          IF (TH_BC_YMIN(N).EQ.0) THEN
            JSTART_TH(N)=2
          ELSE IF (TH_BC_YMIN(N).EQ.1) THEN
            JSTART_TH(N)=1
          ELSE
            JSTART_TH(N)=2
          END IF
        END DO
      ELSE IF (RANKY.eq.NPROCY-1) THEN
        write(*,*) 'U_BC_YMAX: ',U_BC_YMAX
        JSTART=2
        IF (U_BC_YMAX.EQ.0) THEN
          JEND=NY-1
        ELSE IF (U_BC_YMAX.EQ.1) THEN
          JEND=NY
        ELSE
          JEND=NY-1
        END IF

! Set the upper and lower limits of timestepping of the scalar equations
        DO N=1,N_TH
        JSTART_TH(N)=2
        IF (TH_BC_YMAX(N).EQ.0) THEN
          JEND_TH(N)=NY-1
        ELSE IF (TH_BC_YMAX(N).EQ.1) THEN
          JEND_TH(N)=NY
        ELSE
          JEND_TH(N)=NY-1
        END IF
        END DO

      ELSE
! Here, we are on a middle process
        JSTART=2
        JEND=NY
        DO N=1,N_TH
           JSTART_TH(N)=2
           JEND_TH(N)=NY
        END DO
      END IF

      RETURN
      END

      SUBROUTINE APPLY_BC_TH_MPI(MATL,MATD,MATU,VEC,N)
! This subroutine applies the boundary conditions to the
! scalar fields prior to the implicit solve
      INCLUDE 'header'

      INTEGER N
! We first need to check to see which processor we are, if we are
! the upper or lowermost process, then apply boundary conditions
        IF (RANKY.eq.0) THEN
! If we have the lowest plane, apply the boundary conditions
          CALL APPLY_BC_TH_LOWER(MATL,MATD,MATU,VEC,N)
        END IF
        IF (RANKY.eq.NPROCY-1) THEN
! If we have the upper plane, apply the boundary conditions
          CALL APPLY_BC_TH_UPPER(MATL,MATD,MATU,VEC,N)
        END IF
      RETURN
      END

      SUBROUTINE APPLY_BC_TH_MPI_C(MATL_C,MATD_C,MATU_C,VEC_C,N)
! This subroutine applies the boundary conditions to the
! scalar fields prior to the implicit solve
      INCLUDE 'header'

      INTEGER N
! We first need to check to see which processor we are, if we are
! the upper or lowermost process, then apply boundary conditions
        IF (RANK.eq.0) THEN
! If we have the lowest plane, apply the boundary conditions
          CALL APPLY_BC_TH_LOWER_C(MATL_C,MATD_C,MATU_C,VEC_C,N)
        END IF
        IF (RANK.eq.NPROCS-1) THEN
! If we have the upper plane, apply the boundary conditions
          CALL APPLY_BC_TH_UPPER_C(MATL_C,MATD_C,MATU_C,VEC_C,N)
        END IF
      RETURN
      END


      SUBROUTINE APPLY_BC_U2_MPI(MATL,MATD,MATU,VEC)
! This subroutine applies the boundary conditions to the
! velocity field prior to the implicit solve
! Note, MATL, MATD, etc. are dimensioned in header
      INCLUDE 'header'

! We first need to check to see which processor we are, if we are
! the upper or lowermost process, then apply boundary conditions
        IF (RANKY.eq.0) THEN
! If we have the lowest plane, apply the boundary conditions
          CALL APPLY_BC_2_LOWER(MATL,MATD,MATU,VEC)
        END IF
        IF (RANKY.eq.NPROCY-1) THEN
! If we have the highest plane, apply the boundary conditions 
          CALL APPLY_BC_2_UPPER(MATL,MATD,MATU,VEC)
        END IF
      RETURN
      END

      SUBROUTINE APPLY_BC_U2_MPI_C(MATL_C,MATD_C,MATU_C,VEC_C)
! This subroutine applies the boundary conditions to the
! velocity field prior to the implicit solve
! Note, MATL, MATD, etc. are dimensioned in header
      INCLUDE 'header'

! We first need to check to see which processor we are, if we are
! the upper or lowermost process, then apply boundary conditions
        IF (RANK.eq.0) THEN
! If we have the lowest plane, apply the boundary conditions
          CALL APPLY_BC_2_LOWER_C(MATL_C,MATD_C,MATU_C,VEC_C)
        END IF
        IF (RANK.eq.NPROCS-1) THEN
! If we have the highest plane, apply the boundary conditions 
          CALL APPLY_BC_2_UPPER_C(MATL_C,MATD_C,MATU_C,VEC_C)
        END IF
      RETURN
      END


      SUBROUTINE APPLY_BC_U1_MPI(MATL,MATD,MATU,VEC)
! This subroutine applies the boundary conditions to the
! velocity field prior to the implicit solve
! Note, MATL, MATD, etc. are dimensioned in header
      INCLUDE 'header'

! We first need to check to see which processor we are, if we are
! the upper or lowermost process, then apply boundary conditions
        IF (RANKY.eq.0) THEN
! If we have the lowest plane, apply the boundary conditions
          CALL APPLY_BC_1_LOWER(MATL,MATD,MATU,VEC)        
        END IF
        IF (RANKY.eq.NPROCY-1) THEN
! If we have the highest plane, apply the boundary conditions 
          CALL APPLY_BC_1_UPPER(MATL,MATD,MATU,VEC)
        END IF
      RETURN
      END

      SUBROUTINE APPLY_BC_U1_MPI_C(MATL_C,MATD_C,MATU_C,VEC_C)
! This subroutine applies the boundary conditions to the
! velocity field prior to the implicit solve
! Note, MATL, MATD, etc. are dimensioned in header
      INCLUDE 'header'

! We first need to check to see which processor we are, if we are
! the upper or lowermost process, then apply boundary conditions
        IF (RANK.eq.0) THEN
! If we have the lowest plane, apply the boundary conditions
          CALL APPLY_BC_1_LOWER_C(MATL_C,MATD_C,MATU_C,VEC_C)
        END IF
        IF (RANK.eq.NPROCS-1) THEN
! If we have the highest plane, apply the boundary conditions 
          CALL APPLY_BC_1_UPPER_C(MATL_C,MATD_C,MATU_C,VEC_C)
        END IF
      RETURN
      END


      SUBROUTINE APPLY_BC_U3_MPI(MATL,MATD,MATU,VEC)
! This subroutine applies the boundary conditions to the
! velocity field prior to the implicit solve
! Note, MATL, MATD, etc. are dimensioned in header
      INCLUDE 'header'

! We first need to check to see which processor we are, if we are
! the upper or lowermost process, then apply boundary conditions
        IF (RANKY.eq.0) THEN
! If we have the lowest plane, apply the boundary conditions
          CALL APPLY_BC_3_LOWER(MATL,MATD,MATU,VEC)
        END IF
        IF (RANKY.eq.NPROCY-1) THEN
! If we have the highest plane, apply the boundary conditions 
          CALL APPLY_BC_3_UPPER(MATL,MATD,MATU,VEC)
        END IF
      RETURN
      END

      SUBROUTINE APPLY_BC_U3_MPI_C(MATL_C,MATD_C,MATU_C,VEC_C)
! This subroutine applies the boundary conditions to the
! velocity field prior to the implicit solve
! Note, MATL, MATD, etc. are dimensioned in header
      INCLUDE 'header'

! We first need to check to see which processor we are, if we are
! the upper or lowermost process, then apply boundary conditions
        IF (RANK.eq.0) THEN
! If we have the lowest plane, apply the boundary conditions
          CALL APPLY_BC_3_LOWER_C(MATL_C,MATD_C,MATU_C,VEC_C)
        END IF
        IF (RANK.eq.NPROCS-1) THEN
! If we have the highest plane, apply the boundary conditions 
          CALL APPLY_BC_3_UPPER_C(MATL_C,MATD_C,MATU_C,VEC_C)
        END IF
      RETURN
      END


      SUBROUTINE APPLY_BC_REM_DIV_MPI(MATL_C,MATD_C,MATU_C,VEC_C,K)
! This subroutine applies the boundary conditions for the Poisson Eq.
! Note, MATL, MATD, etc. are dimensioned in header
      INCLUDE 'header'

      INTEGER I,K
! We first need to check to see which processor we are, if we are
! the upper or lowermost process, then apply boundary conditions
C Apply the boundary conditions
        IF (RANKY.EQ.0) THEN
        DO I=0,NXP-1
C Use homogeneous dirichlet BCS for kx=kz=0 component at bottom wall
          IF ((K.EQ.0).AND.(I.EQ.0).AND.(RANKZ.EQ.0)) THEN
C Otherwise the matrix will be singular
C Use homogeneous dirichlet BCS for kx=kz=0 component at bottom wall
            MATL_C(I,1)=0.
            MATD_C(I,1)=1.
            MATU_C(I,1)=0.
            VEC_C(I,1)=(0.,0.)
          ELSE
C Use Dirichlet boundary conditions, dp/dz=0 at walls
            MATL_C(I,1)=0.
            MATD_C(I,1)=1.
            MATU_C(I,1)=-1.
            VEC_C(I,1)=(0.,0.)
          END IF
        END DO
        END IF
C Apply the boundary conditions
        IF (RANKY.EQ.NPROCY-1) THEN
          DO I=0,NXP-1
            MATL_C(I,NY)=1.
            MATD_C(I,NY)=-1.
            MATU_C(I,NY)=0.
            VEC_C(I,NY)=(0.,0.)
          END DO
        END IF

        RETURN
        END

      SUBROUTINE APPLY_BC_POISSON_MPI(MATL_C,MATD_C,MATU_C,VEC_C,K)
! This subroutine applies the boundary conditions for the Poisson Eq.
! Note, MATL, MATD, etc. are dimensioned in header
      INCLUDE 'header'

      INTEGER I,K
! We first need to check to see which processor we are, if we are
! the upper or lowermost process, then apply boundary conditions
C Use dirichlet boundary condition at the lower wall to
C prevent the tridiagonal matrix from becomming singular for i,k=0
        IF (RANKY.eq.0) THEN
        DO I=0,NXP-1
          IF ((I.EQ.0).AND.(K.EQ.0).AND.(RANKZ.EQ.0)) THEN
            MATD_C(I,1)=1.
            MATU_C(I,1)=0.
            VEC_C(I,1)=(0.,0.)
          ELSE
! Here, apply Neumann boundary conditions (dp/dz=0) at the walls
            MATD_C(I,1)=1.
            MATU_C(I,1)=-1.
            VEC_C(I,1)=(0.,0.)
          END IF
        END DO
        END IF
C Use dirichlet boundary condition at the lower wall to
C prevent the tridiagonal matrix from becomming singular for i,k=0
        IF (RANKY.eq.NPROCY-1) THEN
        DO I=0,NXP-1
            MATD_C(I,NY)=-1.
            MATL_C(I,NY)=1.
            VEC_C(I,NY)=(0.,0.)
        END DO
        END IF

      RETURN
      END

      SUBROUTINE APPLY_BC_VEL_MPI
! This subroutine applies the boundary conditions for the Poisson Eq.
! Note, MATL, MATD, etc. are dimensioned in header
      INCLUDE 'header'

! Apply Boundary conditions to velocity field
      IF (RANKY.EQ.0) THEN
        CALL APPLY_BC_VEL_LOWER
      END IF
      IF (RANKY.EQ.NPROCY-1) THEN
        CALL APPLY_BC_VEL_UPPER
      END IF
 
      RETURN
      END


      SUBROUTINE APPLY_BC_VEL_PHYS_MPI
! This subroutine applies the boundary conditions for the Poisson Eq.
! Note, MATL, MATD, etc. are dimensioned in header
      INCLUDE 'header'

! Apply Boundary conditions to velocity field
      IF (RANKY.EQ.0) THEN
        CALL APPLY_BC_VEL_PHYS_LOWER
      END IF
      IF (RANKY.EQ.NPROCY-1) THEN
        CALL APPLY_BC_VEL_PHYS_UPPER
      END IF

      RETURN
      END

 
      SUBROUTINE TRANSPOSE_MPI_XZ_TO_XY()
! This subroutine starts with all arrays decomposed in x-z slabs
! and transposes the data so that it is decomposed in x-y slabs
! x-y slabs.
c$$$      include 'header'
c$$$
c$$$      integer i,j,k,N
c$$$
c$$$      real*8 A(0:NX+1,0:NZ+1,0:NY+1) 
c$$$      real*8 buffer(0:NX+1,0:NY+1,0:NZ+1)
c$$$
c$$$      real*8 test1(1:NX,1:NY,1:NZ)
c$$$      real*8 test2(1:NX,1:NY,1:NZ)
c$$$
c$$$       do i=0,NXM
c$$$       do j=1,NY
c$$$         write(100+RANK,*) 'i,j,U1: ',i,j,A(I,0,J)
c$$$       end do
c$$$       end do
c$$$
c$$$      do k=0,NZ-1
c$$$        do j=1,NY
c$$$          do i=0,NX-1
c$$$            buffer(i,j,k)=A(i,k,j)
c$$$            test1(i+1,j,k+1)=A(i,k,j)
c$$$          end do
c$$$        end do
c$$$      end do
c$$$
c$$$      CALL MPI_AllToAll(test1,(NX)*(NY)*(NZ)/NPROCS
c$$$     &   ,MPI_DOUBLE_PRECISION,test2,(NX)*(NY)*(NZ)/NPROCS
c$$$     &       ,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,IERROR)
c$$$      CALL MPI_Barrier(MPI_COMM_WORLD)
c$$$
c$$$! A should now be indexed by A(NX,NZ/nprocs,NY*nprocs)
c$$$
c$$$         n=1
c$$$         do k=0,NZ/NPROCS-1
c$$$         do i=0,NX-1 
c$$$         do j=1,NY
c$$$           A(i,k,j)=test2(i+1,j,k+1)
c$$$         end do
c$$$         end do
c$$$         end do
c$$$  
c$$$         do n=2,NPROCS
c$$$         do k=0,NZ/NPROCS-1
c$$$         do i=0,NX-1 
c$$$         do j=2,NY
c$$$           A(i,k,NY+j-1)=test2(i+1,NY+j,k+1)
c$$$         end do
c$$$         end do
c$$$         end do
c$$$         end do
c$$$
c$$$         NX_t=NX
c$$$         NZ_t=NZ/NPROCS
c$$$         NY_t=NY*NPROCS-(NPROCS-1)
c$$$ 
c$$$! A should now be indexed from A(0:NX_t-1,0:NZ_t/NPROCS-1,1:NY*NPROCS-(NPROCS-1))
c$$$! (The new wall locations are 1 and NY-(NPROCS-1) )
c$$$
c$$$          do i=0,NXM
c$$$          do j=1,NY_T
c$$$            write(200+RANK,*) 'i,j,U1: ',i,j,A(i,0,J)
c$$$          end do
c$$$          end do
c$$$

   
      RETURN
      END

      


      SUBROUTINE get_minimum_mpi(val)
      
      include 'header'

      REAL*8 val,vmin

      call MPI_ALLREDUCE(val,vmin,1,MPI_DOUBLE_PRECISION,
     &        MPI_MIN,MPI_COMM_WORLD,ierror)
      
      val=vmin

      RETURN

      END 



      SUBROUTINE FFT_XZ_TO_FOURIER(V,VV,JMIN,JMAX)

      include 'header'

      REAL*8     V  (0:NX+1,0:NZP+1,0:NY+1)
      COMPLEX*16 VV (0:NXP ,0:NZ +1,0:NY+1)
      COMPLEX*16 TMP(0:NX/2,0:NZP+1,0:NY+1)

      INTEGER I,J,K
      INTEGER JMIN,JMAX 

      !write(100+RANK,'(1E25.15)') V(0:NX-1,0:NZP-1,1)
      !write(100+RANK) V(0:NX-1,0:NZP-1,1)

      ! FFT in X
      DO J=JMIN,JMAX
       CALL RFFTWND_F77_REAL_TO_COMPLEX(FFTW_X_TO_F_PLAN,NZP,
     *    V(0,0,J), 1, NX+2, TMP(0,0,J), 1, NX/2+1)
        DO K=0,NZP-1
          DO I=0,NKX
            TMP(I,K,J)=TMP(I,K,J)/NX
          END DO
          DO I=NKX+1,NX/2
            TMP(I,K,J)=cmplx(0.d0,0.d0)
          END DO
        END DO
      END DO

      !write(110+RANK,'(2E25.15)') TMP(0:NX/2,0:NZP-1,1)
      !write(110+RANK) TMP(0:NX/2,0:NZP-1,1)

      DO J=JMIN,JMAX
         call mpi_alltoall(TMP(0,0,J),1,XY2ZY_1,
     &        VV(0,0,J),1,XY2ZY_2,MPI_COMM_Z,IERROR)
      END DO

      !write(120+RANK,'(2E25.15)') VV(0:NXP-1,0:NZ-1,1)
      !write(120+RANK) VV(0:NXP-1,0:NZ-1,1)
      
      ! FFT in Z
      DO J=JMIN,JMAX
         CALL FFTWND_F77(FFTW_Z_TO_F_PLAN,NXP,
     *        VV(0,0,J), NXP+1, 1, VV(0,0,J), NXP+1, 1)
         DO K=0,NKZ
            DO I=0,NXP
               VV(I,K,J)=VV(I,K,J)/NZ
            END DO
         END DO
c$$$         DO K=NKZ+1,NZ-NKZ-1
c$$$            DO I=0,NXP
c$$$               VV(I,K,J)=cmplx(0.d0,0.d0)
c$$$            END DO
c$$$         END DO
c$$$         DO K=NZ-NKZ,NZ-1
c$$$            DO I=0,NXP
c$$$               VV(I,K,J)=VV(I,K,J)/NZ
c$$$            END DO
c$$$         END DO
         DO K=1,NKZ
            DO I=0,NXP
               VV(I,NKZ+K,J)=VV(I,NZ-1+K-NKZ,J)/NZ
            END DO
         END DO
      END DO

C$$$      ! write(130+RANK,'(2E25.15)') VV(0:NXP-1,0:NZ-1,1)
      !write(130+RANK) VV(0:NXP-1,0:NZ-1,1)

      END SUBROUTINE



      SUBROUTINE FFT_XZ_TO_PHYSICAL(VV,V,JMIN,JMAX) 

      include 'header'

      REAL*8     V  (0:NX+1,0:NZP+1,0:NY+1)
      COMPLEX*16 VV (0:NXP ,0:NZ +1,0:NY+1)
      COMPLEX*16 TMP(0:NX/2,0:NZP+1,0:NY+1)

      INTEGER I,J,K
      INTEGER JMIN,JMAX

      ! FFT in Z
      DO J=JMIN,JMAX
         ! UNPACK
         DO K=NKZ,1,-1
            DO I=0,NXP
               VV(I,NZ-1+K-NKZ,J)=VV(I,NKZ+K,J)
            END DO
         END DO
         DO I=0,NXP
            DO K=NKZ+1,NZ-NKZ-1
               VV(I,K,J)=cmplx(0.d0,0.d0)
            END DO
!            DO K=NZ,NZ+1
!               VV(I,K,J)=cmplx(0.d0,0.d0)
!            END DO
         END DO
         CALL FFTWND_F77(FFTW_Z_TO_P_PLAN,NXP,
     *        VV(0,0,J), NXP+1, 1, VV(0,0,J), NXP+1, 1)
      END DO

      DO J=JMIN,JMAX
         call mpi_alltoall(VV(0,0,J),1,XY2ZY_2,
     &        TMP(0,0,J),1,XY2ZY_1,MPI_COMM_Z,IERROR)
      END DO

      ! FFT in X
      DO J=JMIN,JMAX
       CALL RFFTWND_F77_COMPLEX_TO_REAL(FFTW_X_TO_P_PLAN,NZP,
     *    TMP(0,0,J), 1, NX/2+1, V(0,0,J), 1, NX+2)
      END DO


      END SUBROUTINE


      SUBROUTINE INTEGRATE_Y_VAR(VAR,RES,COMM)

      INCLUDE 'header'

      INTEGER J,COMM
      REAL*8 VAR(0:NY+1),RES

! Integrat the instantaneous mean profile numerically at GY points
      RES=0.
      IF (USE_MPI) THEN
         do j=2,NY
            RES=RES+0.5*(VAR(j)+VAR(j-1))*DY(j)
         end do
         call MPI_ALLREDUCE(MPI_IN_PLACE,RES,1,
     &        MPI_DOUBLE_PRECISION,MPI_SUM,COMM,ierror)
      ELSE
         do j=2,NY
            RES=RES+0.5*(VAR(j)+VAR(j-1))*DY(j)
         end do
      END IF
      RES=RES/LY      

      END SUBROUTINE


      SUBROUTINE END_RUN_MPI(FLAG)
      
      INCLUDE 'header'

      LOGICAL FLAG

      IF (RANK.EQ.0) THEN
         CALL END_RUN(FLAG)
      END IF
      CALL MPI_BCAST(FLAG,1,MPI_LOGICAL,0,MPI_COMM_WORLD,IERROR)
      
      END

