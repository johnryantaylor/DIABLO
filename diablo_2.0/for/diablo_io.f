C This file contains subroutines for inputting and outputting data in
C Diablo as well as all subroutines called directly from diablo.f
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE INITIALIZE
      INCLUDE 'header'

      REAL    VERSION, CURRENT_VERSION
      logical RESET_TIME
      INTEGER I, J, K, N

      OPEN (11,file='input.dat',form='formatted',status='old')      

C Read input file.
C   (Note - if you change the following section of code, update the
C    CURRENT_VERSION number to make obsolete previous input files!)

      CURRENT_VERSION=2.0
      READ(11,*)
      READ(11,*)
      READ(11,*)
      READ(11,*)
      READ(11,*) FLAVOR,   VERSION
      IF (VERSION .NE. CURRENT_VERSION) STOP 'Wrong input data format.'
      READ(11,*)
      READ(11,*) USE_MPI,    LES
      READ(11,*)
      READ(11,*) NU, BETA, LX, LY, LZ
      READ(11,*)
      READ(11,*) NU_V_SCALE
      READ(11,*)
      READ(11,*) NUM_PER_DIR, CREATE_NEW_FLOW
      READ(11,*)
      READ(11,*) N_TIME_STEPS, TIME_LIMIT, DELTA_T, RESET_TIME, 
     &     VARIABLE_DT, CFL, UPDATE_DT
      READ(11,*)
      READ(11,*) VERBOSITY, SAVE_FLOW_INT, SAVE_STATS_INT, MOVIE
      READ(11,*)
! Read in the parameters for the N_TH scalars
      DO N=1,N_TH
        READ(11,*)
        READ(11,*) CREATE_NEW_TH(N)
        READ(11,*)
        READ(11,*) FILTER_TH(N), FILTER_INT(N)
        READ(11,*)
        READ(11,*) RI(N), PR(N)
      END DO

C If we are using MPI, then Initialize the MPI Variables
      IF (USE_MPI) THEN
        CALL INIT_MPI
      ELSE
         RANK=0
         RANKY=0
         RANKZ=0
      END IF

C Check compatibility
      IF ((USE_MPI).AND.(NUM_PER_DIR.eq.3)) THEN
        WRITE(6,*) 'ERROR: periodic.f isnt parallelized in MPI'
        pause
      END IF 

      IF (RANK.EQ.0) THEN
         WRITE(6,*) 
         WRITE(6,*) '             ****** WELCOME TO DIABLO ******'
         WRITE(6,*)
      END IF

C Initialize case-specific packages.
      IF (NUM_PER_DIR.EQ.3) THEN
        CALL INPUT_PER
        CALL CREATE_GRID_PER
        CALL INIT_PER
      ELSEIF (NUM_PER_DIR.EQ.2) THEN
        CALL INPUT_CHAN
        CALL CREATE_GRID_CHAN
        IF (USE_MPI) THEN
          CALL INIT_CHAN_MPI
        ELSE 
          CALL INIT_CHAN
        END IF
        IF (MOVIE) THEN
          CALL INIT_CHAN_MOVIE
        END IF
      ELSEIF (NUM_PER_DIR.EQ.1) THEN
        CALL INPUT_DUCT
        CALL CREATE_GRID_DUCT
        CALL INIT_DUCT
      ELSEIF (NUM_PER_DIR.EQ.0) THEN
        CALL INPUT_CAV 
        CALL CREATE_GRID_CAV
        CALL INIT_CAV
      END IF

C Initialize grid
         IF (RANK.EQ.0) THEN
      WRITE(6,*) 'Note that this code is distributed under the ',
     *           'GNU General Public License.'
      WRITE(6,*) 'No warranty is expressed or implied.'
      WRITE(6,*)
      write(*,*) 'Flavor: ',FLAVOR
      WRITE(6,*) 'Grid size: NX =',NX,', NY =',NY,', NZ =',NZ,'.'
      DO N=1,N_TH
        WRITE(6,*) 'Scalar number: ',N
        WRITE(6,*) '  Richardson number: ',RI(N)
        WRITE(6,*) '  Prandlt number: ',PR(N)
      END DO
      WRITE(6,*) 'NU: ',NU
      WRITE(6,*) 'BETA: ',BETA
      END IF

C Initialize FFT package (includes defining the wavenumber vectors).
      CALL INIT_FFT

C Initialize RKW3 parameters.
      H_BAR(1)=DELTA_T*(8.0/15.0)
      H_BAR(2)=DELTA_T*(2.0/15.0)
      H_BAR(3)=DELTA_T*(5.0/15.0)
      BETA_BAR(1)=1.0
      BETA_BAR(2)=25.0/8.0
      BETA_BAR(3)=9.0/4.0
      ZETA_BAR(1)=0.0
      ZETA_BAR(2)=-17.0/8.0
      ZETA_BAR(3)=-5.0/4.0

C Initialize values for reading of scalars
      NUM_READ_TH=0
      DO N=1,N_TH
        IF (CREATE_NEW_TH(N)) THEN
          NUM_READ_TH=NUM_READ_TH 
        ELSE
          NUM_READ_TH=NUM_READ_TH + 1
          READ_TH_INDEX(NUM_READ_TH)=N
        END IF
      END DO
      IF (NUM_PER_DIR.EQ.2) THEN
        CALL CREATE_TH_CHAN 
      ELSE IF (NUM_PER_DIR.EQ.3) THEN
        CALL CREATE_TH_PER
      END IF 

C Create flow.
      IF (CREATE_NEW_FLOW) THEN
        IF (NUM_PER_DIR.EQ.3) THEN
          CALL CREATE_FLOW_PER
        ELSEIF (NUM_PER_DIR.EQ.2) THEN
          CALL CREATE_FLOW_CHAN
        ELSEIF (NUM_PER_DIR.EQ.1) THEN
          CALL CREATE_FLOW_DUCT
        ELSEIF (NUM_PER_DIR.EQ.0) THEN
          CALL CREATE_FLOW_CAV
        END IF
        IF (RANK.EQ.0) 
     &          write(*,*) 'A new flowfield has been created'
      ELSE
        IF (RANK.EQ.0) 
     &        write(*,*) 'Reading flow...'
        CALL READ_FLOW
        IF (RANK.EQ.0) 
     &       write(*,*) 'Done reading flow'
        
        IF (USE_MPI) THEN
           CALL GHOST_CHAN_MPI
        END IF

C Initialize flow.
        IF (RESET_TIME .OR. CREATE_NEW_FLOW) THEN
           PREVIOUS_TIME_STEP=0
           TIME_STEP=0
           TIME=0
        END IF

        CALL SAVE_STATS(.FALSE.)
      end if

      RETURN
      END



C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE SAVE_STATS(FINAL)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      INCLUDE 'header'
      LOGICAL FINAL

      IF (NUM_PER_DIR.EQ.3) THEN
        CALL SAVE_STATS_PER(FINAL)          
      ELSEIF (NUM_PER_DIR.EQ.2) THEN
        CALL SAVE_STATS_CHAN(FINAL)          
      ELSEIF (NUM_PER_DIR.EQ.1) THEN
        CALL SAVE_STATS_DUCT(FINAL)          
      ELSEIF (NUM_PER_DIR.EQ.0) THEN
        CALL SAVE_STATS_CAV(FINAL)          
      END IF

      RETURN
      END

C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE READ_FLOW
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      INCLUDE 'header'
      CHARACTER*55 FNAME
      CHARACTER*55 FNAME_TH(N_TH)
      INTEGER I, J, K, N, NUM_PER_DIR_T
      
!      FNAME='diablo.start'
      FNAME='start.h5'
      IF (RANK.EQ.0) 
     &     WRITE(6,*)   'Reading flow from ',FNAME
      if (FNAME(len_trim(FNAME)-2:len_trim(FNAME)).eq.".h5") then
         IF (NUM_PER_DIR.NE.2) THEN
            IF (RANK.EQ.0) THEN
            WRITE(6,*) ' READING TO HDF5 NOT IMPLEMENTED '
            WRITE(6,*) ' FOR NUM_PER_DIR DIFFERENT FROM 2'
            END IF
            STOP 
         END IF
#ifdef HDF5
         call mpi_barrier(MPI_COMM_WORLD,ierror)
         call ReadHDF5(FNAME)
#else
         IF (RANK.EQ.0) then     
         write(*,*) ' **** ERROR ******************************'
         write(*,*) ' Program not compiled with HDF5 libraries.'
         end if
         stop 
#endif
      else

         DO N=1,N_TH
            FNAME_TH(N)='diablo_th'
     &           //CHAR(MOD(N,100)/10+48)
     &           //CHAR(MOD(N,10)+48) // '.start'
         END DO

         OPEN(UNIT=10,FILE=FNAME,STATUS="OLD",FORM="UNFORMATTED")
         READ (10) NX_T, NY_T, NZ_T, NUM_PER_DIR_T, TIME, TIME_STEP


         IF (RANK.EQ.0) 
     &        write(*,*) 'NX_T, NY_T, NZ_T: ',NX_T,NY_T,NZ_T

        IF ((NX .NE. NX_T) .OR. (NY .NE. NY_T) .OR. (NZ .NE. NZ_T))
     *     STOP 'Error: old flowfield wrong dimensions. '
        IF (NUM_PER_DIR .NE. NUM_PER_DIR_T)
     *     STOP 'Error: old flowfield wrong NUM_PER_DIR. '

         IF (RANK.EQ.0) 
     &        write(*,*) 'READING FLOW'
         IF (NUM_PER_DIR.EQ.3) THEN

         READ (10) (((CU1(I,K,J),I=0,NKX),K=0,TNKZ),J=0,TNKY),
     *            (((CU2(I,K,J),I=0,NKX),K=0,TNKZ),J=0,TNKY),
     *            (((CU3(I,K,J),I=0,NKX),K=0,TNKZ),J=0,TNKY)
         DO N=1,NUM_READ_TH
!     Specify in input.dat which scalars are to be read
         OPEN(UNIT=11,FILE=FNAME_TH(READ_TH_INDEX(N)),STATUS="OLD"
     &           ,FORM="UNFORMATTED")
         READ (11) NX_T, NY_T, NZ_T, NUM_PER_DIR_T, TIME, TIME_STEP
         READ (11) (((CTH(I,K,J,READ_TH_INDEX(N))
     &           ,I=0,NKX),K=0,TNKZ),J=0,TNKY)
         CLOSE(11)
         END DO

         ELSEIF (NUM_PER_DIR.EQ.2) THEN
            READ (10) (((CU1(I,K,J),I=0,NKX),K=0,TNKZ),J=1,NY),
     *           (((CU2(I,K,J),I=0,NKX),K=0,TNKZ),J=2,NY),
     *           (((CU3(I,K,J),I=0,NKX),K=0,TNKZ),J=1,NY)

            DO N=1,NUM_READ_TH
!     Specify in input.dat which scalars are to be read
               OPEN(UNIT=11,FILE=FNAME_TH(READ_TH_INDEX(N)),STATUS="OLD"
     &              ,FORM="UNFORMATTED")
               READ (11) NX_T, NY_T, NZ_T, NUM_PER_DIR_T, 
     *              TIME, TIME_STEP
               READ (11) (((CTH(I,K,J,READ_TH_INDEX(N))
     &              ,I=0,NKX),K=0,TNKZ),J=1,NY)
               CLOSE(11)
            END DO
         ELSEIF (NUM_PER_DIR.EQ.1) THEN
            READ (10) (((CU1(I,K,J),I=0,NKX),K=1,NZ),J=1,NY),
     *           (((CU2(I,K,J),I=0,NKX),K=1,NZ),J=2,NY),
     *           (((CU3(I,K,J),I=0,NKX),K=2,NZ),J=1,NY)
         ELSEIF (NUM_PER_DIR.EQ.0) THEN
            READ (10) (((U1(I,K,J),I=2,NX),K=1,NZ),J=1,NY),
     *           (((U2(I,K,J),I=1,NX),K=1,NZ),J=2,NY),
     *           (((U3(I,K,J),I=1,NX),K=2,NZ),J=1,NY)
         END IF
         CLOSE(10)
         CLOSE(11)
      end if

C Apply initial boundary conditions, set ghost cells
      IF (USE_MPI) THEN
        call APPLY_BC_VEL_MPI
      ELSE
        call APPLY_BC_VEL_LOWER
        call APPLY_BC_VEL_UPPER
      END IF

      RETURN
      END

C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE SAVE_FLOW(FINAL)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      INCLUDE 'header'
      CHARACTER*55 FNAME
      CHARACTER*55 FNAME_TH(N_TH)
      INTEGER      I, J, K, N
      LOGICAL      FINAL,FLAGE,SAVE_PRESSURE

      SAVE_PRESSURE=.FALSE.
      if (FINAL) then 
         FNAME='end.h5'
         SAVE_PRESSURE=.TRUE.
      else
         FNAME='out.'
     &        //CHAR(MOD(TIME_STEP,1000000)/100000+48)
     &        //CHAR(MOD(TIME_STEP,100000)/10000+48)
     &        //CHAR(MOD(TIME_STEP,10000)/1000+48)
     &        //CHAR(MOD(TIME_STEP,1000)/100+48)
     &        //CHAR(MOD(TIME_STEP,100)/10+48)
     &        //CHAR(MOD(TIME_STEP,10)+48)
     &        //'.h5'
      end if
      if (FNAME(len_trim(FNAME)-2:len_trim(FNAME)).eq.".h5") then
         IF (NUM_PER_DIR.NE.2) THEN
            IF (RANK.EQ.0) THEN
            WRITE(6,*) ' SAVING TO HDF5 NOT IMPLEMENTED '
            WRITE(6,*) ' FOR NUM_PER_DIR DIFFERENT FROM 2'
            END IF
            STOP 
         END IF
#ifdef HDF5
         call mpi_barrier(MPI_COMM_WORLD,ierror)
         call WriteHDF5(FNAME,SAVE_PRESSURE)      
#else
         IF (RANK.EQ.0) then      
         write(*,*) ' **** ERROR ******************************'
         write(*,*) ' Program not compiled with HDF5 libraries.'
         end if
         stop 
#endif
      else
         IF (FINAL) THEN
            FNAME='diablo.res'
            DO N=1,N_TH
               FNAME_TH(N)='diablo_th'
     &              //CHAR(MOD(N,100)/10+48)
     &              //CHAR(MOD(N,10)+48) // '.res'
            END DO
         ELSE
            FNAME='diablo.saved'
!     &        //CHAR(MOD(TIME_STEP,10000)/1000+48)
!     &        //CHAR(MOD(TIME_STEP,1000)/100+48)
!     &        //CHAR(MOD(TIME_STEP,100)/10+48)
!     &        //CHAR(MOD(TIME_STEP,10)+48)
            DO N=1,N_TH
!     FNAME_TH(N)='/media/mybook/chemotaxis/diablo_th'
               FNAME_TH(N)='diablo_th'
     &              //CHAR(MOD(N,100)/10+48)
     &              //CHAR(MOD(N,10)+48) // '.saved'
!     &        //CHAR(MOD(TIME_STEP,10000)/1000+48)
!     &        //CHAR(MOD(TIME_STEP,1000)/100+48)
!     &        //CHAR(MOD(TIME_STEP,100)/10+48)
!     &        //CHAR(MOD(TIME_STEP,10)+48)


            END DO
         END IF
         IF (RANK.EQ.0) 
     &        WRITE(6,*) 'Writing flow to ',FNAME

         OPEN(UNIT=10,FILE=FNAME,STATUS="UNKNOWN",FORM="UNFORMATTED")
         WRITE(10) NX, NY, NZ, NUM_PER_DIR, TIME, TIME_STEP


         IF (NUM_PER_DIR.EQ.3) THEN
            WRITE(10) (((CU1(I,K,J),I=0,NKX),K=0,TNKZ),J=0,TNKY),
     *           (((CU2(I,K,J),I=0,NKX),K=0,TNKZ),J=0,TNKY),
     *           (((CU3(I,K,J),I=0,NKX),K=0,TNKZ),J=0,TNKY)
            DO N=1,N_TH
               OPEN(UNIT=11,FILE=FNAME_TH(N),STATUS="UNKNOWN"
     &              ,FORM="UNFORMATTED")
               WRITE(11) NX, NY, NZ, NUM_PER_DIR, TIME, TIME_STEP
               WRITE(11) (((CTH(I,K,J,N),I=0,NKX),K=0,TNKZ),J=0,TNKY)
               CLOSE(11)
            END DO
         ELSEIF (NUM_PER_DIR.EQ.2) THEN
            WRITE(10) (((CU1(I,K,J),I=0,NKX),K=0,TNKZ),J=1,NY),
     *           (((CU2(I,K,J),I=0,NKX),K=0,TNKZ),J=2,NY),
     *           (((CU3(I,K,J),I=0,NKX),K=0,TNKZ),J=1,NY)
            DO N=1,N_TH
               OPEN(UNIT=11,FILE=FNAME_TH(N),STATUS="UNKNOWN"
     &              ,FORM="UNFORMATTED")
               WRITE(11) NX, NY, NZ, NUM_PER_DIR, TIME, TIME_STEP
               WRITE(11) (((CTH(I,K,J,N),I=0,NKX),K=0,TNKZ),J=0,NY)
               CLOSE(11)
            END DO
         ELSEIF (NUM_PER_DIR.EQ.1) THEN
            WRITE(10) (((CU1(I,K,J),I=0,NKX),K=1,NZ),J=1,NY),
     *           (((CU2(I,K,J),I=0,NKX),K=1,NZ),J=2,NY),
     *           (((CU3(I,K,J),I=0,NKX),K=2,NZ),J=1,NY)
         ELSEIF (NUM_PER_DIR.EQ.0) THEN
            WRITE(10) (((U1(I,K,J),I=2,NX),K=1,NZ),J=1,NY),
     *           (((U2(I,K,J),I=1,NX),K=1,NZ),J=2,NY),
     *           (((U3(I,K,J),I=1,NX),K=2,NZ),J=1,NY)
         END IF
         CLOSE(10)
         CLOSE(11)
      end if
      RETURN
      END

      
      SUBROUTINE END_RUN(FLAG)
      
      INCLUDE 'header'
      
      LOGICAL FLAG,FILE_EXISTS

      FLAG=.FALSE.
      ! Check for the time
      call WALL_TIME(END_TIME)
      if (END_TIME-START_TIME.gt.TIME_LIMIT) THEN
         IF (RANK.EQ.0) 
     &        write(*,*) ' STOP because of wall-time hit!'
         FLAG=.TRUE.
      END IF
      
      INQUIRE(FILE="stop.now", EXIST=FILE_EXISTS)
      IF ( FILE_EXISTS ) THEN
         IF (RANK.EQ.0) 
     &        write(*,*) ' STOP because of stop.now file!'
         FLAG=.TRUE.
      END IF
      
      RETURN
      END 



      subroutine wall_time(wt)
c
c     Return wall-clock time as seconds after Jan. 1, 2016.
c     Support for leap year is not included anymore.
c
c     By using a 'save' statement, the wall-time after the first
c     call to the subroutine could be computed, but that is not
c     intended with the present subroutine (e.g. the history file)
c
      implicit none

      real*8 wt
      integer val(8),i,shift,day

      integer mon(12,2)
      data mon /
     &     31,28,31,30,31,30,31,31,30,31,30,31,
     &     31,29,31,30,31,30,31,31,30,31,30,31/
c
c     Get current date and time
c     val(1) : year
c     val(2) : month
c     val(3) : day
c     val(4) : difference to GMT
c     val(5) : hour
c     val(6) : minute
c     val(7) : second
c     val(8) : 1/1000 second
c
      call date_and_time(values=val)
c
c     Determine leap year
c
      if (mod(val(1),4).eq.0) then
         if (mod(val(1),100).eq.0) then
            if (mod(val(1),400).eq.0) then
               shift=2
            else
               shift=1
            end if
         else
            shift=2
         end if
      else
         shift = 1
      end if
c
c     Construct day of the year
c
      day = val(3)-1
      do i=1,val(2)-1
         day=day+mon(i,shift)
      end do
c
c     And compute wall-clock time
c
      wt = (val(1)-2016)*365*86400+
     &     day*86400+val(5)*3600+val(6)*60+val(7)+dble(val(8)/1000.d0)

      end 



