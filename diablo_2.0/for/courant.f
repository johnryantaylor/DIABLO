      subroutine courant
! This subroutine sets the timestep based on the specified CFL number
! The subroutine should be called with the velocity in physical space

      include 'header'

      real*8 vel
      real*8 dt
      real*8 dt_x,dt_y,dt_z
      real*8 Nmax
      integer i,j,k,n
      integer imin,jmin,kmin

! Set the initial dt to some arbitrary large number
      dt=999.d0

! Set the timestep based on viscosity and diffusivity
      dt=min(dt,0.5d0*min(dx(1),dy(1))/NU)
      do n=1,N_TH
        dt=min(dt,dt*NU/(NU/PR(n)))
      end do
! Make sure that we capture the inertial period (for rotating flows)
      if (I_RO.ne.0.d0) then
        dt=min(dt,2.d0*PI/abs(I_RO)/20.d0)
      end if
! Make sure that we capture the buoyancy period (for stratified flows)
      do n=1,N_TH
      if (RI(n).ne.0.d0) then
        Nmax=sqrt(abs(RI(N)*TH_BC_YMIN_C1(N)))
        dt=min(dt,0.1*2.d0*PI/Nmax)
      end if
      end do

! Use the model velocity to calculate the CFL number
      do j=1,NY
        do k=0,NZP-1
          do i=0,NXM
            dt_x=cfl*dx(i)/abs(U1(i,k,j))
            dt_y=cfl*dy(j)/abs(U2(i,k,j))
            dt_z=cfl*dz(k)/abs(U3(i,k,j))
            dt=min(dt,dt_x,dt_y,dt_z)
          end do
        end do
      end do
      if (USE_MPI) then
         call get_minimum_mpi(dt)
      end if

      if (dt.le.0) then
        IF (RANK.EQ.0)
     &        write(*,*) 'Error: dt<=0 in courant'
! Set DELTA_T to some small default value
        DELTA_T=0.0001d0
      else if (dt.ge.999.) then
!        write(*,*) 'WARNING: DELTA_T > 999, value capped at 999'
        DELTA_T=999.d0
      else
        DELTA_T=dt
      end if
      H_BAR(1)=DELTA_T*(8.d0/15.d0)
      H_BAR(2)=DELTA_T*(2.d0/15.d0)
      H_BAR(3)=DELTA_T*(5.d0/15.d0)

      return
      end

