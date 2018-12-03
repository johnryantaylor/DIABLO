      SUBROUTINE WriteHDF5(FNAME,SAVE_PRESSURE)
      use hdf5

      INCLUDE 'header'

      CHARACTER*55 FNAME
      LOGICAL FINAL,SAVE_PRESSURE

      REAL*8 tmp(NX,NY,NZ)
c
c     HDF5 ------------------------------------------------------
c
c     Dataset names
      character(len=10) :: dname

c     Identifiers
      integer(hid_t) :: file_id, dset_id
      integer(hid_t) :: filspace_id, memspace_id

!     Identifiers
      integer(hid_t) :: gid, selspace_id
      integer(hid_t) :: plist_id_d

c     Dimensions in the memory and in the file
      integer(hsize_t), dimension(3) :: dimsm,dimsf

      integer(hsize_t), dimension(3) :: chunk_dims, count, offset
      integer(hsize_t), dimension(3) :: stride, block, offset_m

      integer :: rHDF5 = 3, arank = 1

      integer(hsize_t),dimension(1)       :: adims
      integer(hid_t)                      :: aid,tspace
      real   , dimension(1)               :: treal(3)
      integer, dimension(1)               :: tint(3)
      character*80                        :: namnbuf
      character*20                        :: sttimec

      integer error, ith

      double precision En(4)

      dimsm(1:3) = (/NX,NY,NZ/)
      dimsf(1:3) = (/NX,NY,NZ/)

!     Flow fields are saved in the fractional grid. We use a basic
!     interpolation
!
!     u_j+1/2=u_j+u_j+1
!
!     on the way back we just invert the relation in order to have
!     exact values.

      WRITE(6,*) 'Writing flow to ',FNAME

      chunk_dims(1) = NX
      chunk_dims(2) = 1
      chunk_dims(3) = NZ

      block(1) = NX
      block(3) = NZ

!     Stride and count for number of rows and columns in each dimension
      stride = 1
      count  = 1

!     Offset determined by the rank of a processor
      offset(1) = 0
      offset(3) = 0

      block(2) =  NY
      offset(2) = 0

!     Initialize interface
      call h5open_f(error)

!     Create the file collectively
      call h5fcreate_f(trim(FNAME), H5F_ACC_TRUNC_F,
     &                 file_id, error)

      arank=1
      adims=3
      call h5screate_simple_f(arank,adims,tspace, error)

      ! -----------------------------
      ! Resolution
      call h5acreate_f(file_id,'Resolution',H5T_STD_I32LE,tspace,
     &                 aid, error)
      tint(1)=NX
      tint(2)=NY
      tint(3)=NZ
      call h5awrite_f(aid,H5T_NATIVE_INTEGER,tint,adims,error)
      call h5aclose_f(aid, error)
      ! -----------------------------
      ! Close dataspace
      call h5sclose_f(tspace, error)

      ! -----------------------------
      ! Date
      adims=20
      call h5screate_simple_f(arank,adims,tspace, error)

      call h5acreate_f(file_id,'Date',H5T_C_S1,tspace,aid,
     &                 error)
      call time_string(sttimec)
      call h5awrite_f(aid,H5T_C_S1,sttimec,adims,error)
      call h5aclose_f(aid, error)
      call h5sclose_f(tspace,error)
      ! -----------------------------

      ! -----------------------------
      ! Extra info
      adims=80
      call h5screate_simple_f(arank,adims,tspace, error)

      call h5acreate_f(file_id,'Info',H5T_C_S1,tspace,aid,
     &                 error)
      namnbuf=' (Put here what you want) '
      call h5awrite_f(aid,H5T_C_S1,namnbuf,adims,error)
      call h5aclose_f(aid, error)
      call h5sclose_f(tspace,error)
      ! -----------------------------

      call h5screate_f(H5S_SCALAR_F,tspace,error)

c$$$      ! -----------------------------
c$$$      ! HDF5-saving version
c$$$      call h5acreate_f(file_id,'Version',H5T_STD_I32LE,tspace,aid,
c$$$     &                 error)
c$$$      call h5awrite_f(aid,H5T_NATIVE_INTEGER,hver,adims,error)
c$$$      call h5aclose_f(aid, error)
c$$$      ! -----------------------------


!     Create the timestep group
      call h5gcreate_f(file_id,'Timestep',gid, error)

      ! -----------------------------
      ! Time
      call h5acreate_f(gid,'Time',H5T_IEEE_F64LE,tspace,
     &                 aid, error)
      call h5awrite_f(aid,H5T_NATIVE_DOUBLE,TIME,adims,error)
      call h5aclose_f(aid, error)
      ! ----------------------------

      call h5sclose_f(tspace,error)

      ! Convert to physical space
      call fft_xz_to_physical(CU1,U1,0,NY+1)
      call fft_xz_to_physical(CU2,U2,0,NY+1)
      call fft_xz_to_physical(CU3,U3,0,NY+1)
      do ith=1,N_TH
         call fft_xz_to_physical(CTH(0,0,0,ith),TH(0,0,0,ith),0,NY+1)
      end do

!     Create property list for the chunked dataset creation
      call h5pcreate_f(H5P_DATASET_CREATE_F, plist_id_d, error)
      call h5pset_chunk_f(plist_id_d, rHDF5, chunk_dims, error)

!     Create the dataspace for ur
      call h5screate_simple_f(rHDF5, dimsf, filspace_id, error)
      call h5screate_simple_f(rHDF5, dimsm, memspace_id,
     +                        error)

      do ith=1,3+N_TH
!     Here it starts the loop--->

         select case(ith)
         case (1)
            call SWAPZY(U1,tmp)
            dname="U"
         case (2)
!     Interpolation to the fractional grid
            call G2GF(U2)
            call SWAPZY(U2,tmp)
            call GF2G(U2)
            dname="V"
         case (3)
            call SWAPZY(U3,tmp)
            dname="W"
         case (4:)
            call SWAPZY(TH(0,0,0,ith-3),tmp)
            dname="TH"//CHAR(ith+45)
         end select

         call h5dcreate_f(gid, trim(dname), H5T_IEEE_F64LE,
     +        filspace_id, dset_id, error, dcpl_id = plist_id_d)

!     Select hyperslab in the file.
!     call h5dget_space_f(dsetur_id, selspace_id, error)
         call h5sselect_hyperslab_f (filspace_id, H5S_SELECT_SET_F,
     &        offset, count, error, stride, block)

         offset_m(1:3)=0
         call h5sselect_hyperslab_f (memspace_id, H5S_SELECT_SET_F,
     +        offset_m, count, error, stride, block)

!     Write the dataset collectively
         call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE,
     +        tmp,
     +        dimsm, error, file_space_id = filspace_id,
     +        mem_space_id = memspace_id)

!     Close dateset
         call h5dclose_f(dset_id, error)
      end do

!     In the case of saving for the pressure as well
      if (SAVE_PRESSURE) then
         call fft_xz_to_physical(CP,P,0,NY+1)

         call SWAPZY(P,tmp)
         dname="P"

         call h5dcreate_f(gid, trim(dname), H5T_IEEE_F64LE,
     +        filspace_id, dset_id, error, dcpl_id = plist_id_d)

!     Select hyperslab in the file.
!     call h5dget_space_f(dsetur_id, selspace_id, error)
         call h5sselect_hyperslab_f (filspace_id, H5S_SELECT_SET_F,
     &        offset, count, error, stride, block)

         call h5sselect_hyperslab_f (memspace_id, H5S_SELECT_SET_F,
     +        offset_m, count, error, stride, block)

!     Write the dataset collectively
         call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE,
     +        tmp,
     +        dimsm, error, file_space_id = filspace_id,
     +        mem_space_id = memspace_id)

!     Close dateset
         call h5dclose_f(dset_id, error)
         call fft_xz_to_fourier(P,CP,0,NY+1)
      end if


!     Close the dataspace for the memory and for the file
      call h5sclose_f(filspace_id, error)
      call h5sclose_f(memspace_id, error)

!     Close the properties for the dataspace creation and the writing
      call h5pclose_f(plist_id_d, error)

      ! Close groups
      call h5gclose_f(gid, error)
      call h5fclose_f(file_id, error)
      call h5close_f(error)

      call fft_xz_to_fourier(U1,CU1,0,NY+1)
      call fft_xz_to_fourier(U2,CU2,0,NY+1)
      call fft_xz_to_fourier(U3,CU3,0,NY+1)
      do ith=1,N_TH
         call fft_xz_to_fourier(TH(0,0,0,ith),CTH(0,0,0,ith),0,NY+1)
      end do

      ! call mpi_finalize(ierror)
      ! stop

      end subroutine WriteHDF5



      SUBROUTINE ReadHDF5(FNAME)
      use hdf5

      INCLUDE 'header'

      CHARACTER*55 FNAME
      LOGICAL FINAL, READ_PRESSURE

      REAL*8 tmp(NX,NY,NZ)
c
c     HDF5 ------------------------------------------------------
c
c     Dataset names
      character(len=10) :: dname

c     Identifiers
      integer(hid_t) :: file_id, dset_id
      integer(hid_t) :: filspace_id, memspace_id

!     Identifiers
      integer(hid_t) :: gid, selspace_id

c     Dimensions in the memory and in the file
      integer(hsize_t), dimension(3) :: dimsm,dimsf

      integer(hsize_t), dimension(3) :: count, offset
      integer(hsize_t), dimension(3) :: stride, block, offset_m

      integer :: rHDF5 = 3

      integer(hsize_t),dimension(1)       :: adims
      integer(hid_t)                      :: aid,tspace
      real   , dimension(1)               :: treal(3)
      integer, dimension(1)               :: tint(3)
      character*80                        :: namnbuf
      character*20                        :: sttimec

      integer error, ith

      double precision En(4)

      dimsm(1:3) = (/NX,NY,NZ/)
      dimsf(1:3) = (/NX,NY,NZ/)

!     Flow fields are saved in the fractional grid. We use a basic
!     interpolation
!
!     u_j+1/2=u_j+u_j+1
!
!     on the way back we just invert the relation in order to have
!     exact values.

      block(1) = NX
      block(3) = NZ

!     Stride and count for number of rows and columns in each dimension
      stride = 1
      count  = 1

!     Offset determined by the rank of a processor
      offset(1) = 0
      offset(3) = 0

      block(2) =  NY
      offset(2) = 0

!     Initialize interface
      call h5open_f(error)

!     Create the file collectively
      call h5fopen_f(trim(FNAME), H5F_ACC_RDONLY_F,
     &                 file_id, error)

      adims=3

!     -----------------------------
!     Resolution
!     -----------------------------
      call h5aopen_by_name_f(file_id,'.','Resolution',aid,error)
      call h5aread_f(aid,H5T_NATIVE_INTEGER,tint,adims,error)
      call h5aclose_f(aid, error)
      ! Check that the resolution is of the same kind
      if ((tint(1).ne.NX) .or.
     &    (tint(2).ne.NY) .or.
     &    (tint(3).ne.NZ)    ) then
         if (RANK.eq.0) then
            write(*,*) ' Error. File and program have ',
     &        'different resolutions. '
            write(*,*) ' Program: ', NX,NY,NZ
            write(*,*) ' File   : ', tint(1:3)
         end if
         stop
      end if


      call h5gopen_f(file_id,"/Timestep",gid, error)
!     -----------------------------
!     Time stamps
!     -----------------------------
      call h5aopen_by_name_f(gid,'.','Time',aid,error)
      call h5aread_f(aid,H5T_NATIVE_DOUBLE,TIME,adims,error)
      call h5aclose_f(aid, error)
!     -----------------------------

!     Dataspace in memory
      call h5screate_simple_f(rHDF5, dimsm, memspace_id,
     +                        error)

      do ith=1,3+N_TH
!     Here it starts the loop--->
         select case(ith)
         case (1)
            dname="U"
         case (2)
            dname="V"
         case (3)
            dname="W"
         case (4:)
            dname="TH"//CHAR(ith+45)
         end select

         call h5dopen_f(gid,trim(dname),dset_id,error)
         call h5dget_space_f(dset_id,filspace_id,error)

!     Select hyperslab in the file.
         call h5sselect_hyperslab_f (filspace_id, H5S_SELECT_SET_F,
     &        offset, count, error, stride, block)

         offset_m(1:3)=0
         call h5sselect_hyperslab_f (memspace_id, H5S_SELECT_SET_F,
     +        offset_m, count, error, stride, block)

!     Write the dataset collectively
         call h5dread_f(dset_id, H5T_NATIVE_DOUBLE,
     +        tmp,
     +        dimsm, error, file_space_id = filspace_id,
     +        mem_space_id = memspace_id)

         select case(ith)
         case (1)
            call SWAPYZ(tmp,U1)
         case (2)
            call SWAPYZ(tmp,U2)
!     Interpolation to the collocated grid
            call GF2G(U2)
         case (3)
            call SWAPYZ(tmp,U3)
         case (4:)
            call SWAPYZ(tmp,TH(0,0,0,ith-3))
         end select

!     Close dateset
         call h5sclose_f(filspace_id, error)
         call h5dclose_f(dset_id, error)
      end do

!     Decide whether to compute the pressure or to read
      call h5lexists_f(gid, 'P', READ_PRESSURE, error)
      if (READ_PRESSURE) then
         dname="P"
         call h5dopen_f(gid,trim(dname),dset_id,error)
         call h5dget_space_f(dset_id,filspace_id,error)

!     Select hyperslab in the file.
         call h5sselect_hyperslab_f (filspace_id, H5S_SELECT_SET_F,
     &        offset, count, error, stride, block)

         offset_m(1:3)=0
         call h5sselect_hyperslab_f (memspace_id, H5S_SELECT_SET_F,
     +        offset_m, count, error, stride, block)

!     Write the dataset collectively
         call h5dread_f(dset_id, H5T_NATIVE_DOUBLE,
     +        tmp,
     +        dimsm, error, file_space_id = filspace_id,
     +        mem_space_id = memspace_id) !, xfer_prp = plist_id_w)

         call SWAPYZ(tmp,P)
!     Close dateset
         call h5sclose_f(filspace_id, error)
         call h5dclose_f(dset_id, error)

         call fft_xz_to_fourier(P,CP,0,NY+1)
      end if

!     Close the dataspace for the memory
      call h5sclose_f(memspace_id, error)

      ! Close groups
      call h5gclose_f(gid, error)
      call h5fclose_f(file_id, error)
      call h5close_f(error)

      IF (VARIABLE_DT) THEN
         CALL COURANT
      END IF

      ! Convert to physical space
      call fft_xz_to_fourier(U1,CU1,0,NY+1)
      call fft_xz_to_fourier(U2,CU2,0,NY+1)
      call fft_xz_to_fourier(U3,CU3,0,NY+1)
      do ith=1,N_TH
         call fft_xz_to_fourier(TH(0,0,0,ith),CTH(0,0,0,ith),0,NY+1)
      end do

      if (.NOT.READ_PRESSURE) then
         IF (NUM_PER_DIR.EQ.3) THEN
!     CALL POISSON_P_PER
            CALL REM_DIV_PER
         ELSEIF (NUM_PER_DIR.EQ.2) THEN
            CALL POISSON_P_CHAN
         ELSEIF (NUM_PER_DIR.EQ.1) THEN
            CALL POISSON_P_DUCT
         ELSEIF (NUM_PER_DIR.EQ.0) THEN
            CALL POISSON_P_CAV
         END IF
      end if

      end subroutine ReadHDF5


      subroutine ReadGridHDF5(FNAME,coord)
      use hdf5

      implicit none

      INCLUDE 'header'

      CHARACTER*55 FNAME

c     Identifiers
      integer(hid_t) :: file_id, dset_id
      integer(hid_t) :: filspace_id, memspace_id
!     Identifiers
      integer(hid_t) :: gid, selspace_id

c     Dataset names
      integer           :: coord
      character(len=10) :: cname
c     Dimensions in the memory and in the file
      integer(hsize_t),dimension(1)  :: dimsm,dimsf
      integer(hsize_t),dimension(1)  :: idims,imaxd

      integer(hsize_t),dimension(1)  :: count, offset
      integer(hsize_t),dimension(1)  :: stride, block, offset_m

      integer           :: error,rHDF5,ith

!     Initialize interface
      call h5open_f(error)

!     Create the file collectively
      call h5fopen_f(trim(FNAME), H5F_ACC_RDONLY_F,
     &                 file_id, error)

c$$$!     IF ONE WANTS TO PUT A CHECK FOR THE RESOLUTION
c$$$!     -----------------------------
c$$$!     Resolution
c$$$!     -----------------------------
c$$$      adims=3
c$$$      call h5aopen_by_name_f(file_id,'.','Resolution',aid,error)
c$$$      call h5aread_f(aid,H5T_NATIVE_INTEGER,tint,adims,error)
c$$$      call h5aclose_f(aid, error)

      select case(coord)
      case(1)
         write(*,*) ' Error 235454. Not implemented yet! '
         stop
      case(2)
         cname='y'

         dimsm=NY+2
         dimsf=NY

!     Stride and count for number of rows and columns in each dimension
         stride = 1
         count  = 1

!     Offset determined by the rank of a processor
         block  =  NY+1
         offset =  0
      case(3)
         write(*,*) ' Error 235455. Not implemented yet! '
         stop
      end select

      call h5gopen_f(file_id,"/grids",gid, error)

!     Dataspace in memory
      rHDF5=1
      call h5screate_simple_f(rHDF5, dimsm, memspace_id,
     +                        error)

      call h5dopen_f(gid,trim(cname),dset_id,error)
      call h5dget_space_f(dset_id,filspace_id,error)

      ! Check for the dimensions
      call h5sget_simple_extent_dims_f(filspace_id,idims,imaxd,error)
      if (idims(1)-1.ne.dimsf(1)) then
         write(*,*) ' Grid file and program do not match. '
         write(*,*) '   gridfile (',trim(cname),'): ',dimsf
         write(*,*) '   program     : '              ,idims-1
         stop
      end if

!     Select hyperslab in the file.
      call h5sselect_hyperslab_f (filspace_id, H5S_SELECT_SET_F,
     &     offset, count, error, stride, block)

      offset_m=1
      call h5sselect_hyperslab_f (memspace_id, H5S_SELECT_SET_F,
     +     offset_m, count, error, stride, block)


!     Write the dataset collectively
      call h5dread_f(dset_id, H5T_NATIVE_DOUBLE,
     +     GY,
     +     dimsm, error, file_space_id = filspace_id,
     +     mem_space_id = memspace_id)

!     Close dateset
      call h5sclose_f(filspace_id, error)
      call h5dclose_f(dset_id, error)

!     Close the dataspace for the memory
      call h5sclose_f(memspace_id, error)

      ! Close groups
      call h5gclose_f(gid, error)
      call h5fclose_f(file_id, error)
      call h5close_f(error)

      ! Calculate the GYF in the interior
      do ith=1,NY
         GYF(ith)=0.5*(GY(ith)+GY(ith+1))
      end do

      ! ###############################
      !    Get the outer ghost cells
      ! ###############################

      GYF(0) = 2.d0*GYF(1)-GYF(2)
      GYF(NY+1)=2.d0*GYF(NY)-GYF(NY-1)

      end subroutine













      subroutine G2GF(var)

      implicit none

      include 'header'

      REAL*8 var(0:NX+1,0:NZ+1,0:NY+1)
      INTEGER X,Z,Y

      do x=0,NX-1
      do z=0,NZ-1
      var(x,z,1)=var(x,z,2)
      var(x,z,NY+1)=var(x,z,NY)
      do y=1,NY
         var(x,z,y)=0.5*(var(x,z,y)+var(x,z,y+1))
      end do
      end do
      end do

      end subroutine



      subroutine GF2G(var)

      implicit none

      include 'header'

      REAL*8 var(0:NX+1,0:NZ+1,0:NY+1)
      INTEGER X,Z,Y
      INTEGER XZBOX

      do x=0,NX-1
         do z=0,NZ-1
            var(x,z,NY+1)=var(x,z,NY)
         end do
      end do

      do x=0,NX-1
         do z=0,NZ-1
            do y=NY,1,-1
               var(x,z,y)=2*var(x,z,y)-var(x,z,y+1)
            end do
         end do
      end do

      ! Impose the values at the boundary as prescribed in the
      ! code in order to have zero mass flux
      var(:,:,1)=-var(:,:,2)
      var(:,:,NY+1)=-var(:,:,NY)

      end subroutine







      subroutine SWAPZY(in,out)

      include 'header'

      REAL*8 in(0:NX+1,0:NZ+1,0:NY+1)
      REAL*8 out(1:NX,1:NY,1:NZ)
      INTEGER X,Z,Y

      out=0.d0
      do x=0,NX-1
         do y=1,NY
            do z=0,NZ-1
               out(x+1,y,z+1)=in(x,z,y)
            end do
         end do
      end do

      end subroutine




      subroutine SWAPYZ(in,out)

      include 'header'

      REAL*8 out(0:NX+1,0:NZ+1,0:NY+1)
      REAL*8 in(1:NX,1:NY,1:NZ)
      INTEGER X,Z,Y

      out=0.d0
      do x=0,NX-1
         do y=1,NY
            do z=0,NZ-1
               out(x,z,y)=in(x+1,y,z+1)
            end do
         end do
      end do

      end subroutine



      subroutine time_string(cdt)
c
c     Construct string in the format '19-DEC-2005 22:47:06'
c
      implicit none

      integer i

      integer val(8)
      character*20 cdt
      character*3 monc

      call date_and_time(values=val)

      if (val(2).eq.1) then
         monc  = 'JAN'
      else if (val(2).eq.2) then
         monc  = 'FEB'
      else if (val(2).eq.3) then
         monc  = 'MAR'
      else if (val(2).eq.4) then
         monc  = 'APR'
      else if (val(2).eq.5) then
         monc  = 'MAY'
      else if (val(2).eq.6) then
         monc  = 'JUN'
      else if (val(2).eq.7) then
         monc  = 'JUL'
      else if (val(2).eq.8) then
         monc  = 'AUG'
      else if (val(2).eq.9) then
         monc  = 'SEP'
      else if (val(2).eq.10) then
         monc  = 'OCT'
      else if (val(2).eq.11) then
         monc  = 'NOV'
      else if (val(2).eq.12) then
         monc  = 'DEC'
      else
         monc  = 'XXX'
      end if

      write(cdt,'(i2,a1,a3,a1,i4,a1,i2,a1,i2,a1,i2)')
     &     val(3),'-',monc,'-',val(1),' ',val(5),':',val(6),':',val(7)
      do i=1,2
         if (cdt(i:i).eq.' ') then
            cdt(i:i)='0'
         end if
      end do
      do i=13,20
         if (cdt(i:i).eq.' ') then
            cdt(i:i)='0'
         end if
      end do

      end subroutine time_string




c$$$
c$$$      SUBROUTINE WriteHDF5_var_real(FNAME)
c$$$      use hdf5
c$$$
c$$$      INCLUDE 'header'
c$$$
c$$$      CHARACTER*55 FNAME
c$$$      LOGICAL FINAL
c$$$
c$$$      REAL*8 tmp(NX,NY,NZ)
c$$$      ! REAL*8 var(0:NX+1,0:NZ+1,0:NY+1)
c$$$c
c$$$c     HDF5 ------------------------------------------------------
c$$$c
c$$$c     Dataset names
c$$$      character(len=10) :: dname
c$$$
c$$$c     Identifiers
c$$$      integer(hid_t) :: file_id, dset_id
c$$$      integer(hid_t) :: filspace_id, memspace_id
c$$$
c$$$!     Identifiers
c$$$      integer(hid_t) :: gid, selspace_id
c$$$      integer(hid_t) :: plist_id_d
c$$$
c$$$c     Dimensions in the memory and in the file
c$$$      integer(hsize_t), dimension(3) :: dimsm,dimsf
c$$$
c$$$      integer(hsize_t), dimension(3) :: chunk_dims, count, offset
c$$$      integer(hsize_t), dimension(3) :: stride, block, offset_m
c$$$
c$$$      integer :: rHDF5 = 3, arank = 1
c$$$
c$$$      integer(hsize_t),dimension(1)       :: adims
c$$$      integer(hid_t)                      :: aid,tspace
c$$$      real   , dimension(1)               :: treal(3)
c$$$      integer, dimension(1)               :: tint(3)
c$$$      character*80                        :: namnbuf
c$$$      character*20                        :: sttimec
c$$$
c$$$      integer error, ith
c$$$
c$$$      double precision En(4)
c$$$
c$$$      dimsm(1:3) = (/NX,NY,NZ/)
c$$$      dimsf(1:3) = (/NX,NY,NZ/)
c$$$
c$$$      write(100,'(1E)') tvar(5,5,:)
c$$$
c$$$!     Flow fields are saved in the fractional grid. We use a basic
c$$$!     interpolation
c$$$!
c$$$!     u_j+1/2=u_j+u_j+1
c$$$!
c$$$!     on the way back we just invert the relation in order to have
c$$$!     exact values.
c$$$!     The interpolation formula is anyway not too important. Anything
c$$$!     could be used. The important is rather to get back exactly the
c$$$!     same when reading the file. At the boundary, we impose Neumann
c$$$!     condition which allows to pertain in the system N
c$$$!
c$$$!     NOTE that inverting the formula above requires a solution of a
c$$$!     tridiagonal system in the vertical direction. This is similar
c$$$!     to the Thomas algoritm in the implemantation, in the sensa that
c$$$!     a pipeline strategy must be used. Thus, the parallel performance
c$$$!     is rather poor
c$$$
c$$$      WRITE(6,*) 'Writing flow to ',FNAME
c$$$
c$$$      chunk_dims(1) = NX
c$$$      chunk_dims(2) = 1
c$$$      chunk_dims(3) = NZ
c$$$
c$$$      block(1) = NX
c$$$      block(3) = NZ
c$$$
c$$$!     Stride and count for number of rows and columns in each dimension
c$$$      stride = 1
c$$$      count  = 1
c$$$
c$$$!     Offset determined by the rank of a processor
c$$$      offset(1) = 0
c$$$      offset(3) = 0
c$$$
c$$$      block(2) =  NY
c$$$      offset(2) = 0
c$$$
c$$$!     Initialize interface
c$$$      call h5open_f(error)
c$$$
c$$$!     Create the file collectively
c$$$      call h5fcreate_f(trim(FNAME), H5F_ACC_TRUNC_F,
c$$$     &                 file_id, error)
c$$$
c$$$!      ! Convert to physical space
c$$$!      call fft_xz_to_physical(CU1,U1,0,NY+1)
c$$$
c$$$!     Create property list for the chunked dataset creation
c$$$      call h5pcreate_f(H5P_DATASET_CREATE_F, plist_id_d, error)
c$$$      call h5pset_chunk_f(plist_id_d, rHDF5, chunk_dims, error)
c$$$
c$$$!     Create the dataspace for ur
c$$$      call h5screate_simple_f(rHDF5, dimsf, filspace_id, error)
c$$$      call h5screate_simple_f(rHDF5, dimsm, memspace_id,
c$$$     +                        error)
c$$$
c$$$      call SWAPZY(tvar,tmp)
c$$$      dname="U"
c$$$
c$$$      call h5gcreate_f(file_id,"/Timestep",gid, error)
c$$$
c$$$      call h5dcreate_f(gid, trim(dname), H5T_IEEE_F64LE,
c$$$     +     filspace_id, dset_id, error, dcpl_id = plist_id_d)
c$$$
c$$$!     Select hyperslab in the file.
c$$$!     call h5dget_space_f(dsetur_id, selspace_id, error)
c$$$      call h5sselect_hyperslab_f (filspace_id, H5S_SELECT_SET_F,
c$$$     &     offset, count, error, stride, block)
c$$$
c$$$      offset_m(1:3)=0
c$$$      call h5sselect_hyperslab_f (memspace_id, H5S_SELECT_SET_F,
c$$$     +     offset_m, count, error, stride, block)
c$$$
c$$$!     Write the dataset collectively
c$$$      call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE,
c$$$     +     tmp,
c$$$     +     dimsm, error, file_space_id = filspace_id,
c$$$     +     mem_space_id = memspace_id)
c$$$
c$$$!     Close dateset
c$$$      call h5dclose_f(dset_id, error)
c$$$
c$$$!     Close the dataspace for the memory and for the file
c$$$      call h5sclose_f(filspace_id, error)
c$$$      call h5sclose_f(memspace_id, error)
c$$$
c$$$!     Close the properties for the dataspace creation and the writing
c$$$      call h5pclose_f(plist_id_d, error)
c$$$
c$$$      ! Close groups
c$$$      call h5gclose_f(gid, error)
c$$$      call h5fclose_f(file_id, error)
c$$$      call h5close_f(error)
c$$$
c$$$      stop
c$$$
c$$$      end subroutine WriteHDF5_var_real
c$$$
c$$$
c$$$
c$$$
c$$$      SUBROUTINE WriteHDF5_var_complex(FNAME)
c$$$      use hdf5
c$$$
c$$$      INCLUDE 'header'
c$$$
c$$$      CHARACTER*55 FNAME
c$$$      LOGICAL FINAL
c$$$
c$$$      REAL*8 tmp(NX,NY,NZ)
c$$$      !REAL*8 var(0:NX+1,0:NZ+1,0:NY+1)
c$$$      !COMPLEX*16 cvar(0:NX/2,0:NZ+1,0:NY+1)
c$$$      !EQUIVALENCE (var,cvar)
c$$$c
c$$$c     HDF5 ------------------------------------------------------
c$$$c
c$$$c     Dataset names
c$$$      character(len=10) :: dname
c$$$
c$$$c     Identifiers
c$$$      integer(hid_t) :: file_id, dset_id
c$$$      integer(hid_t) :: filspace_id, memspace_id
c$$$
c$$$!     Identifiers
c$$$      integer(hid_t) :: gid, selspace_id
c$$$      integer(hid_t) :: plist_id_d
c$$$
c$$$c     Dimensions in the memory and in the file
c$$$      integer(hsize_t), dimension(3) :: dimsm,dimsf
c$$$
c$$$      integer(hsize_t), dimension(3) :: chunk_dims, count, offset
c$$$      integer(hsize_t), dimension(3) :: stride, block, offset_m
c$$$
c$$$      integer :: rHDF5 = 3, arank = 1
c$$$
c$$$      integer(hsize_t),dimension(1)       :: adims
c$$$      integer(hid_t)                      :: aid,tspace
c$$$      real   , dimension(1)               :: treal(3)
c$$$      integer, dimension(1)               :: tint(3)
c$$$      character*80                        :: namnbuf
c$$$      character*20                        :: sttimec
c$$$
c$$$      integer error, ith
c$$$
c$$$      double precision En(4)
c$$$
c$$$      dimsm(1:3) = (/NX,NY,NZ/)
c$$$      dimsf(1:3) = (/NX,NY,NZ/)
c$$$
c$$$      write(100,'(2E)') ctvar(5,5,:)
c$$$
c$$$!     Flow fields are saved in the fractional grid. We use a basic
c$$$!     interpolation
c$$$!
c$$$!     u_j+1/2=u_j+u_j+1
c$$$!
c$$$!     on the way back we just invert the relation in order to have
c$$$!     exact values.
c$$$!     The interpolation formula is anyway not too important. Anything
c$$$!     could be used. The important is rather to get back exactly the
c$$$!     same when reading the file. At the boundary, we impose Neumann
c$$$!     condition which allows to pertain in the system N
c$$$!
c$$$!     NOTE that inverting the formula above requires a solution of a
c$$$!     tridiagonal system in the vertical direction. This is similar
c$$$!     to the Thomas algoritm in the implemantation, in the sensa that
c$$$!     a pipeline strategy must be used. Thus, the parallel performance
c$$$!     is rather poor
c$$$
c$$$      WRITE(6,*) 'Writing flow to ',FNAME
c$$$
c$$$      chunk_dims(1) = NX
c$$$      chunk_dims(2) = 1
c$$$      chunk_dims(3) = NZ
c$$$
c$$$      block(1) = NX
c$$$      block(3) = NZ
c$$$
c$$$!     Stride and count for number of rows and columns in each dimension
c$$$      stride = 1
c$$$      count  = 1
c$$$
c$$$!     Offset determined by the rank of a processor
c$$$      offset(1) = 0
c$$$      offset(3) = 0
c$$$
c$$$      block(2) =  NY
c$$$      offset(2) = 0
c$$$
c$$$!     Initialize interface
c$$$      call h5open_f(error)
c$$$
c$$$!     Create the file collectively
c$$$      call h5fcreate_f(trim(FNAME), H5F_ACC_TRUNC_F,
c$$$     &                 file_id, error)
c$$$
c$$$!      ! Convert to physical space
c$$$      call fft_xz_to_physical(ctvar,tvar,0,NY+1)
c$$$
c$$$!     Create property list for the chunked dataset creation
c$$$      call h5pcreate_f(H5P_DATASET_CREATE_F, plist_id_d, error)
c$$$      call h5pset_chunk_f(plist_id_d, rHDF5, chunk_dims, error)
c$$$
c$$$!     Create the dataspace for ur
c$$$      call h5screate_simple_f(rHDF5, dimsf, filspace_id, error)
c$$$      call h5screate_simple_f(rHDF5, dimsm, memspace_id,
c$$$     +                        error)
c$$$
c$$$      call SWAPZY(tvar,tmp)
c$$$      dname="U"
c$$$
c$$$      call h5gcreate_f(file_id,"/Timestep",gid, error)
c$$$
c$$$      call h5dcreate_f(gid, trim(dname), H5T_IEEE_F64LE,
c$$$     +     filspace_id, dset_id, error, dcpl_id = plist_id_d)
c$$$
c$$$!     Select hyperslab in the file.
c$$$!     call h5dget_space_f(dsetur_id, selspace_id, error)
c$$$      call h5sselect_hyperslab_f (filspace_id, H5S_SELECT_SET_F,
c$$$     &     offset, count, error, stride, block)
c$$$
c$$$      offset_m(1:3)=0
c$$$      call h5sselect_hyperslab_f (memspace_id, H5S_SELECT_SET_F,
c$$$     +     offset_m, count, error, stride, block)
c$$$
c$$$!     Write the dataset collectively
c$$$      call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE,
c$$$     +     tmp,
c$$$     +     dimsm, error, file_space_id = filspace_id,
c$$$     +     mem_space_id = memspace_id)
c$$$
c$$$!     Close dateset
c$$$      call h5dclose_f(dset_id, error)
c$$$
c$$$!     Close the dataspace for the memory and for the file
c$$$      call h5sclose_f(filspace_id, error)
c$$$      call h5sclose_f(memspace_id, error)
c$$$
c$$$!     Close the properties for the dataspace creation and the writing
c$$$      call h5pclose_f(plist_id_d, error)
c$$$
c$$$      ! Close groups
c$$$      call h5gclose_f(gid, error)
c$$$      call h5fclose_f(file_id, error)
c$$$      call h5close_f(error)
c$$$
c$$$      stop
c$$$
c$$$      end subroutine WriteHDF5_var_complex
c$$$
c$$$
c$$$
c$$$
c$$$
c$$$
