module hdf5_f2f

  use hdf5
  implicit none

  integer(HID_T)   :: fileid
  integer(HID_T)   :: xferid
  integer(HID_T)   :: propid
  integer(HID_T)   :: dsetid
  integer(HID_T)   :: fspcid
  integer(HID_T)   :: mspcid
  integer(HID_T)   :: dspcid
  integer(HSIZE_T) :: dimsf(7)
  integer(HSIZE_T) :: ndimsf

integer, parameter :: ipar_out = 0

  interface fh5_write
     module procedure                 &
          fh5_write_integer_scalar,   &
          fh5_write_real_scalar,      &
          fh5_write_character_scalar, &
          fh5_write_real8_scalar,     &
          fh5_write_logical_scalar,   &
          fh5_write_integer_array,    &
          fh5_write_real_array,       &
          fh5_write_character_array,  &
          fh5_write_real8_array,      &
          fh5_write_logical_array
  end interface fh5_write

  interface fh5d_read
     module procedure                  &
          fh5d_read_integer_scalar,    &
          fh5d_read_real_scalar,       &
          fh5d_read_character_scalar,  &
          fh5d_read_real8_scalar,      &
          fh5d_read_logical_scalar,    &
          fh5d_read_integer_array,     &
          fh5d_read_real_array,        &
          fh5d_read_character_array,   &
          fh5d_read_real8_array,       &
          fh5d_read_logical_array
  end interface fh5d_read


contains


  subroutine fh5f_open(locfn, iaccess, hdferr)
    implicit none

    character(*), intent(IN)  :: locfn
    integer,      intent(IN)  :: iaccess
    integer,      intent(OUT) :: hdferr

    integer(HID_T) :: access_id
    integer        :: flags
    integer        :: nulo

    ! turn off default error handling
    call h5eset_auto_f(0, nulo)

    access_id = H5P_DEFAULT_F
    if(iaccess == 1) flags = H5F_ACC_RDONLY_F
    if(iaccess == 2) flags = H5F_ACC_RDWR_F

    call h5fopen_f(locfn, flags, fileid, hdferr, access_id)

    hdferr = int(fileid) + hdferr

    return
  end subroutine fh5f_open



  subroutine fh5f_close(hdferr)
    implicit none

    integer, intent(OUT) :: hdferr

    call h5fclose_f(fileid, hdferr)

    return
  end subroutine fh5f_close



  subroutine fh5f_create(locfn, iaccess, hdferr)

#if defined(OLAM_MPI) && defined(OLAM_PARALLEL_HDF5)
    use misc_coms, only: ipar_out, io6
    use mpi
#endif

    implicit none

    character(*), intent(IN)  :: locfn
    integer,      intent(IN)  :: iaccess
    integer,      intent(OUT) :: hdferr

    integer(HID_T) :: access_id
    integer(HID_T) :: create_id
    integer        :: flags

    create_id = H5P_DEFAULT_F
    access_id = H5P_DEFAULT_F
    if (iaccess == 1) flags = H5F_ACC_TRUNC_F
    if (iaccess == 2) flags = H5F_ACC_EXCL_F

#if defined(OLAM_MPI) && defined(OLAM_PARALLEL_HDF5)
    if (ipar_out == 1) then
       write(io6,*)
       write(io6,*) "Enabling parallel HDF5 output"
       write(io6,*)
       call h5pcreate_f(H5P_FILE_ACCESS_F, access_id, hdferr)
       call h5pset_fapl_mpio_f(access_id, MPI_COMM_WORLD, MPI_INFO_NULL, hdferr)
       ! call h5pset_fapl_mpiposix_f(access_id, MPI_COMM_WORLD, .true., hdferr)
    endif
#endif

    call h5fcreate_f(locfn, flags, fileid, hdferr, create_id, access_id)

#if defined(OLAM_MPI) && defined(OLAM_PARALLEL_HDF5)
    if (ipar_out == 1) then
       call h5pclose_f(access_id, hdferr)
    endif
#endif

    hdferr = int(fileid) + hdferr

    return
  end subroutine fh5f_create



  subroutine fh5_prepare_write(ndims, dims, hdferr, icompress, &
                               mcoords, fcoords, ifsize)
!    use misc_coms, only: ipar_out
    implicit none

    integer, intent(IN)           :: ndims
    integer, intent(IN)           :: dims(:)
    integer, intent(OUT)          :: hdferr
    integer, intent(IN), optional :: icompress
    integer, intent(IN), optional :: mcoords(:), fcoords(:), ifsize

    integer          :: i, j, iop
    integer(HSIZE_T) :: dims_file(ndims)
    integer(HSIZE_T) :: offset1d(1), count1d(1)
    integer(HSIZE_T) :: offset2d(2), count2d(2)
    integer(HSIZE_T) :: offset3d(3), count3d(3)

    ! Output dimensions

    dimsf = 1
    ndimsf = ndims

    do i = 1, ndims
       dimsf(i) = dims(i)
    end do

    dims_file = 1
     do i = 1, ndims
       dims_file(i) = dims(i)
    end do
   
    ! For distributed output, global data size will be different!

    if (present(fcoords) .and. present(ifsize)) then
       dims_file(ndims) = ifsize
    endif

    ! Create a property list for compression/chunking/filters

    call h5pcreate_f(H5P_DATASET_CREATE_F, propid, hdferr) 	 

    ! If no parallel IO, activate HDF5 compression for valid icompress

    if (present(icompress)) then
       if (icompress > 0 .and. icompress < 10 .and. ipar_out /= 1) then
          if (.not. (ndims == 1 .and. dimsf(1) == 1)) then
             call h5pset_chunk_f(propid, ndims, dimsf, hdferr)
             call h5pset_shuffle_f(propid, hdferr)
             call h5pset_deflate_f(propid, icompress, hdferr)	 
          endif
       endif
    endif

    ! Create transfer property list for collective dataset write

    xferid = H5P_DEFAULT_F

#if defined(OLAM_MPI) && defined(OLAM_PARALLEL_HDF5)
    if (ipar_out == 1) then
       call h5pcreate_f(h5p_dataset_xfer_f, xferid, hdferr)
       call h5pset_dxpl_mpio_f(xferid, h5fd_mpio_collective_f, hdferr)
    endif
#endif

    ! Create the global file space for the data

    call h5screate_simple_f(ndims, dims_file, fspcid, hdferr)
    
    ! Create the local memory space for the data

    call h5screate_simple_f(ndims, dimsf, mspcid, hdferr)

    ! If we are only writing a subset of the data, pick the points
    ! that we will output

    if (present(mcoords)) then

       if (size(mcoords) == 0) then

          call h5sselect_none_f(mspcid, hdferr)

       elseif (ndims == 1) then

          offset1d = mcoords(1) - 1
          count1d  = 1
          iop      = H5S_SELECT_SET_F

          do j = 2, size(mcoords)
             if (mcoords(j) /= mcoords(j-1) + 1) then
                call h5sselect_hyperslab_f(mspcid, iop, offset1d, &
                     count1d, hdferr)
                offset1d = mcoords(j) - 1
                count1d  = 1
                iop      = H5S_SELECT_OR_F
             else
                count1d  = count1d + 1
             endif
          enddo

          call h5sselect_hyperslab_f(mspcid, iop, offset1d, count1d, hdferr)

       elseif (ndims == 2) then

          offset2d = (/ 0, mcoords(1)-1 /)
          count2d  = (/ dims(1), 1 /)
          iop      = H5S_SELECT_SET_F
          
          do j = 2, size(mcoords)
             if (mcoords(j) /= mcoords(j-1) + 1) then
                call h5sselect_hyperslab_f(mspcid, iop, offset2d, &
                     count2d, hdferr)
                offset2d(2) = mcoords(j) - 1
                count2d(2)  = 1
                iop         = H5S_SELECT_OR_F
             else
                count2d(2)  = count2d(2) + 1
             endif
          enddo
          
          call h5sselect_hyperslab_f(mspcid, iop, offset2d, count2d, hdferr)

       elseif (ndims == 3) then

          offset3d = (/ 0, 0, mcoords(1)-1 /)
          count3d  = (/ dims(1), dims(2), 1 /)
          iop      = H5S_SELECT_SET_F

          do j = 2, size(mcoords)
             if (mcoords(j) /= mcoords(j-1) + 1) then
                call h5sselect_hyperslab_f(mspcid, iop, offset3d, &
                     count3d, hdferr)
                offset3d(3) = mcoords(j) - 1
                count3d(3)  = 1
                iop         = H5S_SELECT_OR_F
             else
                count3d(3)  = count3d(3) + 1
             endif
          enddo
          
          call h5sselect_hyperslab_f(mspcid, iop, offset3d, count3d, hdferr)

       else

          stop 'ndims > 3 using parallel I/O is not implemented'
          
       endif
             
    endif

    ! Create the local data space for the data

    call h5screate_simple_f(ndims, dims_file, dspcid, hdferr)

    ! For parallel output, pick the points in the output file that
    ! we will be writing to

    if (present(fcoords) .and. present(ifsize)) then

       if (size(fcoords) == 0) then

          call h5sselect_none_f(dspcid, hdferr)

       elseif (ndims == 1) then

          offset1d = fcoords(1) - 1
          count1d  = 1
          iop      = H5S_SELECT_SET_F

          do j = 2, size(fcoords)
             if (fcoords(j) /= fcoords(j-1) + 1) then
                call h5sselect_hyperslab_f(dspcid, iop, offset1d, &
                     count1d, hdferr)
                offset1d = fcoords(j) - 1
                count1d  = 1
                iop      = H5S_SELECT_OR_F
             else
                count1d  = count1d + 1
             endif
          enddo

          call h5sselect_hyperslab_f(dspcid, iop, offset1d, count1d, hdferr)

       elseif (ndims == 2) then

          offset2d = (/ 0, fcoords(1)-1/)
          count2d  = (/ dims(1), 1 /)
          iop      = H5S_SELECT_SET_F
          
          do j = 2, size(fcoords)
             if (fcoords(j) /= fcoords(j-1) + 1) then
                call h5sselect_hyperslab_f(dspcid, iop, offset2d, &
                     count2d, hdferr)
                offset2d(2) = fcoords(j) - 1
                count2d(2)  = 1
                iop         = H5S_SELECT_OR_F
             else
                count2d(2)  = count2d(2) + 1
             endif
          enddo
          
          call h5sselect_hyperslab_f(dspcid, iop, offset2d, count2d, hdferr)

       elseif (ndims == 3) then

          offset3d = (/ 0, 0, fcoords(1)-1 /)
          count3d  = (/ dims(1), dims(2), 1 /)
          iop      = H5S_SELECT_SET_F

          do j = 2, size(fcoords)
             if (fcoords(j) /= fcoords(j-1) + 1) then
                call h5sselect_hyperslab_f(dspcid, iop, offset3d, &
                     count3d, hdferr)
                offset3d(3) = fcoords(j) - 1
                count3d(3)  = 1
                iop         = H5S_SELECT_OR_F
             else
                count3d(3)  = count3d(3) + 1
             endif
          enddo
          
          call h5sselect_hyperslab_f(dspcid, iop, offset3d, count3d, hdferr)

       else

          stop 'ndims > 3 using parallel I/O is not implemented'
          
       endif
             
    endif

    return
  end subroutine fh5_prepare_write



  subroutine fh5_write_integer_array(buf_integer, dname, hdferr)
    implicit none

    integer,      intent(IN)  :: buf_integer(*)
    character(*), intent(IN)  :: dname
    integer,      intent(OUT) :: hdferr

    call h5dcreate_f(fileid, dname, H5T_NATIVE_INTEGER, fspcid, dsetid, hdferr, propid)

    call h5dwrite_f(dsetid, H5T_NATIVE_INTEGER, buf_integer, dimsf, hdferr, mspcid, dspcid, xferid)

    return
  end subroutine fh5_write_integer_array



  subroutine fh5_write_real_array(buf_real, dname, hdferr)
    implicit none
 
    real,         intent(IN)  :: buf_real(*)
    character(*), intent(IN)  :: dname
    integer,      intent(OUT) :: hdferr

    call h5dcreate_f(fileid, dname, H5T_NATIVE_REAL, fspcid, dsetid, hdferr, propid)
    
    call h5dwrite_f(dsetid, H5T_NATIVE_REAL, buf_real, dimsf, hdferr, mspcid, dspcid, xferid)

    return
  end subroutine fh5_write_real_array



  subroutine fh5_write_character_array(buf_character, dname, hdferr)
    implicit none

    character(*), intent(IN)  :: buf_character(*)
    character(*), intent(IN)  :: dname
    integer,      intent(OUT) :: hdferr

    call h5dcreate_f(fileid, dname, H5T_NATIVE_CHARACTER, fspcid, dsetid, hdferr, propid)

    call h5dwrite_f(dsetid, H5T_NATIVE_CHARACTER, buf_character, dimsf, hdferr, mspcid, dspcid, xferid)

    return
  end subroutine fh5_write_character_array



  subroutine fh5_write_real8_array(buf_real8, dname, hdferr)
    implicit none

    real(kind=8), intent(IN)  :: buf_real8(*)
    character(*), intent(IN)  :: dname
    integer,      intent(OUT) :: hdferr

    call h5dcreate_f(fileid, dname, H5T_NATIVE_DOUBLE, fspcid, dsetid, hdferr, propid)

    call h5dwrite_f(dsetid, H5T_NATIVE_DOUBLE, buf_real8, dimsf, hdferr, mspcid, dspcid, xferid)

    return
  end subroutine fh5_write_real8_array



  subroutine fh5_write_logical_array(buf_logical, dname, hdferr)
    implicit none

    logical,      intent(IN)  :: buf_logical(*)
    character(*), intent(IN)  :: dname
    integer,      intent(OUT) :: hdferr
    integer,      allocatable :: buf_integer(:)
    integer                   :: i, size

    call h5dcreate_f(fileid, dname, H5T_NATIVE_INTEGER, fspcid, dsetid, hdferr, propid)

    ! converting logical to integer
    size = 1
    do i = 1, ndimsf
       size = size * dimsf(i)
    enddo

    allocate (buf_integer(size))
    buf_integer = 0

    do i = 1, size
       if(buf_logical(i)) buf_integer(i) = -1
    enddo

    call h5dwrite_f(dsetid, H5T_NATIVE_INTEGER, buf_integer, dimsf, hdferr, mspcid, dspcid, xferid)

    deallocate(buf_integer)

    return
  end subroutine fh5_write_logical_array



  subroutine fh5_write_integer_scalar(buf_integer, dname, hdferr)
    implicit none

    integer,      intent(IN)  :: buf_integer
    character(*), intent(IN)  :: dname
    integer,      intent(OUT) :: hdferr

    call h5dcreate_f(fileid, dname, H5T_NATIVE_INTEGER, fspcid, dsetid, hdferr, propid)

    call h5dwrite_f(dsetid, H5T_NATIVE_INTEGER, buf_integer, dimsf, hdferr, mspcid, dspcid, xferid)

    return
  end subroutine fh5_write_integer_scalar



  subroutine fh5_write_real_scalar(buf_real, dname, hdferr)
    implicit none

    real,         intent(IN)  :: buf_real
    character(*), intent(IN)  :: dname
    integer,      intent(OUT) :: hdferr

    call h5dcreate_f(fileid, dname, H5T_NATIVE_REAL, fspcid, dsetid, hdferr, propid)

    call h5dwrite_f(dsetid, H5T_NATIVE_REAL, buf_real, dimsf, hdferr, mspcid, dspcid, xferid)

    return
  end subroutine fh5_write_real_scalar



  subroutine fh5_write_character_scalar(buf_character, dname, hdferr)
    implicit none

    character(*), intent(IN)  :: buf_character
    character(*), intent(IN)  :: dname
    integer,      intent(OUT) :: hdferr

    call h5dcreate_f(fileid, dname, H5T_NATIVE_CHARACTER, fspcid, dsetid, hdferr, propid)

    call h5dwrite_f(dsetid, H5T_NATIVE_CHARACTER, buf_character, dimsf, hdferr, mspcid, dspcid, xferid)

    return
  end subroutine fh5_write_character_scalar



  subroutine fh5_write_real8_scalar(buf_real8, dname, hdferr)
    implicit none

    real(kind=8), intent(IN)  :: buf_real8
    character(*), intent(IN)  :: dname
    integer,      intent(OUT) :: hdferr

    call h5dcreate_f(fileid, dname, H5T_NATIVE_DOUBLE, fspcid, dsetid, hdferr, propid)

    call h5dwrite_f(dsetid, H5T_NATIVE_DOUBLE, buf_real8, dimsf, hdferr, mspcid, dspcid, xferid)

    return
  end subroutine fh5_write_real8_scalar



  subroutine fh5_write_logical_scalar(buf_logical, dname, hdferr)
    implicit none

    logical,      intent(IN)  :: buf_logical
    character(*), intent(IN)  :: dname
    integer,      intent(OUT) :: hdferr
    integer                   :: buf_integer

    call h5dcreate_f(fileid, dname, H5T_NATIVE_INTEGER, fspcid, dsetid, hdferr, propid)

    buf_integer = 0
    if(buf_logical) buf_integer = -1

    call h5dwrite_f(dsetid, H5T_NATIVE_INTEGER, buf_integer, dimsf, hdferr, mspcid, dspcid, xferid)

    return
  end subroutine fh5_write_logical_scalar



  subroutine fh5_close_write(hdferr)
    implicit none

    integer, intent(OUT) :: hdferr

    call h5sclose_f(fspcid, hdferr)
    call h5sclose_f(dspcid, hdferr)
    call h5sclose_f(mspcid, hdferr)

    call h5dclose_f(dsetid, hdferr)
    call h5pclose_f(xferid, hdferr)
    call h5pclose_f(propid, hdferr)

    return
  end subroutine fh5_close_write



  subroutine fh5d_open(dname, hdferr)
    implicit none

    character(*), intent(IN)  :: dname
    integer,      intent(OUT) :: hdferr

    call h5dopen_f(fileid, dname, dsetid, hdferr)

    if(dsetid < 0) then
       hdferr = int(dsetid)
       return
    end if

    call h5dget_space_f(dsetid, dspcid, hdferr)

    hdferr = int(dspcid) + hdferr

    return
  end subroutine fh5d_open



  subroutine fh5s_get_ndims(ndims)
    implicit none

    integer, intent(OUT) :: ndims

    integer :: hdferr

    call h5sget_simple_extent_ndims_f(dspcid, ndims, hdferr)

    return
  end subroutine fh5s_get_ndims



  subroutine fh5s_get_dims(dims)
    implicit none

    integer, intent(OUT) :: dims(:)
    integer(HSIZE_T)     :: maxdimsc(7), dimsc(7)
    integer              :: ndims, i
    integer              :: hdferr

    call h5sget_simple_extent_ndims_f(dspcid, ndims, hdferr)

    call h5sget_simple_extent_dims_f(dspcid, dimsc, maxdimsc, hdferr)

    do i = 1, ndims
       dims(i) = int(dimsc(i))
    end do

    return
  end subroutine fh5s_get_dims



  subroutine fh5d_close(hdferr)
    implicit none

    integer, intent(OUT) :: hdferr

    call h5sclose_f(dspcid, hdferr)
    call h5dclose_f(dsetid, hdferr)

    return
  end subroutine fh5d_close



  subroutine fh5d_read_integer_array(buf_integer, hdferr)
    implicit none

    integer, intent(INOUT) :: buf_integer(*)
    integer, intent(OUT)   :: hdferr

    call h5dread_f(dsetid, H5T_NATIVE_INTEGER, buf_integer, dimsf, &
         hdferr, file_space_id=dspcid, mem_space_id=mspcid)

    return
  end subroutine fh5d_read_integer_array



  subroutine fh5d_read_real_array(buf_real, hdferr)
    implicit none

    real,    intent(INOUT) :: buf_real(*)
    integer, intent(OUT)   :: hdferr

    call h5dread_f(dsetid, H5T_NATIVE_REAL, buf_real, dimsf, &
         hdferr, file_space_id=dspcid, mem_space_id=mspcid)

    return
  end subroutine fh5d_read_real_array



  subroutine fh5d_read_character_array(buf_character, hdferr)
    implicit none

    character, intent(INOUT) :: buf_character(*)
    integer,   intent(OUT)   :: hdferr

    call h5dread_f(dsetid, H5T_NATIVE_CHARACTER, buf_character, dimsf, &
         hdferr, file_space_id=dspcid, mem_space_id=mspcid)

    return
  end subroutine fh5d_read_character_array



  subroutine fh5d_read_real8_array(buf_real8, hdferr)
    implicit none

    real(KIND=8), intent(INOUT) :: buf_real8(*)
    integer,      intent(OUT)   :: hdferr

    call h5dread_f(dsetid, H5T_NATIVE_DOUBLE, buf_real8, dimsf, &
         hdferr, file_space_id=dspcid, mem_space_id=mspcid)

    return
  end subroutine fh5d_read_real8_array



  subroutine fh5d_read_logical_array(buf_logical, hdferr)
    implicit none

    logical, intent(INOUT) :: buf_logical(*)
    integer, intent(OUT)   :: hdferr
    integer, allocatable   :: buf_integer(:)
    integer                :: i, size

    hdferr = 1

    size = 1
    do i = 1, ndimsf
       size = size * dimsf(i)
    enddo
    allocate (buf_integer(size))

    call h5dread_f(dsetid, H5T_NATIVE_INTEGER, buf_integer, dimsf, &
         hdferr, file_space_id=dspcid, mem_space_id=mspcid)

    ! converting integer to logical
    do i = 1, size
       if(buf_integer(i) == -1) then
          buf_logical(i) = .true.
       else
          buf_logical(i) = .false.
       endif
    enddo

    deallocate(buf_integer)

    return
  end subroutine fh5d_read_logical_array



  subroutine fh5d_read_integer_scalar(buf_integer, hdferr)
    implicit none

    integer, intent(INOUT) :: buf_integer
    integer, intent(OUT)   :: hdferr

    call h5dread_f(dsetid, H5T_NATIVE_INTEGER, buf_integer, dimsf, hdferr)

    return
  end subroutine fh5d_read_integer_scalar



  subroutine fh5d_read_real_scalar(buf_real, hdferr)
    implicit none

    real,    intent(INOUT) :: buf_real
    integer, intent(OUT)   :: hdferr

    call h5dread_f(dsetid, H5T_NATIVE_REAL, buf_real, dimsf, hdferr)

    return
  end subroutine fh5d_read_real_scalar



  subroutine fh5d_read_character_scalar(buf_character, hdferr)
    implicit none

    character, intent(INOUT) :: buf_character
    integer,   intent(OUT)   :: hdferr

    call h5dread_f(dsetid, H5T_NATIVE_CHARACTER, buf_character, dimsf, hdferr)

    return
  end subroutine fh5d_read_character_scalar



  subroutine fh5d_read_real8_scalar(buf_real8, hdferr)
    implicit none

    real(KIND=8), intent(INOUT) :: buf_real8
    integer, intent(OUT)        :: hdferr

    call h5dread_f(dsetid, H5T_NATIVE_DOUBLE, buf_real8, dimsf, hdferr)

    return
  end subroutine fh5d_read_real8_scalar



  subroutine fh5d_read_logical_scalar(buf_logical, hdferr)
    implicit none

    logical, intent(INOUT) :: buf_logical
    integer, intent(OUT)   :: hdferr
    integer                :: buf_integer

    hdferr = 1

    call h5dread_f(dsetid, H5T_NATIVE_INTEGER, buf_integer, dimsf, hdferr)

    buf_logical = .false.
    if(buf_integer == -1) buf_logical = .true.

    return
  end subroutine fh5d_read_logical_scalar



  subroutine fh5_prepare_read(dname, ndims, dims, hdferr, coords)
    implicit none

    character(*),      intent(IN)  :: dname
    integer,           intent(IN)  :: ndims
    integer,           intent(IN)  :: dims(:)
    integer,           intent(OUT) :: hdferr
    integer, optional, intent(IN)  :: coords(:)
    integer                        :: i, j, k, ndims_file, iop
    integer(HSIZE_T)               :: dims_file(7), maxdims_file(7)
    integer(HSIZE_T)               :: offset1d(1), count1d(1)
    integer(HSIZE_T)               :: offset2d(2), count2d(2)
    integer(HSIZE_T)               :: offset3d(3), count3d(3)

    dimsf  = 1
    ndimsf = ndims

    do i = 1, ndims
       dimsf(i) = dims(i)
    enddo

    ! opening the dataset from file

    call h5dopen_f(fileid, dname, dsetid, hdferr)
    if(dsetid < 0) then
       hdferr = int(dsetid)
       return
    end if

    ! getting the dataspace of the dataset (dimension, size, element type etc)

    call h5dget_space_f(dsetid, dspcid, hdferr)
    if(dspcid < 0) then
       hdferr = int(dspcid)
       return
    end if

    ! prepare the memspace to receive the data

    call h5screate_simple_f(ndims, dimsf, mspcid, hdferr)

    ! if coords is present select the points on dataspace

    if (present(coords)) then

       ! report an error if coords not dimensioned correctly

       if (size(coords) /= dims(ndims)) then
          write(*,*) "Error in fh5_prepare_read:"
          stop       "Invalid number of points selected."
       endif

       ! check the dataspace dimensions in the file

       call h5sget_simple_extent_ndims_f(dspcid, ndims_file, hdferr)
       call h5sget_simple_extent_dims_f(dspcid, dims_file, maxdims_file, hdferr)
       
       ! if the number of points we want is less than that in the file,
       ! select the points we want (else just read the entire dataspace)

       if (dims(ndims) < dims_file(ndims)) then

          if (ndims == 1) then

             offset1d = coords(1) - 1
             count1d  = 1
             iop      = H5S_SELECT_SET_F

             do j = 2, size(coords)
                if (coords(j) /= coords(j-1) + 1) then
                   call h5sselect_hyperslab_f(dspcid, iop, offset1d, &
                        count1d, hdferr)
                   offset1d = coords(j) - 1
                   count1d  = 1
                   iop      = H5S_SELECT_OR_F
                else
                   count1d  = count1d + 1
                endif
             enddo

             call h5sselect_hyperslab_f(dspcid, iop, offset1d, count1d, hdferr)

          else if (ndims == 2) then

             offset2d = (/ 0, coords(1)-1/)
             count2d  = (/ dims(1), 1 /)
             iop      = H5S_SELECT_SET_F
             
             do j = 2, size(coords)
                if (coords(j) /= coords(j-1) + 1) then
                   call h5sselect_hyperslab_f(dspcid, iop, offset2d, &
                        count2d, hdferr)
                   offset2d(2) = coords(j) - 1
                   count2d(2)  = 1
                   iop         = H5S_SELECT_OR_F
                else
                   count2d(2)  = count2d(2) + 1
                endif
             enddo

             call h5sselect_hyperslab_f(dspcid, iop, offset2d, count2d, hdferr)

          else if (ndims == 3) then

             offset3d = (/ 0, 0, coords(1)-1 /)
             count3d  = (/ dims(1), dims(2), 1 /)
             iop      = H5S_SELECT_SET_F

             do j = 2, size(coords)
                if (coords(j) /= coords(j-1) + 1) then
                   call h5sselect_hyperslab_f(dspcid, iop, offset3d, &
                        count3d, hdferr)
                   offset3d(3) = coords(j) - 1
                   count3d(3)  = 1
                   iop         = H5S_SELECT_OR_F
                else
                   count3d(3)  = count3d(3) + 1
                endif
             enddo
             
             call h5sselect_hyperslab_f(dspcid, iop, offset3d, count3d, hdferr)

          else

             stop 'ndims > 3 using partial I/O is not implemented'

          endif
       endif
    endif

    return
  end subroutine fh5_prepare_read



  subroutine fh5_close_read(hdferr)
    implicit none

    integer, intent(OUT) :: hdferr

    call h5sclose_f(mspcid, hdferr)
    call h5sclose_f(dspcid, hdferr)
    call h5dclose_f(dsetid, hdferr)

    return
  end subroutine fh5_close_read



end module hdf5_f2f
