!===============================================================================
! OLAM version 4.0

! Portions of this software are copied or derived from the RAMS software
! package.  The following copyright notice pertains to RAMS and its derivatives,
! including OLAM:  

   !----------------------------------------------------------------------------
   ! Copyright (C) 1991-2006  ; All Rights Reserved ; Colorado State University; 
   ! Colorado State University Research Foundation ; ATMET, LLC 

   ! This software is free software; you can redistribute it and/or modify it 
   ! under the terms of the GNU General Public License as published by the Free
   ! Software Foundation; either version 2 of the License, or (at your option)
   ! any later version. 

   ! This software is distributed in the hope that it will be useful, but
   ! WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
   ! or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
   ! for more details.
 
   ! You should have received a copy of the GNU General Public License along
   ! with this program; if not, write to the Free Software Foundation, Inc.,
   ! 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA 
   ! (http://www.gnu.org/licenses/gpl.html) 
   !----------------------------------------------------------------------------

! OLAM was developed at Duke University and the University of Miami, Florida. 
! For additional information, including published references, please contact
! the software authors, Robert L. Walko (rwalko@rsmas.miami.edu)
! or Roni Avissar (ravissar@rsmas.miami.edu).
!===============================================================================
Module hdf5_utils

integer, parameter :: icompress = 0

Contains

subroutine shdf5_open(locfn, access, idelete)

  use hdf5_f2f
  implicit none

  character(*),      intent(in) :: locfn   ! file name
  character(*),      intent(in) :: access  ! File access ('R','W','RW')
  integer, optional, intent(in) :: idelete ! If W, overwrite file if exists?
                                           !    1=yes, 0=no

  integer                       :: hdferr  ! Error flag
  integer                       :: iaccess ! int access flag
  logical                       :: exists  ! File existence

! Check for existence of RAMS file.

  inquire(file=locfn, exist=exists)

! Create a new file or open an existing RAMS file.

  if (access(1:1) == 'R') then

     if (.not. exists) then
        print*, 'shdf5_open:'
        print*, '   Attempt to open a file for reading that does not exist.'
        print*, '   Filename: ', trim(locfn)
        stop    'shdf5_open: no file'
     else
        if (access == 'R ') iaccess = 1
        if (access == 'RW') iaccess = 2
        call fh5f_open(locfn, iaccess, hdferr)
        
        if (hdferr < 0) then
           print*, 'shdf5_open:'
           print*, '   Error opening hdf5 file - error -', hdferr
           print*, '   Filename: ',trim(locfn)
           stop    'shdf5_open: open error'      
        endif
     endif
     
  elseif (access(1:1) == 'W') then

     if (.not. exists) then
        iaccess = 2
        call fh5f_create(locfn, iaccess, hdferr)
     else
        if (.not. present(idelete) ) then
           print*, 'shdf5_open: idelete not specified when access=W'
           stop    'shdf5_open: no idelete'
        endif
      
        if (idelete == 0) then
           print*, 'In shdf5_open:'
           print*, '   Attempt to open an existing file for writing,'
           print*, '      but overwrite is disabled. idelete=', idelete
           print*, '   Filename: ', trim(locfn)
           stop    'shdf5_open'
        else
           iaccess = 1
           call fh5f_create(locfn, iaccess, hdferr)
        endif

     endif
   
     if(hdferr < 0) then
        print*, 'HDF5 file create failed:', hdferr
        print*, 'file name:', trim(locfn), ' ', trim(access), idelete
        stop    'shdf5_open: bad create'
     endif
  endif

  return
end subroutine shdf5_open

!===============================================================================

subroutine shdf5_info(dsetname, ndims, dims)
  use hdf5_f2f

  implicit none

  character(*), intent(in)    :: dsetname ! Dataset name
  integer,      intent(inout) :: dims(:)
  integer,      intent(inout) :: ndims    ! Dataset rank (in file)
  integer                     :: hdferr   ! Error flag

! Open the dataset.

  call fh5d_open(dsetname, hdferr)

  if (hdferr < 0) then
     print*, 'In shdf5_info:'
     print*, 'Variable ', trim(dsetname), ' is not in the currently opened hdf5 file'
     ndims   = 0
     dims(1) = 0
     return
  endif

! Get dataset's dimensions

  call fh5s_get_ndims(ndims)
  call fh5s_get_dims(dims)

  call fh5d_close(hdferr)

  return
end subroutine shdf5_info

!===============================================================================

subroutine shdf5_orec(ndims,dims,dsetname,ivara,rvara,cvara,dvara,lvara, &
                                          ivars,rvars,cvars,dvars,lvars, &
                                          nglobe, lpoints, gpoints)
  use hdf5_f2f
    
  implicit none

  character(*), intent(in) :: dsetname ! Variable label
  integer,      intent(in) :: ndims    ! Number of dimensions or rank
  integer,      intent(in) :: dims(:)  ! Dataset dimensions.

! Array and scalar arguments for different types. Only specify one in each call
  integer,      intent(in), optional :: ivara(*), ivars
  real,         intent(in), optional :: rvara(*), rvars
  character,    intent(in), optional :: cvara(*), cvars
  real(kind=8), intent(in), optional :: dvara(*), dvars
  logical,      intent(in), optional :: lvara(*), lvars

! Optional arrays to determine cells for partial/parallel IO
  integer,      intent(in), optional :: lpoints(:), gpoints(:), nglobe

! Local variables
  integer :: hdferr  ! Error flag

! Check dimensions and set compression chunk size

  if (ndims <=0 .or. minval(dims(1:ndims)) <=0) then
     print*, 'Dimension error in shdf5_orec:', ndims, dims(1:ndims)
     stop    'shdf5_orec: bad dims'
  endif
     
! Prepare memory and options for the write

  call fh5_prepare_write(ndims, dims, hdferr, icompress, &
       mcoords=lpoints, fcoords=gpoints, ifsize=nglobe)

  if (hdferr /= 0) then
     print*, "shdf5_orec: can't prepare requested field:", trim(dsetname)
     return
  endif

! Write the dataset.

      if (present(ivars)) then ; call fh5_write(ivars, dsetname, hdferr)
  elseif (present(rvars)) then ; call fh5_write(rvars, dsetname, hdferr)
  elseif (present(cvars)) then ; call fh5_write(cvars, dsetname, hdferr)
  elseif (present(dvars)) then ; call fh5_write(dvars, dsetname, hdferr)
  elseif (present(lvars)) then ; call fh5_write(lvars, dsetname, hdferr)
  elseif (present(ivara)) then ; call fh5_write(ivara, dsetname, hdferr)
  elseif (present(rvara)) then ; call fh5_write(rvara, dsetname, hdferr)
  elseif (present(cvara)) then ; call fh5_write(cvara, dsetname, hdferr)
  elseif (present(dvara)) then ; call fh5_write(dvara, dsetname, hdferr)
  elseif (present(lvara)) then ; call fh5_write(lvara, dsetname, hdferr)
  else
     print*, 'Incorrect or missing data field argument in shdf5_orec'
     stop    'shdf5_orec: bad data field'
  endif

  if (hdferr /= 0) then
     print*, 'In shdf5_orec: hdf5 write error =', hdferr
     stop    'shdf5_orec: hdf5 write error'
  endif

! Close the dataset, the dataspace for the dataset, and the dataspace properties.

  call fh5_close_write(hdferr)

  return
end subroutine

!===============================================================================

subroutine shdf5_irec(ndims,dims,dsetname,ivara,rvara,cvara,dvara,lvara  &
                                         ,ivars,rvars,cvars,dvars,lvars  &
                                         ,points)
  use hdf5_f2f
  implicit none

  character(*), intent(IN) :: dsetname ! Dataset name
  integer,      intent(IN) :: ndims    ! Number of dimensions or rank
  integer,      intent(IN) :: dims(:)  ! Dataset dimensions

! Array and scalar arguments for different types. Only specify one in each call.
  integer,      intent(OUT), optional :: ivara(*), ivars
  real,         intent(OUT), optional :: rvara(*), rvars
  character,    intent(OUT), optional :: cvara(*), cvars
  real(kind=8), intent(OUT), optional :: dvara(*), dvars
  logical,      intent(OUT), optional :: lvara(*), lvars

! Optional arrays to determine cells for partial/parallel IO
  integer,      intent(IN),  optional :: points(:)

! Local variables
  integer :: hdferr  ! Error flag

! Check dimensions

  if (ndims <=0 .or. minval(dims(1:ndims)) <=0) then
     print*, 'Dimension error in shdf5_irec:', ndims, dims(1:ndims)
     stop    'shdf5_irec: bad dims'
  endif
    
! Prepare file and memory space for the read

  call fh5_prepare_read(dsetname, ndims, dims, hdferr, coords=points)
  if (hdferr < 0) then
     print*,'shdf5_irec: can''t prepare requested field:',trim(dsetname)
     return
  endif

! Read data from hyperslab in the file into the hyperslab in memory.

      if (present(ivars)) then ; call fh5d_read(ivars, hdferr)
  elseif (present(rvars)) then ; call fh5d_read(rvars, hdferr)
  elseif (present(cvars)) then ; call fh5d_read(cvars, hdferr)
  elseif (present(dvars)) then ; call fh5d_read(dvars, hdferr)
  elseif (present(lvars)) then ; call fh5d_read(lvars, hdferr)
  elseif (present(ivara)) then ; call fh5d_read(ivara, hdferr)
  elseif (present(rvara)) then ; call fh5d_read(rvara, hdferr)
  elseif (present(cvara)) then ; call fh5d_read(cvara, hdferr)
  elseif (present(dvara)) then ; call fh5d_read(dvara, hdferr)
  elseif (present(lvara)) then ; call fh5d_read(lvara, hdferr)
  else
     print*,'Incorrect or missing data field argument in shdf5_irec'
     print*, 'field = ', dsetname
     stop    'shdf5_irec: bad data field'
  endif
  
  if (hdferr /= 0) then
     print*, 'shdf5_irec: call fh5d_read: hdf5 error =', hdferr
     print*, 'Error reading ', trim(dsetname)
     print*, 'ndims = ', ndims
     print*, 'dims  = ', dims(1:ndims)
     stop
  endif

! Close the dataset, the dataspace for the dataset, and the memory space.

  call fh5_close_read(hdferr)

  return
end subroutine

!===============================================================================

subroutine shdf5_close()
  use hdf5_f2f
  implicit none

  integer :: hdferr  ! Error flags

! Close hdf file.

  call fh5f_close(hdferr)

  return
end  subroutine

!===============================================================================

subroutine shdf5_io(action,ndims,dims,dsetname,ivara,rvara,cvara,dvara,lvara  &
                                              ,ivars,rvars,cvars,dvars,lvars)
  use hdf5_f2f
  implicit none

  character(*), intent(in)              :: dsetname, action
  integer,      intent(in)              :: ndims,    dims(:)
  integer,      intent(inout), optional :: ivara(*), ivars
  real,         intent(inout), optional :: rvara(*), rvars
  character,    intent(inout), optional :: cvara(*), cvars
  real(kind=8), intent(inout), optional :: dvara(*), dvars
  logical,      intent(inout), optional :: lvara(*), lvars
 
  ! THIS ROUTINE CALLS SHDF5_IREC OR SHDF5_OREC TO READ OR WRITE A VARIABLE
  ! DEPENDING ON WHETHER 'ACTION' EQUALS 'READ' OR 'WRITE'

  if (action == 'READ') then
     
     call shdf5_irec(ndims,dims,dsetname,ivara,rvara,cvara,dvara,lvara,  &
                     ivars,rvars,cvars,dvars,lvars)
     
  elseif (action == 'WRITE') then
     
     call shdf5_orec(ndims,dims,dsetname,ivara,rvara,cvara,dvara,lvara,  &
                     ivars,rvars,cvars,dvars,lvars)
     
  else
     
     print *, "Illegal action in shdf5_io."
     print *, "Action should be 'READ' or 'WRITE'"
     stop     "Ending model run"

  endif
  
end subroutine shdf5_io

end module
