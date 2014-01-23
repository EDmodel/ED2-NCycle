program main
  use hdf5_utils
  implicit none

  character(len=256) :: inprefix =   &
       '/scratch/lustre/dmedvigy/ed2_mods/ED.africa/dd2/dd2.'
  character(len=256) :: outname = 'test_map.h5'
  integer, parameter :: nsites = 30
  integer, parameter, dimension(nsites) :: lat_list = (/   &
       0, 0, 0, 0, 0,  &
       1, 1, 1, 1, 1,  &
       2, 2, 2, 2, 2,  &
       3, 3, 3, 3, 3,  &
       4, 4, 4, 4, 4,  &
       5, 5, 5, 5, 5  &
       /)
  integer, parameter, dimension(nsites) :: lon_list = (/ 1, 3, 5, 7, 9,  &
       1, 3, 5, 7, 9,  &
       1, 3, 5, 7, 9,  &
       1, 3, 5, 7, 9,  &
       1, 3, 5, 7, 9,  &
       1, 3, 5, 7, 9   /)
  integer, parameter :: iyear=1650
  integer, parameter :: imonth=01
  integer, parameter :: npft=20

  !====================================================

  real, dimension(nsites) :: my_lat, my_lon
  character(len=256) :: fname
  character(len=10) :: post1, post2
  real, dimension(npft,1) :: agb_pft_in
  real, dimension(nsites,npft) :: agb_pft_out
  integer, dimension(5) :: idims
  integer :: ndims
  integer :: isite, ipft
  real :: grid_res
  logical :: fexists

  !====================================================

  agb_pft_out = 0.

  do isite = 1, nsites
     write(post1,'(i10)')lat_list(isite)
     write(post2,'(i10)')lon_list(isite)
     write(fname,'(6a,i4.4,a,i2.2,a)')  &
          trim(inprefix), 'lat', trim(adjustl(post1)),  &
          'lon', trim(adjustl(post2)), '-E-',  &
          iyear,'-',imonth,'-00-000000-g01.h5'
     inquire(file=trim(fname),exist=fexists)
     if(.not.fexists)cycle
     call shdf5_open(trim(fname),'R')
     ndims = 2
     idims(1) = npft
     idims(2) = 1
     call shdf5_irec(ndims,idims,'AGB_PFT',rvara=agb_pft_in)
     do ipft = 1, npft
        agb_pft_out(isite,ipft) = agb_pft_in(ipft,1)
     enddo
     ndims = 1
     idims(1) = 1
     call shdf5_irec(ndims,idims,'LONGITUDE',rvars=my_lon(isite))
     call shdf5_irec(ndims,idims,'LATITUDE',rvars=my_lat(isite))
     call shdf5_close()
  enddo

  call shdf5_open(trim(outname),'W',1)
  ndims = 2
  idims(1) = nsites
  idims(2) = npft
  call shdf5_orec(ndims,idims,'AGB_PFT',rvara=agb_pft_out)
  ndims = 1
  idims(1) = nsites
  call shdf5_orec(ndims,idims,'LON',rvara=my_lon)
  call shdf5_orec(ndims,idims,'LAT',rvara=my_lat)
  ndims = 1
  idims(1) = 1
  grid_res=2.5
  call shdf5_orec(ndims,idims,'GRID_RES',rvars=grid_res)
  call shdf5_close()

end program main
