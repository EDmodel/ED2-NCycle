program main
  use hdf5_utils
  implicit none

  integer, parameter :: iy1=1500
  integer, parameter :: iy2=1650
  character(len=256) :: indir='../run/hist/'
  character(len=256), parameter :: outname='dd2.lat2lon9-ts-M-extract.h5'
  integer, parameter :: nens=1
  integer, parameter :: nyear=iy2-iy1+1
  integer, parameter :: ntime=nyear*12

  integer, parameter :: npft = 20

  character(len=256) :: fname
  integer :: iy,im,itype,mylon,mylat,ndims,ilon,ilat
  integer, dimension(5) :: idims
  integer :: k,ierr,yind2
  logical :: lexist
  real, allocatable, dimension(:,:) :: ntext_soil
  logical, dimension(20) :: sclass
  integer :: lonind, latind,yind,hclass
  real, dimension(22), parameter :: slz=(/ 5.,4.5, 4., 3.5, 3., &
       2.5, 2., 1.75, 1.5, 1.25, &
       1., .85, .7, .6, .5,  &
       .4, .3, .2, .15, .1, &
       .05, 0./)
  real, dimension(100) :: area
  integer, dimension(100) :: paco_id, paco_n
  integer :: npatch,ipatch,icohort, iens
  character(len=256) :: cargv
  real :: nep_in,gpp_in,rh_in,precip_in,rshort_in
  real, dimension(npft,1) :: lai_in
  real, dimension(nens,ntime) :: nep_out,gpp_out,rh_out, precip_out, rshort_out
  real, dimension(nens,ntime,npft) :: lai_out
  character(len=10) :: post
  integer :: tind

  nep_out = 0.
  gpp_out = 0.
  rh_out = 0.
  precip_out = 0.
  rshort_out = 0.
  lai_out = 0.

  do iens = 1, nens
     write(post,'(i10)')iens
     post='r'//trim(adjustl(post))
     do iy = iy1,iy2
        yind = iy-iy1+1
        print*,iens,iy
        do im = 1, 12
           tind = (yind-1)*12 + im
           write(fname,'(2a,i4.4,a,i2.2,a)')  &
                trim(indir),  &
                'dd-nodist.lat3lon3-E-',iy,'-',im,'-00-000000-g01.h5'
           call shdf5_open(trim(fname),'R')
           ndims = 1
           idims(1) = 1
           call shdf5_irec(ndims,idims,'MMEAN_NEP',rvars=nep_in)
           call shdf5_irec(ndims,idims,'MMEAN_GPP',rvars=gpp_in)
           call shdf5_irec(ndims,idims,'MMEAN_RH',rvars=rh_in)
           call shdf5_irec(ndims,idims,'MMEAN_PCPG',rvars=precip_in)
           call shdf5_irec(ndims,idims,'MMEAN_RSHORT',rvars=rshort_in)

           ndims = 2
           idims(1) = npft
           idims(2) = 1
           call shdf5_irec(ndims,idims,'MMEAN_LAI_PFT',rvara=lai_in)

           nep_out(iens,tind) = nep_in
           gpp_out(iens,tind) = gpp_in
           rh_out(iens,tind) = rh_in
           precip_out(iens,tind) = precip_in
           rshort_out(iens,tind) = rshort_in
           lai_out(iens,tind,1:npft) = lai_in(1:npft,1)
           print*,iens,tind,precip_out(iens,tind)
           call shdf5_close()
        enddo
     enddo
  enddo

  call shdf5_open(trim(outname),'W',1)
  ndims = 2
  idims(1) = nens
  idims(2) = ntime
  call shdf5_orec(ndims,idims,'NEP',rvara=nep_out)
  call shdf5_orec(ndims,idims,'GPP',rvara=gpp_out)
  call shdf5_orec(ndims,idims,'RH',rvara=rh_out)
  call shdf5_orec(ndims,idims,'PRECIP',rvara=precip_out)
  call shdf5_orec(ndims,idims,'RSHORT',rvara=rshort_out)
  ndims = 3
  idims(1) = nens
  idims(2) = ntime
  idims(3) = npft
  call shdf5_orec(ndims,idims,'LAI_PFT',rvara=lai_out)
  call shdf5_close()

end program main
