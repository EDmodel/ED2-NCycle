program main
  use hdf5_utils
  implicit none

  integer, parameter :: iy1=1500
  integer, parameter :: iy2=1650
  character(len=256) :: indir='../run/hist/'
  character(len=256), parameter :: outname='dd2.lat2lon9-ts-Y-extract.h5'
  integer, parameter :: nens=1
  integer, parameter :: nyear=iy2-iy1+1

  real, dimension(12,nyear) :: gpp_out,nep_out,plresp_out,rh_out,fsw_out,soilw_out,cantmp_out, rshort_out, lai_pine_out, lai_oak_out, lai_out,gpp_pine_out,gpp_oak_out,precip_out, transp_out, evap_out
  real, dimension(12,nyear,7) :: gpp_hclass_out
  real, dimension(nyear) :: agpp_out,anep_out,aplresp_out,arh_out,afsw_out,asoilw_out,acantmp_out,arshort_out,agpp_pine_out,agpp_oak_out,aprecip_out, atransp_out, aevap_out
  character(len=256) :: fname
  integer :: iy,im,itype,mylon,mylat,ndims,ilon,ilat
  integer, dimension(5) :: idims
  integer :: k,ierr,yind2
  real :: mmean_gpp,mmean_nep,mmean_plresp,mmean_rh,mmean_fsw,mmean_cantmp,mmean_rshort,mmean_pcpg, mmean_evap, mmean_transp
  real, dimension(21) :: mmean_soilw
  real, dimension(15) :: mmean_lai
  logical :: lexist
  real, allocatable, dimension(:,:) :: ntext_soil
  logical, dimension(20) :: sclass
  integer :: lonind, latind,yind,hclass
  real, dimension(22), parameter :: slz=(/ 5.,4.5, 4., 3.5, 3., &
       2.5, 2., 1.75, 1.5, 1.25, &
       1., .85, .7, .6, .5,  &
       .4, .3, .2, .15, .1, &
       .05, 0./)
  integer, parameter :: max_patches=100
  real, dimension(1000) :: mmean_gpp_co,nplant,hite
  integer, dimension(1000) :: pft
  real, dimension(max_patches) :: area
  integer, dimension(100) :: paco_id, paco_n
  integer :: npatch,ipatch,icohort, iens
  character(len=256) :: cargv
  real, dimension(:,:,:), allocatable :: agb_in
  real, dimension(nens,nyear) :: agb_out
  character(len=10) :: post
  integer :: npatches, ip
  real, dimension(max_patches) :: stsc, ssc, fsc
  real, dimension(nens,nyear) :: soilc_out

  iy = iy1
  iens = 1
  write(fname,'(2a,i4.4,a)')  &
       trim(indir),  &
       'dd2.lat2lon9-Y-',iy,'-00-00-000000-g01.h5'
  call shdf5_open(trim(fname),'R')
  call shdf5_info('AGB',ndims,idims)
  allocate(agb_in(idims(1),idims(2),idims(3)))
  call shdf5_close()

  agb_out = 0.
  soilc_out = 0.

  do iens = 1, nens
     do iy = iy1,iy2
        yind = iy-iy1+1
        print*,iens,iy
        write(fname,'(2a,i4.4,a)')  &
             trim(indir),  &
             'dd2.lat2lon9-Y-',iy,'-00-00-000000-g01.h5'
        call shdf5_open(trim(fname),'R')
        call shdf5_info('AGB',ndims,idims)
        call shdf5_irec(ndims,idims,'AGB',rvara=agb_in)
        agb_out(iens,yind) = sum(agb_in)

        call shdf5_info('AREA',ndims,idims)
        call shdf5_irec(ndims,idims,'AREA',rvara=area)
        npatches = idims(1)
        call shdf5_irec(ndims,idims,'STRUCTURAL_SOIL_C',rvara=stsc)
        call shdf5_irec(ndims,idims,'SLOW_SOIL_C',rvara=ssc)
        call shdf5_irec(ndims,idims,'FAST_SOIL_C',rvara=fsc)
        soilc_out(iens,yind) = 0.0
        do ip = 1, npatches
           soilc_out(iens,yind) = soilc_out(iens,yind) + area(ip) * &
                (stsc(ip)+ssc(ip)+fsc(ip))
        enddo

        print*,iens,yind,agb_out(iens,yind)
        call shdf5_close()
     enddo
  enddo

  call shdf5_open(trim(outname),'W',1)
  ndims = 2
  idims(1) = nens
  idims(2) = nyear
  call shdf5_orec(ndims,idims,'AGB',rvara=agb_out)
  call shdf5_orec(ndims,idims,'SOILC',rvara=soilc_out)
  call shdf5_close()

end program main
