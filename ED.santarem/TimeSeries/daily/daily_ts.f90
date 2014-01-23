program main

  use hdf5, only: h5open_f, h5close_f
  use hdf5_utils

  implicit none

  integer, parameter :: imonth1 = 2
  integer, parameter :: idate1 = 3
  integer, parameter :: iyear1 = 2012
  integer, parameter :: imonthf = 5
  integer, parameter :: idatef = 17
  integer, parameter :: iyearf = 2012

  character(len=256), parameter ::   &
!       afilpref='../run/hist/ht4-plim-d-v3-e'
       afilpref='../run/output/b'
  character(len=256), parameter :: outname='ainit_sat_daily_ts.txt'
  character(len=256), parameter :: outname2='control-plim-a-v3_yearly_ts.txt'
  character(len=256), parameter :: outname3='control-plim-a-v3_seascy_ts.txt'

  integer, parameter :: nens=10

  integer, parameter :: nyears=iyearf-iyear1+1
  integer :: iyear, m1, m2, imonth,ipoly,k
  character(len=256) :: fname
  integer :: ndims,iy
  integer, dimension(4) :: idims
  real :: ra
  integer :: nzg, npoly,nmonth,d1,d2
  real :: ea_yearly_gpp, ea_sd_gpp, ea_cumul_gpp, ea_sd_cumul_gpp
  !=====================================================================

  real, allocatable, dimension(:,:) :: yearly_gpp, cumul_gpp
  real, dimension(:), allocatable :: slz
  real, dimension(:), allocatable :: mmean_leaf_resp,mmean_root_resp,  &
       mmean_growth_resp,mmean_storage_resp,mmean_vleaf_resp,mmean_plresp
  real, allocatable, dimension(:) :: mmean_gpp, mmean_rh, mmean_nep,  &
       mmean_sens, mmean_rshort, mmean_precip, mmean_temp, mmean_rhum, &
       mmean_can_temp
  real, allocatable, dimension(:,:) :: mmean_lai, mmean_agb
  real, allocatable, dimension(:) :: ymean_rh, ymean_nep, &
       ymean_leaf_resp, ymean_storage_resp, ymean_growth_resp,   &
       ymean_root_resp,ymean_vleaf_resp,ymean_plresp, ymean_sens,   &
       ymean_rshort, ymean_temp, ymean_precip, ymean_rhum
  real, dimension(:,:), allocatable :: mmean_soil_water
  real, dimension(:), allocatable :: tot_soil_water, mmean_evap,   &
       mmean_transp,mmean_runoff
  integer :: iens
  real, allocatable, dimension(:,:) :: ea_mmean_soil_water
  real, allocatable, dimension(:) :: ea_mmean_leaf_resp, ea_mmean_root_resp, &
       ea_mmean_growth_resp, ea_mmean_storage_resp, ea_mmean_vleaf_resp,   &
       ea_mmean_plresp, ea_mmean_gpp, ea_mmean_rh, ea_mmean_nep, &
       ea_ymean_rh, ea_ymean_nep, ea_ymean_leaf_resp, ea_ymean_storage_resp, &
       ea_ymean_growth_resp, ea_ymean_root_resp, ea_ymean_vleaf_resp, &
       ea_mmean_evap, ea_mmean_transp, ea_mmean_runoff, ea_tot_soil_water, &
       ea_ymean_plresp, ea_mmean_sens, ea_mmean_precip, ea_mmean_rshort,   &
       ea_mmean_temp, ea_ymean_sens, ea_ymean_precip, ea_ymean_rshort,   &
       ea_ymean_temp, ea_mmean_rhum, ea_ymean_rhum, ea_msd_gpp, ea_mmean_lai, &
       ea_mmean_agb, ea_mmean_can_temp
  integer :: hdferr,idate,nday,id1
  real, external :: rhovsl
  real :: ea_seascy_gpp, ea_seascy_sd_gpp,ea_seascy_lai, ea_seascy_sd_lai
  real :: ea_seascy_transp, ea_seascy_et, ea_seascy_sens
  real :: ea_seascy_sd_transp, ea_seascy_sd_et, ea_seascy_sd_sens
  real, dimension(12,nens) :: seascy_gpp, seascy_transp, seascy_et,   &
       seascy_sens, seascy_lai
  real, dimension(21,1) :: ch4_conc,soil_water
  !=====================================================================

  ! Initialize HDF5 environment
  call h5open_f(hdferr)

  ! Open output file
  open(12,file=trim(outname),form='formatted',status='replace')

  ! Loop over years
  do iyear = iyear1, iyearf

     ! First and last month to process
     m1 = 1
     m2 = 12

     ! Special case: what if simulation does not begin in January?
     if(iyear == iyear1)m1 = imonth1

     ! Special case: what if simulation does not end in December?
     if(iyear == iyearf)m2 = imonthf

     ! Loop over months
     do imonth = m1, m2

        ! Set the start and end day for this month
        d1 = 1
        d2 = 31

        ! Not all months have 31 days.
        if(imonth == 4 .or. imonth == 6 .or. imonth == 9 .or. imonth == 11)then
           d2 = 30
        elseif(imonth == 2)then
           d2 = 28
           if(mod(iyear,4) == 4)d2 = 29
        endif

        ! Special case: what if the simulation does not start on the 
        ! first day of the month?
        if(iyear == iyear1 .and. imonth == imonth1) d1 = idate1

        ! Special case: what if the simulation does not end on the 
        ! last day of the month?
        if(iyear == iyearf .and. imonth == imonthf) d2 = idatef

        ! Loop over days
        do idate = d1, d2

           ! Create input file name
           write(fname,'(a,a,i4.4,a,i2.2,a,i2.2,a)')trim(afilpref),'-D-',iyear,'-',imonth,'-',idate,'-000000-g01.h5'
           print*,trim(fname)

           ! Open input file
           call shdf5_open(trim(fname),'R')

           ! Get data
           call shdf5_irec(ndims,idims,'DMEAN_SOIL_CH4_CONC_WTR',rvara=ch4_conc)
           call shdf5_irec(ndims,idims,'DMEAN_SOIL_WATER',rvara=soil_water)

           ! Write to output file
           write(12,'(21e14.6)')ch4_conc/soil_water*1000

           ! Close input file
           call shdf5_close()
        enddo
     enddo
  enddo

  ! Close output file
  close(12)

  ! Close HDF5 environment
  call h5close_f(hdferr)

end program main
