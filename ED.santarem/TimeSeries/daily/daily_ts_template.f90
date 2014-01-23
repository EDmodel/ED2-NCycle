program main

  use hdf5, only: h5open_f, h5close_f
  use hdf5_utils

  implicit none

!JL CHANGE DATES
  integer, parameter :: imonth1 = 01
  integer, parameter :: idate1 = 01
  integer, parameter :: iyear1 = 2002
  integer, parameter :: imonthf = 01
  integer, parameter :: idatef = 02
  integer, parameter :: iyearf = 2003

  character(len=256), parameter ::   &
       afilpref='/home/jhlevy/ED2/ED.santarem/run/PFT3'
  character(len=256), parameter :: outname='/home/jhlevy/ED2/ED.santarem/run/PFT3.txt'

  integer :: ndims, hdferr, iyear, m1, m2, imonth, d1, d2, idate
  character(len=256) :: fname
  integer, dimension(3) :: idims
  real, dimension(:,:), allocatable :: lai_pft
 
!JL!   
  real, dimension(:,:), allocatable :: agb_pft
  real, dimension(:), allocatable :: nstorage
  real, dimension(:), allocatable :: mineralized_soil_N
  real, dimension(:), allocatable :: actual_nitrogen_uptake

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

           ! Find array dimensions
            call shdf5_info('LAI_PFT',ndims,idims)
            allocate(lai_pft(idims(1),idims(2)))

          !  Get data
           call shdf5_irec(ndims,idims,'LAI_PFT',rvara=lai_pft)

! JL
           
           call shdf5_info('AGB_PFT',ndims,idims)
           allocate(agb_pft(idims(1),idims(2)))
           call shdf5_irec(ndims,idims,'AGB_PFT',rvara=agb_pft) 

            call shdf5_info('NSTORAGE',ndims,idims)
            allocate(nstorage(idims(1)))
            call shdf5_irec(ndims,idims,'NSTORAGE',rvara=nstorage)

           call shdf5_info('MINERALIZED_SOIL_N',ndims,idims)
            allocate(mineralized_soil_N(idims(1)))
            call shdf5_irec(ndims,idims,'MINERALIZED_SOIL_N',rvara=mineralized_soil_N)

            call shdf5_info('ACTUAL_NITROGEN_UPTAKE',ndims,idims)
            allocate(actual_nitrogen_uptake(idims(1)))
            call shdf5_irec(ndims,idims,'ACTUAL_NITROGEN_UPTAKE',rvara=actual_nitrogen_uptake)


           ! Write to output file (entry for PFT number 8, grid cell number 1)
 !JL          write(12,'(1e14.6)')lai_pft(8,1)
           write(12,'(1e14.6,1e14.6,1e14.6,1e14.6,1e14.6,1e14.6,1e14.6)')lai_pft(3,1),lai_pft(31,1),agb_pft(3,1),agb_pft(31,1),nstorage(1), mineralized_soil_N(1),actual_nitrogen_uptake(1)
          ! write(12,'(1e14.6)')agb_pft(:,1) will give a string of all agb_pft values.. so at each timestep this will give 35 values in an array, each pos          !ition (1-35) will  represent each pft in order
 !JL!          ! Deallocate memory
           deallocate(agb_pft)
           deallocate(nstorage)
           deallocate(mineralized_soil_N)
           deallocate(actual_nitrogen_uptake)
           deallocate(lai_pft)
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

! go to timeseries/daily folder
!1) make
!2)  ./daily_ts_template
!3) the txt file will be created in the run folder or wherever the ouput is specified above
!4) can open it up in emacs
