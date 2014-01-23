program main

  use hdf5, only: h5open_f, h5close_f
  use hdf5_utils

  implicit none

  integer, parameter :: imonth1 = 1
  integer, parameter :: idate1 = 1
  integer, parameter :: iyear1 = 2000
  integer, parameter :: imonthf = 11
  integer, parameter :: idatef = 30
  integer, parameter :: iyearf = 2050

  character(len=256), parameter ::   &
       afilpref='/scratch/gpfs/att/PFT21/PFT21_50newallom'
  character(len=256), parameter :: outname='PFT21_50newallom.txt'

  integer :: ndims, hdferr, iyear, m1, m2, imonth, d1, d2, idate
  character(len=256) :: fname
  integer, dimension(3) :: idims
  real, dimension(:,:), allocatable :: lai_pft
  real, dimension(:), allocatable :: nep
  real, dimension(:,:), allocatable :: agb_pft
  real, dimension(:), allocatable :: gpp
  real, dimension(:), allocatable :: transp
  real, dimension(:), allocatable :: evap
  real, dimension(:), allocatable :: nppdaily
 

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
        d1 = 00
        d2 = 00

        ! Not all months have 31 days.
       ! if(imonth == 4 .or. imonth == 6 .or. imonth == 9 .or. imonth == 11)then
        !   d2 = 30
       ! elseif(imonth == 2)then
        !   d2 = 28
         !  if(mod(iyear,4) == 4)d2 = 29
       ! endif

        ! Special case: what if the simulation does not start on the 
        ! first day of the month?
       ! if(iyear == iyear1 .and. imonth == imonth1) d1 = idate1

        ! Special case: what if the simulation does not end on the 
        ! last day of the month?
       ! if(iyear == iyearf .and. imonth == imonthf) d2 = idatef

        ! Loop over days
        do idate = d1, d2

           ! Create input file name
           write(fname,'(a,a,i4.4,a,i2.2,a,i2.2,a)')trim(afilpref),'-E-',iyear,'-',imonth,'-',idate,'-000000-g01.h5'
           print*,trim(fname)

           ! Open input file
           call shdf5_open(trim(fname),'R')

           ! Find array dimensions
           call shdf5_info('MMEAN_LAI_PFT',ndims,idims)
           allocate(lai_pft(idims(1),idims(2)))

           ! Get data
           call shdf5_irec(ndims,idims,'MMEAN_LAI_PFT',rvara=lai_pft)

           call shdf5_info('MMEAN_NEP',ndims,idims)
           allocate(nep(idims(1)))
           call shdf5_irec(ndims,idims,'MMEAN_NEP',rvara=nep)

           call shdf5_info('AGB_PFT',ndims,idims)
           allocate(agb_pft(idims(1),idims(2)))
           call shdf5_irec(ndims,idims,'AGB_PFT',rvara=agb_pft) 

           call shdf5_info('MMEAN_GPP',ndims,idims)
           allocate(gpp(idims(1)))
           call shdf5_irec(ndims,idims,'MMEAN_GPP',rvara=gpp)

           call shdf5_info('MMEAN_TRANSP',ndims,idims)
           allocate(transp(idims(1)))
           call shdf5_irec(ndims,idims,'MMEAN_TRANSP',rvara=transp)

           call shdf5_info('MMEAN_EVAP',ndims,idims)
           allocate(evap(idims(1)))
           call shdf5_irec(ndims,idims,'MMEAN_EVAP',rvara=evap)

           call shdf5_info('MMEAN_NPPDAILY',ndims,idims)
           allocate(nppdaily(idims(1)))
           call shdf5_irec(ndims,idims,'MMEAN_NPPDAILY',rvara=nppdaily)



           ! Write to output file (entry for PFT number 8, grid cell number 1)
           write(12,'(1e14.6,1e14.6,1e14.6,1e14.6,1e14.6,1e14.6,1e14.6,1e14.6,1e14.6)')lai_pft(21,1),lai_pft(22,1),nep(1),agb_pft(21,1),agb_pft(22,1),gpp(1),transp(1),evap(1),nppdaily(1)

           ! Deallocate memory
           deallocate(lai_pft)
           deallocate(nep)
           deallocate(agb_pft)
           deallocate(gpp)
           deallocate(transp)
           deallocate(evap)
           deallocate(nppdaily)

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
