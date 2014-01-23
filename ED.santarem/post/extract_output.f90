program main
  use hdf5_utils
  implicit none

  integer, parameter :: iy1=1500
  integer, parameter :: iy2=1540
  character(len=256) :: indir='/scratch/lustre/dmedvigy/ed2_mods/spinuptest/r1/'
!  character(len=256) :: outname
  character(len=256) :: outname='ts-Y-extract.h5'
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
  real, dimension(1000) :: mmean_gpp_co,nplant,hite
  integer, dimension(1000) :: pft
  real, dimension(100) :: area
  integer, dimension(100) :: paco_id, paco_n
  integer :: npatch,ipatch,icohort
  character(len=256) :: cargv

!  call get_command_argument(1,cargv)
!  indir = trim(indir)//trim(cargv)
!  outname = trim(cargv)//'-flux-extract.h5'

  gpp_out = 0.
  nep_out = 0.
  plresp_out = 0.
  rh_out = 0.
  fsw_out = 0.
  evap_out = 0.
  transp_out = 0.
  soilw_out = 0.
  cantmp_out = 0.
  rshort_out = 0.
  lai_pine_out = 0.
  lai_oak_out = 0.
  lai_out = 0.
  gpp_pine_out = 0.
  gpp_oak_out = 0.
  precip_out = 0.

  gpp_hclass_out = 0.

  agpp_out = 0.
  anep_out = 0.
  aplresp_out = 0.
  arh_out = 0.
  atransp_out = 0.
  aevap_out = 0.
  afsw_out = 0.
  asoilw_out = 0.
  acantmp_out = 0.
  arshort_out = 0.
  agpp_pine_out = 0.
  agpp_oak_out = 0.
  aprecip_out = 0.

  do iy = iy1,iy2
     yind = iy-iy1+1
     print*,iy
     do im = 1, 12
        write(fname,'(2a,i4.4,a,i2.2,a)')  &
             trim(indir),  &
             '-E-',iy,'-',im,'-00-000000-g01.h5'
        call shdf5_open(trim(fname),'R')
        call shdf5_info('MMEAN_GPP',ndims,idims)
        call shdf5_irec(ndims,idims,'MMEAN_GPP',rvars=mmean_gpp)
        call shdf5_irec(ndims,idims,'MMEAN_NEP',rvars=mmean_nep)
        call shdf5_irec(ndims,idims,'MMEAN_PLRESP',rvars=mmean_plresp)
        call shdf5_irec(ndims,idims,'MMEAN_RH',rvars=mmean_rh)
        call shdf5_irec(ndims,idims,'MMEAN_FSW',rvars=mmean_fsw)
        call shdf5_irec(ndims,idims,'MMEAN_EVAP',rvars=mmean_evap)
        call shdf5_irec(ndims,idims,'MMEAN_TRANSP',rvars=mmean_transp)
        call shdf5_irec(ndims,idims,'MMEAN_CAN_TEMP',rvars=mmean_cantmp)
        call shdf5_irec(ndims,idims,'MMEAN_RSHORT',rvars=mmean_rshort)
        call shdf5_irec(ndims,idims,'MMEAN_PCPG',rvars=mmean_pcpg)
        call shdf5_info('MMEAN_SOIL_WATER',ndims,idims)
        call shdf5_irec(ndims,idims,'MMEAN_SOIL_WATER',rvara=mmean_soilw)
        call shdf5_info('MMEAN_LAI_PFT',ndims,idims)
        call shdf5_irec(ndims,idims,'MMEAN_LAI_PFT',rvara=mmean_lai)
        call shdf5_info('MMEAN_GPP_CO',ndims,idims)
        call shdf5_irec(ndims,idims,'MMEAN_GPP_CO',rvara=mmean_gpp_co)
        call shdf5_irec(ndims,idims,'NPLANT',rvara=nplant)
        call shdf5_irec(ndims,idims,'HITE',rvara=hite)
        call shdf5_irec(ndims,idims,'PFT',ivara=pft)
        call shdf5_info('AREA',ndims,idims)
        call shdf5_irec(ndims,idims,'AREA',rvara=area)
        call shdf5_irec(ndims,idims,'PACO_ID',ivara=paco_id)
        call shdf5_irec(ndims,idims,'PACO_N',ivara=paco_n)
        npatch = idims(1)
        
        do ipatch = 1, npatch
           do icohort = paco_id(ipatch),paco_id(ipatch)+paco_n(ipatch)-1
              if(pft(icohort) == 6)then
                 gpp_pine_out(im,yind) = gpp_pine_out(im,yind) + mmean_gpp_co(icohort) * nplant(icohort) * area(ipatch)/12.
              else
                 gpp_oak_out(im,yind) = gpp_oak_out(im,yind) + mmean_gpp_co(icohort) * nplant(icohort) * area(ipatch)/12.
              endif
           enddo
        enddo
        
        do ipatch = 1, npatch
           do icohort = paco_id(ipatch),paco_id(ipatch)+paco_n(ipatch)-1
              hclass = min(7,1+int(hite(icohort)/5.))
              gpp_hclass_out(im,yind,hclass) = gpp_hclass_out(im,yind,hclass) + mmean_gpp_co(icohort)*nplant(icohort)/12.*area(ipatch)
           enddo
        enddo

        gpp_out(im,yind) = mmean_gpp/12.
        nep_out(im,yind) = mmean_nep/12.
        plresp_out(im,yind) = mmean_plresp/12.
        rh_out(im,yind) = mmean_rh/12.
        fsw_out(im,yind) = mmean_fsw
        transp_out(im,yind) = mmean_transp
        evap_out(im,yind) = mmean_evap
        cantmp_out(im,yind) = mmean_cantmp
        rshort_out(im,yind) = mmean_rshort
        precip_out(im,yind) = mmean_pcpg
        lai_pine_out(im,yind) = mmean_lai(6)
        lai_oak_out(im,yind) = mmean_lai(10)
        lai_out(im,yind) = mmean_lai(6) + mmean_lai(10)
        do k = 11,21
           soilw_out(im,yind) = soilw_out(im,yind) + mmean_soilw(k)*(slz(k)-slz(k+1))/12.
        enddo

        agpp_out(yind) = agpp_out(yind) + mmean_gpp/12.
        anep_out(yind) = anep_out(yind) + mmean_nep/12.
        aplresp_out(yind) = aplresp_out(yind) + mmean_plresp/12.
        arh_out(yind) = arh_out(yind) + mmean_rh/12.
        afsw_out(yind) = afsw_out(yind) + mmean_fsw
        aevap_out(yind) = aevap_out(yind) + mmean_evap
        atransp_out(yind) = atransp_out(yind) + mmean_transp
        acantmp_out(yind) = acantmp_out(yind) + mmean_cantmp
        arshort_out(yind) = arshort_out(yind) + mmean_rshort
        aprecip_out(yind) = aprecip_out(yind) + mmean_pcpg
        do k = 11,21
           asoilw_out(yind) = asoilw_out(yind) + mmean_soilw(k)*(slz(k)-slz(k+1))/12.
        enddo

        do ipatch = 1, npatch
           do icohort = paco_id(ipatch),paco_id(ipatch)+paco_n(ipatch)-1
              if(pft(icohort) == 6)then
                 agpp_pine_out(yind) = agpp_pine_out(yind) + mmean_gpp_co(icohort) * nplant(icohort) * area(ipatch)/12.
              else
                 agpp_oak_out(yind) = agpp_oak_out(yind) + mmean_gpp_co(icohort) * nplant(icohort) * area(ipatch)/12.
              endif
           enddo
        enddo

        
        call shdf5_close()
     enddo
  enddo

  call shdf5_open(trim(outname),'W',1)
  ndims = 3
  idims(1) = 12
  idims(2) = nyear
  idims(3) = 7
  call shdf5_orec(ndims,idims,'MMEAN_GPP_HCLASS',rvara=gpp_hclass_out)
  ndims = 2
  idims(1) = 12
  idims(2) = nyear
  call shdf5_orec(ndims,idims,'MMEAN_GPP',rvara=gpp_out)
  call shdf5_orec(ndims,idims,'MMEAN_PINE_GPP',rvara=gpp_pine_out)
  call shdf5_orec(ndims,idims,'MMEAN_OAK_GPP',rvara=gpp_oak_out)
  call shdf5_orec(ndims,idims,'MMEAN_NEP',rvara=nep_out)
  call shdf5_orec(ndims,idims,'MMEAN_PLRESP',rvara=plresp_out)
  call shdf5_orec(ndims,idims,'MMEAN_RH',rvara=rh_out)
  call shdf5_orec(ndims,idims,'MMEAN_FSW',rvara=fsw_out)
  call shdf5_orec(ndims,idims,'MMEAN_TRANSP',rvara=transp_out)
  call shdf5_orec(ndims,idims,'MMEAN_EVAP',rvara=evap_out)
  call shdf5_orec(ndims,idims,'MMEAN_CANTMP',rvara=cantmp_out)
  call shdf5_orec(ndims,idims,'MMEAN_RSHORT',rvara=rshort_out)
  call shdf5_orec(ndims,idims,'MMEAN_PRECIP',rvara=precip_out)
  call shdf5_orec(ndims,idims,'MMEAN_SOILW',rvara=soilw_out)
  call shdf5_orec(ndims,idims,'MMEAN_LAI_PINE',rvara=lai_pine_out)
  call shdf5_orec(ndims,idims,'MMEAN_LAI_OAK',rvara=lai_oak_out)
  call shdf5_orec(ndims,idims,'MMEAN_LAI',rvara=lai_out)
  ndims = 1
  idims(1) = nyear
  call shdf5_orec(ndims,idims,'YMEAN_GPP',rvara=agpp_out)
  call shdf5_orec(ndims,idims,'YMEAN_PINE_GPP',rvara=agpp_pine_out)
  call shdf5_orec(ndims,idims,'YMEAN_OAK_GPP',rvara=agpp_oak_out)
  call shdf5_orec(ndims,idims,'YMEAN_NEP',rvara=anep_out)
  call shdf5_orec(ndims,idims,'YMEAN_PLRESP',rvara=aplresp_out)
  call shdf5_orec(ndims,idims,'YMEAN_RH',rvara=arh_out)
  call shdf5_orec(ndims,idims,'YMEAN_FSW',rvara=afsw_out)
  call shdf5_orec(ndims,idims,'YMEAN_TRANSP',rvara=atransp_out)
  call shdf5_orec(ndims,idims,'YMEAN_EVAP',rvara=aevap_out)
  call shdf5_orec(ndims,idims,'YMEAN_RSHORT',rvara=arshort_out)
  call shdf5_orec(ndims,idims,'YMEAN_PRECIP',rvara=aprecip_out)
  call shdf5_orec(ndims,idims,'YMEAN_CANTMP',rvara=acantmp_out)
  call shdf5_orec(ndims,idims,'YMEAN_SOILW',rvara=asoilw_out)
  call shdf5_close()

end program main
