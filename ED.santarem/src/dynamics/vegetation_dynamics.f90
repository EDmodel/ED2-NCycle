!==========================================================================================!
!==========================================================================================!
!     This subroutine is the main driver for the longer-term vegetation dynamics.  This    !
! has become a file by itself to reduce the number of sub-routines that are doubled        !
! between ED-2.1 stand alone and the coupled model.                                        !
!------------------------------------------------------------------------------------------!
subroutine vegetation_dynamics(new_month,new_year)
   use grid_coms        , only : ngrids
   use ed_misc_coms     , only : current_time           & ! intent(in)
                               , dtlsm                  & ! intent(in)
                               , frqsum                 & ! intent(in)
                               , ied_init_mode          ! ! intent(in)
   use disturb_coms     , only : include_fire           ! ! intent(in)
   use disturbance_utils, only : apply_disturbances     & ! subroutine
                               , site_disturbance_rates ! ! subroutine
   use fuse_fiss_utils  , only : fuse_patches           ! ! subroutine
   use ed_state_vars    , only : edgrid_g               & ! intent(inout)
                               , edtype,polygontype,sitetype,patchtype                 ! ! variable type
   use growth_balive    , only : dbalive_dt             & ! subroutine
                               , dbalive_dt_eq_0        ! ! subroutine
   use consts_coms      , only : day_sec                & ! intent(in)
                               , yr_day                 ! ! intent(in)
   use mem_polygons     , only : maxpatch               ! ! intent(in)

 !  use ed_state_vars, only : edtype          & ! structure                    !JL
 !                          , polygontype     & ! structure                    !JL
 !                          , sitetype        ! ! structure                    !JL

use pft_coms, only: c2n_leaf, c2n_storage, c2n_recruit, c2n_stem
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   logical     , intent(in)   :: new_month
   logical     , intent(in)   :: new_year
 !  type(edtype)     , target   :: cgrid                                         !JL
   !----- Local variables. ----------------------------------------------------------------!
   type(edtype), pointer      :: cgrid
   real                       :: tfact1
   real                       :: tfact2
   integer                    :: doy
   integer                    :: ip
   integer                    :: isite
   integer                    :: ifm
   real                       :: oldmsn!JL!
   real                        :: newmsn !JL!
   !----- External functions. -------------------------------------------------------------!
   integer     , external     :: julday
   !---------------------------------------------------------------------------------------!

real :: oldpn, oldsn, newpn, newsn, oldtn, newtn,new_balive,new_bdead,new_bstorage,old_balive,old_bdead,old_bstorage
real :: DailyForestFixation, DailyFixationDemand !JL
 

type(polygontype), pointer :: cpoly
type(sitetype), pointer :: csite
type(patchtype),pointer :: cpatch
integer :: ipy, isi, ipa, ico

oldpn=0.;oldtn=0.;oldsn=0.;newpn=0.;newsn=0.;newtn=0.;newmsn=0.;oldmsn=0. !old and newmsn were added JL
old_balive = .0;old_bdead = .0; old_bstorage = .0;  new_balive = .0; new_bdead = .0; new_bstorage = 0. 



   !----- Find the day of year. -----------------------------------------------------------!
   doy = julday(current_time%month, current_time%date, current_time%year)
   !----- Time factor for normalizing daily variables updated on the DTLSM step. ----------!
   tfact1 = dtlsm / day_sec
   !----- Time factor for averaging dailies. ----------------------------------------------!
   tfact2 = 1.0 / yr_day

   !----- Apply events. -------------------------------------------------------------------!

   call prescribed_event(current_time%year,doy)
 
   !---------------------------------------------------------------------------------------!
   !   Loop over all domains.                                                              !
   !---------------------------------------------------------------------------------------!
   do ifm=1,ngrids

      cgrid => edgrid_g(ifm) 

      do ipy=1,cgrid%npolygons
         cpoly => cgrid%polygon(ipy)
         do isi=1,cpoly%nsites
            csite => cpoly%site(isi)
            do ipa=1,csite%npatches
               oldsn = oldsn + (csite%slow_soil_C(ipa)/10.+csite%fast_soil_N(ipa) + csite%structural_soil_C(ipa)/150. + csite%mineralized_soil_N(ipa)) * csite%area(ipa)
               cpatch => csite%patch(ipa)
               oldpn = oldpn + csite%area(ipa) * (csite%repro(1,ipa)/c2n_recruit(1)+csite%repro(2,ipa)/c2n_recruit(2)+csite%repro(3,ipa)/c2n_recruit(3)+csite%repro(4,ipa)/c2n_recruit(4))
              oldmsn = oldmsn + csite%mineralized_soil_N(ipa) * csite%area(ipa)!JL!
             do ico=1,cpatch%ncohorts
      !            oldpn = oldpn + csite%area(ipa) * cpatch%nplant(ico) * (cpatch%balive(ico)/c2n_leaf(cpatch%pft(ico)) + cpatch%bdead(ico)/c2n_stem(cpatch%pft(ico))+cpatch%bstorage(ico)/c2n_storage)
                    oldpn = oldpn + csite%area(ipa) * cpatch%nplant(ico) * (cpatch%balive(ico)/c2n_leaf(cpatch%pft(ico)) + cpatch%bdead(ico)/c2n_stem(cpatch%pft(ico))+cpatch%nstorage(ico)) !jl changed bstorage to nstorage
                    old_balive = old_balive + csite%area(ipa) * cpatch%nplant(ico) * cpatch%balive(ico)/c2n_leaf(cpatch%pft(ico))
                    old_bdead = old_bdead + csite%area(ipa) * cpatch%nplant(ico) * cpatch%bdead(ico)/c2n_stem(cpatch%pft(ico))
                    old_bstorage = old_bstorage + csite%area(ipa) * cpatch%nplant(ico) * cpatch%nstorage(ico)

              enddo
            enddo
         enddo
      enddo
print*,'DN, INITIAL', 'SoilN',oldsn,'PlantN',oldpn, 'ForestN',oldsn+oldpn,'MSN',oldmsn,'BaliveN',old_balive,'BdeadN',old_bdead,'NstorageN',old_bstorage


      !     The following block corresponds to the daily time-step.                        !
      !------------------------------------------------------------------------------------!
      !----- Standardise the fast-scale uptake and respiration, for growth rates. ---------!
 
     call normalize_ed_daily_vars(cgrid, tfact1)

      !----- Update phenology and growth of live tissues. ---------------------------------!
      select case (ied_init_mode)
      case (-8)
         !----- Special case, in which we don't solve the actual vegetation dynamics. -----!
         call phenology_driver_eq_0(cgrid,doy,current_time%month, tfact1)
         call dbalive_dt_eq_0(cgrid,tfact2)
      case default

         call phenology_driver(cgrid,doy,current_time%month, tfact1)
  
         call dbalive_dt(cgrid,tfact2)

      end select
      !------------------------------------------------------------------------------------!

      !------------------------------------------------------------------------------------!
      !     The following block corresponds to the monthly time-step:                      !
      !------------------------------------------------------------------------------------!
      if (new_month) then

         !----- Update the mean workload counter. -----------------------------------------!
         call update_workload(cgrid)
 !print*,'1'
         !----- Update the growth of the structural biomass. ------------------------------!
         call structural_growth(cgrid, current_time%month)


         !----- Solve the reproduction rates. ---------------------------------------------!

 !print*,'2'
            call reproduction(cgrid,current_time%month)

! print*,'3'

         !----- Update the fire disturbance rates. ----------------------------------------!
         if (include_fire /= 0) then
            call fire_frequency(current_time%month,cgrid) 
        end if

! print*,'4'
         !----- Update the disturbance rates. ---------------------------------------------!
         call site_disturbance_rates(current_time%month, current_time%year, cgrid)
      endif
! print*,'5'
      !------  update dmean and mmean values for NPP allocation terms ---------------------!
      call normalize_ed_dailyNPP_vars(cgrid)
 ! print*,'6'    
      !------------------------------------------------------------------------------------!
      !     This should be done every day, but after the longer-scale steps.  We update    !
      ! the carbon and nitrogen pools, and re-set the daily variables.                     !
      !------------------------------------------------------------------------------------!
      call update_C_and_N_pools(cgrid)
 !print*,'7' 
     call zero_ed_daily_vars(cgrid)
 !print*,'8'
      !------------------------------------------------------------------------------------!
 
! print*,'DailyForestFixation', DailyForestFixation, 'DailyFixationDemand', DailyFixationDemand
    !-----Zero Yearly Nitrogen Variables. --------------------------------------------------!
    ! Yearly variables are reset every June. Reset the Annual Nitrogen variables to 0 on 
    ! the first day of June! JL
     if (current_time%month == 6 .and.current_time%date == 2) then 
        cgrid%total_N_Fixation(1) = 0.0
        cgrid%total_DON_loss(1)   = 0.0
        cgrid%total_DIN_loss(1)   = 0.0
        cgrid%total_Ngas_loss(1)  = 0.0
        print*, 'YEARLY ZERO'
     endif



     do ipy=1,cgrid%npolygons
         cpoly => cgrid%polygon(ipy)
         do isi=1,cpoly%nsites
            csite => cpoly%site(isi)
            do ipa=1,csite%npatches
               newsn = newsn + (csite%slow_soil_C(ipa)/10.+csite%fast_soil_N(ipa) + csite%structural_soil_C(ipa)/150. + csite%mineralized_soil_N(ipa)) * csite%area(ipa)
               cpatch => csite%patch(ipa)
               newpn = newpn + csite%area(ipa) * (csite%repro(1,ipa)/c2n_recruit(1)+csite%repro(2,ipa)/c2n_recruit(2)+csite%repro(3,ipa)/c2n_recruit(3)+csite%repro(4,ipa)/c2n_recruit(4))
                  newmsn = newmsn + csite%mineralized_soil_N(ipa) * csite%area(ipa)!JL! 
             do ico=1,cpatch%ncohorts
                    newpn = newpn + csite%area(ipa) * cpatch%nplant(ico) * (cpatch%balive(ico)/c2n_leaf(cpatch%pft(ico)) &
                            + cpatch%bdead(ico)/c2n_stem(cpatch%pft(ico))+cpatch%nstorage(ico))
                  new_balive = new_balive + csite%area(ipa) * cpatch%nplant(ico) * cpatch%balive(ico)/c2n_leaf(cpatch%pft(ico))
                  new_bdead = new_bdead + csite%area(ipa) * cpatch%nplant(ico) * cpatch%bdead(ico)/c2n_stem(cpatch%pft(ico))
                  new_bstorage = new_bstorage + csite%area(ipa) * cpatch%nplant(ico) * cpatch%nstorage(ico)
               enddo
            enddo
         enddo
      enddo
print*,'DN, FINAL','SoilN',newsn,'PlantN',newpn,'ForestN',newsn+newpn,'MSN', newmsn,'Balive',new_balive,'Bdead',new_bdead,'Nstorage',new_bstorage ! 


cgrid%ForestN(1)           = newsn+newpn !JL!
cgrid%PlantN(1)            = newpn !JL!
cgrid%SoilN(1)             = newsn !JL!

         !----- This is actually the yearly time-step, apply the disturbances. ------------!
         if (new_month .and. new_year) then
 
           call apply_disturbances(cgrid)
 !print*,'9'
         end if


      !------------------------------------------------------------------------------------!
      !      Fuse patches last, after all updates have been applied.  This reduces the     !
      ! number of patch variables that actually need to be fused.                          !

      !------------------------------------------------------------------------------------!

               

!print*,'veg fuse patches'
      if(new_year) then
         if (maxpatch >= 0) call fuse_patches(cgrid,ifm)
      end if
 !print*,'10'
      !------------------------------------------------------------------------------------!

      !----- Recalculate the AGB and basal area at the polygon level. ---------------------!
      call update_polygon_derived_props(cgrid)

 !print*,'11'
      call print_C_and_N_budgets(cgrid)
 !print*,'12'
      !------------------------------------------------------------------------------------!

   end do

   return
end subroutine vegetation_dynamics
!==========================================================================================!
!==========================================================================================!
