!==========================================================================================!
!==========================================================================================!
module growth_balive
   !=======================================================================================!
   !=======================================================================================!


   contains



   !=======================================================================================!
   !=======================================================================================!
   !     This subroutine will update the alive biomass, and compute the respiration terms  !
   ! other than leaf respiration.                                                          !
   ! IMPORTANT: The order of the operations here affect the C/N budgets, so don't change   !
   !            the order of the operations unless you really know what you are doing.     !
   !---------------------------------------------------------------------------------------!
   subroutine dbalive_dt(cgrid, tfact)
      use ed_state_vars   , only : edtype                 & ! structure
                                 , polygontype            & ! structure
                                 , sitetype               & ! structure
                                 , patchtype              ! ! structure
      use pft_coms        , only : q                      & ! intent(in)
                                 , qsw                    & ! intent(in)
                                 , plant_N_supply_scale   & ! intent(in)
                                 , c2n_storage            & ! intent(in)
                                 , c2n_leaf               &  ! intent(in)             !JL!   
                                 , c2n_stem               & ! intent(in)              !JL!                
                                 , growth_resp_factor     & ! intent(in)
                                 , storage_turnover_rate  & ! intent(in)
                                 , phenology              ! ! intent(in)
      use physiology_coms , only : N_plant_lim            ! ! intent(in)
      use grid_coms       , only : nzg                    ! ! intent(in)
      use ed_therm_lib    , only : calc_veg_hcap          & ! function
                                 , update_veg_energy_cweh ! ! function
      use allometry       , only : area_indices           & ! subroutine
                                 , ed_biomass             &! ! function
                                 , dbh2bl                  ! !JL
      use mortality       , only : mortality_rates        ! ! subroutine
      use phenology_coms  , only : theta_crit             ! ! intent(in)
      use decomp_coms     , only : FSN_ndays_to_avg                                   !JL!
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type(edtype)     , target     :: cgrid
      real             , intent(in) :: tfact
      !----- Local variables. -------------------------------------------------------------!
      type(polygontype), pointer    :: cpoly
      type(sitetype)   , pointer    :: csite
      type(patchtype)  , pointer    :: cpatch
      integer                       :: ipy
      integer                       :: isi
      integer                       :: ipa
      integer                       :: ico
      integer                       :: ipft
      real                          :: salloc
      real                          :: salloci
      real                          :: bl
      real                          :: br
      real                          :: daily_C_gain
      real                          :: carbon_balance
      real                          :: carbon_balance_pot
      real                          :: carbon_balance_max
      real                          :: balive_in
      real                          :: dndt
      real                          :: old_leaf_hcap
      real                          :: old_wood_hcap
      real                          :: nitrogen_uptake
      real                          :: N_uptake_pot
      real                          :: temp_dep
      real                          :: shadow_N_uptake2 !JL!
      real                          :: Cost_fixer !JL!
      real                          :: Cost_BNF   !JL!
      real                          :: C_cost      !JL!
      real                          :: ObligateCost      !JL!
      real                          :: Nstorage_after_fixation !JL!
      real                          :: N_to_storage !JL!
      real                          :: N_fixed_Extra !JL!
      real                          :: Nstorage_before_fixation !JL!     
      real                          :: N  
      real                          :: i            !JL! 
      real                          :: nstorage_max !JL!
      real                          :: Nstorage_before_resorption !JL!
      real                          :: Nstorage_after_resorption !JL
      real                          :: N_not_resorbed
      real                          :: N_resorption
      real                          :: Density !JL!
      real                          :: total_N_supply !JL!
      real                          ::  DailyForestFixation !JL!
      !------------------------------------------------------------------------------------!


      do ipy = 1,cgrid%npolygons
         cpoly => cgrid%polygon(ipy)

         do isi = 1,cpoly%nsites
            csite => cpoly%site(isi)

            do ipa = 1,csite%npatches
               cpatch => csite%patch(ipa)

               !----- Reset averaged variables. -------------------------------------------!
               csite%total_plant_nitrogen_uptake(ipa) = 0.0
               !----- Loop over cohorts. --------------------------------------------------!
               do ico = 1,cpatch%ncohorts

               !-----JL checks. -----------------------------------------------------------!
               ! print*,'ico',ico,'INIT NPLANT',sum(cpatch%nplant * (cpatch%balive         &
               ! / c2n_leaf(cpatch%pft(ico)) + cpatch%bdead / c2n_stem(cpatch%pft(ico)) +  &
               ! cpatch%nstorage)) * csite%area(ipa)
               !print*,'ico',ico,'INIT BALIVE',sum(cpatch%nplant * (cpatch%balive  /       &
               !c2n_leaf(cpatch%pft(ico)))) * csite%area(ipa)
               !print*,'ico',ico,'INIT NSTORAGE',sum(cpatch%nplant * cpatch%nstorage)      &
               !* csite%area(ipa)
               !print*,'ico',ico,'INIT N UPTAKE',csite%total_plant_nitrogen_uptake(ipa)
               !print*,'ico',ico,'INIT BSTORAGE',sum(cpatch%nplant * cpatch%nstorage       &
               !/ c2n_leaf(cpatch%pft(ico))) * csite%area(ipa)
               !-----JL checks. -----------------------------------------------------------!

                !----- Alias for current PFT. -------------------------------------------!
                  ipft = cpatch%pft(ico)

                  !----- Update the elongation factor. ------------------------------------!
                  select case (phenology(ipft))
                  case (4)
                     cpatch%elongf(ico) = max(0.0, min(1.0, cpatch%paw_avg(ico)            &
                                                          / theta_crit))
                  case default
                     cpatch%elongf(ico) = 1.0

                  end select


                  !----- Initialize cohort nitrogen uptake. -------------------------------!
                  nitrogen_uptake = 0.0
                  N_uptake_pot    = 0.0
                  shadow_N_uptake2= 0.0 !JL!
                  cpatch%shadow_N_uptake(ico) = 0.0 !JL!

                  !----- Set allocation factors. ------------------------------------------!
                  salloc  = 1.0 + qsw(ipft) * cpatch%hite(ico) + q(ipft)
                  salloci = 1.0 / salloc

                  !------------------------------------------------------------------------!
                  !     Compute maintenance costs using actual pools.                      !
                  !------------------------------------------------------------------------!
                  call plant_maintenance(cpatch,ico,cpatch%broot(ico),cpatch%bleaf(ico)    &
                                        ,tfact,daily_C_gain,csite%avg_daily_temp(ipa))




                 !----- Subtract maintenance costs from pools. ---------------------------!
                  cpatch%balive(ico)    = cpatch%balive(ico)                               &
                                        - cpatch%leaf_maintenance(ico)                     &
                                        - cpatch%root_maintenance(ico)
                  cpatch%bleaf(ico)     = cpatch%bleaf(ico)                                &
                                        - cpatch%leaf_maintenance(ico)       
                  cpatch%broot(ico)     = cpatch%broot(ico)                                &
                                        - cpatch%root_maintenance(ico)
                  cpatch%cb(13,ico)     = cpatch%cb(13,ico)                                &
                                        - cpatch%leaf_maintenance(ico)                     &
                                        - cpatch%root_maintenance(ico)
                  cpatch%cb_max(13,ico) = cpatch%cb_max(13,ico)                            &
                                        - cpatch%leaf_maintenance(ico)                     &
                                        - cpatch%root_maintenance(ico)

                  !------------------------------------------------------------------------!
                  !    Storage respriation/turnover_rate.                                  !
                  !    Calculate in same way as leaf and root turnover in kgC/plant/year.  !
                  !------------------------------------------------------------------------!

                  !------------------------------------------------------------------------!
                  !     The commented line is an experimental and arbitrary test, borrowed !
                  ! from maintainence temperature dependency. [[MCD]]                      !
                  !------------------------------------------------------------------------!
                  ! temp_dep = 1.0                                                         &
                  !          / ( 1.0  + exp( 0.4 * (278.15 - csite%avg_daily_temp(ipa))))
                  temp_dep = 1.0
                  !------------------------------------------------------------------------!

                  cpatch%storage_respiration(ico) = cpatch%bstorage(ico)                   &
                                                  * storage_turnover_rate(ipft)            &
                                                  * tfact * temp_dep

                  cpatch%bstorage(ico) = cpatch%bstorage(ico)                              &
                                         - cpatch%storage_respiration(ico)

                  !------------------------------------------------------------------------!
                  !     When storage carbon is lost, allow the associated nitrogen to go   !
                  ! to litter in order to maintain prescribed C2N ratio.                   !
                  !------------------------------------------------------------------------!
                  !csite%fsn_in(ipa) = csite%fsn_in(ipa) + cpatch%storage_respiration(ico)  &
                  !                                     / c2n_storage * cpatch%nplant(ico)
           
                  !Storage respiration cost for nitrogen is not active so
                  ! fsn does not need to be updated here    

                  !------------------------------------------------------------------------!
                  !      Calculate actual, potential and maximum carbon balances.          !
                  !------------------------------------------------------------------------!
                  call plant_carbon_balances(cpatch,ipa,ico,daily_C_gain,carbon_balance    &
                                            ,carbon_balance_pot,carbon_balance_max)
                     cpatch%carbon_balance(ico) = carbon_balance !JL!

                  !------------------------------------------------------------------------!

                  !------------------------------------------------------------------------! 
                  !      Compute respiration rates for coming day [kgC/plant/day].         !
                  !------------------------------------------------------------------------!
                  cpatch%growth_respiration(ico) = max(0.0, daily_C_gain                   &
                                                          * growth_resp_factor(ipft))
                  !------------------------------------------------------------------------!



                  !------------------------------------------------------------------------!
                  !     Find the "virtual" leaf respiration.                               !
                  !------------------------------------------------------------------------!
                  cpatch%vleaf_respiration(ico) = (1.0-cpoly%green_leaf_factor(ipft,isi))  &
                                                * salloci * cpatch%balive(ico)             &
                                                * storage_turnover_rate(ipft)              &
                                                * tfact * temp_dep
                  !------------------------------------------------------------------------!


                  !------------------------------------------------------------------------!
                  !Define maximum Nstorage and update N storage pool now that maintence    !
                  ! costs are known. Resorb nitrogen from leaf turnover. Fast N input      !
                  !through mortality or disturbance does not get resorbed.                 !
                  !The resorption rate is 0.48% (McGroddy et al. 2004).                    !
                  !Once the Nstorage max is reached put the rest of the N into FSN_in      !         
                  !------------------------------------------------------------------------! 

                  !define maximum Nstorage 

                  Nstorage_before_resorption = cpatch%nstorage(ico)

                  nstorage_max = dbh2bl(cpatch%dbh(ico),ipft)/c2n_leaf(ipft) *1.33

                  if (cpatch%nstorage(ico) <  nstorage_max) then                                
                     cpatch%nstorage(ico) = min((cpatch%nstorage(ico) +                     &
                                            ((( cpatch%leaf_maintenance(ico))               &
                                          / c2n_leaf(ipft))*0.48)),nstorage_max)
                  endif
                  !------------------------------------------------------------------------!

                  !------------------------------------------------------------------------! 
                  ! Add unresorbed N that would have exceeded nstorage_max to the fast     !  
                  ! liter pool.The fraction of N that is not resorbed will be added to     !
                  ! fsn_in in the litter subroutine                                        !     
                  !------------------------------------------------------------------------!

                  N_resorption = 0
                  N_not_resorbed = 0

                  Nstorage_after_resorption = cpatch%nstorage(ico)
                  N_resorption = Nstorage_after_resorption - Nstorage_before_resorption
                  N_not_resorbed = (((cpatch%leaf_maintenance(ico)                         &
                                 + cpatch%root_maintenance(ico))/ c2n_leaf(ipft))*0.48)    &
                                 - N_resorption  

                  csite%fsn_in(ipa) = csite%fsn_in(ipa)                                    &
                                    + N_not_resorbed * cpatch%nplant(ico)
                  !------------------------------------------------------------------------!


                  !-----------------------------------------------------------------------!
                  !    Shadow calculation for Nitrogen Uptake                              ! !JL,XXT
                  !------------------------------------------------------------------------!
                  balive_in = cpatch%balive(ico)                                   
                  call shadow_n(csite,ipa,ico,salloc,salloci,carbon_balance   &
                                            ,shadow_N_uptake2                               &
                                            ,cpoly%green_leaf_factor(ipft,isi))
                   cpatch%shadow_N_uptake(ico) = shadow_N_uptake2

                  !------------------------------------------------------------------------!
                 


                  !------------------------------------------------------------------------!
                  !     Do a shadow calculation to see what would have happened if stomata !
                  ! were open.  This is used to calculate potential nitrogen uptake,       !
                  ! N_uptake_pot.                                                          !
                  !------------------------------------------------------------------------!
                  if (N_plant_lim == 1) then
                     call potential_N_uptake(cpatch,ico,salloc,salloci,balive_in           &
                                            ,carbon_balance_pot,N_uptake_pot               &
                                            ,cpoly%green_leaf_factor(ipft,isi))
                 
                  cpatch%carbon_balance_pot(ico) = carbon_balance_pot !JL!
                  cpatch%N_uptake_pot(ico) = N_uptake_pot !JL!

                  end if

                  !------------------------------------------------------------------------!

                  !------------------------------------------------------------------------!
                  !                      Checks !JL!                                       !
                  !------------------------------------------------------------------------! 
                  ! print*, 'CB pot', cpatch%carbon_balance_pot(ico),                      &
                  !                 'CB',cpatch%carbon_balance(ico)
                  !------------------------------------------------------------------------!  

                  !------------------------------------------------------------------------!  
                  !      Calculate plant N limitation factor, N fixation,                  !
                  !      actual nitrogen uptake constrained by soil supply                 !
                  !      and excess carbon assimilated                                     !
                  !------------------------------------------------------------------------! 
                
                   !by XXT,JL Initialize Nitrogen parameters to 0 
 
                  !------------------------------------------------------------------------!  
                  cpatch%excess_carbon(ico) = 0.0 
                  cpatch%excess_carbon_fixer(ico) = 0.0                   
                  cpatch%actual_nitrogen_uptake(ico) = 0.0
                  cpatch%actual_nitrogen_uptake_fixer(ico) = 0.0
                  cpatch%fixation_demand(ico) = 0.0
                  cpatch%N_fixation(ico) = 0.0 
                  cpatch%nitrogen_supply(ico) = 0.0 
                  cpatch%nitrogen_supply_fixer(ico) = 0.0 
                  cpatch%N_fixation(ico) = 0.0 
                  Density = 0.0 
                  !------------------------------------------------------------------------!      

                  !------------------------------------------------------------------------!  
                  !If plants are not Nitrogen limited set fsn = 1                          ! 
                  !If there's no limitation then excess_carbon = 0                         !
                  !------------------------------------------------------------------------!               
              !    if (n_plant_lim == 0 .or. N_uptake_pot <= 0.0 .and. ipft .ne. 31 )  then ! IF #1  
                  if (n_plant_lim == 0 .or.  N_uptake_pot <= 0.0) then 
                     if(ipft .ne. 31) then
                     cpatch%fsn(ico) = 1.0
                     cpatch%fsn_fixer(ico) = 1.0   !JL    

		     cpatch%nitrogen_supply(ico) = (cpatch%broot(ico)                      &
                                                 / sum(cpatch%broot * cpatch%nplant))      &
                                                 * csite%mineralized_soil_N(ipa)

		     cpatch%actual_nitrogen_uptake(ico) = min(cpatch%nitrogen_supply(ico)  &
                                                        ,cpatch%shadow_N_uptake(ico))  
                     endif

                  !------------------------------------------------------------------------! 

                  !------------------------------------------------------------------------! 
                  ! If plants are non-fixers calculate nitrogen supply, actual nitrogen    ! 
                  ! uptake, excess carbon, and nitrogen limitation factor.                 !
                  ! Excess carbon in put into fsn_in in the litter subroutine              !
                  !------------------------------------------------------------------------! 
                  else if (ipft < 30 )then   
                  print*, 'ipft',ipft,'NON FIXER START'

		    cpatch%nitrogen_supply(ico) = (cpatch%broot(ico)                        &
                                                / sum(cpatch%broot * cpatch%nplant))        &
				         	* csite%mineralized_soil_N(ipa)			

		    cpatch%actual_nitrogen_uptake(ico) = min( cpatch%nitrogen_supply(ico)   &
                                                       , cpatch%shadow_N_uptake(ico))           		  
 
                    if ( cpatch%shadow_N_uptake(ico) >  cpatch%actual_nitrogen_uptake(ico)  &
                    .and. cpatch%phenology_status(ico) == 1 ) then
 
                           cpatch%excess_carbon(ico) = ( cpatch%shadow_N_uptake(ico)        &
                                                     - cpatch%actual_nitrogen_uptake(ico))  &
                                                     * c2n_leaf(ipft) 
                    endif

                   cpatch%fsn(ico) =  cpatch%nitrogen_supply(ico)                           &
                                   / (cpatch%nitrogen_supply(ico)                           &
                                   +  cpatch%N_uptake_pot(ico))  
                                  
                  !------------------------------------------------------------------------!                  
            !       print*, 'shadowN',cpatch%shadow_N_uptake(ico),                       &
            !               'NSUPPLY', cpatch%nitrogen_supply(ico),                      &
            !               'actual',cpatch%actual_nitrogen_uptake(ico),                 &
            !               'excess C',cpatch%excess_carbon(ico),                        & 
            !              ' N_uptake_pot', N_uptake_pot ,'CB',carbon_balance,          &
            !             'phenology',  cpatch%phenology_status(ico)
              !             'cpatch%broot(ico)',cpatch%broot(ico),                       &
              !             'all roots', sum(cpatch%broot * cpatch%nplant),              &
              !              'msn', csite%mineralized_soil_N(ipa),                       &
              !             'ipa',ipa
                   print*, 'NON FIXER END', 'fsn', cpatch%fsn(ico) 
                  !------------------------------------------------------------------------! 

                  !------------------------------------------------------------------------! 
                  ! If plants are facultative-fixers calculate nitrogen supply, actual     !
                  ! nitrogen uptake, biological nitrogen fixation, excess carbon, and      !
                  ! nitrogen limitation factor                                             !
                  !------------------------------------------------------------------------! 
                  else if (ipft == 30)then !NITROGEN FIXER Facultative 
                  print*, 'FACULTATIVE FIXER START', 'ipft', ipft
 
                  ! Define carbon costs for BNF

                  Cost_fixer = 0 ! #arbitrary number ( g C gN -1)
                  Cost_BNF = 9.12 !#( 9.12 g C gN -1)  (gutschick, 1981)     
                  C_cost = Cost_fixer + Cost_BNF ! #kg C kgN -1
                  !------------------------------------------------------------------------!
                  !print*, 'C_cost', C_cost
                  !------------------------------------------------------------------------!
           
                  ! Nitrogen fixer supply comes from soil and fixation.
                  ! Fixation rate is based off of Sarah Batterman's field data 
                  ! 0.0003941912 kg N fixed/total tree biomass(above and belowground)/day
                  ! Total tree biomass(kg biomass) =(bstorage + bdead+balive)*2

                  cpatch%nitrogen_supply(ico) =  (cpatch%broot(ico)                        &
                                              / sum(cpatch%broot * cpatch%nplant))         &
                                              * csite%mineralized_soil_N(ipa)            


		     !soil N uptake
                  cpatch%actual_nitrogen_uptake(ico) = min( cpatch%nitrogen_supply(ico)    &
                                                     , cpatch%shadow_N_uptake(ico))     ! XXT		
                  !------------------------------------------------------------------------!                 
                  !print*, ' cpatch%shadow_N_uptake(ico)', cpatch%shadow_N_uptake(ico), &
                  !        ' Soil N supply', cpatch%nitrogen_supply(ico), &
                  !        ' Actual N uptake from soil', cpatch%actual_nitrogen_uptake(ico)
                  !------------------------------------------------------------------------!


                     if ( cpatch%shadow_N_uptake(ico) > cpatch%actual_nitrogen_uptake(ico) &
                     .and. cpatch%phenology_status(ico) == 1 ) then
                 
                    !excess C before fixation, only considering soil N        
                     cpatch%excess_carbon(ico) = (cpatch%shadow_N_uptake(ico)              &
                                               -  cpatch%actual_nitrogen_uptake(ico))      &
                                               * c2n_leaf(ipft)  
                  !------------------------------------------------------------------------!
                  !   print*,' excess C before fixation',cpatch%excess_carbon(ico), &
                  !        'N needed', (cpatch%shadow_N_uptake(ico) -  cpatch%actual_nitrogen_uptake(ico)),&
                  !        'c2n_leaf(ipft)',c2n_leaf(ipft)
                  !------------------------------------------------------------------------!
                     endif
       
                  ! Calculate Fixation Demand
                  ! If the plant is in postive carbon balance and if actual_nitrogen_uptake 
                  ! (from soil) does not satisfy N demand then Fixation demand is the amount
                  ! of nitrogen the plant would fix (kg N) based on the daily carbon_balance
                  ! taking into account the carbon cost of fixation.
                  ! This uses the potential gain in C from doing 
                  ! additonal fixation (numerator) divided by the cost and c2n raio of the 
                  !plant structure that the nitrogen would be going to (depends on phenology)
     
                     if (cpatch%carbon_balance(ico) > 0                                     &
                     .and. cpatch%shadow_N_uptake(ico) > cpatch%actual_nitrogen_uptake(ico) &
                     .and. cpatch%phenology_status(ico) == 1) then                                                                                                      
                    cpatch%fixation_demand(ico) = cpatch%excess_carbon(ico)                 &
                                                 /(C_cost + c2n_leaf(ipft))                                                            


                 ! Calculate Biological Nitrogen Fixation.This comes at a carbon cost      
                 ! which is incorproated into the fixation_demand term. 
                 ! note, BNF maximum rate observed is 0.0003941912. The mean rate is
                 ! 1.784159e-05  at Sarah's field site(Agua Salud) 
               
                     cpatch%N_fixation(ico) =  min(((cpatch%BSTORAGE(ico) +cpatch%BDEAD(ico)&
                                            + cpatch%BALIVE(ico))*2*0.0003941912)           &
                                            , cpatch%fixation_demand(ico))
                     else                               
                           cpatch%N_fixation(ico) = 0
                     endif

                  !------------------------------------------------------------------------!                  
                  !   print*,' cpatch%fixation_demand(ico)', cpatch%fixation_demand(ico),  &
                  !        '(C_cost + c2n_leaf(ipft))',(C_cost + c2n_leaf(ipft)), &
                  !        ' cpatch%N_fixation(ico)',cpatch%N_fixation(ico),  &
                  !        'max daily fixation', ((cpatch%BSTORAGE(ico) +cpatch%BDEAD(ico)+&
                  !          cpatch%BALIVE(ico))*2*0.0003941912) 
                  !------------------------------------------------------------------------!
  
                ! include N acquired from soil and BNF
                  cpatch%actual_nitrogen_uptake(ico) = cpatch%actual_nitrogen_uptake(ico)   &
                                                     +  cpatch%N_fixation(ico) 

                  !------------------------------------------------------------------------!
                  !print*, 'actual N uptake after fixation',  cpatch%actual_nitrogen_uptake(ico), &
                  !     ' N fixation', cpatch%N_fixation(ico)
                  !------------------------------------------------------------------------!                                                                                   
                  ! Recalculate excess C based on additional N acquired by fixation
                  ! This is the carbon cost of fixation + any carbon that the plant cannot
                  ! use becaue total_N_supply < shadow_N_uptake
                     if (cpatch%shadow_N_uptake(ico) >  cpatch%actual_nitrogen_uptake(ico)  &
                     .and. cpatch%phenology_status(ico) == 1 ) then
                         
                      cpatch%excess_carbon(ico) = (cpatch%shadow_N_uptake(ico)              & 
                                                - cpatch%actual_nitrogen_uptake(ico)) 	    &
				            	* c2n_leaf(ipft)  
                     endif
                  !------------------------------------------------------------------------!
                  !   print*, ' excess C after fixation', cpatch%excess_carbon(ico)
                  !------------------------------------------------------------------------!

                  ! Calculate N limitation factor 
                  total_N_supply =  (cpatch%broot(ico)                                     &
                                              / sum(cpatch%broot * cpatch%nplant))         &
                                              * csite%mineralized_soil_N(ipa)              &
                                              + ((cpatch%BSTORAGE(ico)+cpatch%BDEAD(ico)   &
                                              + cpatch%BALIVE(ico))*2*0.0003941912)
                  

                  cpatch%fsn(ico) =  total_N_supply                                        &
                                  / (total_N_supply                                        &
                                  + cpatch% N_uptake_pot(ico))

                  !Reset cpatch%nitrogen_supply(ico) to total N supply because this value
                  ! is used below to calculate how much N it can put into nstorage
                  ! as well as in the structural growth file to calculate how much nitrogen
                  ! a plant has (in addition to its storage) to put towards structural 
                  !growth and reproduction.               

                  cpatch%nitrogen_supply(ico) =   total_N_supply
                  !------------------------------------------------------------------------!
                  ! print*,'FACULTATIVE FIXER END'!, 'fsn', cpatch%fsn(ico) ,'NSUPPLY',     &
                  ! cpatch%nitrogen_supply(ico),' nitrogen_uptake',nitrogen_uptake,         &
                  !'cpatch%shadow_N_uptake(ico)',cpatch%shadow_N_uptake(ico)  ,'FIXATION',  & 
                  !cpatch%N_fixation(ico),'actual N upatake'                                &
                  !,cpatch%actual_nitrogen_uptake(ico),' N_uptake_pot',                     &
                  !cpatch%N_uptake_pot(ico), 'soil supply',                                 &
                  !((cpatch%broot(ico)/ sum(cpatch%broot * cpatch%nplant))                  &
                  ! * csite%mineralized_soil_N(ipa))   

                  !------------------------------------------------------------------------! 
                  ! If plants are obligate-fixers calculate nitrogen supply, actual        !
                  ! nitrogen uptake, biological nitrogen fixation, excess carbon, and      !
                  ! nitrogen limitation factor                                             !
                  !------------------------------------------------------------------------! 

                  else if (ipft == 31)then !NITROGEN FIXER Obligate 

                  !------------------------------------------------------------------------!                                        

                  ! Define carbon costs for BNF


                  Cost_fixer = 0 ! #arbitrary number ( g C gN -1), this should be evaluated   
                  Cost_BNF = 9.12 !#( 9.12 g C gN -1), metabolic cost, (gutschick, 1981)     
                  C_cost = Cost_fixer + Cost_BNF ! #kg C kgN -1



                  !Calculate fixation if plant fixed all day
                  cpatch%N_fixation(ico) =  ((cpatch%BSTORAGE(ico) +cpatch%BDEAD(ico)       &
                                            + cpatch%BALIVE(ico))*2*0.0003941912) 
          


                 ! The fixation comes at a carbon cost. Calculate the carbon cost associated 
                 ! with the daily fixation. First set the cost to 0 for the day and then
                 ! calculate the Carbon used for fixation
  
                  ObligateCost = 0
                  ObligateCost =  cpatch%N_fixation(ico)* C_cost  ! [ kg C/plant]



                 ! Remove the Carbon used for fixation from the plant's carbon balance by
                 ! putting it into the excess C term. This term will be removed from the
                 ! carbon balance before the plant can use the carbon for growth 
                   
                  cpatch%excess_carbon(ico)  = ObligateCost


                  !------------------------------------------------------------------------!  
                 !Determine Nitrogen supply for plant. This will be used for the fsn 
                 !calculation.Nitrogen supply will be the plants access to N from soil
                 !and fixation.
                 
                  if(cpatch%broot(ico)<0)then
                   cpatch%broot(ico) = 0
                  endif              

                  total_N_supply              =  (cpatch%broot(ico)                        &
                                              / sum(cpatch%broot * cpatch%nplant))         &
                                              * csite%mineralized_soil_N(ipa)              &
                                              + ((cpatch%BSTORAGE(ico)+cpatch%BDEAD(ico)   &
                                              + cpatch%BALIVE(ico))*2*0.0003941912)


!print*, ' total_N_supply', total_N_supply,&
!'soil N supply',(cpatch%broot(ico)  / sum(cpatch%broot * cpatch%nplant)) * csite%mineralized_soil_N(ipa)
                  !------------------------------------------------------------------------! 

                  !Determine How much N the plant takes in (soil + fixation) and update
                  ! actual_Nitrogen uptake                 
                  

                    cpatch%actual_nitrogen_uptake(ico) = min(total_N_supply,               &
                                                         cpatch%shadow_N_uptake(ico)) 

                     if(cpatch%shadow_N_uptake(ico) == 0)then
                      cpatch%actual_nitrogen_uptake(ico) = cpatch%N_fixation(ico)
                     endif

                   !------------------------------------------------------------------------! 
!print*, 'soil + fixation uptake determined',                           &
!'cpatch%actual_nitrogen_uptake(ico)',cpatch%actual_nitrogen_uptake(ico), &
!'total_N_supply',total_N_supply
                   !------------------------------------------------------------------------! 
                  
                  !If the soil and fixation do not meet the shadow N uptake add extra C
                  !to excess C pool. The model ran for 17 hrs and did not come across
                  !a situation where this occurs.The code has not been checked but it
                  !should not make much of a difference in the long term 
                 
                  if (cpatch%shadow_N_uptake(ico) >  cpatch%actual_nitrogen_uptake(ico)     &
                  .and. cpatch%phenology_status(ico) == 1 ) then

!print*,    'soil and fixation not enough to meet demand',             &
!'cpatch%excess_carbon(ico)1',cpatch%excess_carbon(ico)  ,             &
!'cpatch%actual_nitrogen_uptake(ico)',cpatch%actual_nitrogen_uptake(ico) 
 
                    cpatch%excess_carbon(ico) = cpatch%excess_carbon(ico) +                  &
                           (cpatch%shadow_N_uptake(ico)- cpatch%actual_nitrogen_uptake(ico)) &
                              *c2n_leaf(ipft)  
    
! print*, 'cpatch%excess_carbon(ico)2',cpatch%excess_carbon(ico),        &
!  'difference', (cpatch%shadow_N_uptake(ico)                            &
!- cpatch%actual_nitrogen_uptake(ico)) *c2n_leaf(ipft)
!stop
                  endif

                  ! If fixation does satisfy N growth demand, the plant gets all of the N 
                  !it fixed and does not take up N from the soil

                   if(cpatch%shadow_N_uptake(ico) <= cpatch%N_fixation(ico)) then 
!print*, 'fixation satisfies demand',                                                     &
!' cpatch%actual_nitrogen_uptake(ico)', cpatch%actual_nitrogen_uptake(ico)
                  endif                                     
           

                  if (cpatch%shadow_N_uptake(ico) < cpatch%N_fixation(ico)                &
                  .and. cpatch%nstorage(ico) ==  nstorage_max) then
! print*,'fixed extra and nstorage = Nstorage_max'  ,                                     &
! 'cpatch%nstorage(ico)',cpatch%nstorage(ico),'nstorage_max',nstorage_max,                &
! 'csite%fsn_in(ipa) before fixation input', csite%fsn_in(ipa)         

                  csite%fsn_in(ipa) = csite%fsn_in(ipa)                                    &
                                   + (cpatch%N_fixation(ico)- cpatch%shadow_N_uptake(ico)) &
                                     * cpatch%nplant(ico)
!print*,'after fsn input',  csite%fsn_in(ipa)

                  endif


                  ! If the plant fixs more N then its demand for growth, add the excess N to
                  ! N storage up until it reaches nstorage_max. Add any remaining N to fsn_in.
                  
                   Nstorage_before_fixation = 0
                   Nstorage_before_fixation =  cpatch%nstorage(ico)

                   if (cpatch%shadow_N_uptake(ico) < cpatch%N_fixation(ico)                &
                   .and. cpatch%nstorage(ico) <  nstorage_max)then
 
!print*,' more N fixed then growth demand',                                                &
!' Nstorage_before_fixation', Nstorage_before_fixation,                                    &
!'nstorage_max',nstorage_max

                       cpatch%nstorage(ico) = min(nstorage_max,cpatch%nstorage(ico)        &
                                   + (cpatch%N_fixation(ico) - cpatch%shadow_N_uptake(ico)))
 
                   !Determine how much N goes to fsn_in and add it to the fast pool
                   Nstorage_after_fixation = 0
                   N_to_storage = 0
                   N_fixed_Extra = 0 

                   Nstorage_after_fixation =  cpatch%nstorage(ico)
                   N_to_storage = Nstorage_after_fixation -  Nstorage_before_fixation
                   N_fixed_Extra = cpatch%N_fixation(ico)- cpatch%shadow_N_uptake(ico)   &
                                  - N_to_storage

!print*,' Nstorage_after fixation', Nstorage_after_fixation,                               &
!'  N_to_storage', N_to_storage, ' N_fixed_Extra', N_fixed_Extra,'N fixed extra *nplant', &
!N_fixed_Extra* cpatch%nplant(ico)

!print*, ' csite%fsn_in(ipa) before fixed input', csite%fsn_in(ipa)
                  csite%fsn_in(ipa) = csite%fsn_in(ipa) +  N_fixed_Extra* cpatch%nplant(ico)
!print* ,'csite%fsn_in(ipa)2',csite%fsn_in(ipa)           
 
                   endif
                                
                  ! Calculate N limitation factor 
                  
                  cpatch%fsn(ico) =  total_N_supply                                       &
                                  / ( total_N_supply                                      &
                                  + cpatch% N_uptake_pot(ico))
!print*,'cpatch%fsn(ico)',cpatch%fsn(ico),'total_N_supply',total_N_supply

                  !Reset cpatch%nitrogen_supply(ico) to total N supply because this value
                  ! is used below to calculate how much N it can put into nstorage
                  ! as well as in the structural growth file to calculate how much nitrogen
                  ! a plant has (in addition to its storage) to put towards structural 
                  !growth and reproduction. 
 
                   cpatch%nitrogen_supply(ico) =   total_N_supply

!print*, 'NSUPPLY END', cpatch%nitrogen_supply(ico), &
! 'actual N uptake end',cpatch%actual_nitrogen_uptake(ico),&
! 'excess C end',cpatch%excess_carbon(ico)

!                   print*,'OBLIGATE FIXER END'

  
                  !------------------------------------------------------------------------! 
                  end if ! Fixers and non-fixers that are n limited 
                  !------------------------------------------------------------------------! 

                 !------------------------------------------------------------------------!
                  !      Allocate plant carbon balance to balive and bstorage.             !!JL,XXT, moved
                  !------------------------------------------------------------------------!
                  balive_in = cpatch%balive(ico)
                
                  cpatch%carbon_balance(ico) =  cpatch%carbon_balance(ico)                 &
                                             - cpatch%excess_carbon(ico)!JL!
                  carbon_balance = carbon_balance - cpatch%excess_carbon(ico)!JL!

                  call alloc_plant_c_balance(csite,ipa,ico,salloc,salloci,carbon_balance   &
                                            ,nitrogen_uptake                               &
                                            ,cpoly%green_leaf_factor(ipft,isi))
                  !Checks
                  !print*,'delta N in Balive', (cpatch%balive(ico) - balive_in)            &
                  !  / c2n_leaf(ipft),'nitrogen uptake',nitrogen_uptake
                  !print*,'ipft',ipft,'nitrogen uptake',nitrogen_uptake
                  !------------------------------------------------------------------------! 

                  !------------------------------------------------------------------------! 
                  ! Fill Nstorage to max_Nstorage if plant has access to more N then it    !
                  !needs for growth                                                        !
                  !------------------------------------------------------------------------!                       

                  ! Storage input will be the min of what it needs to reach its max or what
                  ! it has access to after n uptake for growth
                  ! NOTE: Nstorage fill is based on just Nstorage.. 
                  ! not cpatch%nstorage(ico)+cpatch%bleaf(ico)/c2n_leaf(ipft).
                  ! if a plant looses its leaves to pest attack it can resprout sometimes
                  ! 1 or 2 times so 1.3x crown
  
                   if (cpatch%nitrogen_supply(ico) > cpatch%actual_nitrogen_uptake(ico)    &
                   .and. cpatch%nstorage(ico)< nstorage_max .and. carbon_balance>0) then

                        nitrogen_uptake =  nitrogen_uptake                                 &
                                        + min((nstorage_max-cpatch%nstorage(ico)),         &
                                        (cpatch%nitrogen_supply(ico)                       &
                                        - cpatch%actual_nitrogen_uptake(ico)))
                        cpatch%nstorage(ico) = cpatch%nstorage(ico)                        &
                                             +  min((nstorage_max-cpatch%nstorage(ico))    &
                                             ,(cpatch%nitrogen_supply(ico)                 &
                                             -cpatch%actual_nitrogen_uptake(ico)))
                   endif

                  !------------------------------------------------------------------------! 
                  !------------------------------------------------------------------------! !JL, XXT
                  !  Increment the [kgN/m2] taken up during previous day.                  !
                  !------------------------------------------------------------------------!
                  csite%total_plant_nitrogen_uptake(ipa) =                                 &
                                       csite%total_plant_nitrogen_uptake(ipa)              &
                                       + nitrogen_uptake * cpatch%nplant(ico)  
                  !------------------------------------------------------------------------! 
                  ! Calculate Ecosystem N fixation [kgN/m2]                                !
                  !------------------------------------------------------------------------!


                  DailyForestFixation = DailyForestFixation + cpatch%N_fixation(ico) *     &
                                        csite%area(ipa) * cpatch%nplant(ico)



                  !Checks 
                  !print*,'BALIVE', 'ipa',ipa,'ipft',ipft,'ico',ico,                        &
                  !'cpatch%fixation_demand(ico)' ,cpatch%fixation_demand(ico),              &
                  !'cpatch%N_fixation(ico)',cpatch%fixation_demand(ico)
                  !print*,'ico',ico,'END NPLANT',sum(cpatch%nplant * (cpatch%balive  /      &
                  !c2n_leaf(cpatch%pft(ico)) + cpatch%bdead / c2n_stem(cpatch%pft(ico))     &
                  !   + cpatch%nstorage)) * csite%area(ipa)
                  !print*,'ico',ico,'END BALIVE',sum(cpatch%nplant * (cpatch%balive  /      &
                  !c2n_leaf(cpatch%pft(ico)))) * csite%area(ipa)
                  !print*,'ico',ico,'END NSTORAGE',sum(cpatch%nplant * cpatch%nstorage)     &
                  ! * csite%area(ipa)
                  ! print*,'ico',ico,'END N UPTAKE',csite%total_plant_nitrogen_uptake(ipa)
                  !------------------------------------------------------------------------!


                  !------------------------------------------------------------------------!
                  !      Do mortality --- note that only frost mortality changes daily.    !
                  !------------------------------------------------------------------------!
                  call mortality_rates(cpatch,ipa,ico,csite%avg_daily_temp(ipa)            &
                                      ,csite%age(ipa))
                  dndt = - sum(cpatch%mort_rate(:,ico)) * cpatch%nplant(ico) * tfact

                  !------- Update monthly mortality rate [plants/m2/month]. ---------------!
                  cpatch%monthly_dndt(ico) = cpatch%monthly_dndt(ico) + dndt

              
                   !----- Updating LAI, WPA, and WAI. --------------------------------------!
                  call area_indices(cpatch%nplant(ico),cpatch%bleaf(ico)                   &
                                   ,cpatch%bdead(ico),cpatch%balive(ico),cpatch%dbh(ico)   &
                                   ,cpatch%hite(ico) ,cpatch%pft(ico),cpatch%sla(ico)      &
                                   ,cpatch%lai(ico),cpatch%wpa(ico),cpatch%wai(ico)        &
                                   ,cpatch%crown_area(ico),cpatch%bsapwood(ico))

                  !----- Update above-ground biomass. -------------------------------------!
                  cpatch%agb(ico) = ed_biomass(cpatch%bdead(ico),cpatch%balive(ico)        &
                                              ,cpatch%bleaf(ico),cpatch%pft(ico)           &
                                              ,cpatch%hite(ico),cpatch%bstorage(ico)       &
                                              ,cpatch%bsapwood(ico))

                  !------------------------------------------------------------------------!
                  !     It is likely that biomass has changed, therefore, update           !
                  ! vegetation energy and heat capacity.                                   !
                  !------------------------------------------------------------------------!
                  old_leaf_hcap         = cpatch%leaf_hcap(ico)
                  old_wood_hcap         = cpatch%wood_hcap(ico)

                  call calc_veg_hcap(cpatch%bleaf(ico) ,cpatch%bdead(ico)                  &
                                    ,cpatch%bsapwood(ico),cpatch%nplant(ico)               &
                                    ,cpatch%pft(ico)                                       &
                                    ,cpatch%leaf_hcap(ico),cpatch%wood_hcap(ico))

 
                  call update_veg_energy_cweh(csite,ipa,ico,old_leaf_hcap,old_wood_hcap)

                 !----- Update the stability status. -------------------------------------!
                  call is_resolvable(csite,ipa,ico,cpoly%green_leaf_factor(:,isi))

               end do
                   !----- Update litter. ----------------------------------------------------!
                    call litter(csite,ipa)
                   ! This is where the N that was not able to be resorbed from leaves and fine
                   ! roots (52%) goes into the fsn_in pool JL

                   !------------------------------------------------------------------------!
                     
               !----- Update patch LAI, WAI, height, roughness... -------------------------!
               call update_patch_derived_props(csite,cpoly%lsl(isi),cpoly%met(isi)%prss,ipa)

               !----- Recalculate storage terms (for budget assessment). ------------------!
               call update_budget(csite,cpoly%lsl(isi),ipa,ipa)

               !----- It's a new day, reset average daily temperature. --------------------!
               csite%avg_daily_temp(ipa) = 0.0 
            end do
         end do
      end do

                  cgrid%total_N_Fixation(1)  = cgrid%total_N_Fixation(1) + DailyForestFixation !JL!
                  DailyForestFixation= 0.

      return
   end subroutine dbalive_dt
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This subroutine will compute the respiration terms other than leaf                !
   ! respiration, plus the carbon balance and maintenance costs but without                !
   ! updating the pools.                                                                   !
   !---------------------------------------------------------------------------------------!
   subroutine dbalive_dt_eq_0(cgrid, tfact)
      use ed_state_vars   , only : edtype                 & ! structure
                                 , polygontype            & ! structure
                                 , sitetype               & ! structure
                                 , patchtype              ! ! structure
      use pft_coms        , only : q                      & ! intent(in)
                                 , qsw                    & ! intent(in)
                                 , plant_N_supply_scale   & ! intent(in)
                                 , c2n_storage            & ! intent(in)
                                 , growth_resp_factor     & ! intent(in)
                                 , storage_turnover_rate  & ! intent(in)
                                 , phenology              ! ! intent(in)
      use physiology_coms , only : N_plant_lim            ! ! intent(in)
      use grid_coms       , only : nzg                    ! ! intent(in)
      use ed_therm_lib    , only : calc_veg_hcap          & ! function
                                 , update_veg_energy_cweh ! ! function
      use allometry       , only : area_indices           & ! subroutine
                                 , ed_biomass             ! ! function
      use mortality       , only : mortality_rates        ! ! subroutine
      use phenology_coms  , only : theta_crit             ! ! intent(in)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type(edtype)     , target     :: cgrid
      real             , intent(in) :: tfact
      !----- Local variables. -------------------------------------------------------------!
      type(polygontype), pointer    :: cpoly
      type(sitetype)   , pointer    :: csite
      type(patchtype)  , pointer    :: cpatch
      integer                       :: ipy
      integer                       :: isi
      integer                       :: ipa
      integer                       :: ico
      integer                       :: ipft
      real                          :: salloc
      real                          :: salloci
      real                          :: bl
      real                          :: br
      real                          :: daily_C_gain
      real                          :: carbon_balance
      real                          :: carbon_balance_pot
      real                          :: carbon_balance_max
      real                          :: balive_in
      real                          :: nitrogen_supply
      real                          :: dndt
      real                          :: old_leaf_hcap
      real                          :: old_wood_hcap
      real                          :: nitrogen_uptake
      real                          :: N_uptake_pot
      !------------------------------------------------------------------------------------!


      do ipy = 1,cgrid%npolygons
         cpoly => cgrid%polygon(ipy)

         do isi = 1,cpoly%nsites
            csite => cpoly%site(isi)

            do ipa = 1,csite%npatches
               cpatch => csite%patch(ipa)

               !----- Reset averaged variables. -------------------------------------------!
               csite%total_plant_nitrogen_uptake(ipa) = 0.0

               !----- Loop over cohorts. --------------------------------------------------!
               do ico = 1,cpatch%ncohorts

                  !----- Alias for current PFT. -------------------------------------------!
                  ipft = cpatch%pft(ico)

                  !----- Update the elongation factor. ------------------------------------!
                  select case (phenology(ipft))
                  case (4)
                     cpatch%elongf(ico) = max(0.0, min(1.0, cpatch%paw_avg(ico)/theta_crit))
                  case default
                     cpatch%elongf(ico) = 1.0
                  end select

                  !----- Initialize cohort nitrogen uptake. -------------------------------!
                  nitrogen_uptake = 0.0
                  N_uptake_pot    = 0.0

                  !----- Set allocation factors. ------------------------------------------!
                  salloc  = 1.0 + qsw(ipft) * cpatch%hite(ico) + q(ipft)
                  salloci = 1.0 / salloc
                  
                  !----- Leaf and root biomass. -------------------------------------------!
                  bl = cpatch%bleaf(ico)
                  br = cpatch%broot(ico)

                  !------------------------------------------------------------------------!
                  !     Compute maintenance costs.                                         !
                  !------------------------------------------------------------------------!
                  call plant_maintenance(cpatch,ico,br,bl,tfact,daily_C_gain               &
                                        ,csite%avg_daily_temp(ipa))

                  !----- Subtract maintenance costs from balive. --------------------------!
                  cpatch%cb(13,ico)     = cpatch%cb(13,ico)                                &
                                        - cpatch%leaf_maintenance(ico)                     &
                                        - cpatch%root_maintenance(ico)
                  cpatch%cb_max(13,ico) = cpatch%cb_max(13,ico)                            &
                                        - cpatch%leaf_maintenance(ico)                     &
                                        - cpatch%root_maintenance(ico)

                  !------------------------------------------------------------------------!
                  !      Calculate actual, potential and maximum carbon balances.          !
                  !------------------------------------------------------------------------!
                  call plant_carbon_balances(cpatch,ipa,ico,daily_C_gain,carbon_balance    &
                                            ,carbon_balance_pot,carbon_balance_max)

                  !------------------------------------------------------------------------!
                  !      Compute respiration rates for coming day [kgC/plant/day].         !
                  !------------------------------------------------------------------------!
                  cpatch%growth_respiration(ico)  = max(0.0, daily_C_gain                  &
                                                           * growth_resp_factor(ipft))
                  cpatch%storage_respiration(ico) = cpatch%bstorage(ico)                   &
                                                  * storage_turnover_rate(ipft) * tfact
                  cpatch%vleaf_respiration(ico) =                                          &
                                        (1.0 - cpoly%green_leaf_factor(ipft,isi))          &
                                      / (1.0 + q(ipft) + qsw(ipft) * cpatch%hite(ico))     &
                                      * cpatch%balive(ico) * storage_turnover_rate(ipft)   &
                                      * tfact

                  !------------------------------------------------------------------------!
                  !     Do a shadow calculation to see what would have happened if stomata !
                  ! were open.  This is used to calculate potential nitrogen uptake,       !
                  ! N_uptake_pot.                                                          !
                  !------------------------------------------------------------------------!
                  if (N_plant_lim == 1) then
                     call potential_N_uptake(cpatch,ico,salloc,salloci,balive_in           &
                                            ,carbon_balance_pot,N_uptake_pot               &
                                            ,cpoly%green_leaf_factor(ipft,isi))
                  end if

                  !------------------------------------------------------------------------!
                  !  Increment the [kgN/m2] taken up during previous day.                  !
                  !------------------------------------------------------------------------!
                 ! csite%total_plant_nitrogen_uptake(ipa) =                                 &
                 !                                 csite%total_plant_nitrogen_uptake(ipa)  &
                 !                              + nitrogen_uptake * cpatch%nplant(ico)

                  !----- Calculate plant N limitation factor. -----------------------------!
                  if (n_plant_lim == 0 .or. N_uptake_pot <= 0.0) then
                     cpatch%fsn(ico) = 1.0
                  else
                     nitrogen_supply = plant_N_supply_scale * br                           &
                                     * csite%mineralized_soil_N(ipa)
                     cpatch%fsn(ico) = nitrogen_supply / (nitrogen_supply + N_uptake_pot)
                  end if
                  
                  !------------------------------------------------------------------------!
                  !      Do mortality --- note that only frost mortality changes daily.    !
                  !------------------------------------------------------------------------!
                  call mortality_rates(cpatch,ipa,ico,csite%avg_daily_temp(ipa)            &
                                      ,csite%age(ipa))
               end do

               !----- It's a new day, reset average daily temperature. --------------------!
               csite%avg_daily_temp(ipa) = 0.0 
            end do
         end do
      end do

      return
   end subroutine dbalive_dt_eq_0
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !    This subroutine will transfer some of the stored carbon to balive in order to put  !
   ! the plant back on allometry.                                                          !
   !---------------------------------------------------------------------------------------!
   subroutine transfer_C_from_storage(cpatch,ico,salloc,salloci,nitrogen_uptake            &
                                     ,N_uptake_pot)
      use ed_state_vars , only : patchtype
      use pft_coms      , only : c2n_leaf    & ! intent(in)
                               , c2n_storage & ! intent(in)
                               , c2n_stem    & ! intent(in)
                               , q           & ! intent(in)
                               , qsw         ! ! intent(in)
      use decomp_coms   , only : f_labile    ! ! intent(in)
      use allometry     , only : dbh2bl      ! ! function
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type(patchtype), target        :: cpatch
      integer        , intent(in)    :: ico
      real           , intent(in)    :: salloc
      real           , intent(in)    :: salloci
      real           , intent(inout) :: nitrogen_uptake
      real           , intent(inout) :: N_uptake_pot
      !----- Local variables. -------------------------------------------------------------!
      integer                        :: ipft
      real                           :: off_allometry_cb
      real                           :: increment
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Only do the transfer if leaves exist.                                          !
      !------------------------------------------------------------------------------------!
      if (cpatch%phenology_status(ico) == 2) return
     
      !----- Alias for pft type. ----------------------------------------------------------!
      ipft = cpatch%pft(ico)
     
      !----- Determine how much biomass we need to go back to allometry. ------------------!
      off_allometry_cb = dbh2bl(cpatch%dbh(ico),ipft) * salloc - cpatch%balive(ico)

      !----- If plants have storage, transfer it to balive. -------------------------------!
      increment            = max(0.0,min(max(0.0, off_allometry_cb),cpatch%bstorage(ico)))
      cpatch%balive(ico)   = cpatch%balive(ico)   + increment
      cpatch%bstorage(ico) = cpatch%bstorage(ico) - increment

      !----- Compute sapwood and fine root biomass. ---------------------------------------!
      cpatch%broot(ico)    = q(ipft) * cpatch%balive(ico) * salloci
      cpatch%bsapwood(ico) = qsw(ipft) * cpatch%hite(ico) * cpatch%balive(ico) * salloci

      !------------------------------------------------------------------------------------!
      !      N uptake is required since c2n_leaf < c2n_storage.  Units are kgN/plant/day.  !
      !------------------------------------------------------------------------------------!
      nitrogen_uptake = increment * (        f_labile(ipft)  / c2n_leaf(ipft)              &
                                    + (1.0 - f_labile(ipft)) / c2n_stem(ipft)              &
                                    -  1.0 / c2n_storage)
      N_uptake_pot    = nitrogen_uptake

      return
   end subroutine transfer_C_from_storage
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   subroutine plant_maintenance(cpatch,ico,br,bl,tfact,daily_C_gain,tempk)
      use ed_state_vars, only : patchtype          ! ! structure
      use pft_coms     , only : phenology          & ! intent(in)
                              , root_turnover_rate & ! intent(in)
                              , leaf_turnover_rate ! ! intent(in)
      use consts_coms  , only : umol_2_kgC         & ! intent(in)
                              , day_sec            ! ! intent(in)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type(patchtype), target        :: cpatch
      integer        , intent(in)    :: ico
      real           , intent(in)    :: br
      real           , intent(in)    :: bl
      real           , intent(in)    :: tfact
      real           , intent(in)    :: tempk
      real           , intent(out)   :: daily_C_gain
      !----- Local variables. -------------------------------------------------------------!
      integer                        :: ipft
      real                           :: maintenance_temp_dep
      !------------------------------------------------------------------------------------!

      !------ Alias for plant functional type. --------------------------------------------!
      ipft = cpatch%pft(ico)

      !------ Get the temperature dependence. ---------------------------------------------!
      if (phenology(ipft) == 0) then
         maintenance_temp_dep = 1.0 / (1.0 + exp(0.4 * (278.15 - tempk)))
      else
         maintenance_temp_dep = 1.0
      end if

      !----- Calculate maintenance demand (kgC/plant/year). -------------------------------!
      cpatch%root_maintenance(ico) = root_turnover_rate(ipft) * br * maintenance_temp_dep
      if (phenology(ipft) /= 3)then
         cpatch%leaf_maintenance(ico) = leaf_turnover_rate(ipft) * bl * maintenance_temp_dep
      else
         cpatch%leaf_maintenance(ico) = leaf_turnover_rate(ipft) * bl                      &
                                      * cpatch%turnover_amp(ico) * maintenance_temp_dep
      end if


      !----- Convert units of maintenance to [kgC/plant/day]. -----------------------------!
      cpatch%leaf_maintenance(ico) = cpatch%leaf_maintenance(ico) * tfact
      cpatch%root_maintenance(ico) = cpatch%root_maintenance(ico) * tfact


      !----- Compute daily C uptake [kgC/plant/day]. --------------------------------------!
      if(cpatch%nplant(ico) > tiny(1.0)) then
         daily_C_gain = umol_2_kgC * day_sec * ( cpatch%today_gpp(ico)                     &
                                               - cpatch%today_leaf_resp(ico)               &
                                               - cpatch%today_root_resp(ico))              &
                                             / cpatch%nplant(ico)
      else
         daily_C_gain = 0.0
      end if

      return


   end subroutine plant_maintenance
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   subroutine plant_carbon_balances(cpatch,ipa,ico,daily_C_gain,carbon_balance             &
                                   ,carbon_balance_pot,carbon_balance_max)
      use ed_state_vars, only : patchtype          ! ! structure
      use pft_coms     , only : growth_resp_factor ! ! intent(in)
      use consts_coms  , only : umol_2_kgC         & ! intent(in)
                              , day_sec            ! ! intent(in)
      use ed_misc_coms , only : current_time       ! ! intent(in)
      use ed_max_dims  , only : n_pft              ! ! intent(in)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type(patchtype)          , target      :: cpatch
      integer                  , intent(in)  :: ipa
      integer                  , intent(in)  :: ico
      real                     , intent(in)  :: daily_C_gain
      real                     , intent(out) :: carbon_balance
      real                     , intent(out) :: carbon_balance_pot
      real                     , intent(out) :: carbon_balance_max
      !----- Local variables. -------------------------------------------------------------!
      real                                   :: daily_C_gain_pot
      real                                   :: daily_C_gain_max
      real                                   :: growth_respiration_pot
      real                                   :: growth_respiration_max
      integer                                :: ipft
      !----- Local constants. -------------------------------------------------------------!
      logical                  , parameter   :: print_debug = .false.
      !----- Locally saved variables. -----------------------------------------------------!
      logical, dimension(n_pft), save        :: first_time  = .true.
      !------------------------------------------------------------------------------------!

      !----- Alias for PFT type. ----------------------------------------------------------!
      ipft = cpatch%pft(ico)

      !------ Calculate actual daily carbon balance: kgC/plant/day. -----------------------!
      carbon_balance = daily_C_gain - cpatch%growth_respiration(ico)                       &
                                    - cpatch%vleaf_respiration(ico)

      if (cpatch%nplant(ico) > tiny(1.0)) then

         !---------------------------------------------------------------------------------!
         !      Calculate potential carbon balance (used for nitrogen demand function).    !
         ! [kgC/plant/day].                                                                !
         !---------------------------------------------------------------------------------!
         daily_C_gain_pot       = umol_2_kgC * day_sec * ( cpatch%today_gpp_pot(ico)       &
                                                         - cpatch%today_leaf_resp(ico)     &
                                                         - cpatch%today_root_resp(ico))    &
                                                       / cpatch%nplant(ico)
         growth_respiration_pot = max(0.0, daily_C_gain_pot * growth_resp_factor(ipft))

         carbon_balance_pot = daily_C_gain_pot - cpatch%growth_respiration(ico)            & !JL!
                                               - cpatch%vleaf_respiration(ico)
         
       !----- Calculate maximum carbon balance (used for mortality). --------------------!
         daily_C_gain_max       = umol_2_kgC * day_sec * ( cpatch%today_gpp_max(ico)       &
                                                         - cpatch%today_leaf_resp(ico)     &
                                                         - cpatch%today_root_resp(ico) )   &
                                                       / cpatch%nplant(ico)
         growth_respiration_max = max(0.0, daily_C_gain_max * growth_resp_factor(ipft))
         carbon_balance_max     = daily_C_gain_max - growth_respiration_max                &
                                                   - cpatch%vleaf_respiration(ico)
      else
         carbon_balance_max = 0.0
         carbon_balance_pot = 0.0
      end if

      !----- Carbon balances for mortality. -----------------------------------------------!
      cpatch%cb(13,ico)     = cpatch%cb(13,ico)     + carbon_balance
      cpatch%cb_max(13,ico) = cpatch%cb_max(13,ico) + carbon_balance_max

      if (print_debug) then

         if (first_time(ipft)) then
            first_time(ipft) = .false.
            write (unit=30+ipft,fmt='(a10,15(1x,a12))')                                    &
                '      TIME','       PATCH','      COHORT','      NPLANT','    CB_TODAY'   &
                            ,' GROWTH_RESP','  VLEAF_RESP','   TODAY_GPP','TODAY_GPPMAX'   &
                            ,'  TODAY_LEAF','  TODAY_ROOT',' CBMAX_TODAY','          CB'   &
                            ,'       CBMAX','  LEAF_MAINT','  ROOT_MAINT'
         end if

         write(unit=30+ipft,fmt='(2(i2.2,a1),i4.4,2(1x,i12),13(1x,es12.5))')               &
              current_time%month,'/',current_time%date,'/',current_time%year               &
             ,ipa,ico,cpatch%nplant(ico),carbon_balance,cpatch%growth_respiration(ico)     &
             ,cpatch%vleaf_respiration(ico),cpatch%today_gpp(ico)                          &
             ,cpatch%today_gpp_max(ico),cpatch%today_leaf_resp(ico)                        &
             ,cpatch%today_root_resp(ico),carbon_balance_max,cpatch%cb(13,ico)             &
             ,cpatch%cb_max(13,ico),cpatch%leaf_maintenance(ico)                           &
             ,cpatch%root_maintenance(ico)
      end if

      return

   end subroutine plant_carbon_balances
   !=======================================================================================!
   !=======================================================================================!







   !=======================================================================================!
   !=======================================================================================!
   subroutine alloc_plant_c_balance(csite,ipa,ico,salloc,salloci,carbon_balance            &
                                   ,nitrogen_uptake,green_leaf_factor)
      use ed_state_vars , only : sitetype     & ! structure
                               , patchtype    ! ! structure
      use pft_coms      , only : c2n_storage  & ! intent(in)
                               , c2n_leaf     & ! intent(in)
                               , sla          & ! intent(in)
                               , q            & ! intent(in)
                               , qsw          & ! intent(in)
                               , c2n_stem     ! ! intent(in)
      use decomp_coms   , only : f_labile     ! ! intent(in)
      use allometry     , only : dbh2bl       ! ! function
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type(sitetype) , target        :: csite
      integer        , intent(in)    :: ipa
      integer        , intent(in)    :: ico
      real           , intent(in)    :: salloc
      real           , intent(in)    :: salloci
      real           , intent(in)    :: carbon_balance
      real           , intent(inout) :: nitrogen_uptake
      real           , intent(in)    :: green_leaf_factor
      !----- Local variables. -------------------------------------------------------------!
      type(patchtype), pointer       :: cpatch
      integer                        :: ipft
      real                           :: bl_max
      real                           :: balive_max
      real                           :: bl_pot
      real                           :: increment
      real                           :: old_status
      real                           :: delta_bleaf
      real                           :: delta_broot
      real                           :: delta_bsapwood
      real                           :: available_carbon
      real                           :: f_total
      real                           :: f_bleaf
      real                           :: f_broot
      real                           :: f_bsapwood
      real                           :: f_resp
      real                           :: tr_bleaf
      real                           :: tr_broot
      real                           :: tr_bsapwood
      real                           :: bl
      logical                        :: on_allometry
      logical                        :: time_to_flush

      !------------------------------------------------------------------------------------!

      cpatch => csite%patch(ipa)
      
      ipft = cpatch%pft(ico) 

      !------------------------------------------------------------------------------------!
      !      When plants transit from dormancy to leaf flushing, it is possible that       !
      ! carbon_balance is negative, but the sum of carbon_balance and bstorage is          !
      ! positive. Under this circumstance, we have to allow plants to grow leaves.         !
      !------------------------------------------------------------------------------------!

      increment     = cpatch%bstorage(ico) + carbon_balance
      time_to_flush = (carbon_balance <= 0.0) .and. (increment > 0.0) .and.                &
                      (cpatch%phenology_status(ico) == 1) 

      if (carbon_balance > 0.0 .or. time_to_flush) then  !#1
         if (cpatch%phenology_status(ico) == 1) then   !#2
            !------------------------------------------------------------------------------!
            ! There are leaves, we are not actively dropping leaves and we're off          !
            ! allometry.  Here we will compute the maximum amount that can go to balive    !
            ! pools, and put any excess in storage.                                        !
            !------------------------------------------------------------------------------!
            !  available_carbon = cpatch%bstorage(ico) + carbon_balance
            ! for this section, excess carbon is included and will be subtracted from the
            ! balive portion  when the tr_bleaf ect are calculated
           available_carbon = cpatch%bstorage(ico) + carbon_balance + cpatch%excess_carbon(ico)

            !------------------------------------------------------------------------------!
            !     Maximum bleaf that the allometric relationship would allow.  If the      !
            ! plant is drought stress (elongf < 1), we do not allow the plant to get back  !
            ! to full allometry.                                                           !
            !------------------------------------------------------------------------------!
            bl_max     = dbh2bl(cpatch%dbh(ico),ipft) * green_leaf_factor                  &
                       * cpatch%elongf(ico)
            balive_max = dbh2bl(cpatch%dbh(ico),ipft) * salloc * cpatch%elongf(ico)

            !--- Amount that bleaf, broot, and bsapwood are off allometry -----------------!
            delta_bleaf = max (0.0, bl_max- cpatch%bleaf(ico))
            delta_broot = max (0.0, balive_max * q(ipft) * salloci - cpatch%broot(ico))
            delta_bsapwood = max (0.0, balive_max * qsw(ipft) * cpatch%hite(ico) * salloci &
                                     - cpatch%bsapwood(ico))

            !------------------------------------------------------------------------------!
            ! If the available carbon is less than what we need to get back to allometry.  !
            ! Grow pools in proportion to demand.  If we have enough carbon, we'll put the !
            ! extra into bstorage.                                                         !
            !------------------------------------------------------------------------------!
            
            f_bleaf    = delta_bleaf / bl_max
            f_broot    = delta_broot / (balive_max * q(ipft) * salloci )
            f_bsapwood = delta_bsapwood / (balive_max * qsw(ipft) * cpatch%hite(ico)       &
                       * salloci)
            f_total    = f_bleaf + f_broot + f_bsapwood

            !------------------------------------------------------------------------------!
            !     We only allow transfer from storage to living tissues if there is need   !
            ! to transfer.                                                                 !
            !------------------------------------------------------------------------------!
            if (f_total > 0.0) then
               tr_bleaf    = min( delta_bleaf   , (f_bleaf/f_total)    * available_carbon)
               tr_broot    = min( delta_broot   , (f_broot/f_total)    * available_carbon)
               tr_bsapwood = min( delta_bsapwood, (f_bsapwood/f_total) * available_carbon)
            else
               tr_bleaf    = 0.
               tr_broot    = 0.
               tr_bsapwood = 0.
            end if
            !------------------------------------------------------------------------------!
            !------------------------------------------------------------------------------!
            !XXT
            tr_bleaf = tr_bleaf - cpatch%excess_carbon(ico) * (f_bleaf/f_total)
            tr_broot = tr_broot - cpatch%excess_carbon(ico) * (f_broot/f_total)
            tr_bsapwood = tr_bsapwood - cpatch%excess_carbon(ico) * (f_bsapwood/f_total)

            cpatch%bleaf(ico)    = cpatch%bleaf(ico)    + tr_bleaf
            cpatch%broot(ico)    = cpatch%broot(ico)    + tr_broot
            cpatch%bsapwood(ico) = cpatch%bsapwood(ico) + tr_bsapwood
            !print*,'tr_balive',(tr_bleaf + tr_broot + tr_bsapwood)/c2n_leaf(cpatch%pft(ico))
            !print*,'balive_OLD',cpatch%balive(ico)/c2n_leaf(cpatch%pft(ico))
            cpatch%balive(ico)   = cpatch%bleaf(ico) + cpatch%broot(ico)                   &
                                 + cpatch%bsapwood(ico)
            !print*,'balive_NEW',cpatch%balive(ico)/c2n_leaf(cpatch%pft(ico))
    
            !----- NPP allocation in diff pools in KgC/m2/day. ----------------------------!
            cpatch%today_nppleaf(ico)   = tr_bleaf       * cpatch%nplant(ico)
            cpatch%today_nppfroot(ico)  = tr_broot       * cpatch%nplant(ico)
            cpatch%today_nppsapwood(ico)= tr_bsapwood    * cpatch%nplant(ico)
            cpatch%today_nppdaily(ico)  = carbon_balance * cpatch%nplant(ico)
            
            !------------------------------------------------------------------------------!
            !    Find the amount of carbon used to recover the tissues that were off-      !
            ! -allometry, take that from the carbon balance first, then use some of the    !
            ! storage if needed be.                                                        !
            !------------------------------------------------------------------------------!
            increment = carbon_balance -  tr_bleaf - tr_broot - tr_bsapwood
            cpatch%bstorage(ico) = max(0.0, cpatch%bstorage(ico) + increment) 
            !------------------------------------------------------------------------------!

           ! if (increment <= 0.0)  then
           !      nitrogen_uptake = nitrogen_uptake + max(0.0,min(cpatch%bstorage(ico)+    &
           !                carbon_balance,carbon_balance -increment))/  c2n_leaf(ipft)
                nitrogen_uptake = nitrogen_uptake + (carbon_balance -increment)/ c2n_leaf(ipft)
               !---------------------------------------------------------------------------!
               !    We are using up all of daily C gain and some of bstorage.  First       !
               ! calculate N demand from using daily C gain.                               !
               !---------------------------------------------------------------------------!
            !else
               !---------------------------------------------------------------------------!
               !     N uptake for fraction of daily C gain going to balive.                !
               !---------------------------------------------------------------------------!

             ! nitrogen_uptake = nitrogen_uptake +(carbon_balance - increment)/ c2n_leaf(ipft) 
 
               !----- N uptake for fraction of daily C gain going to bstorage. ------------!
              ! nitrogen_uptake = nitrogen_uptake + increment / c2n_storage 
              !JL. this is no longer necessray since bstorage is no longer linked to nstorage
           !endif

            on_allometry = 2.0 * abs(balive_max - cpatch%balive(ico))                      &
                         / (balive_max + cpatch%balive(ico))          < 1.e-6
            if (cpatch%elongf(ico) == 1.0 .and. on_allometry) then 
               !---------------------------------------------------------------------------!
               !     We're back to allometry, change phenology_status.                     !
               !---------------------------------------------------------------------------!
               cpatch%phenology_status(ico) = 0
            end if
         else !phenology is not 1 #2
            !------------------------------------------------------------------------------!
            !     Put carbon gain into storage.  If we're not actively dropping leaves or  !
            ! off-allometry, this will be used for structural growth at the end of the     !
            ! month.                                                                       !
            !------------------------------------------------------------------------------!
            cpatch%bstorage(ico) = cpatch%bstorage(ico) + carbon_balance
         !  cpatch%nstorage(ico) = cpatch%nstorage(ico) + carbon_balance/c2n_storage 
         !  !JL! plant is not getting N from carbon gain
         !  nitrogen_uptake      = nitrogen_uptake      + carbon_balance / c2n_storage
         !JL. this is no longer necessray since bstorage is no longer linked to nstorage                            
            !----- NPP allocation in diff pools in Kg C/m2/day. ---------------------------!
            cpatch%today_nppleaf(ico)    = 0.0
            cpatch%today_nppfroot(ico)   = 0.0
            cpatch%today_nppsapwood(ico) = 0.0
            cpatch%today_nppdaily(ico)   = carbon_balance * cpatch%nplant(ico)
         end if !phenology branch end #2
 

      else  !#1
         !---------------------------------------------------------------------------------!
         !   Carbon balance is negative, take it out of storage.                           !
         !---------------------------------------------------------------------------------!
         increment =  cpatch%bstorage(ico) + carbon_balance

         if (increment <= 0.0)  then !CB will consume all bstorage #3
                !nitrogen_uptake = nitrogen_uptake +                                       &
                !max(0.0,min(cpatch%bstorage(ico)+carbon_balance,carbon_balance -increment))  
           !----- Use Storage pool first then take out of balive. ------------------------!
            increment            =  - increment
            cpatch%bstorage(ico) = 0.0
           ! cpatch%nstorage(ico) = 0.0
           ! csite%fsn_in(ipa)    = csite%fsn_in(ipa) + cpatch%bstorage(ico) / c2n_storage  &
           !                                         * cpatch%nplant(ico)

           ! csite%fsn_in(ipa)    = csite%fsn_in(ipa) + cpatch%nstorage(ico)  &
           !                                         * cpatch%nplant(ico) !JL!
           !JL! C is being used up in bstorage but not necessarily n storage.
           ! Nstorage does not have to go to 0 unless the plant dies,
           ! therefore fsn_in is not updated here
            if (cpatch%phenology_status(ico) == 0)  then !#4
               !---------------------------------------------------------------------------!
               !     We were on allometry, but now we need to burn carbon and go off-      !
               ! -allometry.                                                               !
               !---------------------------------------------------------------------------!
               cpatch%balive(ico)   = cpatch%balive(ico) - increment
               cpatch%bleaf(ico)    = cpatch%balive(ico) * salloci * green_leaf_factor
               cpatch%broot(ico)    = cpatch%balive(ico) * q(ipft) * salloci
               cpatch%bsapwood(ico) = cpatch%balive(ico) * cpatch%hite(ico) * qsw(ipft)    &
                                    * salloci
               cpatch%phenology_status(ico) = 1
            else
               f_resp = cpatch%today_leaf_resp(ico)                                        &
                      / ( cpatch%today_leaf_resp(ico) + cpatch%today_root_resp(ico) )
               bl     = cpatch%bleaf(ico) - f_resp * (increment)

               if (bl > 0.0) then
                  cpatch%bleaf(ico) = bl
                  cpatch%broot(ico) = cpatch%broot(ico) - (1.0 - f_resp) * increment 
               else
                  cpatch%broot(ico) = cpatch%broot(ico) - (increment - cpatch%bleaf(ico))
                  cpatch%bleaf(ico) = 0.0
                  cpatch%elongf(ico) = 0.0
                  cpatch%phenology_status(ico) = 2
              end if 

               cpatch%balive(ico) = cpatch%bleaf(ico) + cpatch%broot(ico)                  &
                                  + cpatch%bsapwood(ico)    
          end if!#4

            csite%fsn_in(ipa) = csite%fsn_in(ipa) + increment                              &
                              * ( f_labile(ipft) / c2n_leaf(ipft)                          &
                                + (1.0 - f_labile(ipft)) / c2n_stem(ipft) )                &
                              * cpatch%nplant(ico)

         else !#3 !increment bigger then 0, bstorage large enough to satisfy neg CB
            !------ Burn the storage pool.  Dont' forget the nitrogen. --------------------!
            cpatch%bstorage(ico) =  cpatch%bstorage(ico) + carbon_balance
          !  cpatch%nstorage(ico) =  cpatch%nstorage(ico) + carbon_balance/c2n_storage   
          !   csite%fsn_in(ipa)    = csite%fsn_in(ipa) - carbon_balance / c2n_storage        &
          !                       * cpatch%nplant(ico)  
          ! CB is negative but there is enough carbon in bstorage for the plant to use 
          ! and this does not affect nstorge
 
         end if !#3 

         !---- NPP allocation in diff pools in KgC/m2/day. --------------------------------!
         cpatch%today_nppleaf(ico)    = 0.0
         cpatch%today_nppfroot(ico)   = 0.0
         cpatch%today_nppsapwood(ico) = 0.0
         cpatch%today_nppdaily(ico)   = carbon_balance * cpatch%nplant(ico)
      end if !#1

      return
   end subroutine alloc_plant_c_balance
   !=======================================================================================!
   !=======================================================================================!





   !=======================================================================================!
   !=======================================================================================!
   subroutine potential_N_uptake(cpatch,ico,salloc,salloci,balive_in,carbon_balance_pot    &
                                ,N_uptake_pot,green_leaf_factor)
      use ed_state_vars , only : patchtype    ! ! structure
      use pft_coms      , only : c2n_storage  & ! intent(in)
                               , c2n_leaf     & ! intent(in)
                               , c2n_stem     ! ! intent(in)
      use decomp_coms   , only : f_labile     ! ! intent(in)
      use allometry     , only : dbh2bl       ! ! intent(in)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type(patchtype), target        :: cpatch
      integer        , intent(in)    :: ico
      real           , intent(in)    :: salloc
      real           , intent(in)    :: salloci
      real           , intent(in)    :: balive_in
      real           , intent(in)    :: carbon_balance_pot
      real           , intent(inout) :: N_uptake_pot
      real           , intent(in)    :: green_leaf_factor
      !----- Local variables. -------------------------------------------------------------!
      integer                        :: ipft
      real                           :: bl_max
      real                           :: bl_pot
      real                           :: increment
      real                           :: nstorage_max !JL!
      real                           :: balive_max !JL!
      !------------------------------------------------------------------------------------!

      ipft = cpatch%pft(ico) 

      ! N_uptake_pot is 0 when CB is 0!JL
       if (carbon_balance_pot < 0.0 )then 
           N_uptake_pot = 0.0        
       end if                            

      if ( cpatch%phenology_status(ico) == 0 .and. carbon_balance_pot > 0.0 ) then

      !Check 
      !print *, 'Positive carbon balance with plants fully flushed' ,',phenstatus',&
      ! cpatch%phenology_status(ico)

         !----- Positive carbon balance with plants fully flushed. ------------------------!
         ! N_uptake_pot = N_uptake_pot + carbon_balance_pot / c2n_storage! JL not necessary  
        
      elseif (cpatch%phenology_status(ico) == 1) then

         !print *, 'Positive carbon balance plants growing' 

          !JL Calculate how much C will go into balive for full flush, not just bleaf
          balive_max = dbh2bl(cpatch%dbh(ico),ipft) * salloc * cpatch%elongf(ico)
          bl_pot = cpatch%balive(ico) + carbon_balance_pot

        ! bl_max = dbh2bl(cpatch%dbh(ico),ipft) * green_leaf_factor * cpatch%elongf(ico)
        !  bl_pot = cpatch%bleaf(ico) + carbon_balance_pot

         if (bl_pot > balive_max) then
            !------------------------------------------------------------------------------!
            !     This increment would take us over the limit, so we assign all that can   !
            ! go for leaves to them, and put the remainder in storage.                     !
            !------------------------------------------------------------------------------!
            !Increment = amount that would go into storage
            !increment    = carbon_balance_pot - (bl_max-cpatch%bleaf(ico)) !JL
 
            ! N_uptake_pot = N that accompanied the carbon that went into storage
            ! N_uptake_pot = N_uptake_pot + increment / c2n_storage !JL 
            
            ! amount that would go into storage
            increment    = carbon_balance_pot - (balive_max-cpatch%balive(ico)) 
 
           !JL the amount of growth if you have no limitations
            increment    = balive_max-cpatch%balive(ico)  
           !JL the amount o f N you use for the extra growth  
            N_uptake_pot = N_uptake_pot + increment                                        &
                         * ( f_labile(ipft) / c2n_leaf(ipft)                               &
                           + (1.0 - f_labile(ipft)) / c2n_stem(ipft)) 
           !Checks
           !if(ipft ==25 )then
           !print*, 'ipft',ipft,' bl_pot > balive_max', ' increment', increment,           &
           !'N_uptake_addition', increment + ( f_labile(ipft) / c2n_leaf(ipft)  +          &
           ! (1.0 - f_labile(ipft)) / c2n_stem(ipft)),'N uptake pot', N_uptake_pot
           ! endif

           !  if(N_uptake_pot==0)then
           !  stop
           !  endif

         elseif (carbon_balance_pot > 0.0) then

            !------------------------------------------------------------------------------!
            !      This increment did not exceed the limit, put everything in leaves.  We  !
            ! don't compute the uptake if carbon balance is negative, just because there   !
            ! will be no uptake...                                                         !
            !------------------------------------------------------------------------------!
            N_uptake_pot = N_uptake_pot + carbon_balance_pot                               &
                         * ( f_labile(ipft) / c2n_leaf(ipft)                               &
                           + (1.0 - f_labile(ipft)) / c2n_stem(ipft))                    
           !Checks
           !if(ipft ==25)then
           !print*, 'ipft',ipft,' bl_pot > balive_max everything in leaves', ' increment', &
           ! increment,'N_uptake_addition', increment + ( f_labile(ipft) / c2n_leaf(ipft)  &
           !  + (1.0 - f_labile(ipft)) / c2n_stem(ipft)),'N uptake pot', N_uptake_pot
           !endif

           !if(N_uptake_pot==0)then
           !print*,'ALL LEAVES'
           !stop
           !endif

         end if
      end if

     
    return
   end subroutine potential_N_uptake
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   subroutine litter(csite,ipa)

      use ed_state_vars, only : patchtype & ! structure
                              , sitetype  ! ! structure
      use pft_coms     , only : c2n_leaf  & ! intent(in)
                              , c2n_stem  & ! intent(in)
                              , l2n_stem  ! ! intent(in)
      use decomp_coms  , only : f_labile  ! ! intent(in)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type(sitetype)  , target     :: csite
      integer         , intent(in) :: ipa
      !----- Local variables. -------------------------------------------------------------!
      type(patchtype) , pointer    :: cpatch
      integer                      :: ico
      integer                      :: ipft
      real                         :: plant_litter
      real                         :: plant_litter_f
      real                         :: plant_litter_s
      real                         :: plant_litter_f_no_excess_c!JL
       !------------------------------------------------------------------------------------!

      cpatch => csite%patch(ipa)
   
      !------------------------------------------------------------------------------------!
      !      Add fine root and leaf turnover to the litter.                                !
      !------------------------------------------------------------------------------------!
      do ico=1,cpatch%ncohorts
         ipft = cpatch%pft(ico)

         plant_litter   = ( cpatch%leaf_maintenance(ico) + cpatch%root_maintenance(ico) )  &
                        * cpatch%nplant(ico)
         plant_litter_f = plant_litter * f_labile(ipft)                                    & 
                        + (cpatch%excess_carbon(ico)* cpatch%nplant(ico))                 !JL
         plant_litter_f_no_excess_c = plant_litter * f_labile(ipft)                       !JL
         !plant_litter_s = plant_litter - plant_litter_f                                  !JL
         plant_litter_s = plant_litter - plant_litter_f_no_excess_c                       !JL!

         csite%fsc_in(ipa) = csite%fsc_in(ipa) + plant_litter_f    
         csite%fsc_in_no_excess(ipa) =  csite%fsc_in_no_excess(ipa)                        &
                                     + plant_litter_f_no_excess_c                         !JL!

        ! the 48% of leaf and fine root turnover N available for resorption was already 
        !added to nstorage or put in fsn_in. The remaining 52% is put in fsn_in here..  

        ! csite%fsn_in(ipa) = csite%fsn_in(ipa) + plant_litter_f / c2n_leaf(ipft) 
        csite%fsn_in(ipa) = csite%fsn_in(ipa)                                              &
                          + ((plant_litter_f_no_excess_c / c2n_leaf(ipft))*0.52) !JL


         csite%ssc_in(ipa) = csite%ssc_in(ipa) + plant_litter_s
         csite%ssl_in(ipa) = csite%ssl_in(ipa) + plant_litter_s * l2n_stem / c2n_stem(ipft)

end do
      return
   end subroutine litter
   !=======================================================================================!
   !=======================================================================================!
   ! NEW SUBROUTINE., JL,XXT
   !=======================================================================================!
   !=======================================================================================!

subroutine shadow_n(csite,ipa,ico,salloc,salloci,carbon_balance            &
                                   ,nitrogen_uptake,green_leaf_factor)
      use ed_state_vars , only : sitetype     & ! structure
                               , patchtype    ! ! structure
      use pft_coms      , only : c2n_storage  & ! intent(in)
                               , c2n_leaf     & ! intent(in)
                               , sla          & ! intent(in)
                               , q            & ! intent(in)
                               , qsw          & ! intent(in)
                               , c2n_stem     ! ! intent(in)
      use decomp_coms   , only : f_labile     ! ! intent(in)
      use allometry     , only : dbh2bl       ! ! function
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type(sitetype) , target        :: csite
      integer        , intent(in)    :: ipa
      integer        , intent(in)    :: ico
      real           , intent(in)    :: salloc
      real           , intent(in)    :: salloci
      real           , intent(in)    :: carbon_balance
      real           , intent(inout) :: nitrogen_uptake
      real           , intent(in)    :: green_leaf_factor

!----- Local variables. -------------------------------------------------------------!
      type(patchtype), pointer       :: cpatch
      integer                        :: ipft
      real                           :: bl_max
      real                           :: balive_max
      real                           :: bl_pot
      real                           :: increment
      real                           :: old_status
      real                           :: delta_bleaf
      real                           :: delta_broot
      real                           :: delta_bsapwood
      real                           :: available_carbon
      real                           :: f_total
      real                           :: f_bleaf
      real                           :: f_broot
      real                           :: f_bsapwood
      real                           :: f_resp
      real                           :: tr_bleaf
      real                           :: tr_broot
      real                           :: tr_bsapwood
      real                           :: bl
      logical                        :: on_allometry
      logical                        :: time_to_flush
      !JL additions for shadow calculation!      
      real                           :: temp_bleaf
      real                           :: temp_bstorage
      real                           :: temp_broot
      real                           :: temp_bsapwood
      real                           :: temp_balive
      integer                        :: temp_phenology_status
      real                           :: temp_elongf
      real                           :: temp_hite
      real                           :: temp_dbh
      !------------------------------------------------------------------------------------!

      cpatch => csite%patch(ipa)
      
      ipft = cpatch%pft(ico) 

      !Initialize new vairables !
      temp_bleaf = cpatch%bleaf(ico)
      temp_bstorage = cpatch%bstorage(ico)
      temp_broot = cpatch%broot(ico)
      temp_bsapwood= cpatch%bsapwood(ico)
      temp_balive = cpatch%balive(ico)  
      temp_phenology_status = cpatch%phenology_status(ico) 
      temp_elongf=cpatch%elongf(ico)
      temp_dbh=cpatch%dbh(ico)
      temp_hite = cpatch%hite(ico)


      !------------------------------------------------------------------------------------!
      !      When plants transit from dormancy to leaf flushing, it is possible that       !
      ! carbon_balance is negative, but the sum of carbon_balance and bstorage is          !
      ! positive. Under this circumstance, we have to allow plants to grow leaves.         !
      !------------------------------------------------------------------------------------!
      increment     = temp_bstorage + carbon_balance
      time_to_flush = (carbon_balance <= 0.0) .and. (increment > 0.0) .and.                &
                      (temp_phenology_status == 1) 

      if (carbon_balance > 0.0 .or. time_to_flush) then 
         if (temp_phenology_status == 1) then
            !------------------------------------------------------------------------------!
            ! There are leaves, we are not actively dropping leaves and we're off          !
            ! allometry.  Here we will compute the maximum amount that can go to balive    !
            ! pools, and put any excess in storage.                                        !
            !------------------------------------------------------------------------------!
            available_carbon = temp_bstorage + carbon_balance

            !------------------------------------------------------------------------------!
            !     Maximum bleaf that the allometric relationship would allow.  If the      !
            ! plant is drought stress (elongf < 1), we do not allow the plant to get back  !
            ! to full allometry.                                                           !
            !------------------------------------------------------------------------------!
           
           
            bl_max     = dbh2bl(cpatch%dbh(ico),ipft) * green_leaf_factor                  &
                       * temp_elongf
            balive_max = dbh2bl(cpatch%dbh(ico),ipft) * salloc *  temp_elongf

            !--- Amount that bleaf, broot, and bsapwood are off allometry -----------------!
            delta_bleaf = max (0.0, bl_max- temp_bleaf)
            delta_broot = max (0.0, balive_max * q(ipft) * salloci - temp_broot)
            delta_bsapwood = max (0.0, balive_max * qsw(ipft) * temp_hite * salloci &
                                     - temp_bsapwood)

            !------------------------------------------------------------------------------!
            ! If the available carbon is less than what we need to get back to allometry.  !
            ! Grow pools in proportion to demand.  If we have enough carbon, we'll put the !
            ! extra into bstorage.                                                         !
            !------------------------------------------------------------------------------!
            
            f_bleaf    = delta_bleaf / bl_max
            f_broot    = delta_broot / (balive_max * q(ipft) * salloci )
            f_bsapwood = delta_bsapwood / (balive_max * qsw(ipft) * temp_hite      &
                       * salloci)
            f_total    = f_bleaf + f_broot + f_bsapwood

            !------------------------------------------------------------------------------!
            !     We only allow transfer from storage to living tissues if there is need   !
            ! to transfer.                                                                 !
            !------------------------------------------------------------------------------!
            if (f_total > 0.0) then
               tr_bleaf    = min( delta_bleaf   , (f_bleaf/f_total)    * available_carbon)
               tr_broot    = min( delta_broot   , (f_broot/f_total)    * available_carbon)
               tr_bsapwood = min( delta_bsapwood, (f_bsapwood/f_total) * available_carbon)
            else
               tr_bleaf    = 0.
               tr_broot    = 0.
               tr_bsapwood = 0.
            end if
            !------------------------------------------------------------------------------!
            temp_bleaf = temp_bleaf
            temp_bleaf    = temp_bleaf    + tr_bleaf
            temp_broot    = temp_broot    + tr_broot
            temp_bsapwood = temp_bsapwood + tr_bsapwood

            temp_balive   =temp_bleaf + temp_broot                   &
                                 + temp_bsapwood
            
            !----- NPP allocation in diff pools in KgC/m2/day. ----------------------------!
            !cpatch%today_nppleaf   = tr_bleaf       * cpatch%nplant(ico)
            !cpatch%today_nppfroot(ico)  = tr_broot       * cpatch%nplant(ico)
            !cpatch%today_nppsapwood(ico)= tr_bsapwood    * cpatch%nplant(ico)
            !cpatch%today_nppdaily(ico)  = carbon_balance * cpatch%nplant(ico)
            

            !------------------------------------------------------------------------------!
            !    Find the amount of carbon used to recover the tissues that were off-      !
            ! -allometry, take that from the carbon balance first, then use some of the    !
            ! storage if needed be.                                                        !
            !------------------------------------------------------------------------------!
            increment = carbon_balance -  tr_bleaf - tr_broot - tr_bsapwood
            temp_bstorage = max(0.0, temp_bstorage + increment) 
            !------------------------------------------------------------------------------!

           ! if (increment <= 0.0)  then
           !     nitrogen_uptake = nitrogen_uptake +                                       &
                                 !max(0.0,min(cpatch%bstorage(ico)+carbon_balance,         &
                                 ! carbon_balance -increment))/  c2n_leaf(ipft)
            nitrogen_uptake = nitrogen_uptake + (carbon_balance -increment)/ c2n_leaf(ipft)
               !---------------------------------------------------------------------------!
               !    We are using up all of daily C gain and some of bstorage.  First       !
               ! calculate N demand from using daily C gain.                               !
               !---------------------------------------------------------------------------!
               !if (carbon_balance < 0.0) then
               ! nitrogen_uptake = nitrogen_uptake + carbon_balance / c2n_storage
               ! nitrogen_uptake = nitrogen_uptake + increment / c2n_storage 
               ! nitrogen_uptake = nitrogen_uptake                                        &
               !                 + (carbon_balance - increment)                           &
               !                 * ( f_labile(ipft) / c2n_leaf(ipft)                      &
               !                 + (1.0 - f_labile(ipft)) / c2n_stem(ipft)                &
               !                 -  1.0 / c2n_storage)
             
               ! else
               !    nitrogen_uptake = nitrogen_uptake + carbon_balance                       &
               !                   * ( f_labile(ipft) / c2n_leaf(ipft)                       &
               !                     + (1.0 - f_labile(ipft)) / c2n_stem(ipft) )

                  !------------------------------------------------------------------------!
                  !     Now calculate additional N uptake required from transfer of C from !
                  ! storage to balive.                                                     !
                  !------------------------------------------------------------------------!
                ! nitrogen_uptake  = nitrogen_uptake +  ( - 1.0 * increment )              &
                !                   * ( f_labile(ipft)  / c2n_leaf(ipft)                   &
                !                     + (1.0 - f_labile(ipft)) / c2n_stem(ipft)            &
                 !                    -  1.0 / c2n_storage)

               !end if

            !else
               !---------------------------------------------------------------------------!
               !     N uptake for fraction of daily C gain going to balive.                !
               !---------------------------------------------------------------------------!
 
              ! nitrogen_uptake = nitrogen_uptake + (carbon_balance - increment)            &
              !                 * ( f_labile(ipft) / c2n_leaf(ipft)                         &
              !                   + (1.0 - f_labile(ipft)) / c2n_stem(ipft))

               !----- N uptake for fraction of daily C gain going to bstorage. ------------!
           !     nitrogen_uptake = nitrogen_uptake + increment / c2n_storage !JL
           !end if

            on_allometry = 2.0 * abs(balive_max - temp_balive)                      &
                         / (balive_max  + temp_balive)          < 1.e-6
            if (temp_elongf == 1.0 .and. on_allometry) then
               !---------------------------------------------------------------------------!
               !     We're back to allometry, change phenology_status.                     !
               !---------------------------------------------------------------------------!
             temp_phenology_status = 0
            end if
         else
            !------------------------------------------------------------------------------!
            !     Put carbon gain into storage.  If we're not actively dropping leaves or  !
            ! off-allometry, this will be used for structural growth at the end of the     !
            ! month.                                                                       !
            !------------------------------------------------------------------------------!
           temp_bstorage =temp_bstorage + carbon_balance
           ! nitrogen_uptake      = nitrogen_uptake      + carbon_balance / c2n_storage
                         
            !----- NPP allocation in diff pools in Kg C/m2/day. ---------------------------!
            !cpatch%today_nppleaf(ico)    = 0.0
            !cpatch%today_nppfroot(ico)   = 0.0
            !cpatch%today_nppsapwood(ico) = 0.0
            !cpatch%today_nppdaily(ico)   = carbon_balance * cpatch%nplant(ico)
         end if
 

      else
         !---------------------------------------------------------------------------------!
         !   Carbon balance is negative, take it out of storage.                           !
         !---------------------------------------------------------------------------------!
         increment =  temp_bstorage + carbon_balance

         if (increment <= 0.0)  then
            !----- Use Storage pool first then take out of balive. ------------------------!
            increment            =  - increment
            temp_bstorage = 0.0
           ! temp_fsn_in    =temp_fsn_in + temp_bstorage / c2n_storage  &
           !                                          * cpatch%nplant

            if (temp_phenology_status == 0)  then
               !---------------------------------------------------------------------------!
               !     We were on allometry, but now we need to burn carbon and go off-      !
               ! -allometry.                                                               !
               !---------------------------------------------------------------------------!
              temp_balive   = temp_balive - increment
              temp_bleaf    = temp_balive * salloci * green_leaf_factor
              temp_broot   =temp_balive * q(ipft) * salloci
              temp_bsapwood = temp_balive * temp_hite * qsw(ipft)    &
                                    * salloci
               temp_phenology_status = 1
            else
               f_resp = cpatch%today_leaf_resp(ico)                                        &
                      / ( cpatch%today_leaf_resp(ico) + cpatch%today_root_resp(ico) )
               bl     = temp_bleaf - f_resp * (increment)

               if (bl > 0.0) then
                 temp_bleaf = bl
                  temp_broot = temp_broot - (1.0 - f_resp) * increment 
               else
                 temp_broot =temp_broot - (increment - temp_bleaf)
                  temp_bleaf = 0.0
                 temp_elongf = 0.0
                  temp_phenology_status = 2
               end if

               temp_balive= temp_bleaf + temp_broot                &
                                  + temp_bsapwood   
            end if

           !temp_fsn_in = temp_fsn_in + increment                                          &
           !                   * ( f_labile(ipft) / c2n_leaf(ipft)                         &
           !                   + (1.0 - f_labile(ipft)) / c2n_stem(ipft) )                 &
           !                   * cpatch%nplant(ico)

         else
            !------ Burn the storage pool.  Dont' forget the nitrogen. --------------------!
           temp_bstorage =  temp_bstorage + carbon_balance
           ! temp_fsn_in(ipa)    = temp_fsn_in(ipa) - carbon_balance / c2n_storage        &
           !                      * cpatch%nplant(ico)   
         end if

         !---- NPP allocation in diff pools in KgC/m2/day. --------------------------------!
         !cpatch%today_nppleaf(ico)    = 0.0
         !cpatch%today_nppfroot(ico)   = 0.0
         !cpatch%today_nppsapwood(ico) = 0.0
         !cpatch%today_nppdaily(ico)   = carbon_balance * cpatch%nplant(ico)n
      end if

      return
   end subroutine shadow_n
!==========================================================================================!
!==========================================================================================!

end module growth_balive
!==========================================================================================!
!==========================================================================================!
