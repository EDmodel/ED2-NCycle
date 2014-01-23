!==========================================================================================!
!==========================================================================================!
!     This subroutine computes the soil respiration terms (root and heterotrophic).        !
!------------------------------------------------------------------------------------------!
subroutine soil_respiration(csite,ipa,mzg,ntext_soil)

   use ed_state_vars, only : sitetype                 & ! structure
                           , patchtype                ! ! structure
   use soil_coms    , only : soil                     ! ! intent(in)
   use pft_coms     , only : root_respiration_factor  ! ! intent(in)

   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   type(sitetype)                , target     :: csite
   integer                       , intent(in) :: ipa
   integer                       , intent(in) :: mzg
   integer       , dimension(mzg), intent(in) :: ntext_soil
   !----- Local variables. ----------------------------------------------------------------!
   type(patchtype)               , pointer    :: cpatch
   integer                                    :: ico
   integer                                    :: ipft
   real                                       :: r_resp_temp_fac
   real                                       :: Lc
   real                                       :: r_resp
   !----- External functions. -------------------------------------------------------------!
   real                          , external   :: resp_weight
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !      This is the temperature dependence of root respiration.  Same for all cohorts.   !
   !---------------------------------------------------------------------------------------!
   r_resp_temp_fac = 1.0                                                                   &
                   / (1.0 + exp(0.4 * ( 278.15 - csite%soil_tempk(mzg,ipa) ) ) )           &
                   / (1.0 + exp(0.4 * ( csite%soil_tempk(mzg,ipa) - 318.15 ) ) )           &
                   * exp( 10.41 - 3000.0/csite%soil_tempk(mzg,ipa) )

   cpatch => csite%patch(ipa)
   do ico = 1,cpatch%ncohorts
      ipft = cpatch%pft(ico)
      r_resp = root_respiration_factor(ipft) * r_resp_temp_fac * cpatch%broot(ico)        &
             * cpatch%nplant(ico)
      cpatch%root_respiration(ico) = r_resp
      cpatch%mean_root_resp(ico)   = cpatch%mean_root_resp(ico)  + r_resp
      cpatch%today_root_resp(ico)  = cpatch%today_root_resp(ico) + r_resp
   end do

   !----- Compute soil/temperature modulation of heterotrophic respiration. ---------------!
   csite%A_decomp(ipa) = resp_weight(csite%soil_tempk(mzg,ipa),csite%soil_water(mzg,ipa)   &
                                    ,ntext_soil(mzg))

   !----- Compute nitrogen immobilization factor. -----------------------------------------!
   call resp_f_decomp(csite,ipa, Lc)

   !----- Compute heterotrophic respiration. ----------------------------------------------!
   call resp_rh(csite,ipa, Lc)

   !----- Update averaged variables. ------------------------------------------------------!
   csite%today_A_decomp(ipa)  = csite%today_A_decomp(ipa)  + csite%A_decomp(ipa)
   csite%today_Af_decomp(ipa) = csite%today_Af_decomp(ipa)                                 &
                              + csite%A_decomp(ipa) * csite%f_decomp(ipa)
   csite%mean_rh(ipa)         = csite%mean_rh(ipa) + csite%rh(ipa)

   return
end subroutine soil_respiration
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This function computes the respiration limitation factor, which includes limitations !
! due to temperature and moisture.                                                         !
!------------------------------------------------------------------------------------------!
real function resp_weight(soil_tempk,soil_water,nsoil)

   use decomp_coms, only : resp_temperature_increase  & ! intent(in)
                         , resp_opt_water             & ! intent(in)
                         , resp_water_below_opt       & ! intent(in)
                         , resp_water_above_opt       & ! intent(in)
                         , LloydTaylor                ! ! intent(in)
   use soil_coms  , only : soil                       ! ! intent(in)

   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   real   , intent(in) :: soil_tempk
   real   , intent(in) :: soil_water
   integer, intent(in) :: nsoil
   !----- Local variables. ----------------------------------------------------------------!
   real                :: temperature_limitation
   real                :: water_limitation
   real                :: rel_soil_moist
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     Find the temperature dependence.                                                  !
   !---------------------------------------------------------------------------------------!
   if (LloydTaylor) then 
      !----- Use Lloyd and Taylor (1994) temperature dependence. --------------------------!
      temperature_limitation = min( 1.0                                                    &
                                  , resp_temperature_increase                              &
                                  * exp(308.56 * (1./56.02 - 1./(soil_tempk-227.15)) ) )
   else 
      !----- Use original exponential temperature dependence. -----------------------------!
      temperature_limitation = min( 1.0                                                    &
                                  , exp( resp_temperature_increase * (soil_tempk-294.)))
!      temperature_limitation = min( 1.0                                                    &
!                                  , exp( resp_temperature_increase * (soil_tempk-318.15)))
   end if


   !---------------------------------------------------------------------------------------!
   !     Find the relative soil moisture, then the moisture dependence.                    !
   !---------------------------------------------------------------------------------------!
   rel_soil_moist = (soil_water         - soil(nsoil)%soilcp)                              &
                  / (soil(nsoil)%slmsts - soil(nsoil)%soilcp)
   if (rel_soil_moist <= resp_opt_water)then
      water_limitation = exp( (rel_soil_moist - resp_opt_water) * resp_water_below_opt)
   else
      water_limitation = exp( (resp_opt_water - rel_soil_moist) * resp_water_above_opt)
   end if
   
   !----- Compute the weight, which is just the combination of both. ----------------------!
   resp_weight = temperature_limitation * water_limitation
      
   return
end function resp_weight
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This subroutine computes the Nitrogen immobilization factor.                         !
!------------------------------------------------------------------------------------------!
subroutine resp_f_decomp(csite,ipa,Lc)

   use ed_state_vars, only : sitetype               ! ! structure
   use decomp_coms  , only : r_stsc                 & ! intent(in)
                           , N_immobil_supply_scale & ! intent(in)
                           , K1                     & ! intent(in)
                           , n_decomp_lim           ! ! intent(in)
   use pft_coms     , only : c2n_structural         & ! intent(in)
                           , c2n_slow               ! ! intent(in)

   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   type(sitetype), target      :: csite
   integer       , intent(in)  :: ipa
   real          , intent(out) :: Lc
   !----- Local variables. ----------------------------------------------------------------!
   real                        :: N_immobilization_demand
   !---------------------------------------------------------------------------------------!

 
   if (csite%structural_soil_C(ipa) > 0.0) then
      if (csite%structural_soil_L(ipa) == csite%structural_soil_C(ipa)) then
         Lc = 0.049787 ! = exp(-3.0)
      else
         Lc = exp(-3.0 * csite%structural_soil_L(ipa)/csite%structural_soil_C(ipa))
      end if
   else
      Lc=0.0
   end if
   
   if (n_decomp_lim == 1) then
      N_immobilization_demand = csite%A_decomp(ipa) * Lc * K1                              &
                              * csite%structural_soil_C(ipa)                               &
                              * ((1.0 - r_stsc) / c2n_slow - 1.0 / c2n_structural)
      
      csite%f_decomp(ipa)     = N_immobil_supply_scale * csite%mineralized_soil_N(ipa)     &
                              / ( N_immobilization_demand                                  &
                                + N_immobil_supply_scale  * csite%mineralized_soil_N(ipa))
   else
      !----- Option for no plant N limitation. --------------------------------------------!
      csite%f_decomp(ipa)     = 1.0
   end if
  csite%f_decomp(ipa)     = 1.0 
       !csite%f_decomp(ipa) = 1.0 turns off immobilization/ whether decomposition can be
       ! limited by nitrogen. this can also be turned off in ED2IN(n_decomp_lim) !
  
       !write(*,*),'immobilization demand',N_immobilization_demand
 return
end subroutine resp_f_decomp
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This subroutine computes the heterotrophic respiration.                              !
!------------------------------------------------------------------------------------------!
subroutine resp_rh(csite,ipa,Lc)

   use ed_state_vars, only : sitetype       ! ! structure
   use consts_coms  , only : kgCday_2_umols ! ! intent(in)
   use decomp_coms  , only : K1             & ! intent(in)
                           , K2             & ! intent(in)
                           , K3             & ! intent(in)
                           , r_fsc          & ! intent(in)
                           , r_ssc          & ! intent(in)
                           , r_stsc         & ! intent(in)
                           , cwd_frac       ! ! intent(in)

   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   type(sitetype), target      :: csite
   integer       , intent(in)  :: ipa
   real          , intent(in)  :: Lc
   !----- Local variables. ----------------------------------------------------------------!
   real                        :: fast_C_loss
   real                        :: structural_C_loss
   real                        :: slow_C_loss
   !---------------------------------------------------------------------------------------!



   !----- The following variables have units of [kgC/m2/day]. -----------------------------!
   fast_C_loss       = csite%A_decomp(ipa) * K2 * csite%fast_soil_C(ipa)
   structural_C_loss = csite%A_decomp(ipa) * Lc * K1 * csite%structural_soil_C(ipa)        &
                     * csite%f_decomp(ipa)
   slow_C_loss       = csite%A_decomp(ipa) * K3 * csite%slow_soil_C(ipa) 

   !----- The following variables have units of [umol_CO2/m2/s]. --------------------------!
   csite%rh(ipa)     = kgCday_2_umols * ( r_fsc * fast_C_loss + r_stsc * structural_C_loss &
                                        + r_ssc * slow_C_loss)      
                       
   csite%cwd_rh(ipa) = kgCday_2_umols * (r_stsc * structural_C_loss + r_ssc * slow_C_loss) &
                     * cwd_frac

   return
end subroutine resp_rh
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This subroutine will update the soil carbon and nitrogen pools.                      !
!------------------------------------------------------------------------------------------!
subroutine update_C_and_N_pools(cgrid)
   
   use ed_state_vars, only : edtype          & ! structure
                           , polygontype     & ! structure
                           , sitetype        ! ! structure
   use decomp_coms  , only : K1              & ! intent(in)
                           , K2              & ! intent(in)
                           , K3              & ! intent(in)
                           , r_stsc          ! ! intent(in
   use pft_coms     , only : c2n_slow        & ! intent(in)
                           , c2n_structural  ! ! intent(in)
 use  soil_coms     , only : soil &          !the structure stores the soil information  !JL
                           , nslcon &        !An index                                   !JL
                           , dslz &          !soil thickness of each layer               !JL
                           , slz             ! soil depth of each layer                  !JL
  use grid_coms, only: nzg
 

 implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   type(edtype)     , target   :: cgrid
  ! real   , dimension(nzg), intent(in)    :: soil_water                                 !JL
   !----- Local variables. ----------------------------------------------------------------!
   type(polygontype), pointer  :: cpoly
   type(sitetype)   , pointer  :: csite
   integer                     :: ipy
   integer                     :: isi
   integer                     :: ipa
   real                        :: Lc
   real                        :: fast_C_loss
   real                        :: fast_N_loss 
   real                        :: structural_C_loss
   real                        :: structural_L_loss
   real                        :: slow_C_input
   real                        :: slow_C_loss
   real                        :: DON_loss                                               !JL
   real                        :: DON_leached                                            !JL
   real                        :: structural_N_loss                                      !JL
   real                        :: soil_depth                                             !JL
   real                        :: water_density                                          !JL
   real                        :: sec_in_day                                             !JL
   real                        :: N_leached                                              !JL
   real                        :: Nconc                                                  !JL
   real                        :: Avg_soil_water                                         !JL
   real                        :: N_deposition                                           !JL
   real                        :: efficiency_factor                                      !JL
   real                        :: water_input                                            !JL
   integer                     :: k            ! soil layer counter                      !JL              
   integer                     :: nsoil        ! soil class                              !JL
   real                        :: fracw        ! fraction of saturation (m/m)            !JL                   
   real                        :: total_soil_water                                       !JL
   real                        :: total_field_capacity                                   !JL
   real                        :: total_soil_water_capacity                              !JL
   real                        :: MSN_after_plant_uptake                                 !JL
   real                        :: N_loss_total                                           !JL
   real                        :: N_loss_emission                                        !JL
   real                        :: f_soluble                                              !JL
   !---------------------------------------------------------------------------------------!


   polygonloop: do ipy = 1,cgrid%npolygons

      cpoly => cgrid%polygon(ipy)

      siteloop: do isi = 1,cpoly%nsites
         
         csite => cpoly%site(isi)
   
         patchloop: do ipa = 1,csite%npatches

            if (csite%structural_soil_C(ipa) > 0.0) then
               if (csite%structural_soil_L(ipa) == csite%structural_soil_C(ipa)) then
                  Lc = 0.049787 ! = exp(-3.0)
               else
                  Lc = exp( -3.0 * csite%structural_soil_L(ipa)                            &
                          /  csite%structural_soil_C(ipa))
               end if
            else
               Lc=0.0
            end if
      
            !----- Fast pools. ------------------------------------------------------------!
            fast_C_loss = csite%today_A_decomp(ipa) * K2 * csite%fast_soil_C(ipa)
            fast_N_loss = csite%today_A_decomp(ipa) * K2 * csite%fast_soil_N(ipa)

            !----- Structural pools. ------------------------------------------------------!
            structural_C_loss = csite%today_Af_decomp(ipa) * Lc * K1                       &
                              * csite%structural_soil_C(ipa)
            structural_L_loss = csite%today_Af_decomp(ipa) * Lc * K1                       &
                              * csite%structural_soil_L(ipa)
            structural_N_loss =   structural_C_loss/c2n_structural     

            !----- Slow pools. ------------------------------------------------------------!
            slow_C_input = (1.0 - r_stsc) * structural_C_loss           
            slow_C_loss  = csite%today_A_decomp(ipa) * K3 * csite%slow_soil_C(ipa)
            ! NOTE: r_stsc = 1 so slow C and N loss is 0

            !------------------------------------------------------------------------------! 
            !Nitrogen deposition.                                                          !
            !Deposition rate is based on 0.0005kg N/m2/yr (Hedin et al.2009) and will be   ! 
            !added to MSN pool daily, independent of weather.                              !
            !------------------------------------------------------------------------------!
            N_deposition = 0.0009/365 !(kg N/m2/day )
            ! Add an extra 0.0004 kg N/m2/yr to represent asymbiotic BNF
            !N_deposition = 0.0009/365 + 0.0004/365  !(kg N/m2/day )
   
            !------------------------------------------------------------------------------! 
  
            !------------------------------------------------------------------------------! 
            !Hydrologic and Gas Nitrogen loss                                              !
            !Nitrogen is lost from mineralized soil N (MSN) pool when soils are at field   ! 
            !capacity.Only a fraction of the MSN is lost on a given day                    !
            !------------------------------------------------------------------------------!

            !----INITIALIZATION------------------------------------------------------------!
             soil_depth = 1.5        ! (Jackson 1996), most fine roots within top  1.5 m of  
                                     ! soil in tropical forests 
             water_density = 1000    !1000kg water/m3 water
             sec_in_day = 86400      ! (sec/day conversion)
             efficiency_factor= 0.05 ! fraction of N that can be lost by leaching in a day

             Avg_soil_water = sum(csite%soil_water(6:nzg,ipa))/9 !9 soil layers

             Nconc = 0.0 
             total_soil_water = 0.0
             total_soil_water_capacity =  0.0
             total_field_capacity = 0.0

             N_loss_emission = 0.0
             N_loss_total = 0.0
             DON_leached = 0.0
             csite%N_leached(ipa) =    0.0
             csite%N_gas_loss(ipa)  =    0.0

            !------------------------------------------------------------------------------!
 
            !---CALCULATIONS---------------------------------------------------------------!
  
            !1) Amount of mineralized soil nitrogen available for loss through leaching or
               !gas emissions [kg N/m2]      
            MSN_after_plant_uptake =  csite%mineralized_soil_N(ipa)                        &
                                   - csite%total_plant_nitrogen_uptake(ipa)         

            !2) Calculate [N] in soil, kg N/kg water          
            Nconc =  MSN_after_plant_uptake*1/soil_depth                                   &
                  * (1/ Avg_soil_water)*1/water_density  
            !UNITS: kg N/m2  *1 /m * m3soil/m3water *m3 water/kg water = kg N/kg water
       
            !3) Calulate the water input to the soil system  (Precip - evapotranspiration)  
                !(kg water/m2/day)         
               
                ! avg_transp  ! Average transpiration of water vapor [kg/m2/s] 
                              ! (the values look like it is kg/m2 not per sec)   
                ! avg_evap    !Average evaporation from leaf and ground surfaces,CAS [kg/m2/s] 

             water_input = cpoly%daily_pcpg(1) * sec_in_day  - cpoly%daily_transp(1) * sec_in_day      &
                         - cpoly%daily_evap(1) * sec_in_day  
       
            !4) Determine if soil moisture has exceeded field capacity
               ! if soil water exceed field capacity, allow N to leave the system
               !slmsts     ! Soil moisture at saturation[m3/m3]

               !layerloop:do k = cpoly%lsl(isi),nzg
               !loop through soil layers can change to (do k = 2,nzg) indicates to start
               ! from the second bottom layer 
     
               layerloop: do k = 6,nzg  !start at 1.5m depth
               !look up soil type (switched to using SITE level soils [mcd 9/30/08]
               nsoil = cpoly%ntext_soil(k,isi)  
               total_soil_water = total_soil_water + csite%soil_water(k,ipa) * dslz(k)  
               ! weighted for layers
               total_soil_water_capacity = total_soil_water_capacity + soil(nsoil)%slmsts  &
                                         * dslz(k) 
               !weighted for layers   
               total_field_capacity = total_field_capacity + soil(nsoil)%sfldcap * dslz(k)
               !calculate fraction of moisture capacity
               !fracw = csite%soil_water(k,ipa) / soil(nsoil)%slmsts 
               end do layerloop
       
            !5) If soil water has exceeded field capacity and it has rained,
               !calculate max N leached [kg N/m2/day]
               !make sure you don't leach all of the available soil N
               ! [water input is in mm/day = kg h20/m2/day]       
               if (total_soil_water > total_field_capacity  .and. water_input > 0) then
               csite%N_leached(ipa) =  water_input* Nconc
               endif  
             
            ! Just in case, make sure N_leached is never a negative value
               if (csite%N_leached(ipa)<0)then
               csite%N_leached(ipa) = 0
               endif
            !6) Modify leaching by the efficiency factor
                csite%N_leached(ipa) = csite%N_leached(ipa)* efficiency_factor 
            !7) Calculate Gas loss
            !Gas loss is approximately 20% of total N loss (Houlton et al 2006 PNAS)
            ! Generally, less then 2700 mm rain/yr is ~ 20 % but greater is about 50% gas loss  
            !20% based on reported rainfall at BCI
    
             N_loss_total = csite%N_leached(ipa)/.80
             csite%N_gas_loss(ipa) = N_loss_total*.20

            !-----Dissolved Organic Nitrogen loss JL. -------------------------------------!
            ! DON_loss = daily DON that is being set aside and building up until leaching 
            ! occurs. DON_loss_total is the build up of DON in between leaching intervals
            ! DON_leached = DON lost from MSN when soil water has exceeded field capacity
            ! and it has rained


            !DON loss is 1/5 of the total N deposition
             DON_loss    = 0.0001/365 !JL
 
             if (structural_N_loss > DON_loss) then
               structural_N_loss  = structural_N_loss  - DON_loss
               !print*, 'structural_N_loss > DON_loss'

             elseif (structural_N_loss <= DON_loss) then
               DON_loss = structural_N_loss
               structural_N_loss = 0
               !print*, 'structural_N_loss <= DON_loss'
             endif

             csite%interval_DON_loss(ipa) =  csite%interval_DON_loss(ipa) + DON_loss
             !print*, 'csite%interval_DON_loss(ipa)',csite%interval_DON_loss(ipa)

             if (total_soil_water > total_field_capacity  .and. water_input > 0) then
                DON_leached =  csite%interval_DON_loss(ipa)
                csite%interval_DON_loss(ipa) = 0.0
                !print*,'LEACHING CONDITIONS'
             endif  

 
!print*, 'total_don_loss' ,cgrid%total_DON_loss(1),'  cgrid%total_DIN_loss(1)',  cgrid%total_DIN_loss(1), &
!' cgrid%total_Ngas_loss(1)', cgrid%total_Ngas_loss(1)
               cgrid%total_DON_loss(1)   = cgrid%total_DON_loss(1)                         &
                                         + DON_leached* csite%area(ipa)
               cgrid%total_DIN_loss(1)   = cgrid%total_DIN_loss(1)                         &
                                         + csite%N_leached(ipa)* csite%area(ipa)
               cgrid%total_Ngas_loss(1)  = cgrid%total_Ngas_loss(1)                        &
                                         +  csite%N_gas_loss(ipa)* csite%area(ipa)
!print*, 'total_don_loss' ,cgrid%total_DON_loss(1),' DON_leached', DON_leached,'csite%area(ipa)', csite%area(ipa)
            !-------------------------------------------------------------------------------!
            !----- Mineralized pool. -------------------------------------------------------!

            csite%mineralized_N_input(ipa) = fast_N_loss + slow_C_loss / c2n_slow           & 
                                            +  structural_N_loss                            & 
                                            + N_deposition                               

            csite%mineralized_N_loss(ipa)  = csite%total_plant_nitrogen_uptake(ipa)         &
                                           + csite%N_leached(ipa) + csite%N_gas_loss(ipa)   &
                                          -    DON_leached                                  !JL
           
           !print*,'csite%mineralized_N_loss(ipa)',csite%mineralized_N_loss(ipa),           &
           !'csite%total_plant_nitrogen_uptake(ipa)',csite%total_plant_nitrogen_uptake(ipa),&
           !'csite%N_leached(ipa)',csite%N_leached(ipa),'gas loss', csite%N_gas_loss(ipa),  &
           !'DON_leached',DON_leached


           !JL Turn off immoblilization
                                    !   + csite%today_Af_decomp(ipa) * Lc * K1             & 
                                    !  * csite%structural_soil_C(ipa)                      
                                    ! * ( (1.0 - r_stsc) / c2n_slow - 1.0 / c2n_structural)  
 
            !JL added N from respired structural Pool and N deposition into MSN pool. 
            !JL added N_leached, N gas loss, and DON_loss to MSN loss calculation
            !Slow_C_loss = 0           
 
           !------------------------------------------------------------------------------!
            !JL: NOTES ON Immobilization equation above!

            !The First 2 lines calculate the Structural C loss. A fraction of that loss is 
            !respired through heterotrophic respiration (r_stsc), and a fraction goes into 
            !the slow soil !carbon pool (1.0 - r_stsc). Since the slow soil carbon pool has a
            !lower c2n raio(10) then the structural soil carbon pool(150), Nitrogen gets taken
            !out of the mineralized soil N pool to accomodate the carbon input to 
            !ssc(and acts as immobilization). The N needed for the structural soil pool is
            ! (1.0 - r_stsc) / c2n_slow, and the a!mount of N that is released from the structural 
            !soil C during decomposition is 1.0 / c2n_structural. The remaining amount is what is
            ! taken from mineralized soil N and it is therefore recorded as a mineralized soil N loss 

            !I set r_stsc = 1 to make the structural pool decomposition all goes to hetrotrophic
            !respiration, therefore, there is no need to have mineralized_soil_N loss associated with
            ! immobilization. Note, r_stsc is only used in this file for calculations and it should 
            !not impact any other calculations..
            !------------------------------------------------------------------------------!

            !------------------------------------------------------------------------------!
            !      All carbon fluxes have units kgC/m2/day, and we are updating on the     !
            ! daily time step.                                                             !
            !------------------------------------------------------------------------------!
            csite%fast_soil_C(ipa)       = csite%fast_soil_C(ipa) + csite%fsc_in(ipa)      &
                                         - fast_C_loss
            csite%structural_soil_C(ipa) = csite%structural_soil_C(ipa)                    &
                                         + csite%ssc_in(ipa) - structural_C_loss
            csite%structural_soil_L(ipa) = csite%structural_soil_L(ipa)                    &
                                         + csite%ssl_in(ipa) - structural_L_loss
            csite%slow_soil_C(ipa)       = csite%slow_soil_C(ipa) + slow_C_input           &
                                         - slow_C_loss
            
            !------------------------------------------------------------------------------!
            !      All nitrogen fluxes have units kgN/m2/day, and we are updating on the   !
            ! daily time step.                                                             !
            !------------------------------------------------------------------------------!
            csite%fast_soil_N(ipa)        = csite%fast_soil_N(ipa) + csite%fsn_in(ipa)     &
                                          - fast_N_loss
            csite%mineralized_soil_N(ipa) = csite%mineralized_soil_N(ipa)                  &
                                          + csite%mineralized_N_input(ipa)                 &
                                          - csite%mineralized_N_loss(ipa)

            !------------------------------------------------------------------------------!
            !      Force all pools to be either zero or positive.                          !
            !------------------------------------------------------------------------------!
            csite%fast_soil_C(ipa)        = max(0.0,csite%fast_soil_C(ipa))
            csite%structural_soil_C(ipa)  = max(0.0,csite%structural_soil_C(ipa))
            csite%structural_soil_L(ipa)  = max(0.0,csite%structural_soil_L(ipa))
            csite%slow_soil_C(ipa)        = max(0.0,csite%slow_soil_C(ipa))
            csite%fast_soil_N(ipa)        = max(0.0,csite%fast_soil_N(ipa))
            csite%mineralized_soil_N(ipa) = max(0.0,csite%mineralized_soil_N(ipa))
  
         end do patchloop
      end do siteloop
   end do polygonloop
   return
end subroutine update_C_and_N_pools
!==========================================================================================!
!==========================================================================================!
