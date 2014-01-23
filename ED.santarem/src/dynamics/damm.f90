!  This is based off of Davidson et al. (2012) Global Change Biology Paper
!-------------------------------------------------------------------------

Module damm_model
  implicit none

  real, parameter :: enzyme_capacity=5.38e10 / 3600.
  !Davidson refers to this as [alpha_Sx].  Davidson's value is:
  !   5.38e10 (mg C) / (cm3 soil) / (hour)
  !      Divide by 1e6 to get kg C.
  !      Multiply by 1e6 to get m3 soil.
  !      Divide by 3600 to get seconds.

  real, parameter :: activation_energy = 72.26 * 1000. / 8.314
  !Davidson refers to this as [Easx] in kJ / mol.
  !   Here, we multiply by 1000 to get J/mol.
  !   Then, divide by gas constant (8.314 J/mol/K) to get K.

  real, parameter :: Michaelis_constant_substrate = 9.95e-7 * 1.e3
  !This is the Michaelis constant for substrate supply.
  !  Davidson has units of (g C) / (cm3 soil).
  !  We divide by 1000 to get kg C.
  !  We multiply by 1.e6 to get m3 soil.
  !  Final units: (kgC)/(m3 soil)

  real, parameter :: soluble_fraction=4.14e-4 ! Dimensionless
  !This is "p" in Davidson et al. 2012.

  real, parameter :: Michaelis_constant_O2 = 0.121
  !This is the Michaelis constant for O2.
  !  Units are m3/m3.

Contains
  
  subroutine damm_reaction_velocity(temperature, total_substrate,  &
       saturated_water_content, actual_water_content, reaction_velocity)
    implicit none
    
    real, intent(in) :: temperature  ! K
    real, intent(in) :: total_substrate ! kgC/m3
    real, intent(in) :: saturated_water_content ! m3/m3
    real, intent(in) :: actual_water_content ! m3/m3
    real, intent(out) :: reaction_velocity ! kgC/m3/second

    real :: Vmax  ! kgC/m3/s
    real :: soluble_substrate ! kgC/m3; [Sxsoluble] in Davidson
    real :: diffusion_coeff_substrate ! m3/m3; Dliq in Davidson
    real :: substrate_at_enzyme ! kgC/m3; [Sx] in Davidson
    real :: michaelis_menton_substrate  ! [-]
    real :: air_filled_porosity ! [m3/m3]; a in Davidson
    real :: diffusion_coeff_O2 ! [-]; Dgas in Davidson
    real :: O2_at_enzyme ! [O2] in Davidson
    real :: michaelis_menton_O2  ! [-]

    Vmax = enzyme_capacity * exp(-activation_energy / temperature)

    soluble_substrate = soluble_fraction * total_substrate
    diffusion_coeff_substrate = 1./saturated_water_content**3
    substrate_at_enzyme = soluble_substrate * diffusion_coeff_substrate *  &
         actual_water_content**3
    michaelis_menton_substrate = substrate_at_enzyme /   &
         (Michaelis_constant_substrate + substrate_at_enzyme)

    air_filled_porosity = saturated_water_content - actual_water_content
    diffusion_coeff_O2 = 1. / (saturated_water_content)**(4./3.)
    O2_at_enzyme = diffusion_coeff_O2 * 0.209 * air_filled_porosity**(4./3.)
    michaelis_menton_O2 = O2_at_enzyme /   &
         (Michaelis_constant_O2 + O2_at_enzyme)

    print*,Vmax

    reaction_velocity = Vmax * michaelis_menton_substrate *  &
         michaelis_menton_O2

    return
  end subroutine damm_reaction_velocity

end Module damm_model
