  !! ----------------------------------------------------------------------------
  !! ----------------------------------------------------------------------------
  !!  *PROGRAM* crsim_subrs_airborne
  !!  @version crsim 4.0
  !!
  !!  *SRC_FILE*
  !!  crsim/src/crsim_subrs_airborne.f90
  !!
  !!
  !!
  !!  *LAST CHANGES*
  !!
  !! Oct 16 2023  - M.O.      surface backscatter
  !!
  !!  *DESCRIPTION* 
  !!
  !!  This program  contains the subroutines needed for computation of forward
  !!  airborne radar variables
  !!
  !!
  !!  Copyright (C) 2024 Stony Brook University - Brookhaven National Laboratory
  !!  Radar Science Group (//radarscience.weebly.com//)
  !!  Contact email address: mariko.oue@stonybrook.edu pavlos.kollias@stonybrook.edu
  !!
  !!  This source code is free software: you can redistribute it and/or modify
  !!  it under the terms of the GNU General Public License as published by
  !!  the Free Software Foundation, either version 3 of the License, or
  !!  any later version.
  !!
  !!  This program is distributed in the hope that it will be useful,
  !!  but WITHOUT ANY WARRANTY; without even the implied warranty of
  !!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  !!  GNU General Public License for more details.
  !!
  !!  You should have received a copy of the GNU General Public License
  !!  along with this program.  If not, see <http://www.gnu.org/licenses/>.
  !!
  !!----------------------------------------------------------------------------------------------------------------- 
  !!-----------------------------------------------------------------------------------------------------------------
  !!
  subroutine surface_backscatter_ocean(Freq,elev,pulse_len,uu,temp,Zsfc)
  Use crsim_mod
  Use phys_param_mod, only: c,K2,d2r,pi,m999
  Implicit None
  !
!  Type(conf_var),Intent(in)          :: conf
  Real*8, Intent(In)                 :: Freq ! radar frequency [GHz] (user input)
  Real*8, Intent(In)                 :: elev ! elevation angle
  Real*8, Intent(In)                 :: pulse_len ! pulse_length [m] (user input)
  Real*8, Intent(In)                 :: uu ![m/s] surface wind speed
  Real*8, Intent(In)                 :: temp ! [C] surface temperature (SST) 
  !--
  Real*8, Intent(Out)                :: Zsfc ! dBZ
  !--
  !-- salinity. default value = 38. Difference in Zsfc between sal=0 and 40 is 0.03 dB.
  Real*8, parameter                 :: sal=38 ! salinity User can change the value here
  !
  real*8                            :: thdeg ! incident angle [deg] 90-elev 
  real*8                            :: sigma0,s0dB ! surface backscatter sigma0
  real*8                            :: s2,lambda,KK2,frequency,th
  real*8                            :: w0,w1,F2
  !
  real*8                            :: w,tau0,a,b,tau
  real*8                            :: epsT,epsS,eps0,epsinf
  real*8                            :: alpha,beta,delta
  real*8                            :: sigma,sigma25
  complex*8,parameter               :: jj=(0.d0,1.d0)
  complex*8                         :: eps00

    ! for each grid point
    !------------------------------------------------------------------------------------
    !
  if(elev > -15.d0) Goto 999

  if(temp < -990.d0) Goto 999

  ! NRCS using a scattering model 'wu90'

  thdeg = 90 - abs(elev)
  th = thdeg * d2r

  frequency = Freq * 1.d9;
  lambda = c/frequency;
  
  s2 = m999; !for uu>20
  if(uu<7.0) then
     w0 = 0.009;
     w1 = 0.0276;       
     s2 = w0 + w1 * log10(uu);
  endif
  if((uu>7.0) .and. (uu<=20.0)) then
     w0 = -0.084;
     w1 = 0.138;
     s2 = w0 + w1 * log10(uu);
  endif

  if(s2 < 0.003) s2 = 0.003;
  
  ! Calculation of F2 (gamma fanction) using Klein method
  w = frequency * 2 *pi;

  tau0 = 1.768d-11 - 6.086d-13*temp + 1.104d-14*(temp**2.d0) - 8.111d-17*(temp**3.d0);
  b = 1 + 2.282d-5 * sal * temp - 7.638d-4*sal - 7.760d-6*(sal**2.d0) + 1.105d-8*(sal**3.d0);
  tau = tau0 * b;
  epsT = 87.134 - 1.949d-1*temp - 1.276d-2*(temp**2.d0) + 2.491d-4*(temp**3.d0);
  a = 1 + 1.613e-5*sal*temp - 3.656d-3*sal + 3.210d-5*(sal**2.d0) - 4.232d-7*(sal**3.d0);
  epsS = epsT * a;
  delta = 25 - temp;
  beta = 2.033d-2 + 1.266d-4*delta + 2.464d-6*(delta**2.d0) - sal*(1.849d-5 - 2.551d-7*delta + 2.551d-8*(delta**2.d0));
  sigma25 = sal*(0.182521 - 1.46192d-3*sal + 2.09324d-5*(sal**2.d0) - 1.28205d-7*(sal**3.d0));
  sigma = sigma25 * exp(-delta* beta);
  alpha = 0;
  epsinf = 4.9;
  eps0 = 8.854d-12;
  eps00 = epsinf + (epsS - epsinf)/(1 + ((jj*w*tau)**(1-alpha))) - jj * sigma / (w*eps0);
  F2 = abs((1 - sqrt(eps00))/(1 + sqrt(eps00)))**2.d0;

  !--
  sigma0 = F2 / (s2 * (cos(th)**4.d0)) * exp(-(tan(th)**2.d0)/s2); 

  s0dB = 10*log10(sigma0);
  !write(*,*) sigma0,F2,s2,thdeg,cos(th),tan(th)

  !--Convert surface cross section to radar reflectivity factor
  KK2 = K2;! 0.75 for W, 0.92 for all others
  if(Freq>90.0) then
     KK2 = 0.75;
  endif

  Zsfc = (10.d0**(s0dB/10))*(lambda**4.d0)*(10.d0**18.d0)/(pi**5*KK2*pulse_len);
  Zsfc = 10.d0*log10(Zsfc)

999 if((elev > -15.d0) .or. (temp < -990.d0)) then
     Zsfc = m999
  endif

  return
end subroutine surface_backscatter_ocean

!-------------------------------------------------------------------------------
! subroutine surface backscatter over land
!
subroutine surface_backscatter_land(Freq,elev,pulse_len, Zland_coarse,Zland_flat)
  Use crsim_mod
  Use phys_param_mod, only: c,K2,m999,pi
  Implicit None
  Real*8, Intent(In)                 :: Freq ! radar frequency [GHz] (user input) 
  Real*8, Intent(In)                 :: elev ! elevation angle                 
  Real*8, Intent(In)                 :: pulse_len ! pulse_length [m] (user input)
  !--
  Real*8, Intent(Out)                :: Zland_coarse,Zland_flat ! dBZ  
  !--
  real*8                            :: thdeg ! incident angle [deg] 90-elev 
  real*8                            :: sigma0_coarse_dB,sigma0_flat_dB ! surface backscatter sigma0 
  real*8                            :: lambda,KK2,frequency
  !     

  if(elev > 0.d0) Goto 999

  thdeg = 90 - abs(elev)

  frequency = Freq * 1.d9;
  lambda = c/frequency;

  ! Empilical equations for sigma0 as a function of incidence angle
  ! for coarse surface and flat surface
  !For coarse (i.e. wet snow), 
  sigma0_coarse_dB = 6.767*exp(-0.293*thdeg) - 0.127*thdeg - 6.767
  !For flat (i.e. flat road), 
  sigma0_flat_dB = 6.767*exp(-0.293*thdeg) - 0.127*thdeg + 8.233

  !--Convert surface cross section to radar reflectivity factor
  KK2 = K2;! 0.75 for W, 0.92 for all others
  if(Freq>90.0) then
     KK2 = 0.75;
  endif

  Zland_coarse = (10.d0**(sigma0_coarse_dB/10))*(lambda**4.d0)*(10.d0**18.d0)/(pi**5*KK2*pulse_len); 
  Zland_flat = (10.d0**(sigma0_flat_dB/10))*(lambda**4.d0)*(10.d0**18.d0)/(pi**5*KK2*pulse_len); 

  Zland_coarse = 10.d0*log10(Zland_coarse)
  Zland_flat   = 10.d0*log10(Zland_flat)

!  write(*,*) sigma0_coarse_dB,sigma0_flat_dB
999 if(elev > 0.d0)  then
     Zland_coarse = m999
     Zland_flat = m999
  endif

  return
end subroutine surface_backscatter_land


subroutine airborne_radial_velocity(conf,elev,azim,Dopp_airborne)
  Use crsim_mod
  Use phys_param_mod, only: d2r,m999
  Implicit None
  Type(conf_var),Intent(in)          :: conf
  Real*8, Intent(In)                 :: elev
  Real*8, Intent(In)                 :: azim
  Real*8, Intent(Out)                :: Dopp_airborne
  Real*8                             :: uu,vv,ww
  Real*8                             :: azth,elth

  !- airplane speed (u,v,w)
  azth = conf%airborne_azdeg * d2r
  elth = conf%airborne_eldeg * d2r
  ww = conf%airborne_spd * sin(elth)
  uu = conf%airborne_spd * cos(elth) * sin(azth)
  vv = conf%airborne_spd * cos(elth) * cos(azth)

  !- radial velocity at each grid box
  azth = azim * d2r
  elth = elev * d2r
  Dopp_airborne = (uu * cos(azth) + vv * sin(azth)) * cos(elth) + ww * sin(elth)
  Dopp_airborne = Dopp_airborne * (-1.d0)

end subroutine airborne_radial_velocity



  !
  !======================================================!
  !======================================================!
  !=== END Subroutines for airborne radar simulations ===!
  !======================================================!
  !======================================================!
     
