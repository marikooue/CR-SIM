!! ----------------------------------------------------------------------------
!! ----------------------------------------------------------------------------
!!  *PROGRAM* crsim_mod
!!  @version crsim 3.x
!!
!!  *SRC_FILE*
!!  crsim/src/crsim_mod.f90
!!
!!
!!  *LAST CHANGES*
!!
!!  Sep 15 2015  :: A.T.    included horizontal wind variables u and v in the structure env_var.
!!                          Subroutines for allocating,deallocating and nullifying env_var variables
!!                          are accordingly updated.
!!  Sep 17 2015  :: A.T     added a new variable azim in the rmout structure 
!!  Sep 23 2015  :: A.T     included new variables in the structure rmout for Doppler velocity and
!!                          Spectrum Width. Now those variables are computed for both slant and nadir 
!!                          view. Subroutines for allocating,deallocating and nullifying rmout variables
!!                          are accordingly updated.
!!  Jan 12 2016  :: A.T     added new configuration variable thr_mix_ratio (structure conf_var).
!!                          The input mixing  ratio <= than thr_mix_ratio value is set to 0.
!!  Jan 13 2016  :: A.T     included xlat, xlong and tke (turbulence kinetic energy) in the structure
!!                          env_var. Subroutines for allocating,deallocating and nullifying env_var variables 
!!                           are accordingly updated.
!!  Jan 13 2016  :: A.T     added new configuration variables Theta1 and dr (beamwidth and range 
!!                          resolution respectivelly). Also added variables computed with Theta1 and dr,
!!                          sigma_theta and sigma_r, needed for spectrum width computations
!!  Jan 14 2016  :: A.T     included a new variable sw_t in the rmout structure (spectrum width due to 
!!                          the turbulence contribution only).Subroutines for allocating, deallocating and
!!                          nullifying rmout_var variables  are accordingly updated.
!!  Jan 14 2016  :: A.T     included new variables SWt and SWt_tot in the rmout structure (spectrum width due to 
!!                          all contributions computed).Subroutines for allocating, deallocating and nullifying
!!                          rmout_var variables  are accordingly updated.
!!  Jan 14 2016  :: A.T     included new wind shear variables Ku,Kv and Kw to the structure env_var. 
!!                          Subroutines for allocating, deallocating and nullifying env_var variables are
!!                          accordingly updated.
!!  Jan 15 2016  :: A.T     included new variables sw_s and sw_v in the rmout structure (spectrum width due to 
!!                          wind shear in radar volume and cross wind, respectively).Subroutines for allocating,
!!                          deallocating and nullifying  rmout_var variables  are accordingly updated.
!!  FEb 15 2016  :: A.T     included new variable "range" in the rmout structure. Subroutines for allocating, 
!!                          deallocating and nullifying  rmout_var variables are accordingly updated.
!!  FEb 15 2016  :: A.T     included new configuration parameter ZMIN, a coefficient in expression for radar
!!                          sensitivity limitation with range 
!!  FEb 17 2016  :: A.T     included new configuration parameter ceiloID for cloud lidar measurements
!!  FEb 17 2016  :: A.T     included new structure lout_var containing the cloud lidar output variables.
!!                          Also written are new subroutines for  allocating, deallocating and nullifying lout variables. 
!!  FEb 17 2016  :: A.T     included new variable diff_back_phase in the polarimetric output. Affected structures are
!!                          mout_var and mrad_var, as well as all the subroutines for allocating, deallocating and
!!                          nullifying this structures.
!!  Apr 29 2016  :: A.T     Included new configuration parameter mplID for ice cloud vertical measurements
!!  Apr 29 2016  :: A.T     Included new structure mpl_var containing the mpl lidar output variables.
!!                          Also written are new subroutines for  allocating, deallocating and nullifying mpl variables.
!!  May 06 2016  :: A.T.    Included new configuration parameters aeroID, aero_tau and aero_lidar_ratio for optional aerosol
!!                          profile in MPL measurements
!!  May 09 2016  :: A.T.    Included new structure aero_var containing the mpl aerosol lidar output variables
!!  MAY 12 2016  :: A.T.    Introduced  a new MPL variable,rayleigh_back, for the molecular backscatter [m-1 sr-1]
!!  MAR 23 2017  :: M.O.    Introduced ARSCL and MWR LWP parameters
!!  JUL 21 2017  :: M.O     Included to read RAMS outout 
!!  OCT 11 2017  :: K. YU   Added variable "pi"
!!  OCT 30 2017  :: D.W.    Added hydro50_var, allocate_hydro50_var, initialize_hydro50_var and deallocate_hydro50_var for P3 scheme (added by DW)
!!  Apr      2018  ::M.O    Incorporated SAM warm bin microphysics (MP_PHYSICS=70)
!!  Jun 17   2018  ::M.O    Incorporated SAM Morrison 2-moment microphysics (MP_PHYSICS=75)
!!  MAY  27  2018  ::M.O    Incorporated spectra generator for bulk moment microphysics (spectra_var)
!!                          Added suproutines of processing_ds, compute_doppler_spectra, spectra_generator,
!!                          spectra_unfolding,Turbulence_Convolution,interp1,hsmethod.
!!                          Modified processing*, GetPolarimetricInfofromLUT*, 
!!                          compute_polarim_vars, compute_rad_forward_vars
!!
!!  
!!  Nov 2018       :: A.T.  Introduced parameter maxNumOfHydroTypes equal to the maximal number of Hydrometeor Types
!!                          allowed (currrently 8 for RAMS, MP=40). Corrected 
!!                          dimension for those variables whose dimension was not updated (i.e. was equal to less than maxNumOfHydroTypes)
!!  Aug 2019       :: A.T.  Some numerical issue corrected (changes in value in conversion from real to integer and similar) 
!!  Dec 2019       :: A.T   Included new variable RHOhvc in the mrad_var structure. Subroutines for allocating, 
!!                          deallocating and nullifying  mrad_var variables are accordingly updated.
!!  Dec 2019       :: A.T   Included new variables RHOhv and RHOhv_tot in the mout_var structure. Subroutines for allocating, 
!!                          deallocating and nullifying  mout_var variables are accordingly updated.
!!  Oct 2023       :: M.O.  Included airborne radar simulation. airborne_var. env_var and conf_var are also modified to
!!                          include airborne-simulation related variables.

!!  *DESCRIPTION* 
!!
!!  This module contains a number of subroutines for allocating, deallocating  and
!!  nullifying (or setting to the missed value) data types used in crsim. Those
!!  data types are also defined in this module.
!!
!!
!!  Copyright (C) 2014 McGill Clouds Research Group (//radarscience.weebly.com//)
!!  Contact email address: aleksandra.tatarevic@mcgill.ca pavlos.kollias@mcgill.ca
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
module crsim_mod
Use phys_param_mod, ONLY: maxNumOfHydroTypes
Implicit None


#ifdef __PRELOAD_LUT__


#ifndef __LUT_FILE_LIST_SIZE__
#define __LUT_FILE_LIST_SIZE__ 10240
#endif

Type lut_file_type

  Character(len=1024) :: filenameLUTS

  real*8                      :: llam! mm
  real*8                      :: ltemp ! K
  real*8                      :: lref_re,lref_im
  integer                     :: lnd,lnelev

  real*8,Dimension(:),Allocatable              :: elev
  real*8,Dimension(:),Allocatable              :: ldiam ! mm
  real*8,Dimension(:),Allocatable              :: aoverb
  complex*16,Dimension(:,:,:,:),Allocatable    :: sb ! mm
  complex*16,Dimension(:,:),Allocatable        :: sf11,sf22 ! mm

End Type



Type(lut_file_type), Dimension(__LUT_FILE_LIST_SIZE__) :: lut_files
Integer*4                           :: size_lut_files
Integer*4                           :: lut_search_flag




Type mpl_lut_file_type

  Character(len=1024) :: filenameLUTS

  real                       :: tempK,density,wavelength !  K,gr/cm^3 and um
  real                       :: ni_r,ni_c ! real and imaginary part of refr. index
  integer                    :: ndiam ! number of diams stored

  real,  Dimension(:),Allocatable   :: ldiam, qext, qback ! um, m^2,m^2

End Type

Type(mpl_lut_file_type), Dimension(maxNumOfHydroTypes) :: mpl_lut_files 
!Type(mpl_lut_file_type), Dimension(8) :: mpl_lut_files
!
Type lut_ceilo_file_type

  real                       :: tempK,wavelength !  K and um
  real                       :: ni_r,ni_c ! real and imaginary part of refr. index
  integer                    :: ndiam ! number of diams stored

  real,  Dimension(:),Allocatable   :: ldiam, qext, qback ! um, m^2,m^2

End Type

Type(lut_ceilo_file_type)   :: lut_ceilo_file


#endif






Type conf_var

Integer            ::  OMPThreadNum                      ! The number of OpenMP threads used in the code. This number is used for OMP_SET_NUM_THREADS().
Character(len=365) ::  WRFInputFile                      ! the name of the WRF netcdf input file 
Character(len=365) ::  WRFmpInputFile                    ! the name of the spectral WRF netcdf input file with explicit microphysics 
                                                         ! for MP_PHYSICS=20; ="NULL" for MP_PHYSICS=10; RAMS header text file for MP_PHYSICS=40 
Character(len=365) ::  WRF_MP_PHYSICS_20_InputFiles(4)   ! names of the additional 4 ASCII files with diameters,
                                                         ! densities, masses and fall velocities for MP_PHYSICS = 20, =NULL otherwise

Integer            ::  ix_start,ix_end      ! the starting and ending indices of domain in x-direction (E-W)
Integer            ::  iy_start,iy_end      ! the starting and ending indices of domain in y-direction (S-N)
Integer            ::  iz_start,iz_end      ! the starting and ending indices of domain in z-direction (vertical)
Integer            ::  it                   ! time step of the WRF simoulation 

Integer            ::  MP_PHYSICS           ! =10 or =20  (WRF microphysics Morrison 10, JFan 20), 9 for Milbrandt&Yau with snow spherical,
                                            !                                                     901 for Milbrandt&Yau with snow not spherical
Integer            ::  MP10HailOption       ! Morrison (MP_PHYSICS=10) hail option instead of graupel (graupel=0, hail=1)
Integer            ::  snow_spherical       ! =1 for Milbrandt&Yau with snow spherical, any other int for snow not spherical
Character(len=365) ::  OutFile              ! the name of the main output file
Integer            ::  ixc,iyc              ! horizontal indices of the radar position
Real*8             ::  zc                   ! [m]  height at whivh radar is positioned  

Character(len=365) ::  path_back_scatt_libs ! the path of LUTS for backscattering
Character(len=365) ::  file_ice_parameter   ! the precalculated of ice PSD parameters added by DW 2017/10/30 for P3
Character(len=20)  ::  hydro_luts(maxNumOfHydroTypes)        ! names of LUTS backscattering species 1-cloud, 2-rain, 3-ice, 4-snow 5-graupel, 6-hail ; 7-drizzle (RAMS) ; 8-aggregates (RAMS), 8 is MAX NUM
!Character(len=20)  ::  hydro_luts(8)
Real*8             ::  thr_mix_ratio(maxNumOfHydroTypes)     ! threshold values for the input mixing ratio. Input mixing ratios lower or equal than the specified threshold are set to zero.
Integer            ::  horientID(maxNumOfHydroTypes)         ! distribution of particle orientation =1 for fully chaotic orientation,=2 for random orientation in
                                           ! horizontal plane and=3 for two-dimensional axisymmetric Gaussian distribution of
                                           ! orientation with zero mean and sigma standard deviation. 
Real*8             ::  sigma(maxNumOfHydroTypes)             ! [degree] -The width in degrees for horientID==3. If horientID/=3, the values of sigma are read but disregarded.
Real*8             ::  freq                 ! [GHz] frequency 
Real*8             ::  elev                 ! [degrees] elevation, fixed if positive, and if negative then 
                                            !           elevation of each pixel (ix,iy,iz) is relative to the radar origin 
                                            !           given by indices (ixc,iyc)
Integer            :: radID                 ! =1 to turn off the computation of polarimetric variables, the output would consist of reflectivity, doppler, spect. width and attenuation, /=1 the full polarimetric output is written 
Integer            :: nht                   ! number of hydrometeor classes in the WRF microphysical package, =5 for  MP_PHYSICS=20,10 and =6 for  MP_PHYSICS=9,30, 8 for MP_PHYSICS=40
!
Real*8             ::  Theta1               ! [deg] radar beamwidth i.e.  one-way angular resolution 
Real*8             ::  dr                   ! [m] radar range resolution
!
Real*8             :: sigma_theta           ! sigma_theta= (Theta1) / (4 sqrt(ln2))
Real*8             :: sigma_r               ! sigma_r=0.35 dr
Real*8             :: ZMIN                  ! coefficient in expression for radar sensitivity limitation with range DBZ_MIN[dBZ]=ZMIN + 20 log10 (range [km]) 
Integer            :: ceiloID               ! whether or not to include cloud lidar (celiometer) measurements =1 yes, 0 otherwise
Integer            :: mplID                 ! whether or not to include micropulse lidar measurements (=1 fr=353 nm, =2 fr=532 nm). If mplID<=0 do not include mpl simmulations
Real*8             :: mpl_wavel             ! wavelength of micropulse lidar (353 nm for mplID=1 or 532 for mplID=2). 
Integer            :: aeroID                ! Whether (aeroID=1)or not(aeroID=0) to include an average aerosol profile in MPL simulation
Real*8             :: aero_tau              ! Value of desired  path integrated extinction to normalize the aerosol profile, if negative, do not normalize
Real*8             :: aero_lidar_ratio      ! fixed value for aerosol lidar ratio
Integer            :: spectraID             ! =1 to simulate Doppler spectra; /=1 to not include Doppler spectra
Integer            :: arsclID               ! Whether (arsclID=1)or not(arsclID=0) to include ARSCL products after radar&lidar simulations (needs mplID>0 and ceiloID=1)
Integer            :: mwrID                 ! Whether (mwrID=1)or not(mwrID=0) to include MWR LWP estimations after radar&lidar simulations 
Real*8             :: mwr_view              ! MWR field of view in degree used when mwrID=1
Real*8             :: mwr_alt               ! Altitude [m] of MWR location used when mwrID=1

Integer            :: InpProfile_flag       ! Check if use Input profile data file (third file of the second argument) used for SAM (MP_PHYSICS=70)
Character(len=365) :: InputProfile          ! the name of the profiling data file for SAM (netcdf) 
!
! airborne added by oue Oct 2023
Integer            :: airborne              ! ==0 for radar looking up  (no airborne) and /= 0 otherwise
Real*8             :: pulse_len             ! pulse length for airborne simulation
Real*8             :: airborne_spd          ! speed of air craft m/s
Real*8             :: airborne_azdeg        ! azimuth direction of moving aircraft deg from north
Real*8             :: airborne_eldeg        ! elevation angle of moving aircraft deg from horizon

!
End Type conf_var
!!
!!
Type env_var
!
Integer   :: nx
Integer   :: ny
Integer   :: nz
!
Real*8                                    :: dx     ! [m] x-resolution 
Real*8                                    :: dy     ! [m] y-resolution 
Real*8,Dimension(:),Allocatable           :: x      ! [m] extent of the scene in E-W direction
Real*8,Dimension(:),Allocatable           :: y      ! [m] extent of the scene in S-N direction
Real*8,Dimension(:,:),Allocatable         :: xlat   ! [deg]  latitude, north is positive
Real*8,Dimension(:,:),Allocatable         :: xlong  ! [deg]  longitude, west is positive
Real*8,Dimension(:,:,:),Allocatable       :: z      ! [m] geopotential height
Real*8,Dimension(:,:,:),Allocatable       :: u      ! [ms^-1]  horizontal u wind component
Real*8,Dimension(:,:,:),Allocatable       :: v      ! [ms^-1]  horizontal  v wind component
Real*8,Dimension(:,:,:),Allocatable       :: w      ! [ms^-1]  vertical w wind component, positive up (+z)
Real*8,Dimension(:,:,:),Allocatable       :: Ku     ! [s^-1]  wind shear along x (for u wind component)
Real*8,Dimension(:,:,:),Allocatable       :: Kv     ! [s^-1]  wind shear along y (for  v wind component)
Real*8,Dimension(:,:,:),Allocatable       :: Kw     ! [s^-1]  wind shear in vertical (for w wind component)

Real*8,Dimension(:,:,:),Allocatable       :: rho_d  ! [kg/m^3] dry air density
Real*8,Dimension(:,:,:),Allocatable       :: qvapor ! [kg/kg]  Water vapor mixing ratio
Real*8,Dimension(:,:,:),Allocatable       :: temp   ! [C]      temperature
Real*8,Dimension(:,:,:),Allocatable       :: press  ! [mb]     pressure 
Real*8,Dimension(:,:,:),Allocatable       :: tke    ! [m^2/s^2]  turbulence kinetic energy

Real*8,Dimension(:,:),Allocatable       :: sfctemp! [degC]  surface temperature (sea surface temperature over ocean)
Real*8,Dimension(:,:),Allocatable       :: sfcwspd! [m/s]  surface wind speed
!
End Type env_var
!!
Type hydro_var
!
integer                 :: nx
integer                 :: ny
integer                 :: nz
integer                 :: nht
!
real*8, Dimension(:,:,:,:),Allocatable    :: qhydro  ! hydrometeor mixing ratio [kg/kg] 1-cloud, 2-rain, 3-ice, 4-snow 5-graupel
real*8, Dimension(:,:,:,:),Allocatable    :: qnhydro ! hydro number concentrat. [1/kg] 
!
End Type hydro_var
!!-------------------------------
! Type hydro50_var is added by DW 2017/10/30 for P3
Type hydro50_var
!
integer                 :: nx
integer                 :: ny
integer                 :: nz
integer                 :: nht
!
real*8, Dimension(:,:,:,:),Allocatable    :: qhydro  ! hydrometeor mixing ratio [kg/kg] 1-cloud, 2-rain, 3-small ice, 4-unrimed ice, 5-graupel, 6-partially rimed ice
real*8, Dimension(:,:,:,:),Allocatable    :: qnhydro ! hydro number concentrat. [1/kg] 
real*8, Dimension(:,:,:),Allocatable      :: qir     ! Rimed ice mass mixing ratio. [kg/kg]
real*8, Dimension(:,:,:),Allocatable      :: qib     ! Rimed ice volume mixing ratio [m3/kg]
!
End Type hydro50_var
!!
Type hydro20_var
!
integer                 :: nx
integer                 :: ny
integer                 :: nz
integer                 :: nbins
!
! ff1,ff2,ff5 are bin number concentration in 1/m^3 multiplied by bin mass and
! divided by density of dry air
real*8,Dimension(:,:,:,:),Allocatable      :: ff1 ! [kg/kg] cloud/rain bin output [nx,ny,nz,nbins]
real*8,Dimension(:,:,:,:),Allocatable      :: ff5 ! [kg/kg] ice/snow bin output [nx,ny,nz,nbins]
real*8,Dimension(:,:,:,:),Allocatable      :: ff6 ! [kg/kg] graupel bin output [nx,ny,nz,nbins] 
!!
End Type hydro20_var
!!
Type hydro70_var
!
integer                 :: nx
integer                 :: ny
integer                 :: nz
integer                 :: nbins, nht
!
! fn1 is number concentration 1/m3 for liquid,fm1 is mass kg/kg for liquid
real*8,Dimension(:,:,:,:),Allocatable  :: fm1 ! [kg/kg] cloud/rain mass bin output [nx,ny,nz,nbins]
real*8,Dimension(:,:,:,:),Allocatable  :: fn1 ! [1/m3] cloud/rain number bin output [nx,ny,nz,nbins]
real*8,Dimension(:,:,:,:),Allocatable  :: qhydro  ! hydrometeor mixing ratio [kg/kg] 1-cloud, 2-rain,
real*8,Dimension(:,:,:,:),Allocatable  :: qnhydro ! hydro number concentrat. [1/m3] 
real*8,Dimension(:,:,:),Allocatable    :: naero   ! [1/m3] Aerosol Number concentration
!!
End Type hydro70_var
!!
!!
Type scatt_type_var
integer                                        :: nx
integer                                        :: ny
integer                                        :: nz
integer                                        :: nbins ! number of bins
integer                                        :: ihtf  ! the hydro category 1=cloud+rain, 2=ice+snow,3=graupel+hail
integer                                        :: id1,id2 ! starting and ending bin indices for each scatt. type
real*8,Dimension(:), Pointer :: diam  ! [m] diameters [nbins]
real*8,Dimension(:), Pointer :: rho   ! [kg/m^3] density   [nbins]
real*8,Dimension(:), Pointer :: mass  ! [kg] mass  [nbins]
real*8,Dimension(:), Pointer :: fvel  ! [m/s] bin fall velocity [nbins]
real*8,Dimension(:,:,:), Pointer :: qq    ! [kg/kg] mixing ratio [nx,ny,nz]
real*8,Dimension(:,:,:,:), Pointer :: N   ! [kg/kg] bin mass concentration [nx,ny,nz,nbins]
real*8,Dimension(:,:,:,:), Pointer :: NN  ! [1/m^3] bin number concentration [nx,ny,nz,nbins] for SAM
real*8,Dimension(:), Pointer :: diam1  ! [m] diameters [nbins+1] boundary, for SAM
real*8,Dimension(:), Pointer :: mass1  ! [/kg] mass [nbins+1] boundary, for SAM
real*8,Dimension(:,:,:,:), Pointer :: diam2 ! [m] diameters [nbins*nx*ny*nz] changing wigh N&NN
real*8,Dimension(:,:,:,:), Pointer :: rho2  ! [kg/m^3] density [nbins*nx*ny*nz] changing wigh N&NN
real*8,Dimension(:,:,:,:), Pointer :: fvel2 ! [m/s] bin fall velocity [nbins*nx*ny*nz] changing wigh N&NN
End Type scatt_type_var 
!!
!!
Type mout_var
Integer   :: nx
Integer   :: ny
Integer   :: nz
Integer   :: nht
!
Real*8,Dimension(:,:,:,:),Allocatable :: Zvv   ! [mm^6/m^3] Reflectivity at vv polarization per hydrometeor
Real*8,Dimension(:,:,:,:),Allocatable :: Zvh   ! [mm^6/m^3] Reflectivity at vh polarization per hydrometeor
Real*8,Dimension(:,:,:,:),Allocatable :: RHOhv ! [-] Cross-correlation coefficient  per hydrometeor
Real*8,Dimension(:,:,:,:),Allocatable :: Zdr   ! [-]        Differential reflectivity per hydrometeor
Real*8,Dimension(:,:,:,:),Allocatable :: LDRh  ! [-]        Linear Depolarization Ratio per hydrometeor
Real*8,Dimension(:,:,:,:),Allocatable :: Kdp   ! [deg/km]   Specific Differencial Phase per hydrometeor
Real*8,Dimension(:,:,:,:),Allocatable :: Adp   ! [dB/km]    Differencial Attenuation per hydrometeor
Real*8,Dimension(:,:,:,:),Allocatable :: Av    ! [dB/km]    Specific Vertical Attenuation per hydrometeor
Real*8,Dimension(:,:,:,:),Allocatable :: diff_back_phase ! [deg]  Differential Backscattered Phase per hydrometeor
!
Real*8,Dimension(:,:,:),Allocatable :: Zvv_tot   ! [mm^6/m^3] Total Reflectivity at vv polarization
Real*8,Dimension(:,:,:),Allocatable :: Zvh_tot   ! [mm^6/m^3] Total Reflectivity at vh polarization
Real*8,Dimension(:,:,:),Allocatable :: RHOhv_tot ! [-]        Total Cross-correaltion coefficient
Real*8,Dimension(:,:,:),Allocatable :: Zdr_tot   ! [-]        Total Differential reflectivity
Real*8,Dimension(:,:,:),Allocatable :: LDRh_tot  ! [-]        Total Linear Depolarization Ratio
!
Real*8,Dimension(:,:,:),Allocatable :: Kdp_tot   ! [deg/km]   Total Specific Differencial Phase
Real*8,Dimension(:,:,:),Allocatable :: Adp_tot   ! [dB/km]    Total Differencial Attenuation
Real*8,Dimension(:,:,:),Allocatable :: Av_tot    ! [dB/km]    Total Specific Vertical Attenuation
Real*8,Dimension(:,:,:),Allocatable :: diff_back_phase_tot ! [deg]  Differential Backscattered Phase
!
End type mout_var
!!
Type rmout_var
Integer   :: nx
Integer   :: ny
Integer   :: nz
Integer   :: nht
!
Real*8,Dimension(:,:,:,:),Allocatable :: Zhh     ! [mm^6/m^3] Reflectivity at hh polarization per hydrometeor
Real*8,Dimension(:,:,:,:),Allocatable :: Dopp    ! [m/s]      Radial Doppler Velocity per hydrometeor class
Real*8,Dimension(:,:,:,:),Allocatable :: dDVh    ! [m/s]      Spectrum width per hydrometeor class due to diff hydrometeor fall velocities at diff sizes
Real*8,Dimension(:,:,:,:),Allocatable :: SWt     ! [m/s]      Total (final) Spectrum width per hydrometeor class due to different contributions (turbulence, hydr. term. velocity ...)
Real*8,Dimension(:,:,:,:),Allocatable :: DVh     ! [m/s]      Reflectifity Weighted Velocity per hydrometeor (def:air vert veloc=0, elev=90 deg)
Real*8,Dimension(:,:,:,:),Allocatable :: Dopp90  ! [m/s]    Doppler Velocity per hydrometeor class at elev=90
Real*8,Dimension(:,:,:,:),Allocatable :: dDVh90  ! [m/s]    Spectrum width per hydrometeor class at elev=90
Real*8,Dimension(:,:,:,:),Allocatable :: Ah      ! [dB/km]    Specific Horizontal Attenuation per hydrometeor
!
Real*8,Dimension(:,:,:),Allocatable :: Zhh_tot    ! [mm^6/m^3] Total Reflectivity at hh polarization
Real*8,Dimension(:,:,:),Allocatable :: Dopp_tot   ! [m/s]      Total hydrom. Radial Doppler Velocity 
Real*8,Dimension(:,:,:),Allocatable :: dDVh_tot   ! [m/s]      Total hydrom. Spectrum width  due to diff hydrometeor fall velocities at diff sizes
Real*8,Dimension(:,:,:),Allocatable :: SWt_tot    ! [m/s]      Total hydrom. Spectrum width due to different contributions (turbulence, hydr. term. velocity ...)
Real*8,Dimension(:,:,:),Allocatable :: DVh_tot    ! [m/s]      Total hydrom. Reflectifity Weighted Velocity (def: w=0 m/s; elev=90deg)
Real*8,Dimension(:,:,:),Allocatable :: Dopp90_tot ! [m/s]    Total hydrom. Vertical Doppler Velocity (elev=90)
Real*8,Dimension(:,:,:),Allocatable :: dDVh90_tot ! [m/s]    Total hydrom. Spectrum width (elev=90)

Real*8,Dimension(:,:,:),Allocatable :: Ah_tot    ! [dB/km]    Total Specific Horizontal Attenuation
!
Real*8,Dimension(:,:,:),Allocatable :: elev      ! [degree] elevation
Real*8,Dimension(:,:,:),Allocatable :: azim      ! [degree] azimuth (convention E=0; W=180, N=90, and S=270 degrees)
Real*8,Dimension(:,:,:),Allocatable :: range     ! [m] range
Real*8,Dimension(:,:,:),Allocatable :: sw_t      ! [m/s]  spectrum width due to turbulence contribution only   
Real*8,Dimension(:,:,:),Allocatable :: sw_s      ! [m/s]  spectrum width due to wind shear in radar volume
Real*8,Dimension(:,:,:),Allocatable :: sw_v      ! [m/s]  spectrum width due to cross wind
  
!
End type rmout_var
!!
Type lout_var ! output lidar variables 4th dimension added by oue gor RAMS
Integer                             :: nx,ny,nz,nht 
Real*8,Dimension(:,:,:,:),Allocatable :: ceilo_back_true  ! lidar true backscatter [m-1 sr-1]
Real*8,Dimension(:,:,:,:),Allocatable :: ceilo_ext        ! lidar extinction coefficient [m-1 sr-1]]
Real*8,Dimension(:,:,:,:),Allocatable :: ceilo_back_obs   ! lidar observed backscatter (true back + attenuation) [m-1 sr-1]]
Real*8,Dimension(:,:,:,:),Allocatable :: lidar_ratio      ! lidar ratio [-] (=ceilo_ext/ceilo_back_true)
Real*8,Dimension(:,:),  Allocatable :: ceilo_first_cloud_base ! the height of the first cloud base [m]
Real*8,Dimension(:,:,:),Allocatable :: ceilo_back_true_tot  ! lidar true backscatter [m-1 sr-1]
Real*8,Dimension(:,:,:),Allocatable :: ceilo_ext_tot        ! lidar extinction coefficient [m-1 sr-1]
Real*8,Dimension(:,:,:),Allocatable :: ceilo_back_obs_tot   ! lidar observed backscatter (true back tot+ attenuation) [m-1 sr-1]]
Real*8,Dimension(:,:,:),Allocatable :: lidar_ratio_tot ! lidar ratio [-] (=ceilo_ext_tot/ceilo_back_true_tot)
End type lout_var
!
Type mpl_var ! output mpl variables 
!
Integer                               :: nx,ny,nz,nht
Real*8                                :: wavel          ! MPLL wavelenght [nm]
!
Real*8,Dimension(:,:,:,:),Allocatable :: back_true      ! mpl true backscatter [m-1 sr-1]
Real*8,Dimension(:,:,:,:),Allocatable :: ext            ! mpl extinction coefficient [m-1 sr-1]]
Real*8,Dimension(:,:,:,:),Allocatable :: back_obs       ! mpl  observed backscatter (true back + attenuation) [m-1 sr-1]]
Real*8,Dimension(:,:,:,:),Allocatable :: lidar_ratio    ! lidar ratio [-] (=ceilo_ext/ceilo_back_true)
!
Real*8,Dimension(:,:,:),Allocatable :: rayleigh_back      ! mpl molecular backscatter [m-1 sr-1]
Real*8,Dimension(:,:,:),Allocatable :: back_true_tot      ! mpl true backscatter [m-1 sr-1]
Real*8,Dimension(:,:,:),Allocatable :: ext_tot            ! mpl extinction coefficient [m-1 sr-1]]
Real*8,Dimension(:,:,:),Allocatable :: back_obs_tot       ! mpl  observed backscatter (true back + attenuation) [m-1 sr-1]]
Real*8,Dimension(:,:,:),Allocatable :: lidar_ratio_tot    ! lidar ratio [-] (=ext/back_true)

End type mpl_var
!!
!!
Type aero_var ! output mpl aero variables 
!
Integer                               :: nx,ny,nz
!
Real*8,Dimension(:,:,:),Allocatable :: back_true      ! mpl true backscatter [m-1 sr-1]
Real*8,Dimension(:,:,:),Allocatable :: ext            ! mpl extinction coefficient [m-1 sr-1]]
Real*8,Dimension(:,:,:),Allocatable :: back_obs       ! mpl  observed backscatter (true back + attenuation) [m-1 sr-1]]
Real*8,Dimension(:,:,:),Allocatable :: lidar_ratio    ! lidar ratio [-] (=ext/back_true)
!
End type aero_var
!!
!!
Type rmrad_var ! variables used for each grid point and directly computed from LUTS, in function of number of hydrometeor species
Integer   :: nht
Real*8,Dimension(:),Allocatable     :: zhh   ! [mm^6/m^3] 
Real*8,Dimension(:),Allocatable     :: dvh   ! [mm^6/m^3 * m/s]
Real*8,Dimension(:),Allocatable     :: d_dvh ! [mm^6/m^3 *(m/s)^2]
Real*8,Dimension(:),Allocatable     :: Dopp  ! [mm^6/m^3 * m/s]
Real*8,Dimension(:),Allocatable     :: Ah    ! [dB/km]
End type rmrad_var
!
Type mrad_var ! variables used for each grid point and directly computed from LUTS, in function of number of hydrometeor species
Integer   :: nht
Real*8,Dimension(:),Allocatable     :: zvv   ! [mm^6/m^3]
Real*8,Dimension(:),Allocatable     :: zvh   ! [mm^6/m^3]
Real*8,Dimension(:),Allocatable     :: RHOhvc   ! [mm^6/m^3]
Real*8,Dimension(:),Allocatable     :: Kdp   ! [deg/km]
Real*8,Dimension(:),Allocatable     :: Adp   ! [dB/km]
Real*8,Dimension(:),Allocatable     :: Av    ! [dB/km]
Real*8,Dimension(:),Allocatable     :: diff_back_phase ! [deg]
End type mrad_var
!


Type spectra_var
Integer   :: nx,ny,nz,nht
Integer   :: NFFT ! number of velocity bins and also number of FFT points
integer   :: Nave ! number of spectral averages
real*8    :: TimeSampling  ! integration time in secs
real*8    :: PRF           ! Pulse Repetition Frequency Hz
real*8    :: Lambda        ! Wavelength m
real*8    :: VNyquist      ! Nyquist velocity (m/sec)
real*8    :: NOISE_1km     ! Linear units of noise power in dBZ at 1-km

real*8,Dimension(:,:,:,:,:),Allocatable     :: zhh_spectra
real*8,Dimension(:,:,:,:,:),Allocatable     :: zvh_spectra
real*8,Dimension(:,:,:,:,:),Allocatable     :: zvv_spectra
real*8,Dimension(:,:,:,:),Allocatable       :: zhh_spectra_tot
real*8,Dimension(:,:,:,:),Allocatable       :: zvh_spectra_tot
real*8,Dimension(:,:,:,:),Allocatable       :: zvv_spectra_tot
real*8,Dimension(:),Allocatable       :: vel_bins
End type spectra_var


Type airborne_var
Integer   :: nx,ny,nz
real*8,Dimension(:,:),Allocatable     :: Zsfc_ocean
real*8,Dimension(:,:),Allocatable     :: Zsfc_land_coarse
real*8,Dimension(:,:),Allocatable     :: Zsfc_land_flat
real*8,Dimension(:,:,:),Allocatable     :: Dopp_airborne
End type airborne_var



Contains
subroutine allocate_env_var(str)
Implicit None
Type(env_var),Intent(InOut)            :: str
integer                                :: nx,ny,nz
!
nx=str%nx
ny=str%ny
nz=str%nz
!
Allocate(str%x(nx))
Allocate(str%y(ny))
Allocate(str%xlat(nx,ny))
Allocate(str%xlong(nx,ny))
Allocate(str%z(nx,ny,nz))
Allocate(str%u(nx,ny,nz))
Allocate(str%v(nx,ny,nz))
Allocate(str%w(nx,ny,nz))
Allocate(str%Ku(nx,ny,nz))
Allocate(str%Kv(nx,ny,nz))
Allocate(str%Kw(nx,ny,nz))
Allocate(str%rho_d(nx,ny,nz))
Allocate(str%qvapor(nx,ny,nz))
Allocate(str%temp(nx,ny,nz))
Allocate(str%press(nx,ny,nz))
Allocate(str%tke(nx,ny,nz))
Allocate(str%sfctemp(nx,ny))
Allocate(str%sfcwspd(nx,ny))
!
return
end subroutine allocate_env_var
!
subroutine nullify_env_var(str)
Use phys_param_mod, ONLY: zero
Implicit None
Type(env_var),Intent(InOut)            :: str
!
str%u=zero
str%v=zero
str%w=zero
str%Ku=zero
str%Kv=zero
str%Kw=zero
str%rho_d=zero
str%qvapor=zero
str%temp=zero
str%press=zero
str%tke=zero
str%sfctemp=zero
str%sfcwspd=zero
!
str%x=zero
str%y=zero
str%xlat=zero
str%xlong=zero
str%z=zero
!
return
end subroutine nullify_env_var
!
subroutine deallocate_env_var(str)
Implicit None
Type(env_var),Intent(InOut)            :: str
!
Deallocate(str%u)
Deallocate(str%v)
Deallocate(str%w)
!
Deallocate(str%Ku)
Deallocate(str%Kv)
Deallocate(str%Kw)
!
Deallocate(str%rho_d)
Deallocate(str%qvapor)
Deallocate(str%temp)
Deallocate(str%press)
Deallocate(str%sfctemp)
Deallocate(str%sfcwspd)
!
Deallocate(str%x)
Deallocate(str%y)
Deallocate(str%xlat)
Deallocate(str%xlong)
Deallocate(str%z)
!
str%nx=0
str%ny=0
str%nz=0
!
return
end subroutine deallocate_env_var
!!
!!
subroutine allocate_hydro_var(str)
Implicit None
Type(hydro_var),Intent(InOut)        :: str
Integer                              :: nx,ny,nz,nht
!
nx=str%nx
ny=str%ny
nz=str%nz
nht=str%nht
!
Allocate(str%qhydro(nx,ny,nz,nht))
Allocate(str%qnhydro(nx,ny,nz,nht))
!
return
end subroutine allocate_hydro_var
!!
subroutine initialize_hydro_var(str)
Use phys_param_mod, ONLY: zero
Implicit None
Type(hydro_var),Intent(InOut)        :: str
!
!
str%qhydro=zero
str%qnhydro=zero
!
return
end subroutine initialize_hydro_var
!!
subroutine deallocate_hydro_var(str)
Implicit None
Type(hydro_var),Intent(InOut)        :: str
!
Deallocate(str%qhydro)
Deallocate(str%qnhydro)
!
str%nx=0
str%ny=0
str%nz=0
str%nht=0
!
return
end subroutine deallocate_hydro_var
!!
!subroutine allocate_hydro50_var(str) is added by DW 2017/10/30 for P3
subroutine allocate_hydro50_var(str)
Implicit None
Type(hydro50_var),Intent(InOut)      :: str
Integer                              :: nx,ny,nz,nht
!
nx=str%nx
ny=str%ny
nz=str%nz
nht=str%nht
!
Allocate(str%qhydro(nx,ny,nz,nht))
Allocate(str%qnhydro(nx,ny,nz,nht))
Allocate(str%qir(nx,ny,nz))
Allocate(str%qib(nx,ny,nz))
!
return
end subroutine allocate_hydro50_var
!!subroutine initialize_hydro50_var(str) is added by DW 2017/10/30 for P3
subroutine initialize_hydro50_var(str)
Use phys_param_mod, ONLY: zero
Implicit None
Type(hydro50_var),Intent(InOut)      :: str
!
!
str%qhydro=zero
str%qnhydro=zero
str%qir=zero
str%qib=zero
!
return
end subroutine initialize_hydro50_var
!!subroutine deallocate_hydro50_var(str) is added by DW 2017/10/30 for P3
subroutine deallocate_hydro50_var(str)
Implicit None
Type(hydro50_var),Intent(InOut)        :: str
!
Deallocate(str%qhydro)
Deallocate(str%qnhydro)
Deallocate(str%qir)
Deallocate(str%qib)
!
str%nx=0
str%ny=0
str%nz=0
str%nht=0
!
return
end subroutine deallocate_hydro50_var
!!
subroutine allocate_hydro20_var(str)
Implicit None
Type(hydro20_var),Intent(InOut)        :: str
Integer                                :: nx,ny,nz,nbins
!
nx=str%nx
ny=str%ny
nz=str%nz
nbins=str%nbins
!
Allocate(str%ff1(nx,ny,nz,nbins))
Allocate(str%ff5(nx,ny,nz,nbins))
Allocate(str%ff6(nx,ny,nz,nbins))
!
return
end subroutine allocate_hydro20_var
!!
subroutine initialize_hydro20_var(str)
Use phys_param_mod, ONLY: zero
Implicit None
Type(hydro20_var),Intent(InOut)        :: str
!
str%ff1=zero
str%ff5=zero
str%ff6=zero
!
return
end subroutine initialize_hydro20_var
!!
subroutine deallocate_hydro20_var(str)
Implicit None
Type(hydro20_var),Intent(InOut)        :: str
!
Deallocate(str%ff1)
Deallocate(str%ff5)
Deallocate(str%ff6)
!
str%nx=0
str%ny=0
str%nz=0
str%nbins=0
!
return
end subroutine deallocate_hydro20_var
!!
subroutine allocate_hydro70_var(str)
Implicit None
Type(hydro70_var),Intent(InOut)        :: str
Integer                                :: nx,ny,nz,nbins,nht
!
nx=str%nx
ny=str%ny
nz=str%nz
nbins=str%nbins
nht=str%nht
!
Allocate(str%fn1(nx,ny,nz,nbins))
Allocate(str%fm1(nx,ny,nz,nbins))
Allocate(str%qhydro(nx,ny,nz,nht))
Allocate(str%qnhydro(nx,ny,nz,nht))
Allocate(str%naero(nx,ny,nz))
!
return
end subroutine allocate_hydro70_var
!!
subroutine initialize_hydro70_var(str)
Use phys_param_mod, ONLY: zero
Implicit None
Type(hydro70_var),Intent(InOut)        :: str
!
str%fn1=zero
str%fm1=zero
str%qhydro=zero
str%qnhydro=zero
str%naero=zero
!
return
end subroutine initialize_hydro70_var
!!
subroutine deallocate_hydro70_var(str)
Implicit None
Type(hydro70_var),Intent(InOut)        :: str
!
Deallocate(str%fn1)
Deallocate(str%fm1)
Deallocate(str%qhydro)
Deallocate(str%qnhydro)
Deallocate(str%naero)
!
str%nx=0
str%ny=0
str%nz=0
str%nbins=0
!
return
end subroutine deallocate_hydro70_var
!!
!!
subroutine allocate_rmout_var(str)
Implicit None
Type(rmout_var),Intent(InOut)           :: str
Integer                                 :: nx,ny,nz,nht
!
nx=str%nx
ny=str%ny
nz=str%nz
nht=str%nht
!
Allocate(str%Zhh(nx,ny,nz,nht))

Allocate(str%Dopp(nx,ny,nz,nht))
Allocate(str%dDVh(nx,ny,nz,nht))
Allocate(str%SWt(nx,ny,nz,nht))
Allocate(str%DVh(nx,ny,nz,nht))
Allocate(str%Dopp90(nx,ny,nz,nht))
Allocate(str%dDVh90(nx,ny,nz,nht))
Allocate(str%Ah(nx,ny,nz,nht))
!
Allocate(str%Zhh_tot(nx,ny,nz))
Allocate(str%Dopp_tot(nx,ny,nz))
Allocate(str%dDVh_tot(nx,ny,nz))
Allocate(str%SWt_tot(nx,ny,nz))
Allocate(str%DVh_tot(nx,ny,nz))
Allocate(str%Dopp90_tot(nx,ny,nz))
Allocate(str%dDVh90_tot(nx,ny,nz))
Allocate(str%Ah_tot(nx,ny,nz))
!
Allocate(str%elev(nx,ny,nz))
Allocate(str%azim(nx,ny,nz))
Allocate(str%range(nx,ny,nz))
Allocate(str%sw_t(nx,ny,nz))
Allocate(str%sw_s(nx,ny,nz))
Allocate(str%sw_v(nx,ny,nz))
!
return
end subroutine allocate_rmout_var
!
subroutine m999_rmout_var(str)
Use phys_param_mod, ONLY: m999
Implicit None
Type(rmout_var),Intent(InOut)            :: str
!
str%Zhh=m999
str%Dopp=m999
str%dDVh=m999
str%SWt=m999
str%DVh=m999
str%Dopp90=m999
str%dDVh90=m999
str%Ah=m999
!
str%Zhh_tot=m999
str%Dopp_tot=m999
str%dDVh_tot=m999
str%SWt_tot=m999
str%DVh_tot=m999
str%Dopp90_tot=m999
str%dDVh90_tot=m999
str%Ah_tot=m999
!
str%elev=m999
str%azim=m999
str%range=m999
str%sw_t=m999
str%sw_s=m999
str%sw_v=m999
!
return
end subroutine m999_rmout_var
!
subroutine deallocate_rmout_var(str)
Implicit None
Type(rmout_var),Intent(InOut)            :: str

Deallocate(str%Zhh)
Deallocate(str%Dopp)
Deallocate(str%dDVh)
Deallocate(str%SWt)
Deallocate(str%DVh)
Deallocate(str%Dopp90)
Deallocate(str%dDVh90)
Deallocate(str%Ah)
!
Deallocate(str%Zhh_tot)
Deallocate(str%Dopp_tot)
Deallocate(str%dDVh_tot)
Deallocate(str%SWt_tot)
Deallocate(str%DVh_tot)
Deallocate(str%Dopp90_tot)
Deallocate(str%dDVh90_tot)
Deallocate(str%Ah_tot)
!
Deallocate(str%elev)
Deallocate(str%azim)
Deallocate(str%range)
Deallocate(str%sw_t)
Deallocate(str%sw_s)
Deallocate(str%sw_v)
!
return
end subroutine deallocate_rmout_var
!!
!!
subroutine allocate_mout_var(str)
Implicit None
Type(mout_var),Intent(InOut)            :: str
Integer                                 :: nx,ny,nz,nht
!
nx=str%nx
ny=str%ny
nz=str%nz
nht=str%nht
!
Allocate(str%Zvv(nx,ny,nz,nht))
Allocate(str%Zvh(nx,ny,nz,nht))
Allocate(str%RHOhv(nx,ny,nz,nht))
Allocate(str%Zdr(nx,ny,nz,nht))
Allocate(str%LDRh(nx,ny,nz,nht))
!
Allocate(str%Kdp(nx,ny,nz,nht))
Allocate(str%Adp(nx,ny,nz,nht))
Allocate(str%Av(nx,ny,nz,nht))
Allocate(str%diff_back_phase(nx,ny,nz,nht))
!
Allocate(str%Zvv_tot(nx,ny,nz))
Allocate(str%Zvh_tot(nx,ny,nz))
Allocate(str%RHOhv_tot(nx,ny,nz))
Allocate(str%Zdr_tot(nx,ny,nz))
Allocate(str%LDRh_tot(nx,ny,nz))
!
Allocate(str%Kdp_tot(nx,ny,nz))
Allocate(str%Adp_tot(nx,ny,nz))
Allocate(str%Av_tot(nx,ny,nz))
Allocate(str%diff_back_phase_tot(nx,ny,nz))
!
return
end subroutine allocate_mout_var
!
subroutine m999_mout_var(str)
use phys_param_mod, ONLY: m999
Implicit None
Type(mout_var),Intent(InOut)            :: str
!
str%Zvv=m999
str%Zvh=m999
str%RHOhv=m999
str%Zdr=m999
str%LDRh=m999
str%Kdp=m999
str%Adp=m999
str%Av=m999
str%diff_back_phase=m999
!
str%Zvv_tot=m999
str%Zvh_tot=m999
str%RHOhv_tot=m999
str%Zdr_tot=m999
str%LDRh_tot=m999
str%Kdp_tot=m999
str%Adp_tot=m999
str%Av_tot=m999
str%diff_back_phase_tot=m999
!
return
end subroutine m999_mout_var
!
subroutine deallocate_mout_var(str)
Implicit None
Type(mout_var),Intent(InOut)            :: str

Deallocate(str%Zvv)
Deallocate(str%Zvh)
Deallocate(str%RHOhv)
Deallocate(str%Zdr)
Deallocate(str%LDRh)
Deallocate(str%Kdp)
Deallocate(str%Adp)
Deallocate(str%Av)
!
Deallocate(str%Zvv_tot)
Deallocate(str%Zvh_tot)
Deallocate(str%RHOhv_tot)
Deallocate(str%Zdr_tot)
Deallocate(str%LDRh_tot)
Deallocate(str%diff_back_phase)
Deallocate(str%Kdp_tot)
Deallocate(str%Adp_tot)
Deallocate(str%Av_tot)
Deallocate(str%diff_back_phase_tot)
!
return
end subroutine deallocate_mout_var
!!
!!
subroutine allocate_lout_var(str)!4th dimension added by oue for RAMS
Implicit None
Type(lout_var),Intent(InOut)   :: str
Integer                         :: nx,ny,nz,nht
!
nx=str%nx
ny=str%ny
nz=str%nz
nht=str%nht
!
Allocate(str%ceilo_back_true(nx,ny,nz,nht))
Allocate(str%ceilo_ext(nx,ny,nz,nht))
Allocate(str%ceilo_back_obs(nx,ny,nz,nht))
Allocate(str%lidar_ratio(nx,ny,nz,nht))
Allocate(str%ceilo_first_cloud_base(nx,ny))
Allocate(str%ceilo_back_true_tot(nx,ny,nz))
Allocate(str%ceilo_ext_tot(nx,ny,nz))
Allocate(str%ceilo_back_obs_tot(nx,ny,nz))
Allocate(str%lidar_ratio_tot(nx,ny,nz))
!
return
end subroutine allocate_lout_var
!
subroutine m999_lout_var(str)
Use phys_param_mod, ONLY: zero
Implicit None
Type(lout_var),Intent(InOut)   :: str
!
str%ceilo_back_true=zero
str%ceilo_ext=zero
str%ceilo_back_obs=zero
str%lidar_ratio=zero
str%ceilo_first_cloud_base=zero
str%ceilo_back_true_tot=zero
str%ceilo_ext_tot=zero
str%ceilo_back_obs_tot=zero
str%lidar_ratio_tot=zero
!
return
end subroutine m999_lout_var
!
subroutine deallocate_lout_var(str)
Implicit None
Type(lout_var),Intent(InOut)   :: str
!
Deallocate(str%ceilo_back_true)
Deallocate(str%ceilo_ext)
Deallocate(str%ceilo_back_obs)
Deallocate(str%lidar_ratio)
Deallocate(str%ceilo_first_cloud_base)
Deallocate(str%ceilo_back_true_tot)
Deallocate(str%ceilo_ext_tot)
Deallocate(str%ceilo_back_obs_tot)
Deallocate(str%lidar_ratio_tot)
!
str%nx=0; str%ny=0; str%nz=0
!
return
end subroutine deallocate_lout_var
!!
subroutine allocate_mpl_var(str)
Implicit None
Type(mpl_var),Intent(InOut)   :: str
Integer                         :: nx,ny,nz,nht
!
nx=str%nx
ny=str%ny
nz=str%nz
nht=str%nht
!
Allocate(str%back_true(nx,ny,nz,nht))
Allocate(str%ext(nx,ny,nz,nht))
Allocate(str%back_obs(nx,ny,nz,nht))
Allocate(str%lidar_ratio(nx,ny,nz,nht))
!
Allocate(str%rayleigh_back(nx,ny,nz))
Allocate(str%back_true_tot(nx,ny,nz))
Allocate(str%ext_tot(nx,ny,nz))
Allocate(str%back_obs_tot(nx,ny,nz))
Allocate(str%lidar_ratio_tot(nx,ny,nz))

return
end subroutine allocate_mpl_var
!
subroutine m999_mpl_var(str)
Use phys_param_mod, ONLY: zero
Implicit None
Type(mpl_var),Intent(InOut)   :: str
!
str%back_true=zero
str%ext=zero
str%back_obs=zero
str%lidar_ratio=zero
!
str%rayleigh_back=zero
str%back_true_tot=zero
str%ext_tot=zero
str%back_obs_tot=zero
str%lidar_ratio_tot=zero
!
return
end subroutine m999_mpl_var
!
subroutine deallocate_mpl_var(str)
Implicit None
Type(mpl_var),Intent(InOut)   :: str
!
Deallocate(str%back_true)
Deallocate(str%ext)
Deallocate(str%back_obs)
Deallocate(str%lidar_ratio)
!
Deallocate(str%rayleigh_back)
Deallocate(str%back_true_tot)
Deallocate(str%ext_tot)
Deallocate(str%back_obs_tot)
Deallocate(str%lidar_ratio_tot)
!
str%nx=0; str%ny=0; str%nz=0 ; str%nht=0
!
return
end subroutine deallocate_mpl_var
!!
!!
subroutine allocate_aero_var(str)
Implicit None
Type(aero_var),Intent(InOut)   :: str
Integer                         :: nx,ny,nz
!
nx=str%nx
ny=str%ny
nz=str%nz
!
Allocate(str%back_true(nx,ny,nz))
Allocate(str%ext(nx,ny,nz))
Allocate(str%back_obs(nx,ny,nz))
Allocate(str%lidar_ratio(nx,ny,nz))
!
return
end subroutine allocate_aero_var
!
subroutine m999_aero_var(str)
Use phys_param_mod, ONLY: zero
Implicit None
Type(aero_var),Intent(InOut)   :: str
!
str%back_true=zero
str%ext=zero
str%back_obs=zero
str%lidar_ratio=zero
!
return
end subroutine m999_aero_var
!
subroutine deallocate_aero_var(str)
Implicit None
Type(aero_var),Intent(InOut)   :: str
!
Deallocate(str%back_true)
Deallocate(str%ext)
Deallocate(str%back_obs)
Deallocate(str%lidar_ratio)
!
str%nx=0; str%ny=0; str%nz=0 
!
return
end subroutine deallocate_aero_var
!!
!!
subroutine allocate_rmrad_var(str)
Implicit None
Type(rmrad_var),Intent(InOut)   :: str
Integer                         :: nht
!
nht=str%nht
!
Allocate(str%zhh(nht))
Allocate(str%dvh(nht))
Allocate(str%d_dvh(nht))
Allocate(str%Dopp(nht))
Allocate(str%Ah(nht))
!
return
end subroutine allocate_rmrad_var
!!
subroutine nullify_rmrad_var(str)
Use phys_param_mod, ONLY: zero
Implicit None
Type(rmrad_var),Intent(InOut)   :: str
!
str%zhh=zero
str%dvh=zero
str%d_dvh=zero
str%Dopp=zero
str%Ah=zero
!
return
end subroutine nullify_rmrad_var
!
subroutine deallocate_rmrad_var(str)
Implicit None
Type(rmrad_var),Intent(InOut)   :: str
!
Deallocate(str%zhh)
Deallocate(str%dvh)
Deallocate(str%d_dvh)
Deallocate(str%Dopp)
Deallocate(str%Ah)
!
str%nht=0
!
return
end subroutine deallocate_rmrad_var
!
subroutine allocate_mrad_var(str)
Implicit None 
Type(mrad_var),Intent(InOut)   :: str
Integer                        :: nht
!
nht=str%nht
!
Allocate(str%zvv(nht))
Allocate(str%zvh(nht))
Allocate(str%RHOhvc(nht))
Allocate(str%Kdp(nht))
Allocate(str%Adp(nht))
Allocate(str%Av(nht))
Allocate(str%diff_back_phase(nht))
!
return
end subroutine allocate_mrad_var
!!
!!
subroutine nullify_mrad_var(str)
Use phys_param_mod, ONLY: zero
Implicit None
Type(mrad_var),Intent(InOut)   :: str
!
str%zvv=zero
str%zvh=zero
str%RHOhvc=zero
str%Kdp=zero
str%Adp=zero
str%Av=zero
str%diff_back_phase=zero
!
return
end subroutine nullify_mrad_var
!
subroutine deallocate_mrad_var(str)
Implicit None
Type(mrad_var),Intent(InOut)   :: str
!
Deallocate(str%zvv)
Deallocate(str%zvh)
Deallocate(str%RHOhvc)
Deallocate(str%Kdp)
Deallocate(str%Adp)
Deallocate(str%Av)
Deallocate(str%diff_back_phase)
!
str%nht=0
!
return
end subroutine deallocate_mrad_var 


subroutine allocate_spectra_var(str,radID)
Use phys_param_mod, ONLY: zero
Implicit None
Type(spectra_var),Intent(InOut)        :: str
Integer, Intent(In) :: radID
Integer                                :: nx,ny,nz,nbins,nht
!
nx=str%nx
ny=str%ny
nz=str%nz
nbins=str%NFFT
nht=str%nht
!
Allocate(str%vel_bins(nbins))
Allocate(str%zhh_spectra(nx,ny,nz,nbins,nht))
Allocate(str%zhh_spectra_tot(nx,ny,nz,nbins))
if(radID /=1)then
  Allocate(str%zvh_spectra(nx,ny,nz,nbins,nht))
  Allocate(str%zvv_spectra(nx,ny,nz,nbins,nht))
  Allocate(str%zvh_spectra_tot(nx,ny,nz,nbins))
  Allocate(str%zvv_spectra_tot(nx,ny,nz,nbins))
end if
! initialize allocated variables
str%vel_bins=zero
str%zhh_spectra=zero
str%zhh_spectra_tot=zero
if(radID /=1)then
  str%zvh_spectra=zero
  str%zvv_spectra=zero
  str%zvh_spectra_tot=zero
  str%zvv_spectra_tot=zero
end if
!
return
end subroutine allocate_spectra_var
!!
subroutine deallocate_spectra_var(str,radID)
Use phys_param_mod, ONLY: zero,izero
Implicit None
Type(spectra_var),Intent(InOut)        :: str
Integer, Intent(In) :: radID
!
Deallocate(str%zhh_spectra)
Deallocate(str%zhh_spectra_tot)
Deallocate(str%vel_bins)
if(radID /=1)then
  Deallocate(str%zvh_spectra)
  Deallocate(str%zvv_spectra)
  Deallocate(str%zvh_spectra_tot)
  Deallocate(str%zvv_spectra_tot)
end if
!
str%nx=0
str%ny=0
str%nz=0
str%NFFT=0
str%Lambda = zero
str%VNyquist = zero
str%Nave = izero

!
return
end subroutine deallocate_spectra_var
!!

subroutine allocate_airborne_var(str)
Use phys_param_mod, ONLY: zero
Implicit None
Type(airborne_var),Intent(InOut)        :: str
Integer                                :: nx,ny,nz
!
nx=str%nx
ny=str%ny
nz=str%nz
!
Allocate(str%Zsfc_ocean(nx,ny))
Allocate(str%Zsfc_land_coarse(nx,ny))
Allocate(str%Zsfc_land_flat(nx,ny))
Allocate(str%Dopp_airborne(nx,ny,nz))
! initialize allocated variables
str%Zsfc_ocean=zero
str%Zsfc_land_coarse=zero
str%Zsfc_land_flat=zero
str%Dopp_airborne=zero
!
return
end subroutine allocate_airborne_var
!!
subroutine deallocate_airborne_var(str)
Use phys_param_mod, ONLY: zero
Implicit None
Type(airborne_var),Intent(InOut)        :: str
!
Deallocate(str%Zsfc_ocean)
Deallocate(str%Zsfc_land_coarse)
Deallocate(str%Zsfc_land_flat)
Deallocate(str%Dopp_airborne)
!
str%nx=0
str%ny=0
str%nz=0
!
return
end subroutine deallocate_airborne_var
!!

end module crsim_mod
