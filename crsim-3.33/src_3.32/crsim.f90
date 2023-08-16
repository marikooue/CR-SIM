  !!-----------------------------------------------------------------------------------------------------------------
  !!-----------------------------------------------------------------------------------------------------------------
  !! *PROGRAM* : crsim
  !!  @version crsim 3.x
  !!
  !!  *SRC_FILE*
  !!  /crsim/src/crsim.f90
  !!
  !!
  !!  *LAST CHANGES* 
  !! 
  !!  Sep 15- Sep 23  2015   - A.T.   changes introduced in order to compute the radial Doppler velocity and spectrum width
  !!  Oct 20 2015            - A.T.   corrected issue with dlog10(0.d0) appearing
  !!  Jan 05 2016            - A.T.   bug when computing the Zdr_tot and LDRh_tot is now corrected.Unrealistic Zdr_tot and
  !!                                  LDRh_tot were due to incorrectly placed line for changing of units from linear to logaritmic
  !!                                  for rmout%Zhh_tot
  !!  Jan 11 2016            - A.T.   control variable hmask is removed as it is now obsolate  
  !!  Jan 11 2016            - A.T.   included threshold value Zthr (= -100 dBZ)  of minimum reflectivity such that:
  !!                                     1. if Zhh>Zthr => Zhh,Ah and all Doopler-related variables exist
  !!                                     2. if Zvv>Zthr => Zvv, Av exist
  !!                                     3. if both Zhh,Zvv > Zthr => all other polarimetric variables exist (Zvh,Zdr,LDR,Kdp,Zdr)
  !! Jan 12 2016             - A.T.   included new parameters in Configuration file, the minimum values of the input mixing ratios
  !!                                  for species. Mixing ratio values smaller or equal than specified thresholds are set to 0. In the
  !!                                  case of MP_PHYSICS=20, in the pixels where computed mixing ratio is <= specified threshold,
  !!                                  the input bin total number distributions are set to zeros.
  !! Jan 13 2016             - A.T.   Zvh was not reported in logaritmic but linear units. It is now corrected and all the  reflecivity 
  !!                                  related variables are reported in logaritmic units.                                   
  !! Jan 14-JAn 15  2016     - A.T.   Included computation of spectrum width contributions due to turnulence, wind shear in radar
  !!                                  volume and cross-wind.
  !! Feb 17 2016             - A.T.   Icluded computation od cloud lidar (ceilometer) measurements
  !! Feb 18 2016             - A.T.   Corrected bug in sign when computing Doppler velocity at elevations different than 90 deg. All 
  !!                                  velocties are now positive away from the radar (radar notation).
  !!
  !! March 26 2016           -A.T.    Subroutine get_hydro_vars is modified to work for MP_PHYSICS=9. For this  reason get_hydro_vars is
  !!                                  splitted out to 2 subroutines: get_hydro10_vars and get_hydro09_vars. In the case of Morrisson scheme,
  !!                                  a possibility of the input cloud liquid total concentration /=0 is taken into account. 
  !! May 2  2016             -A.T.    Included computation od micropulse lidar (MPL) measurements for cloud and ice hydrometeors
  !! May 9  2016             -A.T.    Included computation of aerosol average profile in MPL measurements
  !! May 12 2016             -A.T.    Included computation of molecular backscatter in MPL simulated measurements  
  !!
  !! Sep 19 2016             -M.O.    Incorporated Thompson microphysics
  !! Mar 23 2017             -M.O.    Included post processing: estimation of ARSCL products and MWR LWP
  !! JUL 17 2017             -M.O     Incorporated  ICON 2-moment microphysics. 
  !! JUL 21 2017             -M.O     Incorporated  RAMS 2-moment microphysics.(8 categories) 
  !! JUN 30 2017             -K. YU   Implemented "__PRELOAD_LUT__" mode and related subroutines are implemented in crsim_subrs.f90. 
  !!                                  OpenMP parallelization is applied.
  !! 
  !! Oct 30 2017            - D.W.    Incorporated P3 microphysics (added by DW)
  !! Mar 13 2018            -M.O.     Corrected bug when computing the ARSCL radar_zhh_min 
  !! Apr    2018            -M.O      Incorporated SAM warm bin microphysics (MP_PHYSICS=70)
  !! Jun 17 2018            -M.O      Incorporated SAM Morrison 2-moment microphysics (MP_PHYSICS=75)
  !! MAY 27 2018            -M.O      Incorporated Doppler spectra generator for bulk moment microphysics
  !!                                  Added an option to use hail instead of graupel for MP_PHYSICS=10 (WRF Morrison 2 moment)
  !!                                  Modified ARSCL's radar_zhh_min: Use rmout%range, but use env%z when vertically-pointing
  !!                                  measurements (conf%elev>89.0) 
  !! Sep 03 2018            -M.O.     Fixed ww=env%v(ix,iy,iz)--> ww=env%w(ix,iy,iz) in crsim.f90.
  !! Nov    2018            -A.T.     Introduced indentation.
  !! Nov    2018            -A.T.     Module phys_param_mod introduced here and throghout the CR-SIM code
  !! Nov    2018            -A.T.     Modified parts of the code with Doppler spectra generator for bulk moment microphysics 
  !! Nov    2018            -A.T.     The main code cleaned, description updated
  !! Dec    2018            -A.T.     Introduced progress messages in procentage
  !! Aug    2019            -A.T.     Conversion from real to double  precision variables everywhere
  !! Jul,Aug   2019       -M.O,A.T.   Work on introducing airborne scanning observations    
  !! Dec    2019            -A.T.     Added computation of cross-correlation coefficient rho_hv 
  !!
  !!
  !!  *DESCRIPTION* 
  !!
  !!  The CR-SIM uses the inputs from the high resolution Weather Research and
  !!  Forecasting Model (WRF), for microphysical schemes MP_PHYSICS=8 (Thompson
  !!  scheme), MP_PHYSICS=10 (Morrison double-moment scheme), MP_PHYSICS=9 (Milbrandt
  !!  and Yau double-moment scheme), MP_PHYSICS=20 (bin explicit scheme), and
  !!  MP_PHYSICS=50 (Predicted particle properties scheme) and computes “idealized”
  !!  forward modeled scanning (or vertical-pointing) radar observations (Figure 1)
  !!  and profiling lidar observables. The latest version can use the inputs from the
  !!  ICOsahedral Non-hydrostatic general circulation model (ICON, MP_PHYSICS=30) and
  !!  the Regional Atmospheric Modeling System (RAMS, MP_PHYSICS=40), and System for
  !!  Atmospheric Modeling (SAM, MP_PHYSICS=70 for bin microphysics, MP_PHYSICS=75 for
  !!  Morrison double-moment scheme).
  !!  
  !!  The CR-SIM employs the T-matrix method for computation of scattering
  !!  characteristics such as cloud water, cloud ice, rain, snow, graupel, hail and 
  !!  allows the specifications of the following radar frequencies for scattering
  !!  calculations: 3 GHz, 5.5 GHz, 9.5 GHz, 35 GHz and 94 GHz. 
  !!
  !!  Given the particles size distributions (explicit or reconstructed from the
  !!  provided PSD moments), polarimetric radar variables can be calculated if the
  !!  scattering amplitudes are known. The complex scattering amplitudes are
  !!  pre-computed and stored as look-up tables (LUTs) for equally spaced particle
  !!  sizes using the Mishchenko's T-matrix code for a non-spherical particle at a
  !!  fixed orientation (Mishchenko, 2000) and for elevation angles from 0° to 90°
  !!  with a spacing of 1°, for specified radar frequencies, temperatures and
  !!  different possibilities of particles densities and aspect ratios. 
  !!
  !!  It is worth noting that some particle characteristics needed to simulate
  !!  polarimetric variables such as the shape or the statistical properties of the
  !!  particle orientations are not explicitly specified in the model. Thus certain
  !!  assumptions have to be made, based on much available information as possible.
  !!  For example, the mean canting angles of all hydrometeors are assumed to be 0°.
  !!  This assumption enables the use of Ryzkov et al. (2010) formulations for
  !!  computation of polarimetric variables and especially the use of simple
  !!  expressions for the angular moments, which in turn make possible that the width
  !!  of 2D Gaussian distribution of canting angle that is different for various
  !!  hydrometeor species can be specified. 
  !!
  !!
  !!  *INPUT*
  !!
  !!  The I/O files are specified in th ecommand line while the 
  !!  configuration parameters are specified in the file "PARAMETERS"
  !!  For additional information, please see CR-SIM User Guide.
  !!
  !!  *OUTPUT* 
  !!
  !!  A set of netcdf output files is generated. The main output file is the one
  !!  that refers to the CR-SIM output for the total hydrometeor content.
  !!  The additional output files will be created for each hydrometeor class. 
  !!  Per default, they will have extensions “_cloud”, ”_rain”, ”_ice”, “_snow”,
  !!  “_graupel”, “_hail”, "drizzle", "aggregare", "smallice", "unrimed_ice", 
  !!  and parimedice" before the “.nc”.
  !!  The polarimetric variables simulated include reflectivity at vertical and
  !!  horizontal polarization, differential reflectivity, specific differential phase,
  !!  specific attenuation at horizontal and vertical polarization, specific
  !!  differential attenuation, and linear depolarization ratio. Also simulated are
  !!  the vertical reflectivity weighted velocity, mean Doppler velocity and spectrum
  !!  width. The computation of fall velocities is consistent to the computational
  !!  method used in the specific WRF microphysical package.
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
  !!  You should have received a copy of the GNU  General Public License
  !!  along with this program.  If not, see <http://www.gnu.org/licenses/>.
  !!
  !!----------------------------------------------------------------------------------------------------------------- 
  !!-----------------------------------------------------------------------------------------------------------------
  !!
  program crsim
  Use crsim_mod
  Use wrf_var_mod
  Use postprocess_mod  !added by oue 2017.03.23
  Use crsim_luts_mod !--Added by oue for preloading LUTs of all densities
  use omp_lib
  Use phys_param_mod, ONLY: pi, N_layers, m999,d2r, Zthr ! AT -- Added here and throghout the CR-SIM code 
  !
  Implicit None
  
  Type(conf_var)                 :: conf
  
  Type(wrf_var)                  :: wrf
  Type(wrf_var_mp09)             :: mp09
  Type(wrf_var_mp08)             :: mp08 !added by oue 2016/09/19
  Type(wrf_var_mp10)             :: mp10
  Type(wrf_var_mp20)             :: mp20
  Type(wrf_var_mp30)             :: mp30 !added by oue 2017/07/17 for ICON
  Type(wrf_var_mp40)             :: mp40 !added by oue 2017/07/21 for RAMS
  Type(wrf_var_mp50)             :: mp50 ! added by DW 2017/10/30 for P3
  Type(wrf_var_mp70)             :: mp70 ! added by oue for SAM warm bin
  Type(wrf_var_mp75)             :: mp75 ! added by oue for SAM morr
  !
  Type(env_var)                  :: env 
  Type(hydro_var)                :: hydro
  Type(hydro20_var)              :: hydro20
  Type(hydro50_var)              :: hydro50 ! added by DW 2017/10/30 for P3
  Type(hydro70_var)              :: hydro70 ! added by oue for SAM warm bin
  !
  Type(mout_var)                 :: mout
  !Type(mrad_var)                 :: mrad
  !
  Type(rmout_var)                 :: rmout
  !Type(rmrad_var)                 :: rmrad
  !
  Type(lout_var)                  :: lout
  Type(mpl_var)                   :: mpl
  Type(aero_var)                  :: aero
  Type(arscl_var)                 :: arscl !added by oue 2017.03.23
  Type(mwr_var)                   :: mwr !added by oue 2017.03.23
  !
  Type(scatt_type_var),Dimension(:),Allocatable       :: scatt_type 
  !
  integer                        :: nsc ! number of scatt. species in LUTS 
  integer                        :: isc
  Character(len=256)             :: ConfigInpFile
  Integer*4                      :: ix,iy,iz,iht,ihtc
  integer                        :: nx,ny,nz,nt
  integer                        :: status
  !
  real*8                         :: elev,elevx
  real*8                         :: azim
  real*8                         :: rr  ! radar range [m] 
  !
  real*8                         :: uu,vv ! horizontal wind componenets in m/s
  real*8                         :: ww ! vertical wind component in m/s
  real*8                         :: tke ! turbulence kinetic energy m^2/s^2
  real*8,Dimension(:,:,:), Allocatable    :: qqc ! [kg/m^3] water content 
  !
  Character(len=256)             :: path
  Character(len=100)             :: filename
  Character(len=10)              :: extention
  Character(len=10)              :: hstring
  integer                        :: nn,nnf
  Character(len=356)             :: WRFmpInputFile
  Character(len=356)             :: OutFileName
  Character(len=340)             :: LutFileName
  !
  real*8                         :: qsum   ! [kg/kg] total mixing ratio at specific position
  real*8                         :: qq     ! [kg/kg] mixing ratio at specific position
  integer                        :: nkr    ! total number of bins in each category (=33 iper default for MP_PHYSICS=20 and 70)
  integer                        :: ii
  !
  real*8                         :: zhh,zvv,zvh,RHOhvc,dvh,d_dvh,Dopp,Kdp,Adp,Ah,Av,diff_back_phase
  real*8                         :: sw_t2  ! [m2/s2] ! squared spectrum width due to turbulence only
  real*8                         :: sw_h2  ! [m2/s2] ! squared spectrum width due to diff. hydrometeor terminal vel. at diff. sizes
  real*8                         :: sw_s2  ! [m2/s2] ! squared spectrum width due to wind shear in radar volume
  real*8                         :: sw_v2  ! [m2/s2] ! squared spectrum width due to cross wind 
  real*8                         :: work,work1
  integer                        :: nht
  !
  Real*8                             :: ceilo_back_true ! true (unatenuated) lidar backscatter [m sr]^-1
  Real*8                             :: ceilo_ext       ! lidar extinction coefficient [m sr]^-1
  Real*8,Dimension(:), Allocatable   :: ceilo_back_obs  ! observed lidar backscatter [m sr]^-1
  Real*8                             :: ceilo_first_cloud_base  ! the height of the first cloud base [m]
  !
  Real*8                             :: mpl_back_true ! true (unatenuated) lidar backscatter [m sr]^-1
  Real*8                             :: mpl_ext       ! lidar extinction coefficient [m sr]^-1
  Real*8,Dimension(:), Allocatable   :: mpl_back_obs  ! observed lidar backscatter [m sr]^-1
  !-Used for ARSCL simulation, added by oue 2017.03.23
  Real*8,Dimension(:,:,:), Allocatable   :: mpl_back_hydro  ! MPL hydrometeor backscatter [m sr]^-1 
  Real*8,Dimension(:,:,:), Allocatable   :: mpl_back_aero  ! aerosol backscatter [m sr]^-1 
  Real*8,Dimension(:), Allocatable   :: radar_zhh_min  ! radar limitation sensitivity [dBZ] 
  Integer*2,Dimension(:), Allocatable  :: cloud_detection_flag,cloud_mask ! ARSCL cloud mask
  Real*8,Dimension(:), Allocatable  :: cloud_layer_base,cloud_layer_top ! ARSCL cloud layers
  !- Used for MWR LWP
  Real*8,Dimension(:,:,:), Allocatable   :: lwc  ! Liquid water content [kg m^-3]
  !--
  !
  integer :: clck_counts_beg, clck_counts_end, clck_rate
  Real*8 :: wstart, wend, wend2, wstart3, wend3, wstart4, wend4
  
  !--!--Added by oue for preloading LUTs of all densities
  Integer :: nden,iden
  Real*8,Dimension(:),Allocatable :: lden
  !--
  !--!--Added by oue for Doppler spectra simulation
  Type(spectra_var)               :: spectra ! added by oue for Doppler spectra simulation
  Integer                         :: ibin,nfft
  real*8                          :: sw_dyn ! [m/s] ! spectra broadening due to dynamics (sw_t+sw_s+sw_v) 
  real*8                          :: w_r ! [m/s] ! radial component of wind field 
  real*8                          :: dist_from_radar ! [m] distance from radar 
  Real*8,Dimension(:),Allocatable ::  spectra_bins,zhh_spectra,zvh_spectra,zvv_spectra ![mm^6 m^-3]
  !--
  
  Integer  :: psum_xy_tot,psum_xy,percentage,percentage_saved
  
  !#ifndef __MAX_TH_NUM__
  !#define __MAX_TH_NUM__ 16
  !#endif
  Integer :: MAX_TH_NUM
  !Integer, PARAMETER :: MAX_TH_NUM=__MAX_TH_NUM__
  Integer :: tid
  
  !Type(rmrad_var), Dimension(MAX_TH_NUM) :: rmrad_vec
  !Type(mrad_var), Dimension(MAX_TH_NUM) :: mrad_vec
  Type(rmrad_var), Dimension(:), Allocatable :: rmrad_vec
  Type(mrad_var),  Dimension(:), Allocatable :: mrad_vec
  !wstart = omp_get_wtime();
  !
  !
    call system_clock ( clck_counts_beg, clck_rate )
    !
    ! get command line arguments
    Call get_command_line_arguments(status)  
    If (status /=0) Then
      write(*,*) 'Error when geting the command line arguments'
      Call Exit(1)
    Endif
    !
    !
    ! read conf file
    call ReadConfParameters(trim(ConfigInpFile),conf)
    !
    ! OMP Thread related setting
#ifdef __PRELOAD_LUT__
      !write(*,*) 'conf%OMPThreadNum : ', conf%OMPThreadNum
  
      
      if(conf%OMPThreadNum < 1) then
        conf%OMPThreadNum = 1
      endif
!!$    conf%OMPThreadNum = min(conf%OMPThreadNum, omp_get_max_threads())
  
#else
      conf%OMPThreadNum = 1
      write(*,*) 'Only 1 OpenMP thread is supported due to __PRELOAD_LUT__ is not used!!!'
#endif
  
    ! reset the conf%OMPThreadNum where  __PRELOAD_LUT_ not supported
    if ( (conf%MP_PHYSICS==20) .or. (conf%MP_PHYSICS==50) .or. (conf%MP_PHYSICS==70) ) then
      MAX_TH_NUM = conf%OMPThreadNum
      if (MAX_TH_NUM /= 1) then 
        write(*,*) 'Number of OpenMP threads reset to 1 for selected MP_PHYSICS'
        conf%OMPThreadNum = 1
      endif 
    endif
    !
    !
    MAX_TH_NUM = conf%OMPThreadNum
    !$ call omp_set_num_threads(MAX_TH_NUM) ! AT April 2019
    !$ write(*,*) 'OpenMP Thread Number: ', MAX_TH_NUM ! AT April 2019

    allocate(rmrad_vec(MAX_TH_NUM), mrad_vec(MAX_TH_NUM))
  
    !
    ! read dimensions from the input wrf netcdf file
    if(conf%MP_PHYSICS==40) then !read RAMS added by oue to read RAMS
      call ReadInpRAMS_dim(trim(conf%WRFInputFile),wrf,status)
    else if ((conf%MP_PHYSICS==70) .or. (conf%MP_PHYSICS==75)) then!read SAM: added by oue to read SAM
      call ReadInpSAM_dim(trim(conf%WRFInputFile),wrf,status)
    else
      call ReadInpWRF_dim(trim(conf%WRFInputFile),wrf,status)
    endif
    ! 
    !------------------------------------
    ! save dimensions
    nt=wrf%nt ; nx=wrf%nx ; ny=wrf%ny ; nz=wrf%nz
    nsc=conf%nht
    !------------------------------------
    !! define dimensions for the env structure
    env%nx=wrf%nx ; env%ny=wrf%ny ; env%nz=wrf%nz
    !
    if (MIN(conf%ix_end,conf%ix_start)>0) env%nx=conf%ix_end-conf%ix_start + 1
    if (MIN(conf%iy_end,conf%iy_start)>0) env%ny=conf%iy_end-conf%iy_start + 1
    if (MIN(conf%iz_end,conf%iz_start)>0) env%nz=conf%iz_end-conf%iz_start + 1
    !----------------------------------------
    !! allocate env vars 
    call allocate_env_var(env)
    call nullify_env_var(env)
    !
    !! get env vars from the WRF input file
    call allocate_wrf_var(wrf)
    call initialize_wrf_var(wrf)
    if(conf%MP_PHYSICS==40) then !read RAMS modified by oue to read RAMS 
      call ReadInpRAMS_var(trim(conf%WRFInputFile),trim(conf%WRFmpInputFile),wrf,status)
      call get_env_vars_rams(conf,wrf,env)
    else if ((conf%MP_PHYSICS==70) .or. (conf%MP_PHYSICS==75)) then!read SAM: modified by oue to read SAM
      if((conf%MP_PHYSICS==70) .and. (conf%InpProfile_flag/=2)) conf%InputProfile = "NULL"
      if((conf%MP_PHYSICS==75) .and. (conf%InpProfile_flag==0)) conf%InputProfile = "NULL"
      if((conf%MP_PHYSICS==75) .and. (conf%InpProfile_flag==1)) conf%InputProfile = conf%WRFmpInputFile
      call ReadInpSAM_var(trim(conf%WRFInputFile),trim(conf%InputProfile),wrf,status)
      call get_env_vars_sam(conf,wrf,env)
    else
      call ReadInpWRF_var(trim(conf%WRFInputFile),wrf,status)
      call get_env_vars(conf,wrf,env)
    endif
    call deallocate_wrf_var(wrf) 
    !
    ! at this point all needed env variables are extracted and provided on domain defined by user
    ! grid points correspond to wrf points where centre of the mass is defined
  
    ! now we read microphysical data
    ! conf%MP_PHYSICS==9
    if (conf%MP_PHYSICS==9) then ! milbrandt-yau 2-moment scheme
      WRFmpInputFile=conf%WRFInputFile
      mp09%nx=nx
      mp09%ny=ny
      mp09%nz=nz
      mp09%nt=nt
      !
      call allocate_wrf_var_mp09(mp09)
      call initialize_wrf_var_mp09(mp09)
      call ReadInpWRF_MP_PHYSICS_09(Trim(WRFmpInputFile),mp09,status)
      !
      write(*,*) '-----------------'
      write(*,*) '---mix ratio---'
      write(*,*) 'cloud',minval(mp09%qcloud),maxval(mp09%qcloud)
      write(*,*) 'rain',minval(mp09%qrain),maxval(mp09%qrain)
      write(*,*) 'ice',minval(mp09%qice),maxval(mp09%qice)
      write(*,*) 'snow',minval(mp09%qsnow),maxval(mp09%qsnow)
      write(*,*) 'graupel',minval(mp09%qgraup),maxval(mp09%qgraup)
      write(*,*) 'hail',minval(mp09%qhail),maxval(mp09%qhail)
  
      write(*,*) '---concentr---'
      write(*,*) 'cloud',minval(mp09%qncloud),maxval(mp09%qncloud)
      write(*,*) 'rain',minval(mp09%qnrain),maxval(mp09%qnrain)
      write(*,*) 'ice',minval(mp09%qnice),maxval(mp09%qnice)
      write(*,*) 'snow',minval(mp09%qnsnow),maxval(mp09%qnsnow)
      write(*,*) 'graupel',minval(mp09%qngraup),maxval(mp09%qngraup)
      write(*,*) 'hail',minval(mp09%qnhail),maxval(mp09%qnhail)
  
      write(*,*) '-----------------'
  
      hydro%nx=env%nx
      hydro%ny=env%ny
      hydro%nz=env%nz
      hydro%nht=nsc
      !----------------------------------------
      !----------------------------------------
      call allocate_hydro_var(hydro)
      call initialize_hydro_var(hydro)
      !
      call get_hydro09_vars(conf,mp09,hydro)
      call deallocate_wrf_var_mp09(mp09)
      !
      write(*,*) '--------------------------------------'
      write(*,*) '--------------------------------------'
      write(*,*) 'cloud',minval(hydro%qhydro(:,:,:,1)),maxval(hydro%qhydro(:,:,:,1))
      write(*,*) 'rain ',minval(hydro%qhydro(:,:,:,2)),maxval(hydro%qhydro(:,:,:,2))
      write(*,*) 'ice  ',minval(hydro%qhydro(:,:,:,3)),maxval(hydro%qhydro(:,:,:,3))
      write(*,*) 'snow ',minval(hydro%qhydro(:,:,:,4)),maxval(hydro%qhydro(:,:,:,4))
      write(*,*) 'graup',minval(hydro%qhydro(:,:,:,5)),maxval(hydro%qhydro(:,:,:,5))
      write(*,*) 'hail',minval(hydro%qhydro(:,:,:,6)),maxval(hydro%qhydro(:,:,:,6))
      write(*,*) '--------------------------------------'
      !
  
      !write(*,*) 'cloud',minval(hydro%qhydro(:,:,:,1), mask=hydro%qhydro(:,:,:,1)>0.d0),maxval(hydro%qhydro(:,:,:,1))
      !write(*,*) 'rain ',minval(hydro%qhydro(:,:,:,2), mask=hydro%qhydro(:,:,:,2)>0.d0),maxval(hydro%qhydro(:,:,:,2))
      !write(*,*) 'ice  ',minval(hydro%qhydro(:,:,:,3),mask=hydro%qhydro(:,:,:,3)>0.d0),maxval(hydro%qhydro(:,:,:,3))
      !write(*,*) 'snow ',minval(hydro%qhydro(:,:,:,4),mask=hydro%qhydro(:,:,:,4)>0.d0),maxval(hydro%qhydro(:,:,:,4))
      !write(*,*) 'graup',minval(hydro%qhydro(:,:,:,5),mask=hydro%qhydro(:,:,:,5)>0.d0),maxval(hydro%qhydro(:,:,:,5))
      !write(*,*) 'hail',minval(hydro%qhydro(:,:,:,6),mask=hydro%qhydro(:,:,:,6)>0.d0),maxval(hydro%qhydro(:,:,:,6))
      !
      !write(*,*) '--------------------------------------'
      !write(*,*) '--------------------------------------'
      write(*,*) 'cloud ',minval(hydro%qnhydro(:,:,:,1)),maxval(hydro%qnhydro(:,:,:,1))
      write(*,*) 'rain ',minval(hydro%qnhydro(:,:,:,2)),maxval(hydro%qnhydro(:,:,:,2))
      write(*,*) 'ice  ',minval(hydro%qnhydro(:,:,:,3)),maxval(hydro%qnhydro(:,:,:,3))
      write(*,*) 'snow ',minval(hydro%qnhydro(:,:,:,4)),maxval(hydro%qnhydro(:,:,:,4))
      write(*,*) 'graup',minval(hydro%qnhydro(:,:,:,5)),maxval(hydro%qnhydro(:,:,:,5))
      write(*,*) 'hail',minval(hydro%qnhydro(:,:,:,6)),maxval(hydro%qnhydro(:,:,:,6))
      !
      !write(*,*) '--------------------------------------'
      !write(*,*)'cloud',minval(hydro%qnhydro(:,:,:,1),mask=hydro%qnhydro(:,:,:,1)>0.d0),maxval(hydro%qnhydro(:,:,:,1))
      !write(*,*) 'rain',minval(hydro%qnhydro(:,:,:,2),mask=hydro%qnhydro(:,:,:,2)>0.d0),maxval(hydro%qnhydro(:,:,:,2))
      !write(*,*) 'ice',minval(hydro%qnhydro(:,:,:,3),mask=hydro%qnhydro(:,:,:,3)>0.d0),maxval(hydro%qnhydro(:,:,:,3))
      !write(*,*) 'snow',minval(hydro%qnhydro(:,:,:,4),mask=hydro%qnhydro(:,:,:,4)>0.d0),maxval(hydro%qnhydro(:,:,:,4))
      !write(*,*) 'graup',minval(hydro%qnhydro(:,:,:,5),mask=hydro%qnhydro(:,:,:,5)>0.d0),maxval(hydro%qnhydro(:,:,:,5))
      !write(*,*) 'hail',minval(hydro%qnhydro(:,:,:,6),mask=hydro%qnhydro(:,:,:,6)>0.d0),maxval(hydro%qnhydro(:,:,:,6))
  
      !
      !write(*,*) '--------------------------------------'
      !write(*,*) '--------------------------------------'
      !!pause
  
    endif
    !-------------------------------------------------------------------------------------------------------------------
    !-------------------------------------------------------------------------------------------------------------------
    !
    ! conf%MP_PHYSICS==10
    if (conf%MP_PHYSICS==10) then ! morrison 2-moment scheme
      WRFmpInputFile=conf%WRFInputFile
      mp10%nx=nx
      mp10%ny=ny
      mp10%nz=nz
      mp10%nt=nt
      !
      call allocate_wrf_var_mp10(mp10)
      call initialize_wrf_var_mp10(mp10)
      call ReadInpWRF_MP_PHYSICS_10(Trim(WRFmpInputFile),mp10,status)
      !
      hydro%nx=env%nx
      hydro%ny=env%ny
      hydro%nz=env%nz
      hydro%nht=nsc  
      !----------------------------------------
      call allocate_hydro_var(hydro)
      call initialize_hydro_var(hydro)
      !
      call get_hydro10_vars(conf,mp10,hydro)
      call deallocate_wrf_var_mp10(mp10)
      !
      write(*,*) '--------------------------------------'
      write(*,*) '---mix ratio--------------------------'
      write(*,*) 'cloud',minval(hydro%qhydro(:,:,:,1)),maxval(hydro%qhydro(:,:,:,1))
      write(*,*) 'rain ',minval(hydro%qhydro(:,:,:,2)),maxval(hydro%qhydro(:,:,:,2))
      write(*,*) 'ice  ',minval(hydro%qhydro(:,:,:,3)),maxval(hydro%qhydro(:,:,:,3))
      write(*,*) 'snow ',minval(hydro%qhydro(:,:,:,4)),maxval(hydro%qhydro(:,:,:,4))
      write(*,*) 'graup',minval(hydro%qhydro(:,:,:,5)),maxval(hydro%qhydro(:,:,:,5))
      write(*,*) '--------------------------------------'
      !
  
      !write(*,*) 'cloud',minval(hydro%qhydro(:,:,:,1), mask=hydro%qhydro(:,:,:,1)>0.d0),maxval(hydro%qhydro(:,:,:,1))
      !write(*,*) 'rain ',minval(hydro%qhydro(:,:,:,2), mask=hydro%qhydro(:,:,:,2)>0.d0),maxval(hydro%qhydro(:,:,:,2))
      !write(*,*) 'ice  ',minval(hydro%qhydro(:,:,:,3), mask=hydro%qhydro(:,:,:,3)>0.d0),maxval(hydro%qhydro(:,:,:,3))
      !write(*,*) 'snow ',minval(hydro%qhydro(:,:,:,4), mask=hydro%qhydro(:,:,:,4)>0.d0),maxval(hydro%qhydro(:,:,:,4))
      !write(*,*) 'graup',minval(hydro%qhydro(:,:,:,5), mask=hydro%qhydro(:,:,:,5)>0.d0),maxval(hydro%qhydro(:,:,:,5))
      !
      !write(*,*) '--------------------------------------'
      write(*,*) '---concentration-----------------------'
      write(*,*) 'rain ',minval(hydro%qnhydro(:,:,:,2)),maxval(hydro%qnhydro(:,:,:,2))
      write(*,*) 'ice  ',minval(hydro%qnhydro(:,:,:,3)),maxval(hydro%qnhydro(:,:,:,3))
      write(*,*) 'snow ',minval(hydro%qnhydro(:,:,:,4)),maxval(hydro%qnhydro(:,:,:,4))
      write(*,*) 'graup',minval(hydro%qnhydro(:,:,:,5)),maxval(hydro%qnhydro(:,:,:,5))
      !
      !write(*,*) '--------------------------------------'
      
      !write(*,*) 'rain ',minval(hydro%qnhydro(:,:,:,2),mask=hydro%qnhydro(:,:,:,2)>0.d0),maxval(hydro%qnhydro(:,:,:,2))
      !write(*,*) 'ice  ',minval(hydro%qnhydro(:,:,:,3),mask=hydro%qnhydro(:,:,:,3)>0.d0),maxval(hydro%qnhydro(:,:,:,3))
      !write(*,*) 'snow ',minval(hydro%qnhydro(:,:,:,4),mask=hydro%qnhydro(:,:,:,4)>0.d0),maxval(hydro%qnhydro(:,:,:,4))
      !write(*,*) 'graup',minval(hydro%qnhydro(:,:,:,5),mask=hydro%qnhydro(:,:,:,5)>0.d0),maxval(hydro%qnhydro(:,:,:,5))
      !
      write(*,*) '--------------------------------------'
      write(*,*) '--------------------------------------'
      
      !pause
    endif
    !
    !-------------------------------------------------------------------------------------------------------------------
    !!! Added by oue, 2016/09/19
    ! conf%MP_PHYSICS==8
    if (conf%MP_PHYSICS==8) then ! Thompson microphysics
      WRFmpInputFile=conf%WRFInputFile
      mp08%nx=nx
      mp08%ny=ny
      mp08%nz=nz
      mp08%nt=nt
      !
      call allocate_wrf_var_mp08(mp08) !from wrf_var_mod.f90
      call initialize_wrf_var_mp08(mp08) !from wrf_var_mod.f90
      call ReadInpWRF_MP_PHYSICS_08(Trim(WRFmpInputFile),mp08,status) ! from ReadInpWRFFile.f90
      !
      hydro%nx=env%nx
      hydro%ny=env%ny
      hydro%nz=env%nz
      hydro%nht=nsc  
      !----------------------------------------
      call allocate_hydro_var(hydro) !from wrf_var_mod.f90
      call initialize_hydro_var(hydro) !from wrf_var_mod.f90
      !
      call get_hydro08_vars(conf,mp08,hydro) !from crsim_subrs.f90
      call deallocate_wrf_var_mp08(mp08) !from wrf_var_mod.f90
      !
      write(*,*) '--------------------------------------'
      write(*,*) '---mix ratio--------------------------'
      write(*,*) 'cloud',minval(hydro%qhydro(:,:,:,1)),maxval(hydro%qhydro(:,:,:,1))
      write(*,*) 'rain ',minval(hydro%qhydro(:,:,:,2)),maxval(hydro%qhydro(:,:,:,2))
      write(*,*) 'ice  ',minval(hydro%qhydro(:,:,:,3)),maxval(hydro%qhydro(:,:,:,3))
      write(*,*) 'snow ',minval(hydro%qhydro(:,:,:,4)),maxval(hydro%qhydro(:,:,:,4))
      write(*,*) 'graup',minval(hydro%qhydro(:,:,:,5)),maxval(hydro%qhydro(:,:,:,5))
      write(*,*) '--------------------------------------'
      !
      
      write(*,*) '---concentration----------------------'
      write(*,*) 'rain ',minval(hydro%qnhydro(:,:,:,2)),maxval(hydro%qnhydro(:,:,:,2))
      write(*,*) 'ice  ',minval(hydro%qnhydro(:,:,:,3)),maxval(hydro%qnhydro(:,:,:,3))
      write(*,*) 'snow ',minval(hydro%qnhydro(:,:,:,4)),maxval(hydro%qnhydro(:,:,:,4))
      write(*,*) 'graup',minval(hydro%qnhydro(:,:,:,5)),maxval(hydro%qnhydro(:,:,:,5))
      
      write(*,*) '--------------------------------------'
      write(*,*) '--------------------------------------'
      
    endif
    !!!-- added by oue
    !------------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------------
    !!! Added by oue, 2017/07/17 for ICON
    ! conf%MP_PHYSICS==30
    if (conf%MP_PHYSICS==30) then ! ICON 2 moment
      WRFmpInputFile=conf%WRFInputFile
      mp30%nx=nx
      mp30%ny=ny
      mp30%nz=nz
      mp30%nt=nt
      !
      call allocate_wrf_var_mp30(mp30) !from wrf_var_mod.f90
      call initialize_wrf_var_mp30(mp30) !from wrf_var_mod.f90
      call ReadInpWRF_MP_PHYSICS_30(Trim(WRFmpInputFile),mp30,status) ! from ReadInpWRFFile.f90
      !
      hydro%nx=env%nx
      hydro%ny=env%ny
      hydro%nz=env%nz
      hydro%nht=nsc  
      !----------------------------------------
      call allocate_hydro_var(hydro) !from wrf_var_mod.f90
      call initialize_hydro_var(hydro) !from wrf_var_mod.f90
      !
      call get_hydro30_vars(conf,mp30,hydro) !from crsim_subrs.f90
      call deallocate_wrf_var_mp30(mp30) !from wrf_var_mod.f90
      !
      write(*,*) '--------------------------------------'
      write(*,*) '---mixing ratio-----------------------'
      write(*,*) 'cloud',minval(hydro%qhydro(:,:,:,1)),maxval(hydro%qhydro(:,:,:,1))
      write(*,*) 'rain ',minval(hydro%qhydro(:,:,:,2)),maxval(hydro%qhydro(:,:,:,2))
      write(*,*) 'ice  ',minval(hydro%qhydro(:,:,:,3)),maxval(hydro%qhydro(:,:,:,3))
      write(*,*) 'snow ',minval(hydro%qhydro(:,:,:,4)),maxval(hydro%qhydro(:,:,:,4))
      write(*,*) 'graup',minval(hydro%qhydro(:,:,:,5)),maxval(hydro%qhydro(:,:,:,5))
      write(*,*) 'hail ',minval(hydro%qhydro(:,:,:,6)),maxval(hydro%qhydro(:,:,:,6))
      write(*,*) '--------------------------------------'
      write(*,*) '---concentration----------------------'
      write(*,*) 'cloud',minval(hydro%qnhydro(:,:,:,1)),maxval(hydro%qnhydro(:,:,:,1))
      write(*,*) 'rain ',minval(hydro%qnhydro(:,:,:,2)),maxval(hydro%qnhydro(:,:,:,2))
      write(*,*) 'ice  ',minval(hydro%qnhydro(:,:,:,3)),maxval(hydro%qnhydro(:,:,:,3))
      write(*,*) 'snow ',minval(hydro%qnhydro(:,:,:,4)),maxval(hydro%qnhydro(:,:,:,4))
      write(*,*) 'graup',minval(hydro%qnhydro(:,:,:,5)),maxval(hydro%qnhydro(:,:,:,5))
      write(*,*) 'hail ',minval(hydro%qnhydro(:,:,:,6)),maxval(hydro%qnhydro(:,:,:,6))
      write(*,*) '--------------------------------------'
      write(*,*) '--------------------------------------'
      
    endif
    !!!-- added by oue
    !---------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------------
    !!! Added by oue, 2017/07/21 for RAMS
    ! conf%MP_PHYSICS==40
    if (conf%MP_PHYSICS==40) then ! RAMS 2 moment
      WRFmpInputFile=conf%WRFInputFile
      mp40%nx=nx
      mp40%ny=ny
      mp40%nz=nz
      mp40%nt=nt
      !
      call allocate_wrf_var_mp40(mp40) !from wrf_var_mod.f90
      call initialize_wrf_var_mp40(mp40) !from wrf_var_mod.f90
      call ReadInpRAMS_MP_PHYSICS_40(Trim(WRFmpInputFile),mp40,status) ! from ReadInpRAMSFile.f90
      !
      hydro%nx=env%nx
      hydro%ny=env%ny
      hydro%nz=env%nz
      hydro%nht=nsc  
      !----------------------------------------
      call allocate_hydro_var(hydro) !from wrf_var_mod.f90
      call initialize_hydro_var(hydro) !from wrf_var_mod.f90
      !
      call get_hydro40_vars(conf,mp40,hydro) !from crsim_surs.f90
      call deallocate_wrf_var_mp40(mp40) !from wrf_var_mod.f90
      !
      write(*,*) '--------------------------------------'
      write(*,*) '---mixing ratio-----------------------'
      write(*,*) 'cloud',minval(hydro%qhydro(:,:,:,1)),maxval(hydro%qhydro(:,:,:,1))
      write(*,*) 'rain ',minval(hydro%qhydro(:,:,:,2)),maxval(hydro%qhydro(:,:,:,2))
      write(*,*) 'ice  ',minval(hydro%qhydro(:,:,:,3)),maxval(hydro%qhydro(:,:,:,3))
      write(*,*) 'snow ',minval(hydro%qhydro(:,:,:,4)),maxval(hydro%qhydro(:,:,:,4))
      write(*,*) 'graup',minval(hydro%qhydro(:,:,:,5)),maxval(hydro%qhydro(:,:,:,5))
      write(*,*) 'hail ',minval(hydro%qhydro(:,:,:,6)),maxval(hydro%qhydro(:,:,:,6))
      write(*,*) 'drizzle',minval(hydro%qhydro(:,:,:,7)),maxval(hydro%qhydro(:,:,:,7))
      write(*,*) 'aggregate',minval(hydro%qhydro(:,:,:,8)),maxval(hydro%qhydro(:,:,:,8))
      write(*,*) '--------------------------------------'
      write(*,*) '---concentration----------------------'
      write(*,*) 'cloud',minval(hydro%qnhydro(:,:,:,1)),maxval(hydro%qnhydro(:,:,:,1))
      write(*,*) 'rain ',minval(hydro%qnhydro(:,:,:,2)),maxval(hydro%qnhydro(:,:,:,2))
      write(*,*) 'ice  ',minval(hydro%qnhydro(:,:,:,3)),maxval(hydro%qnhydro(:,:,:,3))
      write(*,*) 'snow ',minval(hydro%qnhydro(:,:,:,4)),maxval(hydro%qnhydro(:,:,:,4))
      write(*,*) 'graup',minval(hydro%qnhydro(:,:,:,5)),maxval(hydro%qnhydro(:,:,:,5))
      write(*,*) 'hail ',minval(hydro%qnhydro(:,:,:,6)),maxval(hydro%qnhydro(:,:,:,6))
      write(*,*) 'drizzle',minval(hydro%qnhydro(:,:,:,7)),maxval(hydro%qnhydro(:,:,:,7))
      write(*,*) 'aggregate',minval(hydro%qnhydro(:,:,:,8)),maxval(hydro%qnhydro(:,:,:,8))
      write(*,*) '--------------------------------------'
      write(*,*) '--------------------------------------'
      
    endif
    !!!-- added by oue
    !---------------------------------------------------------------------------------
    
  
    !!-----------------------------------------------------------
    ! | if (conf%MP_PHYSICS==50) then | added by DW 2017/10/30 for P3
    !
    if (conf%MP_PHYSICS==50) then ! added by DW
      WRFmpInputFile=conf%WRFInputFile ! added by DW
      mp50%nx=nx ! added by DW
      mp50%ny=ny ! added by DW
      mp50%nz=nz ! added by DW
      mp50%nt=nt ! added by DW
      !
      call allocate_wrf_var_mp50(mp50) ! added by DW
      call initialize_wrf_var_mp50(mp50)  ! added by DW
      call ReadInpWRF_MP_PHYSICS_50(Trim(WRFmpInputFile),mp50,status) ! added by DW
      !
      hydro50%nx=env%nx ! added by DW
      hydro50%ny=env%ny ! added by DW
      hydro50%nz=env%nz ! added by DW
      hydro50%nht=nsc  ! added by DW
      !
      call allocate_hydro50_var(hydro50) ! added by DW
      call initialize_hydro50_var(hydro50) ! added by DW
      call get_hydro50_vars(conf,mp50,hydro50)     ! added by DW
      !
      call deallocate_wrf_var_mp50(mp50)! added by DW
    endif  ! added by DW
    !
    !--------------------------------------------------------------
    
    !-------------------------------------------------------------------------------------------------------------------
    ! SAM Morrison microphysics added by oue 2018/06
    ! conf%MP_PHYSICS==75
    if (conf%MP_PHYSICS==75) then !sam morrison 2-moment scheme
      WRFmpInputFile=conf%WRFInputFile
      mp75%nx=nx
      mp75%ny=ny
      mp75%nz=nz
      mp75%nt=nt
      !
      call allocate_wrf_var_mp75(mp75)
      call initialize_wrf_var_mp75(mp75)
      write(*,*)'ReadInpSAM_MP_PHYSICS_75'
      call ReadInpSAM_MP_PHYSICS_75(Trim(WRFmpInputFile),mp75,status)
      !
      hydro%nx=env%nx
      hydro%ny=env%ny
      hydro%nz=env%nz
      hydro%nht=nsc  
      !----------------------------------------
      call allocate_hydro_var(hydro)
      call initialize_hydro_var(hydro)
      !
      call get_hydro75_vars(conf,env,mp75,hydro)
      call deallocate_wrf_var_mp75(mp75)
      !
      write(*,*) '--------------------------------------'
      write(*,*) '---mix ratio--------------------------'
      write(*,*) 'cloud',minval(hydro%qhydro(:,:,:,1)),maxval(hydro%qhydro(:,:,:,1))
      write(*,*) 'rain ',minval(hydro%qhydro(:,:,:,2)),maxval(hydro%qhydro(:,:,:,2))
      write(*,*) 'ice  ',minval(hydro%qhydro(:,:,:,3)),maxval(hydro%qhydro(:,:,:,3))
      write(*,*) 'snow ',minval(hydro%qhydro(:,:,:,4)),maxval(hydro%qhydro(:,:,:,4))
      write(*,*) 'graup',minval(hydro%qhydro(:,:,:,5)),maxval(hydro%qhydro(:,:,:,5))
      write(*,*) '--------------------------------------'
      !
      write(*,*) '---concentration-----------------------'
      write(*,*) 'cloud',minval(hydro%qnhydro(:,:,:,1)),maxval(hydro%qnhydro(:,:,:,1))
      write(*,*) 'rain ',minval(hydro%qnhydro(:,:,:,2)),maxval(hydro%qnhydro(:,:,:,2))
      write(*,*) 'ice  ',minval(hydro%qnhydro(:,:,:,3)),maxval(hydro%qnhydro(:,:,:,3))
      write(*,*) 'snow ',minval(hydro%qnhydro(:,:,:,4)),maxval(hydro%qnhydro(:,:,:,4))
      write(*,*) 'graup',minval(hydro%qnhydro(:,:,:,5)),maxval(hydro%qnhydro(:,:,:,5))
      write(*,*) '--------------------------------------'
      write(*,*) '--------------------------------------'
      !pause
    endif
    !
    !--------------------------------------------------------------
    
    
    !--------------------------------------------------------------
    !conf%MP_PHYSICS==20
    if (conf%MP_PHYSICS==20) then ! bin-explicit Jiwen Fan
      !--------------------------------
    
      Allocate(scatt_type(nsc))
      do isc=1,nsc
      !
        call hydro_info_jfan(isc,nkr,scatt_type(isc)%ihtf,&
                             scatt_type(isc)%id1,scatt_type(isc)%id2,scatt_type(isc)%nbins)
        ! 
        Allocate(scatt_type(isc)%diam(scatt_type(isc)%nbins))
        Allocate(scatt_type(isc)%diam1(scatt_type(isc)%nbins+1)) ! used for spectra simulation
        Allocate(scatt_type(isc)%rho(scatt_type(isc)%nbins))
        Allocate(scatt_type(isc)%mass(scatt_type(isc)%nbins))
        Allocate(scatt_type(isc)%fvel(scatt_type(isc)%nbins))
        !
        scatt_type(isc)%diam1(:)=-999.d0
      enddo
      ! get diameters, densities,mass and fall velocities
      call hydro_jfan(conf,nkr,nsc,scatt_type)
      !
      !---------------------------------
      ! reaf the bin distribution netcf files
      WRFmpInputFile=conf%WRFmpInputFile
      mp20%nx=nx
      mp20%ny=ny
      mp20%nz=nz
      mp20%nt=nt
      mp20%nbins=nkr ! (33)
      !
      call allocate_wrf_var_mp20(mp20)
      call initialize_wrf_var_mp20(mp20)
      write(*,*) 'Reading: ',Trim(WRFmpInputFile)
      call ReadInpWRF_MP_PHYSICS_20(Trim(WRFmpInputFile),mp20,status)
      !
      do ii=1,nt
        write(*,*) '-----------------------'
        write(*,*) 'time step=',ii
        write(*,*) 'ff1',minval(mp20%ff1(:,:,:,ii,:),mask=mp20%ff1(:,:,:,ii,:)>0.e0),maxval(mp20%ff1(:,:,:,ii,:))
        write(*,*) 'ff5',minval(mp20%ff5(:,:,:,ii,:),mask=mp20%ff5(:,:,:,ii,:)>0.e0),maxval(mp20%ff5(:,:,:,ii,:))
        write(*,*) 'ff6',minval(mp20%ff6(:,:,:,ii,:),mask=mp20%ff6(:,:,:,ii,:)>0.e0),maxval(mp20%ff6(:,:,:,ii,:))
        write(*,*) '-----------------------'
      enddo
      !
      hydro20%nx=env%nx
      hydro20%ny=env%ny
      hydro20%nz=env%nz
      hydro20%nbins=mp20%nbins
      !----------------------------------------
      call allocate_hydro20_var(hydro20)
      call initialize_hydro20_var(hydro20)
      !
      call get_hydro20_vars(conf,mp20,hydro20)
      !
      call deallocate_wrf_var_mp20(mp20)
      !
      write(*,*) '--------------------------------------'
      write(*,*) '--------------------------------------'
      write(*,*) 'ff1',minval(hydro20%ff1(:,:,:,:),mask=hydro20%ff1>0.d0),maxval(hydro20%ff1(:,:,:,:))
      write(*,*) 'ff5',minval(hydro20%ff5(:,:,:,:),mask=hydro20%ff5>0.d0),maxval(hydro20%ff5(:,:,:,:))
      write(*,*) 'ff6',minval(hydro20%ff6(:,:,:,:),mask=hydro20%ff6>0.d0),maxval(hydro20%ff6(:,:,:,:))
      
      write(*,*) '--------------------------------------'
      !
      ! now split out ff1 to cloud and rain, and ff5 to ice and snow
      !
      do isc=1,nsc
      !
        scatt_type(isc)%nx=hydro20%nx
        scatt_type(isc)%ny=hydro20%ny
        scatt_type(isc)%nz=hydro20%nz
        !
        Allocate(scatt_type(isc)%qq(scatt_type(isc)%nx,scatt_type(isc)%ny,scatt_type(isc)%nz))
        Allocate(scatt_type(isc)%N(scatt_type(isc)%nx,scatt_type(isc)%ny,scatt_type(isc)%nz,scatt_type(isc)%nbins))
        !
        scatt_type(isc)%qq=0.d0
        scatt_type(isc)%N=0.d0
        !
        !
        if ((isc==1).or.(isc==2)) then
          !write(*,*) isc,scatt_type(isc)%id1,scatt_type(isc)%id2
          scatt_type(isc)%N(:,:,:,:)=hydro20%ff1(:,:,:,scatt_type(isc)%id1:scatt_type(isc)%id2)
          scatt_type(isc)%qq=Sum(scatt_type(isc)%N,dim=4)  ! total mixing ratio per hydrom type in kg/kg = Sum ( bin mass mix. ratio[kg/kg] )
        endif
        !
        if ((isc==3).or.(isc==4)) then
          !write(*,*) isc,scatt_type(isc)%id1,scatt_type(isc)%id2
          scatt_type(isc)%N(:,:,:,:)=hydro20%ff5(:,:,:,scatt_type(isc)%id1:scatt_type(isc)%id2)
          scatt_type(isc)%qq=Sum(scatt_type(isc)%N,dim=4)  ! total mixing ratio per hydrom type in kg/kg = Sum ( bin mass mix. ratio[kg/kg] )   
        endif
        !
        if (isc==5) then
          !write(*,*) isc,scatt_type(isc)%id1,scatt_type(isc)%id2
          scatt_type(isc)%N(:,:,:,:)=hydro20%ff6(:,:,:,scatt_type(isc)%id1:scatt_type(isc)%id2)
          scatt_type(isc)%qq=Sum(scatt_type(isc)%N,dim=4)  ! total mixing ratio per hydrom type in kg/kg = Sum ( bin mass mix. ratio[kg/kg] )
        endif
        !
      enddo ! isc
      !
      call deallocate_hydro20_var(hydro20)
    !!
    endif ! if (conf%MP_PHYSICS==20) 
    !--------------------------------------------------------------
  
  
  
    !--------------------------------------------------------------
    !---SAM warm bin; added by oue----------------------------------
    !conf%MP_PHYSICS==70
    if (conf%MP_PHYSICS==70) then 
    !--------------------------------
      Allocate(scatt_type(nsc))
      do isc=1,nsc
      !
        call hydro_info_samsbm(isc,nkr,scatt_type(isc)%ihtf,&
                               scatt_type(isc)%id1,scatt_type(isc)%id2,scatt_type(isc)%nbins)
        Allocate(scatt_type(isc)%mass(scatt_type(isc)%nbins))
        Allocate(scatt_type(isc)%diam1(scatt_type(isc)%nbins + 1))
        Allocate(scatt_type(isc)%mass1(scatt_type(isc)%nbins + 1))
      !
      enddo
      ! get diameters, densities,mass and fall velocities
      WRFmpInputFile=conf%WRFmpInputFile
      write(*,*) 'Reading: ',Trim(WRFmpInputFile)
      call hydro_samsbm(conf,nsc,nkr,scatt_type)
      !
      !---------------------------------
      ! read the bin distribution netcf files
      mp70%nx=nx
      mp70%ny=ny
      mp70%nz=nz
      mp70%nt=nt
      mp70%nbins=nkr ! (36)
      !
      call allocate_wrf_var_mp70(mp70)
      call initialize_wrf_var_mp70(mp70)
      call ReadInpSAM_MP_PHYSICS_70(Trim(conf%WRFInputFile),nkr,mp70,status)
      !
      do ii=1,nt
        write(*,*) '-----------------------'
        write(*,*) 'time step=',ii, ',within entire domain'
        write(*,*) 'Mass(g/kg)',minval(mp70%fm1(:,:,:,ii,:),mask=mp70%fm1(:,:,:,ii,:)>0.d0),maxval(mp70%fm1(:,:,:,ii,:))
        write(*,*) 'Number(/cm3)',minval(mp70%fn1(:,:,:,ii,:),mask=mp70%fn1(:,:,:,ii,:)>0.d0),maxval(mp70%fn1(:,:,:,ii,:))
        write(*,*) '-----------------------'
      enddo
      !
      hydro70%nx=env%nx
      hydro70%ny=env%ny
      hydro70%nz=env%nz
      hydro70%nbins=mp70%nbins
      hydro70%nht=nsc
      !----------------------------------------
      call allocate_hydro70_var(hydro70)
      call initialize_hydro70_var(hydro70)
      !
      call get_hydro70_vars(conf,mp70,hydro70)
      !
      call deallocate_wrf_var_mp70(mp70)
      !
      write(*,*) '--------------------------------------'
      write(*,*) 'Within selected domain'
      write(*,*) 'Mass(kg/kg)',minval(hydro70%fm1(:,:,:,:),mask=hydro70%fm1>0.d0),maxval(hydro70%fm1(:,:,:,:))
      write(*,*) 'Number(/m3)',minval(hydro70%fn1(:,:,:,:),mask=hydro70%fn1>0.d0),maxval(hydro70%fn1(:,:,:,:))
      write(*,*) '--------------------------------------'
      !
      ! now split out fn1 to cloud and rain, 
      !
      do isc=1,nsc
        !
        scatt_type(isc)%nx=hydro70%nx
        scatt_type(isc)%ny=hydro70%ny
        scatt_type(isc)%nz=hydro70%nz
        !
        Allocate(scatt_type(isc)%qq(scatt_type(isc)%nx,scatt_type(isc)%ny,scatt_type(isc)%nz))
        Allocate(scatt_type(isc)%N(scatt_type(isc)%nx,scatt_type(isc)%ny,scatt_type(isc)%nz,scatt_type(isc)%nbins))
        Allocate(scatt_type(isc)%NN(scatt_type(isc)%nx,scatt_type(isc)%ny,scatt_type(isc)%nz,scatt_type(isc)%nbins))
        Allocate(scatt_type(isc)%diam2(scatt_type(isc)%nx,scatt_type(isc)%ny,scatt_type(isc)%nz,scatt_type(isc)%nbins))
        Allocate(scatt_type(isc)%rho2(scatt_type(isc)%nx,scatt_type(isc)%ny,scatt_type(isc)%nz,scatt_type(isc)%nbins))
        Allocate(scatt_type(isc)%fvel2(scatt_type(isc)%nx,scatt_type(isc)%ny,scatt_type(isc)%nz,scatt_type(isc)%nbins))
        !
        scatt_type(isc)%qq=0.d0
        scatt_type(isc)%N=0.d0
        scatt_type(isc)%NN=0.d0
        !
        if ((isc==1).or.(isc==2)) then
          scatt_type(isc)%N(:,:,:,:)=hydro70%fm1(:,:,:,scatt_type(isc)%id1:scatt_type(isc)%id2)
          scatt_type(isc)%NN(:,:,:,:)=hydro70%fn1(:,:,:,scatt_type(isc)%id1:scatt_type(isc)%id2)
          if (isc==1) scatt_type(isc)%qq(:,:,:)=hydro70%qhydro(:,:,:,1)
          if (isc==2) scatt_type(isc)%qq(:,:,:)=hydro70%qhydro(:,:,:,2)
        endif
      !
      enddo ! isc
      !
      call deallocate_hydro70_var(hydro70)
      !!
    endif ! if (conf%MP_PHYSICS==70) 
  
    !-------------------------------------------------------------
  
  
    !============================================================
    !-- Initialize for Doppler spectra simulation
    call initialize_spectra_var(spectra,conf%freq,conf%ZMIN) ! parameter settings 
    if(conf%spectraID==1)then
      spectra%nx=env%nx
      spectra%ny=env%ny
      spectra%nz=env%nz
      spectra%nht=nsc
      call allocate_spectra_var(spectra,conf%radID)
      nfft = spectra%NFFT
    else
      nfft = 1
      spectra%NFFT=1
    end if
    !Allocate(spectra_bins(nfft),zhh_spectra(nfft),zvh_spectra(nfft),zvv_spectra(nfft))
    !spectra_bins=0.d0; zhh_spectra=0.d0; zvv_spectra=0.d0; zvh_spectra=0.d0
    !============================================================
    !-------------------------------------------------------------
    !
    ! define dimensions of the CR-SIM output (rmout) 
    ! and dimensions of other variables needed (mout) 
    !
    !-----------------------------------------
    if ( (conf%MP_PHYSICS==9) .or. &
         (conf%MP_PHYSICS==10) .or. (conf%MP_PHYSICS==8) .or. &
         (conf%MP_PHYSICS==30) .or. (conf%MP_PHYSICS==40).or. & ! modified by oue 2016/09/19, 2017/07/17 ICON, RAMS
         (conf%MP_PHYSICS==75)) then ! modified by oue 2018/06 SAM
      rmout%nx=hydro%nx
      rmout%ny=hydro%ny
      rmout%nz=hydro%nz
      rmout%nht=hydro%nht
      !
      IF (conf%radID/=1) THEN ! if radID/=1 include polarimet. vars
        mout%nx=hydro%nx
        mout%ny=hydro%ny
        mout%nz=hydro%nz
        mout%nht=hydro%nht
      ENDIF
      !
    endif
    !-----------------------------------------
    !-----------------------------------------
    if ((conf%MP_PHYSICS==20).or.(conf%MP_PHYSICS==70)) then !added by oue for SAM warmbin
      rmout%nx=scatt_type(1)%nx
      rmout%ny=scatt_type(1)%ny
      rmout%nz=scatt_type(1)%nz
      rmout%nht=nsc
      
        IF (conf%radID/=1) THEN ! if radID/=1 include polarimet. vars
          mout%nx=scatt_type(1)%nx
          mout%ny=scatt_type(1)%ny
          mout%nz=scatt_type(1)%nz
          mout%nht=nsc
        ENDIF
    endif
    !-----------------------------------------
    !---------------------------------------------------
    !
    if (conf%MP_PHYSICS==50) then ! added by DW 2017/10/30 for P3
      rmout%nx=hydro50%nx  ! added by DW 
      rmout%ny=hydro50%ny  ! added by DW
      rmout%nz=hydro50%nz  ! added by DW
      rmout%nht=hydro50%nht ! added by DW
      !
        IF (conf%radID/=1) THEN ! added by DW
          mout%nx=hydro50%nx ! added by DW
          mout%ny=hydro50%ny ! added by DW
          mout%nz=hydro50%nz ! added by DW
          mout%nht=hydro50%nht ! added by DW
        ENDIF ! added by DW
      !
    endif ! added by DW
    !---------------------------------------------------
    !-----------------------------------------
    !-----------------------------------------
    if (conf%ceiloID==1) then ! include ceilo measurements
      lout%nx=rmout%nx
      lout%ny=rmout%ny
      lout%nz=rmout%nz
      lout%nht=1
      if(conf%MP_PHYSICS==40) lout%nht=2 !cloud,ice,drizzle
    endif
    !-----------------------------------------
    if (conf%mplID>0) then ! include MPL measurements
      mpl%nx=rmout%nx
      mpl%ny=rmout%ny
      mpl%nz=rmout%nz
      mpl%nht=2 ! cloud and ice
        if(conf%MP_PHYSICS==40) mpl%nht=3 !cloud,ice,drizzle
        mpl%wavel=conf%mpl_wavel
        if (conf%aeroID==1) then ! include aerosol profile
          aero%nx=rmout%nx
          aero%ny=rmout%ny
          aero%nz=rmout%nz
        endif
    endif
    !-----------------------------------------
    !== Dimensions for post processing
    if (conf%arsclID==1) then ! include ARSCL products added by oue 2017.03.23
      arscl%nx=rmout%nx
      arscl%ny=rmout%ny
      arscl%nz=rmout%nz
      arscl%n_layers=N_layers
    endif
    if (conf%mwrID==1) then ! include MWR LWP added by oue 2017.03.23
      mwr%nx=rmout%nx
      mwr%ny=rmout%ny
      mwr%nz=rmout%nz
    endif
    !==================================
    !-----------------------------------------
    !___________________________________________________________________________________________
    !
    ! Allocate rmout and mout vars
    call allocate_rmout_var(rmout)
    call m999_rmout_var(rmout)
    IF (conf%radID/=1) THEN
      call allocate_mout_var(mout)
      call m999_mout_var(mout)
    ENDIF
    !
    ! Allocate lout vars
    IF (conf%ceiloID==1) then
      call allocate_lout_var(lout)
      call m999_lout_var(lout)
    ENDIF
    !
    ! Allocate mpl vars
    IF (conf%mplID>0) then
      call allocate_mpl_var(mpl)
      call m999_mpl_var(mpl)
      IF (conf%aeroID==1) THEN
        call allocate_aero_var(aero)
        call m999_aero_var(aero)
      ENDIF
    ENDIF
    !
    !== Vars for post processing==
    ! Allocate arscl vars added by oue 2017.03.23
    IF (conf%arsclID==1) then
      call allocate_arscl_var(arscl)
      call m999_arscl_var(arscl)
      Allocate(mpl_back_hydro(arscl%nx,arscl%ny,arscl%nz),mpl_back_aero(arscl%nx,arscl%ny,arscl%nz)) 
    ENDIF
    !
    ! Allocate mwr vars added by oue 2017.03.23
    IF (conf%mwrID==1) then
      call allocate_mwr_var(mwr)
      call m999_mwr_var(mwr)
    ENDIF
    !==============================
    !
    DO isc=1, MAX_TH_NUM
      rmrad_vec(isc)%nht=rmout%nht
      call allocate_rmrad_var(rmrad_vec(isc))
      IF (conf%radID/=1) THEN
        mrad_vec(isc)%nht=mout%nht
        call allocate_mrad_var(mrad_vec(isc))
      ENDIF
    ENDDO
    !
    !
    !----------------------------------
    ! 
    !define average aerosol profile as in Spinhirne,1993
    IF (conf%mplID>0) THEN
      call compute_molecular_backscatter(conf,env,mpl)
      IF (conf%aeroID==1) THEN
        call spinhirne_aero_model(conf,env,aero)
      ENDIF
    ENDIF
    !----------------------------------
    ! elevation is fixed if conf%elev >= 0
    if (conf%elev >= 0.d0 ) rmout%elev=conf%elev

    ! elevation is fixed if  -90.d0 <= conf%elev <= +90  ! alex
    !if ( (conf%elev >= -90.d0 ) .and. (conf%elev <= 90.d0) ) then
    !  rmout%elev = conf%elev
    !endif
    !
    rmout%azim=-999.d0
  
  
    !wend = omp_get_wtime()
    !write(*,*) 'Part 1 Time : ',  wend - wstart
  
  
    !wstart = omp_get_wtime()
  
  
    ! preload LUT file data
  
#ifdef __PRELOAD_LUT__
  
      !$ wstart = omp_get_wtime()  ! AT April 2019
  
      size_lut_files = 0
      do ix=1,rmout%nx
        do iy=1,rmout%ny
          do iz=1,rmout%nz
            elev=conf%elev
            rr=-999.d0
  
  
            if (conf%elev < 0.d0) then
            ! AUG2019
            !if ( (conf%elev < -90.d0) .or. (conf%elev > 90.d0) )  then !AT 
              ! compute elevation and range of the pixel with coord x(ix),y(iy)
              ! and z(iz) relative to the radar origin
              call determine_elevation_and_range(env%x(ix),env%y(iy),env%z(ix,iy,iz),&
                                                 conf%ixc,conf%iyc,conf%zc,env%dx,env%dy,elev,rr)  
              rmout%elev(ix,iy,iz)=elev
              rmout%range(ix,iy,iz)=rr
            endif
  
            qq=0.d0
            do iht=1,nsc
            !write(*,*) '(ix,iy,iz,iht) : (',ix,',',iy,',',iz,',',iht,')'
              lut_search_flag = 0
              LutFileName = ''
  
              if (conf%MP_PHYSICS== 9) qq=hydro%qhydro(ix,iy,iz,iht)
              if (conf%MP_PHYSICS==10) qq=hydro%qhydro(ix,iy,iz,iht)
              if (conf%MP_PHYSICS==20) qq=scatt_type(iht)%qq(ix,iy,iz)
              if (conf%MP_PHYSICS== 8) qq=hydro%qhydro(ix,iy,iz,iht) !Added by oue 2016/09/19
              if (conf%MP_PHYSICS==30) qq=hydro%qhydro(ix,iy,iz,iht) !Added by oue 2017/07/17 ICON
              if (conf%MP_PHYSICS==40) qq=hydro%qhydro(ix,iy,iz,iht) !Added by oue 2017/07/21 RAMS
              if (conf%MP_PHYSICS==50) qq=hydro50%qhydro(ix,iy,iz,iht) ! P3 
              if (conf%MP_PHYSICS==70) qq=scatt_type(iht)%qq(ix,iy,iz) ! for SAM warm bin
              if (conf%MP_PHYSICS==75) qq=hydro%qhydro(ix,iy,iz,iht)   ! for SAM morrison
  
              nden=-1
              Allocate(lden(1))
              lden(1)=-1.d0
              !--Added by oue for preloading LUTs of all densities
              if((iht==1) .or. (iht==7)) then
                if ( allocated(lden)) deallocate(lden)
                nden=n_lden_cld
                Allocate(lden(nden))
                lden=lden_cld  
              endif
              if(iht==2) then
                if ( allocated(lden)) deallocate(lden)
                nden=n_lden_rain
                Allocate(lden(nden))
                lden=lden_rain  
              endif
              if(conf%MP_PHYSICS==50) then !oue  P3
                 if(iht==3) then
                   if ( allocated(lden)) deallocate(lden)
                   nden=n_lden_smallice
                   Allocate(lden(nden))
                   lden=lden_smallice  
                 endif
                 if(iht==4) then
                    if ( allocated(lden)) deallocate(lden)
                    nden=n_lden_unrice
                    Allocate(lden(nden))
                    lden=lden_unrice 
                 endif
                 if(iht==5) then
                    if ( allocated(lden)) deallocate(lden)
                    nden=n_lden_graupP3
                    Allocate(lden(nden))
                    lden=lden_graupP3  
                 endif
                 if (iht==6) then
                    if ( allocated(lden)) deallocate(lden)
                    nden=n_lden_parice
                    Allocate(lden(nden))
                    lden=lden_parice  
                 endif
              else !oue  P3
                if(iht==3) then
                  if ( allocated(lden)) deallocate(lden)
                  nden=n_lden_ice
                  Allocate(lden(nden))
                  lden=lden_ice 
                endif
                if((iht==4) .or. (iht==8)) then
                  if ( allocated(lden)) deallocate(lden)
                  nden=n_lden_snow
                  Allocate(lden(nden))
                  lden=lden_snow 
                endif
                if((iht==5).or. (iht==6)) then
                  if ( allocated(lden)) deallocate(lden)
                  nden=n_lden_graup
                  Allocate(lden(nden))
                  lden=lden_graup  
                endif
              endif !-oue  P3
    
              If (nden==-1) then
                write(*,*) 'Problem with nden'
                exit
              Endif
              If (lden(1) < 0.d0) then
                write(*,*) 'Problem with lden'
                exit
              Endif
    
              !-oue
              do iden=1,nden !--Added do loop by oue for preloading LUTs of all densities
                if (qq > 0.d0) then
                  ! AUG2019
                  elevx = elev ! radar looking up 
                  !if (conf%airborne /= 0) elevx = MAX(0.d0, 90.d0 - dabs(elev))  ! radar looking down AT 
                  call make_lutfilename(iht, conf, elevx,env%temp(ix,iy,iz), lden(iden), LutFileName) !for preloading LUTs of all densities
                  !  
                  search: do ii=0,size_lut_files
                    !write(*,*) 'ii : ', ii, ', size_lut_files : ', size_lut_files, ', lut_search_flag', lut_search_flag,  &
                    !    ', String cmp : ', (Trim(LutFileName) .EQ. Trim(lut_files(ii+1)%filenameLUTS))
                    !write(*,*) 'LutFileName : ', Trim(LutFileName)
                    !write(*,*) 'lut_files(ii+1)%filenameLUTS : ', Trim(lut_files(ii+1)%filenameLUTS)
                    if( Trim(LutFileName) .EQ. Trim(lut_files(ii+1)%filenameLUTS) ) then
                      lut_search_flag = lut_search_flag + 1
                      exit search
                    endif
                  enddo search
                  ! 
                  if(lut_search_flag == 0) then
                    call preload_luts(Trim(LutFileName), lut_files(size_lut_files+1))
                    !write(*,*) size_lut_files + 1, lut_files(size_lut_files+1)%filenameLUTS
                    size_lut_files = size_lut_files + 1
                  endif
                  !
                endif
              enddo !-oue
              Deallocate(lden) !-oue
    
    
            enddo
          enddo
        enddo 
      enddo
    
    
      write(*,*) 'loaded lut files : ', size_lut_files
      !$ wend = omp_get_wtime()  ! AT April 2019
      !$ write(*,*) 'LUT files preloading time : ', wend - wstart ! AT April 2019
    
      if(conf%MP_PHYSICS==70) then !added by our for SAM warm bin
        call make_mpl_lutfilename(1, conf, LutFileName)
        call preload_mpl_luts(LutFileName, mpl_lut_files(1))
        write(*,*) 'Loaded MPL LUT File : ', Trim(LutFileName)
      else
        call make_mpl_lutfilename(1, conf, LutFileName)
        call preload_mpl_luts(LutFileName, mpl_lut_files(1))
        write(*,*) 'Loaded MPL LUT File : ', Trim(LutFileName)
        call make_mpl_lutfilename(3, conf, LutFileName)
        call preload_mpl_luts(LutFileName, mpl_lut_files(3))
        write(*,*) 'Loaded MPL LUT File : ', Trim(LutFileName)
     
        if(hydro%nht > 6) then !for RAMS
          call make_mpl_lutfilename(7, conf, LutFileName)
          call preload_mpl_luts(LutFileName, mpl_lut_files(7))
          write(*,*) 'Loaded MPL LUT File : ', Trim(LutFileName)
        endif
      endif
  
      call preload_lut_ceilo(conf)
  
  
      flush(6)
#endif
  
  
  
    Allocate(mpl_back_obs(mpl%nz))
    !
  
    psum_xy=0
    psum_xy_tot=rmout%nx*rmout%nx
    !
    ! MAIN LOOP
    !
    ihtc=-1
    percentage_saved=0
  
#ifdef __PRELOAD_LUT__
      !$OMP PARALLEL DO DEFAULT(SHARED) &
      !$OMP& FIRSTPRIVATE(elev) &
      !$OMP& FIRSTPRIVATE(elevx) &
      !$OMP& PRIVATE(ix,iy,iz,tid,rr,azim,uu,vv,ww,tke,qsum,isc,sw_t2,sw_s2,sw_v2,iht,qq) &
      !$OMP& PRIVATE(zhh,zvv,zvh,RHOhvc,dvh,d_dvh,Dopp,Kdp,Adp,Ah,Av) &
      !$OMP& PRIVATE(diff_back_phase,ceilo_back_true,ceilo_ext,mpl_back_true,mpl_ext,sw_h2,ihtc,work,work1) &
      !$OMP& PRIVATE(wstart3,wstart4,wend2,wend3,wend4) &
      !$OMP& PRIVATE(ibin,sw_dyn,w_r,dist_from_radar) &
      !$OMP& PRIVATE(spectra_bins,zhh_spectra,zvh_spectra,zvv_spectra)
      do ix=1,rmout%nx
#else
      do ix=1,rmout%nx
#endif
      
        if (MAX_TH_NUM > 1) tid=1  ! AT April 2019
        !$  tid = omp_get_thread_num() !AT April 2019
        !$  tid = tid+1 !AT April 2019
  
        do iy=1,rmout%ny
  
          psum_xy = psum_xy + 1
          percentage=nint(real(psum_xy)/real(psum_xy_tot) * 100.e0) 
          IF (TID .EQ. 1) THEN
            if ( (mod(percentage,5) ==0).and.(percentage /=percentage_saved) ) then
              percentage_saved=percentage
              write(*,*) 'Number of Processed profiles in percents: ', percentage
            endif
          ENDIF
  
          !wstart2 = omp_get_wtime()
  
          do iz=1,rmout%nz

    Allocate(spectra_bins(nfft),zhh_spectra(nfft),zvh_spectra(nfft),zvv_spectra(nfft))
    spectra_bins=0.d0; zhh_spectra=0.d0; zvv_spectra=0.d0; zvh_spectra=0.d0
            ! 
            elev=conf%elev
            rr=m999
            !
            ! determine elevation
            ! AUG2019
            !if ( (conf%elev < -90.d0) .or. (conf%elev > 90.d0) )  then !AT 
            if (conf%elev < 0.d0) then
              call determine_elevation_and_range(env%x(ix),env%y(iy),env%z(ix,iy,iz),&
                                                 conf%ixc,conf%iyc,conf%zc,env%dx,env%dy,elev,rr)
              rmout%elev(ix,iy,iz)=elev
              rmout%range(ix,iy,iz)=rr
            endif
            ! 
            if (elev >= 0.d0) then
            ! AUG2019
            !if ( (elev >= -90.d0 ) .and. (elev <= 90.d0) ) then  !AT 
              !wstart3 = omp_get_wtime()
              if (rr<0.d0) then
                call determine_range(env%x(ix),env%y(iy),env%z(ix,iy,iz),conf%ixc,conf%iyc,conf%zc,env%dx,env%dy,rr)
                rmout%range(ix,iy,iz)=rr
              endif
              !
              ! determine azimuth
              call determine_azimuth(env%x(ix),env%y(iy),conf%ixc,conf%iyc,env%dx,env%dy,rmout%elev(ix,iy,iz),azim)
              rmout%azim(ix,iy,iz)=azim    
              !
              ! AUG2019
              ! at this point elev and rmout%elev(ix,iy,iz) are equal 
              ! define elevx (the lut elevation angle from 0 to 90
              elevx = elev ! radar looking up 
              !if (conf%airborne /= 0) elevx = MAX(0.d0, 90.d0 - dabs(elev)) ! radar looking down AT
              ! 
              uu=env%u(ix,iy,iz)
              vv=env%v(ix,iy,iz)
              ww=env%w(ix,iy,iz)    !ww=env%v(ix,iy,iz) modified by oue 2018/09/03
              tke=env%tke(ix,iy,iz)
              !
              qsum=0.d0
              !
              if ((conf%MP_PHYSICS==9 ).or.(conf%MP_PHYSICS==10).or.(conf%MP_PHYSICS==8) .or.&
                  (conf%MP_PHYSICS==30).or.(conf%MP_PHYSICS==40).or.(conf%MP_PHYSICS==75)) then !Modified by oue 2016/09/19, 2017/07/17, 2018/06/17 ICON,RAMS,SAM
                qsum=Sum(hydro%qhydro(ix,iy,iz,1:nsc))  
              elseif (conf%MP_PHYSICS==50) then ! added by DW for P3
                qsum=Sum(hydro50%qhydro(ix,iy,iz,1:nsc)) ! added by DW for P3
              ! 
              elseif (conf%MP_PHYSICS==20) then
                !!$OMP PARALLEL DO PRIVATE(isc)
                do isc=1,nsc
                  ! removing N from pixels in which the mixing ratio <= specifed threshold
                  work=scatt_type(isc)%qq(ix,iy,iz)
                  if (work<conf%thr_mix_ratio(isc)) then
                    scatt_type(isc)%qq(ix,iy,iz)=0.d0
                    scatt_type(isc)%N(ix,iy,iz,:)=0.d0
                  endif 
                  ! 
                  qsum=qsum+scatt_type(isc)%qq(ix,iy,iz) ! total mixing ratio
                enddo  
                !!$OMP END PARALLEL DO
              !endif
      
              elseif (conf%MP_PHYSICS==70) then
                !!$OMP PARALLEL DO PRIVATE(isc)
                do isc=1,nsc
                  qsum=qsum+scatt_type(isc)%qq(ix,iy,iz) ! total mixing ratio
                enddo  
                !!$OMP END PARALLEL DO
              !endif
              else
                write(*,*) 'Problem with qsum==0'
                Exit
              endif
              !
              !wend3 = omp_get_wtime()
              !write(*,*) 'Part 2 - 1 - 1 Time : ',  wend3 - wstart3
  
              IF (qsum > 0.d0) then   ! IF QSUM > 0 
                !wstart3 = omp_get_wtime()
                call nullify_rmrad_var(rmrad_vec(tid))
                IF (conf%radID/=1)  call nullify_mrad_var(mrad_vec(tid))
                ! 
                call determine_sw_contrib_terms(rr,elev,azim,uu,vv,ww,tke,conf%sigma_r,conf%sigma_theta,&
                                                env%Ku(ix,iy,iz),env%Kv(ix,iy,iz),env%Kw(ix,iy,iz),&
                                                sw_t2,sw_s2,sw_v2)
    
                if ((tke > 0.d0).and.(sw_t2 > 0.d0)) rmout%sw_t(ix,iy,iz)=dsqrt(sw_t2)     
                if (sw_s2 > 0.d0) rmout%sw_s(ix,iy,iz)=dsqrt(sw_s2)
                if (sw_v2 > 0.d0) rmout%sw_v(ix,iy,iz)=dsqrt(sw_v2)
     
                !!!! Total spectral broadening estimation except hydrometeor fall velocity, used for Doppler spectra
                sw_dyn = 0.2d0
                if((sw_t2+sw_s2+sw_v2) > 0.d0) sw_dyn = dsqrt(sw_t2+sw_s2+sw_v2)
                !!!! Radial component of wind field, used for Doppler spectra
                ! AT check AUG2019
                w_r = (env%u(ix,iy,iz)*dcos(azim*d2r) + env%v(ix,iy,iz)*dsin(azim*d2r)) * dcos(elev*d2r) + &
                       env%w(ix,iy,iz)*dsin(elev*d2r)
                dist_from_radar = rr
                !
                !!!! Replace range with height from radar for Dopper spectra for vertically-pointing measurements
                !AUG2019
!                AT
!                if (   ( (conf%airborne == 0) .and. (elev>89.d0) .and. (elev<91.d0) ) & !if radar looking up
!                & .or. ( (conf%airborne /= 0) .and. (elev>-1.d0) .and. (elev<1.d0) ) ) & !or it looks down
!                & then !
                if((conf%elev>89.d0) .and. (conf%elev<91.d0)) then 
                  dist_from_radar = dabs( env%z(ix,iy,iz) - conf%zc ) 
                endif
                !wend3 = omp_get_wtime()
                !write(*,*) 'Part 2 - 1 - 2 Time : ',  wend3 - wstart3
                !wstart3 = omp_get_wtime() AT April 2019
                !
                !
                qq=0.d0
                do iht=1,nsc ! nsc == 5 in small example
    
                  if (conf%MP_PHYSICS== 9) qq=hydro%qhydro(ix,iy,iz,iht)
                  if (conf%MP_PHYSICS==10) qq=hydro%qhydro(ix,iy,iz,iht)
                  if (conf%MP_PHYSICS==20) qq=scatt_type(iht)%qq(ix,iy,iz)
                  if (conf%MP_PHYSICS== 8) qq=hydro%qhydro(ix,iy,iz,iht) !Added by oue 2016/09/19
                  if (conf%MP_PHYSICS==30) qq=hydro%qhydro(ix,iy,iz,iht) !Added by oue 2017/07/17 ICON
                  if (conf%MP_PHYSICS==40) qq=hydro%qhydro(ix,iy,iz,iht) !Added by oue 2017/07/21 RAMS
                  if (conf%MP_PHYSICS==50) qq=hydro50%qhydro(ix,iy,iz,iht) ! added by DW 2017/10/30 P3
                  if (conf%MP_PHYSICS==70) qq=scatt_type(iht)%qq(ix,iy,iz) ! added by oue for SAM warm bin
                  if (conf%MP_PHYSICS==75) qq=hydro%qhydro(ix,iy,iz,iht)  ! added by oue for SAM morr
        
                  !================================================================!
                  !== Doppler spectra
                  !spectra_bins = 0.d0 ; zhh_spectra = 0.d0 ; zvh_spectra = 0.d0 ; zvv_spectra = 0.d0 
                  !================================================================!
  
                  if (qq>0.d0) then 
               
                    !wstart4 = omp_get_wtime()
                    !write(*,*) '------------------' 
                    !write(*,*) 'iht,iz',iht,iz,hydro%qnhydro(ix,iy,iz,iht)        
                    !write(*,*) '------------------'
                    !----------------------------------------------------------------------------------
 
                    if ((conf%MP_PHYSICS== 9).or.(conf%MP_PHYSICS==10).or.(conf%MP_PHYSICS== 8) .or. &
                        (conf%MP_PHYSICS==30).or.(conf%MP_PHYSICS==40).or.(conf%MP_PHYSICS==75)) then !Modified by oue 2016/09/19, 2017/07/17 2018/06/17 ICON RAMS SAM
                      !AUG2019
                      ! AT NOTE: two parameteres relative for airborne obs enter
                      ! into processing subrs: conf% airborne and elevx
                      call processing(iht, conf,elevx,env%w(ix,iy,iz),&
                                      env%temp(ix,iy,iz),env%rho_d(ix,iy,iz), env%rho_d(ix,iy,1), &
                                      hydro%qhydro(ix,iy,iz,iht),hydro%qnhydro(ix,iy,iz,iht),&
                                      spectra%VNyquist,spectra%NOISE_1km,spectra%NFFT,spectra%Nave,&
                                      dist_from_radar,w_r,sw_dyn,&
                                      zhh,zvv,zvh,RHOhvc,&
                                      dvh,d_dvh,Dopp,&
                                      Kdp,Adp,Ah,Av,&
                                      diff_back_phase,&
                                      ceilo_back_true,ceilo_ext,&
                                      mpl_back_true,mpl_ext,&
                                      spectra_bins,zhh_spectra,zvh_spectra,zvv_spectra) 
                    endif
                    !
                    if (conf%MP_PHYSICS==50) then ! added by DW 2017/10/30 P3
                      call processing_P3(iht,conf,elevx,env%w(ix,iy,iz),&
                                         env%temp(ix,iy,iz),env%rho_d(ix,iy,iz),env%rho_d(ix,iy,1),& 
                                         hydro50%qhydro(ix,iy,iz,iht),hydro50%qnhydro(ix,iy,iz,iht),&
                                         hydro50%qir(ix,iy,iz),hydro50%qib(ix,iy,iz),&
                                         spectra%VNyquist,spectra%NOISE_1km,spectra%NFFT,spectra%Nave,&
                                         dist_from_radar,w_r,sw_dyn,&
                                         zhh,zvv,zvh,RHOhvc,&
                                         dvh,d_dvh,Dopp,&
                                         Kdp,Adp,Ah,Av,&
                                         diff_back_phase,&
                                         ceilo_back_true,ceilo_ext,&
                                         mpl_back_true,mpl_ext,&
                                         spectra_bins,zhh_spectra,zvh_spectra,zvv_spectra) 
                    endif ! added by DW
                    !
                    if (conf%MP_PHYSICS==20) then 
                      call processing_sbm(iht,conf,elevx,&
                                          env%w(ix,iy,iz),env%temp(ix,iy,iz),scatt_type(iht)%nbins,&
                                          (scatt_type(iht)%N(ix,iy,iz,:)/scatt_type(iht)%mass)*env%rho_d(ix,iy,iz),&  ! N in 1/m^3
                                          scatt_type(iht)%diam,scatt_type(iht)%rho,&  ! diam in m, rho in kg/m^3
                                          scatt_type(iht)%fvel,&  ! m/s
                                          spectra%VNyquist,spectra%NOISE_1km,spectra%NFFT,spectra%Nave,&
                                          dist_from_radar,w_r,sw_dyn,scatt_type(iht)%diam1,&
                                          zhh,zvv,zvh,RHOhvc,&
                                          dvh,d_dvh,Dopp,&
                                          Kdp,Adp,Ah,Av,&
                                          diff_back_phase,&
                                          ceilo_back_true,ceilo_ext,&
                                          mpl_back_true,mpl_ext,&
                                          spectra_bins,zhh_spectra,zvh_spectra,zvv_spectra) 
                    endif
                    ! 
                    if (conf%MP_PHYSICS==70) then !Added in April 2018 for SAM warmbin
                      call hydro_samsbm_2mom(iht,scatt_type(iht)%nbins,env%rho_d(ix,iy,iz),&
                                             scatt_type(iht)%diam1,scatt_type(iht)%mass1,&
                                             scatt_type(iht)%NN(ix,iy,iz,:),scatt_type(iht)%N(ix,iy,iz,:),&
                                             scatt_type(iht)%diam2(ix,iy,iz,:),scatt_type(iht)%rho2(ix,iy,iz,:),&
                                             scatt_type(iht)%fvel2(ix,iy,iz,:))
                      call processing_sbm(iht,conf,elevx,&
                                          env%w(ix,iy,iz),env%temp(ix,iy,iz),scatt_type(iht)%nbins,&
                                          scatt_type(iht)%NN(ix,iy,iz,:),&  ! N in 1/m^3
                                          scatt_type(iht)%diam2(ix,iy,iz,:),scatt_type(iht)%rho2(ix,iy,iz,:),& ! diam in m, rho in kg/m^3
                                          scatt_type(iht)%fvel2(ix,iy,iz,:),&  ! m/s
                                          spectra%VNyquist,spectra%NOISE_1km,spectra%NFFT,spectra%Nave,&
                                          dist_from_radar,w_r,sw_dyn,scatt_type(iht)%diam1,&
                                          zhh,zvv,zvh,RHOhvc,&
                                          dvh,d_dvh,Dopp,&
                                          Kdp,Adp,Ah,Av,&
                                          diff_back_phase,&
                                          ceilo_back_true,ceilo_ext,&
                                          mpl_back_true,mpl_ext,&
                                          spectra_bins,zhh_spectra,zvh_spectra,zvv_spectra) 
  
                    endif
                    ! 
                    !wend4 = omp_get_wtime()
                    !write(*,*) 'Part 2 - 1 - 3 - 1 Time : ',  wend4 - wstart4
                    !wstart4 = omp_get_wtime()
                    !
                    rmrad_vec(tid)%zhh(iht)=zhh
                    rmrad_vec(tid)%Dopp(iht)= Dopp ! at el=90, no cont of azimuth 
                    rmrad_vec(tid)%d_dvh(iht)=d_dvh ! at el=90, no cont of azimuth 
                    rmrad_vec(tid)%dvh(iht)=dvh ! rwv, def: w=0 m/s, el=90
                    rmrad_vec(tid)%Ah(iht)=Ah
                    !
                    IF (conf%radID/=1) THEN  ! POLARIMETRIC 
                      mrad_vec(tid)%zvv(iht)=zvv
                      mrad_vec(tid)%zvh(iht)=zvh
                      mrad_vec(tid)%RHOhvc(iht)=RHOhvc
                      mrad_vec(tid)%Kdp(iht)=Kdp
                      mrad_vec(tid)%Adp(iht)=Adp         
                      mrad_vec(tid)%Av(iht)=Av
                      mrad_vec(tid)%diff_back_phase(iht)=diff_back_phase
                    ENDIF
                    !----------------------------------------------------------------------------------
                    !
                    if (rmrad_vec(tid)%zhh(iht)>Zthr) then ! if reflectivity larger than -100 dBZ  
                      rmout%Zhh(ix,iy,iz,iht) =rmrad_vec(tid)%zhh(iht)
                      rmout%Dopp(ix,iy,iz,iht)=(uu*dcos(azim*d2r) + &
                                                vv*dsin(azim*d2r) ) * dcos(elev*d2r) + &
                                                dsin(elev*d2r)*rmrad_vec(tid)%Dopp(iht)/rmrad_vec(tid)%zhh(iht)
                      sw_h2=(dsin(elev*d2r)*dsin(elev*d2r)) * rmrad_vec(tid)%d_dvh(iht)/rmrad_vec(tid)%zhh(iht)
                      rmout%dDVh(ix,iy,iz,iht)=dsqrt( sw_h2 ) 
                      rmout%SWt(ix,iy,iz,iht)=dsqrt( sw_h2 + sw_t2 + sw_s2 + sw_v2)
                      rmout%DVh(ix,iy,iz,iht)=rmrad_vec(tid)%dvh(iht)/rmrad_vec(tid)%zhh(iht) ! refl. weigh. velocity (def: w=0, elev=90 deg)
                      ! 
                      rmout%Dopp90(ix,iy,iz,iht)=rmrad_vec(tid)%Dopp(iht)/rmrad_vec(tid)%zhh(iht)
                      rmout%dDVh90(ix,iy,iz,iht)=dsqrt(rmrad_vec(tid)%d_dvh(iht)/rmrad_vec(tid)%zhh(iht))
                      rmout%Ah(ix,iy,iz,iht)  =rmrad_vec(tid)%Ah(iht)
                      !
                    endif
                    !
                    ihtc=-1
                    !------------------------------------------------------
                    ! CEILO
                    ! for lidar variables we don't apply the reflectivity
                    ! threshold -100 dBZ
                    if ((conf%ceiloID==1) .and. ((iht==1).or.(iht==7))) then
                      If (iht==1) ihtc=1   ! cloud  
                      If (iht==7) ihtc=2   ! drizzle
                      !
                      lout%ceilo_back_true(ix,iy,iz,ihtc)=ceilo_back_true
                      lout%ceilo_ext(ix,iy,iz,ihtc)=ceilo_ext
                      lout%lidar_ratio(ix,iy,iz,ihtc)=ceilo_ext/ceilo_back_true ! sr
                    endif
                    !
                    ihtc=-1
                    !------------------------------------------------------
                    ! MPL
                    !for  mpl lidar variables we don't apply the reflectivity
                    !threshold -100 dBZ 
                    if (conf%mplID>0) Then
                      if ((iht==1).or.(iht==3).or.(iht==7)) Then
                        !
                        If (iht==1) ihtc=1   ! cloud
                        If (iht==3) ihtc=2   ! ice
                        If (iht==7) ihtc=3   ! drizzle
                        ! 
                        mpl%back_true(ix,iy,iz,ihtc)=mpl_back_true
                        mpl%ext(ix,iy,iz,ihtc)=mpl_ext
                        mpl%lidar_ratio(ix,iy,iz,ihtc)=mpl_ext/mpl_back_true ! sr
                      endif
                    endif
                    !------------------------------------------------------
                    !------------------------------------------------------
                    ! Doppler spectra
                    if (conf%spectraID==1) Then
                      if(iht==1) spectra%vel_bins = spectra_bins
                      if (rmrad_vec(tid)%zhh(iht)>Zthr) then ! if total reflectivity larger than -100 dBZ 
                        spectra%zhh_spectra(ix,iy,iz,:,iht) = zhh_spectra(:)
                        !write(*,*) iht, 10.0*dlog10(Maxval(zhh_spectra))
                        if (conf%radID/=1) Then !polarimetry
                          spectra%zvh_spectra(ix,iy,iz,:,iht) = zvh_spectra(:)
                          spectra%zvv_spectra(ix,iy,iz,:,iht) = zvv_spectra(:)
                        end if
                      end if
                    end if
                    !------------------------------------------------------
                    !------------------------------------------------------
                    !
                    !----------------------------------------------------------------------------------
                    !----------------------------------------------------------------------------------
                    IF (conf%radID/=1) THEN  ! POLARIMETRIC
                      !
                      if (mrad_vec(tid)%zvv(iht)>Zthr)  then
                        mout%Zvv(ix,iy,iz,iht) =mrad_vec(tid)%zvv(iht)
                        mout%Av(ix,iy,iz,iht)  =mrad_vec(tid)%Av(iht)
                        mout%diff_back_phase(ix,iy,iz,iht)  =mrad_vec(tid)%diff_back_phase(iht)
                        !
                        if (rmrad_vec(tid)%zhh(iht)>Zthr) then
                          mout%Zvh(ix,iy,iz,iht) =mrad_vec(tid)%zvh(iht)
                          mout%RHOhv(ix,iy,iz,iht)=mrad_vec(tid)%RHOhvc(iht)/&
                            dsqrt(rmrad_vec(tid)%zhh(iht)*mrad_vec(tid)%zvv(iht))
                          mout%Zdr(ix,iy,iz,iht)  = rmrad_vec(tid)%zhh(iht)/mrad_vec(tid)%zvv(iht)
                          mout%LDRh(ix,iy,iz,iht) = mrad_vec(tid)%zvh(iht)/rmrad_vec(tid)%zhh(iht)
                          !
                          mout%Kdp(ix,iy,iz,iht) =mrad_vec(tid)%Kdp(iht)
                          mout%Adp(ix,iy,iz,iht) =mrad_vec(tid)%Adp(iht)
                        endif
                        !
                      endif
                      !
                    ENDIF ! IF  POLARIMETRIC 
                   !----------------------------------------------------------------------------------
                   !----------------------------------------------------------------------------------

                   !wend4 = omp_get_wtime()
                   !write(*,*) 'Part 2 - 1 - 3 - 2 Time : ',  wend4 - wstart4
                  endif ! IF qq > 0
                enddo  ! iht
                !wend3 = omp_get_wtime()
                !write(*,*) 'Part 2 - 1 - 3 Time : ',  wend3 - wstart3
                !wstart3 = omp_get_wtime()
                !
                !
                !---------------------------------------------------------------------------------
                ! TOTAL HYDROMETEOR 
                work=0.d0; work1=0.d0
                !
                work=Sum(rmrad_vec(tid)%zhh,Mask=rmrad_vec(tid)%zhh>Zthr) 
                if (work>Zthr) then
                  rmout%Zhh_tot(ix,iy,iz)=work
                  rmout%Dopp_tot(ix,iy,iz)=( uu*dcos(azim*d2r) + vv*dsin(azim*d2r) ) * dcos(elev*d2r)  +&
                                             dsin(elev*d2r) * Sum(rmrad_vec(tid)%Dopp)/rmout%Zhh_tot(ix,iy,iz)
                  sw_h2=dsin(elev*d2r)*dsin(elev*d2r) * (Sum(rmrad_vec(tid)%d_dvh)/rmout%Zhh_tot(ix,iy,iz))
                  rmout%dDVh_tot(ix,iy,iz)=dsqrt( sw_h2)
                  rmout%SWt_tot(ix,iy,iz)=dsqrt(sw_h2 + sw_t2 + sw_s2 + sw_v2)
                  rmout%DVh_tot(ix,iy,iz)=Sum(rmrad_vec(tid)%dvh)/rmout%Zhh_tot(ix,iy,iz)
                  !
                  rmout%Dopp90_tot(ix,iy,iz)=Sum(rmrad_vec(tid)%Dopp)/rmout%Zhh_tot(ix,iy,iz)
                  rmout%dDVh90_tot(ix,iy,iz)=dsqrt(Sum(rmrad_vec(tid)%d_dvh)/rmout%Zhh_tot(ix,iy,iz))
                  !
                  rmout%Ah_tot(ix,iy,iz)  =Sum(rmrad_vec(tid)%Ah)
                endif
                !
                !---------------------------------------------------------------------------------
                If (conf%ceiloID>0) Then ! Ceilometer measurements
                  work=Sum(lout%ceilo_back_true(ix,iy,iz,:))
                  lout%ceilo_back_true_tot(ix,iy,iz)=work
                  lout%ceilo_ext_tot(ix,iy,iz)=Sum(lout%ceilo_ext(ix,iy,iz,:))
                  lout%lidar_ratio_tot(ix,iy,iz)=lout%ceilo_ext_tot(ix,iy,iz)/work
                Endif
                !---------------------------------------------------------------------------------
                If (conf%mplID>0) Then ! MPL measurements
                  If (conf%aeroID==1) Then ! If aerosol
                    work=Sum(mpl%back_true(ix,iy,iz,:)) + aero%back_true(ix,iy,iz) + mpl%rayleigh_back(ix,iy,iz)
                    mpl%back_true_tot(ix,iy,iz)=work
                    mpl%ext_tot(ix,iy,iz)=Sum(mpl%ext(ix,iy,iz,:)) + aero%ext(ix,iy,iz)
                    mpl%lidar_ratio_tot(ix,iy,iz)=mpl%ext_tot(ix,iy,iz)/work
                    !-- To store mpl_backscatter for hydro and aerosol added by oue 2017.03.23
                    If (conf%arsclID==1) Then 
                      mpl_back_hydro(ix,iy,iz) = Sum(mpl%back_true(ix,iy,iz,:))
                      mpl_back_aero(ix,iy,iz) = aero%back_true(ix,iy,iz)
                    Endif
                    !--
                  Else ! if aeroID/=1 cloud,ice only + molecular backscattering
                    work=Sum(mpl%back_true(ix,iy,iz,:)) + mpl%rayleigh_back(ix,iy,iz) 
                    mpl%back_true_tot(ix,iy,iz)=work
                    mpl%ext_tot(ix,iy,iz)=Sum(mpl%ext(ix,iy,iz,:))  
                    mpl%lidar_ratio_tot(ix,iy,iz)=mpl%ext_tot(ix,iy,iz)/work
                    !-- To store mpl_backscatter for hydro and aerosol added by oue 2017.03.23
                    If (conf%arsclID==1) Then 
                      mpl_back_hydro(ix,iy,iz) = Sum(mpl%back_true(ix,iy,iz,:))
                      mpl_back_aero(ix,iy,iz) = 0.d0
                    Endif
                    !--
                  EndIF
                EndIf
                !---------------------------------------------------------------------------------
                ! Doppler spectra
                if (conf%spectraID==1) Then
                  !work=Sum(rmrad_vec(tid)%zhh,Mask=rmrad_vec(tid)%zhh>Zthr) 
                  !if (work>Zthr) then ! if total reflectivity larger than -100 dBZ  
                  !$OMP PARALLEL DO PRIVATE(ibin,work)
                  do ibin=1,NFFT !spectra%NFFT
                    work=Sum(spectra%zhh_spectra(ix,iy,iz,ibin,:))
                    if(work>0.d0)then
                      spectra%zhh_spectra_tot(ix,iy,iz,ibin) = 10.d0 * dlog10(work)
                    else
                      spectra%zhh_spectra_tot(ix,iy,iz,ibin) = m999
                    end if
                    !
                    if (conf%radID/=1) Then !polarimetry
                      !
                      work=Sum(spectra%zvh_spectra(ix,iy,iz,ibin,:)) 
                      if(work>0.d0)then
                        spectra%zvh_spectra_tot(ix,iy,iz,ibin) = 10.d0 * dlog10(work)
                      else
                        spectra%zvh_spectra_tot(ix,iy,iz,ibin) = m999
                      end if
                      !
                      work=Sum(spectra%zvv_spectra(ix,iy,iz,ibin,:))
                      if(work>0.d0)then
                        spectra%zvv_spectra_tot(ix,iy,iz,ibin) = 10.d0 * dlog10(work)
                      else
                        spectra%zvv_spectra_tot(ix,iy,iz,ibin) = m999
                      end if
                      !
                    end if
                    !
                    !$OMP PARALLEL DO PRIVATE(iht)
                    do iht=1,nsc
                      !write(*,*) spectra%zhh_spectra(ix,iy,iz,ibin,iht)
                      if(spectra%zhh_spectra(ix,iy,iz,ibin,iht)>0.d0)then
                        spectra%zhh_spectra(ix,iy,iz,ibin,iht) = 10.d0 * dlog10(spectra%zhh_spectra(ix,iy,iz,ibin,iht))
                      else
                        spectra%zhh_spectra(ix,iy,iz,ibin,iht) = m999
                      end if
                      if (conf%radID/=1) Then !polarimetry
                        if(spectra%zvh_spectra(ix,iy,iz,ibin,iht)>0.d0)then
                          spectra%zvh_spectra(ix,iy,iz,ibin,iht) = 10.d0 * dlog10(spectra%zvh_spectra(ix,iy,iz,ibin,iht))
                        else
                          spectra%zvh_spectra(ix,iy,iz,ibin,iht) = m999 
                        end if
                        if(spectra%zvv_spectra(ix,iy,iz,ibin,iht)>0.d0)then
                          spectra%zvv_spectra(ix,iy,iz,ibin,iht) = 10.d0 * dlog10(spectra%zvv_spectra(ix,iy,iz,ibin,iht))
                        else
                          spectra%zvv_spectra(ix,iy,iz,ibin,iht) = m999
                        end if
                      end if
                    end do
                    !$OMP END PARALLEL DO
                  end do
                  !$OMP END PARALLEL DO
                end if
                !---------------------------------------------------------------------------------
                !------------------------------------------------------
                !
                IF (conf%radID/=1) THEN !POLARIMETRIC
                  ! 
                  work1=Sum(mrad_vec(tid)%zvv, Mask=mrad_vec(tid)%zvv>Zthr)
                  !
                  if (work1>Zthr) then
                    mout%Zvv_tot(ix,iy,iz)=work1
                    mout%Av_tot(ix,iy,iz)  =Sum(mrad_vec(tid)%Av)
                    mout%diff_back_phase_tot(ix,iy,iz)=Sum(mrad_vec(tid)%diff_back_phase)
                    !
                    if (work>Zthr) then
                      mout%Zvh_tot(ix,iy,iz)=Sum(mrad_vec(tid)%zvh)
                      mout%Zdr_tot(ix,iy,iz)=rmout%Zhh_tot(ix,iy,iz)/mout%Zvv_tot(ix,iy,iz)
                      mout%LDRh_tot(ix,iy,iz)=mout%Zvh_tot(ix,iy,iz)/rmout%Zhh_tot(ix,iy,iz)
                      mout%Kdp_tot(ix,iy,iz) =Sum(mrad_vec(tid)%Kdp)
                      mout%Adp_tot(ix,iy,iz) =Sum(mrad_vec(tid)%Adp)
                      if ( (mout%Zvv_tot(ix,iy,iz) >0.e0) .and. (rmout%Zhh_tot(ix,iy,iz) >0.d0)) Then
                        mout%RHOhv_tot(ix,iy,iz) =Sum(mrad_vec(tid)%RHOhvc)/&
                                         dsqrt(mout%Zvv_tot(ix,iy,iz)*rmout%Zhh_tot(ix,iy,iz))
                      endif
                    endif
                    !
                  endif
                  !
                  ! units from linear to logaritmic
                  if (mout%Zvv_tot(ix,iy,iz)>0.d0) & 
                                   mout%Zvv_tot(ix,iy,iz)=10.d0 * dlog10( mout%Zvv_tot(ix,iy,iz) )
                  if (mout%Zdr_tot(ix,iy,iz)>0.d0) &
                                   mout%Zdr_tot(ix,iy,iz)=10.d0 * dlog10( mout%Zdr_tot(ix,iy,iz) )
                  if (mout%LDRh_tot(ix,iy,iz)>0.d0) &
                                   mout%LDRh_tot(ix,iy,iz)=10.d0 * dlog10( mout%LDRh_tot(ix,iy,iz) )
                  if (mout%Zvh_tot(ix,iy,iz)>0.d0) &
                                   mout%Zvh_tot(ix,iy,iz)=10.d0 * dlog10(mout%Zvh_tot(ix,iy,iz) )
   
                ENDIF
                ! --------------------------------------------------------------------------------
                !
                ! units from linear to logaritmic
                if (rmout%Zhh_tot(ix,iy,iz)>0.d0) then
                  rmout%Zhh_tot(ix,iy,iz)=10.d0 * dlog10(rmout%Zhh_tot(ix,iy,iz))
                endif
      
                !  change units: linear to logarithmic for species     
                !!$OMP PARALLEL DO PRIVATE(iht)
                do iht=1,nsc
                  if (rmout%Zhh(ix,iy,iz,iht) > 0.d0) &
                      rmout%Zhh(ix,iy,iz,iht)=10.d0 * dlog10( rmout%Zhh(ix,iy,iz,iht) )
                  !
                  IF (conf%radID/=1) THEN !POLARIMETRIC
                    if (mout%Zvv(ix,iy,iz,iht) > 0.d0) &
                        mout%Zvv(ix,iy,iz,iht)=10.d0 * dlog10( mout%Zvv(ix,iy,iz,iht) )
                    if (mout%Zvh(ix,iy,iz,iht) > 0.d0) &
                        mout%Zvh(ix,iy,iz,iht)=10.d0 * dlog10( mout%Zvh(ix,iy,iz,iht) )
                    if (mout%Zdr(ix,iy,iz,iht)> 0.d0) &
                        mout%Zdr(ix,iy,iz,iht)=10.d0 * dlog10( mout%Zdr(ix,iy,iz,iht) )
                    if (mout%LDRh(ix,iy,iz,iht)> 0.d0) &
                        mout%LDRh(ix,iy,iz,iht)=10.d0 * dlog10( mout%LDRh(ix,iy,iz,iht) ) 
                  ENDIF                 
                  !
                enddo ! iht
                !!$OMP END PARALLEL DO
                !-------------------------------------------------------------------------------------
                !wend3 = omp_get_wtime()
                !write(*,*) 'Part 2 - 1 - 4 Time : ',  wend3 - wstart3
              ELSE ! for the cases QSUM==0 
                if (conf%spectraID==1) Then
                  spectra%zhh_spectra_tot(ix,iy,iz,:) = -999.d0
                  spectra%zhh_spectra(ix,iy,iz,:,:) = -999.d0
                  if (conf%radID/=1) Then !polarimetry
                    spectra%zvh_spectra_tot(ix,iy,iz,:) = -999.d0
                    spectra%zvv_spectra_tot(ix,iy,iz,:) = -999.d0
                    spectra%zvh_spectra(ix,iy,iz,:,:) = -999.d0
                    spectra%zvv_spectra(ix,iy,iz,:,:) = -999.d0
                  end if
                end if
              ENDIF ! If (Sum(hydro%qhydro(ix,iy,iz,:) > 0.d0) 
            ELSE ! for elev <0 
              if (conf%spectraID==1) Then
                spectra%zhh_spectra_tot(ix,iy,iz,:) = -999.d0
                spectra%zhh_spectra(ix,iy,iz,:,:) = -999.d0
                if (conf%radID/=1) Then !polarimetry
                  spectra%zvh_spectra_tot(ix,iy,iz,:) = -999.d0
                  spectra%zvv_spectra_tot(ix,iy,iz,:) = -999.d0
                  spectra%zvh_spectra(ix,iy,iz,:,:) = -999.d0
                  spectra%zvv_spectra(ix,iy,iz,:,:) = -999.d0
                end if
              end if
            ENDIF ! if(elev >= 0.d0)
  
 Deallocate(spectra_bins,zhh_spectra,zvh_spectra,zvv_spectra)   
          enddo ! iz
          !call system_clock ( clck_counts_end1, clck_rate1 )
          !write (*, *) 'time=', (clck_counts_end1 - clck_counts_beg1) / real (clck_rate1)
          !write(*,*) ''
    
          !wend2 = omp_get_wtime()
          !write(*,*) 'Part 2 - 1 Time : ',  wend2 - wstart2
        enddo  ! end do iy

#ifdef __PRELOAD_LUT__
      enddo   ! enddo ix
      !$OMP END PARALLEL DO
#else
    enddo  ! enddo ix
#endif
  
  
  
    !Deallocate(spectra_bins,zhh_spectra,zvh_spectra,zvv_spectra)
  
  
    do ix=1,rmout%nx
      do iy=1,rmout%ny
        !---------------------------------------------------------------------------------------------------------------------
        ! CEILOMETER -> account for attenuation due to hydrometeors midified by oue for RAMS 
  
        if (conf%ceiloID==1) then
          work=MaxVal(lout%ceilo_back_true_tot(ix,iy,:))
          IF (work>0.d0) then
            Allocate(ceilo_back_obs(lout%nz))
            ceilo_back_obs=0.d0
            call compute_attenuated_lidar_signal(lout%nz,env%z(ix,iy,:),lout%ceilo_back_true_tot(ix,iy,:),&
                                                 lout%ceilo_ext_tot(ix,iy,:),ceilo_back_obs)
            where (lout%ceilo_back_true_tot(ix,iy,:) > 0.d0) lout%ceilo_back_obs_tot(ix,iy,:)=ceilo_back_obs
            where (lout%ceilo_back_true_tot(ix,iy,:) == 0.d0) lout%ceilo_back_obs_tot(ix,iy,:)=m999
            Deallocate(ceilo_back_obs)
            call estimate_ceilo_first_cloud_base(lout%nz,env%z(ix,iy,:),lout%ceilo_back_obs_tot(ix,iy,:),ceilo_first_cloud_base)
            lout%ceilo_first_cloud_base(ix,iy)=ceilo_first_cloud_base
          ENDIF
        endif
        !---------------------------------------------------------------------------------------------------------------------
        ! MPL -> account for attenuation due to hydrometeors (cloud+ice considered only
        ! + aerosol if conf%aeroID==1)
        !
        ! 
        IF (conf%mplID>0) THEN
          do iht=1,mpl%nht
            if  (MaxVal(mpl%back_true(ix,iy,:,iht))>0.d0) then
              mpl_back_obs=0.d0
              call compute_attenuated_lidar_signal(mpl%nz,env%z(ix,iy,:),mpl%back_true(ix,iy,:,iht),&
                                                   mpl%ext(ix,iy,:,iht),mpl_back_obs)
              where (mpl%back_true(ix,iy,:,iht) > 0.d0) mpl%back_obs(ix,iy,:,iht)=mpl_back_obs
            endif
          enddo
          !
          ! up to this point the mpl "tot" vars are defined only in pixels where the total hydrometeor content>0.
          ! 
          If (conf%aeroID==1) Then ! aerosol case)
            !!$OMP PARALLEL DO PRIVATE(iz,work)
            do iz=1,mpl%nz
              work=mpl%back_true_tot(ix,iy,iz)
              if (work == 0.d0) then ! the case if no hydrometeor is present in the column
                mpl%back_true_tot(ix,iy,iz)=aero%back_true(ix,iy,iz) + mpl%rayleigh_back(ix,iy,iz)
                mpl%ext_tot(ix,iy,iz)=aero%ext(ix,iy,iz)
                mpl%lidar_ratio_tot(ix,iy,iz)=mpl%ext_tot(ix,iy,iz)/mpl%back_true_tot(ix,iy,iz) 
              endif
            enddo
            !!$OMP END PARALLEL DO
            call compute_attenuated_lidar_signal(mpl%nz,env%z(ix,iy,:),mpl%back_true_tot(ix,iy,:),mpl%ext_tot(ix,iy,:),mpl_back_obs)
            mpl%back_obs_tot(ix,iy,:)=mpl_back_obs
          Else  ! back_obs_tot  for no aerosol case
            !!$OMP PARALLEL DO PRIVATE(iz,work)
            do iz=1,mpl%nz
              work=mpl%back_true_tot(ix,iy,iz)
              if (work == 0.d0) then ! the case if no hydrometeor is present in the column
                mpl%back_true_tot(ix,iy,iz)=mpl%rayleigh_back(ix,iy,iz)
                mpl%ext_tot(ix,iy,iz)=0.d0 ! atmospheric moleculs do not have extinction
                mpl%lidar_ratio_tot(ix,iy,iz)=mpl%ext_tot(ix,iy,iz)/mpl%back_true_tot(ix,iy,iz)
              endif
            enddo
            !!$OMP END PARALLEL DO
            call compute_attenuated_lidar_signal(mpl%nz,env%z(ix,iy,:),mpl%back_true_tot(ix,iy,:),mpl%ext_tot(ix,iy,:),mpl_back_obs)
            mpl%back_obs_tot(ix,iy,:)=mpl_back_obs
            !
          Endif
          !
          ! missing values
          do iht=1,mpl%nht ! mpl%nht == 2 in small example
            where (mpl%back_true(ix,iy,:,iht) == 0.d0) mpl%back_true(ix,iy,:,iht)=m999
            where (mpl%ext(ix,iy,:,iht) == 0.d0) mpl%ext(ix,iy,:,iht)=m999
            where (mpl%back_obs(ix,iy,:,iht) == 0.d0) mpl%back_obs(ix,iy,:,iht)=m999
            where (mpl%lidar_ratio(ix,iy,:,iht) == 0.d0) mpl%lidar_ratio(ix,iy,:,iht)=m999
          enddo
          !
          where (mpl%rayleigh_back(ix,iy,:) == 0.d0) mpl%rayleigh_back(ix,iy,:)=m999
          where (mpl%back_true_tot(ix,iy,:) == 0.d0) mpl%back_true_tot(ix,iy,:)=m999
          where (mpl%ext_tot(ix,iy,:) == 0.d0) mpl%ext_tot(ix,iy,:)=m999
          where (mpl%back_obs_tot(ix,iy,:) == 0.d0) mpl%back_obs_tot(ix,iy,:)=m999
          where (mpl%lidar_ratio_tot(ix,iy,:) == 0.d0) mpl%lidar_ratio_tot(ix,iy,:)=m999
        END IF ! if conf%mplID>0
    
      enddo
    enddo
    !---------------------------------------------------------------------------------------------------------------------
    !---------------------------------------------------------------------------------------------------------------------  
    !
    !---------------------------------------------------------------------------------------------------------------------
    
    Deallocate(mpl_back_obs)
  
    !===========================================================================================!
    !=== Post processing to produce observational products using simulated variables ===========!
    !=== added by oue 2017.03.23 
    !--------------------------------------------------------------------------------------------
    !--------------------------------------------------------------------------------------------
    ! ARSCL products 
    If ((conf%arsclID==1) .and. (conf%mplID>0) .and. (conf%ceiloID==1)) Then 
      do ix=1,arscl%nx
        do iy=1,arscl%ny
          work=MaxVal(mpl_back_hydro(ix,iy,:))
          Allocate(radar_zhh_min(arscl%nz))
          Allocate(cloud_detection_flag(arscl%nz),cloud_mask(arscl%nz))
          Allocate(cloud_layer_base(arscl%n_layers),cloud_layer_top(arscl%n_layers))
          if (work > 0.d0) then
            radar_zhh_min(:)=conf%ZMIN + 20.d0 * dlog10(1.d-3 * rmout%range(ix,iy,:))
            !!!! Replace range with height from radar for vertical-pointing measurements
            if((conf%elev>89.d0) .and. (conf%elev<91.d0)) then
            !AUG2019
            ! AT
            !if (   ( (conf%airborne == 0) .and. (elev>89.d0) .and.(elev<91.d0) ) & !if radar looking up
            !& .or. ( (conf%airborne /= 0) .and. (elev>-1.d0) .and. (elev<1.d0) ) ) & !or it looks down
            !& then 
              !dist_from_radar = dabs( env%z(ix,iy,iz) - conf%zc )  AT inserted only in version 3.3.0 but not correct
              !radar_zhh_min(:)=conf%ZMIN + 20.d0 * dlog10(1.d-3 * dist_from_radar) !  so removed for the next version (3.3.1)
              radar_zhh_min(:)=conf%ZMIN + 20.d0 * dlog10(1.d-3 * (env%z(ix,iy,:)-env%z(ix,iy,1)))
            endif
            !
            call estimate_arscl_cloudmask(arscl%nz,env%z(ix,iy,:),arscl%n_layers,&
                                          rmout%Zhh_tot(ix,iy,:),radar_zhh_min,rmout%Ah_tot(ix,iy,:),&
                                          mpl_back_hydro(ix,iy,:), mpl%back_obs_tot(ix,iy,:),&
                                          mpl_back_aero(ix,iy,:), mpl%rayleigh_back(ix,iy,:),&
                                          lout%ceilo_first_cloud_base(ix,iy),&
                                          cloud_detection_flag,cloud_mask,cloud_layer_base,cloud_layer_top)
            arscl%cloud_detection_flag(ix,iy,:)=cloud_detection_flag
            arscl%cloud_mask(ix,iy,:)=cloud_mask
            arscl%cloud_layer_base(ix,iy,:)=cloud_layer_base
            arscl%cloud_layer_top(ix,iy,:)=cloud_layer_top
          else
            arscl%cloud_detection_flag(ix,iy,:)=int2(1)
            arscl%cloud_mask(ix,iy,:)=int2(0)
            arscl%cloud_layer_base(ix,iy,:)=m999
            arscl%cloud_layer_top(ix,iy,:)=m999
          endif
          Deallocate(radar_zhh_min)
          Deallocate(cloud_detection_flag,cloud_mask)
          Deallocate(cloud_layer_base,cloud_layer_top)
        enddo
      enddo
    endif
    if (conf%arsclID==1) then
      Deallocate(mpl_back_hydro,mpl_back_aero) 
    endif
    !--------------------------------------------------------------------------------------------
    !--------------------------------------------------------------------------------------------     
    ! MWR LWP taking account of field of view
    If (conf%mwrID==1) Then
      Allocate(lwc(mwr%nx,mwr%ny,mwr%nz))
      !- Liquid water content (cloud + rain) ! [kg/m3] water content 
      if (conf%MP_PHYSICS== 8)  lwc=(hydro%qhydro(:,:,:,1)+hydro%qhydro(:,:,:,2))*env%rho_d 
      if (conf%MP_PHYSICS== 9)  lwc=(hydro%qhydro(:,:,:,1)+hydro%qhydro(:,:,:,2))*env%rho_d 
      if (conf%MP_PHYSICS== 10) lwc=(hydro%qhydro(:,:,:,1)+hydro%qhydro(:,:,:,2))*env%rho_d 
      if (conf%MP_PHYSICS== 30) lwc=(hydro%qhydro(:,:,:,1)+hydro%qhydro(:,:,:,2))*env%rho_d 
      if (conf%MP_PHYSICS== 40) lwc=(hydro%qhydro(:,:,:,1)+hydro%qhydro(:,:,:,2)&
                                    +hydro%qhydro(:,:,:,7))*env%rho_d 
      if (conf%MP_PHYSICS== 50) lwc=(hydro50%qhydro(:,:,:,1)+hydro50%qhydro(:,:,:,2))*env%rho_d 
      if (conf%MP_PHYSICS== 20) lwc=(scatt_type(1)%qq+scatt_type(2)%qq)*env%rho_d 
      if (conf%MP_PHYSICS== 70) lwc=(scatt_type(1)%qq+scatt_type(2)%qq)*env%rho_d 
      if (conf%MP_PHYSICS== 75) lwc=(hydro%qhydro(:,:,:,1)+hydro%qhydro(:,:,:,2))*env%rho_d 
  
      call estimate_mwr_lwp(mwr%nx,mwr%ny,mwr%nz,env%z,env%dx,env%dy,conf%mwr_view,conf%mwr_alt,&
                            lwc,mwr%true_lwp,mwr%mwr_lwp,mwr%n_samples)
      Deallocate(lwc)
    endif
    !--------------------------------------------------------------------------------------------   
    !== End post processing
    !===========================================================================================!
  
  
    DO isc=1, MAX_TH_NUM
      call deallocate_rmrad_var(rmrad_vec(isc))
      IF (conf%radID/=1)  call deallocate_mrad_var(mrad_vec(isc))
    ENDDO
  
    deallocate(rmrad_vec, mrad_vec)
  
    !---------------------------------------------------
    !------------------------------------------------------------------------------------
    ! write the main output file and the output files for each hydrometeor type
    !
    call get_filename(Trim(conf%WRFInputFile),filename)
    conf%WRFInputFile=Trim(filename)
    !
    call get_path(Trim(conf%OutFile),path)
    call get_filename(Trim(conf%OutFile),filename)
    call get_extention(trim(filename),extention)
    nnf=len(Trim(filename))
    nn=len(Trim(extention))
    filename=Trim(filename(1:nnf-nn-1))
    
    nht=rmout%nht! or =nsc
    
    do iht=0,nht
      if (iht==0) hstring="total"
      if (iht==1) hstring="cloud"
      if (iht==2) hstring="rain"
      if (iht==3) hstring="ice"
      if (iht==4) hstring="snow"
      if (iht==5) hstring="graupel"
      if (iht==6) hstring="hail"
      if (iht==7) hstring="drizzle"
      if (iht==8) hstring="aggregate"
      if ((iht==3).and.(conf%MP_PHYSICS==50)) hstring="smallice"   ! added by DW 2017/10/30 P3 
      if ((iht==4).and.(conf%MP_PHYSICS==50)) hstring="unrimedice" ! added by DW
      if ((iht==5).and.(conf%MP_PHYSICS==50)) hstring="graupel"    ! added by DW
      if ((iht==6).and.(conf%MP_PHYSICS==50)) hstring="parimedice" ! added by DW
      !
      if (iht==0) then
        OutFileName=Trim(conf%OutFile)
      else
        OutFileName=Trim(Adjustl(path))//Trim(Adjustl(filename))//"_"//Trim(Adjustl(hstring))//"."//Trim(Adjustl(extention))
      endif
      !
      write(*,*) 'Writing:',Trim(OutFileName)
      Allocate (qqc(rmout%nx,rmout%ny,rmout%nz))
      if ((iht>=1) .and. (iht<=nht)) then
        if (conf%MP_PHYSICS== 8)  qqc=hydro%qhydro(:,:,:,iht) !added by oue 2016/09/19
        if (conf%MP_PHYSICS== 9)  qqc=hydro%qhydro(:,:,:,iht)
        if (conf%MP_PHYSICS==10)  qqc=hydro%qhydro(:,:,:,iht)
        if (conf%MP_PHYSICS==20)  qqc=scatt_type(iht)%qq
        if (conf%MP_PHYSICS==30)  qqc=hydro%qhydro(:,:,:,iht) !added by oue 2017/07/17 ICON
        if (conf%MP_PHYSICS==40)  qqc=hydro%qhydro(:,:,:,iht) !added by oue 2017/07/21 RAMS
        if (conf%MP_PHYSICS==50)  qqc=hydro50%qhydro(:,:,:,iht) ! added by DW 2017/10/30 P3
        if (conf%MP_PHYSICS==70)  qqc=scatt_type(iht)%qq !added by oue SAM 
        if (conf%MP_PHYSICS==75)  qqc=hydro%qhydro(:,:,:,iht) !added by oue SAM morr
      else if (iht==0) then
        qqc=0.d0
        do isc=1,nht
          if (conf%MP_PHYSICS== 8) qqc=qqc+hydro%qhydro(:,:,:,isc) !added by oue 2016/09/19
          if (conf%MP_PHYSICS== 9) qqc=qqc+hydro%qhydro(:,:,:,isc)
          if (conf%MP_PHYSICS==10) qqc=qqc+hydro%qhydro(:,:,:,isc)
          if (conf%MP_PHYSICS==20) qqc=qqc+scatt_type(isc)%qq
          if (conf%MP_PHYSICS==30) qqc=qqc+hydro%qhydro(:,:,:,isc) !added by oue 2017/07/17 ICON
          if (conf%MP_PHYSICS==40) qqc=qqc+hydro%qhydro(:,:,:,isc) !added by oue 2017/07/21 RAMS
          if (conf%MP_PHYSICS==50) qqc=qqc+hydro50%qhydro(:,:,:,isc) ! added by DW 2017/10/30 P3
          if (conf%MP_PHYSICS==70) qqc=qqc+scatt_type(isc)%qq !added by oue SAM 
          if (conf%MP_PHYSICS==75) qqc=qqc+hydro%qhydro(:,:,:,isc) !added by oue SAM morr
        enddo
      endif
      qqc=qqc*env%rho_d ! [kg/m3] water content
  
      IF (conf%radID==1) THEN
        write(*,*) ''
        write(*,*) 'hydrometeor ',Trim(hstring)
        if (iht==0) then
          write(*,*) 'Min,Max Zh',MinVal(rmout%Zhh_tot,Mask=rmout%Zhh_tot>-888.d0),MaxVal(rmout%Zhh_tot)
          write(*,*) ''
        else
          write(*,*) 'Min,Max Zh',MinVal(rmout%Zhh(:,:,:,iht),Mask=rmout%Zhh(:,:,:,iht)>-888.d0),MaxVal(rmout%Zhh(:,:,:,iht))
          write(*,*) ''
        endif 
        
        call WriteOutNetcdf_rmp(iht,Trim(OutFileName),rmout,qqc,env,conf,lout,mpl,aero,arscl,mwr,spectra,status)!arscl&mwr are added by oue 2017.03.23, spectra added by oue in May 2018
      ELSE 
        write(*,*) ''
        write(*,*) 'hydrometeor ',Trim(hstring)
        if (iht==0) then
          write(*,*) 'Min,Max Zh',MinVal(rmout%Zhh_tot,Mask=rmout%Zhh_tot>-888.d0),MaxVal(rmout%Zhh_tot)
          !write(*,*) 'Min,Max Zdr',MinVal(mout%Zdr_tot,Mask=mout%Zdr_tot>-888.d0),MaxVal(mout%Zdr_tot)
          write(*,*) ''
        else
          write(*,*) 'Min,Max Zh',MinVal(rmout%Zhh(:,:,:,iht),Mask=rmout%Zhh(:,:,:,iht)>-888.d0),MaxVal(rmout%Zhh(:,:,:,iht))
          !write(*,*) 'Min,Max Zdr',MinVal(mout%Zdr(:,:,:,iht),Mask=mout%Zdr(:,:,:,iht)>-888.d0),MaxVal(mout%Zdr(:,:,:,iht))
          write(*,*) ''
        endif
    
        call WriteOutNetcdf_mp(iht,Trim(OutFileName),rmout,mout,qqc,env,conf,lout,mpl,aero,arscl,mwr,spectra,status) !arscl&mwr are added by oue 2017.03.23, spectra added by oue in May 2018
      ENDIF
      Deallocate(qqc)
    enddo ! iht
    !
    call deallocate_rmout_var(rmout)
    IF (conf%radID/=1) call deallocate_mout_var(mout)
    call deallocate_env_var(env)
    !
    !deallocate lout vars
    IF (conf%ceiloID==1) call deallocate_lout_var(lout)
    !
    !deallocate mpl (and aero) vars
    IF (conf%mplID>0) THEN
       call deallocate_mpl_var(mpl)
       IF (conf%aeroID==1) call deallocate_aero_var(aero)
    ENDIF
    !
    !== Vars for Doppler spectra ==
    !deallocate spectra vars
    IF (conf%spectraID==1) THEN
      call deallocate_spectra_var(spectra,conf%radID)
    ENDIF
    !
    !== Vars for post processing====
    !deallocate arscl vars added by oue 2017.03.23
    IF (conf%arsclID==1) THEN
      call deallocate_arscl_var(arscl)
    ENDIF
    !deallocate mwr vars added by oue 2017.03.25
    IF (conf%mwrID==1) THEN
      call deallocate_mwr_var(mwr)
    ENDIF
    !==============================
    !
    !
    if (conf%MP_PHYSICS==08) call deallocate_hydro_var(hydro) !added by oue 2016/09/19
    if (conf%MP_PHYSICS==09) call deallocate_hydro_var(hydro)
    if (conf%MP_PHYSICS==10) call deallocate_hydro_var(hydro)
    if (conf%MP_PHYSICS==30) call deallocate_hydro_var(hydro) !added by oue 2017/07/17 ICON
    if (conf%MP_PHYSICS==40) call deallocate_hydro_var(hydro) !added by oue 2017/07/21 RAMS
    if (conf%MP_PHYSICS==50) call deallocate_hydro50_var(hydro50) ! added by DW 2017/10/30 P3
    if (conf%MP_PHYSICS==75) call deallocate_hydro_var(hydro)  !added by oue SAM morr
    if (conf%MP_PHYSICS==20) then
      do iht=1,nht
        deallocate(scatt_type(iht)%diam)
        deallocate(scatt_type(iht)%rho)
        deallocate(scatt_type(iht)%mass)
        deallocate(scatt_type(iht)%fvel)
        deallocate(scatt_type(iht)%qq)
        deallocate(scatt_type(iht)%N)
      enddo
    endif
    if (conf%MP_PHYSICS==70) then
      do iht=1,nht
        deallocate(scatt_type(iht)%diam2)
        deallocate(scatt_type(iht)%rho2)
        deallocate(scatt_type(iht)%mass)
        deallocate(scatt_type(iht)%fvel2)
        deallocate(scatt_type(iht)%qq)
        deallocate(scatt_type(iht)%N)
        deallocate(scatt_type(iht)%NN)
        deallocate(scatt_type(iht)%diam1)
        deallocate(scatt_type(iht)%mass1)
      enddo
    endif
    !
    call system_clock ( clck_counts_end, clck_rate )
  
    write (*, *) ''
    write (*, *) 'tot. run time in sec', (clck_counts_end - clck_counts_beg) / real (clck_rate)
    write (*, *) ''
    write (*, *) 'Finished!'
    write (*, *) ''
  
    Contains
  
    subroutine get_command_line_arguments(status)
      Use crsim_mod, ONLY: conf_var
      Implicit None
  
      Integer, Intent(Out)             :: status
      !
      Integer                          :: nargs 
      Integer                          :: status2
      Character(len=768)               :: arg_str
      Character(len=256)               :: error_str,error_str2
      !
      error_str='Usage: crsim config_file inp1_file,inp2_file out_file'
  
      status=0
      status2=0
      !
      conf%InpProfile_flag = 0 !added by oue for SAM profile data (third file)
      !
      error_str2=''
      !
      nargs=COMMAND_ARGUMENT_COUNT()
      !
      If (nargs /=3) Then
        status=1
        error_str2='ERROR: Wrong number of arguments'
      Endif
      !
      ! read the name of configuration file
      If (status == 0) Then
        CALL GET_COMMAND_ARGUMENT(1,arg_str)
        ConfigInpFile = Trim(Adjustl(arg_str))
      Endif
      !
      ! read the input file(s)
      If (status == 0) Then
        CALL GET_COMMAND_ARGUMENT(2,arg_str)
        Call get_nargs_sim_input(arg_str,',',nn)
        If ((nn > 3).or. (nn < 1)) Then
          status=1
          error_str2='Error in inputfile list'
        Else
          Call get_arg_sim_input(arg_str,',',1,conf%WRFInputFile,status2)
          conf%WRFmpInputFile=Trim(Adjustl(conf%WRFInputFile))
          if (nn==2) Then
            Call get_arg_sim_input(arg_str,',',2,conf%WRFmpInputFile,status2)
            write(*,*) 'conf%WRFmpInputFile=',trim(conf%WRFmpInputFile)
            conf%InpProfile_flag = 1
          endif
            if (nn==3) Then !added by our for SAM to read profile data
              Call get_arg_sim_input(arg_str,',',2,conf%WRFmpInputFile,status2)
              Call get_arg_sim_input(arg_str,',',3,conf%InputProfile,status2)
              conf%InpProfile_flag = 2
            endif
          Endif
        EndIf
        ! 
        ! read the output file
        If (status == 0) Then
          CALL GET_COMMAND_ARGUMENT(3,arg_str)
          conf%OutFile=Trim(Adjustl(arg_str))
        Endif
        !
        If (status /= 0) Then
          write(*,*) Trim(error_str)
          write(*,*) Trim(error_str2)
        Endif
        ! 
        write(*,*) '*****************'
        write(*,*) 'ConfigInpFile',trim(ConfigInpFile)
        write(*,*) 'conf%WRFInputFile=',trim(conf%WRFInputFile) 
        write(*,*) 'conf%WRFmpInputFile=',trim(conf%WRFmpInputFile)
        write(*,*) 'conf%OutFile=',trim(conf%OutFile) 
        write(*,*) '*****************'
   
      return
    end subroutine get_command_line_arguments
  
  end program crsim
  
