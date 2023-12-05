  !! ----------------------------------------------------------------------------
  !! ----------------------------------------------------------------------------
  !!  *PROGRAM* crsim_subrs 
  !!  @version crsim 3.x
  !!
  !!  *SRC_FILE*
  !!  crsim/src/crsim_subrs.f90
  !!
  !!
  !!
  !!  *LAST CHANGES*
  !!
  !! Sep 15 2015  - A.T.      included horizontal wind in subroutine get_env_vars
  !! Sep 17 2015  - A.T.      written a new subroutine "determine_azimuth".The
  !!                          azimuth is needed in expression of forward radial 
  !!                          Doppler velocity 
  !! Jan 04 2016  - A.T.      added few lines at the end of subroutine get_hydro20_vars 
  !!                          that remove negative values of the input bin concentrations
  !!                          for MP_PHYSICS==20 
  !! Jan  06 2016 - A.T.      change of sign for all velocities. In previous versions
  !!                          velosities  are defines as positive downard i.e. towards
  !!                          the radar. Now we follow radar convenction where velocity
  !!                          is positive upward,i.e. away from the radar. Note that the
  !!                          sign convenction is also applied for reflectivity weighted velocity. 
  !! Jan 12 2016  - A.T.      modified subroutine get_hydro_vars. Input mixing ratios and 
  !!                          total number distributions are set to zero if mixing ratio lower
  !!                          than the input threshold value and if the input N <0. 
  !! Jan 12 2016  - A.T.      modified subroutine get_env_vars in order to include xlat, xlong
  !!                          and tke variables
  !! Jan 12 2016  - A.T.      added new subroutine "determine_sw_contrib_terms" for computation
  !!                          of the spectrum width due to turbulence
  !! Jan 14 2016  - A.T.      included computations for wind shear in get_env_vars.
  !! Jan 15 2016  - A.T.      included contrubutions of wind shear and cross-wind to the spectrum 
  !!                          width in subr. determine_sw_contrib_terms 
  !! Feb 17 2016  - A.T       new sbroutine added: GetCloudLidarMeasurements for cloud lidar 
  !!                          (ceilometer) simulated measurements
  !! Feb 18 2016  - A.T.      Corrected bug in sign when computing Doppler velocity at elevations 
  !!                          different than 90 deg. All velocties are kept positive away from
  !!                          the radar (radar notation). Sign for reflectivity weighted velocity
  !!                          is as for the fall velocity, positive downward or toward the radar 
  !!                          (i.e. related change from Jan 06 is nullified)
  !! March 26 2016 -A.T.      Subroutine get_hydro_vars is modified to work for MP_PHYSICS=9. For this 
  !!                          reason get_hydro_vars is splitted out to 2 subroutines: get_hydro10_vars and
  !!                          get_hydro09_vars. In the case of Morrisson scheme,a possibility of the input 
  !!                          cloud liquid total concentration /=0 is taken into account. 
  !! March 26 2016 -A.T.      In all soubrutines GetPolarimetricInfofromLUT*, when selecting the LUTs densities
  !!                          for hydrometeor isc the line "if (isc==5)" is replaced by "if (isc>=5)" in order
  !!                          to account not only graupel but hail also       
  !! April 18 2016 - A.T      corrected bug in expression for Kazim in subroutine determine_sw_contrib_terms (removed 
  !!                          term dcos(elev*d2r) from the expression)
  !! May  03-05 2016 -A.T     Modified subroutines iprocessing and processing_sbm in order to introduce MPL
  !!                          simulated measurements. Written new soubroutines: GetCloud_MPL_Measurements,
  !!                          GetCloud_MPL_Measurements_sbm  and read_mpl_luts.
  !! May   09 2016   -A.T     Added aerosol MPL measurement. The average aerosol model implemented, as in 
  !!                          Spinhirne D. James, 1993:  Micro Pulse Lidar.IEEE Trans. Geo. Rcem. Sens.,31, 48-55. 
  !!                          New sabroutines added: spinhirne_aero_model and normalize_spinhirne_aero_model
  !! May   09 2016   -A.T     Corrected units for lidar ratio and for extinction coefficient. Lidar ratio 
  !!                          is now in sr and ext. coeff in m^-1 (ano not m^-1 sr^-1).  
  !! May   10 2016   -A.T     Corrected bug in Mie Backscattering particle cross-sections read from LUTs; they are now multiplied by 4 pi.
  !!                          Actually, the backscattering efficiencies computed by Bohren-Huffman Mie scattering 
  !!                          BHMIE code have to be 4 pi larger, as in ftp://ftp.astro.princeton.edu/draine/scat/bhmie/bhmie.f
  !! May  12  2016   -A.T     Included a new subroutine compute_molecular_backscatter computing the molecular backscatter [m-1 sr-1]
  !!
  !! Sep  19  2016   -M.O     Incorporated subroutines for Thompson microphysics. 
  !! May  22  2017   -M.O     Fixed N0 for rain&ice in Thompson microphysics (considering rho_d). 
  !! JUL  17  2017   -M.O     Incorporated subroutines for ICON 2-moment microphysics. 
  !! JUL  21  2017   -M.O     Incorporated subroutines for RAMS 2-moment microphysics.(8 categories) 
  !! JUN  30  2017   -K. YU   Added "__PRELOAD_LUT__" mode.
  !!                          subroutine GetCloudLidarMeasurements_preload(), make_lutfilename(), preload_luts(), preload_mpl_luts(), 
  !!                          and preload_lut_ceilo() are implemented for "__PRELOAD_LUT__" mode.
  !!
  !! Oct  30  2017   -D.W.     Incorporated subroutines for P3 microphysics, adding processing_P3, hydro_info_P3, cloud_P3, rain_P3, 
  !!                          ice_P3, fvel_P3, GetPolarimetricInfofromLUT_P3, get_hydro50_vars, access_lookup_table, 
  !!                          find_lookupTable_indices_1 subroutines. (added by DW)
  !! Apr    2018     -M.O     Incorporated SAM warm bin microphysics (MP_PHTSICS=70)
  !! Jun 17 2018     -M.O     Incorporated SAM Morrison 2-moment microphysics (MP_PHYSICS=75)
  !! MAY  27  2018   -M.O     Incorporated spectra generator for bulk moment microphysics
  !!                          Added suproutines of processing_ds, compute_doppler_spectra, spectra_generator,
  !!                          spectra_unfolding,Turbulence_Convolution,interp1,hsmethod.
  !!                          Modified processing*, GetPolarimetricInfofromLUT*, 
  !!                          compute_polarim_vars, compute_rad_forward_vars
  !!                          Added an option to use hail instead of graupel for MP_PHYSICS=10 (WRF Morrison 2 moment)!!!
  !! Nov  03 2018    -M.O     Incorporated a retrieval of particle size distribution parameters which are used in the Thompson scheme's
  !!                          reflectivity calculation in WRF instead of those described in Thompson et al. (2008 & 2004)
  !! Dec 07 2018     -A.T     Removed line 'if (conf%aero_tau>0.d0) Deallocate(ext_new)' which deallocates already
  !!                          deallocated veriable ext_new in subroutine spinhirne_aero_model
  !! Jul 09 2019     -M.O.    Modified cloud_P3 (Nconst) and rainP3 (mu is fixed) to adapt it for WRFv4.1
  !!                          A bug in hydro_icon (rain part) was fixed.
  !! Aug 19 2019     -A.T.    Modified computation of the specific attenuation at horizontal and vertical polarizations to include
  !!                          the dependencies of Ah and Av on the angular moments.
  !! Aug 20 2019     -A.T     Resolved compilation warning for possible change of value in conversion elsewhere
  !! Aug 27 2019     -A.T     Conversion from real to dounle precission
  !! Sep    2019     -M.O.    Temporary incorporated calculations for negative elevation angles.
  !! Dec    2019     -A.T     Introduced computation of cross-correlation coefficient RHOhv. Affected subroutine are 
  !!                          compute_polarim_vars(),GetPolarimetricInfoFromLUT*(), and processing*().
  !! Apr    2020     -M.O     Incorporated CM1 Morrison 2-moment microphysics (MP_PHYSICS=80)
  !! Dec    2021     -M.O     Morrified cloud_morrison to incorporate different values for Nconst for SAM (MP-75)
  !!
  !!  *DESCRIPTION* 
  !!
  !!  This program  contains the subroutines needed for computation of forward
  !!  radar variables
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
  subroutine processing(isc,conf,elev,ww,temp,rho_d,rho_ds,qhydro,qnhydro,&
                        spectra_VNyquist,spectra_NOISE_1km,NFFT,spectra_Nave,&
                        range_m,w_r,sw_dyn,&
                        Zhh,Zvv,Zvh,RHOhvc,DVh,dDVh,Dopp,Kdp,Adp,Ah,Av,diff_back_phase,&
                        ceilo_back_true,ceilo_ext,mpl_back_true,mpl_ext,&
                        spectra_bins,zhh_spectra,zvh_spectra,zvv_spectra)
  Use crsim_mod
  Use phys_param_mod, ONLY: m999
  Implicit None
  !
  Integer,Intent(In)                 :: isc
  Type(conf_var),Intent(in)          :: conf
  Real*8, Intent(In)                 :: elev 
  Real*8, Intent(In)                 :: ww ![m/s] vertical air velocity
  Real*8, Intent(In)                 :: temp ! [C] temperature 
  Real*8, Intent(In)                 :: rho_d   ! [kg/m^3]  air density
  Real*8, Intent(In)                 :: rho_ds  ! [kg/m^3]  air density at the first atm. level
  Real*8, Intent(In)                 :: qhydro  ! [kg/kg]  mixing ratio 
  Real*8, Intent(In)                 :: qnhydro ! [1/kg]  total concentration 
  !-- input for Doppler spectra
  Real*8, Intent(In)                 :: spectra_VNyquist,spectra_NOISE_1km
  Integer,Intent(In)                 :: NFFT ! size of spectra_bins. If spectraID/=1, this is 1
  Integer,Intent(In)                 :: spectra_Nave
  Real*8, Intent(In)                 :: range_m ! distance from radar [m], used for Doppler spectra simulation
  Real*8, Intent(In)                 :: w_r     ! radial component of wind field [m/s], used for Doppler spectra 
  Real*8, Intent(In)                 :: sw_dyn  ! spectra broadening due to dynamics [m/s], used for Doppler spectra 
  !--
  Real*8, Intent(Out)                :: Zhh,Zvv,Zvh,RHOhvc ! mm^6/m^3
  Real*8, Intent(Out)                :: DVh  ! mm^6/m^3 m/s   assuming ww=0 m/s
  Real*8, Intent(Out)                :: dDVh  ! mm^6/m^3 (m/s)^2  assuming ww=0 m/s
  Real*8, Intent(Out)                :: Dopp  ! mm^6/m^3 m/sa  ww from WRF
  Real*8, Intent(Out)                :: Kdp   ! deg/km
  Real*8, Intent(Out)                :: Adp   ! dB/km
  Real*8, Intent(Out)                :: Ah    ! dB/km
  Real*8, Intent(Out)                :: Av    ! dB/km
  Real*8, Intent(Out)                :: diff_back_phase  ! deg
  Real*8, Intent(Out)                :: ceilo_back_true ! true (unatenuated) ceilo lidar backscatter [m sr]^-1
  Real*8, Intent(Out)                :: ceilo_ext       ! ceilo lidar extinction coefficient [m]^-1
  Real*8, Intent(Out)                :: mpl_back_true ! mpl true (unatenuated) lidar backscatter [m sr]^-1
  Real*8, Intent(Out)                :: mpl_ext       ! mpl lidar extinction coefficient [m]^-1
  !-- Output for Doppler spectra
  Real*8, Intent(Out)                :: spectra_bins(nfft)  ! spectra velocity bin [m/s]
  Real*8, Intent(Out)                :: zhh_spectra(nfft),zvh_spectra(nfft),zvv_spectra(nfft) !output 1D Doppler spectra [mm^6 m^-3 / (m/s)]
  !--
  !
  integer                           :: nd  
  real*8,dimension(:),Allocatable   :: NN    ! 1/m^3
  real*8,dimension(:),Allocatable   :: fvel ! m/s
  real*8,dimension(:),Allocatable   :: diam,ddiam  ! m
  real*8,dimension(:),Allocatable   :: rho   ! kg/m^3
  real*8,dimension(:),Allocatable   :: zhh_d,zvh_d,zvv_d ! mm^6, used for spectrum generation
  !
  real*8                            :: dmin,dmax ! um 
  real*8                            :: Nconst !added Nconst, Dec. 11, 2021 
    !
    ! for each grid point
    !------------------------------------------------------------------------------------
    !
    Call hydro_info(isc,dmin,dmax,nd) ! number of diams, min,max diam for which the scatt info are stored
    !
    Allocate(NN(nd),fvel(nd),diam(nd),rho(nd))
    NN=0.d0 ; fvel=0.d0 ; diam=0.d0; rho=0.d0 ; Nconst=-999.d0 !added Nconst, Dec. 11, 2021
    !!
    IF ((conf%MP_PHYSICS==10) .or. (conf%MP_PHYSICS==75) .or. (conf%MP_PHYSICS==80) ) THEN !add CM1 morrison by oue Apr 2020
      if (isc==1) then
         !added Nconst, Dec. 11, 2021
         Nconst = 250.d0 ! 1/cm^3
         if ((conf%MP_PHYSICS==75)) then
            Nconst = 100.d0 ! 1/cm^3 
         endif
        !call cloud_morrison(nd,dmin,dmax,temp,rho_d,qhydro,qnhydro,diam,rho,NN,fvel)   
        call cloud_morrison(nd,dmin,dmax,temp,rho_d,qhydro,qnhydro,Nconst,diam,rho,NN,fvel) !added Nconst, Dec. 11, 2021   
      else
        call hydro_morrison(conf%MP10HailOption,isc,nd,dmin,dmax,temp,rho_d,qhydro,qnhydro,diam,rho,NN,fvel) 
      endif
    ENDIF
    !!
    !! Added by oue 2016/09/19 ---
    IF (conf%MP_PHYSICS==8) THEN
      call hydro_thompson(isc,nd,dmin,dmax,temp,rho_d,qhydro,qnhydro,diam,rho,NN,fvel) 
    ENDIF
    !!--- added by oue
    !!
    IF (conf%MP_PHYSICS==9) THEN
      call hydro_milbrandt_yau(isc,conf%snow_spherical,nd,dmin,dmax,temp,rho_d,rho_ds,qhydro,qnhydro,diam,rho,NN,fvel)
    ENDIF
    !!
    !! Added by oue 2017/07/17 --- for ICON
    IF (conf%MP_PHYSICS==30) THEN
      call hydro_icon(isc,nd,dmin,dmax,temp,rho_d,qhydro,qnhydro,diam,rho,NN,fvel) 
    ENDIF
    !!--- added by oue
    !! Added by oue 2017/07/21 --- for RAMS
    IF (conf%MP_PHYSICS==40) THEN
      call hydro_rams(isc,nd,dmin,dmax,temp,rho_d,qhydro,qnhydro,diam,rho,NN,fvel) 
    ENDIF
    !!--- added by oue
    Allocate(zhh_d(nd),zvh_d(nd),zvv_d(nd)) ! used for Doppler spectra simulation
    zhh_d=0.d0 ; zvh_d=0.d0 ; zvv_d=0.d0  
    !------------------------------------------------------------------------------------
    IF ((conf%MP_PHYSICS==9).and.(isc==4).and.(conf%snow_spherical/=1)) THEN
      ! when hydrometeor density changes with size
      call GetPolarimetricInfofromLUT(isc,conf,elev,ww,temp,nd,diam,NN,rho,fvel,Zhh,Zvv,Zvh,RHOhvc,&
                                      DVh,dDVh,Dopp,Kdp,Adp,Ah,Av,diff_back_phase,zhh_d,zvh_d,zvv_d)
    ELSE IF ((conf%MP_PHYSICS==8).and.(isc==4)) THEN
      ! when hydrometeor density changes with size
      call GetPolarimetricInfofromLUT(isc,conf,elev,ww,temp,nd,diam,NN,rho,fvel,Zhh,Zvv,Zvh,RHOhvc,&
                                      DVh,dDVh,Dopp,Kdp,Adp,Ah,Av,diff_back_phase,zhh_d,zvh_d,zvv_d)
    ELSE ! density of hydrometeor doesn't change with size
      call GetPolarimetricInfofromLUT_cdws(isc,conf,elev,ww,temp,nd,diam,NN,rho(1),fvel,Zhh,Zvv,Zvh,RHOhvc,&
                                      DVh,dDVh,Dopp,Kdp,Adp,Ah,Av,diff_back_phase,zhh_d,zvh_d,zvv_d)
    ENDIF
    !
    !=======================================================!
    !-- Simulate Doppler spectrum---------------------------!
    !call spectrum subroutine
    spectra_bins = 0.d0 ; zhh_spectra = 0.d0 ; zvh_spectra = 0.d0 ; zvv_spectra = 0.d0 
    if(conf%spectraID==1)then
      Allocate(ddiam(nd+1))
      ddiam = m999 !!!diameter bin sizes [m]
      call processing_ds(isc,conf,spectra_VNyquist,spectra_NOISE_1km,NFFT,spectra_Nave,&
                         &range_m,elev,w_r,sw_dyn,nd,NN,diam,ddiam,fvel,zhh_d,zvh_d,zvv_d,&
                         &spectra_bins,zhh_spectra,zvh_spectra,zvv_spectra)
      Deallocate(ddiam)
    end if
    !---End Simulate Doppler spectrum---------------------------!
    !=======================================================!
    !   
    Deallocate(zhh_d,zvh_d,zvv_d)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !------------------------------------------------------------------------------------
    !
    ! --- introduce cloud ceilo lidar measurements ----
    ceilo_back_true=0.d0 ; ceilo_ext=0.d0
    if ( ((isc==1).or.(isc==7)) .and. (conf%ceiloID==1)) then
#ifdef __PRELOAD_LUT__
      call GetCloudLidarMeasurements_preload(nd,diam,NN,ceilo_back_true,ceilo_ext)
#else
      call GetCloudLidarMeasurements(conf,nd,diam,NN,ceilo_back_true,ceilo_ext)
#endif
    endif
    !----------------------------------------------
    ! --- introduce MPL measurements for liquid and ice clouds ----
    mpl_back_true=0.d0 ; mpl_ext=0.d0
    IF  (conf%mplID>0) THEN
      if ( (isc==1) .or. (isc==3).or.(isc==7)) then
        call GetCloud_MPL_Measurements(isc,conf,nd,diam,NN,mpl_back_true,mpl_ext)
      endif
    ENDIF
    !----------------------------------------------
    !
    Deallocate(NN,fvel,diam,rho)
    !
  return
  end subroutine processing
  !
  !-------------------------------------------------------------------------------
  ! subroutine processing_P3 is added by DW 2017/10/30 
  !
  subroutine processing_P3(isc,conf,elev,ww,temp,rho_d,rho_ds,qhydro,qnhydro,&
                           qir,qib,&
                           spectra_VNyquist,spectra_NOISE_1km,NFFT,spectra_Nave,&
                           range_m,w_r,sw_dyn,&
                           Zhh,Zvv,Zvh,RHOhvc,&
                           DVh,dDVh,Dopp,Kdp,Adp,Ah,Av,diff_back_phase,&
                           ceilo_back_true,ceilo_ext,mpl_back_true,mpl_ext,&
                           spectra_bins,zhh_spectra,zvh_spectra,zvv_spectra)
  ! 
  Use crsim_mod
  Use phys_param_mod, ONLY: Rd, eps
  Implicit None
  !
  Integer,Intent(In)                 :: isc
  Type(conf_var),Intent(in)          :: conf
  Real*8, Intent(In)                 :: elev
  Real*8, Intent(In)                 :: ww      ! [m/s] vertical air velocity
  Real*8, Intent(In)                 :: temp    ! [C] temperature
  Real*8, Intent(In)                 :: rho_d   ! [kg/m^3] air density
  Real*8, Intent(In)                 :: rho_ds  ! [kg/m^3] air density at the first atm. level
  character(len=365)                 :: lookup_file_1
  !
  Real*8, Intent(InOut)              :: qhydro  ! [kg/kg]  mixing ratio 
  Real*8, Intent(InOut)              :: qnhydro ! [1/kg]   total concentration 
  Real*8, Intent(In)                 :: qir     ! [kg/kg]  Rime mixing ratio 
  Real*8, Intent(In)                 :: qib     ! [m3/kg]  Rime ice volume mixing ratio 
  !-- input for Doppler spectra
  Real*8, Intent(In)                 :: spectra_VNyquist,spectra_NOISE_1km
  Integer,Intent(In)                 :: spectra_Nave
  Integer,Intent(In)                 :: NFFT ! size of spectra_bins. If spectraID/=1, this is 1
  Real*8, Intent(In)                 :: range_m ! distance from radar [m], used for Doppler spectra simulation
  Real*8, Intent(In)                 :: w_r     ! radial component of wind field [m/s], used for Doppler spectra 
  Real*8, Intent(In)                 :: sw_dyn  ! spectra broadening due to dynamics [m/s], used for Doppler spectra 
  !--
  Real*8, Intent(Out)                :: Zhh,Zvv,Zvh,RHOhvc ! mm^6/m^3
  Real*8, Intent(Out)                :: DVh   ! mm^6/m^3 m/s   assuming ww=0 m/s
  Real*8, Intent(Out)                :: dDVh  ! mm^6/m^3 (m/s)^2  assuming ww=0 m/s
  Real*8, Intent(Out)                :: Dopp  ! mm^6/m^3 m/sa  ww from WRF
  Real*8, Intent(Out)                :: Kdp   ! deg/km
  Real*8, Intent(Out)                :: Adp   ! dB/km
  Real*8, Intent(Out)                :: Ah    ! dB/km
  Real*8, Intent(Out)                :: Av    ! dB/km
  Real*8, Intent(Out)                :: diff_back_phase  ! deg
  Real*8, Intent(Out)                :: ceilo_back_true  ! true (unatenuated) ceilo lidar backscatter [m sr]^-1
  Real*8, Intent(Out)                :: ceilo_ext        ! ceilo lidar extinction coefficient [m]^-1
  Real*8, Intent(Out)                :: mpl_back_true    ! mpl true (unatenuated) lidar backscatter [m sr]^-1
  Real*8, Intent(Out)                :: mpl_ext          ! mpl lidar extinction coefficient [m]^-1
  !-- Output for Doppler spectra
  Real*8, Intent(Out)                :: spectra_bins(nfft)  ! spectra velocity bin [m/s]
  Real*8, Intent(Out)                :: zhh_spectra(nfft),zvh_spectra(nfft),zvv_spectra(nfft) !output 1D Doppler spectra [mm^6 m^-3 / (m/s)]
  !--
  !
  integer                           :: nd
  real*8,dimension(:),Allocatable   :: NN    ! 1/m^3
  real*8,dimension(:),Allocatable   :: fvel  ! m/s
  real*8,dimension(:),Allocatable   :: diam,ddiam  ! m
  real*8,dimension(:),Allocatable   :: rho   ! kg/m^3
  !
  real*8                            :: dmin,dmax ! um
  real*8                            :: rrho,fr   
  real*8,dimension(:),Allocatable   :: zhh_d,zvh_d,zvv_d ! mm^6, used for spectrum generation
  !
    Call hydro_info_P3(isc,dmin,dmax,nd) ! number of diams, min,max diam for which the scatt info are stored
    !
    Allocate(NN(nd),fvel(nd),diam(nd),rho(nd))
    NN=0.d0; fvel=0.d0; diam=0.d0; rho=0.d0
    !
    !-calculate fr: rimed fraction
    if ((qir==0.d0) .or. (qhydro==0.d0)) then 
      fr = 0.d0
    else
      fr = qir/qhydro 
    endif
    !
    !-calculate bulk densty
    if ((qir==0.d0) .or. (qib==0.d0)) then
      rrho = 0.d0
    else
      rrho = qir/qib  !kg/m3
    endif
    !
    lookup_file_1=Trim(conf%file_ice_parameter)
    !
    if (isc==1) then
      call cloud_P3(nd,dmin,dmax,temp,rho_d,qhydro,qnhydro,diam,rho,NN,fvel)
    elseif (isc==2) then
      call rain_P3(nd,dmin,dmax,temp,rho_d,qhydro,qnhydro,diam,rho,NN,fvel)
    elseif (isc>=3) then
      call ice_P3(isc,nd,dmin,dmax,temp,rho_d,qhydro,qnhydro,qir,fr,rrho,diam,rho,NN,fvel,lookup_file_1)
    endif
    !
    Allocate(zhh_d(nd),zvh_d(nd),zvv_d(nd)) ! used for Doppler spectra simulation
    zhh_d=0.d0 ; zvh_d=0.d0 ; zvv_d=0.d0 ; 
    !
    IF (isc>3) THEN
      ! density of hydrometeor change with size
      call GetPolarimetricInfofromLUT_P3(isc,conf,elev,ww,temp,nd,diam,NN,rho,fvel,Zhh,Zvv,Zvh,RHOhvc,&
                                         DVh,dDVh,Dopp,Kdp,Adp,Ah,Av,diff_back_phase,zhh_d,zvh_d,zvv_d)   
    ELSE ! density of hydrometeor doesn't change with size
      call GetPolarimetricInfofromLUT_cdws(isc,conf,elev,ww,temp,nd,diam,NN,rho(1),fvel,Zhh,Zvv,Zvh,RHOhvc,&
                                         DVh,dDVh,Dopp,Kdp,Adp,Ah,Av,diff_back_phase,zhh_d,zvh_d,zvv_d)
    ENDIF
     
    !=======================================================!
    !=======================================================!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !-- Simulate Doppler spectrum---------------------------!
    !call spectrum subroutine for P3
    spectra_bins = 0.d0 ; zhh_spectra = 0.d0 ; zvh_spectra = 0.d0 ; zvv_spectra = 0.d0 
    if(conf%spectraID==1)then
      Allocate(ddiam(nd+1))
      ddiam = -999.d0 !!!diameter bin sizes [m]
      call processing_ds(isc,conf,spectra_VNyquist,spectra_NOISE_1km,NFFT,spectra_Nave,&
                        &range_m,elev,w_r,sw_dyn,nd,NN,diam,ddiam,fvel,zhh_d,zvh_d,zvv_d,&
                        &spectra_bins,zhh_spectra,zvh_spectra,zvv_spectra)
      Deallocate(ddiam)
    end if
    Deallocate(zhh_d,zvh_d,zvv_d)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !=======================================================!
     
    !---introduce cloud ceilo lidar measurements---
    ceilo_back_true=0.d0 ; ceilo_ext=0.d0
    if ( (isc==1) .and. (conf%ceiloID==1)) then
      !oue PRELOAD P3
#ifdef __PRELOAD_LUT__ 
      call GetCloudLidarMeasurements_preload(nd,diam,NN,ceilo_back_true,ceilo_ext)
#else 
      call GetCloudLidarMeasurements(conf,nd,diam,NN,ceilo_back_true,ceilo_ext)
#endif
    endif
    ! --- introduce MPL measurements for liquid and ice clouds ----
    mpl_back_true=0.d0 ; mpl_ext=0.d0
    IF  (conf%mplID>0) THEN
      if ( (isc==1) .or. (isc==3) ) then 
        call GetCloud_MPL_Measurements(isc,conf,nd,diam,NN,mpl_back_true,mpl_ext)
      endif
    ENDIF
    !
    Deallocate(NN,fvel,diam,rho)
    !
  return
  end subroutine processing_P3
  !
  !--------------------------------------------------------------------------------
  subroutine processing_sbm(isc,conf,elev,ww,temp,nd,&
                            NN,diam,rho,fvel,&
                            spectra_VNyquist,spectra_NOISE_1km,NFFT,spectra_Nave,&
                            range_m,w_r,sw_dyn,ddiam,&
                            Zhh,Zvv,Zvh,RHOhvc,DVh,dDVh,Dopp,Kdp,Adp,Ah,Av,diff_back_phase,&
                            ceilo_back_true,ceilo_ext,mpl_back_true,mpl_ext,&
                            spectra_bins,zhh_spectra,zvh_spectra,zvv_spectra)
  Use crsim_mod
  Implicit None
  !
  Integer,Intent(In)                 :: isc
  Type(conf_var),Intent(in)          :: conf
  Real*8, Intent(In)                 :: elev
  Real*8, Intent(In)                 :: ww ![m/s] vertical air velocity
  Real*8, Intent(In)                 :: temp
  Integer,Intent(In)                 :: nd
  Real*8, Intent(In)                 :: NN(nd)    !  1/m^3
  Real*8, Intent(In)                 :: diam(nd)  !  m
  Real*8, Intent(In)                 :: rho(nd)   !  kg/m^3
  Real*8, Intent(In)                 :: fvel(nd)  !  m/s
  !-- input for Doppler spectra
  Real*8, Intent(In)                 :: spectra_VNyquist,spectra_NOISE_1km
  Integer,Intent(In)                 :: spectra_Nave
  Integer,Intent(In)                 :: NFFT ! size of spectra_bins. If spectraID/=1, this is 1
  Real*8, Intent(In)                 :: range_m ! distance from radar [m], used for Doppler spectra simulation
  Real*8, Intent(In)                 :: w_r     ! radial component of wind field [m/s], used for Doppler spectra 
  Real*8, Intent(In)                 :: sw_dyn  ! spectra broadening due to dynamics [m/s], used for Doppler spectra 
  Real*8, Intent(In)                 :: ddiam(nd+1) ! diameter bin boundaries
  !--
  Real*8, Intent(Out)                :: Zhh,Zvv,Zvh,RHOhvc ! mm^6/m^3
  Real*8, Intent(Out)                :: DVh   ! mm^6/m^3 m/s  assuming ww=0 m/s
  Real*8, Intent(Out)                :: dDVh  ! mm^6/m^3 (m/s)^2 assuming ww=0 m/s
  Real*8, Intent(Out)                :: Dopp  ! mm^6/m^3 m/s  ww from WRF
  Real*8, Intent(Out)                :: Kdp   ! deg/km
  Real*8, Intent(Out)                :: Adp   ! dB/km
  Real*8, Intent(Out)                :: Ah    ! dB/km
  Real*8, Intent(Out)                :: Av    ! dB/km
  Real*8, Intent(Out)                :: diff_back_phase ! deg
  Real*8, Intent(Out)                :: ceilo_back_true ! ceilo true (unatenuated) lidar backscatter [m sr]^-1
  Real*8, Intent(Out)                :: ceilo_ext       ! ceilo lidar extinction coefficient [m]^-1
  Real*8, Intent(Out)                :: mpl_back_true ! mpl true (unatenuated) lidar backscatter [m sr]^-1
  Real*8, Intent(Out)                :: mpl_ext       ! mpl lidar extinction coefficient [m]^-1
  !-- Output for Doppler spectra
  Real*8, Intent(Out)                :: spectra_bins(nfft)  ! spectra velocity bin [m/s]
  Real*8, Intent(Out)                :: zhh_spectra(nfft),zvh_spectra(nfft),zvv_spectra(nfft) !output 1D Doppler spectra [mm^6 m^-3 / (m/s)]
  !--
  real*8,dimension(:),Allocatable   :: zhh_d,zvh_d,zvv_d ! mm^6, used for spectrum generation
    !
    ! for each grid point
    !------------------------------------------------------------------------------------
    Allocate(zhh_d(nd),zvh_d(nd),zvv_d(nd)) ! used for Doppler spectra simulation
    zhh_d=0.d0 ; zvh_d=0.d0 ; zvv_d=0.d0 ; 
    !
    call GetPolarimetricInfofromLUT_sbm(isc,conf,elev,ww,temp,nd,diam,NN,rho,fvel,Zhh,Zvv,Zvh,RHOhvc,&
                                        DVh,dDVh,Dopp,Kdp,Adp,Ah,Av,diff_back_phase,zhh_d,zvh_d,zvv_d)
     
    !=======================================================!
    !=======================================================!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !-- Simulate Doppler spectrum---------------------------!
    !call spectrum subroutine for bin model
    spectra_bins = 0.d0 ; zhh_spectra = 0.d0 ; zvh_spectra = 0.d0 ; zvv_spectra = 0.d0 
    if(conf%spectraID==1)then
      call processing_ds(isc,conf,spectra_VNyquist,spectra_NOISE_1km,NFFT,spectra_Nave,&
                        &range_m,elev,w_r,sw_dyn,nd,NN,diam,ddiam,fvel,zhh_d,zvh_d,zvv_d,&
                        &spectra_bins,zhh_spectra,zvh_spectra,zvv_spectra)
    end if
    Deallocate(zhh_d,zvh_d,zvv_d)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !=======================================================!
     
    ! --- introduce CEILO cloud lidar measurements ----
    ceilo_back_true=0.d0 ; ceilo_ext=0.d0
#ifdef __PRELOAD_LUT__
    if ( ((isc==1).or.(isc==7)) .and. (conf%ceiloID==1)) &
      call GetCloudLidarMeasurements_preload(nd,diam,NN,ceilo_back_true,ceilo_ext)
#else
    if ( ((isc==1).or.(isc==7)) .and. (conf%ceiloID==1)) &
      call GetCloudLidarMeasurements(conf,nd,diam,NN,ceilo_back_true,ceilo_ext)
#endif
    !---------------------------------------------
    ! --- introduce MPL measurements for liquid and ice clouds ----
    mpl_back_true=0.d0 ; mpl_ext=0.d0
    IF (conf%mplID>0) THEN
      if ( (isc==1) .or. (isc==3).or.(isc==7)) then
        call GetCloud_MPL_Measurements_sbm(isc,conf,nd,diam,rho,NN,mpl_back_true,mpl_ext)
      endif
    ENDIF
    !----------------------------------------------
    ! 
  return
  end subroutine processing_sbm
  !   
  subroutine hydro_info(isc,dmin,dmax,nd)
  Implicit None
  Integer,Intent(In)                 :: isc
  Real*8,Intent(Out)                 :: dmin,dmax ! um
  integer,Intent(Out)                :: nd
  !
  Real*8                             :: ddiam ! um
    !
    ! -dmin,dmax,ddiam are the those from luts
    !
    if((isc==1) .or. (isc==7)) then ! cloud&drizzle modified by oue 2017/07/21 for RAMS
      dmin=1.d0 ; dmax=250.d0 ; ddiam=1.d0
      !nd=nint((dmax-dmin)/ddiam) + 1
      nd=250
    endif
    !
    if(isc==2) then ! rain
      dmin=100.d0 ; dmax=9000.d0 ; ddiam=20.d0
      !nd=nint((dmax-dmin)/ddiam) + 1
      nd=446
    endif
    !
    if(isc==3) then ! ice
      dmin=1.d0 ; dmax=1496.d0 ; ddiam=5.d0
      !nd=nint((dmax-dmin)/ddiam) + 1
      nd=300 
    endif
    !
    if((isc==4) .or. (isc==8)) then ! snow and aggregates modified by oue 2017/07/21 for RAMS
      dmin=100.d0 ; dmax=50000.d0; ddiam=100.d0
      !nd=nint((dmax-dmin)/ddiam) + 1
      nd=500
    endif
    !
    if(isc==5) then ! graupel
      dmin=5.d0 ; dmax=50005.d0 ; ddiam=100.d0
      !nd=nint((dmax-dmin)/ddiam) + 1
      nd=501
    endif
    !
    if(isc==6) then ! hail
      dmin=5.d0 ; dmax=50005.d0 ; ddiam=100.d0
      !nd=nint((dmax-dmin)/ddiam) + 1
      nd=501
    endif
    !
  return
  end subroutine hydro_info
  !
  !----------------------------------
  ! subroutine hydro_info_P3 is added by DW
  subroutine hydro_info_P3(isc,dmin,dmax,nd)
  !
  Implicit None
  Integer,Intent(In)                 :: isc
  Real*8,Intent(Out)                 :: dmin,dmax ! um
  integer,Intent(Out)                :: nd
  !
  Real*8                             :: ddiam ! um
    !
    ! -dmin,dmax,ddiam are the those from luts
    !
    if(isc==1) then ! cloud
      dmin=1.d0 ; dmax=250.d0 ; ddiam=1.d0
      nd=250
    endif
    !
    if(isc==2) then ! rain
      dmin=100.d0 ; dmax=9000.d0 ; ddiam=20.d0
      nd=446
    endif
    !
    if(isc==3) then ! small ice 
      dmin=2.d0 ; dmax=1502.d0 ; ddiam=5.d0
      nd=301
    endif
    !
    if(isc==4) then ! unrimed ice  
      dmin=2.d0 ; dmax=20002.d0 ; ddiam=100.d0
      nd=201
    endif
    !
    if(isc==5) then ! graupel  
      dmin=2.d0 ; dmax=50002.d0 ; ddiam=100.d0  
      nd=501
    endif
    !
    if(isc==6) then ! partially rimed ice   
      dmin=2.d0 ; dmax=50002.d0 ; ddiam=100.d0
      nd=501
    endif

  !   
  return
  end subroutine hydro_info_P3
  !
  !--------------------------
  subroutine hydro_info_jfan(isc,nkr,ihtf,id1,id2,nd)
  Implicit None
  Integer,Intent(In)                 :: isc     ! scattering type 
  integer,Intent(Out)                :: nkr     ! total number of bins in each category =33
  integer,Intent(Out)                :: ihtf    ! the hydro category 1=cloud+rain,2=ice+snow, 3 graupel
  integer,Intent(Out)                :: id1,id2 ! starting and ending bin indices for each scatt. type    
  integer,Intent(Out)                :: nd      ! number of diameters
  !
  Integer                            :: ict    ! bin number to devide droplets and rain, and also ice and snow
                                               ! here ict is the last bin number in droplet category, and also 
                                               ! in ice category, and 100 um is the threshold size
   
    !
    !-----------------------------
    ! values from model simulation
    ict=17
    nkr=33
    !-----------------------------
    
    if (isc==1) then ! cloud
      id1=1 ; id2=ict
      nd=id2-id1+1
      ihtf=1 
    endif
    !
    if (isc==2) then ! rain
      id1=ict+1 ; id2=nkr  
      nd=id2-id1+1
      ihtf=1
    endif
    !
    if (isc==3) then ! ice
         id1=1 ; id2=ict
         nd=id2-id1+1
         ihtf=2
     endif
     !
     if (isc==4) then ! snow
       id1=ict+1 ; id2=nkr
       nd=id2-id1+1
       ihtf=2
     endif
     !
     if (isc==5) then ! graupel
       id1=1 ; id2=nkr
       nd=id2-id1+1
       ihtf=3
     endif 
     !
  return
  end subroutine hydro_info_jfan
  !
  !--------------------------
  subroutine hydro_info_samsbm(isc,nkr,ihtf,id1,id2,nd)
  Implicit None
  Integer,Intent(In)                 :: isc     ! scattering type 1: cloud, 2: rain
  integer,Intent(Out)                :: nkr     ! total number of bins in each category =33
  integer,Intent(Out)                :: ihtf    ! the hydro category 1=cloud+rain
  integer,Intent(Out)                :: id1,id2 ! starting and ending bin indices for each scatt. type    
  integer,Intent(Out)                :: nd      ! number of diameters
  !
  Integer                            :: ict    ! bin number to devide droplets and rain, and also ice and snow
                                               ! here ict is the last bin number in droplet category, and also 
                                               ! in ice category, and 100 um is the threshold size
   
    !
    !-----------------------------
    ! values from model simulation
    ict=12 ! index of the last cloud bin
    nkr=36 ! index of the last rain bin
    !-----------------------------
   
    if (isc==1) then ! cloud
      id1=1 ; id2=ict
      nd=id2-id1+1
      ihtf=1 
    endif
    !
    if (isc==2) then ! rain
      id1=ict+1 ; id2=nkr  
      nd=id2-id1+1
      ihtf=1
    endif
     !
  return
  end subroutine hydro_info_samsbm
  !
  !subroutine cloud_morrison(nd,dmin,dmax,temp,rho_d,qhydro,qnhydro,diam,rho,NN,fvel)
  subroutine cloud_morrison(nd,dmin,dmax,temp,rho_d,qhydro,qnhydro,Nconst,diam,rho,NN,fvel) !added Nconst, Dec. 11, 2021
  use crsim_mod
  Use phys_param_mod, ONLY :pi,grav, rhow, rhow_25C,oneOverThree
  Implicit None
  !
  integer,Intent(in)      :: nd
  real*8,Intent(in)       :: dmin,dmax ! um
  real*8,Intent(in)       :: temp      ! C
  real*8,Intent(in)       :: rho_d     ! kg/m^3
  real*8,Intent(in)       :: qhydro    ! kg/kg
  Real*8,Intent(in)       :: qnhydro   ! 1/kg
  Real*8,Intent(in)       :: Nconst    ! 1/cm^3 !added Nconst, Dec. 11, 2021
  real*8,Intent(out)      :: diam(nd)  ! m
  real*8,Intent(out)      :: rho(nd)   ! kg/m^3
  real*8,Intent(Out)      :: NN(nd)    ! 1/m^3
  real*8,Intent(Out)      :: fvel(nd)  ! m/s
  !
  real*8                           :: rho_surf ! kg/m^3
  real*8                           :: visc ! Viscosity (kg/m/s)
!  real*8,parameter                 :: Nconst=250.d0 ! 1/cm^3 ! commented Dec. 11, 2021
  real*8                           :: av,bv
  real*8                           :: mu,Nc
  real*8,Dimension(:),Allocatable  :: lambda,N0
  !
  integer                            :: ir
  real*8                             :: ddiam  ! m
    ! --------------------------------------
    do ir=1,nd
      diam(ir)=dmin*1.d-4+dble(ir-1)*(dmax-dmin)/dble(nd-1)*1.d-4 ! diameter in cm
    enddo
    ddiam=diam(10)-diam(9) ! dr in cm
    !
    ddiam=ddiam*1.d-2 ! dr in m
    diam=diam*1.d-2 ! diameter in m 
    !---------------------------------------
    rho_surf=85000.d0/(287.04d0*273.15d0)   
    visc=1.72d0 *1.d-5 *(393.d0/(temp+393.d0))*((temp+273.d0)/273.d0)
    ! parametrs for fall velocity 
    av = grav*rhow/(18.d0*visc)  ! 
    bv=2.d0
    !-----------------------------
    if (qnhydro>0) then
     Nc=qnhydro
    else
     Nc=(Nconst*1.d+6)/rho_d ! now in 1/kg
    endif
    !-----------------------------
    !
    mu=0.0005714d0*  (Nc  * 1.d-6 * rho_d)  +  0.2714d0
    mu=1.d0/(mu*mu) - 1d0
    mu=MAX(mu,2.d0)
    mu=MIN(mu,10.d0)
    !
    Allocate(lambda(nd),N0(nd))
    !
    lambda=(pi*rhow_25C*Nc)/(6.d0*qhydro) * dgamma(mu+4.d0)/dgamma(mu+1.d0)
    lambda=lambda**oneOverThree
    !
    N0=(Nc * lambda**(mu+1.d0))/(dgamma(mu+1.d0))
    !
    do ir=1,nd
      fvel(ir) = av * (diam(ir))**bv  ! m/s
      NN(ir)=(N0(ir) *(diam(ir))**mu)  * dexp(-lambda(ir)*diam(ir))
      NN(ir)=NN(ir)*ddiam   ! in m-3
    enddo
    !
    Deallocate(lambda,N0)
    !
  return
  endsubroutine cloud_morrison
  !
  !
  subroutine hydro_morrison(mpphys,isc,nd,dmin,dmax,temp,rho_d,qhydro,qnhydro,diam,rho,NN,fvel)
  use crsim_mod
  Use phys_param_mod, ONLY: pi,oneOverThree, Rd,T0K,p0_850
  Implicit None
  !
  integer, intent(In)     :: isc,mpphys
  integer,Intent(in)      :: nd
  real*8,Intent(in)       :: dmin,dmax ! um
  real*8,Intent(in)       :: temp      ! C
  real*8,Intent(in)       :: rho_d     ! kg/m^3
  real*8,Intent(in)       :: qhydro    ! kg/kg
  real*8,Intent(in)       :: qnhydro   ! 1/kg 
  real*8,Intent(Out)      :: diam(nd)  ! m
  real*8,Intent(Out)      :: rho(nd)   ! kg/m^3
  real*8,Intent(Out)      :: NN(nd)    ! 1/m^3
  real*8,Intent(Out)      :: fvel(nd)  ! m/s
  !
  real*8                               :: rho_surf ! kg/m^3
  real*8                               :: av,bv,k,fc,fck
  real*8,Dimension(:),Allocatable      :: lambda,N0
  !
  integer                              :: ir
  real*8                               :: ddiam
    !
    ! --------------------------------------
    rho_surf=p0_850/(Rd*T0K)   ! 
    fc=rho_surf/rho_d
    !  --------------------------------------
    !
    if(isc==2) then ! rain
      rho=997.d0 
      k=0.54d0 ; av=841.99667d0 ; bv=0.8d0
    endif
    !
    if(isc==3) then ! ice
      rho=500.d0
      k=0.35d0 ; av=700.d0 ; bv=1.d0
    endif
    !
    if(isc==4) then ! snow
      rho=100.d0 
      k=0.54d0 ; av=11.72d0 ; bv=0.41d0
    endif
    !
    if(mpphys == 0)then
      if(isc==5) then ! graupel
        rho=400.d0
        k=0.54d0 ; av=19.3d0 ; bv=0.37d0
      endif
    endif
    if(mpphys == 1) then
      if(isc==5) then ! hail
        rho=900.d0
        k=0.5d0 ; av=206.89d0 ; bv=0.6384d0
      endif
    endif
    !
    !---------------------------------------
    !--------------------------------------
    !
    do ir=1,nd
      diam(ir)=dmin*1.d-4+dble(ir-1)*(dmax-dmin)/dble(nd-1)*1.d-4 ! diameter in cm
    enddo
    ddiam=diam(10)-diam(9) ! ddiam in cm
    !
    ddiam=ddiam*1.d-2 ! ddiam in m
    diam=diam*1.d-2 ! diameter in m 
    !
    !---------------------------------------
    !
    !
    Allocate(lambda(nd),N0(nd))
    !
    lambda=(pi * rho * qnhydro/qhydro)**oneOverThree ! 1/m
    N0=qnhydro *lambda * rho_d  ! 1/m^4
    !
    !
    fck=fc**k
    do ir=1,nd
      fvel(ir) = av * (diam(ir))**bv * fck  ! m/s
      NN(ir)=N0(ir) * dexp(-lambda(ir)*diam(ir))
      NN(ir)=NN(ir)*ddiam   ! in m-3
    enddo
    !
    Deallocate(lambda,N0)
    !
    go to 234
      if (isc==2) write(*,*) 'rain',qhydro,qnhydro
      if (isc==3) write(*,*) 'ice',qhydro,qnhydro
      if (isc==4) write(*,*) 'snow',qhydro,qnhydro
      if (isc==5) write(*,*) 'graup',qhydro,qnhydro
      ! computed wcont/model_wcont
      if (isc==2) write(*,*) 'rain',(pi/6.d0*Sum(rho*NN*diam*diam*diam))/ (qhydro*rho_d)
      if (isc==3) write(*,*) 'ice',(pi/6.d0*Sum(rho*NN*diam*diam*diam))/ (qhydro*rho_d)
      if (isc==4) write(*,*) 'snow',(pi/6.d0*Sum(rho*NN*diam*diam*diam))/ (qhydro*rho_d)
      if (isc==5) write(*,*) 'graup',(pi/6.d0*Sum(rho*NN*diam*diam*diam))/ (qhydro*rho_d)
      ! computed NN/NN_model
      if (isc==2) write(*,*) 'rain',Sum(NN)/(qnhydro*rho_d)
      if (isc==3) write(*,*) 'ice', Sum(NN)/(qnhydro*rho_d)
      if (isc==4) write(*,*) 'snow',Sum(NN)/(qnhydro*rho_d)
      if (isc==5) write(*,*) 'graup',Sum(NN)/(qnhydro*rho_d)
      !
      !pause
    234 continue
  !
  return
  end subroutine hydro_morrison
  !
  !
  !
  !---------------------------------------------------------------------------------------
  !!! Added by oue, 2016/09/19
  subroutine hydro_thompson(isc,nd,dmin,dmax,temp,rho_d,qhydro,qnhydro,diam,rho,NN,fvel)
  use crsim_mod
  Use phys_param_mod, ONLY: pi,oneOverThree,Rd
  Implicit None
  !-- WRFv3.7 and Thompson et al. (2008)
  integer, intent(In)     :: isc
  integer,Intent(in)      :: nd
  real*8,Intent(in)       :: dmin,dmax ! um
  real*8,Intent(in)       :: temp      ! C
  real*8,Intent(in)       :: rho_d     ! kg/m^3
  real*8,Intent(in)       :: qhydro    ! kg/kg
  real*8,Intent(in)       :: qnhydro   ! 1/kg 
  real*8,Intent(Out)      :: diam(nd)  ! m
  real*8,Intent(Out)      :: rho(nd)   ! kg/m^3
  real*8,Intent(Out)      :: NN(nd)    ! 1/m^3
  real*8,Intent(Out)      :: fvel(nd)  ! m/s
  !
  real*8                               :: rho_hydro,rho_surf ! kg/m^3
  real*8                               :: av,bv,k,fc,fck,fd
  real*8                               :: lambda,N0,mu
  !
  integer                              :: ir,j
  real*8                               :: ddiam
  !
  real*8                               :: am,bm
  real*8                               :: Nt_c !/cm^3
  real*8                               :: rc ! water content kg/m^3
  real*8                               :: qn ! number concentration, identical to qnhydro
  real*8                               :: gammln
  ! Parameters used for cloud and graupel size distributions
  real*8                               :: cce1,cce2,cce3,ocg1,ocg2,ccg1,ccg2,ccg3
  real*8                               :: X,Y,TMP,SER,STP
  ! Parameters used for graupel size distribution
  real*8                               :: N0_exp,lam_exp ! used for graupel
  ! Parameters used for cloud size distribution
  real*8, parameter                    :: N1 = 9.d9
  real*8, parameter                    :: N2 = 2.d6
  real*8, parameter                    :: qr0 = 1.d-4
  ! Parameters used for snow size distribition
  real*8, parameter                    :: Kap0 = 490.6d0
  real*8, parameter                    :: Kap1 = 17.46d0
  real*8, parameter                    :: Lam0 = 20.78d0
  real*8, parameter                    :: Lam1 = 3.29d0
  real*8                               :: Tc !temperature in C used for snow
  real*8                               :: cse1,a_,b_,loga_,M0,M2,M3,Mrat,slam1,slam2
  real*8, Dimension(10), parameter       :: &
       sa = (/ 5.065339d0, -0.062659d0, -3.032362d0, 0.029469d0, -0.000285d0, &
               0.31255d0,   0.000204d0,  0.003199d0, 0.0d0,      -0.015952d0/)
  real*8, Dimension(10), parameter       :: &
       sb = (/ 0.476221d0, -0.015896d0,  0.165977d0, 0.007468d0, -0.000141d0, &
               0.060366d0,  0.000079d0,  0.000594d0, 0.0d0,      -0.003577d0/)
  double precision, dimension(6), parameter:: &
       COF = (/76.18009172947146D0, -86.50532032941677D0, &
               24.01409824083091D0, -1.231739572450155D0, &
                 .1208650973866179D-2, -.5395239384953D-5/)
  !  parameters for graupel size distribution
  real*8  :: mvd_r,rg,N0_min,xslw1,ygra1,zans1
  LOGICAL :: L_qr
 

    ! init 
    qn=0.d0
    am=0.d0
    rc=0.d0
    rg=0.d0
    ! --------------------------------------
    rho_surf=101325.d0/(Rd*298.d0)
    fc=rho_surf/rho_d
    !  --------------------------------------
    !
    !
    if (isc==1) then ! cloud
      rho_hydro=1000.d0
      k=0.5d0 ; av=0.316946d8 ; bv=2.0d0 ; fd=0.d0 ! coefficients for fall velocity
      Nt_c = 100.d6 ! in m^-3 = 100 cm^-3
      mu = MIN(15.d0, (1000.d6/Nt_c + 2.d0))
      am = pi*rho_hydro/6.d0 ;   bm=3.d0
      rc = qhydro * rho_d ! in kg m^-3
      qn = Nt_c * 1.d6 / rho_d !in kg^-1
    endif
    !
    if(isc==2) then ! rain
      rho_hydro=1000.d0 
      k=0.5d0 ; av=4854.0d0 ; bv=1.0d0 ; fd=195.d0
      mu = 0.d0
      am =pi*rho_hydro/6.d0 ;   bm=3.d0
      qn = qnhydro
    endif
    !
    if(isc==3) then ! ice
      rho_hydro=890.d0
      k=0.5d0 ; av=1847.5d0 ; bv=1.d0 ; fd=0.d0
      mu =0.d0
      am = pi*rho_hydro/6.d0 ;   bm=3.d0
      qn = qnhydro
    endif
    !
    if(isc==4) then ! snow
      rho_hydro=100.d0 
      k=0.5d0 ; av=40.d0 ; bv=0.55d0 ; fd=100.d0
      mu = 0.6357d0
      am = 0.069d0 ;   bm = 2.d0
      rc = qhydro * rho_d   
    endif
    !
    if(isc==5) then ! graupel
      rho_hydro=500.d0
      k=0.5d0 ; av=442.d0 ; bv=0.89d0 ; fd=0.d0
      mu =0.d0
      am =pi*rho_hydro/6.d0 ;   bm=3.d0
      rc = qhydro * rho_d !fixed 2017.05
    endif
    !
    !--------------------------------------
    ! Diameter
    do ir=1,nd
      diam(ir)=dmin*1.d-4+dble(ir-1)*(dmax-dmin)/dble(nd-1)*1.d-4 ! diameter in cm
    enddo
    ddiam=diam(10)-diam(9) ! ddiam in cm
    !
    ddiam=ddiam*1.d-2 ! ddiam in m
    diam=diam*1.d-2 ! diameter in m 
    !
    !---------------------------------------
    !
    rho=rho_hydro
    if(isc==4) then ! snow
      do ir=1,nd
        rho(ir)=0.069d0 * 6.d0 / pi / diam(ir)
      enddo
    endif
    !
    !---------------------------------------
    ! Fall velocity
    fck=fc**k
    do ir=1,nd
      fvel(ir) = av * (diam(ir))**bv * fck * dexp((-1.d0)*fd*diam(ir))  ! m/s
    enddo
    !---------------------------------------
    !
    !
    !lambda=(pi * rho_hydro * qn/qhydro)**bb ! 1/m
    !gamma_alphap1= dgamma(mu+1)
    !gfac=dgamma(mu+bm) / gamma_alphap1
    !lambda=(gfac * (am*qn)/(qhydro))**bb  ! 1/m
  
  
    if(isc==1) then !cloud
      cce1=mu+1.d0
      cce2=bm+mu+1.d0
      !-- get ccg1 from cce1
      X=cce1
      Y=X
      TMP=X+5.5D0
      TMP=(X+0.5D0)*dlog(TMP)-TMP
      SER=1.000000000190015D0
      STP = 2.5066282746310005D0
      do j=1,6
        Y=Y+1.D0
        SER=SER+COF(j)/Y
      enddo
      gammln =TMP+DLOG(STP*SER/X)
      ccg1 = dexp(gammln)
      ocg1=1.d0/ccg1
      !--       
      !-- get ccg2 from cce2
      X=cce2
      Y=X
      TMP=X+5.5D0
      TMP=(X+0.5D0)*dlog(TMP)-TMP
      SER=1.000000000190015D0
      STP = 2.5066282746310005D0
      do j=1,6
        Y=Y+1.D0
        SER=SER+COF(j)/Y
      enddo
      gammln =TMP+DLOG(STP*SER/X)
      ccg2 = dexp(gammln)
      !--
      lambda = 1.0D-6 * (Nt_c*am* ccg2 * ocg1 / rc)**oneOverThree !/m^1
      N0 = 1.0D-18 * Nt_c*ocg1 * lambda**cce1
      do ir=1,nd
        NN(ir) = N0* (diam(ir)*1.0D6)**mu * DEXP(-lambda*(diam(ir)*1.0D6))*ddiam * 1.0D24
      enddo
    endif   
    !
    if((isc==2) .or. (isc==3)) then !rain, ice
      cce2 = mu+1.d0
      cce3 = bm+mu+1.d0
      !-- get ccg2 from cce2
      X=cce2
      Y=X
      TMP=X+5.5D0
      TMP=(X+0.5D0)*dlog(TMP)-TMP
      SER=1.000000000190015D0
      STP = 2.5066282746310005D0
      do j=1,6
        Y=Y+1.D0
        SER=SER+COF(j)/Y
      end do
      gammln =TMP+dlog(STP*SER/X)
      ccg2 = dexp(gammln)
      !-- get ccg3 from cce3
      X=cce3
      Y=X
      TMP=X+5.5D0
      TMP=(X+0.5D0)*dlog(TMP)-TMP
      SER=1.000000000190015D0
      STP = 2.5066282746310005D0
      do j=1,6
        Y=Y+1.D0
        SER=SER+COF(j)/Y
      end do
      gammln =TMP+dlog(STP*SER/X)
      ccg3 = dexp(gammln)
      !--
      lambda = (am*ccg3/ccg2 * qn/qhydro)**oneOverThree
      N0 = qn*rho_d/ccg2 * lambda**cce2 !fixed: rho_d 2017.05.22
      do ir=1,nd
        NN(ir) = N0*diam(ir)**mu *DEXP(-lambda*diam(ir))*ddiam
      enddo
    endif
    !
    if(isc==4) then !snow
      cse1 = bm+1.d0
      Tc = temp
      if(Tc >0) then
        Tc = -0.01d0
      endif
      M2 = rc/am *1.0d0 !(bm)th moment (here, bm=3)
      loga_ = sa(1) + sa(2)*Tc + sa(3)*cse1 &
            + sa(4)*Tc*cse1 + sa(5)*Tc*Tc &
            + sa(6)*cse1*cse1 + sa(7)*Tc*Tc*cse1 &
            + sa(8)*Tc*cse1*cse1 + sa(9)*Tc*Tc*Tc &
            + sa(10)*cse1*cse1*cse1
      a_ = 10.d0**loga_
      b_ = sb(1)+sb(2)*Tc+sb(3)*cse1 + sb(4)*Tc*cse1 &
         + sb(5)*Tc*Tc + sb(6)*cse1*cse1 &
         + sb(7)*Tc*Tc*cse1 + sb(8)*Tc*cse1*cse1 &
         + sb(9)*Tc*Tc*Tc+sb(10)*cse1*cse1*cse1
      M3 = a_ * M2**b_ !(bm+1)th moment 
      Mrat = M2*(M2/M3)*(M2/M3)*(M2/M3)
      M0   = (M2/M3)**mu
      slam1 = M2 / M3 * Lam0
      slam2 = M2 / M3 * Lam1
      do ir=1,nd
        NN(ir) = Mrat*(Kap0*DEXP(-slam1*diam(ir)) &
               + Kap1*M0*diam(ir)**mu * DEXP(-slam2*diam(ir)))*ddiam
      enddo
    endif
    !
    if(isc==5) then !graupel
      cce1=bm+1.d0
      cce2=mu+1.d0
      cce3=bm+mu+1.d0
      !-- get ccg1 from cce1
      X=cce1
      Y=X
      TMP=X+5.5D0
      TMP=(X+0.5D0)*dlog(TMP)-TMP
      SER=1.000000000190015D0
      STP = 2.5066282746310005D0
      do j=1,6
        Y=Y+1.D0
        SER=SER+COF(j)/Y
      enddo
      gammln =TMP+dlog(STP*SER/X)
      ccg1 = dexp(gammln)
      ocg1=1.d0/ccg1
      !-- get ccg2 from cce2
      X=cce2
      Y=X
      TMP=X+5.5D0
      TMP=(X+0.5D0)*dlog(TMP)-TMP
      SER=1.000000000190015D0
      STP = 2.5066282746310005D0
      do j=1,6
        Y=Y+1.D0
        SER=SER+COF(j)/Y
      enddo
      gammln =TMP+dlog(STP*SER/X)
      ccg2 = dexp(gammln)
      ocg2=1.d0/ccg2
      !-- get ccg3 from cce3
      X=cce3
      Y=X
      TMP=X+5.5D0
      TMP=(X+0.5D0)*dlog(TMP)-TMP
      SER=1.000000000190015D0
      STP = 2.5066282746310005D0
      do j=1,6
        Y=Y+1.D0
        SER=SER+COF(j)/Y
      enddo
      gammln =TMP+dlog(STP*SER/X)
      ccg3 = dexp(gammln)
      !--

      !N0_exp=1.D6 ! in 1/m^-4
      !lam_exp = (N0_exp*am*ccg1/rc)**(1./cce1)
      !lambda = lam_exp * (ccg3*ocg2*ocg1)**oneOverThree
      !N0 = max(1.D4,min(200./qhydro,5.D6))

      !+---+-----------------------------------------------------------------+
      !Calculate y-intercept, slope values for graupel. used for Z calculation
      !in WRF
      !+---+-----------------------------------------------------------------+
      !mvd_r and L_qr depends on rain mixing ratio when rain mixing ratio > 1E-6 in WRF. 
      !However, due to the CRSIM structure, rain mixing ratio is assumed to be <1E-6
      !for graupel size distribution
      mvd_r = 50.d-6
      L_qr = .false.
      !if (qhydro .gt. 1.E-6) then
      rg = rc
      !   L_qg = .true.
      !else
      !   rg = 1.E-12
      !   L_qg = .false.
      !endif
      N0_min = 3.d6
      if (L_qr .and. mvd_r.gt.100.d-6) then
        xslw1 = 4.01d0 + dlog10(mvd_r)
      else
        xslw1 = 0.01d0
      endif
      ygra1 = 4.31d0 + dlog10(max(5.d-5, rg))
      zans1 = 3.1d0 + (100.d0/(300.d0*xslw1*ygra1/(10.d0/xslw1+1.d0+0.25d0*ygra1)+30.d0+10.d0*ygra1))
      N0_exp = 10.d0**(zans1)
      N0_exp = MAX(1.d4, MIN(N0_exp, 3.d6))
      N0_min = MIN(N0_exp, N0_min)
      N0_exp = N0_min
      lam_exp = (N0_exp*am*ccg1/rg)**(1.d0/cce1)
      lambda = lam_exp * (ccg3*ocg2*ocg1)**(1.d0/bm)
      N0 = N0_exp/(ccg2*lam_exp) * lambda**(1.d0/cce2)
      !+---+-----------------------------------------------------------------+
      !
      do ir=1,nd
        NN(ir) = N0* diam(ir)**mu * dexp(-lambda*diam(ir))*ddiam
      enddo
    endif
    !
    !
    go to 234
      if (isc==1) write(*,*) 'cloud',qhydro,qn
      if (isc==2) write(*,*) 'rain',qhydro,qnhydro
      if (isc==3) write(*,*) 'ice',qhydro,qnhydro
      if (isc==4) write(*,*) 'snow',qhydro,qnhydro
      if (isc==5) write(*,*) 'graup',qhydro,qnhydro
      ! computed wcont/model_wcont
      if (isc==1) write(*,*) 'cloud',(pi/6.d0*Sum(rho*NN*diam*diam*diam))/ (qhydro*rho_d)
      if (isc==2) write(*,*) 'rain',(pi/6.d0*Sum(rho*NN*diam*diam*diam))/ (qhydro*rho_d)
      if (isc==3) write(*,*) 'ice',(pi/6.d0*Sum(rho*NN*diam*diam*diam))/ (qhydro*rho_d)
      if (isc==4) write(*,*) 'snow',(pi/6.d0*Sum(rho*NN*diam*diam*diam))/ (qhydro*rho_d)
      if (isc==5) write(*,*) 'graup',(pi/6.d0*Sum(rho*NN*diam*diam*diam))/ (qhydro*rho_d)
      ! computed NN/NN_model
      if (isc==1) write(*,*) 'cloud',Sum(NN)/(qn*rho_d)
      if (isc==2) write(*,*) 'rain',Sum(NN)/(qnhydro*rho_d)
      if (isc==3) write(*,*) 'ice', Sum(NN)/(qnhydro*rho_d)
      if (isc==4) write(*,*) 'snow',Sum(NN)/(qnhydro*rho_d)
      if (isc==5) write(*,*) 'graup',Sum(NN)/(qnhydro*rho_d)
      !
      !pause
    234 continue
    !
  return
  end subroutine hydro_thompson
  !!
  !!-- addedd by oue
  !--------------------------------------------------------------------
  !
  subroutine hydro_milbrandt_yau(isc,snow_spherical,nd,dmin,dmax,temp,rho_d,rho_ds,qhydro,qnhydro,diam,rho,NN,fvel)
  use crsim_mod
  Use phys_param_mod, ONLY: pi
  Implicit None
  !
  integer, intent(In)     :: isc
  integer, intent(In)     :: snow_spherical ! =1 for snow spherical
  integer,Intent(in)      :: nd
  real*8,Intent(in)       :: dmin,dmax ! um
  real*8,Intent(in)       :: temp      ! C
  real*8,Intent(in)       :: rho_d     ! kg/m^3
  real*8,Intent(in)       :: rho_ds    ! kg/m^3
  real*8,Intent(in)       :: qhydro    ! kg/kg
  real*8,Intent(in)       :: qnhydro   ! 1/kg 
  real*8,Intent(Out)      :: diam(nd)  ! m
  real*8,Intent(Out)      :: rho(nd)   ! kg/m^3
  real*8,Intent(Out)      :: NN(nd)    ! 1/m^3
  real*8,Intent(Out)      :: fvel(nd)  ! m/s
  !
  real*8                               :: rho_surf ! kg/m^3
  real*8                               :: av,bv,k,fc,fck
  real*8,Dimension(:),Allocatable      :: lambda
  !
  integer                              :: ir
  real*8                               :: ddiam
  !
  real*8                               :: piov6  ! pi/6
  real*8                               :: nu ! the dispersion parameter of the gener. gamma distr. funct
  real*8                               :: alpha ! the shape parameter
  real*8                               :: alphap1 ! alpha+1
  real*8                               :: gamma_alphap1 ! gamma(alpha+1)
  real*8                               :: gfac  ! gamma(alpha +1  + bm/nu) / gamma(alpha+1)
  real*8                               :: am,bm ! coefficients in the mass-size relationship m = am * D^bm
  real*8                               :: rhoh  ! [kg/m^3]  hydrometeor bulk density 
  real*8                               :: ibm   ! 1/bm
  real*8                               :: a_rhoh, b_rhoh ! coefficients in relations rhoh = a_rhoh * D^(b_rhoh), 
                                                         ! or  rhoh = 6/pi am D^(bm-3) where rhoh is the bulk density of a sphere 
                                                         ! with the mass equivalent to the mass of a particle with m=am D^m
  
    ! generalized gamma distribution function 
      ! N(D)=Nt nu/G(1+alpha) lambda^[nu(1+alpha)]  D^[nu(1+alpha)-1]  exp[ -(lambda D)^nu]
    ! lambda=the slope; nu=the dispersion parameter; alpha the shape parameter where
    ! alpha== (mu + 1)/nu -1  and mu is the tail of the distribution
    !
  
    piov6=pi/6.d0
    a_rhoh=0.d0; b_rhoh=0.d0
    am=0.d0 ; bm=0.d0
    nu=0.d0 ; alpha=0.d0
    k=0.d0 ; av=0.d0 ; bv=0.d0
    !
    ! --------------------------------------
    rho_surf=rho_ds
    fc=rho_surf/rho_d
    !  --------------------------------------
    !!
    if(isc==1) then ! cloud
      nu=3.d0 ; alpha=1.d0 
      rhoh=1000.d0
      am=piov6*rhoh ; bm=3.d0
      k=0.d0 ; av=0.d0 ; bv=0.d0
    end if
    !
    if(isc==2) then ! rain
      nu=1.d0 ; alpha=2.d0
      rhoh=1000.d0
      am=piov6*rhoh ; bm=3.d0
      k=0.5d0 ; av=149.1d0 ; bv=0.5d0
    end if
    !
    if(isc==3) then ! ice
      nu=1.d0 ; alpha=0.d0
      rhoh=500.d0 ! not used actually in the code 
      am=440.d0 ; bm=3.d0 ! as bm=3 ice is spherical and density doesn't change with size, and is equal to 6/pi * am
      rhoh=am/piov6
      k=0.5d0 ; av=71.34d0 ; bv=0.6635d0
    endif
    !
    if(isc==4) then ! snow
      nu=1.d0 ; alpha=0.d0
      rhoh=100.d0
      if (snow_spherical==1) then
        am=piov6*rhoh ; bm=3.d0     ! for snow spherical
      else  
        am=0.1597d0 ; bm=2.078d0   ! snow not spherical  
        ! PROBABLY ERROR IN WRF CODE, am =0.1597d0 is too small,and gives rho=64.6 kg/m^3 fr D=3 mm and rho=23.4 kg/m^3 fr D=9 mm 
        !  assumed here  am=15.97, and then  rho= 646 kg/m^3 for D=3 mm and rho=234 kg/m^3 for D=9 mm.
  
        am=1.597d0  ! AT assumed correct value 
        a_rhoh=am/piov6; b_rhoh=bm-3.d0 ! rhoh = 6/pi am D^(bm-3) -> the bulk density of a sphere with the equivalent mass
                                        ! rhoh = a_rhoh * D^(b_rhoh)  
   
      endif
      k=0.5d0 ; av=11.72d0 ; bv=0.41d0
    endif
    !
    if(isc==5) then ! graupel
      nu=1.d0 ; alpha=0.d0
      rhoh=400.d0
      am=piov6*rhoh ; bm=3.d0 
      k=0.5d0 ; av=19.3d0 ; bv=0.37d0
    endif
    !
    if(isc==6) then ! hail
      nu=1.d0 ; alpha=3.d0
      rhoh=900.d0
      am=piov6*rhoh ; bm=3.d0
      k=0.5d0 ; av=206.89d0 ; bv=0.6384d0
    endif
    !
    !---------------------------------------
    !
    !--------------------------------------
    !
    do ir=1,nd
      diam(ir)=dmin*1.d-4+dble(ir-1)*(dmax-dmin)/dble(nd-1)*1.d-4 ! diameter in cm
    enddo
    ddiam=diam(10)-diam(9) ! ddiam in cm
    !
    ddiam=ddiam*1.d-2 ! ddiam in m
    diam=diam*1.d-2 ! diameter in m 
    !
    !---------------------------------------
    !
    !
    ! rho does not depend on the size so in fact lambda doesn't need to have
    ! dimension nd
    rho(:)=rhoh  ! if the hydrometeor density is the same for all sizes
    if ((isc==4).and.(snow_spherical/=1)) rho(:)=a_rhoh*diam(:)**b_rhoh ! the hydr. density changes with the size
   
    Allocate(lambda(nd))
    !
    alphap1=1.d0+alpha
    gamma_alphap1= dgamma(alphap1)
    gfac=dgamma(alphap1+bm/nu) / gamma_alphap1
    !
    ibm=1.d0/bm
    lambda=(gfac * (am*qnhydro)/(qhydro))**ibm  ! 1/m
    !
    fck=fc**k
    do ir=1,nd
      NN(ir)=qnhydro*rho_d * nu/gamma_alphap1 * lambda(ir)**(nu*alphap1) * &
             diam(ir)**(nu*alphap1-1.d0)*dexp(-(lambda(ir)*diam(ir))**nu) ! 1/m^4
      NN(ir)=NN(ir)*ddiam   ! in m-3
      fvel(ir) = av * (diam(ir))**bv * fck  ! m/s
    enddo
  
    Deallocate(lambda)
    !
  return
  end subroutine  hydro_milbrandt_yau
  !!
  !!
  !---------------------------------------------------------------------------------------
  !!! Added by oue, 2017/07/17 for ICON
  subroutine hydro_icon(isc,nd,dmin,dmax,temp,rho_d,qhydro,qnhydro,diam,rho,NN,fvel)
  use crsim_mod
  use phys_param_mod, ONLY: pi
  Implicit None
  !-- ICON two moment, Seifert and Beheng (2006)
  integer, intent(In)     :: isc
  integer,Intent(in)      :: nd
  real*8,Intent(in)       :: dmin,dmax ! um
  real*8,Intent(in)       :: temp      ! C
  real*8,Intent(in)       :: rho_d     ! kg/m^3
  real*8,Intent(in)       :: qhydro    ! kg/kg
  real*8,Intent(in)       :: qnhydro   ! 1/kg 
  real*8,Intent(Out)      :: diam(nd)  ! m
  real*8,Intent(Out)      :: rho(nd)   ! kg/m^3
  real*8,Intent(Out)      :: NN(nd)    ! 1/m^3
  real*8,Intent(Out)      :: fvel(nd)  ! m/s
  !
  real*8                               :: rho_surf,rho_hydro ! kg/m^3
  real*8                               :: av,bv,k,fc !fall velocity coefficients
  real*8                               :: ag,bg !D-M coefficients
  real*8                               :: lambda,N0,mu,nu,mue
  real*8                               :: Dm !mean volime diameter [m]
  real*8                               :: Mm !mean particle mass [kg]
  integer                              :: ir
  real*8                               :: ddiam
  real*8,Allocatable,Dimension(:)      :: mass,dmass
  !
  real*8                               :: am,bm,ibm
  real*8                               :: rc ! water content kg/m^3
  real*8                               :: qn ! number concentration, identical to qnhydro
  real*8                               :: gamma_alphap1,gfac,gam_lam,gam_n0
  !  particle_rain_coeffs
  real*8, parameter                     ::        alfa = 9.292000d0
  real*8, parameter                     ::        beta = 9.623000d0
  real*8, parameter                     ::        gama = 6.222d+2
  real*8, parameter                     ::        cmu0 = 6.000000d0
  real*8, parameter                     ::        cmu1 = 3.000d+1
  real*8, parameter                     ::        cmu2 = 1.000d+3
  real*8, parameter                     ::        cmu3 = 1.100d-3 != D_br
  real*8, parameter                     ::        cmu4 = 1.000000d0
  real*8, parameter                     ::        cmu5 = 2.d0
    !
    Allocate( mass(nd),dmass(nd))
    !
    ! --------------------------------------
    rho_surf=1.225d0 
    fc=rho_surf/rho_d
    !  --------------------------------------
    !
    if(isc==1) then ! cloud
      av=3.75d5 ; bv=2.0d0/3.0d0 ; k = 1.0d0 ! coefficients for fall velocity
      ag=1.24d-01 ; bg=0.333333d0  ! coefficients for D-Mass relationship
      mu = 1.d0 ; nu = 1.d0;
    endif
    !
    if(isc==2) then ! rain
      av=159.0d0 ; bv=0.266d0 ; k=0.5d0
      ag=1.24d-01 ; bg=0.333333d0  ! coefficients for D-Mass relationship
      mu = 1.d0/3.d0 ; nu = 0.d0 !nu = (-1) * 2.d0/3.d0
    endif
    !
    if(isc==3) then ! ice
      av=317.0d0 ; bv=0.363d0 ; k=0.5d0
      ag=0.835d0 ; bg=0.39d0  ! coefficients for D-Mass relationship
     !ag=0.217d0 ; bg=0.302d0  ! coefficients for D-Mass relationship
     mu =1.d0/3.d0 ;  nu = 0.d0 !nu = 1.d0
    endif
    !
    if(isc==4) then ! snow
      av=27.7d0 ; bv=0.216d0 ; k=0.5d0
      ag=5.13d0 ; bg=0.5d0  ! coefficients for D-Mass relationship
      !ag=8.156d0 ; bg=0.526d0  ! coefficients for D-Mass relationship
        mu =0.5d0 ;  nu = 0.d0 !nu = 1.d0
    endif
    !
    if(isc==5) then ! graupel
      av=40.d0 ; bv=0.23d0 ; k=0.5d0
      ag=1.42d-01 ; bg=0.314d0  ! coefficients for D-Mass relationship
      !ag=0.19d0 ; bg=0.323d0  ! coefficients for D-Mass relationship
      mu =1.d0/3.d0 ; nu = 1.d0
    endif
    if(isc==6) then ! hail
      av=40.d0 ; bv=0.23d0 ; k=0.5d0
      ag=0.1366d0 ; bg=0.333333d0  ! coefficients for D-Mass relationship
      mu =1.d0/3.d0 ; nu = 1.d0
    endif
    !
    rc = qhydro * rho_d 
    qn = qnhydro * rho_d 
    !--------------------------------------
    ! Mean particle mass, Mean volume diameter
    Mm = qhydro / qnhydro !kg
    Dm = ag * (Mm **bg) !m
    !---------------------------------------
    ! particle density kg/m^3
    if((isc==1) .or. (isc==2)) then
      rho_hydro = 1000.d0
    endif
    !
    if(isc>=3) then
      rho_hydro = Mm * 6.d0 / pi / (Dm ** 3.d0)
    endif
    am =pi*rho_hydro/6.d0 ;   bm=3.d0
    rho = rho_hydro
    !--------------------------------------
    ! Diameter
    do ir=1,nd
      diam(ir)=dmin*1.d-4+dble(ir-1)*(dmax-dmin)/dble(nd-1)*1.d-4 ! diameter in cm
    enddo
    ddiam=diam(10)-diam(9) ! ddiam in cm
    !
    ddiam=ddiam*1.d-2 ! ddiam in m
    diam=diam*1.d-2 ! diameter in m 
    !
    !---------------------------------------
    !
    !---------------------------------------
    ! Fall velocity
    do ir=1,nd
      mass(ir) = (diam(ir) / ag) ** (1.d0/bg)
      fvel(ir) = av * (mass(ir)**bv) * (fc ** k)  ! m/s
    enddo
    do ir=1,nd
      if(ir==1) then
        dmass(ir) = mass(ir+1)-mass(ir)
      else
        dmass(ir) = mass(ir)-mass(ir-1)
      endif
    enddo
    !---------------------------------------
    gam_lam=dgamma((nu+1.d0)/mu)/dgamma((nu+2.d0)/mu);
    gam_n0=1.d0/dgamma((nu+1.d0)/mu);
    !
    if(isc==1) then !cloud
      lambda = (gam_lam * Mm) ** (-mu) !lambda = (0.5d0 * Mm) ** (-mu)
      N0 = mu * qn * gam_n0 * (lambda ** ((nu+1.d0)/mu)) !N0 = mu * qn * (lambda ** 2.d0)
      do ir=1,nd
        NN(ir) = N0* (mass(ir)**nu) * DEXP(-lambda*(mass(ir)**mu)) * dmass(ir)
      enddo
    endif   
    ! 
    !
    if(isc==2) then !rain
      mue = 0.0
      if(Dm .le. cmu3) then
        mue = cmu0*DTANH((4.d0*cmu2*(Dm-cmu3))**cmu5) + cmu4 !Dm [mm]
      else
        mue = cmu1*DTANH((1.d0*cmu2*(Dm-cmu3))**cmu5) + cmu4
      endif
        gamma_alphap1= dgamma(mue+1.d0)
        gfac=dgamma(mue+4.d0) / gamma_alphap1
        ibm=1.d0/3.d0
        lambda = (am*(mue+3.d0)*(mue+2.d0)*(mue+1.d0)/Mm) ** ibm 
        N0 = qn * (lambda ** (mue + 1.d0)) / gamma_alphap1
        do ir=1,nd
          NN(ir) = N0* (diam(ir)**mue) * DEXP(-lambda*diam(ir)) *ddiam
          !if((mass(ir)<2.6d-10) .or. (mass(ir)>5.d-6)) NN(ir)=0 ! commented 20190709
        enddo
    endif
    !
    if(isc>=3) then !ice, snow, graupel, hail
      lambda = (gam_lam * Mm) ** (-mu) 
      N0 = mu * qn * gam_n0 * (lambda ** ((nu+1.d0)/mu)) 
      !lambda = (0.0166667d0 * Mm) ** (-mu) !ice 
      !N0 = mu * qn * 0.5d0 * (lambda ** 3.d0)
      !lambda = (0.1666667d0 * Mm) ** (-mu) !snow
      !N0 = mu * qn * (lambda ** 2.d0) 
      !lambda = (0.297619d-02 * Mm) ** (-mu) !graupel, hail
      !N0 = mu * qn * 0.00833d0 * (lambda ** 6.d0)
      do ir=1,nd
        NN(ir) = N0* (mass(ir)**nu) * DEXP(-lambda*(mass(ir)**mu)) * dmass(ir)
      enddo
    endif
    !
    !
    Deallocate(mass,dmass)
    !
    go to 234
      !if (isc==1) write(*,*) 'cloud',qhydro,qnhydro
      !if (isc==2) write(*,*) 'rain',qhydro,qnhydro
      !if (isc==3) write(*,*) 'ice',qhydro,qnhydro
      !if (isc==4) write(*,*) 'snow',qhydro,qnhydro
      !if (isc==5) write(*,*) 'graup',qhydro,qnhydro
      !if (isc==6) write(*,*) 'ghail',qhydro,qnhydro
      if (isc==1) write(*,*) 'cloud',(pi/6.d0*Sum(rho*NN*diam*diam*diam))/ (qhydro*rho_d)
      if (isc==2) write(*,*) 'rain',(pi/6.d0*Sum(rho*NN*diam*diam*diam))/ (qhydro*rho_d)
      if (isc==3) write(*,*) 'ice',(pi/6.d0*Sum(rho*NN*diam*diam*diam))/ (qhydro*rho_d)
      if (isc==4) write(*,*) 'snow',(pi/6.d0*Sum(rho*NN*diam*diam*diam))/ (qhydro*rho_d)
      if (isc==5) write(*,*) 'graup',(pi/6.d0*Sum(rho*NN*diam*diam*diam))/ (qhydro*rho_d)
      if (isc==6) write(*,*) 'hail',(pi/6.d0*Sum(rho*NN*diam*diam*diam))/ (qhydro*rho_d)
      if (isc==1) write(*,*) 'cloud',Sum(NN)/(qnhydro*rho_d)
      if (isc==2) write(*,*) 'rain',Sum(NN)/(qnhydro*rho_d)
      if (isc==3) write(*,*) 'ice', Sum(NN)/(qnhydro*rho_d)
      if (isc==4) write(*,*) 'snow',Sum(NN)/(qnhydro*rho_d)
      if (isc==5) write(*,*) 'graup',Sum(NN)/(qnhydro*rho_d)
      if (isc==6) write(*,*) 'hail',Sum(NN)/(qnhydro*rho_d)
      !
      !pause
    234 continue
    !
  return
  end subroutine hydro_icon
  !!-- addedd by oue
  !--------------------------------------------------------------------
  !
  !---------------------------------------------------------------------------------------
  !!! Added by oue, 2017/07/21 for RAMS
  subroutine hydro_rams(isc,nd,dmin,dmax,temp,rho_d,qhydro,qnhydro,diam,rho,NN,fvel)
  use crsim_mod
  Use phys_param_mod, ONLY: pi
  Implicit None
  !-- RAMS, 2 moment
  integer, intent(In)     :: isc
  integer,Intent(in)      :: nd
  real*8,Intent(in)       :: dmin,dmax ! um
  real*8,Intent(in)       :: temp      ! C
  real*8,Intent(in)       :: rho_d     ! kg/m^3
  real*8,Intent(in)       :: qhydro    ! kg/kg
  real*8,Intent(in)       :: qnhydro   ! 1/kg 
  real*8,Intent(Out)      :: diam(nd)  ! m
  real*8,Intent(Out)      :: rho(nd)   ! kg/m^3
  real*8,Intent(Out)      :: NN(nd)    ! 1/m^3
  real*8,Intent(Out)      :: fvel(nd)  ! m/s
  !
  real*8                               :: rho_hydro ! kg/m^3
  real*8                               :: av,bv,fc !fall velocity coefficients
  real*8                               :: am,bm,ibm !D-M coefficients
  real*8                               :: nu
  real*8                               :: Dm !mean volime diameter [m]
  integer                              :: ir
  real*8                               :: ddiam,diam_min,diam_max
  !
  real*8                               :: qn ! number concentration,
  real*8                               :: gam_n0
    !
    ! --------------------------------------
    fc=(0.7d0/rho_d) ** 0.362d0
    qn = qnhydro * rho_d 
    !  --------------------------------------
    !
    !
    if(isc==1) then ! cloud
      am=524.d0 ; bm=3.d0  ! coefficients for D-Mass relationship
      nu = 4.d0
      diam_min=2.d-6 ; diam_max=5.d-5 ! m
    endif
    if(isc==7) then ! drizzle
      av=3.173d3 ; bv=2.0d0  ! coefficients for fall velocity
      am=524.d0 ; bm=3.d0  ! coefficients for Dn
      nu = 4.d0
      diam_min=6.5d-5 ; diam_max=1.d-4 ! m
    endif
    !
    if(isc==2) then ! rain
      av=149.d0 ; bv=0.5d0  ! coefficients for fall velocity
      am=524.d0 ; bm=3.d0  ! coefficients for Dn
      nu = 2.d0
      diam_min=1.d-4 ; diam_max=5.d-3 ! m
    endif
    !
    if(isc==3) then ! ice
      av=5.77d5 ; bv=1.88d0  ! coefficients for fall velocity
      am=110.8d0 ; bm=2.91d0  ! coefficients for Dn
      nu = 2.d0
      diam_min=1.5d-5 ; diam_max=1.25d-4 ! m
    endif
     !
    if(isc==4) then ! snow
      av=188.146d0 ; bv=0.933d0  ! coefficients for fall velocity
      am=2.74d-3 ; bm=1.74d0  ! coefficients for Dn
      nu = 2.d0;
      diam_min=1.d-4 ; diam_max=1.d-2 ! m
    endif
    if(isc==8) then ! aggregates
      av=3.084d0 ; bv=0.2d0  ! coefficients for fall velocity
      am=0.496d0 ; bm=2.4d0  ! coefficients for Dn
      nu = 2.d0;
      diam_min=1.d-4 ; diam_max=1.d-2 ! m
    endif
    !
    if(isc==5) then ! graupel
      av=93.3d0 ; bv=0.5d0  ! coefficients for fall velocity
      am=157.d0 ; bm=3.d0  ! coefficients for Dn
      nu = 2.d0;
      diam_min=1.d-4 ; diam_max=5.d-3 ! m
    endif
    if(isc==6) then ! hail
      av=161.d0 ; bv=0.5d0  ! coefficients for fall velocity
      am=471.d0 ; bm=3.d0  ! coefficients for Dn
      nu = 2.d0;
      diam_min=8.d-4 ; diam_max=1.d-2 ! m
    endif
    !
    !--------------------------------------
    ! Mean particle mass, base diameter
    ibm=1.d0/bm
    Dm = (qhydro*dgamma(nu)/qnhydro/am/dgamma(nu+bm)) ** ibm !m
    !---------------------------------------
    ! particle density kg/m^3
    if((isc==1) .or. (isc==2) .or. (isc==7)) then
      rho_hydro = 1000.d0
    endif
    if ((isc==3) .or. (isc==4) .or. (isc==5) .or. (isc==6) .or. (isc==8)) then
      rho_hydro = am * 6.d0 * (Dm ** (bm-3.d0)) / pi
    endif
    rho = rho_hydro
    !--------------------------------------
    ! Diameter
    do ir=1,nd
      diam(ir)=dmin*1.d-4+dble(ir-1)*(dmax-dmin)/dble(nd-1)*1.d-4 ! diameter in cm
    enddo
    ddiam=diam(10)-diam(9) ! ddiam in cm
    !
    ddiam=ddiam*1.d-2 ! ddiam in m
    diam=diam*1.d-2 ! diameter in m 
    !---------------------------------------
    ! Fall velocity
    fvel=0.0d0 ! cloud
    if(isc>1) then
      do ir=1,nd
        fvel(ir) = av * (diam(ir)**bv) * fc   ! m/s
      enddo
    endif
    !---------------------------------------
    gam_n0=1.d0/dgamma(nu)
    !
    ! size distribution
    do ir=1,nd
      NN(ir) = qn*gam_n0 * (diam(ir)**(nu-1))/(Dm**nu) * DEXP(-(diam(ir)/Dm)) * ddiam
      if(diam(ir) < diam_min) NN(ir)=0.d0
      if(diam(ir) > diam_max) NN(ir)=0.d0
    enddo
    !
    go to 234
      if (isc==1) write(*,*) 'cloud',(pi/6.d0*Sum(rho*NN*diam*diam*diam))/ (qhydro*rho_d)
      if (isc==7) write(*,*) 'drizzle',(pi/6.d0*Sum(rho*NN*diam*diam*diam))/ (qhydro*rho_d)
      if (isc==2) write(*,*) 'rain',(pi/6.d0*Sum(rho*NN*diam*diam*diam))/ (qhydro*rho_d)
      if (isc==3) write(*,*) 'ice',(pi/6.d0*Sum(rho*NN*diam*diam*diam))/ (qhydro*rho_d)
      if (isc==4) write(*,*) 'snow',(pi/6.d0*Sum(rho*NN*diam*diam*diam))/ (qhydro*rho_d)
      if (isc==8) write(*,*) 'aggrregates',(pi/6.d0*Sum(rho*NN*diam*diam*diam))/ (qhydro*rho_d)
      if (isc==5) write(*,*) 'graup',(pi/6.d0*Sum(rho*NN*diam*diam*diam))/ (qhydro*rho_d)
      if (isc==6) write(*,*) 'hail',(pi/6.d0*Sum(rho*NN*diam*diam*diam))/ (qhydro*rho_d)
      if (isc==1) write(*,*) 'cloud',Sum(NN)/(qnhydro*rho_d)
      if (isc==7) write(*,*) 'drizzle',Sum(NN)/(qnhydro*rho_d)
      if (isc==2) write(*,*) 'rain',Sum(NN)/(qnhydro*rho_d)
      if (isc==3) write(*,*) 'ice', Sum(NN)/(qnhydro*rho_d)
      if (isc==4) write(*,*) 'snow',Sum(NN)/(qnhydro*rho_d)
      if (isc==8) write(*,*) 'aggregates',Sum(NN)/(qnhydro*rho_d)
      if (isc==5) write(*,*) 'graup',Sum(NN)/(qnhydro*rho_d)
      if (isc==6) write(*,*) 'hail',Sum(NN)/(qnhydro*rho_d)
      !
      !pause
    234 continue
    !
    !
  return
  end subroutine hydro_rams
  !!-- addedd by oue
  !--------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------------
  ! subroutine cloud_P3 is added by DW 2017/10/30 for P3
  !
  subroutine cloud_P3(nd,dmin,dmax,temp,rho_d,qhydro,qnhydro,diam,rho,NN,fvel)
  Use crsim_mod
  Use phys_param_mod, ONLY: pi,grav,oneOverThree, rhow, rhow_25C
  Implicit None
  !
  integer,Intent(in)      :: nd
  real*8,Intent(in)       :: dmin,dmax ! um
  real*8,Intent(in)       :: temp      ! C
  real*8,Intent(inout)    :: qhydro    ! kg/kg
  real*8,Intent(in)       :: qnhydro   ! 1/kg
  real*8,Intent(in)       :: rho_d     ! kg/m3
  real*8,Intent(out)      :: diam(nd)  ! m
  real*8,Intent(out)      :: rho(nd)   ! kg/m3
  real*8,Intent(Out)      :: NN(nd)    ! 1/m3
  real*8,Intent(Out)      :: fvel(nd)  ! m/s
  !
  real*8                           :: visc           ! Viscosity kg/m/s
  real*8,parameter                 :: Nconst=200.d0  ! 1/cm3 ! change 400 to 200, consistent with WRFv4.1, 20190709
  real*8                           :: av,bv
  real*8                           :: mu
  real*8,Dimension(:),Allocatable  :: lambda,N0
  real*8                           :: Nc
  !
  integer                            :: ir
  real*8                             :: ddiam  ! m
    !
    do ir=1,nd
      diam(ir)=dmin*1.d-4+dble(ir-1)*(dmax-dmin)/dble(nd-1)*1.d-4 ! diameter in cm
    enddo
    ddiam=diam(10)-diam(9) ! dr in cm
    ddiam=ddiam*1.d-2 ! dr in m
    diam=diam*1.d-2 ! diameter in m
    !---------------------------------------
    !
    visc=1.496d0*1.d-6*(temp+273.15d0)**1.5d0/(temp+393.15d0)
    !
    ! parametrs for fall velocity 
    av=grav*rhow/(18.d0*visc)  !'a' parameter for droplet fallspeed (Stokes'law)
    bv=2.d0  
    !-----------------------------
    if (qnhydro>0.d0) then
      Nc=qnhydro
    else
      Nc=(Nconst*1.d+6)/rho_d ! in 1/kg
    endif
    !-----------------------------
    mu=0.0005714d0*(Nc*1.d-6*rho_d)+0.2714d0
    mu=1.d0/(mu*mu)-1.d0
    ! 
    mu=MAX(mu,2.d0)
    mu=MIN(mu,15.d0)
    !
    Allocate(lambda(nd),N0(nd))
    !
    lambda=(pi*rhow_25C*Nc)/(6.d0*qhydro)*(mu+3.d0)*(mu+2.d0)*(mu+1.d0)
    lambda=lambda**oneOverThree
    !
    N0=(Nc*lambda**(mu+1.d0))/(dgamma(mu+1.d0))
    !
    do ir=1,nd
      fvel(ir) = av*(diam(ir))**bv  ! m/s 
      NN(ir)=(N0(ir)*(diam(ir))**mu)*dexp(-lambda(ir)*diam(ir))
      NN(ir)=NN(ir)*ddiam   ! in m-3
    enddo
    !
    Deallocate(lambda,N0)
  return
  end subroutine cloud_P3
  !---------------------------------------
  ! subroutine rain_P3 added by DW 2017/10/30 for P3
  subroutine rain_P3(nd,dmin,dmax,temp,rho_d,qhydro,qnhydro,diam,rho,NN,fvel)
  Use crsim_mod
  Use phys_param_mod, ONLY: pi,grav, rhow,rhow_25C, Rd,p0, T0K, oneOverThree
  Implicit None
  !
  integer,Intent(in)      :: nd
  real*8,Intent(in)       :: dmin,dmax ! um
  real*8,Intent(in)       :: temp      ! C
  real*8,Intent(inout)    :: qhydro    ! kg/kg
  real*8,Intent(in)       :: qnhydro   ! 1/kg  
  real*8,Intent(in)       :: rho_d 
  real*8,Intent(Out)      :: diam(nd)   ! m
  real*8,Intent(Out)      :: rho(nd)   ! kg/m^3
  real*8,Intent(Out)      :: NN(nd)    ! 1/m^3
  real*8,Intent(Out)      :: fvel(nd)  ! m/s
  !
  real*8                               :: rho_surf ! kg/m^3
  real*8                               :: k,fc,fck
  real*8,Dimension(:),Allocatable      :: lambda,N0,amg
  real*8,Dimension(150)                :: mu_r_table
  !real*8                               :: rdumii
  real*8                               :: dum, inv_dum,lamold
  real*8                               :: mu
  real*8                               :: initlamr
  !
  integer                              :: ir,i,ii !,dumii
  real*8                               :: ddiam
  real*8                               :: lambdae
  real*8,parameter                     :: mu_r_constant=0.d0 ! apply constant mu, consistent with WRFv4.1, 20190709
  !
  !
    do ir=1,nd
      diam(ir)=dmin*1.d-4+dble(ir-1)*(dmax-dmin)/dble(nd-1)*1.d-4 ! diameter in cm
    enddo
    ddiam=diam(10)-diam(9) ! ddiam in cm
    ddiam=ddiam*1.d-2 ! ddiam in m
    diam=diam*1.d-2 ! diameter in m 
    !
    k=0.54d0
    rho_surf=p0/(Rd*T0K)
    fc=rho_surf/rho_d
    !
    ! Generate lookup table for rain shape parameter mu (150x1)
    !
    do i = 1,150   ! loop over lookup table values
      initlamr = 1.d0/((dble(i)*2.d0)*1.d-6+250.d-6)
      ! mu-lambda relationship is from Cao et al.(2008), eq.(7)
      ! first guess
      mu = 0.d0
      do ii=1,50
        lambdae = initlamr*(dgamma(mu+4.d0)/(6.d0*dgamma(mu+1.d0)))**oneOverThree
        ! new estimate for mu based on lambda
        ! set max lambda in formula for mu to 20 mm-1, according to Cao et al.
        dum = min(20.d0,lambdae*1.d-3)
        
        !== apply constant mu, consistent with WRFv4.1, 20190709
        mu = mu_r_constant 
        !!mu  = max(0.d0,-0.0201d0*dum**2.d0+0.902d0*dum-1.718d0) ! commented out to fix mu=0, consistent with WRFv4.1, 20190709 
    
        ! if lambda is converged within 0.1%, then exit loop
        if (ii.ge.2) then
          if (dabs((lamold-lambdae)/lambdae).lt.0.001d0) goto 111
        endif
        lamold = lambdae
      enddo

      111 continue
      mu_r_table(i) = mu
    enddo
    !
    inv_dum = (qhydro/(pi*rhow*qnhydro))**oneOverThree

!==== apply constant mu, consistent with WRFv4.1, 20190709
    mu = mu_r_constant 
!    !
!    if (inv_dum.lt.282.d-6) then
!      mu = 8.282d0
!    elseif (inv_dum.ge.282.d-6 .and. inv_dum.lt.502.d-6) then
!      ! interpolate
!      rdumii  = (inv_dum-250.d-6)*1.d+6*0.5d0
!      rdumii  = max(rdumii,1.d0)
!      rdumii  = min(rdumii,150.d0)
!      dumii   = int(rdumii)
!      dumii   = min(149,dumii)
!      mu      = mu_r_table(dumii)+(mu_r_table(dumii+1)-mu_r_table(dumii))*(rdumii-real(dumii))
!    elseif (inv_dum.ge.502.d-6) then
!      mu  = 0.d0
!    endif
!====
    !
    Allocate(lambda(nd),N0(nd),amg(nd))
    !
    lambda = pi/6.d0*rhow_25C*qnhydro*(mu+3.d0)*(mu+2.d0)*(mu+1.d0)/(qhydro)
    lambda = lambda**oneOverThree
    !
    lambda = min((mu+1.d0)*1.d+5,lambda)
    lambda = max((mu+1.d0)*1250.d0,lambda)
    !
    N0=(qnhydro*lambda**(mu+1.d0))/(dgamma(mu+1.d0))
    !
    do ir=1,nd
      amg(ir) = pi/6.d0*rhow*diam(ir)**3.d0      ! mass in kg
      amg(ir) = amg(ir)*1000.d0                  ! kg to g
      if (diam(ir)*1.d+6.le.134.43d0)      then
        fvel(ir) = 4.5795d+3*amg(ir)**(2.d0/3.d0)
      elseif (diam(ir)*1.d+6.lt.1511.64d0) then
        fvel(ir) = 4.962d+1*amg(ir)**(1.d0/3.d0)
      elseif (diam(ir)*1.d+6.lt.3477.84d0) then
        fvel(ir) = 1.732d+1*amg(ir)**(1.d0/6.d0)
      else
        fvel(ir) = 9.17d0
      endif
      !
      fck=fc**k
      fvel(ir)=fvel(ir)*fck  ! air density modification of particle fall speeds
      !
      NN(ir)=(N0(ir)*(diam(ir))**mu)*dexp(-lambda(ir)*diam(ir))
      NN(ir)=NN(ir)*ddiam   ! in m-3
    enddo
    !
    Deallocate(lambda,N0,amg)
    !
  return
  end subroutine rain_P3
  !-------------------------------
  ! subroutine ice_P3 added by DW 2017/10/30 for P3
  subroutine ice_P3(isc,nd,dmin,dmax,temp,rho_d,qhydro,qnhydro,qir,fr,rrho,diam,rho,NN,fvel,lookup_file_1)
  Use crsim_mod
  use phys_param_mod, ONLY: pi,Rd,T0K,p0_600
  Implicit None
  !
  INTERFACE 
    DOUBLE PRECISION FUNCTION gammq(a,x)
      Real*8,Intent(in) :: a,x
      Real*8            :: gammcf,gamser,gln
    END FUNCTION gammq
  END INTERFACE
  !
  Integer,Intent(in)        :: nd,isc
  character(*), intent(in)  :: lookup_file_1
  !
  integer,parameter       :: isize   = 50
  integer,parameter       :: densize = 10
  integer,parameter       :: rimsize = 4
  integer,parameter       :: tabsize = 8  ! number of quantities used from lookup table
  real*8,dimension(densize,rimsize,isize,tabsize) :: itab  ! ice lookup table values
  !
  real*8,Intent(in)       :: dmin,dmax ! um
  real*8,Intent(in)       :: temp      ! C
  real*8,Intent(inout)    :: qhydro    ! kg/kg
  real*8,Intent(inout)    :: qnhydro   ! 1/kg   
  real*8,Intent(in)       :: fr        ! rime fraction 
  real*8,Intent(in)       :: qir       ! kg/kg rime mixing ratio
  real*8,Intent(in)       :: rrho      ! kg/m3 bulk density
  real*8,Intent(in)       :: rho_d     ! kg/m3 air density
  !
  real*8,Intent(Out)      :: diam(nd)      ! m  
  real*8,Intent(Out)      :: rho(nd)       ! kg/m^3
  real*8,Intent(Out)      :: NN(nd)        ! 1/m^3
  real*8,Intent(Out)      :: fvel(nd)      ! m/s
  !
  real*8                  :: rho_surf, a_rho, b_rho ! kg/m^3
  real*8                  :: k,fc
  real*8                  :: lambda,N0
  real*8                  :: mu,mu_a
  !
  integer                 :: ir,i,ii,jj
  integer                 :: duma,dumb,dumc
  real*8                  :: ddiam
  real*8                  :: f1pr02,f1pr03,f1pr04,f1pr05,&
                                           f1pr06,f1pr07,f1pr08,fvel1
  !
  real*8                  :: ds,cs,cgp,dg,csr,dsr,cs1,ds1,bas,aas,bag,aag,bas1,aas1
  real*8                  :: dum,dumm,dum1,dum2,dum3,dum4,dum11,dum22,dum33
  real*8                  :: dcrit,dcrits,dcritr
  real*8                  :: m1,m2,m3
  real*8                  :: a0,b0,c0,c1,c2,del0
  real*8                  :: sum1
    !
    k=0.54d0
    !
    do ir=1,nd
      diam(ir)=dmin*1.d-4+dble(ir-1)*(dmax-dmin)/dble(nd-1)*1.d-4 ! diameter in cm
    enddo
    ddiam=diam(10)-diam(9) ! ddiam in cm
    !
    ddiam=ddiam*1.d-2 ! ddiam in m
    diam=diam*1.d-2   ! diameter in m 
    !------------------------------------------
    ! from dmin to dcrit    (small ice)
    ! from dcrit to dcrits  (unrimed ice)
    ! from dcrits to dcritr (graupel)
    ! from dcritr to dmax   (partially rimed ice)
    !
    rho_surf=p0_600/(Rd*T0K)              ! air density at the reference level(kg m-3)
    mu_a=1.496d-6*253.15d0**1.5d0/(253.15d0+120.d0)/rho_surf   ! viscosity of air
    fc=rho_surf/rho_d
    !
    ! m-D parameters for unrimed ice:
    ! Brown and Francis (1995)
    ds = 1.9d0
    cs = 0.0121d0
    ! m-D parameter for fully rimed ice (spherical)
    dg  = 3.d0
    bag = 2.d0
    aag = pi*0.25d0
    ! aggreagtes of side planes, bullets, etc.
    bas = 1.88d0
    aas = 0.2285d0*100.d0**bas/(100.d0**2.d0)
    ! parameters for surface roughness of ice particle used for fallspeed (mitchell and heymsfield 2005)
    del0 = 5.83d0
    c0   = 0.6d0
    c1   = 4.d0/(del0**2.d0*c0**0.5d0)
    c2   = del0**2.d0/4.d0
    !
    dcrit = (pi/(6.d0*cs)*900.d0)**(1.d0/(ds-3.d0))
    !
    ! set up m-D relationship for solid ice with D < Dcrit
    cs1  = pi*900.d0/6.d0
    ds1  = 3.d0
    !
    ! read in ice PSD parameter lookuptable
    !
    open(unit=10,file=lookup_file_1, status='old')
    !
    do jj = 1,densize
      do ii = 1,rimsize
        do i = 1,isize
          read(10,*) dum,dum,dum,itab(jj,ii,i,1),itab(jj,ii,i,2),itab(jj,ii,i,3),&
                     itab(jj,ii,i,4),itab(jj,ii,i,5),itab(jj,ii,i,6),itab(jj,ii,i,7),&
                     itab(jj,ii,i,8)
        enddo
      enddo
    enddo
    close(10)
    !
    call find_lookupTable_indices_1(duma,dumb,dumc,dum11,dum22,dum33,&
                                    isize,rimsize,densize,qhydro,qnhydro,qir,rrho)
    !
    call access_lookup_table(dumb,dumc,duma,2,densize,rimsize,isize,tabsize,itab,f1pr02)
    call access_lookup_table(dumb,dumc,duma,3,densize,rimsize,isize,tabsize,itab,f1pr03)
    call access_lookup_table(dumb,dumc,duma,4,densize,rimsize,isize,tabsize,itab,f1pr04)
    call access_lookup_table(dumb,dumc,duma,5,densize,rimsize,isize,tabsize,itab,f1pr05)
    call access_lookup_table(dumb,dumc,duma,6,densize,rimsize,isize,tabsize,itab,f1pr06)
    call access_lookup_table(dumb,dumc,duma,7,densize,rimsize,isize,tabsize,itab,f1pr07)
    call access_lookup_table(dumb,dumc,duma,8,densize,rimsize,isize,tabsize,itab,f1pr08)
    !
    lambda    = f1pr02
    mu        = f1pr03
    dcrits    = f1pr04
    dcritr    = f1pr05
    csr       = f1pr06
    dsr       = f1pr07
    cgp       = f1pr08
    !
    ! adapted N0 to real Q 
    dum1 = lambda**(-ds1-mu-1.d0)*dgamma(mu+ds1+1.d0)*(1.d0-gammq(mu+ds1+1.d0,dcrit*lambda))
    dum2 = lambda**(-ds-mu-1.d0)*dgamma(mu+ds+1.d0)*(gammq(mu+ds+1.d0,dcrit*lambda))
    dumm = lambda**(-ds-mu-1.d0)*dgamma(mu+ds+1.d0)*(gammq(mu+ds+1.d0,dcrits*lambda))
    dum2 = dum2-dumm
    dum3 = lambda**(-dg-mu-1.d0)*dgamma(mu+dg+1.d0)*(gammq(mu+dg+1.d0,dcrits*lambda))
    dumm = lambda**(-dg-mu-1.d0)*dgamma(mu+dg+1.d0)*(gammq(mu+dg+1.d0,dcritr*lambda))
    dum3 = dum3-dumm
    dum4 = lambda**(-dsr-mu-1.d0)*dgamma(mu+dsr+1.d0)*(gammq(mu+dsr+1.d0,dcritr*lambda))
    !
    N0 = qhydro/(cs1*dum1+cs*dum2+cgp*dum3+csr*dum4)
    !
    if(isc==3) then  !-small ice
      !
      sum1 = 0.d0
      !
      do ir=1,nd
        !
        if (diam(ir).le.dcrit) then
          NN(ir)=N0*diam(ir)**mu*dexp(-lambda*diam(ir))
          NN(ir)=NN(ir)*ddiam   ! in m-3 
          if (NN(ir).le.1) then 
            NN(ir)=0.d0
          endif
        else 
          NN(ir)=0.d0
        endif
        !
        !-recalculate Q for small ice
        sum1 = sum1+NN(ir)*cs1*diam(ir)**ds1
      enddo ! loop ir
      qhydro = sum1/rho_d
      !
    else if(isc==4) then !-unrimed ice
      !
      sum1 = 0.d0
      !
      do ir=1,nd
        if (diam(ir).gt.dcrit.and.diam(ir).le.dcrits) then
          NN(ir)=N0*diam(ir)**mu*dexp(-lambda*diam(ir))
          NN(ir)=NN(ir)*ddiam
          if (NN(ir).le.1) then 
            NN(ir)=0.d0
          endif
        else
          NN(ir)=0.d0
        endif   
        !
        !-recalculate Q for unrimed ice
        sum1 = sum1+NN(ir)*cs*diam(ir)**ds
        !
      enddo ! loop ir  
      qhydro = sum1/rho_d
      ! 
    else if(isc==5) then ! graupel
      !
      sum1 = 0.d0
      !
      do ir=1,nd
        if (diam(ir).gt.dcrits.and.diam(ir).le.dcritr) then
          NN(ir)=N0*diam(ir)**mu*dexp(-lambda*diam(ir))
          NN(ir)=NN(ir)*ddiam
          if (NN(ir).le.1) then 
            NN(ir)=0.d0
          endif
        else
          NN(ir)=0.d0
        endif  
        !
        !-recalculate Q 
        sum1 = sum1+NN(ir)*cgp*diam(ir)**dg  
        !
      enddo ! loop ir   
      qhydro = sum1/rho_d
      !
    else if(isc==6) then ! partically
      !
      sum1 = 0.d0
      !
      do ir=1,nd
        if (diam(ir).gt.dcritr) then
          NN(ir)=N0*diam(ir)**mu*dexp(-lambda*diam(ir))
          NN(ir)=NN(ir)*ddiam
        else
          NN(ir)=0.d0
        endif
        !-recalculate Q 
        sum1 = sum1+NN(ir)*csr*diam(ir)**dsr
        !
      enddo ! loop ir 
      qhydro = sum1/rho_d
      !
    endif
   
    ! fall speed
    ! neglect turbulent correction for aggregates for now correction for turbulence
    a0 = 0.d0
    b0 = 0.d0
    !
    if(isc==3) then
      do ir=1,nd
        if (diam(ir).le.dcrit) then
          cs1  = pi*900.d0/6.d0
          ds1  = 3.d0
          bas1 = 2.d0
          aas1 = pi/4.d0
          rho(ir) = 917.d0 
          call fvel_P3(cs1,ds1,bas1,aas1,c1,c2,a0,b0,mu_a,rho_surf,diam(ir),fvel1)
          fvel(ir) = fvel1 
        else 
          fvel(ir)=0.d0
          rho(ir) =0.d0
        endif
      enddo        
    endif 
    !
    if(isc==4) then
      do ir=1,nd
        if (diam(ir).gt.dcrit.and.diam(ir).le.dcrits) then
          cs1  = cs
          ds1  = ds
          bas1 = bas
          aas1 = aas
          a_rho = cs1*6.d0/pi 
          b_rho = ds1-3.d0
          rho(ir) = a_rho*diam(ir)**b_rho
          call fvel_P3(cs1,ds1,bas1,aas1,c1,c2,a0,b0,mu_a,rho_surf,diam(ir),fvel1)
          fvel(ir) = fvel1  
        else
          fvel(ir)=0.d0
          rho(ir) =0.d0
        endif
      enddo
    endif
    !
    if(isc==5) then
      do ir=1,nd
        if (diam(ir).gt.dcrits.and.diam(ir).le.dcritr) then         
          cs1  = cgp
          ds1  = dg
          bas1 = bag
          aas1 = aag
          a_rho = cs1*6.d0/pi
          b_rho = ds1-3.d0
          rho(ir) = a_rho*diam(ir)**b_rho
          call fvel_P3(cs1,ds1,bas1,aas1,c1,c2,a0,b0,mu_a,rho_surf,diam(ir),fvel1)        
          fvel(ir) = fvel1     
        else
          fvel(ir)=0.d0
          rho(ir) =0.d0
        endif  
      enddo  
    endif
    !
    if(isc==6) then
      do ir=1,nd
        if (diam(ir).gt.dcritr) then
          cs1   = csr
          ds1   = dsr
          a_rho = cs1*6.d0/pi
          b_rho = ds1-3.d0
          rho(ir) = a_rho*diam(ir)**b_rho
          if (fr.eq.1.d0) then
            aas1 = aas
            bas1 = bas
          else ! for projected area, keep bas1 constant, but modify aas1 according to rimed fraction
            bas1 = bas
            dum1 = aas*diam(ir)**bas
            dum2 = aag*diam(ir)**bag
            m1   = cs1*diam(ir)**ds1
            m2   = cs*diam(ir)**ds
            m3   = cgp**diam(ir)**dg
            ! linearly interpolate based on particle mass
            dum3 = dum1+(m1-m2)*(dum2-dum1)/(m3-m2)
            aas1 = dum3/(diam(ir)**bas)
          endif
          call fvel_P3(cs1,ds1,bas1,aas1,c1,c2,a0,b0,mu_a,rho_surf,diam(ir),fvel1)
          fvel(ir) = fvel1
        else
          fvel(ir)=0.d0
          rho(ir) =0.d0
        endif
      enddo   
    endif
    !
  return
  end subroutine ice_P3
  !-------------------------------------------
  ! subroutine fvel_P3 is added by DW 2017/10/30 for P3
  subroutine fvel_P3(cs1,ds1,bas1,aas1,c1,c2,a0,b0,mu_a,rho_surf,diam,fvel)
  Use phys_param_mod, ONLY: grav
  
  ! fall speed !
  Implicit None
  Real*8,Intent(in)      :: cs1,ds1,bas1,aas1,c1,c2,a0,b0,mu_a,rho_surf,diam
  Real*8,Intent(out)     :: fvel
  Real*8                 :: xx,b1,a1
  ! 
    xx = 2.d0*cs1*grav*rho_surf*diam**(ds1+2.d0-bas1)/(aas1*(mu_a*rho_surf)**2.d0)
    ! drag terms
    b1 = c1*xx**0.5d0/(2.d0*((1.d0+c1*xx**0.5d0)**0.5d0-1.d0)*&
         (1.d0+c1*xx**0.5d0)**0.5d0)-a0*b0*xx**&
         b0/(c2*((1.d0+c1*xx**0.5d0)**0.5d0-1.d0)**2.d0)
    a1 = (c2*((1.d0+c1*xx**0.5d0)**0.5d0-1.d0)**2.d0-a0*xx**b0)/xx**b1
    ! velocity in terms of drag terms
    fvel = a1*mu_a**(1.d0-2.d0*b1)*(2.d0*cs1*grav/(rho_surf*aas1))&
           **b1*diam**(b1*(ds1-bas1+2.d0)-1.d0)
     
  return
  end subroutine fvel_P3
  !
  !-----------------------------------------------------------------
  subroutine hydro_jfan(conf,nkr,nsc,scatt_type)
  Use crsim_mod
  Use phys_param_mod, ONLY : m999
  Implicit None
  Type(conf_var),Intent(in)           :: conf
  Integer,Intent(In)                  :: nkr ! total number of bins in each category
  Integer,Intent(In)                  :: nsc ! total number of scattering types
  Type(scatt_type_var),Intent(InOut)  :: scatt_type(nsc)
  !
  Integer, parameter  :: id=1
  Integer, parameter  :: nf=4  ! number of jfan input files with microph. info  
  Integer             :: nht=3 ! number of hydro category 1=cloud+rain,2=ice+snow,3=graupel+hail  
  !
  Integer               :: iht      ! the hydro category 1=cloud+rain, 2=ice+snow,3=graupel+hail
  Integer               :: ir1, ir2 ! the first and the last bin number of the category
  Integer               :: isc      ! counter for scattering types (1=cloud,2=rain,3=ice, 4=snow,5=graupel)
  Integer               :: ifile
  !
  Real*8,Allocatable,Dimension(:)      :: RO1BL,RO3BL,RO4BL,RO5BL
  Real*8,Allocatable,Dimension(:,:)    :: RO2BL
  Real*8,Allocatable,Dimension(:,:,:)  :: work
  !
  !
  !WRF_MP_PHYSICS_20_InputFiles
  !InpFileName(1)='bulkradii.asc_s_0_03_0_9'
  !InpFileName(2)='bulkdens.asc_s_0_03_0_9'
  !InpFileName(3)='masses.asc'
  !InpFileName(4)='termvels.asc'
  !   
    nht=scatt_type(nsc)%ihtf ! or max scatt_type(isc)%ihtf  , =3
    ! read additional input jfan files with masses,densities,fall velocities and diams
    Allocate(work(nht,nkr,nf))
    !
    Allocate( RO1BL(nkr),RO3BL(nkr),RO4BL(nkr),RO5BL(nkr))
    Allocate( RO2BL(nkr,nht))
    !
    RO1BL=m999 ; RO3BL=m999; RO4BL=m999; RO5BL=m999
    RO2BL=m999
    work=m999
    !
    do ifile=1,nf
      Open(unit=id,File=Trim(conf%WRF_MP_PHYSICS_20_InputFiles(ifile)))
      read(id,*) RO1BL,RO2BL,RO3BL,RO4BL,RO5BL
      work(1,:,ifile)=RO1BL(:)
      work(2,:,ifile)=RO3BL(:)
      work(3,:,ifile)=RO4BL(:)
      close(id)
    enddo ! ifile
    !
    Deallocate(RO1BL,RO3BL,RO4BL,RO5BL,RO2BL)
    !
    !!
    !-------------------------------------------------------------------------------
    !! INFO FOR FILES 
    !! work(iht,ir,ifile)
    !! iht=1 liquid, iht=2 ice and snow, iht=3 graupel and hail
    !! ir =1 -> 33 number of radius bins
    !! ifile=1  bin radius (cm), =2 density(gr/cm3), =3 mass (gr), =4 term vel in cm/s
    !! ------------------------------
    !
    ! init
    ir1=-1
    ir2=-1
    iht=-1
    ! 
    do isc=1,nsc
      !
      if (isc==1) then 
        iht=scatt_type(isc)%ihtf 
        ir1=scatt_type(isc)%id1
        ir2=scatt_type(isc)%id2
      endif
      !
      if (isc==2) then
        iht=scatt_type(isc)%ihtf
        ir1=scatt_type(isc)%id1
        ir2=scatt_type(isc)%id2
      endif
      !
      if (isc==3) then
        iht=scatt_type(isc)%ihtf
        ir1=scatt_type(isc)%id1
        ir2=scatt_type(isc)%id2
      endif
      !
      if (isc==4) then
        iht=scatt_type(isc)%ihtf
        ir1=scatt_type(isc)%id1
        ir2=scatt_type(isc)%id2
      endif
      !
      if (isc==5) then
        iht=scatt_type(isc)%ihtf
        ir1=scatt_type(isc)%id1
        ir2=scatt_type(isc)%id2
      endif
      !
      if ((ir1 == -1) .or. (ir2 == -1) .or. (iht == -1)) then
      write(*,*) 'Problem in subroutine hydro_jfan. Stopping now.'
      call Exit(1)
      endif
      !
      scatt_type(isc)%diam=2.d0*work(iht,ir1:ir2,1)*1.d-2 ! cm -> m
      scatt_type(isc)%rho=work(iht,ir1:ir2,2)*1.d+3  ! gr/cm^3 -> kg/m^3
      scatt_type(isc)%mass=work(iht,ir1:ir2,3)*1.d-3 ! gr -> kg
      scatt_type(isc)%fvel=work(iht,ir1:ir2,4)*1.d-2 ! cm/s ->m/s
      !
    enddo ! isc
    !
    Deallocate(work)
    !
  return
  end subroutine hydro_jfan
  !
  !-----------------------------------------------------------------
   
  subroutine hydro_samsbm(conf,nsc,nkr,scatt_type)
  Use crsim_mod
  Use phys_param_mod, ONLY: m999
  Implicit None
  Type(conf_var),Intent(in)           :: conf
  Integer,Intent(In)                  :: nsc ! total number of scattering types
  Integer,Intent(InOut)               :: nkr ! total number of bins in each category
  Type(scatt_type_var),Intent(InOut)  :: scatt_type(nsc)
  !
  Integer, parameter  :: id=1
  !
  Integer               :: iht      ! the hydro category 1=cloud+rain, 
  Integer               :: ir1, ir2 ! the first and the last bin number of the category
  Integer               :: isc      ! counter for scattering types (1=cloud,2=rain)
  !
  Real*8,Allocatable,Dimension(:)      :: RO1BL,RO3BL,RO4BL,RO2BL
  character(len=220)             :: inputline='#'!
  Integer                 :: n, nnkr, ikr ! number of recorded in file, number of bins in file
  Integer                 :: i
    !
    !SAM warm bin information
    !conf%WRFmpInputFile='micro_taubin_info.txt'
    
    !nht=scatt_type(nsc)%ihtf ! or max scatt_type(isc)%ihtf  , =1 (cloud and rain)
    
    ! read additional input SAM files with masses, diameter, 
    !
    !
    Open(unit=id,File=Trim(conf%WRFmpInputFile))
    !- check number of bins
    n = 0
    read (id, '()') !skip header
    read (id, '()') !skip header
    do
      !read (id, *, end=100) dummy1, dummy2, z  ! go to 100 if file ends
      read(id,'(a)',end=100) inputline ! go to 100 if file ends
      n = n + 1
    end do
    100 continue
    rewind (id)  ! go back to top of file
    nnkr = (n - 2)/2
    if (nkr /= nnkr) then
      write(*,*) 'Worning: the number of bins in txt file and nkr set are inconsitent.',nnkr,nkr
    endif
    Allocate( RO1BL(nkr+1),RO2BL(nkr+1),RO3BL(nkr+1),RO4BL(nkr+1))
    RO1BL=m999 ; RO3BL=m999; RO4BL=m999; RO2BL=m999
    ! ===read again ===
    write(*,*) 'Index= ; ', 'Mass(g/kg)= ; ', 'Diameter(cm)= '
    read (id, '()') !skip header
    read (id, '()') !skip header
    ikr=1
    A1: do i = 1, n
      read(id,'(a)') inputline
      if (inputline(2:6).eq.'-----') then
        read(inputline(8:46),*)  RO1BL(ikr),RO2BL(ikr),RO3BL(ikr),RO4BL(ikr)
        write(*,*) ikr,RO2BL(ikr),RO3BL(ikr)
        ikr=ikr+1
      endif
      if (ikr > (nkr+1)) then
        exit A1
      endif
    enddo A1
     
    close(id)
    
    ! get diameter each bin
    !-------------------------------------------------------------------------------
    !! INFO FOR FILES 
    !! iht=1 cloud +  rain 
    !! ir =1 -> 33 number of radius bins
    !! RO1BL= bin index, RO2BL=mass(g), RO3BL=diameter(cm),RO4BL=radius(cm)
    !! ------------------------------
    
    do isc=1,nsc
      ! isc=1: cloud, isc=2: rain
      iht=scatt_type(isc)%ihtf 
      ir1=scatt_type(isc)%id1
      ir2=scatt_type(isc)%id2
      if (isc==1) write(*,*) 'cloud indices:',ir1,ir2
      if (isc==2) write(*,*) 'rain indices:',ir1,ir2
      scatt_type(isc)%diam1=RO3BL(ir1:ir2)*1.d-2 ! cm -> m
      scatt_type(isc)%mass1=RO2BL(ir1:ir2)*1.d-3 ! gr -> kg
      !
    enddo ! isc
    !
    Deallocate(RO1BL,RO3BL,RO4BL,RO2BL)
  return
  end subroutine hydro_samsbm
     !
  !--------------------------
  subroutine hydro_samsbm_2mom(isc,nd,rhod,diam1,mass1,N,M,diam,rho,fvel)
  !-- Computation of diameter and fall velocity from mass, because SAM bin is 2 moments (mass and number)
  Use phys_param_mod, ONLY: pi
  Implicit None
  Integer,Intent(In)                 :: isc     ! scattering type 1: cloud, 2: rain
  integer,Intent(In)                :: nd      ! number of diameters
  Real*8,Intent(In)                ::  rhod    ! air density kg/m^3
  Real*8,Intent(In)                ::  diam1(nd+1)    ! boundary diameters in m (nd+1)
  Real*8,Intent(In)                ::  mass1(nd+1)    ! boundary mass in kg (nd+1)
  Real*8,Intent(In)                ::  N(nd)    ! number 1/m3 (nd)
  Real*8,Intent(In)                ::  M(nd)    ! mixing ratio kg/kg (nd)
  Real*8,Intent(Out)                ::  diam(nd)    !  diameters in m nd+1
  Real*8,Intent(Out)                ::  rho(nd)    !  particle bulk density kg/m3 (nd)
  Real*8,Intent(Out)                ::  fvel(nd)    !  particle fall velocity m/s (nd)
  Real*8                            :: am,alp,bet
  Integer  :: i
    !
    ! Initialize rho
    rho = 1.d0*1.d+3  ! gr/cm^3  -> kg/m^3
    !
    if ((isc==1).or.(isc==2)) then ! cloud
      rho = 1.d0*1.d+3  ! gr/cm^3  -> kg/m^3
      am = 0.d0
      bet = 0.d0
      alp = 0.d0
      do i=1,nd
        !diam(i) = ( M(i) / ( N(i) + eps ) * rhod / rho(i) * 3.0 / 4.0 / pi ) ** ( 1.0 / 3.0 ) * 2.0
        diam(i) = ( M(i) /  N(i)  * rhod / rho(i) * 3.d0 / 4.d0 / pi ) ** ( 1.d0 / 3.d0 ) * 2.d0
        diam(i) = MAX( diam1(i), MIN( diam1(i+1), diam(i) ))
        ! Set terminal velocity, store max Courant number
        !!! coefficients of velocity formula
        IF ( diam(i) >= 0.313d-5 .AND. diam(i) <= 0.100d-3 ) THEN !cloud(1-12) + rain(1-3)
          bet = 2.d0 / 3.d0
          alp = 0.45795d+06
        ELSE IF (  diam(i) > 0.100d-3 .AND. diam(i) <= 0.800d-3 ) THEN !rain(4-)
          bet = 1.d0 / 3.d0
          alp = 4.962d+03
        ELSE IF ( diam(i) > 0.800d-3 .AND. diam(i) <= 0.403d-2 ) THEN !rain
          bet = 1.d0 / 6.d0
          alp = 1.732d+03
        ELSE IF ( diam(i) > 0.403d-2 ) THEN !rain
          bet = 0.d0
          alp = 917.d0
        ENDIF
        if ( min( M(i), (N(i)/rhod*1000.d0) ) > 0.d0 ) THEN !  mass (g/g); number (#/g) 
          am = MAX( mass1(i), MIN( mass1(i+1), (M(i)/N(i)*rhod) ) ) * 1.d3 !mass for 1 particle in g
          fvel(i) =  alp * am ** bet * 0.01d0  ! cm/s -> m/s (factor 0.01) 
        else
          fvel(i) = 0.0d0
        endif
      enddo
     
    endif !isc
    !
  return
  end subroutine hydro_samsbm_2mom
  !
  subroutine GetPolarimetricInfofromLUT(isc,conf,inp_elev,ww,temp,nd,diam,NN,rho,fvel,Zhh,Zvv,Zvh,RHOhvc,&
                                        DVh,dDVh,Dopp,Kdp,Adp,Ah,Av,diff_back_phase,zhh_d,zvh_d,zvv_d)
  Use crsim_mod
  Use crsim_luts_mod
  Use phys_param_mod, ONLY: pi
  !
  Implicit None
  Integer, intent(In)        :: isc
  Type(conf_var),Intent(in)  :: conf
  real*8,Intent(in)          :: inp_elev ! degrees 
  Real*8, Intent(In)         :: ww       ![m/s] vertical air velocity
  real*8,Intent(in)          :: temp     !K
  integer,Intent(in)         :: nd
  real*8,Intent(in)          :: diam(nd) ! m
  real*8,Intent(in)          :: NN(nd)   ! 1/m^3
  real*8,Intent(in)          :: rho(nd)  ! kg/m^3 
  real*8,Intent(in)          :: fvel(nd) ! m/s
  real*8,Intent(Out)         :: Zhh,Zvv,Zvh,RHOhvc ! mm^6/m^3
  real*8,Intent(Out)         :: DVh  ! mm^6/m^3 * m/s     ! assuming ww=0
  real*8,Intent(Out)         :: dDVh ! mm^6/m^3 * (m/s)^2 ! assuming ww=0 
  real*8,Intent(Out)         :: Dopp ! mm^6/m^3 * m/s     ! ww from the WRF
  real*8,Intent(Out)         :: Kdp ! deg/km
  real*8,Intent(Out)         :: Adp ! dB/km
  real*8,Intent(Out)         :: Ah  ! dB/km
  real*8,Intent(Out)         :: Av  ! dB/km
  real*8,Intent(Out)         :: diff_back_phase ! deg
  real*8,Intent(Out)          :: zhh_d(nd),zvh_d(nd),zvv_d(nd)   ! mm^6 for each size
  !
  !
  real*8                      :: fff,ttt,elev
  character                   :: frq_str*4,t_str*5,el_str*4,rho_str*5
  character                   :: LutFileName*340
  !
  integer                     :: ir,ll
  ! 
  integer                          :: lnd
  real*8                           :: ldmin,ldmax,lddiam
  real*8,Dimension(:),Allocatable  :: ldiam ! m
  !
  integer                          :: nden
  real*8,Dimension(:),Allocatable  :: lden ! kg/m^3
  real*8                           :: den ! kg/m^3
  !
  integer                     :: ii, ii_start
  integer                     :: j1,j2,jmax
  integer                     :: j,nj,jj
  integer                     :: k,k1,k2
  integer                     :: ir1,ir2,nir
  
  complex*16,Dimension(:),Allocatable       :: lsb11,lsb22
  complex*16,Dimension(:),Allocatable       :: lsf11,lsf22
  real*8,Dimension(:),Allocatable           :: laoverb1,ldiam1
  integer,Dimension(:),Allocatable          :: mind
  !
  complex*16,Dimension(:),Allocatable       :: sb11,sb22
  complex*16,Dimension(:),Allocatable       :: sf11,sf22
  real*8,Dimension(:),Allocatable           :: aoverb,diam_check
  
  integer,parameter           :: id=1
  !
  real*8               :: lam_mm ! radar wavelength mm
  integer              :: horientID ! choice for orientation distribution
  real*8               :: sigma ! used only of horientID=3 
  !
  real*8               :: diff,diff0
  !
  real*8               :: abs_elev !absolute value for negative elevation by oue 20190911
  !
    !
    j2=0
    j1=1
    !---------------------------------------------
    !---------------------------------------------
    ! changed inp_elev to abs(inp_elev) for negative elevation by oue 20190911
    abs_elev = abs(inp_elev)
    !---------------------------------------------
    !----------------------------------------------
    ! orientation distribution user parameters
    horientID=conf%horientID(isc)
    if (horientID==3) sigma=conf%sigma(isc)
    !
    lam_mm=299.7925d0/conf%freq ! mm
    !---------------------------------------------
    !----------------------------------------------
    !! - get the frequency
    call get_lut_str_var(n_lfreq,lfreq,conf%freq,fff)
    !---------------------------------------------
    ! get the temperature
    !-------------------------------------
    if (isc==1) then! CLOUD
      call get_lut_str_var(n_ltemp_cld,ltemp_cld,temp+273.15d0,ttt)
    end if !- if cloud
    !------------------------------------
    if (isc==2) then! RAIN
      call get_lut_str_var(n_ltemp_rain,ltemp_rain,temp+273.15d0,ttt)
    end if !- if rain
    !!------------------------------------
    if (isc>2) then ! ICE,SNOW,GRAUPEL,HAIL
      ttt=ltemp_ice(1)
    endif
    if (isc==7) then! DRizzle 
      call get_lut_str_var(n_ltemp_cld,ltemp_cld,temp+273.15d0,ttt)
    end if !- if driz
    !
    !--------------------------------------------
    ! get the elevation
    if (isc>1) then
       if((inp_elev<=90.d0) .and. (inp_elev>=0.d0)) then
       call get_lut_str_var(n_lelev,lelev,inp_elev,elev)
       endif
       if((inp_elev>=-90.d0) .and. (inp_elev<0.d0)) then
       call get_lut_str_var(n_lelev,lelev,abs_elev,elev) ! changed inp_elev to abs(inp_elev) for negative elevation by oue 20190911
       endif
    else
      elev=lelev(n_lelev)
    endif
    if (isc==7) then
      elev=lelev(n_lelev)
    endif
    !
    !----------------------------------------------
    call get_lut_str(fff,ttt,elev,frq_str,t_str,el_str)
    !-------------------------------------------
    ! select the LUTs densities for hydrometeor isc !modified by oue 2017/07/21 for RAMS to include drizzle(isc=7) and aggregates(isc=8)
    if ((isc<=2) .or. (isc==7)) then ! cloud, rain, drizzle
      nden=n_lden_cld 
      Allocate(lden(nden))
      lden=lden_cld
    endif
    if (isc==3) then ! ice
      nden=n_lden_ice
      Allocate(lden(nden))
      lden=lden_ice
    endif
    if ((isc==4) .or. (isc==8)) then ! snow, aggregates
      nden=n_lden_snow
      Allocate(lden(nden))
      lden=lden_snow
    endif
    if ((isc==5) .or. (isc==6)) then ! graupel,hail
      nden=n_lden_graup
      Allocate(lden(nden))
      lden=lden_graup
    endif
    !-------------------------------------------
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! LUT diameters
    call  hydro_info(isc,ldmin,ldmax,lnd)
    !
    Allocate(ldiam(lnd))
    
    do ir=1,lnd
      ldiam(ir)=ldmin*1.d-3+dble(ir-1)*(ldmax-ldmin)/dble(lnd-1)*1.d-3 ! diameter in mm
    enddo
    lddiam=ldiam(10)-ldiam(9) ! ddiam in mm
    lddiam=lddiam*1.d-3 ! lddiam in m
    ldiam=ldiam*1.d-3 ! diameter in m 
    !
    !-------------------------------------------------------------------------------------------------
    ! Allocate variables from LUTS 
    Allocate(diam_check(nd),aoverb(nd),sb11(nd),sb22(nd),sf11(nd),sf22(nd))
    diam_check=0.d0; aoverb=0.d0; sb11=(0.d0,0.d0); sb22=(0.d0,0.d0); sf11=(0.d0,0.d0); sf22=(0.d0,0.d0)
    !-------------------------------------------------------------------------------------------------
    j2=0
    
    211 continue 
     
    ii_start=j2+1
    
     
    M1:do ii=ii_start,nd  ! for each hydrometeor bin diameter 
      !--------------------------------------------
      ! for bin density rho(ii) find the closest density from LUTs
      call get_lut_str_var(nden,lden,rho(ii),den)
      !for the bin density den determine the string rho_str
      call get_lut_den_str(den,rho_str)
      !----------------------------------------------------------------
      !
      ! now we can construct the LUT filename
      LutFileName=Trim(adjustl(conf%hydro_luts(isc)))//&
                 '_fr'//frq_str//'GHz_'//'t'//t_str//'_rho'//rho_str//'_el'//el_str//'.dat'
     
      LutFileName=Trim(Adjustl(conf%path_back_scatt_libs))//&
                  Trim(Adjustl(conf%hydro_luts(isc)))//'/'//&
                  Trim(Adjustl(LutFileName))
      !----------------------------------------------
      
      ! find the bins interval j1->j2 having the same density
      j1=ii
      j2=Min(ii+1,nd)
      jmax=j2      
      if (rho(jmax)==rho(j1)) then
        k1=j2
        k2=nd
        A1: do k=k1,k2
          if (rho(k)==rho(j1)) then
            j2=k
          else
            exit A1
          endif
        enddo A1
      endif
      !-------
      !---- find min,max diameter from luts covering that range   
      ir1=1 ; ir2=lnd  
      B1:do ir=1,lnd
        if (ldiam(ir)<=diam(j1)) then
          ir1=ir
        else
          exit B1
        endif
      enddo B1
      !
      B2:do ir=ir1,lnd
        if (ldiam(ir)<=diam(j2)) then 
          ir2=ir
        else
          exit B2
        endif
      enddo B2
      !
      nir=ir2-ir1+1
      !----------------------------------------------------------
      exit M1
    enddo M1
    !
    !
    jj=1
    ! now we read lookup tables and extract nir look-up variables in desired diameter range
    ! from ldiam(ir1) to ldiam(ir2)
    !-----------------------------------------------------------------------------------
    Allocate(ldiam1(nir),laoverb1(nir),lsb11(nir),lsb22(nir),lsf11(nir),lsf22(nir))
    ldiam1=0.d0; laoverb1=0.d0; lsb11=(0.d0,0.d0); lsb22=(0.d0,0.d0); lsf11=(0.d0,0.d0); lsf22=(0.d0,0.d0)
     
    !oue PRELOAD P3       
    !call read_luts(Trim(LutFileName),inp_elev,ir1,ir2,nir,ldiam1,laoverb1,lsb11,lsb22,lsf11,lsf22)
    ! changed inp_elev to abs(inp_elev) for negative elevation by oue 20190911
#ifdef __PRELOAD_LUT__
    call load_luts(Trim(LutFileName),abs_elev,ir1,ir2,nir,ldiam1,laoverb1,lsb11,lsb22,lsf11,lsf22)
#else
    call read_luts(Trim(LutFileName),abs_elev,ir1,ir2,nir,ldiam1,laoverb1,lsb11,lsb22,lsf11,lsf22)
#endif
!#ifdef __PRELOAD_LUT__
!    call load_luts(Trim(LutFileName),inp_elev,ir1,ir2,nir,ldiam1,laoverb1,lsb11,lsb22,lsf11,lsf22)
!#else
!    call read_luts(Trim(LutFileName),inp_elev,ir1,ir2,nir,ldiam1,laoverb1,lsb11,lsb22,lsf11,lsf22)
!#endif
     
    ldiam1=ldiam1*1.d-3  ! mm -> m
     
    nj=j2-j1+1
    Allocate(mind(nj))
    mind=-1
    !------------------------
    ii=0
    do j=j1,j2
      ii=ii+1
      D1:do ir=1,nir
        if(ldiam1(ir)<=diam(j)) then
          jj=ir 
        else
          exit D1
        endif
      enddo D1
      diff0=1.d+17
      do ll=max(jj-1,1),min(jj+1,nir)
        diff=dabs(ldiam1(ll)-diam(j))
        if (diff<diff0) then
          diff0=diff
          jj=ll
        endif
      enddo
      mind(ii)=jj
    enddo ! do j=j1,j2      
    !------------------------
    ii=0
    do j=j1,j2
      ii=ii+1
      diam_check(j)=ldiam1(mind(ii))
      aoverb(j)=laoverb1(mind(ii))
      sb11(j)=lsb11(mind(ii)) ! fa
      sb22(j)=-lsb22(mind(ii))! fb
      !
      sf11(j)=lsf11(mind(ii))
      sf22(j)=lsf22(mind(ii))
      !
    enddo ! do j
    Deallocate(mind)
    !   
    Deallocate(ldiam1,laoverb1,lsb11,lsb22,lsf11,lsf22)
    ! 
    IF (j2/=nd) go to 211
     
    Deallocate(ldiam)
    Deallocate(lden)
    ! at this point we have stored all information from LUTs at desired sizes diam(nd) 
    ! and we can compute the following polarimetric vars
    !
    !write(*,*) 'compute polarimetric variables'
    
    IF (conf%radID/=1) THEN  ! POLARIMETRIC
      call compute_polarim_vars(inp_elev,horientID,sigma,nd,lam_mm,sb11,sb22,sf11,sf22,NN,fvel,ww,&
                                Zhh,Zvv,Zvh,RHOhvc,DVh,Dopp,dDVh,Kdp,Ah,Av,Adp,diff_back_phase,zhh_d,zvh_d,zvv_d)
    else
      Zvv=0.d0; Zvh=0.d0; Kdp=0.d0; Adp=0.d0; Av=0.d0; Adp=0.d0; zvh_d=0.0; zvv_d=0.0;
      call compute_rad_forward_vars(inp_elev,horientID,sigma,nd,lam_mm,sb11,sb22,sf11,sf22,NN,fvel,ww,&
                                    Zhh,DVh,Dopp,dDVh,Ah,zhh_d)
    endif
    !----------------------------------------------
    Deallocate(diam_check,aoverb,sb11,sb22,sf11,sf22)
    !
  !------------------------------------------------------------------------------------------------------
  return
  end subroutine GetPolarimetricInfofromLUT
     !
  subroutine GetPolarimetricInfofromLUT_cdws(isc,conf,inp_elev,ww,temp,nd,diam,NN,rho_const,fvel,Zhh,Zvv,Zvh,RHOhvc,&
                                             DVh,dDVh,Dopp,Kdp,Adp,Ah,Av,diff_back_phase,zhh_d,zvh_d,zvv_d)
  Use crsim_mod
  Use crsim_luts_mod
  !
  ! for constant hydrometeor density with size
  Implicit None
  Integer, intent(In)        :: isc
  Type(conf_var),Intent(in)  :: conf
  real*8,Intent(in)          :: inp_elev ! degrees 
  Real*8, Intent(In)         :: ww       ![m/s] vertical air velocity
  real*8,Intent(in)          :: temp     !K
  integer,Intent(in)         :: nd
  real*8,Intent(in)          :: diam(nd)    ! m
  real*8,Intent(in)          :: NN(nd)      ! 1/m^3
  real*8,Intent(in)          :: rho_const   ! kg/m^3 
  real*8,Intent(in)          :: fvel(nd)    ! m/s
  real*8,Intent(Out)         :: Zhh,Zvv,Zvh,RHOhvc ! mm^6/m^3
  real*8,Intent(Out)         :: DVh  ! mm^6/m^3 * m/s     ! assuming ww=0
  real*8,Intent(Out)         :: dDVh ! mm^6/m^3 * (m/s)^2 ! assuming ww=0 
  real*8,Intent(Out)         :: Dopp ! mm^6/m^3 * m/s     ! ww from the WRF
  real*8,Intent(Out)         :: Kdp  ! deg/km
  real*8,Intent(Out)         :: Adp  ! dB/km
  real*8,Intent(Out)         :: Ah   ! dB/km
  real*8,Intent(Out)         :: Av   ! dB/km
  real*8,Intent(Out)         :: diff_back_phase ! [deg]
  real*8,Intent(Out)          :: zhh_d(nd),zvh_d(nd),zvv_d(nd)   ! mm^6 for each size
  !
  real*8                      :: fff,ttt,elev
  character                   :: frq_str*4,t_str*5,el_str*4,rho_str*5
  character                   :: LutFileName*340
  !
  integer                          :: lnd
  real*8                           :: ldmin,ldmax
  !
  integer                          :: nden
  real*8,Dimension(:),Allocatable  :: lden ! kg/m^3
  real*8                           :: den  ! kg/m^3
  !
  complex*16,Dimension(:),Allocatable       :: sb11,sb22
  complex*16,Dimension(:),Allocatable       :: sf11,sf22
  real*8,Dimension(:),Allocatable           :: aoverb,diam_check
  
  integer,parameter           :: id=1
  !
  real*8               :: lam_mm ! radar wavelength mm
  integer              :: horientID ! choice for orientation distribution
  real*8               :: sigma ! used only of horientID=3 
  !
  real*8               :: abs_elev !absolute value for negative elevation by oue 20190911
  !                                                                                                         
    !
  !---------------------------------------------
  !---------------------------------------------
  ! changed inp_elev to abs(inp_elev) for negative elevation 20190911
  abs_elev = abs(inp_elev)
  !---------------------------------------------
  !----------------------------------------------
    ! orientation distribution user parameters
    horientID=conf%horientID(isc)
    if (horientID==3) sigma=conf%sigma(isc)
    !
    lam_mm=299.7925d0/conf%freq ! mm
    !---------------------------------------------
    !----------------------------------------------
    !! - get the frequency
    call get_lut_str_var(n_lfreq,lfreq,conf%freq,fff)
    !write(*,*) 'LUTs: Using frequency', fff
    !---------------------------------------------
    ! get the temperature
    !-------------------------------------
    if (isc==1) then! CLOUD
      call get_lut_str_var(n_ltemp_cld,ltemp_cld,temp+273.15d0,ttt)
    end if !- if cloud
    !------------------------------------
    if (isc==2) then! RAIN
      call get_lut_str_var(n_ltemp_rain,ltemp_rain,temp+273.15d0,ttt)
    end if !- if rain
    !!------------------------------------
    if (isc>2) then ! ICE,SNOW,GRAUPEL,HAIL (aggregates for RAMS)
      ttt=ltemp_ice(1)
    endif
    !-------------------------------------
    if (isc==7) then ! Drizzle added by oue 2017/07/21 for RAMS
      call get_lut_str_var(n_ltemp_cld,ltemp_cld,temp+273.15d0,ttt)
    end if !- if drizzle
    !
    !write(*,*) 'LUTs: Using temperature', ttt
    !--------------------------------------------
    !
    ! get the elevation
    if (isc>1) then
      if((inp_elev<=90.d0) .and. (inp_elev>=0.d0)) then
         call get_lut_str_var(n_lelev,lelev,inp_elev,elev)
      endif
      if((inp_elev>=-90.d0) .and. (inp_elev<0.d0)) then
      call get_lut_str_var(n_lelev,lelev,abs_elev,elev)  ! changed inp_elev to abs(inp_elev) for negative elevation by oue 20190911
      endif
    else
      elev=lelev(n_lelev)
    endif
    if (isc==7) then !Drizzle added by oue 2017/07/21 for RAMS
      elev=lelev(n_lelev)
    endif
    !write(*,*) 'fff, ttt, elev, inp_elev : ', fff, ttt, elev, inp_elev
    !
    !----------------------------------------------
    call get_lut_str(fff,ttt,elev,frq_str,t_str,el_str)
    !-------------------------------------------
    !
    ! hydrometeor density
    ! select the LUTs densities for hydrometeor isc !modified by oue 2017/07/21 for RAMS to include drizzle(isc=7) and aggregates(isc=8)
    if ((isc<=2) .or. (isc==7)) then ! cloud, rain, drizzle
      nden=n_lden_cld
      Allocate(lden(nden))
      lden=lden_cld
    endif
    if (isc==3) then ! ice
      nden=n_lden_ice
      Allocate(lden(nden))
      lden=lden_ice
    endif
    if ((isc==4) .or. (isc==8)) then ! snow, aggregates
      nden=n_lden_snow
      Allocate(lden(nden))
      lden=lden_snow
    endif
    if ((isc==5) .or. (isc==6)) then ! graupel,hail
      nden=n_lden_graup
      Allocate(lden(nden))
      lden=lden_graup
    endif
    !-------------------------------------------
    !
    ! for rho_const find the closest density from LUTs
    call get_lut_str_var(nden,lden,rho_const,den)
    !for the bin density den determine the string rho_str
    call get_lut_den_str(den,rho_str)
    !
    Deallocate(lden)
    !----------------------------------------------------------------
    !
    ! now we can construct the LUT filename
    !
    LutFileName=Trim(adjustl(conf%hydro_luts(isc)))//&
               '_fr'//frq_str//'GHz_'//'t'//t_str//'_rho'//rho_str//'_el'//el_str//'.dat'
    LutFileName=Trim(Adjustl(conf%path_back_scatt_libs))//&
                Trim(Adjustl(conf%hydro_luts(isc)))//'/'//&
                Trim(Adjustl(LutFileName))
    !write(*,*) 'LutFileName : ', rho_const,LutFileName
    !
    !----------------------------------------------
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! LUT diameters
    call  hydro_info(isc,ldmin,ldmax,lnd)
    !
    !-------------------------------------------------------------------------------------------------
    ! Allocate variables from LUTS 
    Allocate(diam_check(nd),aoverb(nd),sb11(nd),sb22(nd),sf11(nd),sf22(nd))
    diam_check=0.d0; aoverb=0.d0; sb11=(0.d0,0.d0); sb22=(0.d0,0.d0)
    sf11=(0.d0,0.d0); sf22=(0.d0,0.d0)
    !-------------------------------------------------------------------------------------------------
    ! changed inp_elev to abs(inp_elev) for negative elevation by oue 20190911
#ifdef __PRELOAD_LUT__
    call load_luts(Trim(LutFileName),abs_elev,1,nd,nd,diam_check,aoverb,sb11,sb22,sf11,sf22)
#else
    call read_luts(Trim(LutFileName),abs_elev,1,nd,nd,diam_check,aoverb,sb11,sb22,sf11,sf22)
#endif
!#ifdef __PRELOAD_LUT__
!    call load_luts(Trim(LutFileName),inp_elev,1,nd,nd,diam_check,aoverb,sb11,sb22,sf11,sf22)
!#else
!    call read_luts(Trim(LutFileName),inp_elev,1,nd,nd,diam_check,aoverb,sb11,sb22,sf11,sf22)
!#endif
    diam_check=diam_check*1.d-3  ! mm -> m
    sb22=-sb22
    
    !write(*,*) 'compute polarimetric variables'
    !write(*,*) 'compute polarimetric variables',lam_mm,NN(1),fvel(1),ww
    
    IF (conf%radID/=1) THEN  ! POLARIMETRIC
    call compute_polarim_vars(inp_elev,horientID,sigma,nd,lam_mm,sb11,sb22,sf11,sf22,NN,fvel,ww,&
                               Zhh,Zvv,Zvh,RHOhvc,DVh,Dopp,dDVh,Kdp,Ah,Av,Adp,diff_back_phase,zhh_d,zvh_d,zvv_d)
    else
    Zvv=0.d0; Zvh=0.d0; Kdp=0.d0; Adp=0.d0; Av=0.d0; Adp=0.d0; zvh_d=0.0; zvv_d=0.0
    call compute_rad_forward_vars(inp_elev,horientID,sigma,nd,lam_mm,sb11,sb22,sf11,sf22,NN,fvel,ww,&
                              Zhh,DVh,Dopp,dDVh,Ah,zhh_d)
    endif
     
    !write(*,*) 'zhh : ', Zhh
    
    Deallocate(diam_check,aoverb,sb11,sb22,sf11,sf22)
     
    !
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
  return
  end subroutine GetPolarimetricInfofromLUT_cdws
  !
  !
  subroutine GetPolarimetricInfofromLUT_sbm(isc,conf,inp_elev,ww,temp,nd,diam,NN,rho,fvel,Zhh,Zvv,Zvh,RHOhvc,&
                                           DVh,dDVh,Dopp,Kdp,Adp,Ah,Av,diff_back_phase,zhh_d,zvh_d,zvv_d)
  Use crsim_mod
  Use crsim_luts_mod
  Use phys_param_mod, ONLY: pi
  !
  ! for SBM - 2 possibilities 1) density const with sze and 2) density variable with size
  !
  Implicit None
  !
  Integer, intent(In)        :: isc
  Type(conf_var),Intent(in)  :: conf
  real*8,Intent(in)          :: inp_elev ! degrees 
  Real*8, Intent(In)         :: ww       ![m/s] vertical air velocity
  real*8,Intent(in)          :: temp     !K
  integer,Intent(in)         :: nd
  real*8,Intent(in)          :: diam(nd) ! m
  real*8,Intent(in)          :: NN(nd)   ! 1/m^3
  real*8,Intent(in)          :: rho(nd)  ! kg/m^3 
  real*8,Intent(in)          :: fvel(nd) ! m/s
  real*8,Intent(Out)         :: Zhh,Zvv,Zvh,RHOhvc ! mm^6/m^3
  real*8,Intent(Out)         :: DVh  ! mm^6/m^3 * m/s     ! assuming ww=0
  real*8,Intent(Out)         :: dDVh ! mm^6/m^3 * (m/s)^2 ! assuming ww=0 
  real*8,Intent(Out)         :: Dopp ! mm^6/m^3 * m/s     ! ww from the WRF
  real*8,Intent(Out)         :: Kdp ! deg/km
  real*8,Intent(Out)         :: Adp ! dB/km
  real*8,Intent(Out)         :: Ah  ! dB/km
  real*8,Intent(Out)         :: Av  ! dB/km
  real*8,Intent(Out)         :: diff_back_phase ! deg
  real*8,Intent(Out)          :: zhh_d(nd),zvh_d(nd),zvv_d(nd)   ! mm^6 for each size
  !
  real*8                      :: fff,ttt,elev
  character                   :: frq_str*4,t_str*5,el_str*4,rho_str*5
  character                   :: LutFileName*340
  !
  integer                     :: ir,ll
  ! 
  integer                          :: lnd
  real*8                           :: ldmin,ldmax
  real*8,Dimension(:),Allocatable  :: ldiam ! m
  !
  integer                          :: nden
  real*8,Dimension(:),Allocatable  :: lden ! kg/m^3
  real*8                           :: den ! kg/m^3
  !
  integer                     :: ii, ii_start
  integer                     :: j1,j2,jmax
  integer                     :: j,nj,jj
  integer                     :: k,k1,k2
  integer                     :: ir1,ir2,nir
  
  complex*16,Dimension(:),Allocatable       :: lsb11,lsb22
  complex*16,Dimension(:),Allocatable       :: lsf11,lsf22
  real*8,Dimension(:),Allocatable           :: laoverb1,ldiam1
  integer,Dimension(:),Allocatable          :: mind
  !
  complex*16,Dimension(:),Allocatable       :: sb11,sb22
  complex*16,Dimension(:),Allocatable       :: sf11,sf22
  real*8,Dimension(:),Allocatable           :: aoverb,diam_check
  
  integer,parameter           :: id=1
  !
  real*8               :: lam_mm ! radar wavelength mm
  integer              :: horientID ! choice for orientation distribution
  real*8               :: sigma ! used only of horientID=3 
  !
  real*8               :: diff,diff0
  !
  real*8               :: abs_elev !absolute value for negative elevation by oue 20190911
  !                                                                                                         
    !
    jj=1
    j1=1
    j2=0
    !---------------------------------------------
    !---------------------------------------------
    ! changed inp_elev to abs(inp_elev) for negative elevation by oue 20190911
    abs_elev = abs(inp_elev)
    !----------------------------------------------
    !---------------------------------------------
    ! orientation distribution user parameters
    horientID=conf%horientID(isc)
    if (horientID==3) sigma=conf%sigma(isc)
    !
    lam_mm=299.7925d0/conf%freq ! mm
    !---------------------------------------------
    !----------------------------------------------
    !! - get the frequency
    call get_lut_str_var(n_lfreq,lfreq,conf%freq,fff)
    !write(*,*) 'LUTs: Using frequency', fff
    !---------------------------------------------
    ! get the temperature
    !-------------------------------------
    if (isc==1) then! CLOUD
      call get_lut_str_var(n_ltemp_cld,ltemp_cld,temp+273.15d0,ttt)
    end if !- if cloud
    !------------------------------------
    if (isc==2) then! RAIN
      call get_lut_str_var(n_ltemp_rain,ltemp_rain,temp+273.15d0,ttt)
    end if !- if rain
    !!------------------------------------
    if (isc>2) then ! ICE,SNOW,GRAUPEL,HAIL
      ttt=ltemp_ice(1)
    endif
    !
    !write(*,*) 'LUTs: Using temperature', ttt
    !--------------------------------------------
    !
    ! get the elevation
    if (isc>1) then
       if((inp_elev<=90.d0) .and. (inp_elev>=0.d0)) then
          call get_lut_str_var(n_lelev,lelev,inp_elev,elev)
       endif
       if((inp_elev>=-90.d0) .and. (inp_elev<0.d0)) then
       call get_lut_str_var(n_lelev,lelev,abs_elev,elev) ! changed inp_elev to abs(inp_elev) for negative elevation by oue 20190911
       endif
    else
      elev=lelev(n_lelev)
    endif
    !
    !----------------------------------------------
    call get_lut_str(fff,ttt,elev,frq_str,t_str,el_str)
    !-------------------------------------------
    !-------------------------------------------
    ! select the LUTs densities for hydrometeor isc
    if (isc<=2) then ! cloud, rain
      nden=n_lden_cld
      Allocate(lden(nden))
      lden=lden_cld
    endif
    if (isc==3) then ! ice
      nden=n_lden_ice
      Allocate(lden(nden))
      lden=lden_ice
    endif
    if (isc==4) then ! snow
      nden=n_lden_snow
      Allocate(lden(nden))
      lden=lden_snow
    endif
    if (isc>=5) then ! graupel
      nden=n_lden_graup
      Allocate(lden(nden))
      lden=lden_graup
    endif
    !-------------------------------------------
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! LUT diameters
    call  hydro_info(isc,ldmin,ldmax,lnd)
    !
    !-------------------------------------------------------------------------------------------------
    ! Allocate variables from LUTS at input bin sizes
    Allocate(diam_check(nd),aoverb(nd),sb11(nd),sb22(nd),sf11(nd),sf22(nd))
    diam_check=0.d0; aoverb=0.d0; sb11=(0.d0,0.d0); sb22=(0.d0,0.d0); sf11=(0.d0,0.d0); sf22=(0.d0,0.d0)
    !-------------------------------------------------------------------------------------------------
    ! 
    ! IF DENSITY CONSTANT WITH SIZE
    IF (MinVal(rho)==MaxVal(rho)) THEN 
      !
      !-------------------------------------------------------------------------------------------------
      !--------------------------------------------
      ! for bin density rho(ii) find the closest density from LUTs
      call get_lut_str_var(nden,lden,rho(1),den)
      !for the bin density den determine the string rho_str
      call get_lut_den_str(den,rho_str)
      !----------------------------------------------------------------
      !
      ! now we can construct the LUT filename
      !
      LutFileName=Trim(adjustl(conf%hydro_luts(isc)))//&
                 '_fr'//frq_str//'GHz_'//'t'//t_str//'_rho'//rho_str//'_el'//el_str//'.dat'
    
      LutFileName=Trim(Adjustl(conf%path_back_scatt_libs))//&
                  Trim(Adjustl(conf%hydro_luts(isc)))//'/'//&
                  Trim(Adjustl(LutFileName))
      !----------------------------------------------
      ! read LUTS
      Allocate(ldiam1(lnd),laoverb1(lnd),lsb11(lnd),lsb22(lnd),lsf11(lnd),lsf22(lnd))
      ldiam1=0.d0; laoverb1=0.d0; lsb11=(0.d0,0.d0); lsb22=(0.d0,0.d0);
      lsf11=(0.d0,0.d0); lsf22=(0.d0,0.d0)
      
      !call read_luts(Trim(LutFileName),inp_elev,1,lnd,lnd,ldiam1,laoverb1,lsb11,lsb22,lsf11,lsf22)
      ! changed inp_elev to abs(inp_elev) 20190911
#ifdef __PRELOAD_LUT__
      call load_luts(Trim(LutFileName),abs_elev,1,lnd,lnd,ldiam1,laoverb1,lsb11,lsb22,lsf11,lsf22)
#else
      call read_luts(Trim(LutFileName),abs_elev,1,lnd,lnd,ldiam1,laoverb1,lsb11,lsb22,lsf11,lsf22)
#endif
!#ifdef __PRELOAD_LUT__
!      call load_luts(Trim(LutFileName),inp_elev,1,lnd,lnd,ldiam1,laoverb1,lsb11,lsb22,lsf11,lsf22)
!#else
!      call read_luts(Trim(LutFileName),inp_elev,1,lnd,lnd,ldiam1,laoverb1,lsb11,lsb22,lsf11,lsf22)
!#endif
      ldiam1=ldiam1*1.d-3  ! mm -> m
      !
      ! find the closest LUT diametrs
      Allocate(mind(nd))
      mind=-1
      !------------------------
      ii=0
      do j=1,nd
        ii=ii+1
        AD1:do ir=1,lnd
          if(ldiam1(ir)<=diam(j)) then
            jj=ir
          else
            exit AD1
          endif
        enddo AD1
        diff0=1.d+17
        do ll=max(jj-1,1),min(jj+1,lnd)
          diff=dabs(ldiam1(ll)-diam(j))
          if (diff<diff0) then
            diff0=diff
            jj=ll
          endif
        enddo
        mind(ii)=jj
      enddo ! do j=1,nd    
      !------------------------
      ! get the LUT info at the LUT diameters closest to the input diams
      !------------------------
      ii=0
      do j=1,nd
        ii=ii+1
        diam_check(j)=ldiam1(mind(ii))
        aoverb(j)=laoverb1(mind(ii))
        !
        sb11(j)=lsb11(mind(ii)) ! fa
        sb22(j)=-lsb22(mind(ii))! fb
        !
        sf11(j)=lsf11(mind(ii))
        sf22(j)=lsf22(mind(ii))
        !
      enddo ! do j
      !
      Deallocate(mind)
      !
      Deallocate(ldiam1,laoverb1,lsb11,lsb22,lsf11,lsf22)
      !
      !
    ELSE ! IF DENSITY IS NOT CONSTANT WITH SIZE
      !--------------------------------------------------------------------------------------
      !
      ! reconstruct lut bin diams
      Allocate(ldiam(lnd))
      do ir=1,lnd
        ldiam(ir)=ldmin*1.d-3+dble(ir-1)*(ldmax-ldmin)/dble(lnd-1)*1.d-3 ! diameter in mm
      enddo
      ldiam=ldiam*1.d-3 ! diameter in m 
     
      j2=0
     
      211 continue
     
      ii_start=j2+1
      
      M1:do ii=ii_start,nd  ! for each hydrometeor bin diameter 
  
        !!---------------------------------------------------
        !--------------------------------------------
        ! for bin density rho(ii) find the closest density from LUTs
        call get_lut_str_var(nden,lden,rho(ii),den)
        !for the bin density den determine the string rho_str
        call get_lut_den_str(den,rho_str)
        !----------------------------------------------------------------
        !
        ! now we can construct the LUT filename
        
        LutFileName=Trim(adjustl(conf%hydro_luts(isc)))//&
               '_fr'//frq_str//'GHz_'//'t'//t_str//'_rho'//rho_str//'_el'//el_str//'.dat'
    
        LutFileName=Trim(Adjustl(conf%path_back_scatt_libs))//&
                Trim(Adjustl(conf%hydro_luts(isc)))//'/'//&
                Trim(Adjustl(LutFileName))
        !----------------------------------------------
        ! find the bins interval j1->j2 having the same density
        
        j1=ii
        j2=Min(ii+1,nd)
        jmax=j2
        if (rho(jmax)==rho(j1)) then
          k1=j2
          k2=nd
          A1: do k=k1,k2
            if (rho(k)==rho(j1)) then
              j2=k
            else
              exit A1
            endif
          enddo A1
        endif
        !-------
        !---- find min,max diameter from luts covering that range 
        ir1=1 ; ir2=lnd
        B1:do ir=1,lnd
          if (ldiam(ir)<=diam(j1)) then
            ir1=ir
          else
            exit B1
          endif
        enddo B1
        !
        B2:do ir=ir1,lnd
          if (ldiam(ir)<=diam(j2)) then
            ir2=ir
          else
            exit B2
          endif
        enddo B2
        !
        nir=ir2-ir1+1
        !----------------------------------------------------------
        exit M1
      enddo M1
      !
      !
      ! now we read lookup tables and extract nir look-up variables in desired
      ! diameter range
      ! from ldiam(ir1) to ldiam(ir2)
      !-----------------------------------------------------------------------------------
      Allocate(ldiam1(nir),laoverb1(nir),lsb11(nir),lsb22(nir),lsf11(nir),lsf22(nir))
      ldiam1=0.d0; laoverb1=0.d0
      lsb11=(0.d0,0.d0); lsb22=(0.d0,0.d0)
      lsf11=(0.d0,0.d0); lsf22=(0.d0,0.d0)
     
      !call read_luts(Trim(LutFileName),inp_elev,ir1,ir2,nir,ldiam1,laoverb1,lsb11,lsb22,lsf11,lsf22)
      ! changed inp_elev to abs(inp_elev) for negative elevation by oue 20190911
#ifdef __PRELOAD_LUT__
      call load_luts(Trim(LutFileName),abs_elev,ir1,ir2,nir,ldiam1,laoverb1,lsb11,lsb22,lsf11,lsf22)
#else
      call read_luts(Trim(LutFileName),abs_elev,ir1,ir2,nir,ldiam1,laoverb1,lsb11,lsb22,lsf11,lsf22)
#endif
!#ifdef __PRELOAD_LUT__
!      call load_luts(Trim(LutFileName),inp_elev,ir1,ir2,nir,ldiam1,laoverb1,lsb11,lsb22,lsf11,lsf22)
!#else
!      call read_luts(Trim(LutFileName),inp_elev,ir1,ir2,nir,ldiam1,laoverb1,lsb11,lsb22,lsf11,lsf22)
!#endif
      ldiam1=ldiam1*1.d-3  ! mm -> m
     
   
      nj=j2-j1+1
      Allocate(mind(nj))
      mind=-1
      !------------------------
      ii=0
      do j=j1,j2
        ii=ii+1
        D1:do ir=1,nir
          if(ldiam1(ir)<=diam(j)) then
            jj=ir
          else
            exit D1
          endif
        enddo D1
        diff0=1.d+17
        do ll=max(jj-1,1),min(jj+1,nir)
          diff=dabs(ldiam1(ll)-diam(j))
          if (diff<diff0) then
            diff0=diff
            jj=ll
          endif
        enddo
    
        mind(ii)=jj
      enddo ! do j=j1,j2      
      !------------------------
      ii=0
      do j=j1,j2
        ii=ii+1
        diam_check(j)=ldiam1(mind(ii))
        aoverb(j)=laoverb1(mind(ii))
        ! 
        sb11(j)=lsb11(mind(ii)) ! fa
        sb22(j)=-lsb22(mind(ii))! fb
        !
        sf11(j)=lsf11(mind(ii))
        sf22(j)=lsf22(mind(ii))
        !
      enddo ! do j
      Deallocate(mind)
     
      Deallocate(ldiam1,laoverb1,lsb11,lsb22,lsf11,lsf22)
     
      IF (j2/=nd) go to 211
     
     
      Deallocate(ldiam)
    ENDIF ! -------------------------------------------------------
    !
    Deallocate(lden)
     
    !!-----------------------------------------------------------------------------
    ! at this point we have stored all information from LUTs at desired sizes diam(nd) 
    ! and we can compute the following polarimetric vars
    !
    ! write(*,*) 'compute polarimetric variables'
     
    IF (conf%radID/=1) THEN  ! POLARIMETRIC
      call compute_polarim_vars(inp_elev,horientID,sigma,nd,lam_mm,sb11,sb22,sf11,sf22,NN,fvel,ww,&
                           Zhh,Zvv,Zvh,RHOhvc,DVh,Dopp,dDVh,Kdp,Ah,Av,Adp,diff_back_phase,zhh_d,zvh_d,zvv_d)
    ELSE
      Zvv=0.d0; Zvh=0.d0; Kdp=0.d0; Adp=0.d0; Av=0.d0; Adp=0.d0; zvh_d=0.0; zvv_d=0.0
      call compute_rad_forward_vars(inp_elev,horientID,sigma,nd,lam_mm,sb11,sb22,sf11,sf22,NN,fvel,ww,&
                           Zhh,DVh,Dopp,dDVh,Ah,zhh_d)
    ENDIF
    !----------------------------------------------
    
    Deallocate(diam_check,aoverb,sb11,sb22,sf11,sf22)
     
  return
  end subroutine GetPolarimetricInfofromLUT_sbm
  !
  !-----------------------------
  ! subroutine GetPolarimetricInfofromLUT_P3 is added by DW 2017/10/30 for P3
  subroutine GetPolarimetricInfofromLUT_P3(isc,conf,inp_elev,ww,temp,nd,diam,NN,rho,fvel,Zhh,Zvv,Zvh,RHOhvc,&
                                           DVh,dDVh,Dopp,Kdp,Adp,Ah,Av,diff_back_phase,zhh_d,zvh_d,zvv_d)
  Use crsim_mod
  Use crsim_luts_mod
  Use phys_param_mod, ONLY: pi
  !
  Implicit None
  Integer, intent(In)        :: isc
  Type(conf_var),Intent(in)  :: conf
  real*8,Intent(in)          :: inp_elev    ! degrees 
  Real*8, Intent(In)         :: ww          ! m/s vertical air velocity
  real*8,Intent(in)          :: temp        ! K
  integer,Intent(in)         :: nd
  real*8,Intent(in)          :: diam(nd)    ! m
  real*8,Intent(in)          :: NN(nd)      ! 1/m^3
  real*8,Intent(in)          :: rho(nd)     ! kg/m^3 
  real*8,Intent(in)          :: fvel(nd)    ! m/s
  real*8,Intent(Out)         :: Zhh,Zvv,Zvh,RHOhvc ! mm^6/m^3
  real*8,Intent(Out)         :: DVh         ! mm^6/m^3 * m/s     ! assuming ww=0
  real*8,Intent(Out)         :: dDVh        ! mm^6/m^3 * (m/s)^2 ! assuming ww=0 
  real*8,Intent(Out)         :: Dopp        ! mm^6/m^3 * m/s     ! ww from the WRF
  real*8,Intent(Out)         :: Kdp         ! deg/km
  real*8,Intent(Out)         :: Adp         ! dB/km
  real*8,Intent(Out)         :: Ah          ! dB/km
  real*8,Intent(Out)         :: Av          ! dB/km
  real*8,Intent(Out)         :: diff_back_phase ! deg
  real*8,Intent(Out)          :: zhh_d(nd),zvh_d(nd),zvv_d(nd)   ! mm^6 for each size
  real*8                     :: fff,ttt,elev
  character                  :: frq_str*4,t_str*5,el_str*4,rho_str*5
  character                  :: LutFileName*340
  !
  integer                    :: ir,ll
  ! 
  integer                          :: lnd
  real*8                           :: ldmin,ldmax,lddiam
  real*8,Dimension(:),Allocatable  :: ldiam ! m
  !
  integer                          :: nden
  real*8,Dimension(:),Allocatable  :: lden ! kg/m^3
  real*8                           :: den  ! kg/m^3
  !
  integer                     :: ii, ii_start
  integer                     :: j1,j2,jmax
  integer                     :: j,nj,jj
  integer                     :: k,k1,k2
  integer                     :: ir1,ir2,nir
  !
  complex*16,Dimension(:),Allocatable       :: lsb11,lsb22
  complex*16,Dimension(:),Allocatable       :: lsf11,lsf22
  real*8,Dimension(:),Allocatable           :: laoverb1,ldiam1
  integer,Dimension(:),Allocatable          :: mind
  !
  complex*16,Dimension(:),Allocatable       :: sb11,sb22
  complex*16,Dimension(:),Allocatable       :: sf11,sf22
  real*8,Dimension(:),Allocatable           :: aoverb,diam_check
  !
  integer,parameter                         :: id=1
  !
  real*8               :: lam_mm    ! radar wavelength mm
  integer              :: horientID ! choice for orientation distribution
  real*8               :: sigma     ! used only of horientID=3 
  real*8               :: diff,diff0
  !
  real*8               :: abs_elev !absolute value for negative elevation by oue 20190911
  !
    jj=1
    j2=0
    j1=1
    !---------------------------------------------
    !---------------------------------------------
    ! changed inp_elev to abs(inp_elev) for negative elevation by oue 20190911
    abs_elev = abs(inp_elev)
    !---------------------------------------------
    !---------------------------------------------    
    ! orientation distribution user parameters
    horientID=conf%horientID(isc)
    if (horientID==3) sigma=conf%sigma(isc)
    !
    lam_mm=299.7925d0/conf%freq ! mm
    !
    ! get the frequency
    call get_lut_str_var(n_lfreq,lfreq,conf%freq,fff)
    !
    ! get the temperature
    ttt=ltemp_smallice(1)
    !
    ! get the elevation
    if((inp_elev<=90.d0) .and. (inp_elev>=0.d0)) then
    call get_lut_str_var(n_lelev,lelev,inp_elev,elev)
    end if
    if((inp_elev>=-90.d0) .and. (inp_elev<0.d0)) then
    call get_lut_str_var(n_lelev,lelev,abs_elev,elev) ! changed inp_elev to abs(inp_elev) for negative elevation by oue 20190911
    endif
    !
    call get_lut_str(fff,ttt,elev,frq_str,t_str,el_str)
    !
    ! select the LUTs densities for hydrometeor isc
    if (isc==3) then ! small ice
      nden=n_lden_smallice
      Allocate(lden(nden))
      lden=lden_smallice
    endif
    if (isc==4) then ! unrimed ice
      nden=n_lden_unrice
      Allocate(lden(nden))
      lden=lden_unrice
    endif
    if (isc==5) then ! graupel
      nden=n_lden_graupP3   
      Allocate(lden(nden))
      lden=lden_graupP3
    endif
    if (isc==6) then ! partially rimed ice
      nden=n_lden_parice
      Allocate(lden(nden))
      lden=lden_parice
    endif
    !
    ! LUT diameters
    call  hydro_info_P3(isc,ldmin,ldmax,lnd)
    !
    Allocate(ldiam(lnd))
    !
    do ir=1,lnd
      ldiam(ir)=ldmin*1.d-3+dble(ir-1)*(ldmax-ldmin)/dble(lnd-1)*1.d-3 ! diameter in mm
    enddo
    lddiam=ldiam(2)-ldiam(1) ! ddiam in mm
    lddiam=lddiam*1.d-3 ! lddiam in m
    ldiam=ldiam*1.d-3 ! diameter in m 
    !
    ! Allocate variables from LUTS 
    Allocate(diam_check(nd),aoverb(nd),sb11(nd),sb22(nd),sf11(nd),sf22(nd))
    diam_check=0.d0; aoverb=0.d0
    sb11=(0.d0,0.d0); sb22=(0.d0,0.d0)
    sf11=(0.d0,0.d0); sf22=(0.d0,0.d0)
    !
    j2=0
    211 continue 
    ii_start=j2+1
    !
    M1:do ii=ii_start,nd  ! for each hydrometeor bin diameter 
      !
      ! for bin density rho(ii) find the closest density from LUTs
      call get_lut_str_var(nden,lden,rho(ii),den)
      !
      call get_lut_den_str(den,rho_str)
      !
      ! now we can construct the LUT filename
      LutFileName=Trim(adjustl(conf%hydro_luts(isc)))//&
                 '_fr'//frq_str//'GHz_'//'t'//t_str//'_rho'//rho_str//'_el'//el_str//'.dat'
      !
      LutFileName=Trim(Adjustl(conf%path_back_scatt_libs))//&
                  Trim(Adjustl(conf%hydro_luts(isc)))//'/'//&
                  Trim(Adjustl(LutFileName))
      !
      ! find the bins interval j1->j2 having the same density
      j1=ii
      j2=Min(ii+1,nd)
      jmax=j2      
      if (rho(jmax)==rho(j1)) then
        k1=j2
        k2=nd
        A1: do k=k1,k2
          if (rho(k)==rho(j1)) then
            j2=k
          else
            exit A1
          endif
        enddo A1
      endif
      !---- find min,max diameter from luts covering that range   
      ir1=1 ; ir2=lnd  
      B1:do ir=1,lnd
        if (ldiam(ir)<=diam(j1)) then
          ir1=ir
        else
          exit B1
        endif
      enddo B1
      !
      B2:do ir=ir1,lnd
        if (ldiam(ir)<=diam(j2)) then 
          ir2=ir
        else
          exit B2
        endif
      enddo B2
      !
      nir=ir2-ir1+1
      exit M1
    enddo M1
    !
    ! now we read lookup tables and extract nir look-up variables in desired
    ! diameter range from ldiam(ir1) to ldiam(ir2)
    !------------------------------------------------------------------------------
    Allocate(ldiam1(nir),laoverb1(nir),lsb11(nir),lsb22(nir),lsf11(nir),lsf22(nir))
    ldiam1=0.d0; laoverb1=0.d0
    lsb11=(0.d0,0.d0); lsb22=(0.d0,0.d0)
    lsf11=(0.d0,0.d0); lsf22=(0.d0,0.d0)
    !
    !oue PRELOAD P3
    ! changed inp_elev to abs(inp_elev) for negative elevation by oue 20190911
#ifdef __PRELOAD_LUT__
    call load_luts(Trim(LutFileName),abs_elev,ir1,ir2,nir,ldiam1,laoverb1,lsb11,lsb22,lsf11,lsf22)
#else
    call read_luts(Trim(LutFileName),abs_elev,ir1,ir2,nir,ldiam1,laoverb1,lsb11,lsb22,lsf11,lsf22)
#endif
!#ifdef __PRELOAD_LUT__
!    call load_luts(Trim(LutFileName),inp_elev,ir1,ir2,nir,ldiam1,laoverb1,lsb11,lsb22,lsf11,lsf22)
!#else
!    call read_luts(Trim(LutFileName),inp_elev,ir1,ir2,nir,ldiam1,laoverb1,lsb11,lsb22,lsf11,lsf22)
!#endif
     
    ldiam1=ldiam1*1.d-3  ! mm -> m
    nj=j2-j1+1
    Allocate(mind(nj))
    mind=-1
    ii=0
    do j=j1,j2
      ii=ii+1
      D1:do ir=1,nir
        if(ldiam1(ir)<=diam(j)) then
          jj=ir 
        else
          exit D1
        endif
      enddo D1
      diff0=1.d+17
      do ll=max(jj-1,1),min(jj+1,nir)
        diff=dabs(ldiam1(ll)-diam(j))
        if (diff<diff0) then
          diff0=diff
          jj=ll
        endif
      enddo
      mind(ii)=jj
    enddo ! do j=j1,j2      
  
    ii=0
    do j=j1,j2
      ii=ii+1
      diam_check(j)=ldiam1(mind(ii))
      aoverb(j)=laoverb1(mind(ii))
      !
      sb11(j)=lsb11(mind(ii)) ! fa
      sb22(j)=-lsb22(mind(ii))! fb
      !
      sf11(j)=lsf11(mind(ii))
      sf22(j)=lsf22(mind(ii))
      !
    enddo ! do j
    Deallocate(mind)
    !
    Deallocate(ldiam1,laoverb1,lsb11,lsb22,lsf11,lsf22)
   
    IF (j2/=nd) go to 211
    Deallocate(ldiam)
    Deallocate(lden)
    
    ! at this point we have stored all information from LUTs at desired sizes diam(nd) 
    ! and we can compute the following polarimetric vars
    IF (conf%radID/=1) THEN  ! POLARIMETRIC
      call compute_polarim_vars(inp_elev,horientID,sigma,nd,lam_mm,sb11,sb22,sf11,sf22,NN,fvel,ww,&
                           Zhh,Zvv,Zvh,RHOhvc,DVh,Dopp,dDVh,Kdp,Ah,Av,Adp,diff_back_phase,zhh_d,zvh_d,zvv_d)
    else
      Zvv=0.d0; Zvh=0.d0; Kdp=0.d0; Adp=0.d0; Av=0.d0; Adp=0.d0; zvh_d=0.0; zvv_d=0.0
      call compute_rad_forward_vars(inp_elev,horientID,sigma,nd,lam_mm,sb11,sb22,sf11,sf22,NN,fvel,ww,&
                                    Zhh,DVh,Dopp,dDVh,Ah,zhh_d)
    endif
    !
    Deallocate(diam_check,aoverb,sb11,sb22,sf11,sf22)
    !
  return
  end subroutine GetPolarimetricInfofromLUT_P3
  !
  !--------------------------------------------
  
     
#ifdef __PRELOAD_LUT__
     
  subroutine preload_lut_ceilo(conf)
  Use crsim_mod
  Use phys_param_mod, ONLY: pi
  Implicit None
  
  character                  :: LutFileName*340
  Type(conf_var),Intent(in)  :: conf
  integer, parameter         :: id=1
  integer                    :: ii
  !Type(lut_ceilo_file_type)  :: lut_file
  
    LutFileName='ceilo/cld_ceilo_905nm_p25.dat'
    !
    LutFileName=Trim(Adjustl(conf%path_back_scatt_libs))//&
                Trim(Adjustl(LutFileName))
    !write(*,*) "ceilo LUT filename : ", LutFileName
    !
    open(unit=id,file=LutFileName)
      read(id,*) lut_ceilo_file%tempK, lut_ceilo_file%wavelength  ! K and um
      read(id,*) lut_ceilo_file%ni_r, lut_ceilo_file%ni_c
      read(id,*) lut_ceilo_file%ndiam
      !
      Allocate(lut_ceilo_file%ldiam(lut_ceilo_file%ndiam), &
      lut_ceilo_file%qext(lut_ceilo_file%ndiam), &
      lut_ceilo_file%qback(lut_ceilo_file%ndiam))
      !
      do ii=1,lut_ceilo_file%ndiam
        read(id,*) lut_ceilo_file%ldiam(ii),lut_ceilo_file%qext(ii),lut_ceilo_file%qback(ii) ! um, m^2,m^2
      enddo
     
    close(id)
     
    lut_ceilo_file%qback = real(lut_ceilo_file%qback*(4.d0*pi)) ! m^2
    ! 
  return
  end subroutine preload_lut_ceilo
  ! 
  subroutine GetCloudLidarMeasurements_preload(nd,diam,NN,ceilo_back_true,ceilo_ext)
  Use crsim_mod
  Use phys_param_mod, ONLY: pi
  Implicit None
  
  integer,Intent(in)         :: nd
  real*8,Intent(in)          :: diam(nd) ! m
  real*8,Intent(in)          :: NN(nd)   ! 1/m^3
  real*8,Intent(Out)         :: ceilo_back_true ! [m sr]^-1
  real*8,Intent(Out)         :: ceilo_ext ! [m]^-1 
  !
  real*8,Dimension(:),Allocatable   :: qext_new, qback_new !  m^2,m^2
  ! see ftp://ftp.astro.princeton.edu/draine/scat/bhmie/callbhmie.f 
    !
    Allocate(qext_new(nd),qback_new(nd))
    qext_new=0.d0 ; qback_new=0.d0
    call linpol(dble(lut_ceilo_file%ldiam),dble(lut_ceilo_file%qext),lut_ceilo_file%ndiam,diam*1.d+6,qext_new,nd)
    call linpol(dble(lut_ceilo_file%ldiam),dble(lut_ceilo_file%qback),lut_ceilo_file%ndiam,diam*1.d+6,qback_new,nd)
    !
    ceilo_back_true=(SUM(NN*qback_new))/(4.d0*pi) ! m-1 sr-1
    ceilo_ext=(SUM(NN*qext_new)) ! m-1
    !
    Deallocate(qext_new,qback_new)
    ! 
  return
  end subroutine GetCloudLidarMeasurements_preload
  !  
#endif
  !   
  subroutine GetCloudLidarMeasurements(conf,nd,diam,NN,ceilo_back_true,ceilo_ext)
  Use crsim_mod
  Use phys_param_mod, ONLY: pi
  Implicit None
  Type(conf_var),Intent(in)  :: conf
  integer,Intent(in)         :: nd
  real*8,Intent(in)          :: diam(nd) ! m
  real*8,Intent(in)          :: NN(nd)   ! 1/m^3
  real*8,Intent(Out)         :: ceilo_back_true ! [m sr]^-1
  real*8,Intent(Out)         :: ceilo_ext ! [m]^-1 
  !
  character                  :: LutFileName*340
  integer, parameter         :: id=1
  real                       :: tempK,wavelength !  K and um
  real                       :: ni_r,ni_c ! real and imaginary part of refr. index
  integer                    :: ndiam ! number of diams stored
  !
  real,  Dimension(:),Allocatable   :: ldiam, qext, qback ! um, m^2,m^2
  real*8,Dimension(:),Allocatable   :: qext_new, qback_new !  m^2,m^2
  integer                           :: ii
    !
    !we construct the LUT lidar filename
    LutFileName='ceilo/cld_ceilo_905nm_p25.dat'
    !
    LutFileName=Trim(Adjustl(conf%path_back_scatt_libs))//&
                Trim(Adjustl(LutFileName))
    !
    open(unit=id,file=LutFileName)
      read(id,*) tempK,wavelength  ! K and um
      read(id,*) ni_r,ni_c
      read(id,*) ndiam
      !
      Allocate(ldiam(ndiam),qext(ndiam),qback(ndiam))
      !
      do ii=1,ndiam
        read(id,*) ldiam(ii),qext(ii),qback(ii) ! um, m^2,m^2
      enddo
      !
    close(id)
    !
    qback=qback*real(4.d0*pi) ! m^2
    ! see ftp://ftp.astro.princeton.edu/draine/scat/bhmie/callbhmie.f 
    !
    Allocate(qext_new(nd),qback_new(nd))
    qext_new=0.d0 ; qback_new=0.d0
    call linpol(dble(ldiam),dble(qext),ndiam,diam*1.d+6,qext_new,nd)
    call linpol(dble(ldiam),dble(qback),ndiam,diam*1.d+6,qback_new,nd)
    Deallocate(ldiam,qext,qback)
    !
    ceilo_back_true=(SUM(NN*qback_new))/(4.d0*pi) ! m-1 sr-1
    ceilo_ext=(SUM(NN*qext_new)) ! m-1
    !
    Deallocate(qext_new,qback_new)
    !
  return
  end subroutine GetCloudLidarMeasurements
  !    
  subroutine GetCloud_MPL_Measurements_sbm(isc,conf,nd,diam,rho,NN,mpl_back_true,mpl_ext)
  Use crsim_mod
  Implicit None
  Integer,Intent(In)         :: isc
  Type(conf_var),Intent(in)  :: conf
  integer,Intent(in)         :: nd
  real*8,Intent(in)          :: diam(nd) ! m
  real*8,Intent(in)          :: rho(nd)  ! kg/m^3  hydrometeor density 
  real*8,Intent(in)          :: NN(nd)   ! 1/m^3
  real*8,Intent(Out)         :: mpl_back_true ! [m sr]^-1
  real*8,Intent(Out)         :: mpl_ext ! [m]^-1 
  !
  character                  :: LutFileName*340
  character                  :: hydroID*3,fr_string*3,den_string*3,temp_string*6
  !
  integer,Parameter           :: nden_ice=10
  Real*8                      :: den_ice_mpl(nden_ice)
  !
  Integer                     :: ii,j,irho
  Real*8                      :: var,w1,w2
  Real*8,Dimension(:),Allocatable    :: rho_new
  Real*8,Dimension(:,:),Allocatable  :: lrho
  Integer                     :: ii_start,ii_end
  Integer                     :: n1,n2,nd_new,nd_new_sum
  Real*8                      :: rho_ref, mpl_back_true_temp,mpl_ext_temp
  Real*8,Dimension(:),Allocatable    :: diam_new,NN_new
  !! luts
  data den_ice_mpl/100.d0,200.d0,300.d0,400.d0,500.d0,600.d0,700.d0,800.d0,900.d0,917.d0/
    !  
    !we construct the LUT lidar filename
    !
    fr_string='353'
    if (conf%mplID==2) fr_string='532'
    !!
    If (isc==1) Then ! CLOUD
      !-----------------------------------------------------------------------------------------
      !-----------------------------------------------------------------------------------------
      hydroID='cld'
      den_string='000'
      temp_string='_p20_d'
      !
      LutFileName='mpl/'//Trim(hydroID)//'_lidar_'//&
                          Trim(fr_string)//Trim(temp_string)//Trim(den_string)//'.dat'
      !
      LutFileName=Trim(Adjustl(conf%path_back_scatt_libs))//&
                  Trim(Adjustl(LutFileName))
      !
      !call read_mpl_luts(Trim(LutFileName),nd,diam,NN,mpl_back_true,mpl_ext)
#ifdef __PRELOAD_LUT__
      call load_mpl_luts(isc,nd,diam,NN,mpl_back_true,mpl_ext)
#else
      call read_mpl_luts(Trim(LutFileName),nd,diam,NN,mpl_back_true,mpl_ext)
#endif
      !
      !-----------------------------------------------------------------------------------------
      !-----------------------------------------------------------------------------------------
    Else If (isc==3) Then ! ICE
      !-----------------------------------------------------------------------------------------
      !-----------------------------------------------------------------------------------------
      hydroID='ice'
      temp_string='_m30_d'
      !
      !--------------------------------------------------------------------
      ! den_string
      !
      Allocate(rho_new(nd)) ; rho_new=0.d0
      Allocate(lrho(nden_ice,2)) ; lrho=0.d0
      !
      !---------------------------------------------------------------------
      !create limits for LUTs densities
      do j=1,nden_ice
        If (j==1) Then
          w1=den_ice_mpl(1)
          w2= den_ice_mpl(1) + 0.5d0*(den_ice_mpl(2)-den_ice_mpl(1))
        ElseIf( j==nden_ice) Then
          w1=den_ice_mpl(j) - 0.5d0*(den_ice_mpl(j)-den_ice_mpl(j-1))
          w2=den_ice_mpl(j)
        Else
          w1=den_ice_mpl(j) - 0.5d0*(den_ice_mpl(j)-den_ice_mpl(j-1))
          w2=den_ice_mpl(j) + 0.5d0*(den_ice_mpl(j+1)-den_ice_mpl(j))
        EndIf
        
        lrho(j,1)=w1
        lrho(j,2)=w2
      enddo
      !---------------------------------------------------------------------
      ! create array rho_new with closest densities from luts instead of actual densities
      do ii=1,nd
        var=rho(ii)
        irho=-1
        if ( var<=lrho(1,2)) irho=1
        if ( var>lrho(nden_ice,1)) irho=nden_ice
        !
        if (irho>0) go to 111
        A1:do j=2,nden_ice-1
          if ( (var>lrho(j,1)) .and. (var<=lrho(j,2)) ) then
            irho=j
            exit A1
          endif
        enddo A1!! j
        111 continue
        rho_new(ii)=den_ice_mpl(irho)
      enddo ! ii
      !---------------------------------------------------------------------     
      !
      ii_start=1; ii_end=nd
      rho_ref=rho_new(ii_start)
      n1=-1 ; n2=-1
      
      mpl_back_true_temp=0.d0 ; mpl_ext_temp=0.d0
      mpl_back_true=0.d0 ; mpl_ext=0.d0
      !
      nd_new_sum=0
      
      112 continue
      B1:do ii=ii_start,ii_end
        if (rho_new(ii)/=rho_ref) then
          n1=ii_start ; n2=max(ii_start,ii-1)
          ii_start=ii
          exit B1
        else  
          if (ii_start==ii_end) then
            n1=ii ; n2=ii
          else
            n1=ii_start ; n2=ii_end
          endif      
        endif
      enddo B1
      !------------------------------------------------------------------------------------------------------------
      nd_new=n2-n1+1
      allocate(diam_new(nd_new),NN_new(nd_new))
      diam_new(1:nd_new)=diam(n1:n2)
      NN_new(1:nd_new)=NN(n1:n2)
   
      write(den_string,'(i3)') int(rho_ref) ! K
     
      LutFileName='mpl/'//Trim(hydroID)//'_lidar_'//Trim(fr_string)//Trim(temp_string)//Trim(den_string)//'.dat'
      LutFileName=Trim(Adjustl(conf%path_back_scatt_libs))//&
                  Trim(Adjustl(LutFileName))
         
      !call read_mpl_luts(Trim(LutFileName),nd_new,diam_new,NN_new,mpl_back_true_temp,mpl_ext_temp)
#ifdef __PRELOAD_LUT__
      call load_mpl_luts(isc,nd_new,diam_new,NN_new,mpl_back_true_temp,mpl_ext_temp)
#else
      call read_mpl_luts(Trim(LutFileName),nd_new,diam_new,NN_new,mpl_back_true_temp,mpl_ext_temp)
#endif
      mpl_back_true=mpl_back_true+mpl_back_true_temp
      mpl_ext=mpl_ext+mpl_ext_temp
      deallocate(diam_new,NN_new)
      !------------------------------------------------------------------------------------------------------------
      rho_ref=rho_new((ii_start))
      nd_new_sum=nd_new_sum+nd_new
     
      if (nd_new_sum==nd) go to 113
      if (ii_start<=ii_end)  go to 112
     
      113 continue
      Deallocate(rho_new)
      Deallocate(lrho)
      !-----------------------------------------------------------------------------------------
      !-----------------------------------------------------------------------------------------
    else
      write(*,*) 'Problem: subr.GetCloud_MPL_Measurements_sbm'
      stop
    endif
    !
  return
  end subroutine GetCloud_MPL_Measurements_sbm
  !
#ifdef __PRELOAD_LUT__
  !   
  subroutine make_mpl_lutfilename(isc, conf, LutFileName)
  use crsim_mod
  use crsim_luts_mod
  Implicit None
  
  Integer, intent(In)        :: isc
  Type(conf_var),Intent(in)  :: conf
  
  Character(len=340),Intent(out) :: LutFileName
  
  character                  :: hydroID*3,fr_string*3,den_string*3,temp_string*6
    ! 
    fr_string='353'
    if (conf%mplID==2) fr_string='532'
    !
    temp_string='_m30_d'
    !
    if ((isc==1) .or. (isc==7)) then
      hydroID='cld'
      den_string='000'
      temp_string='_p20_d'
    endif
    !
    if (isc==3) then
      hydroID='ice'
      ! if (conf%MP_PHYSICS==9)  den_string='840'
      if (conf%MP_PHYSICS==9)  den_string='800' !changed by oue because mpl lut does not includes 840 2017/12/08
      if (conf%MP_PHYSICS==10) den_string='500'
      if (conf%MP_PHYSICS==8)  den_string='900' !Added by oue 2016/09/16: Note that Thopmson scheme assumes ice density of 890
      if (conf%MP_PHYSICS==30)  den_string='400' !Added by oue 2017/07/17:
      if (conf%MP_PHYSICS==40)  den_string='400' !Added by oue 2017/07/21:
      if (conf%MP_PHYSICS==50)  den_string='900' !Added by DW for P3 
      if (conf%MP_PHYSICS==75) den_string='500' !Added by oue for SAM morr 2018/06
      if (conf%MP_PHYSICS==80) den_string='500' !Added by oue for CM1 morr Apr 2020
    endif
    !
    LutFileName='mpl/'//Trim(hydroID)//'_lidar_'//Trim(fr_string)//Trim(temp_string)//Trim(den_string)//'.dat'
    !
    LutFileName=Trim(Adjustl(conf%path_back_scatt_libs))//&
                Trim(Adjustl(LutFileName))
    ! 
  return
  end subroutine make_mpl_lutfilename
  !   
  subroutine preload_mpl_luts(LutFileName, mpl_lut_file)
  use crsim_mod
  Use phys_param_mod, ONLY: pi
  Implicit None
  
  Character(len=*),Intent(In)       :: LutFileName
  Type(mpl_lut_file_type),Intent(InOut)   :: mpl_lut_file
  
  integer,parameter                  :: id=1
  integer :: ii
  
    open(unit=id,file=Trim(LutFileName))
      read(id,*) mpl_lut_file%tempK, mpl_lut_file%density, mpl_lut_file%wavelength  ! K,gr/cm^3  and um
      read(id,*) mpl_lut_file%ni_r, mpl_lut_file%ni_c
      read(id,*) mpl_lut_file%ndiam
      !
      Allocate(mpl_lut_file%ldiam(mpl_lut_file%ndiam),mpl_lut_file%qext(mpl_lut_file%ndiam),mpl_lut_file%qback(mpl_lut_file%ndiam))
      !
      do ii=1,mpl_lut_file%ndiam
        read(id,*) mpl_lut_file%ldiam(ii),mpl_lut_file%qext(ii),mpl_lut_file%qback(ii) ! um, m^2,m^2
      enddo
      !
    close(id)

    mpl_lut_file%qback = mpl_lut_file%qback*real(4.d0*pi) ! m^2
     
  return
  end subroutine preload_mpl_luts
  !   
   
  subroutine load_mpl_luts(isc,nd,diam,NN,mpl_back_true,mpl_ext)
  use crsim_mod
  use phys_param_mod, ONLY: pi
  Implicit None
  
  Integer,Intent(In)           :: isc
  Integer,Intent(In)           :: nd
  Real*8,Intent(In)            :: diam(nd) ! m
  real*8,Intent(in)            :: NN(nd)   ! 1/m^3
  Real*8,Intent(Out)           :: mpl_back_true ! m-1 sr-1
  Real*8,Intent(Out)           :: mpl_ext ! m-1
  !
  real*8,Dimension(:),Allocatable   :: qext_new, qback_new !  m^2,m^2
  !!--------------------------------------------------------------
  ! see ftp://ftp.astro.princeton.edu/draine/scat/bhmie/callbhmie.f 
  !!--------------------------------------------------------------
    !
    Allocate(qext_new(nd),qback_new(nd))
    qext_new=0.d0 ; qback_new=0.d0
    call linpol(dble(mpl_lut_files(isc)%ldiam),dble(mpl_lut_files(isc)%qext),mpl_lut_files(isc)%ndiam,diam*1.d+6,qext_new,nd)
    call linpol(dble(mpl_lut_files(isc)%ldiam),dble(mpl_lut_files(isc)%qback),mpl_lut_files(isc)%ndiam,diam*1.d+6,qback_new,nd)
    !
    mpl_back_true=(SUM(NN*qback_new))/(4.d0*pi) ! m-1 sr-1
    mpl_ext=(SUM(NN*qext_new))! m-1 
    !
    Deallocate(qext_new,qback_new)!
    ! 
  return
  end subroutine load_mpl_luts
  !   
#endif
  !     
  subroutine read_mpl_luts(LutFileName,nd,diam,NN,mpl_back_true,mpl_ext)
  use crsim_mod
  use phys_param_mod, ONLY: pi
  Implicit None
  Character(len=*),Intent(In)  :: LutFileName
  Integer,Intent(In)           :: nd
  Real*8,Intent(In)            :: diam(nd) ! m
  real*8,Intent(in)            :: NN(nd)   ! 1/m^3
  Real*8,Intent(Out)           :: mpl_back_true ! m-1 sr-1
  Real*8,Intent(Out)           :: mpl_ext ! m-1
  !
  integer, parameter         :: id=1
  real                       :: tempK,density,wavelength !  K,gr/cm^3 and um
  real                       :: ni_r,ni_c ! real and imaginary part of refr. index
  integer                    :: ndiam ! number of diams stored
  !
  real,  Dimension(:),Allocatable   :: ldiam, qext, qback ! um, m^2,m^2
  real*8,Dimension(:),Allocatable   :: qext_new, qback_new !  m^2,m^2
  integer                           :: ii
    !
    !--------------------------------------------------------------
    open(unit=id,file=Trim(LutFileName))
      read(id,*) tempK,density,wavelength  ! K,gr/cm^3  and um
      read(id,*) ni_r,ni_c
      read(id,*) ndiam
      !
      Allocate(ldiam(ndiam),qext(ndiam),qback(ndiam))
      !
      do ii=1,ndiam
        read(id,*) ldiam(ii),qext(ii),qback(ii) ! um, m^2,m^2
      enddo
      !
    close(id)
    !!--------------------------------------------------------------
    ! see ftp://ftp.astro.princeton.edu/draine/scat/bhmie/callbhmie.f 
    qback=qback*real(4.d0*pi) ! m^2
    !!
    Allocate(qext_new(nd),qback_new(nd))
    qext_new=0.d0 ; qback_new=0.d0
    call linpol(dble(ldiam),dble(qext),ndiam,diam*1.d+6,qext_new,nd)
    call linpol(dble(ldiam),dble(qback),ndiam,diam*1.d+6,qback_new,nd)
    Deallocate(ldiam,qext,qback)
    !
    mpl_back_true=(SUM(NN*qback_new))/(4.d0*pi) ! m-1 sr-1
    mpl_ext=(SUM(NN*qext_new))! m-1 
    !
    Deallocate(qext_new,qback_new)!
    !!
  return
  end subroutine read_mpl_luts
  !
  subroutine GetCloud_MPL_Measurements(isc,conf,nd,diam,NN,mpl_back_true,mpl_ext)
  Use crsim_mod
  Implicit None
  Integer,Intent(In)         :: isc
  Type(conf_var),Intent(in)  :: conf
  integer,Intent(in)         :: nd
  real*8,Intent(in)          :: diam(nd) ! m
  real*8,Intent(in)          :: NN(nd)   ! 1/m^3
  real*8,Intent(Out)         :: mpl_back_true ! [m sr]^-1
  real*8,Intent(Out)         :: mpl_ext ! [m]^-1 
  !
!#ifdef __PRELOAD_LUT__    ! AT April 2019 commented line
  character                  :: LutFileName*340
  character                  :: hydroID*3,fr_string*3,den_string*3,temp_string*6
    !
    LutFileName=''  ! AT April 2019  
    !
    !we construct the LUT lidar filename
    !
    fr_string='353'
    if (conf%mplID==2) fr_string='532'
    !
    temp_string='_m30_d'
    !
    if ((isc==1).or.(isc==7)) then
      hydroID='cld'
      den_string='000'
      temp_string='_p20_d'
    endif
    !
    if (isc==3) then
      hydroID='ice'
      ! if (conf%MP_PHYSICS==9)  den_string='840'
      if (conf%MP_PHYSICS==9)  den_string='800' !changed by oue because mpl lut does not includes 840 2017/12/08
      if (conf%MP_PHYSICS==10) den_string='500'
      if (conf%MP_PHYSICS==8)  den_string='900' !Added by oue 2016/09/16: Note that Thopmson scheme assumes ice density of 890
      if (conf%MP_PHYSICS==30)  den_string='400' !Added by oue 2017/07/17
      if (conf%MP_PHYSICS==40)  den_string='400' !Added by oue 2017/07/21
      if (conf%MP_PHYSICS==50)  den_string='900' !Added by DW for P3 
      if (conf%MP_PHYSICS==75) den_string='500' !Added by oue 2018/06 for SAM morr
      if (conf%MP_PHYSICS==80) den_string='500' !Added by oue for CM1 morr Apr 2020
    endif
    !
    ! AT April 2019 : the definition of LutFileName moved below (inside of ifdef  --PRELOAD_LUT__ loop)
    ! ifdef __PRELOAD_LUT__
!    !write(*,*) Trim(den_string)
!    LutFileName='mpl/'//Trim(hydroID)//'_lidar_'//Trim(fr_string)//Trim(temp_string)//Trim(den_string)//'.dat'
!    !write(*,*) Trim(LutFileName)
!    !
!    LutFileName=Trim(Adjustl(conf%path_back_scatt_libs))//&
!                Trim(Adjustl(LutFileName))
    !
!#endif   ! AT April 2019 commented loop
 
#ifdef __PRELOAD_LUT__
    call load_mpl_luts(isc,nd,diam,NN,mpl_back_true,mpl_ext)
#else
    LutFileName='mpl/'//Trim(hydroID)//'_lidar_'//Trim(fr_string)//Trim(temp_string)//Trim(den_string)//'.dat'
    LutFileName=Trim(Adjustl(conf%path_back_scatt_libs))//&
                Trim(Adjustl(LutFileName))
    call read_mpl_luts(Trim(LutFileName),nd,diam,NN,mpl_back_true,mpl_ext)
#endif
    !
  return
  end subroutine GetCloud_MPL_Measurements
  !!
  subroutine compute_molecular_backscatter(conf,env,mpl)
  Use crsim_mod
  Implicit None
  !!
  Type(conf_var),Intent(in)    :: conf
  Type(env_var),Intent(in)     :: env
  Type(mpl_var),Intent(InOut)  :: mpl
  !
  real*8             :: lambda_fac
  real*8,Dimension(:,:,:),Allocatable    :: tp_fac
    !
    Allocate(tp_fac(mpl%nx,mpl%ny,mpl%nz))
    tp_fac=0.d0
     
    lambda_fac=mpl%wavel*1.d-9 ! wavelength (nm->m)
    lambda_fac=lambda_fac**4.0117
    !
    tp_fac=(env%press*1.d+2)/(env%temp+273.15d0)  ! Pa/K
    !
    mpl%rayleigh_back=(1.d-2*tp_fac)  * ( 1.d0/lambda_fac ) 
    mpl%rayleigh_back=mpl%rayleigh_back* 2.938d0 * 1.d-32  ! m^-1 sr^-1
    !
    ! 8 pi/3 is backstering for pure Rayleigh scattering
    !mpl%rayleigh_back=mpl%rayleigh_back*(8.d0*pi/3.d0)  ! m^-1 sr^-1 !commented out by oue
    !
    Deallocate(tp_fac)
    !
  return
  end subroutine compute_molecular_backscatter
  !!
  !!
  subroutine compute_polarim_vars(elev,horientID,sigma,nd,lam_mm,fa,fb,fa0,fb0,NN,fvel,ww,&
                                  Zhh,Zvv,Zvh,RHOhvc,DVh,Dopp,dDVh,Kdp,Ah,Av,Adp,diff_back_phase,zhh_d,zvh_d,zvv_d)
  use crsim_mod
  use phys_param_mod, ONLY: pi,d2r, K2
  Implicit None
  
  real*8,Intent(In)         :: elev    ! degree
  integer,Intent(In)        :: horientID ! choice for the orientation distribution
  real*8,Intent(In)         :: sigma ! sigma for horientID=3
  integer,Intent(In)        :: nd
  real*8,Intent(In)         :: lam_mm ! radar wavelength mm
  complex*16,Intent(In)     :: fa(nd),fb(nd)   ! backward amplitudes mm
  complex*16,Intent(In)     :: fa0(nd),fb0(nd) ! forward amplitudes mm
  real*8,Intent(In)         :: NN(nd)   ! 1/m^3
  real*8,Intent(In)         :: fvel(nd) ! m/s
  real*8,Intent(In)         :: ww ! m/s air vertical velocity
  real*8,Intent(Out)        :: Zhh,Zvv,Zvh,RHOhvc,DVh,Dopp,dDVh,Kdp,Ah,Av,Adp,diff_back_phase
  real*8,Intent(Out)        :: zhh_d(nd),zvh_d(nd),zvv_d(nd) ! used for spectrum generation
  !
  real*8,Dimension(:),Allocatable      :: zhh1,zvv1,zvh1,kdp1,sigma_ext_h,sigma_ext_v,rhohv1c
  real*8,Dimension(:),Allocatable      :: w1,w2
  real*8                    :: Cpol_mm ! mm^4
  !
  real*8                    :: aa1,aa2,aa3,aa4,aa5,aa6,aa7
    !
  !
  real*8                    :: abs_elev ! absolute value of elev for negative elevation by oue 20190911
  !
  !---------------------------------
  !---------------------------------
  ! changed elev to abs(elev) for negative elevation by oue 20190911
  abs_elev = abs(elev)
  !---------------------------------
  !---------------------------------
    Cpol_mm=(lam_mm*lam_mm*lam_mm*lam_mm)/(pi*pi*pi*pi*pi * K2 ) ! mm^4

    if((elev<=90.d0) .and. (elev>=0.d0)) then
    call get_orientation_coeff(horientID,elev,sigma,aa1,aa2,aa3,aa4,aa5,aa6,aa7)
    end if
    if((elev>=-90.d0) .and. (elev<0.d0)) then
       call get_orientation_coeff(horientID,abs_elev,sigma,aa1,aa2,aa3,aa4,aa5,aa6,aa7) ! changed elev to abs(elev) for negative elevation by oue 20190911
    endif
    Allocate(zhh1(nd),zvv1(nd),zvh1(nd),kdp1(nd),sigma_ext_h(nd),sigma_ext_v(nd),rhohv1c(nd))
    !
    zhh1=0.d0 ; zvv1=0.d0; zvh1=0.d0; kdp1=0.d0; sigma_ext_h=0.d0 
    sigma_ext_v=0.d0
    !
    !-----------------------------------------------------------------------------------------------------------
    zhh1=dble(Cpol_mm*4.d0*pi*(fb*dconjg(fb)-2.d0*dble(real(dconjg(fb)*(fb-fa)*aa2))+((fb-fa)*dconjg(fb-fa))*aa4)) ! mm^6
    zvv1=dble(Cpol_mm*4.d0*pi*(fb*dconjg(fb)-2.d0*dble(real(dconjg(fb)*(fb-fa)*aa1))+((fb-fa)*dconjg(fb-fa))*aa3)) ! mm^6
    zvh1=dble(Cpol_mm*4.d0*pi*(((fb-fa)*dconjg(fb-fa))*aa5)) ! mm^6
    rhohv1c=dble(Cpol_mm*4.d0*pi*(  fb*dconjg(fb) + (fb-fa)*dconjg(fb-fa)*aa5 - dconjg(fb)*(fb-fa)*aa1 &
                                                 - fb*(dconjg(fb)-dconjg(fa))*aa2 ) )  ! mm^6
    !
    kdp1=lam_mm*dble(real(fb0-fa0))*aa7 ! rad/m
    !sigma_ext_h=2.d0*lam_mm*dimag(fb0) ! mm^2
    !sigma_ext_v=2.d0*lam_mm*dimag(fa0) ! mm^2
    ! AT Aug 19, 2019
    sigma_ext_h=2.d0*lam_mm * ( dimag(fb0) - aa2*dimag(fb0-fa0) ) ! mm^2
    sigma_ext_v=2.d0*lam_mm * ( dimag(fb0) - aa1*dimag(fb0-fa0) ) ! mm^2
    !
    
    Allocate(w1(nd),w2(nd))
    w1=0.d0 ; w2=0.d0
    !
    !w1=lam_mm*NN*1.d-9*dble(fb*dconjg(fa)) 
    !w2=lam_mm*NN*1.d-9*dimag(fb*dconjg(fa))
     
    w1=lam_mm*NN*1.d-9*dble(fb*dconjg(fb) + ((fb-fa)*dconjg(fb-fa))*aa3  &
       -(dconjg(fb)*(fb-fa))*aa1 - (fb*(dconjg(fb)-dconjg(fa)))*aa2)
    w2=lam_mm*NN*1.d-9*dimag(fb*dconjg(fb) + ((fb-fa)*dconjg(fb-fa))*aa3 &
       - (dconjg(fb)*(fb-fa))*aa1- (fb*(dconjg(fb)-dconjg(fa)))*aa2)
     
    diff_back_phase=datan2(Sum(w2),Sum(w1))
    diff_back_phase=diff_back_phase*180.d0/pi ! deg
    !
    Deallocate(w1,w2)
    !-----------------------------------------------------------------------------------------------------------
    !
    !do i=1, nd
    !  write(*,*) 'zhh1 * NN (',i,')', zhh1(i)*NN(i), ', NN : ', NN(i), ', zhh1 : ', zhh1(i), ', fb : ', fb(i), ', fa : ', fa(i)
    !enddo
    ! 
    Zhh=Sum(zhh1*NN)  ! mm^6/m^3
    Zvv=Sum(zvv1*NN)  ! mm^6/m^3 
    Zvh=Sum(zvh1*NN)  ! mm^6/m^3 
    RHOhvc=Sum(rhohv1c*NN)! mm^6/m^3 

    !
    DVh=Sum(zhh1*NN*fvel)   ! mm^6/m^3 * m/s
    ! AT use conf%airborne for Dopp here
    if ((elev>=-90.d0) .and. (elev < 0.d0)) DVh=Sum(zhh1*NN*fvel*(-1.d0)) ! modified for negative elevation by oue 20190911
    Dopp=(ww*Zhh-DVh) !  mm^6/m^3 * m/s, positive upward or away from the radar
    !
    if (Zhh>0.d0) then  ! modified for negative elevation by oue 20190911
       dDVh=Sum(zhh1*NN*(fvel-DVh/Zhh)*(fvel-DVh/Zhh)) ! mm^6/m^3*(m/s)^2
       if (elev < 0.d0) dDVh=Sum(zhh1*NN*(fvel*(-1.d0)-DVh/Zhh)*(fvel*(-1.d0)-DVh/Zhh))
    end if
    !if (Zhh>0.d0) dDVh=Sum(zhh1*NN*(fvel-DVh/Zhh)*(fvel-DVh/Zhh)) ! mm^6/m^3*(m/s)^2
    !
    !
    Kdp=180.d0/pi*Sum(kdp1*NN)*1.d-3  ! deg/km  specific differencial phase
    Ah=(10.d0/dlog(10.d0))*1.d-6*Sum(sigma_ext_h*NN) * 1.d+3  ! dB/km  Specific horizontal attenuation
    Av=(10.d0/dlog(10.d0))*1.d-6*Sum(sigma_ext_v*NN) * 1.d+3  ! dB/km  Specific vertical attenuation
    Adp=(10.d0/dlog(10.d0))*1.d-6*Sum((sigma_ext_h-sigma_ext_v)*NN) * 1.d+3  ! dB/km differencial attenuation
    !
    zhh_d = zhh1; zvh_d = zvh1; zvv_d = zvv1
    !
    Deallocate(zhh1,zvv1,zvh1,kdp1,sigma_ext_h,sigma_ext_v,rhohv1c)
    !
  return
  end subroutine compute_polarim_vars
  !
  !
  subroutine compute_rad_forward_vars(elev,horientID,sigma,nd,lam_mm,fa,fb,fa0,fb0,NN,fvel,ww,Zhh,DVh,Dopp,dDVh,Ah,zhh_d)
  use crsim_mod
  use phys_param_mod, ONLY: pi, d2r, K2
  Implicit None
  
  real*8,Intent(In)         :: elev    ! degree
  integer,Intent(In)        :: horientID ! choice of orientation distribution  
  real*8,Intent(In)         :: sigma   ! sigma for horientID=3
  integer,Intent(In)        :: nd
  real*8,Intent(In)         :: lam_mm ! radar wavelength mm
  complex*16,Intent(In)     :: fa(nd),fb(nd)   ! backward amplitudes mm
  complex*16,Intent(In)     :: fa0(nd),fb0(nd) ! forward amplitudes mm
  real*8,Intent(In)         :: NN(nd)   ! 1/m^3
  real*8,Intent(In)         :: fvel(nd) ! m/s
  real*8,Intent(In)         :: ww ! m/s air vertical velocity
  real*8,Intent(Out)        :: Zhh,DVh,Dopp,dDVh,Ah
  real*8,Intent(Out)        :: zhh_d(nd) ! used for spectrum generation
  !!
  real*8,Dimension(:),Allocatable  :: zhh1,sigma_ext_h
  real*8                    :: Cpol_mm ! mm^4
  !
  real*8                    :: aa2,aa4
  !
  real*8                    :: abs_elev ! absolute value of elev for negative elevation by oue 20190911
  !
  !--------------------------------------
  !--------------------------------------
  ! changed elev to abs(elev) for negative elevation by oue 20190911
  abs_elev = abs(elev)
  !--------------------------------------
  !--------------------------------------
  !
    Cpol_mm=(lam_mm*lam_mm*lam_mm*lam_mm)/(pi*pi*pi*pi*pi * K2 ) ! mm^4
    !
    aa2=0.d0; aa4=0.d0;
    if((elev<=90.d0) .and. (elev>=0.d0)) then
       call get_orientation_coeff_24(horientID,elev,sigma,aa2,aa4)
    end if
    if((elev>=-90.d0) .and. (elev<0.d0)) then
       call get_orientation_coeff_24(horientID,abs_elev,sigma,aa2,aa4)! changed elev to abs(elev) for negative elevation by oue 20190911
    end if
    !
    Allocate(zhh1(nd),sigma_ext_h(nd))
    zhh1=0.d0 ; sigma_ext_h=0.d0 
    !-----------------------------------------------------------------------------------------------------------
    zhh1=dble(Cpol_mm*4.d0*pi*(fb*dconjg(fb)-2.d0*dble(real(dconjg(fb)*(fb-fa)*aa2))+ ((fb-fa)*dconjg(fb-fa))*aa4)) ! mm^6
    ! AT Aug 19, 2019
    !sigma_ext_h=2.d0*lam_mm*dimag(fb0) ! mm^2
    sigma_ext_h=2.d0*lam_mm * ( dimag(fb0) - aa2*dimag(fb0-fa0) ) ! mm^2

    !-----------------------------------------------------------------------------------------------------------
    !
    Zhh=Sum(zhh1*NN)  ! mm^6/m^3
    DVh=Sum(zhh1*NN*fvel)   ! mm^6/m^3 * m/s
    ! AT use conf%airborne for Dopp
    if ((elev>=-90.d0) .and. (elev < 0.d0)) DVh=Sum(zhh1*NN*fvel*(-1.d0)) ! modified for negative elevation by oue 20190911
    Dopp=(ww*Zhh-DVh)
    if (Zhh>0.d0) then ! modified for negative elevation by oue 20190911
       dDVh=Sum(zhh1*NN*(fvel-DVh/Zhh)*(fvel-DVh/Zhh)) ! mm^6/m^3*(m/s)^2
       if ((elev>=-90.d0) .and. (elev < 0.d0)) dDVh=Sum(zhh1*NN*(fvel*(-1.d0)-DVh/Zhh)*(fvel*(-1.d0)-DVh/Zhh))
    end if
    !if (Zhh>0.d0) dDVh=Sum(zhh1*NN*(fvel-DVh/Zhh)*(fvel-DVh/Zhh)) ! mm^6/m^3*(m/s)^2
    Ah=(10.d0/dlog(10.d0))*1.d-6*Sum(sigma_ext_h*NN) * 1.d+3  ! dB/km  Specific horizontal attenuation
    !
    zhh_d = zhh1
    !
    Deallocate(zhh1,sigma_ext_h)
    !
  return
  end subroutine compute_rad_forward_vars
  !
  subroutine get_orientation_coeff(orid,theta,sigma,aa1,aa2,aa3,aa4,aa5,aa6,aa7)
  use crsim_mod
  use phys_param_mod, ONLY: pi
  Implicit None
  Integer, Intent(In)   :: orid
  real*8,Intent(In)     :: theta
  real*8,Intent(In)     :: sigma
  real*8,Intent(Out)    :: aa1,aa2,aa3,aa4,aa5,aa6,aa7
  !
  real*8  :: beta,sigmar,rr
    !
    if (orid==1) then ! fully chaotic orientation
      aa1=1.d0/3.d0
      aa2=aa1
      aa3=1.d0/5.d0
      aa4=aa3
      aa5=1.d0/15.d0
      aa6=0.d0
      aa7=0.d0
    endif
    !
    if (orid==2) then ! random orientation in the horizontal plane
      beta=theta*pi/180.d0  ! elev. angle in radians
      aa1=0.5d0*dsin(beta)*dsin(beta)
      aa2=0.5d0
      aa3=3.d0/8.d0*dsin(beta)*dsin(beta)
      aa4=3.d0/8.d0
      aa5=1.d0/8.d0*dsin(beta)*dsin(beta)
      aa6=0.d0
      aa7=-0.5*dcos(beta)*dcos(beta)
    endif
    !
    if (orid==3) then ! 2-d axisymmetric gaussian distribution 
      sigmar=sigma*pi/180.d0
      rr=dexp(-2.d0*sigmar*sigmar)
      
      aa1=0.25d0*(1.d0+rr)*(1.d0+rr)
      aa2=0.25d0*(1.d0-rr*rr)
      aa3=(3.d0/8.d0 + 0.5d0*rr + 1.d0/8.d0*rr**4)**2
      aa4=(3.d0/8.d0 - 0.5d0*rr + 1.d0/8.d0*rr**4) * (3.d0/8.d0 + 0.5d0*rr +1.d0/8.d0*rr**4)
      aa5=1.d0/8.d0*(3.d0/8.d0 + 0.5d0*rr + 1.d0/8.d0*rr**4)*(1.d0-rr**4)
      aa6=0.d0
      aa7=0.5d0*rr*(1.d0+rr)
      
    endif
    !   
  return
  end subroutine get_orientation_coeff
  !
  !
  subroutine get_orientation_coeff_24(orid,theta,sigma,aa2,aa4)
  use crsim_mod
  use phys_param_mod, ONLY: pi
  !
  Implicit None
  Integer, Intent(In)   :: orid
  real*8,Intent(In)     :: theta
  real*8,Intent(In)     :: sigma
  real*8,Intent(Out)    :: aa2,aa4
  !
  real*8  :: beta,sigmar,rr
    !
    if (orid==1) then ! fully chaotic orientation
      aa2=1.d0/3.d0
      aa4=1.d0/5.d0
    endif
    !
    if (orid==2) then ! random orientation in the horizontal plane
      beta=theta*pi/180.d0  ! elev. angle in radians
      aa2=0.5d0
      aa4=3.d0/8.d0
    endif
    !
    if (orid==3) then ! 2-d axisymmetric gaussian distribution 
      sigmar=sigma*pi/180.d0
      rr=dexp(-2.d0*sigmar*sigmar)
      aa2=0.25d0*(1.d0-rr*rr)
      aa4=(3.d0/8.d0 - 0.5d0*rr + 1.d0/8.d0*rr**4) * (3.d0/8.d0 + 0.5d0*rr+1.d0/8.d0*rr**4)
    endif
    !
  return
  end subroutine get_orientation_coeff_24
  !    
#ifdef __PRELOAD_LUT__
  !subroutine make_lutfilename(isc, conf, inp_elev, temp, LutFileName)
  subroutine make_lutfilename(isc, conf, inp_elev, temp, rho_hydro,LutFileName)!by oue 2017/07/17 for preloading LUTs of all densities
  use crsim_mod
  use crsim_luts_mod
  Implicit None
  !
  Integer, intent(In)        :: isc
  Type(conf_var),Intent(in)  :: conf
  real*8,Intent(in)          :: inp_elev ! degrees 
  real*8,Intent(in)          :: temp     !K
  real*8,Intent(in)          :: rho_hydro   ! kg/m^3 !by oue 2017/07/17 for preloading LUTs of all densities
  Character(len=340),Intent(out) :: LutFileName
  ! 
  real*8                      :: fff,ttt,elev
  character                   :: frq_str*4,t_str*5,el_str*4,rho_str*5
  integer                          :: nden
  real*8,Dimension(:),Allocatable  :: lden ! kg/m^3
  real*8                           :: den  ! kg/m^3
  real*8                           :: rho_const   ! kg/m^3 
    ! 
    ! - get the frequency
    call get_lut_str_var(n_lfreq,lfreq,conf%freq,fff)
    !---------------------------------------------
    ! get the temperature
    !-------------------------------------
    if (isc==1) then! CLOUD
      call get_lut_str_var(n_ltemp_cld,ltemp_cld,temp+273.15d0,ttt)
    end if !- if cloud
    !------------------------------------
    if (isc==2) then! RAIN
      call get_lut_str_var(n_ltemp_rain,ltemp_rain,temp+273.15d0,ttt)
    end if !- if rain
    !!------------------------------------
    if (isc>2) then ! ICE,SNOW,GRAUPEL,HAIL
      ttt=ltemp_ice(1)
    endif
    !-------------------------------------
    if (isc==7) then ! Drizzle added by oue 2017/07/21 for RAMS
      call get_lut_str_var(n_ltemp_cld,ltemp_cld,temp+273.15d0,ttt)
    end if !- if drizzle
    !
    if(conf%MP_PHYSICS==50) then !by oue PRELOAD P3
      if (isc==3) ttt=ltemp_smallice(1) ! small ice
      if (isc==4) ttt=ltemp_unrice(1) ! unrimed ice
      if (isc==5) ttt=ltemp_graupP3(1) ! graupel P3
      if (isc==6) ttt=ltemp_parice(1) ! partily-rimed ice
    endif
    ! 
    !write(*,*) 'LUTs: Using temperature', ttt
    !--------------------------------------------
    !
    ! get the elevation
    if (isc>1) then
      call get_lut_str_var(n_lelev,lelev,inp_elev,elev)
    else
      elev=lelev(n_lelev)
    endif
    if (isc==7) then !Drizzle added by oue 2017/07/21 for RAMS
      elev=lelev(n_lelev)
    endif
    ! 
    !----------------------------------------------
    call get_lut_str(fff,ttt,elev,frq_str,t_str,el_str)
    !-------------------------------------------
    !
    ! hydrometeor density
    ! select the LUTs densities for hydrometeor isc !modified by oue 2017/07/21 for RAMS to include drizzle(isc=7) and aggregates(isc=8)
    if ((isc<=2) .or. (isc==7)) then ! cloud, rain, drizzle
      nden=n_lden_cld
      Allocate(lden(nden))
      lden=lden_cld
    endif
    if(conf%MP_PHYSICS==50) then !by oue PRELOAD P3 
      if (isc==3) then ! small ice 
        nden=n_lden_smallice
        Allocate(lden(nden))
        lden=lden_smallice
      endif
      if (isc==4) then ! unrimed ice  
        nden=n_lden_unrice
        Allocate(lden(nden))
        lden=lden_unrice
      endif
      if (isc==5) then ! graupel P3 
        nden=n_lden_graupP3
        Allocate(lden(nden))
        lden=lden_graupP3
      endif
      if (isc==6) then ! partily-rimed ice
        nden=n_lden_parice
        Allocate(lden(nden))
        lden=lden_parice
      endif 
    ELSE  !by oue PRELOAD P3 
      if (isc==3) then ! ice
        nden=n_lden_ice
        Allocate(lden(nden))
        lden=lden_ice
      endif
      if ((isc==4) .or. (isc==8)) then ! snow, aggregates
        nden=n_lden_snow
        Allocate(lden(nden))
        lden=lden_snow
      endif
      if ((isc==5) .or. (isc==6)) then ! graupel,hail
        nden=n_lden_graup
        Allocate(lden(nden))
        lden=lden_graup
      endif
      !  
    ENDIF !by oue PRELOAD P3 
    ! 
    !-------------------------------------------
    rho_const = rho_hydro !added by oue for preloading LUTs of all densities
    !-- from hydro_morrison & cloud_morrison ------
    !IF (conf%MP_PHYSICS==10) THEN
    !if(isc==1) then ! cloud
    !   rho_const=997.d0
    !endif
    !if(isc==2) then ! rain
    !   rho_const=997.d0
    !endif
    !if(isc==3) then ! ice
    !   rho_const=500.d0
    !endif
    !if(isc==4) then ! snow
    !   rho_const=100.d0
    !endif
    !if(isc==5) then ! graupel
    !   rho_const=400.d0
    !endif
    !ENDIF
    !-- from hydro_icon
    !IF (conf%MP_PHYSICS==30) THEN
    !if(isc==1) then ! cloud
    !   rho_const=1000.d0
    !endif
    !if(isc==2) then ! rain
    !   rho_const=1000.d0
    !endif
    !if(isc==3) then ! ice
    !   rho_const=400.d0
    !endif
    !if(isc==4) then ! snow
    !   rho_const=100.d0
    !endif
    !if(isc==5) then ! graupel
    !   rho_const=400.d0
    !endif
    !if(isc==6) then ! graupel
    !   rho_const=900.d0
    !endif
    !rho_const = rho_hydro
    !ENDIF
    !
    !
    ! for rho_const find the closest density from LUTs
    call get_lut_str_var(nden,lden,rho_const,den)
    !for the bin density den determine the string rho_str
    call get_lut_den_str(den,rho_str)
    !
    Deallocate(lden)
    !----------------------------------------------------------------
    !
    ! now we can construct the LUT filename
    !
    LutFileName=Trim(adjustl(conf%hydro_luts(isc)))//&
               '_fr'//frq_str//'GHz_'//'t'//t_str//'_rho'//rho_str//'_el'//el_str//'.dat'
    LutFileName=Trim(Adjustl(conf%path_back_scatt_libs))//&
                Trim(Adjustl(conf%hydro_luts(isc)))//'/'//&
                Trim(Adjustl(LutFileName))
    !  
  return
  end subroutine make_lutfilename
  !   
  subroutine preload_luts(LutFileName, lut_file)
  use crsim_mod
  Implicit None
  !  
  Character(len=*),Intent(In)       :: LutFileName
  Type(lut_file_type),Intent(InOut)   :: lut_file
  integer,parameter                  :: id=1
  real*8                             :: work
  integer      :: ir,ielev
    ! 
    lut_file%filenameLUTS = Trim(LutFileName)
   
    open(unit=id,file=lut_file%filenameLUTS)
      read(id,*) lut_file%llam  ! mm wavelength
      read(id,*) lut_file%ltemp ! K temperature
      read(id,*) lut_file%lref_re, lut_file%lref_im ! real and imaginary part of the refractive index
      read(id,*) lut_file%lnd, lut_file%lnelev ! number of diameters and elevations stored
      read(id,*)
      !
      Allocate(lut_file%elev(lut_file%lnelev),lut_file%aoverb(lut_file%lnd),lut_file%ldiam(lut_file%lnd))
      Allocate(lut_file%sb(2,2,lut_file%lnd,lut_file%lnelev))
      Allocate(lut_file%sf11(lut_file%lnd,lut_file%lnelev),lut_file%sf22(lut_file%lnd,lut_file%lnelev))
      !
      do ielev=1,lut_file%lnelev
        read(id,*) lut_file%elev(ielev), work ! incident angle,- elevation =90-elev
        do ir=1,lut_file%lnd
          read(id,*) lut_file%ldiam(ir),lut_file%aoverb(ir)  ! radius [mm] and (1/aspect ratio) 
          read(id,*) lut_file%sf11(ir,ielev),lut_file%sf22(ir,ielev)
          read(id,*) lut_file%sb(1,1,ir,ielev),lut_file%sb(1,2,ir,ielev),lut_file%sb(2,1,ir,ielev),lut_file%sb(2,2,ir,ielev) ! mm
        enddo ! ir
      enddo ! ielev
    close(id)
    !  
  return
  end subroutine preload_luts
  !     
  !   
  subroutine load_luts(LutFileName,elev0,ir1,ir2,nir,out_diam,out_aoverb,out_sb11,out_sb22,out_sf11,out_sf22)
  use crsim_mod
  Implicit None
  
  Character(len=*),Intent(In)       :: LutFileName
  Real*8,Intent(In)                 :: elev0 ! desired elevation
  Integer,Intent(In)                :: ir1,ir2
  Integer,Intent(In)                :: nir
  Real*8,Intent(Out)                :: out_diam(nir) ! mm
  Real*8,Intent(Out)                :: out_aoverb(nir)
  Complex*16,Intent(Out)            :: out_sb11(nir),out_sb22(nir)  ![mm]
  Complex*16,Intent(Out)            :: out_sf11(nir),out_sf22(nir)  ![mm]
  !
  integer      :: ir,ielev,ielev0
  integer      :: ii
  real*8       :: diff,diff0
  integer :: i
    ! 
    lut_search_flag = 0
    search: do i=1,size_lut_files
      !write(*,*) 'ii : ', ii, ', size_lut_files : ', size_lut_files, ',
      !lut_search_flag', lut_search_flag,  &
      !   ', String cmp : ', (Trim(LutFileName) .EQ. Trim(lut_files(ii+1)%filenameLUTS))
      !write(*,*) 'LutFileName : ', Trim(LutFileName)
      !write(*,*) 'lut_files(ii+1)%filenameLUTS : ',
      !Trim(lut_files(ii+1)%filenameLUTS)
      !
      if( Trim(LutFileName) .EQ. Trim(lut_files(i)%filenameLUTS) ) then
        lut_search_flag = lut_search_flag + 1
        exit search
      endif
    enddo search
    ! 
    if(lut_search_flag == 0) then
      write(*,*) '-----------------------------------------'
      write(*,*) 'Searching preloaded LUT file ERROR!!!'
      write(*,*) 'not preloaded file',Trim(LutFileName)
      write(*,*) 'all loaded files from 1 to', size_lut_files
      do i=1,size_lut_files
        write(*,*) 'file',i, Trim(lut_files(i)%filenameLUTS)
      enddo
      write(*,*)
      write(*,*) '-----------------------------------------'
      stop
      return
    endif
    ! 
    !
    ! find elevation elev0
    !
    if (lut_files(i)%lnelev==1) ielev0=1
    ! 
    if (lut_files(i)%lnelev>1) then
      ielev0=-1
      diff0=1.d+11
      do ielev=1,lut_files(i)%lnelev
        diff=dabs(lut_files(i)%elev(ielev)-(90.d0-elev0))
        if (diff<diff0) then
          diff0=diff
          ielev0=ielev
        endif
      enddo
    endif
    !
    !------------------------------
    ! save lut info in desired range
    ii=0
    do ir=ir1,ir2
      ii=ii+1
      out_diam(ii)=2.d0*lut_files(i)%ldiam(ir)
      out_aoverb(ii)=lut_files(i)%aoverb(ir)
      out_sf11(ii)=lut_files(i)%sf11(ir,ielev0)
      out_sf22(ii)=lut_files(i)%sf22(ir,ielev0)
      out_sb11(ii)=lut_files(i)%sb(1,1,ir,ielev0)
      out_sb22(ii)=lut_files(i)%sb(2,2,ir,ielev0)
    enddo
     
    !do ii=1, nir
    !  write(*,*) 'out_sb22 : ', out_sb22(ii), ', out_sb11 : ', out_sb11(ii)
    !enddo
    !  
  return
  end subroutine load_luts
  ! 
#endif
  ! 
  subroutine read_luts(LutFileName,elev0,ir1,ir2,nir,out_diam,out_aoverb,out_sb11,out_sb22,out_sf11,out_sf22)
  use crsim_mod
  Implicit None
  !
  Character(len=*),Intent(In)       :: LutFileName
  Real*8,Intent(In)                 :: elev0 ! desired elevation
  Integer,Intent(In)                :: ir1,ir2
  Integer,Intent(In)                :: nir
  Real*8,Intent(Out)                :: out_diam(nir) ! mm
  Real*8,Intent(Out)                :: out_aoverb(nir)
  Complex*16,Intent(Out)            :: out_sb11(nir),out_sb22(nir)  ![mm]
  Complex*16,Intent(Out)            :: out_sf11(nir),out_sf22(nir)  ![mm]
  !
  integer,parameter                  :: id=1
  real*8                             :: work
  real*8,Dimension(:),Allocatable              :: elev
  real*8,Dimension(:),Allocatable              :: ldiam ! mm
  real*8,Dimension(:),Allocatable              :: aoverb
  complex*16,Dimension(:,:,:,:),Allocatable    :: sb ! mm
  complex*16,Dimension(:,:),Allocatable        :: sf11,sf22 ! mm
  !
  real*8                      :: llam! mm
  real*8                      :: ltemp ! K
  real*8                      :: lref_re,lref_im
  integer                     :: lnd,lnelev
  integer      :: ir,ielev,ielev0
  integer      :: ii
  real*8       :: diff,diff0
    !
    !write(*,*) 'read_luts() : LutFileName : ', LutFileName
    ! 
    open(unit=id,file=Trim(LutFileName))
      read(id,*) llam  ! mm wavelength
      read(id,*) ltemp ! K temperature
      read(id,*) lref_re, lref_im ! real and imaginary part of the refractive index
      read(id,*) lnd, lnelev ! number of diameters and elevations stored
      read(id,*)
      ! 
      Allocate(elev(lnelev),aoverb(lnd),ldiam(lnd))
      Allocate(sb(2,2,lnd,lnelev))
      Allocate(sf11(lnd,lnelev),sf22(lnd,lnelev))
      ! 
      do ielev=1,lnelev
        !
        read(id,*) elev(ielev), work ! incident angle,- elevation =90-elev
        do ir=1,lnd
          !
          read(id,*) ldiam(ir),aoverb(ir)  ! radius [mm] and (1/aspect ratio) 
          read(id,*) sf11(ir,ielev),sf22(ir,ielev)
          read(id,*) sb(1,1,ir,ielev),sb(1,2,ir,ielev),sb(2,1,ir,ielev),sb(2,2,ir,ielev) ! mm
        enddo ! ir
        !
      enddo ! ielev
    close(id)
    !
    ! find elevation elev0
    !
    if (lnelev==1) ielev0=1
    !  
    if (lnelev>1) then
      ielev0=-1
      diff0=1.d+11
      do ielev=1,lnelev
        diff=dabs(elev(ielev)-(90.d0-elev0))
        if (diff<diff0) then
          diff0=diff
          ielev0=ielev
        endif
      enddo
    endif
    !
    !------------------------------
    ! save lut info in desired range
    ii=0
    do ir=ir1,ir2
      ii=ii+1
      out_diam(ii)=2.d0*ldiam(ir)
      out_aoverb(ii)=aoverb(ir)
      out_sf11(ii)=sf11(ir,ielev0)
      out_sf22(ii)=sf22(ir,ielev0)
      out_sb11(ii)=sb(1,1,ir,ielev0)
      out_sb22(ii)=sb(2,2,ir,ielev0)
    enddo
    ! 
    !do ii=1, nir
    !  write(*,*) 'out_sb22 : ', out_sb22(ii), ', out_sb11 : ', out_sb11(ii)
    !enddo
    !
    !--------------------------------
    !
    Deallocate(aoverb,ldiam)
    Deallocate(elev,sf11,sf22,sb)
    !
  return
  end subroutine read_luts
  !
  subroutine get_env_vars(conf,wrf,env)
  Use wrf_var_mod
  Use crsim_mod
  Use phys_param_mod, ONLY: Rd, eps, m999, T0K
  Implicit None
  !
  Type(conf_var),Intent(In)                  :: conf
  Type(wrf_var),Intent(InOut)                :: wrf
  Type(env_var),Intent(InOut)                :: env
  !
  Integer                                    :: ix1,ix2, iy1,iy2,iz1,iz2,it
  Real*8,Dimension(:),Allocatable            :: xtrack,ytrack
  Real*8,Dimension(:,:,:),Allocatable        :: gheight,ww
  Real*8,Dimension(:,:,:),Allocatable        :: uu,vv
  Real*8,Dimension(:,:,:),Allocatable        :: Ku,Kv,Kw
  Real*8                                     :: hres
  Integer                                    :: ix,iy,iz
  Integer                                    :: lowestiz
  Real*8                                     :: sfcheight
    !-------------------------------------------------------------
    ! 
    ix1=1 ; ix2=env%nx
    iy1=1 ; iy2=env%ny
    iz1=1 ; iz2=env%nz
    it = 1
    !
    If (wrf%nx/=env%nx) Then
      ix1=conf%ix_start ; ix2=conf%ix_end
    EndIf
    !
    If (wrf%ny/=env%ny) Then
      iy1=conf%iy_start ; iy2=conf%iy_end
    EndIf
    !
    If (wrf%nz/=env%nz) Then
      iz1=conf%iz_start ; iz2=conf%iz_end
    EndIf
    !
    If (wrf%nt>1) it=conf%it
    !
    ! ------------------------------------------------------------------------------
    ! ------------------------------------------------------------------------------
    ! get meteor. fields from wrf vars
    !
    wrf%press=(wrf%pb+wrf%p)  ! Pa
    wrf%temp=(wrf%theta+300.d0)*(wrf%press/100000.d0)**(287.d0/1004.d0) ! K 
    wrf%geop_height=(wrf%phb+wrf%ph)/9.81d0 ! m
    !
    !------------------------------------------------------------------------------------
    !-------------------------------------------------------------------------------
    !
    write(*,*) 'Info: Getting x-domain, y-domain and height'
    ! x, y , z
    Allocate(xtrack(wrf%nx),ytrack(wrf%ny),gheight(wrf%nx,wrf%ny,wrf%nz),&
            Kw(wrf%nx,wrf%ny,wrf%nz),Ku(wrf%nx,wrf%ny,wrf%nz),Kv(wrf%nx,wrf%ny,wrf%nz),&
            uu(wrf%nx,wrf%ny,wrf%nz),vv(wrf%nx,wrf%ny,wrf%nz),ww(wrf%nx,wrf%ny,wrf%nz))
    !
    hres=1.d0/wrf%dx(it)
    env%dx=hres
    do ix=1,wrf%nx
      xtrack(ix)=dble(ix-1)*hres
    enddo
    !
    hres=1.d0/wrf%dy(it)
    env%dy=hres
    do iy=1,wrf%ny
      ytrack(iy)=dble(iy-1)*hres
    enddo
    !
    !----------------
    do iz=1,wrf%nz
      gheight(:,:,iz)=0.5d0*(wrf%geop_height(:,:,iz,it)+wrf%geop_height(:,:,iz+1,it))
      ww(:,:,iz)=0.5d0*(wrf%w(:,:,iz,it)+wrf%w(:,:,iz+1,it))
      Kw(:,:,iz)= (wrf%w(:,:,iz+1,it)-wrf%w(:,:,iz,it)) / (wrf%geop_height(:,:,iz+1,it)-wrf%geop_height(:,:,iz,it))
      if(wrf%geop_height(1,1,1,1) > wrf%geop_height(1,1,wrf%nz,1)) Kw(:,:,iz)=Kw(:,:,iz)*(-1.d0) 
    enddo
    !
    if(wrf%geop_height(1,1,1,1) > wrf%geop_height(1,1,wrf%nz,1)) write(*,*) 'Decreasing height'
    !----------------
    do ix=1,wrf%nx
      uu(ix,:,:)=0.5d0*(wrf%u(ix,:,:,it)+wrf%u(ix+1,:,:,it))
      Ku(ix,:,:)=wrf%dx(it) * (wrf%u(ix+1,:,:,it)-wrf%u(ix,:,:,it))
    enddo
    !
    do iy=1,wrf%ny
      vv(:,iy,:)=0.5d0*(wrf%v(:,iy,:,it)+wrf%v(:,iy+1,:,it))
      Kv(:,iy,:)=wrf%dy(it) * (wrf%v(:,iy+1,:,it)-wrf%v(:,iy,:,it))
    enddo
    !
    !-----------------
    !
    env%x(1:env%nx)=xtrack(ix1:ix2)
    env%y(1:env%ny)=ytrack(iy1:iy2)
    env%z(1:env%nx,1:env%ny,1:env%nz)=gheight(ix1:ix2,iy1:iy2,iz1:iz2)
    env%w(1:env%nx,1:env%ny,1:env%nz)=ww(ix1:ix2,iy1:iy2,iz1:iz2)
    !
    env%u(1:env%nx,1:env%ny,1:env%nz)=uu(ix1:ix2,iy1:iy2,iz1:iz2)
    env%v(1:env%nx,1:env%ny,1:env%nz)=vv(ix1:ix2,iy1:iy2,iz1:iz2)
    !
    env%xlat(1:env%nx,1:env%ny)=wrf%xlat(ix1:ix2,iy1:iy2,it) ! deg
    env%xlong(1:env%nx,1:env%ny)=wrf%xlong(ix1:ix2,iy1:iy2,it) ! deg
    !
    env%Kw(1:env%nx,1:env%ny,1:env%nz)=Kw(ix1:ix2,iy1:iy2,iz1:iz2)
    env%Ku(1:env%nx,1:env%ny,1:env%nz)=Ku(ix1:ix2,iy1:iy2,iz1:iz2)
    env%Kv(1:env%nx,1:env%ny,1:env%nz)=Kv(ix1:ix2,iy1:iy2,iz1:iz2)
    
    !-------------------------              
    ! get surface temperature and wind speed
    env%sfctemp(1:env%nx,1:env%ny) = wrf%sfctemp(ix1:ix2,iy1:iy2,it) -T0K;
    lowestiz = 1;
    if(wrf%geop_height(1,1,1,1) > wrf%geop_height(1,1,wrf%nz,1)) lowestiz = wrf%nz;
    do ix=ix1,ix2
       do iy=iy1,iy2
          env%sfcwspd(ix-ix1+1,iy-iy1+1)=sqrt((wrf%sfcuu(ix,iy,it)**2.d0)+(wrf%sfcvv(ix,iy,it)**2.d0)) ;
          !mask land
          sfcheight = wrf%topo(ix,iy)
          if(wrf%topo(ix,iy)<-880.0) sfcheight = gheight(ix,iy,lowestiz)
          if(sfcheight>0.0) then
             env%sfctemp(ix-ix1+1,iy-iy1+1)=m999
             env%sfcwspd(ix-ix1+1,iy-iy1+1)=m999
          end if
       end do
    end do
    if(MaxVal(env%sfctemp)<0.d0) then ! if no surface data, use the lowest layer from 3D vaiables
       env%sfctemp(1:env%nx,1:env%ny)=wrf%temp(ix1:ix2,iy1:iy2,lowestiz,it) -T0K;
       do ix=ix1,ix2
          do iy=iy1,iy2
             !mask land
             sfcheight = wrf%topo(ix,iy)
             if(wrf%topo(ix,iy)<-880.0) sfcheight = gheight(ix,iy,lowestiz)
             if(sfcheight>0.0) then
                env%sfctemp(ix-ix1+1,iy-iy1+1)=m999
             end if
          end do
       end do
    end if
    if(MaxVal(env%sfcwspd)<=-888.d0) then ! if no surface data, use the lowest layer from 3D vaiables
       lowestiz = 1;
       if(wrf%geop_height(1,1,1,1) > wrf%geop_height(1,1,wrf%nz,1)) lowestiz = wrf%nz;
       do ix=ix1,ix2
          do iy=iy1,iy2
             env%sfcwspd(ix-ix1+1,iy-iy1+1)=sqrt((uu(ix,iy,lowestiz)**2.d0)+(vv(ix,iy,lowestiz)**2.d0)) ;
             !mask land
             sfcheight = wrf%topo(ix,iy)
             if(wrf%topo(ix,iy)<-880.0) sfcheight = gheight(ix,iy,lowestiz)
             if(sfcheight>0.0) then
                env%sfcwspd(ix-ix1+1,iy-iy1+1)=m999
             end if
          end do
       end do
    end if
    !---------------------------------------------------------------------------------
    Deallocate(xtrack,ytrack,gheight,ww,uu,vv,Kw,Ku,Kv)
    !---------------------------------------------------------------------------------
    !---------------------------------------------------------------------------------
    ! 
    write(*,*) 'Info: Getting press,temp,rho_d,qvapor,tke'
    env%press(1:env%nx,1:env%ny,1:env%nz)=wrf%press(ix1:ix2,iy1:iy2,iz1:iz2,it) ! Pa 
    env%temp(1:env%nx,1:env%ny,1:env%nz)=wrf%temp(ix1:ix2,iy1:iy2,iz1:iz2,it)   ! K
    env%qvapor(1:env%nx,1:env%ny,1:env%nz)=wrf%qvapor(ix1:ix2,iy1:iy2,iz1:iz2,it) ! kg/kg
    env%tke(1:env%nx,1:env%ny,1:env%nz)=wrf%tke(ix1:ix2,iy1:iy2,iz1:iz2,it) ! m^2/s^2
   !
    ! 
    if (MaxVal(wrf%rho_d)<-888.d0) then
      ! if ALT not present in the input file
      env%rho_d=(env%press*env%qvapor)/(eps+env%qvapor) ! e water vapor pressure in Pa 
      env%rho_d=(env%press-env%rho_d)/(Rd * env%temp)   !  kg/m^3    
    else 
      !if ALT present in the input file
      env%rho_d(1:env%nx,1:env%ny,1:env%nz)=wrf%rho_d(ix1:ix2,iy1:iy2,iz1:iz2,it) ! m^3/kg
      env%rho_d=1.d0/env%rho_d !  kg/m^3
    endif
    !
    ! convert T from K to C
    env%temp=env%temp-273.15d0 ! C
    ! convert press frpm Pa to mb
    env%press=env%press*1.d-2  ! mb
    !
    write(*,*) 'press [mb]',MinVal(env%press),MaxVal(env%press)
    write(*,*) 'temp [C]',MinVal(env%temp),MaxVal(env%temp)
    write(*,*) 'rho_d [kg/m^3]',MinVal(env%rho_d),MaxVal(env%rho_d)
    write(*,*) 'tke [m^2/s^2]',MinVal(env%tke),MaxVal(env%tke)
    !
  return
  end subroutine get_env_vars
  !!
  !!
  subroutine get_hydro10_vars(conf,mp10,hydro)
  Use wrf_var_mod
  Use crsim_mod
  Implicit None
  !
  Type(conf_var),Intent(In)                  :: conf
  Type(wrf_var_mp10),Intent(In)              :: mp10
  Type(hydro_var),Intent(InOut)              :: hydro
  !
  Integer                                    :: ix1,ix2, iy1,iy2,iz1,iz2,it
  Integer                                    :: morrID
  real*8                                     :: thr
  Integer                                    :: iz
    !
    ! default - the whole scene and it==1
    ix1=1 ; ix2=mp10%nx
    iy1=1 ; iy2=mp10%ny
    iz1=1 ; iz2=mp10%nz
    it = 1
    !
    ! reconstruct scene from configuration parameters (defined by user)
    If (hydro%nx/=mp10%nx) Then
      ix1=conf%ix_start ; ix2=conf%ix_end
    EndIf
    !
    If (hydro%ny/=mp10%ny) Then
      iy1=conf%iy_start ; iy2=conf%iy_end
    EndIf
    !
    If (hydro%nz/=mp10%nz) Then
      iz1=conf%iz_start ; iz2=conf%iz_end
    EndIf
    !
    If (mp10%nt>1) it=conf%it
    !
    hydro%qhydro(:,:,:,1)=mp10%qcloud(ix1:ix2,iy1:iy2,iz1:iz2,it)
    hydro%qhydro(:,:,:,2)=mp10%qrain(ix1:ix2,iy1:iy2,iz1:iz2,it)
    hydro%qhydro(:,:,:,3)=mp10%qice(ix1:ix2,iy1:iy2,iz1:iz2,it)
    hydro%qhydro(:,:,:,4)=mp10%qsnow(ix1:ix2,iy1:iy2,iz1:iz2,it)
    hydro%qhydro(:,:,:,5)=mp10%qgraup(ix1:ix2,iy1:iy2,iz1:iz2,it)
    !
    hydro%qnhydro(:,:,:,1)=mp10%qncloud(ix1:ix2,iy1:iy2,iz1:iz2,it)
    hydro%qnhydro(:,:,:,2)=mp10%qnrain(ix1:ix2,iy1:iy2,iz1:iz2,it)
    hydro%qnhydro(:,:,:,3)=mp10%qnice(ix1:ix2,iy1:iy2,iz1:iz2,it)
    hydro%qnhydro(:,:,:,4)=mp10%qnsnow(ix1:ix2,iy1:iy2,iz1:iz2,it)
    hydro%qnhydro(:,:,:,5)=mp10%qngraup(ix1:ix2,iy1:iy2,iz1:iz2,it)
    !
    write(*,*) 'Info: removing values below the given threshold'
    thr=1.d-50
    !-----------------------------------------------------------------------------
    morrID=0
    if (MaxVal(hydro%qnhydro(:,:,:,1))<0.d0) morrID=1
    
    do iz=1,hydro%nz
      !
      thr=conf%thr_mix_ratio(1)
      IF (morrID==0) then
        where(hydro%qhydro (:,:,iz,1) < thr)
          hydro%qhydro (:,:,iz,1) = 0.d0
          hydro%qnhydro(:,:,iz,1) = 0.d0
        endwhere
        where(hydro%qnhydro(:,:,iz,1) <= 0.d0)
          hydro%qhydro (:,:,iz,1) =0.d0
          hydro%qnhydro(:,:,iz,1) =0.d0
        endwhere
      ELSE
        where(hydro%qhydro (:,:,iz,1) < thr)
          hydro%qhydro (:,:,iz,1) = 0.d0
        endwhere
      ENDIF
      !
      !!----------------------------------
      thr=conf%thr_mix_ratio(2)
      where(hydro%qhydro (:,:,iz,2) < thr)
        hydro%qhydro (:,:,iz,2) = 0.d0
        hydro%qnhydro(:,:,iz,2) = 0.d0
      endwhere
      where(hydro%qnhydro(:,:,iz,2) <= 0.d0)
        hydro%qhydro (:,:,iz,2) =0.d0
        hydro%qnhydro(:,:,iz,2) =0.d0
      endwhere
      !!----------------------------------
    
      thr=conf%thr_mix_ratio(3)
      where(hydro%qhydro (:,:,iz,3) < thr)
        hydro%qhydro (:,:,iz,3) = 0.d0
        hydro%qnhydro(:,:,iz,3) = 0.d0
      endwhere
      where(hydro%qnhydro(:,:,iz,3) <= 0.d0)
        hydro%qhydro (:,:,iz,3) = 0.d0
        hydro%qnhydro(:,:,iz,3) = 0.d0
      endwhere
      !!----------------------------------
      !!----------------------------------
      thr=conf%thr_mix_ratio(4)
      where(hydro%qhydro (:,:,iz,4) < thr)
        hydro%qhydro (:,:,iz,4) = 0.d0
        hydro%qnhydro(:,:,iz,4) = 0.d0
      endwhere
      where(hydro%qnhydro(:,:,iz,4) <=0.d0)
        hydro%qhydro (:,:,iz,4) = 0.d0
        hydro%qnhydro(:,:,iz,4) = 0.d0
      endwhere
      !!----------------------------------
      !!----------------------------------
      thr=conf%thr_mix_ratio(5)
      where(hydro%qhydro (:,:,iz,5) < thr)
        hydro%qhydro (:,:,iz,5) = 0.d0
        hydro%qnhydro(:,:,iz,5) = 0.d0
      endwhere
      where(hydro%qnhydro(:,:,iz,5) <= 0.d0)
        hydro%qhydro (:,:,iz,5) = 0.d0
        hydro%qnhydro(:,:,iz,5) = 0.d0
      endwhere
      !!----------------------------------
      !!
    enddo ! iz
    !  
  return
  end subroutine  get_hydro10_vars
  !     
  !-------------------------------------------------------------------------------------------------------------------
  !!! Added by oue, 2016/09/19 , fixed 2017/05/26
  subroutine get_hydro08_vars(conf,mp08,hydro)
  Use wrf_var_mod
  Use crsim_mod
  Implicit None
  !
  Type(conf_var),Intent(In)                  :: conf
  Type(wrf_var_mp08),Intent(In)              :: mp08
  Type(hydro_var),Intent(InOut)              :: hydro
  !
  Integer                                    :: ix1,ix2, iy1,iy2,iz1,iz2,it
  Integer                                    :: morrID
  real*8                                     :: thr
  Integer                                    :: iz,ih
    !
    ! default - the whole scene and it==1
    ix1=1 ; ix2=mp08%nx
    iy1=1 ; iy2=mp08%ny
    iz1=1 ; iz2=mp08%nz
    it = 1
    !
    ! reconstruct scene from configuration parameters (defined by user)
    If (hydro%nx/=mp08%nx) Then
      ix1=conf%ix_start ; ix2=conf%ix_end
    EndIf
    !
    If (hydro%ny/=mp08%ny) Then
      iy1=conf%iy_start ; iy2=conf%iy_end
    EndIf
    !
    If (hydro%nz/=mp08%nz) Then
      iz1=conf%iz_start ; iz2=conf%iz_end
    EndIf
    !
    If (mp08%nt>1) it=conf%it
    !
    hydro%qhydro(:,:,:,1)=mp08%qcloud(ix1:ix2,iy1:iy2,iz1:iz2,it)
    hydro%qhydro(:,:,:,2)=mp08%qrain(ix1:ix2,iy1:iy2,iz1:iz2,it)
    hydro%qhydro(:,:,:,3)=mp08%qice(ix1:ix2,iy1:iy2,iz1:iz2,it)
    hydro%qhydro(:,:,:,4)=mp08%qsnow(ix1:ix2,iy1:iy2,iz1:iz2,it)
    hydro%qhydro(:,:,:,5)=mp08%qgraup(ix1:ix2,iy1:iy2,iz1:iz2,it)
    
    !
    hydro%qnhydro(:,:,:,1)=mp08%qncloud(ix1:ix2,iy1:iy2,iz1:iz2,it)
    hydro%qnhydro(:,:,:,2)=mp08%qnrain(ix1:ix2,iy1:iy2,iz1:iz2,it)
    hydro%qnhydro(:,:,:,3)=mp08%qnice(ix1:ix2,iy1:iy2,iz1:iz2,it)
    hydro%qnhydro(:,:,:,4)=mp08%qnsnow(ix1:ix2,iy1:iy2,iz1:iz2,it)
    hydro%qnhydro(:,:,:,5)=mp08%qngraup(ix1:ix2,iy1:iy2,iz1:iz2,it)
    !
    write(*,*) 'Info: removing values below the given threshold'
    thr=1.d-50
    !-----------------------------------------------------------------------------
    ! 
    morrID=0
    if (MaxVal(hydro%qnhydro(:,:,:,1))<0.d0) morrID=1
    do iz=1,hydro%nz
      !
      do ih=1,5
        morrID=0
        if (MaxVal(hydro%qnhydro(:,:,:,ih))<0.d0) morrID=1
        thr=conf%thr_mix_ratio(ih)
        IF (morrID==0) then
          where(hydro%qhydro (:,:,iz,ih) < thr)
            hydro%qhydro (:,:,iz,ih) = 0.d0
            hydro%qnhydro(:,:,iz,ih) = 0.d0
          endwhere
          where(hydro%qnhydro(:,:,iz,ih) <= 0.d0)
            !hydro%qhydro (:,:,iz,ih) =0.d0 !fixed by oue Apr. 2020
            hydro%qnhydro(:,:,iz,ih) =0.d0
          endwhere
        ELSE
          where(hydro%qhydro (:,:,iz,ih) < thr)
            hydro%qhydro (:,:,iz,ih) = 0.d0
          endwhere
        ENDIF
        !!----------------------------------
      end do ! ih loop
      !!
    enddo ! iz
    ! 
  return
  end subroutine  get_hydro08_vars
  !!-- addedd by oue
  !--------------------------------------------------------------------
  subroutine get_hydro09_vars(conf,mp09,hydro)
  Use wrf_var_mod
  Use crsim_mod
  Implicit None
  !
  Type(conf_var),Intent(In)                  :: conf
  Type(wrf_var_mp09),Intent(In)              :: mp09
  Type(hydro_var),Intent(InOut)              :: hydro
  !
  Integer                                    :: ix1,ix2, iy1,iy2,iz1,iz2,it
  real*8                                     :: thr
   
  Integer                                    :: iz
    !
    ! default - the whole scene and it==1
    ix1=1 ; ix2=mp09%nx
    iy1=1 ; iy2=mp09%ny
    iz1=1 ; iz2=mp09%nz
    it = 1
    !
    ! reconstruct scene from configuration parameters (defined by user)
    If (hydro%nx/=mp09%nx) Then
      ix1=conf%ix_start ; ix2=conf%ix_end
    EndIf
    !
    If (hydro%ny/=mp09%ny) Then
      iy1=conf%iy_start ; iy2=conf%iy_end
    EndIf
    !
    If (hydro%nz/=mp09%nz) Then
      iz1=conf%iz_start ; iz2=conf%iz_end
    EndIf
    !
    If (mp09%nt>1) it=conf%it
    ! 
    hydro%qhydro(:,:,:,1)=mp09%qcloud(ix1:ix2,iy1:iy2,iz1:iz2,it)
    hydro%qhydro(:,:,:,2)=mp09%qrain(ix1:ix2,iy1:iy2,iz1:iz2,it)
    hydro%qhydro(:,:,:,3)=mp09%qice(ix1:ix2,iy1:iy2,iz1:iz2,it)
    hydro%qhydro(:,:,:,4)=mp09%qsnow(ix1:ix2,iy1:iy2,iz1:iz2,it)
    hydro%qhydro(:,:,:,5)=mp09%qgraup(ix1:ix2,iy1:iy2,iz1:iz2,it)
    hydro%qhydro(:,:,:,6)=mp09%qhail(ix1:ix2,iy1:iy2,iz1:iz2,it)
    !
    hydro%qnhydro(:,:,:,1)=mp09%qncloud(ix1:ix2,iy1:iy2,iz1:iz2,it)
    hydro%qnhydro(:,:,:,2)=mp09%qnrain(ix1:ix2,iy1:iy2,iz1:iz2,it)
    hydro%qnhydro(:,:,:,3)=mp09%qnice(ix1:ix2,iy1:iy2,iz1:iz2,it)
    hydro%qnhydro(:,:,:,4)=mp09%qnsnow(ix1:ix2,iy1:iy2,iz1:iz2,it)
    hydro%qnhydro(:,:,:,5)=mp09%qngraup(ix1:ix2,iy1:iy2,iz1:iz2,it)
    hydro%qnhydro(:,:,:,6)=mp09%qnhail(ix1:ix2,iy1:iy2,iz1:iz2,it)
    !
    write(*,*) 'Info: removing values below the given threshold'
    thr=1.d-50
    !-----------------------------------------------------------------------------
    do iz=1,hydro%nz
      !
      thr=conf%thr_mix_ratio(1)
      where(hydro%qhydro (:,:,iz,1) < thr)
        hydro%qhydro (:,:,iz,1) = 0.d0
        hydro%qnhydro(:,:,iz,1) = 0.d0
      endwhere
      where(hydro%qnhydro(:,:,iz,1) <= 0.d0)
        hydro%qhydro (:,:,iz,1) =0.d0
        hydro%qnhydro(:,:,iz,1) =0.d0
      endwhere
      !!----------------------------------
      thr=conf%thr_mix_ratio(2)
      where(hydro%qhydro (:,:,iz,2) < thr) 
        hydro%qhydro (:,:,iz,2) = 0.d0
        hydro%qnhydro(:,:,iz,2) = 0.d0
      endwhere
      where(hydro%qnhydro(:,:,iz,2) <= 0.d0)
        hydro%qhydro (:,:,iz,2) =0.d0 
        hydro%qnhydro(:,:,iz,2) =0.d0
      endwhere
      !!----------------------------------
      thr=conf%thr_mix_ratio(3)
      where(hydro%qhydro (:,:,iz,3) < thr)
        hydro%qhydro (:,:,iz,3) = 0.d0
        hydro%qnhydro(:,:,iz,3) = 0.d0
      endwhere
      where(hydro%qnhydro(:,:,iz,3) <= 0.d0)
        hydro%qhydro (:,:,iz,3) = 0.d0
        hydro%qnhydro(:,:,iz,3) = 0.d0
      endwhere
      !!----------------------------------
      !!----------------------------------
      thr=conf%thr_mix_ratio(4)
      where(hydro%qhydro (:,:,iz,4) < thr)
        hydro%qhydro (:,:,iz,4) = 0.d0
        hydro%qnhydro(:,:,iz,4) = 0.d0
      endwhere
      where(hydro%qnhydro(:,:,iz,4) <=0.d0)
        hydro%qhydro (:,:,iz,4) = 0.d0
        hydro%qnhydro(:,:,iz,4) = 0.d0
      endwhere
      !!----------------------------------
      !!----------------------------------
      thr=conf%thr_mix_ratio(5)
      where(hydro%qhydro (:,:,iz,5) < thr)
        hydro%qhydro (:,:,iz,5) = 0.d0
        hydro%qnhydro(:,:,iz,5) = 0.d0
      endwhere
      where(hydro%qnhydro(:,:,iz,5) <= 0.d0)
        hydro%qhydro (:,:,iz,5) = 0.d0
        hydro%qnhydro(:,:,iz,5) = 0.d0
      endwhere
      !!----------------------------------
      thr=conf%thr_mix_ratio(6)
      where(hydro%qhydro (:,:,iz,6) < thr)
        hydro%qhydro (:,:,iz,6) = 0.d0
        hydro%qnhydro(:,:,iz,6) = 0.d0
      endwhere
      where(hydro%qnhydro(:,:,iz,6) <= 0.d0)
        hydro%qhydro (:,:,iz,6) = 0.d0
        hydro%qnhydro(:,:,iz,6) = 0.d0
      endwhere
    
    enddo ! iz
    !!
    !go to 102
    write(*,*) '--------------------------------------'
    write(*,*) '-mixing ratio-----------------------------------'
    write(*,*) 'cloud',minval(hydro%qhydro(:,:,:,1)),maxval(hydro%qhydro(:,:,:,1))
    write(*,*) 'rain ',minval(hydro%qhydro(:,:,:,2)),maxval(hydro%qhydro(:,:,:,2))
    write(*,*) 'ice  ',minval(hydro%qhydro(:,:,:,3)),maxval(hydro%qhydro(:,:,:,3))
    write(*,*) 'snow ',minval(hydro%qhydro(:,:,:,4)),maxval(hydro%qhydro(:,:,:,4))
    write(*,*) 'graup',minval(hydro%qhydro(:,:,:,5)),maxval(hydro%qhydro(:,:,:,5))
    write(*,*) 'hail',minval(hydro%qhydro(:,:,:,6)),maxval(hydro%qhydro(:,:,:,6))
    
    write(*,*) '--------------------------------------'
    !
    write(*,*) '--------------------------------------'
    write(*,*) '---concentration-------------------------------'
    write(*,*) 'cloud',minval(hydro%qnhydro(:,:,:,1)),maxval(hydro%qnhydro(:,:,:,1))
    write(*,*) 'rain ',minval(hydro%qnhydro(:,:,:,2)),maxval(hydro%qnhydro(:,:,:,2))
    write(*,*) 'ice  ',minval(hydro%qnhydro(:,:,:,3)),maxval(hydro%qnhydro(:,:,:,3))
    write(*,*) 'snow ',minval(hydro%qnhydro(:,:,:,4)),maxval(hydro%qnhydro(:,:,:,4))
    write(*,*) 'graup',minval(hydro%qnhydro(:,:,:,5)),maxval(hydro%qnhydro(:,:,:,5))
    write(*,*) 'hail',minval(hydro%qnhydro(:,:,:,6)),maxval(hydro%qnhydro(:,:,:,6))
    !
    write(*,*) '--------------------------------------'
    
    !102 continue
    !
  return
  end subroutine get_hydro09_vars
  !!
  !-------------------------------------------------------------------------------------------------------------------
  !!! Added by oue, 2017/07/17 for ICON
  subroutine get_hydro30_vars(conf,mp30,hydro)
  Use wrf_var_mod
  Use crsim_mod
  Implicit None
  !
  Type(conf_var),Intent(In)                  :: conf
  Type(wrf_var_mp30),Intent(In)              :: mp30
  Type(hydro_var),Intent(InOut)              :: hydro
  !
  Integer                                    :: ix1,ix2, iy1,iy2,iz1,iz2,it
  Integer                                    :: morrID
  real*8                                     :: thr
  
  Integer                                    :: iz,ih
    !
    ! default - the whole scene and it==1
    ix1=1 ; ix2=mp30%nx
    iy1=1 ; iy2=mp30%ny
    iz1=1 ; iz2=mp30%nz
    it = 1
    !
    ! reconstruct scene from configuration parameters (defined by user)
    If (hydro%nx/=mp30%nx) Then
      ix1=conf%ix_start ; ix2=conf%ix_end
    EndIf
    !
    If (hydro%ny/=mp30%ny) Then
      iy1=conf%iy_start ; iy2=conf%iy_end
    EndIf
    !
    If (hydro%nz/=mp30%nz) Then
      iz1=conf%iz_start ; iz2=conf%iz_end
    EndIf
    !
    If (mp30%nt>1) it=conf%it
    !
    hydro%qhydro(:,:,:,1)=mp30%qcloud(ix1:ix2,iy1:iy2,iz1:iz2,it)
    hydro%qhydro(:,:,:,2)=mp30%qrain(ix1:ix2,iy1:iy2,iz1:iz2,it)
    hydro%qhydro(:,:,:,3)=mp30%qice(ix1:ix2,iy1:iy2,iz1:iz2,it)
    hydro%qhydro(:,:,:,4)=mp30%qsnow(ix1:ix2,iy1:iy2,iz1:iz2,it)
    hydro%qhydro(:,:,:,5)=mp30%qgraup(ix1:ix2,iy1:iy2,iz1:iz2,it)
    hydro%qhydro(:,:,:,6)=mp30%qhail(ix1:ix2,iy1:iy2,iz1:iz2,it)
    !
    hydro%qnhydro(:,:,:,1)=mp30%qncloud(ix1:ix2,iy1:iy2,iz1:iz2,it)
    hydro%qnhydro(:,:,:,2)=mp30%qnrain(ix1:ix2,iy1:iy2,iz1:iz2,it)
    hydro%qnhydro(:,:,:,3)=mp30%qnice(ix1:ix2,iy1:iy2,iz1:iz2,it)
    hydro%qnhydro(:,:,:,4)=mp30%qnsnow(ix1:ix2,iy1:iy2,iz1:iz2,it)
    hydro%qnhydro(:,:,:,5)=mp30%qngraup(ix1:ix2,iy1:iy2,iz1:iz2,it)
    hydro%qnhydro(:,:,:,6)=mp30%qnhail(ix1:ix2,iy1:iy2,iz1:iz2,it)
    !
    write(*,*) 'Info: removing values below the given threshold'
    thr=1.d-50
    !-----------------------------------------------------------------------------
    morrID=0
    if (MaxVal(hydro%qnhydro(:,:,:,1))<0.d0) morrID=1
    ! 
    do iz=1,hydro%nz
      do ih=1,6
        !
        morrID=0
        if (MaxVal(hydro%qnhydro(:,:,:,ih))<0.d0) morrID=1
        thr=conf%thr_mix_ratio(ih)
        IF (morrID==0) then
          where(hydro%qhydro (:,:,iz,ih) < thr)
            hydro%qhydro (:,:,iz,ih) = 0.d0
            hydro%qnhydro(:,:,iz,ih) = 0.d0
          endwhere
          where(hydro%qnhydro(:,:,iz,ih) <= 0.d0)
            hydro%qhydro (:,:,iz,ih) =0.d0
            hydro%qnhydro(:,:,iz,ih) =0.d0
          endwhere
        ELSE
          where(hydro%qhydro (:,:,iz,ih) < thr)
            hydro%qhydro (:,:,iz,ih) = 0.d0
          endwhere
        ENDIF
        !
      end do ! ih loop
      !
    enddo ! iz
    !
  return
  end subroutine  get_hydro30_vars
  !!-- addedd by oue
  !--------------------------------------------------------------------
  !!! Added by oue, 2017/07/21 for RAMS
  subroutine get_hydro40_vars(conf,mp40,hydro)
  Use wrf_var_mod
  Use crsim_mod
  Implicit None
  !
  Type(conf_var),Intent(In)                  :: conf
  Type(wrf_var_mp40),Intent(In)              :: mp40
  Type(hydro_var),Intent(InOut)              :: hydro
  !
  Integer                                    :: ix1,ix2, iy1,iy2,iz1,iz2,it
  Integer                                    :: morrID
  real*8                                     :: thr
  Integer                                    :: iz,ih
    !
    ! default - the whole scene and it==1
    ix1=1 ; ix2=mp40%nx
    iy1=1 ; iy2=mp40%ny
    iz1=1 ; iz2=mp40%nz
    it = 1
    !
    ! reconstruct scene from configuration parameters (defined by user)
    If (hydro%nx/=mp40%nx) Then
      ix1=conf%ix_start ; ix2=conf%ix_end
    EndIf
    !
    If (hydro%ny/=mp40%ny) Then
      iy1=conf%iy_start ; iy2=conf%iy_end
    EndIf
    !
    If (hydro%nz/=mp40%nz) Then
      iz1=conf%iz_start ; iz2=conf%iz_end
    EndIf
    !
    If (mp40%nt>1) it=conf%it
    !
    hydro%qhydro(:,:,:,1)=mp40%qcloud(ix1:ix2,iy1:iy2,iz1:iz2,it)
    hydro%qhydro(:,:,:,2)=mp40%qrain(ix1:ix2,iy1:iy2,iz1:iz2,it)
    hydro%qhydro(:,:,:,3)=mp40%qice(ix1:ix2,iy1:iy2,iz1:iz2,it)
    hydro%qhydro(:,:,:,4)=mp40%qsnow(ix1:ix2,iy1:iy2,iz1:iz2,it)
    hydro%qhydro(:,:,:,5)=mp40%qgraup(ix1:ix2,iy1:iy2,iz1:iz2,it)
    hydro%qhydro(:,:,:,6)=mp40%qhail(ix1:ix2,iy1:iy2,iz1:iz2,it)
    hydro%qhydro(:,:,:,7)=mp40%qdrzl(ix1:ix2,iy1:iy2,iz1:iz2,it)
    hydro%qhydro(:,:,:,8)=mp40%qaggr(ix1:ix2,iy1:iy2,iz1:iz2,it)
    !
    hydro%qnhydro(:,:,:,1)=mp40%qncloud(ix1:ix2,iy1:iy2,iz1:iz2,it)
    hydro%qnhydro(:,:,:,2)=mp40%qnrain(ix1:ix2,iy1:iy2,iz1:iz2,it)
    hydro%qnhydro(:,:,:,3)=mp40%qnice(ix1:ix2,iy1:iy2,iz1:iz2,it)
    hydro%qnhydro(:,:,:,4)=mp40%qnsnow(ix1:ix2,iy1:iy2,iz1:iz2,it)
    hydro%qnhydro(:,:,:,5)=mp40%qngraup(ix1:ix2,iy1:iy2,iz1:iz2,it)
    hydro%qnhydro(:,:,:,6)=mp40%qnhail(ix1:ix2,iy1:iy2,iz1:iz2,it)
    hydro%qnhydro(:,:,:,7)=mp40%qndrzl(ix1:ix2,iy1:iy2,iz1:iz2,it)
    hydro%qnhydro(:,:,:,8)=mp40%qnaggr(ix1:ix2,iy1:iy2,iz1:iz2,it)
    !
    !
    write(*,*) 'Info: removing values below the given threshold'
    thr=1.d-50
    !-----------------------------------------------------------------------------
    morrID=0
    if (MaxVal(hydro%qnhydro(:,:,:,1))<0.d0) morrID=1
    !
    do iz=1,hydro%nz
      !
      do ih=1,8
        !----------------------------------
        morrID=0
        if (MaxVal(hydro%qnhydro(:,:,:,ih))<0.d0) morrID=1
        thr=conf%thr_mix_ratio(ih)
        IF (morrID==0) then
          where(hydro%qhydro (:,:,iz,ih) < thr)
            hydro%qhydro (:,:,iz,ih) = 0.d0
            hydro%qnhydro(:,:,iz,ih) = 0.d0
          endwhere
          where(hydro%qnhydro(:,:,iz,ih) <= 0.d0)
            hydro%qhydro (:,:,iz,ih) =0.d0
            hydro%qnhydro(:,:,iz,ih) =0.d0
          endwhere
        ELSE
          where(hydro%qhydro (:,:,iz,ih) < thr)
            hydro%qhydro (:,:,iz,ih) = 0.d0
          endwhere
        ENDIF
        !----------------------------------
      end do ! ih loop
      !
    enddo ! iz
    !  
  return
  end subroutine  get_hydro40_vars
  !!-- addedd by oue
  !--------------------------------------------------------------------
  !
  ! subroutine get_hydro50_vars is added by DW 2017/10/30 for P3
  subroutine get_hydro50_vars(conf,mp50,hydro50)
  !
  Use wrf_var_mod
  Use crsim_mod
  Implicit None
  !
  Type(conf_var),Intent(In)                  :: conf
  Type(wrf_var_mp50),Intent(In)              :: mp50
  Type(hydro50_var),Intent(InOut)            :: hydro50
  !
  Integer                                    :: ix1,ix2,iy1,iy2,iz1,iz2,it
  Integer                                    :: morrID
  real*8                                     :: thr
  Integer                                    :: iz
    !
    ! default - the whole scene and it==1
    ix1=1 ; ix2=mp50%nx
    iy1=1 ; iy2=mp50%ny
    iz1=1 ; iz2=mp50%nz
    it=1
    !
    ! reconstruct scene from configuration parameters (defined by user)
    If (hydro50%nx/=mp50%nx) Then
      ix1=conf%ix_start ; ix2=conf%ix_end
    EndIf
    !
    If (hydro50%ny/=mp50%ny) Then
      iy1=conf%iy_start ; iy2=conf%iy_end
    EndIf
    !
    If (hydro50%nz/=mp50%nz) Then
      iz1=conf%iz_start ; iz2=conf%iz_end
    EndIf
    !
    If (mp50%nt>1) it=conf%it
    !
    hydro50%qhydro(:,:,:,1)=mp50%qcloud(ix1:ix2,iy1:iy2,iz1:iz2,it)
    hydro50%qhydro(:,:,:,2)=mp50%qrain(ix1:ix2,iy1:iy2,iz1:iz2,it)
    hydro50%qhydro(:,:,:,3)=mp50%qice(ix1:ix2,iy1:iy2,iz1:iz2,it)
    hydro50%qhydro(:,:,:,4)=mp50%qice(ix1:ix2,iy1:iy2,iz1:iz2,it)
    hydro50%qhydro(:,:,:,5)=mp50%qice(ix1:ix2,iy1:iy2,iz1:iz2,it)
    hydro50%qhydro(:,:,:,6)=mp50%qice(ix1:ix2,iy1:iy2,iz1:iz2,it)
    !
    hydro50%qnhydro(:,:,:,1)=mp50%qncloud(ix1:ix2,iy1:iy2,iz1:iz2,it)
    hydro50%qnhydro(:,:,:,2)=mp50%qnrain(ix1:ix2,iy1:iy2,iz1:iz2,it)
    hydro50%qnhydro(:,:,:,3)=mp50%qnice(ix1:ix2,iy1:iy2,iz1:iz2,it)
    hydro50%qnhydro(:,:,:,4)=mp50%qnice(ix1:ix2,iy1:iy2,iz1:iz2,it)
    hydro50%qnhydro(:,:,:,5)=mp50%qnice(ix1:ix2,iy1:iy2,iz1:iz2,it)
    hydro50%qnhydro(:,:,:,6)=mp50%qnice(ix1:ix2,iy1:iy2,iz1:iz2,it)
    !
    hydro50%qir(:,:,:)=mp50%qir(ix1:ix2,iy1:iy2,iz1:iz2,it)
    hydro50%qib(:,:,:)=mp50%qib(ix1:ix2,iy1:iy2,iz1:iz2,it)
    !
    write(*,*) 'Info: removing values below the given threshold'
    thr=1.d-50
    !
    morrID=0
    if (MaxVal(hydro50%qnhydro(:,:,:,1))<0.d0) morrID=1
    !
    do iz=1,hydro50%nz
      !
      thr=conf%thr_mix_ratio(1)
      IF (morrID==0) then
      where(hydro50%qhydro (:,:,iz,1) < thr)
        hydro50%qhydro (:,:,iz,1) = 0.d0
        hydro50%qnhydro(:,:,iz,1) = 0.d0
      endwhere
      where(hydro50%qnhydro(:,:,iz,1) <= 0.d0)
        hydro50%qhydro (:,:,iz,1) =0.d0
        hydro50%qnhydro(:,:,iz,1) =0.d0
      endwhere
      ELSE
      where(hydro50%qhydro (:,:,iz,1) < thr)
        hydro50%qhydro (:,:,iz,1) = 0.d0
      endwhere
      ENDIF
      !
      thr=conf%thr_mix_ratio(2)
      !
      where(hydro50%qhydro (:,:,iz,2) < thr)
        hydro50%qhydro (:,:,iz,2) = 0.d0
        hydro50%qnhydro(:,:,iz,2) = 0.d0
      endwhere
      where(hydro50%qnhydro(:,:,iz,2) <= 0.d0)
        hydro50%qhydro (:,:,iz,2) =0.d0
        hydro50%qnhydro(:,:,iz,2) =0.d0
      endwhere
      !
      thr=conf%thr_mix_ratio(3)
      where(hydro50%qhydro (:,:,iz,3) < thr)
        hydro50%qhydro (:,:,iz,3) = 0.d0
        hydro50%qnhydro(:,:,iz,3) = 0.d0
      endwhere
      where(hydro50%qnhydro(:,:,iz,3) <= 0.d0)
        hydro50%qhydro (:,:,iz,3) = 0.d0
        hydro50%qnhydro(:,:,iz,3) = 0.d0
      endwhere
      !!----------------------------------
      thr=conf%thr_mix_ratio(4)
      where(hydro50%qhydro (:,:,iz,4) < thr)
        hydro50%qhydro (:,:,iz,4) = 0.d0
        hydro50%qnhydro(:,:,iz,4) = 0.d0
      endwhere
      where(hydro50%qnhydro(:,:,iz,4) <= 0.d0)
        hydro50%qhydro (:,:,iz,4) = 0.d0
        hydro50%qnhydro(:,:,iz,4) = 0.d0
      endwhere
      !!----------------------------------
      thr=conf%thr_mix_ratio(5)
      where(hydro50%qhydro (:,:,iz,5) < thr)
        hydro50%qhydro (:,:,iz,5) = 0.d0
        hydro50%qnhydro(:,:,iz,5) = 0.d0
      endwhere
      where(hydro50%qnhydro(:,:,iz,5) <= 0.d0)
        hydro50%qhydro (:,:,iz,5) = 0.d0
        hydro50%qnhydro(:,:,iz,5) = 0.d0
      endwhere
      !!----------------------------------
      thr=conf%thr_mix_ratio(6)
      where(hydro50%qhydro (:,:,iz,6) < thr)
        hydro50%qhydro (:,:,iz,6) = 0.d0
        hydro50%qnhydro(:,:,iz,6) = 0.d0
      endwhere
      where(hydro50%qnhydro(:,:,iz,6) <= 0.d0)
        hydro50%qhydro (:,:,iz,6) = 0.d0
        hydro50%qnhydro(:,:,iz,6) = 0.d0
      endwhere
      !!----------------------------------
      thr=conf%thr_mix_ratio(4)
      where(hydro50%qib (:,:,iz) < thr)
        hydro50%qib (:,:,iz) = 0.d0
      endwhere
      where(hydro50%qib(:,:,iz) <= 0.d0)
        hydro50%qib (:,:,iz) = 0.d0
      endwhere
      !!----------------------------------
      thr=conf%thr_mix_ratio(4)
      where(hydro50%qir(:,:,iz) < thr)
        hydro50%qir(:,:,iz) = 0.d0
      endwhere
      where(hydro50%qir(:,:,iz) <= 0.d0)
        hydro50%qir(:,:,iz) = 0.d0
      endwhere
      !!
    enddo ! iz
    !
  return
  end subroutine  get_hydro50_vars
  !
  !-----------------------------------------------------------------------------
  ! 
  subroutine get_hydro20_vars(conf,mp20,hydro20)
  Use wrf_var_mod
  Use crsim_mod
  Implicit None
  !
  Type(conf_var),Intent(In)                  :: conf
  Type(wrf_var_mp20),Intent(In)              :: mp20
  Type(hydro20_var),Intent(InOut)            :: hydro20
  !
  Integer                                    :: ix1,ix2, iy1,iy2,iz1,iz2,it
    !
    ! default - the whole scene and it==1
    ix1=1 ; ix2=mp20%nx
    iy1=1 ; iy2=mp20%ny
    iz1=1 ; iz2=mp20%nz
    it = 1
    !
    ! reconstruct scene from configuration parameters (defined by user)
    If (hydro20%nx/=mp20%nx) Then
      ix1=conf%ix_start ; ix2=conf%ix_end
    EndIf
    !
    If (hydro20%ny/=mp20%ny) Then
      iy1=conf%iy_start ; iy2=conf%iy_end
    EndIf
    !
    If (hydro20%nz/=mp20%nz) Then
      iz1=conf%iz_start ; iz2=conf%iz_end
    EndIf
    !
    If (mp20%nt>1) it=conf%it
    !-----------------------------------------------------------------
    hydro20%ff1(:,:,:,:)=mp20%ff1(ix1:ix2,iy1:iy2,iz1:iz2,it,:)
    hydro20%ff5(:,:,:,:)=mp20%ff5(ix1:ix2,iy1:iy2,iz1:iz2,it,:)
    hydro20%ff6(:,:,:,:)=mp20%ff6(ix1:ix2,iy1:iy2,iz1:iz2,it,:)
    !-----------------------------------------------------------------
    !
    !-----------------------------------------------------
    if (minval(hydro20%ff1)<0.d0) then
      write(*,*) 'Info: removing negative values for ff1'
      where(hydro20%ff1<0.d0) hydro20%ff1=0.d0
    endif
    !
    if (minval(hydro20%ff5)<0.d0) then
      write(*,*) 'Info: removing negative values for ff5'
      where(hydro20%ff5<0.d0) hydro20%ff5=0.d0
    endif
    !
    if (minval(hydro20%ff6)<0.d0) then
      write(*,*) 'Info: removing negative values for ff6'
      where(hydro20%ff6<0.d0) hydro20%ff6=0.d0
    endif
    !
  return
  end subroutine get_hydro20_vars
  !
  !-----------------------------------------------------------------------------
  ! 
  ! SAM warm bin by oue
  subroutine get_hydro70_vars(conf,mp70,hydro70)
  Use wrf_var_mod
  Use crsim_mod
  Implicit None
  !
  Type(conf_var),Intent(In)                  :: conf
  Type(wrf_var_mp70),Intent(In)              :: mp70
  Type(hydro70_var),Intent(InOut)            :: hydro70
  !
  Real*8                                     :: thr
  Integer                                    :: ix1,ix2, iy1,iy2,iz1,iz2,it
    !
    ! default - the whole scene and it==1
    ix1=1 ; ix2=mp70%nx
    iy1=1 ; iy2=mp70%ny
    iz1=1 ; iz2=mp70%nz
    it = 1
    !
    ! reconstruct scene from configuration parameters (defined by user)
    If (hydro70%nx/=mp70%nx) Then
      ix1=conf%ix_start ; ix2=conf%ix_end
    EndIf
    !
    If (hydro70%ny/=mp70%ny) Then
      iy1=conf%iy_start ; iy2=conf%iy_end
    EndIf
    !
    If (hydro70%nz/=mp70%nz) Then
      iz1=conf%iz_start ; iz2=conf%iz_end
    EndIf
    !
    If (mp70%nt>1) it=conf%it
    
    write(*,*) it,iz1,iz2,hydro70%nz,ix1,ix2,iy1,iy2
    !-----------------------------------------------------------------
    hydro70%fm1(:,:,:,:)=mp70%fm1(ix1:ix2,iy1:iy2,iz1:iz2,it,:) * 1.d-3 !g/kg -> kg/kg
    hydro70%fn1(:,:,:,:)=mp70%fn1(ix1:ix2,iy1:iy2,iz1:iz2,it,:) * 1.d6 !1/cm3 -> 1/m3
    ! 
    hydro70%qhydro(:,:,:,1)=mp70%qcloud(ix1:ix2,iy1:iy2,iz1:iz2,it) * 1.d-3 !g/kg -> kg/kg
    hydro70%qhydro(:,:,:,2)=mp70%qrain(ix1:ix2,iy1:iy2,iz1:iz2,it) * 1.d-3 !g/kg -> kg/kg
    !
    hydro70%qnhydro(:,:,:,1)=mp70%qncloud(ix1:ix2,iy1:iy2,iz1:iz2,it) * 1.d6 !1/cm3 -> 1/m3
    hydro70%qnhydro(:,:,:,2)=mp70%qnrain(ix1:ix2,iy1:iy2,iz1:iz2,it) * 1.d6 !1/cm3 -> 1/m3
    !
    hydro70%naero(:,:,:)=mp70%naero(ix1:ix2,iy1:iy2,iz1:iz2,it) * 1.d6 !1/cm3 -> 1/m3
    !-----------------------------------------------------------------
    if (minval(hydro70%fm1)<0.d0) then
      write(*,*) 'Info: removing negative values for mass bin'
      where(hydro70%fm1<0.d0) hydro70%fm1=0.d0
    endif
    !
    if (minval(hydro70%fn1)<0.d0) then
      write(*,*) 'Info: removing negative values for number bin'
      where(hydro70%fn1<0.d0) hydro70%fn1=0.d0
    endif
    !
    thr=conf%thr_mix_ratio(1)
    where(hydro70%qhydro (:,:,:,1) < thr)
      hydro70%qhydro (:,:,:,1) = 0.d0
      hydro70%qnhydro(:,:,:,1) = 0.d0
    endwhere
    thr=conf%thr_mix_ratio(2)
    where(hydro70%qhydro (:,:,:,2) < thr)
      hydro70%qhydro (:,:,:,2) = 0.d0
      hydro70%qnhydro(:,:,:,2) = 0.d0
    endwhere
    !
  return
  end subroutine get_hydro70_vars
 !
 ! SAM morrison by oue
  subroutine get_hydro75_vars(conf,env,mp75,hydro)
  Use wrf_var_mod
  Use crsim_mod
  Implicit None
  !
  Type(conf_var),Intent(In)                  :: conf
  Type(env_var),Intent(In)                  :: env
  Type(wrf_var_mp75),Intent(In)              :: mp75
  Type(hydro_var),Intent(InOut)              :: hydro
  !
  Integer                                    :: ix1,ix2, iy1,iy2,iz1,iz2,it
  Integer                                    :: morrID
  real*8                                     :: thr
  Integer                                    :: iz
    !
    ! default - the whole scene and it==1
    ix1=1 ; ix2=mp75%nx
    iy1=1 ; iy2=mp75%ny
    iz1=1 ; iz2=mp75%nz
    it = 1
    !
    ! reconstruct scene from configuration parameters (defined by user)
    If (hydro%nx/=mp75%nx) Then
      ix1=conf%ix_start ; ix2=conf%ix_end
    EndIf
    !
    If (hydro%ny/=mp75%ny) Then
      iy1=conf%iy_start ; iy2=conf%iy_end
    EndIf
    !
    If (hydro%nz/=mp75%nz) Then
      iz1=conf%iz_start ; iz2=conf%iz_end
    EndIf
    !
    If (mp75%nt>1) it=conf%it
    !
    ! change unit [g/kg] -> [kg/kg]
    hydro%qhydro(:,:,:,1)=mp75%qcloud(ix1:ix2,iy1:iy2,iz1:iz2,it) * 0.001
    hydro%qhydro(:,:,:,2)=mp75%qrain(ix1:ix2,iy1:iy2,iz1:iz2,it) * 0.001
    hydro%qhydro(:,:,:,3)=mp75%qice(ix1:ix2,iy1:iy2,iz1:iz2,it) * 0.001
    hydro%qhydro(:,:,:,4)=mp75%qsnow(ix1:ix2,iy1:iy2,iz1:iz2,it) * 0.001
    hydro%qhydro(:,:,:,5)=mp75%qgraup(ix1:ix2,iy1:iy2,iz1:iz2,it) * 0.001
    !
    ! change unit [/cm^3] -> [/kg]
    hydro%qnhydro(:,:,:,1)=mp75%qncloud(ix1:ix2,iy1:iy2,iz1:iz2,it) / env%rho_d(:,:,:) * 1.d6
    hydro%qnhydro(:,:,:,2)=mp75%qnrain(ix1:ix2,iy1:iy2,iz1:iz2,it) / env%rho_d(:,:,:) * 1.d6
    hydro%qnhydro(:,:,:,3)=mp75%qnice(ix1:ix2,iy1:iy2,iz1:iz2,it) / env%rho_d(:,:,:) * 1.d6
    hydro%qnhydro(:,:,:,4)=mp75%qnsnow(ix1:ix2,iy1:iy2,iz1:iz2,it) / env%rho_d(:,:,:) * 1.d6
    hydro%qnhydro(:,:,:,5)=mp75%qngraup(ix1:ix2,iy1:iy2,iz1:iz2,it) / env%rho_d(:,:,:) * 1.d6
    !
    !!
    write(*,*) 'Info: removing values below the given threshold'
    thr=1.d-50
    !-----------------------------------------------------------------------------
    
    morrID=0
    if (MaxVal(hydro%qnhydro(:,:,:,1))<0.d0) morrID=1
    !  
    do iz=1,hydro%nz
      !
      thr=conf%thr_mix_ratio(1)
      IF (morrID==0) then
        where(hydro%qhydro (:,:,iz,1) < thr)
          hydro%qhydro (:,:,iz,1) = 0.d0
          hydro%qnhydro(:,:,iz,1) = 0.d0
        endwhere
        where(hydro%qnhydro(:,:,iz,1) <= 0.d0)
          hydro%qhydro (:,:,iz,1) =0.d0
          hydro%qnhydro(:,:,iz,1) =0.d0
        endwhere
      ELSE
        where(hydro%qhydro (:,:,iz,1) < thr)
          hydro%qhydro (:,:,iz,1) = 0.d0
        endwhere
      ENDIF
      !
      !!----------------------------------
      thr=conf%thr_mix_ratio(2)
      where(hydro%qhydro (:,:,iz,2) < thr)
        hydro%qhydro (:,:,iz,2) = 0.d0
        hydro%qnhydro(:,:,iz,2) = 0.d0
      endwhere
      where(hydro%qnhydro(:,:,iz,2) <= 0.d0)
        hydro%qhydro (:,:,iz,2) =0.d0
        hydro%qnhydro(:,:,iz,2) =0.d0
      endwhere
      !!----------------------------------
      thr=conf%thr_mix_ratio(3)
      where(hydro%qhydro (:,:,iz,3) < thr)
        hydro%qhydro (:,:,iz,3) = 0.d0
        hydro%qnhydro(:,:,iz,3) = 0.d0
      endwhere
      where(hydro%qnhydro(:,:,iz,3) <= 0.d0)
        hydro%qhydro (:,:,iz,3) = 0.d0
        hydro%qnhydro(:,:,iz,3) = 0.d0
      endwhere
      !!----------------------------------
      !!----------------------------------
      thr=conf%thr_mix_ratio(4)
      where(hydro%qhydro (:,:,iz,4) < thr)
        hydro%qhydro (:,:,iz,4) = 0.d0
        hydro%qnhydro(:,:,iz,4) = 0.d0
      endwhere
      where(hydro%qnhydro(:,:,iz,4) <=0.d0)
        hydro%qhydro (:,:,iz,4) = 0.d0
        hydro%qnhydro(:,:,iz,4) = 0.d0
      endwhere
      !!----------------------------------
      !!----------------------------------
      thr=conf%thr_mix_ratio(5)
      where(hydro%qhydro (:,:,iz,5) < thr)
        hydro%qhydro (:,:,iz,5) = 0.d0
        hydro%qnhydro(:,:,iz,5) = 0.d0
      endwhere
      where(hydro%qnhydro(:,:,iz,5) <= 0.d0)
        hydro%qhydro (:,:,iz,5) = 0.d0
        hydro%qnhydro(:,:,iz,5) = 0.d0
      endwhere
      !----------------------------------
    enddo ! iz
  return
  end subroutine  get_hydro75_vars
  !
 ! CM1 morrison by oue Apr 2020
  subroutine get_hydro80_vars(conf,env,mp80,hydro)
  Use wrf_var_mod
  Use crsim_mod
  Implicit None
  !
  Type(conf_var),Intent(In)                  :: conf
  Type(env_var),Intent(In)                  :: env
  Type(wrf_var_mp80),Intent(In)              :: mp80
  Type(hydro_var),Intent(InOut)              :: hydro
  !
  Integer                                    :: ix1,ix2, iy1,iy2,iz1,iz2,it
  Integer                                    :: morrID
  real*8                                     :: thr
  Integer                                    :: iz
    !
    ! default - the whole scene and it==1
    ix1=1 ; ix2=mp80%nx
    iy1=1 ; iy2=mp80%ny
    iz1=1 ; iz2=mp80%nz
    it = 1
    !
    ! reconstruct scene from configuration parameters (defined by user)
    If (hydro%nx/=mp80%nx) Then
      ix1=conf%ix_start ; ix2=conf%ix_end
    EndIf
    !
    If (hydro%ny/=mp80%ny) Then
      iy1=conf%iy_start ; iy2=conf%iy_end
    EndIf
    !
    If (hydro%nz/=mp80%nz) Then
      iz1=conf%iz_start ; iz2=conf%iz_end
    EndIf
    !
    If (mp80%nt>1) it=conf%it
    !
    hydro%qhydro(:,:,:,1)=mp80%qcloud(ix1:ix2,iy1:iy2,iz1:iz2,it)
    hydro%qhydro(:,:,:,2)=mp80%qrain(ix1:ix2,iy1:iy2,iz1:iz2,it)
    hydro%qhydro(:,:,:,3)=mp80%qice(ix1:ix2,iy1:iy2,iz1:iz2,it)
    hydro%qhydro(:,:,:,4)=mp80%qsnow(ix1:ix2,iy1:iy2,iz1:iz2,it)
    hydro%qhydro(:,:,:,5)=mp80%qgraup(ix1:ix2,iy1:iy2,iz1:iz2,it)
    !
    hydro%qnhydro(:,:,:,1)=mp80%qncloud(ix1:ix2,iy1:iy2,iz1:iz2,it)
    hydro%qnhydro(:,:,:,2)=mp80%qnrain(ix1:ix2,iy1:iy2,iz1:iz2,it)
    hydro%qnhydro(:,:,:,3)=mp80%qnice(ix1:ix2,iy1:iy2,iz1:iz2,it)
    hydro%qnhydro(:,:,:,4)=mp80%qnsnow(ix1:ix2,iy1:iy2,iz1:iz2,it)
    hydro%qnhydro(:,:,:,5)=mp80%qngraup(ix1:ix2,iy1:iy2,iz1:iz2,it)
    !
    !!
    write(*,*) 'Info: removing values below the given threshold'
    thr=1.d-50
    !-----------------------------------------------------------------------------
    
    morrID=0
    if (MaxVal(hydro%qnhydro(:,:,:,1))<0.d0) morrID=1
    !  
    do iz=1,hydro%nz
      !
      thr=conf%thr_mix_ratio(1)
      IF (morrID==0) then
        where(hydro%qhydro (:,:,iz,1) < thr)
          hydro%qhydro (:,:,iz,1) = 0.d0
          hydro%qnhydro(:,:,iz,1) = 0.d0
        endwhere
        where(hydro%qnhydro(:,:,iz,1) <= 0.d0)
          hydro%qhydro (:,:,iz,1) =0.d0
          hydro%qnhydro(:,:,iz,1) =0.d0
        endwhere
      ELSE
        where(hydro%qhydro (:,:,iz,1) < thr)
          hydro%qhydro (:,:,iz,1) = 0.d0
        endwhere
      ENDIF
      !
      !!----------------------------------
      thr=conf%thr_mix_ratio(2)
      where(hydro%qhydro (:,:,iz,2) < thr)
        hydro%qhydro (:,:,iz,2) = 0.d0
        hydro%qnhydro(:,:,iz,2) = 0.d0
      endwhere
      where(hydro%qnhydro(:,:,iz,2) <= 0.d0)
        hydro%qhydro (:,:,iz,2) =0.d0
        hydro%qnhydro(:,:,iz,2) =0.d0
      endwhere
      !!----------------------------------
      thr=conf%thr_mix_ratio(3)
      where(hydro%qhydro (:,:,iz,3) < thr)
        hydro%qhydro (:,:,iz,3) = 0.d0
        hydro%qnhydro(:,:,iz,3) = 0.d0
      endwhere
      where(hydro%qnhydro(:,:,iz,3) <= 0.d0)
        hydro%qhydro (:,:,iz,3) = 0.d0
        hydro%qnhydro(:,:,iz,3) = 0.d0
      endwhere
      !!----------------------------------
      !!----------------------------------
      thr=conf%thr_mix_ratio(4)
      where(hydro%qhydro (:,:,iz,4) < thr)
        hydro%qhydro (:,:,iz,4) = 0.d0
        hydro%qnhydro(:,:,iz,4) = 0.d0
      endwhere
      where(hydro%qnhydro(:,:,iz,4) <=0.d0)
        hydro%qhydro (:,:,iz,4) = 0.d0
        hydro%qnhydro(:,:,iz,4) = 0.d0
      endwhere
      !!----------------------------------
      !!----------------------------------
      thr=conf%thr_mix_ratio(5)
      where(hydro%qhydro (:,:,iz,5) < thr)
        hydro%qhydro (:,:,iz,5) = 0.d0
        hydro%qnhydro(:,:,iz,5) = 0.d0
      endwhere
      where(hydro%qnhydro(:,:,iz,5) <= 0.d0)
        hydro%qhydro (:,:,iz,5) = 0.d0
        hydro%qnhydro(:,:,iz,5) = 0.d0
      endwhere
      !----------------------------------
    enddo ! iz
  return
  end subroutine  get_hydro80_vars
  !
  !
  subroutine determine_elevation_and_range(x,y,z,ixc,iyc,zc,dxc,dyc,elev,rr)
  use crsim_mod
  use phys_param_mod, ONLY: r2d
  Implicit None
  !
  Real*8, Intent(In)       :: x,y,z    ! [m]
  Integer,Intent(In)       :: ixc,iyc  ! indices of xc,yc
  Real*8, Intent(In)       :: zc       ! [m]
  Real*8, Intent(In)       :: dxc,dyc  ! [m] resolution in x and y (constant)
  Real*8,Intent(Out)       :: elev     ! [degree]
  Real*8,Intent(Out)       :: rr       ! [m] radar range gate
  !
  real*8                   :: xc,yc
  real*8                   :: d_x,d_y,d_z
    !
    xc=dble(ixc-1)*dxc 
    yc=dble(iyc-1)*dyc 
    !
    d_x=x-xc
    d_y=y-yc
    d_z=z-zc
    !
    ! AUG2019
    !-- added for airborne radar by oue 20190911
    if (ixc == -999)  d_x = 0.d0
    if (iyc == -999)  d_y = 0.d0
    !---
    !
    rr=dsqrt(d_x*d_x + d_y*d_y)

    if (rr>0.d0) then
      elev=datan(d_z/rr) * r2d ! the sign depends on d_z>0 or d_z<0
    else
       if (d_z>=0.d0) then  ! negative elevation by oue 20190911
          elev=90.d0
       else ! negative elevation by oue 20190911
          elev=-90.d0 ! negative elevation by oue 20190911
       end if ! this is the airborne case  AT ! AUG2019
    endif
    !
    ! output radar range gate
    rr=dsqrt(d_x*d_x + d_y*d_y + d_z*d_z)
    !
  return
  end subroutine determine_elevation_and_range
  !
  subroutine determine_range(x,y,z,ixc,iyc,zc,dxc,dyc,rr)
  use crsim_mod
  use phys_param_mod, ONLY: r2d
  Implicit None
  !
  Real*8, Intent(In)       :: x,y,z    ! [m]
  Integer,Intent(In)       :: ixc,iyc  ! indices of xc,yc
  Real*8, Intent(In)       :: zc       ! [m]
  Real*8, Intent(In)       :: dxc,dyc  ! [m] resolution in x and y (constant)
  Real*8,Intent(Out)       :: rr       ! [m] radar range gate
  !
  real*8                   :: xc,yc
  real*8                   :: d_x,d_y,d_z
    !
    xc=dble(ixc-1)*dxc
    yc=dble(iyc-1)*dyc
    !
    d_x=x-xc
    d_y=y-yc
    d_z=z-zc
    !
    ! AUG2019
    !-- added for airborne radar by oue 20190911  
    if (ixc == -999)  d_x = 0.d0
    if (iyc == -999)  d_y = 0.d0
    !---
    ! output radar range gate
    rr=dsqrt(d_x*d_x + d_y*d_y + d_z*d_z)
    !
  return
  end subroutine determine_range
  !
  subroutine determine_azimuth(x,y,ixc,iyc,dxc,dyc,elev,azim)
  use crsim_mod
  use phys_param_mod, ONLY: r2d
  Implicit None
  !
  Real*8, Intent(In)       :: x,y    ! [m]
  Integer,Intent(In)       :: ixc,iyc  ! indices of xc,yc
  Real*8, Intent(In)       :: dxc,dyc  ! [m] resolution in x and y (constant)
  Real*8, Intent(In)       :: elev     ! [degree]  elevation
  Real*8,Intent(Out)       :: azim     ! [degree]  E=0, W=180; N=90; S=270
  !
  real*8                   :: xc,yc
  real*8                   :: d_x,d_y
    !
    ! 
  !if (elev==90.d0) then
  if ((elev==90.d0) .or. (elev==-90.d0)) then ! changed elev to abs(elev) by oue 20190911
      azim=0.d0
      go to 123
    endif

! AUG2019
!    AT
!    !if not airborne
!    If ( (ixc /= -999) .or. (iyc /= -999) ) then
!      if ( (elev >= 89.5d0) .and. (elev <= 89.5d0) ) then 
!        azim=0.d0
!        go to 123
!      endif
!    Else ! if airborne  
!      if ( (elev >= -0.5d0 ) .and. (elev <= 0.5d0) ) then
!        azim=0.d0
!        go to 123
!      endif
!    Endif
    !
    xc=dble(ixc-1)*dxc
    yc=dble(iyc-1)*dyc
    !
    d_x=x-xc
    d_y=y-yc
    !
! AUG2019
!AT -- added for airborne radar by oue 20190911
    !if (ixc == -999)  d_x = 0.d0
    !if (iyc == -999)  d_y = 0.d0
    !
    !if (ixc /= -999) then
    if (ixc == -999) then
      if (d_y < 0.d0) azim=270.d0
      if (d_y >= 0.d0) azim=90.d0
      go to 123
    end if
    !
    !if (iyc /= -999) then
    if (iyc == -999) then
      if (d_x < 0.d0) azim=180.d0
      if (d_x >= 0.d0) azim=0.d0
      go to 123
    end if
    !
!------
    ! 
    if (d_x==0.d0) then
      if (d_y>=0) azim=90.d0
      if (d_y<0)  azim=270.d0
      go to 123
    endif
    !
    if (d_y==0.d0) then
      if (d_x>=0) azim=0.d0
      if (d_x<0)  azim=180.d0
      go to 123
    endif
    !
    azim=datan(dabs(d_y)/dabs(d_x)) * r2d
    !
    if ( (d_x<0.d0) .and. (d_y>0.d0) ) azim=180.d0-azim  ! II quadr.
    if ( (d_x<0.d0) .and. (d_y<0.d0) ) azim=180.d0+azim  ! III quadr.
    if ( (d_x>0.d0) .and. (d_y<0.d0) ) azim=360.d0-azim  ! IV quadr.
    !
    123 continue
    
    if (azim<0.d0) then 
      write(*,*)  'problem determne_azimuth, azim=',azim
      stop
    endif
    !
    if (azim>360.d0) then
      write(*,*)  'problem determne_azimuth, azim=',azim
      stop
    endif
    !
  return
  end subroutine determine_azimuth
  !
  subroutine determine_sw_contrib_terms(rr,elev,azim,uu,vv,ww,tke,sigma_r,sigma_theta, Kx,Ky,Kz,sw_t2,sw_s2,sw_v2)
  use crsim_mod
  use phys_param_mod, ONLY: d2r,twoOverThree
  Implicit None
  
  Real*8, Intent(In)       :: rr       ! [m] radar range gate
  Real*8, Intent(In)       :: elev     ! [deg]  elevatiom
  Real*8, Intent(In)       :: azim     ! [deg]  azimuth
  Real*8, Intent(In)       :: uu,vv,ww ! [m/s] wind components in cartesian coordinates
  Real*8, Intent(In)       :: tke      ! [m^2/s^2] turbulence kinetic energy or eddy dissipation rate 
  Real*8, Intent(In)       :: sigma_r 
  Real*8, Intent(In)       :: sigma_theta
  Real*8, Intent(In)       :: Kx,Ky,Kz ! [1/s] wind shear in cartesian coordinates
  Real*8, Intent(Out)      :: sw_t2    ! [m^2/s^2] spectrum width due to turbulence
  Real*8, Intent(Out)      :: sw_s2    ! [m^2/s^2] spectrum width due to wind shear in radar volume
  Real*8, Intent(Out)      :: sw_v2    ! [m^2/s^2] spectrum width due to cross wind
  !
  ! reference:
  ! Zhang et al,2009:Estimate of Eddy Dissipation Rate Using Spectrum Width Observed by the Hong
  ! Kong TDWR Radar,34th Conference on Radar Meteorology. P6.9. 
  !
  Real*8,Parameter  :: aconst=1.6d0 ! a universal dimensionless constant A  between 1.53 and 1.68 in experssion : E(K)= A (tke)^2/3 K (-5/3) 
  real*8            :: aa,cconst,dd,expr
  real*8            :: Kr,Kelev,Kazim !  [1/s]  wind shear terms in radar coordinates
  real*8            :: vcross  
     
    aa=aconst**1.5d0
    cconst=1.d0/0.72
    !
    dd=sigma_r/(rr*sigma_theta) 
    !
    ! turbulence term
    !
    sw_t2=0.d0
    !
    if (tke>0.d0) then
      if (dd<=1.d0) then
        sw_t2 = (tke * rr * sigma_theta * aa * cconst)**twoOverThree
      else
        expr= 11.d0/15.d0 + 4.d0/15.d0  *  (rr*rr  *  (sigma_theta*sigma_theta)  /  (sigma_r*sigma_r) )
        sw_t2 = sigma_r**twoOverThree * 1.35d0 * aconst * tke**twoOverThree * expr
      endif
    endif
    !
    ! wind shear term
    !
    Kr    =    Kx*dcos(azim*d2r)*dcos(elev*d2r)     + Ky*dsin(azim*d2r)*dcos(elev*d2r)    + Kz*dsin(elev*d2r)
    Kelev =    Kx*dcos(azim*d2r)*dsin(elev*d2r)     + Ky*dsin(azim*d2r)*dsin(elev*d2r)    - Kz*dcos(elev*d2r)
    !Kazim =  - Kx*dsin(azim*d2r)*dcos(elev*d2r)     + Ky*dcos(azim*d2r)*dcos(elev*d2r)    
    Kazim =  - Kx*dsin(azim*d2r)                    + Ky*dsin(azim*d2r)
    !
    sw_s2= (sigma_r*Kr)*(sigma_r*Kr)  +   (sigma_theta*Kelev)*(sigma_theta*Kelev)   +  (sigma_theta*Kazim)*(sigma_theta*Kazim)
    ! 
    ! cross wind term
    !
    vcross= uu*dcos(azim*d2r)*dsin(elev*d2r) + vv*dsin(azim*d2r)*dsin(elev*d2r) + ww*dcos(elev*d2r)
    sw_v2= (vcross * sigma_theta) * (vcross * sigma_theta)
    !
  return
  end subroutine determine_sw_contrib_terms
  !
  subroutine compute_attenuated_lidar_signal(nz,zz,back_true,ext,back_obs)
  use phys_param_mod, ONLY: m999
  ! Lidar assumed to be at the first grid level
  !-modifiedd to work either height increment or decrement, by oue for ICON
  Implicit None
  !!
  integer, Intent(in)   :: nz
  Real*8, Intent(In)    :: zz(nz) ! heigh [m]
  Real*8, Intent(In)    :: back_true(nz)  ! true lidar backscatter [m-1 sr-1] 
  Real*8, Intent(In)    :: ext(nz)        ! lidar extinction [m-1 sr-1] 
  Real*8, Intent(OuT)   :: back_obs(nz)   ! attenuated lidar backscatter [m-1 sr-1]
  !
  Integer                 :: iz
  Real*8                  :: loc_tau,tot_tau
  Integer                 :: izstart,izmin,izmax,idz
    !
    tot_tau=0.d0
    loc_tau=0.d0
    !
    !-Check if height is increment or decrement added by oue for ICON
    izstart=1 ; izmin=2 ; izmax=nz ; idz= 1
    if(zz(1) > zz(nz)) then
      izstart=nz ; izmin=nz-1 ; izmax=1 ; idz = -1
    endif
    !
    if (back_true(izstart)>0.d0) then
      loc_tau=Max(0.d0,ext(izstart)) * zz(izstart) ! Np/m *m =Np
      tot_tau=tot_tau +loc_tau
      back_obs(izstart)=back_true(izstart) * dExp(-2.d0*tot_tau)
    endif
    !
    do iz=izmin,izmax,idz
      if (back_true(iz)>0.d0) then
        loc_tau=Max(0.d0,ext(iz)) * (zz(iz)-zz(iz-idz)) ! Np/m *m =Np
        tot_tau=tot_tau +loc_tau
        back_obs(iz)=back_true(iz) * dExp(-2.d0*tot_tau)
      endif
    enddo
    !
  return
  end subroutine compute_attenuated_lidar_signal
  !
  subroutine estimate_ceilo_first_cloud_base(nz,zz,back_obs,first_cloud_base)
  !-modifiedd to work either height increment or decrement, by oue for ICON
  Implicit None
  !
  integer, Intent(In)          :: nz
  Real*8,Intent(In)            :: zz(nz) ! height [m]
  Real*8,Intent(In)            :: back_obs(nz) ! ceilo observed backscatter [m sr]^-1
  Real*8,Intent(Out)           :: first_cloud_base  ! height of the first cloud base [m]
  !
  integer                      :: iz
  integer                      :: izb,izt  ! vertical indices of the begining and the end of the first cloud layer
  integer                      :: iz1,iz2,jj
  Real*8                       :: max_back
  integer                      :: min_num_cl_pix ! min number of cloud pixels in vertical to be considered as cloud 
  Integer                 :: izstart,izmin,izmax,idz
  
    !-Check if height is increment or decrement added by oue for ICON
    izstart=1 ; izmin=2 ; izmax=nz ; idz= 1
    if(zz(1) > zz(nz)) then
      izstart=nz ; izmin=nz-1 ; izmax=1 ; idz = -1
    endif
    !
    min_num_cl_pix=2
    first_cloud_base=-999.d0
    !
    iz1=1 ; iz2=nz
    iz1=izmin ; iz2=izmax
    ! 
    112 continue
    !   
    izb=-1; izt=-1
    !     
    A1: do iz =iz1, iz2, idz
      if (back_obs(iz)>0.d0) then
        izb=iz
        exit A1
      endif
    enddo A1
    ! 
    if (izb==izmax) go to 111
    if (izb<0)   go to 111
    ! 
    A2: do iz=iz1+idz,iz2,idz
      if (back_obs(iz)<0.d0) then
        izt=iz-idz
        exit A2
      endif
    enddo A2
    ! 
    if ( (izb>0) .and. (izt<0) .and. (back_obs(izmax) > 0.d0) )  izt=izmax
    ! 
    if(zz(1) < zz(nz)) then
      if(izt-izb+1<min_num_cl_pix) then
        iz1=izt+1 ; iz2=nz  
        go to 112 
      endif
    else
      if(izb-izt+1<min_num_cl_pix) then
        iz1=izt+idz ; iz2=izmax
        go to 112 
      endif
    endif
    !--------------------------------------
    if ((izb>0) .and. izt>0) Then
      if(izb < izt) then
        max_back=MAxVal(back_obs(izb:izt))
      else
        max_back=MAxVal(back_obs(izt:izb))
      endif
      jj=-1
      A3: do iz=izb,izt,idz
        if (back_obs(iz)==max_back) then
          jj=iz
          exit A3
        endif
      enddo A3
    
      ! if found max is at cloud top (and this is possible if cloud consists of 2 pixels in
      ! vertical), don't consider this layer 
      if(zz(1) <= zz(nz)) then
        if( (jj>0) .and. (jj >= izt)) then
          iz1=izt+1 ; iz2=nz
          go to 112
        endif
      else
        if( (jj>0) .and. (jj <= izt)) then
          iz1=izt+idz ; iz2=izmax
          go to 112
        endif
      endif
    
      if ((jj>0).and.(jj<=nz)) first_cloud_base=zz(jj)  ! to output indices: dble(jj) 
      if ((jj<0).or.(jj>nz)) then
        write(*,*) 'problem subroutine estimate_ceilo_first_cloud_base, jj',jj
        stop
      endif
    !
    endif ! if ((iz1>0) .and. iz2>0) Then
    !--------------------------------------
    !
    111 continue
    ! 
  return
  end subroutine estimate_ceilo_first_cloud_base
  !
  subroutine spinhirne_aero_model(conf,env,aero)
  Use crsim_mod
  Implicit None
  !
  Type(conf_var),Intent(In)          :: conf
  Type(env_var),Intent(In)           :: env
  Type(aero_var),Intent(InOut)       :: aero
  !
  Real*8,parameter      :: sigma_0 = 0.025d0 ! km^-1
  Real*8,parameter      :: a=0.4d0, b= 1.6d0 ! km
  Real*8,parameter      :: ap= 0.2981d0, bp=2.5d0 ! km
  Real*8,parameter      :: f=1.5d-7 ! km^-1
  !
  integer               :: ix,iy
  Real*8,Dimension(:),Allocatable  :: aero_back_obs,ext_new
    !
    ! Ref: Spinhirne D. James, 1993:  Micro Pulse Lidar.IEEE Trans. Geo. Rem. Sens.,31, 48-55.
    !
    !aerosol lidar ratio (sr)
    aero%lidar_ratio=conf%aero_lidar_ratio
    !
    ! extinction profile ! km^-1
    aero%ext =  sigma_0 * ((1+a)*(1+a)) * (  dexp(env%z/(b*1.d+3)) / &
                (  (a + dexp(env%z/(b*1.d+3))) *  (a + dexp(env%z/(b*1.d+3))) ) ) 
    aero%ext =  aero%ext + f*((1+ap)*(1+ap)) * (  dexp(env%z/(bp*1d+3)) / &
                ( (a+dexp(env%z/(bp*1d+3))) * (a+dexp(env%z/(bp*1d+3))) ) )
    !
    !  add boundary layer aerosols
    where (env%z<=1000.d0) 
      aero%ext = aero%ext + 0.05d0 ! km^-1
    end where
    !
    ! units transformation km^-1 ->  m^-1
    aero%ext = (aero%ext*1.d-3) ! m^-1
    !
    !-----------------------
    Allocate(aero_back_obs(aero%nz))
    !
    do ix=1,env%nx
      do iy=1,env%ny
    
        ! normalization of the aerosol extinction profile m^-1
        !
        if (conf%aero_tau>0.d0) then
          Allocate(ext_new(aero%nz))
          ext_new=0.d0
          call normalize_spinhirne_aero_model(conf%aero_tau,env%nz,env%z(ix,iy,:),aero%ext(ix,iy,:),ext_new)
          aero%ext(ix,iy,:)=ext_new(:)
          Deallocate(ext_new)
        endif
        !
        ! aerosol true backscatter aero%back_true
        aero%back_true(ix,iy,:)=aero%ext(ix,iy,:)/aero%lidar_ratio(ix,iy,:) !  m-1 sr-1
        !
        ! attenuated aerosol backscatter aero%back_obs
        aero_back_obs=0.d0   
        call compute_attenuated_lidar_signal(env%nz,env%z(ix,iy,:),aero%back_true(ix,iy,:),aero%ext(ix,iy,:),aero_back_obs)
        aero%back_obs(ix,iy,:)=aero_back_obs
        !
      enddo ! iy
    enddo ! ix
    !
    Deallocate(aero_back_obs)
    !
  return
  end subroutine spinhirne_aero_model
  !
  subroutine normalize_spinhirne_aero_model(aero_tau,nz,zz,ext,ext_new)
  !-modifiedd to work either height increment or decrement, by oue for ICON
  Implicit None
  !
  Real*8,Intent(In)     :: aero_tau
  integer, Intent(in)   :: nz
  Real*8, Intent(In)    :: zz(nz) ! heigh [m]
  Real*8, Intent(In)    :: ext(nz)        ! lidar extinction [m-1] 
  Real*8, Intent(OuT)   :: ext_new(nz)    ! normalized idar extinction profile [m-1] 
  !
  Integer                 :: iz
  Real*8                  :: loc_tau,tot_tau,norm_coeff
  Integer                 :: izstart,izmin,izmax,idz
  
    !-Check if height is increment or decrement added by oue for ICON
    izstart=1 ; izmin=2 ; izmax=nz ; idz= 1
    if(zz(1) > zz(nz)) then
      izstart=nz ; izmin=nz-1 ; izmax=1 ; idz = -1
    endif
    !
    ext_new=0.d0
    tot_tau=0.d0
    loc_tau=0.d0
    !
    if (ext(izstart)>0.d0) then
      loc_tau=Max(0.d0,ext(izstart)) * zz(izstart) ! Np/m *m =Np
      tot_tau=tot_tau +loc_tau
    endif
    !
    do iz=izmin,izmax,idz
      if (ext(iz)>0.d0) then
        loc_tau=Max(0.d0,ext(iz)) * (zz(iz)-zz(iz-idz)) ! Np/m *m =Np
        tot_tau=tot_tau +loc_tau
      endif
    enddo
    !
    norm_coeff=aero_tau/tot_tau
    ext_new=ext*norm_coeff
    !
  return
  end subroutine normalize_spinhirne_aero_model
  !
  subroutine get_path(long_name,path)
  Implicit none
  !
  !Returns the path given a full filename (path+filename)
  !
  character(len=*),intent(in)  :: long_name
  character(len=*),intent(out) :: path
  !
  integer ibegstr
  integer iendstr
  integer i
  integer n
  logical endstrfound,ispath
    !
    n=len(long_name)
    !
    ! find the actual length of the string
    i=0
    endstrfound=.false.
    do while(.not.endstrfound)
      if(long_name(n-i:n).ne.' ') then
        iendstr=n-i+1
        endstrfound=.true.
      endif
      i=i+1
    enddo
    !
    ! Is there a '/' ?
    !
    ispath=.false.
    do i=1,min(iendstr,n)
      if (long_name(i:i)=='/') then
        ispath=.true.
      endif
    enddo
    !
    if (ispath) then
      ! find the last '/'
      i=min(iendstr,n)
      do while(long_name(i:i).ne.'/')
        i=i-1
      enddo
      ibegstr=i+1
      !
      path=long_name(1:ibegstr-1)
      !
    else
      path='./'
    endif
    !
  return
  end subroutine get_path
  !
  subroutine get_filename(long_name,filename)
  Implicit none
  !
  !Returns the filename given a full filename (path+filename)
  !
  character(len=*),intent(in)  :: long_name
  character(len=*),intent(out) :: filename
  !
  integer ibegstr
  integer iendstr
  integer i
  integer n
  logical endstrfound
  logical ispath
    !
    n=len(long_name)
    !
    ! find the actual length of the string
    i=0
    endstrfound=.false.
    do while(.not.endstrfound)
      if(long_name(n-i:n).ne.' ') then
        iendstr=n-i+1
        endstrfound=.true.
      endif
      i=i+1
    enddo
    !
    iendstr=min(iendstr,n)
    !
    ! Is there a '/' ?
    ispath=.false.
    do i=1,min(iendstr,n)
      if (long_name(i:i)=='/') then
        ispath=.true.
      endif
    enddo
    !
    if (ispath) then
      ! find the last '/'
      i=iendstr
      do while(long_name(i:i).ne.'/')
        i=i-1
      enddo
      ibegstr=min(i+1,n)
      !
      !path=long_name(1:ibegstr-1)
      filename=long_name(ibegstr:iendstr)
    else
      filename=long_name(1:iendstr)
    endif
    !
  return
  end subroutine get_filename
  !
  subroutine get_extention(long_name,extention)
  !
  Implicit none
  ! Returns the extention (dsfg.dfgdfg.dat ==> dat)
  ! extention=' ' if there is no extention
  !
  character(len=*),intent(in)  :: long_name
  character(len=*),intent(out) :: extention
  !
  integer ibegstr
  integer iendstr
  integer i
  integer n
  logical endstrfound
  logical ispath
    !
    n=len(long_name)
    !
    ! find the actual length of the string
    i=0
    endstrfound=.false.
    do while(.not.endstrfound)
      if(long_name(n-i:n).ne.' ') then
        iendstr=n-i+1
        endstrfound=.true.
      endif
    i=i+1
    enddo
    !
    iendstr=min(iendstr,n)
    !
    ! Is there a '.' ?
    ispath=.false.
    do i=1,min(iendstr,n)
      if (long_name(i:i)=='.') then
        ispath=.true.
      endif
    enddo
    !
    if (ispath) then
      ! find the last '.'
      i=iendstr
      do while(long_name(i:i).ne.'.')
        i=i-1
      enddo
      ibegstr=min(i+1,n)
      !  path=long_name(1:ibegstr-1)
      extention=long_name(ibegstr:iendstr)
    else
      extention=' '
    endif
    !
  return
  end subroutine get_extention
  !
  subroutine get_lut_str_var(n,lvar,cvar,var)
  Implicit None
  Integer,Intent(in)   :: n
  Real*8,Intent(in)    :: lvar(n)
  Real*8,Intent(in)    :: cvar
  Real*8,Intent(out)   :: var
  !
  integer              :: ii
  real*8               :: diff0,diff
  integer              :: i
  
    ii=-1
    diff0=1.d+23
    !
    do i=1,n
      diff=dabs(lvar(i)-cvar)
      if (diff<diff0) then
        diff0=diff
        ii=i
      endif
    enddo
    !
    var=lvar(ii)
    !
  return
  end subroutine get_lut_str_var
  !
  subroutine get_lut_str(fff,ttt,elev,frq_str,t_str,el_str)
  Implicit None
  Real*8,Intent(In)          :: fff  ! freq GHz
  Real*8,Intent(In)          :: ttt  ! temp K
  Real*8,Intent(In)          :: elev ! elev deg
  Character,Intent(Out)      :: frq_str*4,t_str*5,el_str*4
  !
  character                  :: frq_str2*10,el_str2*10
    !
    if  (fff<10.d0) then
      write(frq_str2,'(f4.1)') fff
      frq_str='0'//Adjustl(Trim(frq_str2))
    else
      write(frq_str2,'(f5.1)') fff
      frq_str=Adjustl(Trim(frq_str2))
    endif
    !
    write(t_str,'(f5.1)') ttt ! K
    !
    if  (elev<10.d0) then
      write(el_str2,'(f4.1)') elev
      el_str='0'//Adjustl(Trim(el_str2))
    else
      write(el_str2,'(f5.1)') elev
      el_str=Adjustl(Trim(el_str2))
    endif
    !
  return
  end subroutine get_lut_str
  !
  subroutine get_lut_den_str(den,rho_str)
  Implicit None
  Real*8,Intent(In)          :: den ! kg/m^3 hydrometeor density 
  Character,Intent(Out)      :: rho_str*5
  !
  character                  :: rho_str2*10
    !
    if  (den<1.d0) then
      write(rho_str2,'(f7.4)') den
      rho_str=Adjustl(Trim(rho_str2))
    else if ( (den>=1.d0) .and.(den<10.d0))  then
      write(rho_str2,'(f6.3)') den
      rho_str=Adjustl(Trim(rho_str2))
    else if ( (den>=10.d0) .and.(den<100.d0))  then
      write(rho_str2,'(f6.2)') den
      rho_str=Adjustl(Trim(rho_str2))
    else
      write(rho_str2,'(f7.1)') den
      rho_str=Adjustl(Trim(rho_str2))
    endif
    !
  return
  end subroutine get_lut_den_str
  !
  subroutine linpol ( x, y, imax, xp, yp, jmax )
  Implicit None
  !     
  ! Does a linear interpolation on the data given by the arrays
  ! x(1..imax) (x-coordinate) and
  ! y(1..imax) (y-coordinate)
  !     
  ! Data are interpolated at the x-values given by the array
  !xp(1..jmax) (x-coordinates of the interpolated values)
  !     
  ! the interpolated y-values are returned in the array  yp(1..jmax)
  ! 
  integer, Intent(In)       :: imax
  real*8, Intent(In)        :: x(imax),y(imax)
  integer, Intent(In)       :: jmax
  real*8, Intent(In)        :: xp(jmax)
  real*8, Intent(Out)       ::yp(jmax)
  integer                   :: i,j,imin
    !     
    imin=1
    !
    do 1 j = 1, jmax
      if (xp(j) .lt. x(1)) then
        i=1
        yp(j)=(y(i+1)-y(i))/(x(i+1)-x(i))*(xp(j)-x(i)) + y(i)
      else if (xp(j) .gt. x(imax)) then
        i=imax-1
        yp(j)=(y(i+1)-y(i))/(x(i+1)-x(i))*(xp(j)-x(i)) + y(i)
      else if (xp(j) .eq. x(imax)) then
        yp(j)=y(imax)
      else
        do 2 i = imin, imax-1
          if (xp(j) .eq. x(i)) then
            imin=i
            yp(j)=y(i)
            goto 1
          else if ((xp(j) .gt. x(i)).and.(xp(j) .le. x(i+1))) then
            imin=i
            yp(j)=(y(i+1)-y(i))/(x(i+1)-x(i))*(xp(j)-x(i)) + y(i)
            goto 1
          endif
        2 continue
      endif
    1 continue
    !     
  return
  end subroutine  linpol
  !   
  !--------------------------------------------------------------------------------
  ! SUBROUTINE access_lookup_table is added by DW 2017/10/30 for P3
  SUBROUTINE access_lookup_table(dumjj,dumii,dumi,index,densize,rimsize,isize,tabsize,itab,proc)
  !
  implicit none
  Integer,Intent(In)    :: dumjj,dumii,dumi,index,densize,rimsize,isize,tabsize
  Real*8,Intent(Out)    :: proc
  real*8,Intent(In),dimension(densize,rimsize,isize,tabsize) :: itab
    !
    proc = itab(dumjj,dumii,dumi,index)
    !
  END SUBROUTINE access_lookup_table
  !
  !------------------------------------------------------------------------------
  !subroutine find_lookupTable_indices_1 is added by DW 2017/10/30 for P3
  subroutine find_lookupTable_indices_1(dumi,dumjj,dumii,dum1,dum4,dum5,&
                            isize,rimsize,densize,qitot,nitot,qirim,rhop)
  !
  ! Finds indices in lookup table
  !
  implicit none
  !
  integer, intent(out) :: dumi,dumjj,dumii
  real*8,  intent(out) :: dum1,dum4,dum5
  integer, intent(in)  :: isize,rimsize,densize
  real*8,  intent(in)  :: qitot,nitot,qirim,rhop
  
    ! find index for qi (normalized ice mass mixing ratio = qitot/nitot)
    dum1 = (dlog10(qitot/nitot)+18.d0)/(0.1d0*dlog10(261.7d0))-10.d0
    dumi = int(dum1)
    !
    ! set limits (to make sure the calculated index doesn't exceed
    ! range of lookup table)
    dum1 = min(dum1,dble(isize))
    dum1 = max(dum1,1.d0)
    dumi = max(1,dumi)
    dumi = min(isize-1,dumi)
    !
    ! find index for rime mass fraction
    dum4  = (qirim/qitot)*3.d0 + 1.d0
    dumii = int(dum4)
    !
    ! set limits
    dum4  = min(dum4,dble(rimsize))
    dum4  = max(dum4,1.d0)
    dumii = max(1,dumii)
    dumii = min(rimsize-1,dumii)
    !
    ! find index for bulk rime density
    ! (account for uneven spacing in lookup table for density)
    if (rhop.le.650.d0) then
      dum5 = (rhop-50.d0)*0.005d0 + 1.d0
    else
      dum5 = (rhop-650.d0)*0.004d0 + 4.d0
    endif
    dumjj = int(dum5)
    ! set limits
    dum5  = min(dum5,dble(densize))
    dum5  = max(dum5,1.d0)
    dumjj = max(1,dumjj)
    dumjj = min(densize-1,dumjj)
    !
  end subroutine find_lookupTable_indices_1
  !
  !-----------------------------------------------------------
  ! Real function gammq(a,x) is added by DW 2017/10/30 for P3
  Double  precision function gammq(a,x) 
  IMPLICIT NONE
  Real*8,Intent(in)     :: a,x
  ! USES gcf,gser
  ! Returns the incomplete gamma function Q(a,x) = 1-P(a,x)
  Real*8                ::  gammcf,gamser
  
    !if (x.lt.0..or.a.le.0)
    if (x.lt.a+1.d0) then
      call gser(gamser,a,x)
      gammq=1.d0-gamser
    else
      call gcf(gammcf,a,x)
      gammq=gammcf
    end if
    !
  !return
  end function gammq
  !
  !-------------------------------------------------------------
  ! subroutine gser(gamser,a,x) is added by DW 2017/10/30 for P3
  subroutine gser(gamser,a,x)       
  Real*8,Intent(In)        :: a,x
  Real*8,Intent(Out)       :: gamser
  Integer,parameter :: itmax=100
  Real*8,parameter  :: eps=3.d-7
  Integer           :: n
  Real*8            :: ap,del,gln,summ
    !
    gln = dlog(dgamma(a))
    !
    if (x.le.0.d0) then
      gamser = 0.d0
      return
    end if
    !
    ap=a
    summ=1.d0/a
    del=summ
      do n=1,itmax
        ap=ap+1.d0
        del=del*x/ap
        summ=summ+del
        if (dabs(del).lt.dabs(summ)*eps) goto 1
      end do
    1 gamser=summ*dexp(-x+a*dlog(x)-gln)
  return
  end
  !------------------------------------------------------------------
  ! subroutine gcf(gammcf,a,x) is added by DW 2017/10/30 for P3
  subroutine gcf(gammcf,a,x)       
  real*8,Intent(In)        :: a,x
  real*8,Intent(Out)       :: gammcf
  integer,parameter :: itmax=100
  real*8,parameter  :: eps=3.d-7
  real*8,parameter  :: fpmin=1.d-30
  integer           :: i
  real*8            :: an,b,c,d,del,h,gln
    !
    gln=dlog(dgamma(a))
    b=x+1.d0-a
    c=1.d0/fpmin
    d=1.d0/b
    h=d
    do i=1,itmax
      an=-i*(i-a)
      b=b+2.d0
      d=an*d+b
      if(dabs(d).lt.fpmin) d=fpmin
      c=b+an/c
      if(dabs(c).lt.fpmin) c=fpmin
      d=1.d0/d
      del=d*c
      h = h*del
      if(dabs(del-1.d0).lt.eps)goto 1
    end do
    1 gammcf=dexp(-x+a*dlog(x)-gln)*h
    !
  return
  end
  !-------------------------------------
  subroutine get_arg_sim_input(instring,delim,iwant,outstring,status)
  Implicit None
  !
  ! Scans through instring and returns the 
  ! substring seperated by occurances of delim 
  ! specified by iwant
  ! status =0 ==> success
  ! status =1 ==> iwant is out of range 
  !
  character(len=*),intent(in)    :: instring
  character(len=1),intent(in)    :: delim
  integer,intent(in)             :: iwant
  !
  character(len=*)               :: outstring
  integer,intent(out)            :: status
  !
  integer                          :: i,nmax,n
  integer,dimension(:),allocatable :: delim_pos
  character(len=1000)              :: temp_string
    !
    status=0
    !
    !-----------------------------
    ! How many delims are there ?
    !-----------------------------
    !
    nmax=len(trim(adjustl(instring)))
    if (nmax.gt.len(temp_string)) then
      status=1
      goto 10
    endif
    !
    temp_string=trim(adjustl(instring))
    outstring=' '
    !
    n=0
    do i=1,nmax
      if (temp_string(i:i)==delim) then
        n=n+1
      endif
    enddo
    !
    if ((iwant.lt.1).or.(iwant.gt.n+1)) then
      status=1
    else
      if (n==0) then
        outstring=trim(adjustl(temp_string))
      else
        allocate(delim_pos(n))
        n=0
        do i=1,nmax
          if (temp_string(i:i)==delim) then
            n=n+1
            delim_pos(n)=i
            temp_string(i:i)=' '
          endif
        enddo
        !
        if (iwant==1) then
          outstring=temp_string(1:delim_pos(1))
        else if (iwant==n+1) then
          outstring=temp_string(delim_pos(n)+1:nmax) !outstring=temp_string(delim_pos(n):nmax)
        else
          outstring=temp_string(delim_pos(iwant-1)+1:delim_pos(iwant)) !outstring=temp_string(delim_pos(iwant-1):delim_pos(iwant))
        endif
        deallocate(delim_pos)
      endif
      !
    end if
    !
  10  return
  end subroutine get_arg_sim_input
  
  subroutine get_nargs_sim_input(instring,delim,n)
  !
  ! Scans through instring and returns the number 
  ! of sub-strings seperated by `delim'
  !
  character(len=*),intent(in)    :: instring
  integer,intent(out)            :: n
  character(len=1)               :: delim
  !
  integer                        :: nmax,i
    !
    nmax=len(trim(adjustl(instring)))
    !
    n=1
    do i=1,nmax
      if (instring(i:i)==delim) then
        n=n+1
      endif
    enddo
    !
  return
  end subroutine get_nargs_sim_input
  !    
  !======================================================!
  !=== Subroutines for Doppler spectra simulations ======!
  !======================================================!
  !======================================================!
  subroutine processing_ds(isc,conf,spectra_VNyquist,spectra_NOISE_1km,spectra_NFFT,spectra_Nave,&
                          &range_m,elev,ww,sw_dyn,nd,N_d,diam,ddiam_in,fvel_d,zhh_d,zvh_d,zvv_d,&
                          &spectra_bins,zhh_spectra,zvh_spectra,zvv_spectra)
  Use crsim_mod
  Implicit None
  !
  Integer,Intent(In)              :: isc
  Type(conf_var),Intent(in)       :: conf
  Real*8, Intent(In)              :: spectra_VNyquist,spectra_NOISE_1km
  Integer, Intent(In)             :: spectra_NFFT,spectra_Nave
  Real*8, Intent(In)              :: range_m ![m] distance from radar
  Real*8, Intent(In)              :: elev ![deg] elevation angle
  Real*8, Intent(In)              :: ww ![m/s] vertical air velocity
  Real*8, Intent(In)              :: sw_dyn ! [m/s] spectrum broadening due to dynamics 
  Integer, Intent(In)             :: nd   ! number of diameter bins
  Real*8, Intent(In)              :: diam(nd) ! diameters [m]
  Real*8, Intent(In)              :: ddiam_in(nd+1) ! diameter bin boundary from data [m]
  Real*8, Intent(In)              :: N_d(nd) ! number density [m^-3]
  Real*8, Intent(In)              :: fvel_d(nd) ! particle fall velocitat each diameter bin
  Real*8, Intent(In)              :: zhh_d(nd),zvh_d(nd),zvv_d(nd) ! [mm^6] each each diameter
  Real*8, Intent(Out)             :: spectra_bins(spectra_NFFT),zhh_spectra(spectra_NFFT)
  Real*8, Intent(Out)             :: zvh_spectra(spectra_NFFT),zvv_spectra(spectra_NFFT)
  !
  Integer                         :: i
  real*8,dimension(:),Allocatable :: spectra_example,spectra_example_unfold,velo_example,velo_example_unfold
  Real*8,dimension(:),Allocatable:: ddiam ! diameter bin sizes [m]
    !
    ! for each grid point
    !------------------------------------------------------------------------------------
    !
    Allocate(ddiam(nd+1))
    if(MaxVal(ddiam_in) <= 0.d0)then
      ddiam(1)=diam(1) - 0.5d0*(diam(2)-diam(1))
      do i=2,nd
        if(diam(i)<diam(i-1)) then
          if(i>2)  ddiam(i)=diam(i-1) + 0.5d0*(diam(i-1)-diam(i-2))
          if(i<=2) ddiam(i)=diam(1) + 0.5d0*dabs(diam(2)-diam(1))
        else
          ddiam(i)=0.5d0*(diam(i) + diam(i-1)) 
        endif
      enddo
      ddiam(nd+1)=diam(nd) + 0.5d0*dabs(diam(nd)-diam(nd-1))
    else
      do i=1,nd+1
        ddiam(i)=ddiam_in(i)
      enddo
    end if
    !
    Allocate(spectra_example(spectra_NFFT),spectra_example_unfold(spectra_NFFT))
    Allocate(velo_example(spectra_NFFT),velo_example_unfold(spectra_NFFT))
    ! 
    spectra_example=0.d0
    spectra_example_unfold=0.d0
    velo_example=0.d0
    velo_example_unfold=0.d0
    !
    call compute_doppler_spectrum(conf%freq,spectra_VNyquist,spectra_NOISE_1km,spectra_NFFT,spectra_Nave,&
                                  range_m,elev,ww,sw_dyn,&
                                  nd,diam,ddiam,fvel_d,N_d,zhh_d,&
                                  spectra_example,spectra_example_unfold,velo_example,velo_example_unfold)
    spectra_bins = velo_example
    zhh_spectra = spectra_example
    ! 
    if(conf%radID/=1)then !polarimetry
      call compute_doppler_spectrum(conf%freq,spectra_VNyquist,spectra_NOISE_1km,spectra_NFFT,spectra_Nave,&
                                    range_m,elev,ww,sw_dyn,&
                                    nd,diam,ddiam,fvel_d,N_d,zvh_d,&
                                    spectra_example,spectra_example_unfold,velo_example,velo_example_unfold)
      !spectra_bins = velo_example
      zvh_spectra = spectra_example
      call compute_doppler_spectrum(conf%freq,spectra_VNyquist,spectra_NOISE_1km,spectra_NFFT,spectra_Nave,&
                                    range_m,elev,ww,sw_dyn,&
                                    nd,diam,ddiam,fvel_d,N_d,zvv_d,&
                                    spectra_example,spectra_example_unfold,velo_example,velo_example_unfold)
      !spectra_bins = velo_example_unfold
      !zvv_spectra = spectra_example_unfold
      !spectra_bins = velo_example
      zvv_spectra = spectra_example
    else
      zvh_spectra = 0.d0
      zvv_spectra = 0.d0
    end if
    !
    Deallocate(spectra_example,spectra_example_unfold,velo_example,velo_example_unfold)
    Deallocate(ddiam)
    !
  return
  end subroutine processing_ds
  !
  subroutine compute_doppler_spectrum(Frequency,VNyquist,NOISE_1km,NFFT,Nave,range_m,elev,ww,sw_dyn,&
             nd,diam,ddiam,fvel,N_d,zhh_d,&
             spectra_example,spectra_example_unfold,velo_example,velo_example_unfold)
  ! DESCRIPTION
  ! LES Radar Doppler spectra simulator 
  !   The McGill Radar Doppler Spectra Simulator (MRDSS) uses the  explicit (bin) microphysics 
  ! schemes input from high resolution Large Eddies Simulations (LES) and Cloud Resolving
  ! Models (CRMs) and outputs forward modeled radar Doppler spectra moments from profiling radars.
  !   The main features of the MRDSS are the following:
  !       - Accounts for gaseous attenuation estimated using Rosenkranz (1998)
  !       - Account for liquid scattering and attenuation using Mie (1908)
  !       - Provides output at two vertical resolutions: model and radar (specified)
  !       - The model resolution output assumes uniform beam filling conditions
  !       - The radar resolution output accounts for non-uniform beam filling conditions 
  !       - Range of input radar frequencies from 0.1-300 GHz
  !       - Offers a complete radar instrument model (Kollias et al., 2014)
  !       - Estimates the kinematic broadening of the Doppler spectrum using model output
  !       - Includes range-dependent radar receiver noise and signal integration
  !       - Post-processing of the radar Doppler spectrum and moment estimation
  !    Copyright (C) 2014 Pavlos Kollias
  !    Contact email address: pavlos.kollias@stonybrook.edu!
  use crsim_mod
  use phys_param_mod, ONLY: pi
  implicit none
  
  interface
    subroutine Turbulence_Convolution(s, xi, yi, nx, xo, yo)
      integer, intent(in) :: nx
      real*8, intent(in) :: s
      real*8, dimension(nx), intent(in) :: xi, yi
      real*8, dimension(:), allocatable, intent(out) :: xo, yo
    end subroutine Turbulence_Convolution
  end interface
  
  
  Real*8,Intent(in) :: Frequency !!! Radar Frequency GHz
  Real*8,Intent(in) :: VNyquist,NOISE_1km !!!parameters set in crsim_mod.f90 
  Integer,Intent(in):: NFFT,Nave !!!parameters set in crsim_mod.f90
  Real*8,Intent(in) :: range_m !!! Distance from radar in meters
  Real*8,Intent(in) :: elev !!! elevation angle degrees
  Real*8,Intent(in) :: ww !!! vertical velocity 
  Real*8,Intent(in) :: sw_dyn !!! spectra bloadning due to dynamics 
  Integer,Intent(in):: nd       !!! Number of diameter bins
  Real*8,Intent(in) :: diam(nd) !!! diameter [m]
  Real*8,Intent(in) :: ddiam(nd+1) !!! diameter bin size [m]
  Real*8,Intent(in) :: fvel(nd)  !!! particle fall velocity m/s 
  Real*8,Intent(in) :: N_d(nd) !!! number concentration /m^3 for each diameter
  Real*8,Intent(in) :: zhh_d(nd) !, zvv_d(nd)  !!! backscatter for each diameter mm^6
  Real*8,Intent(Out) :: spectra_example(NFFT)   !!!! spectrogram at LES resolution
  Real*8,Intent(Out) :: spectra_example_unfold(NFFT)
  Real*8,Intent(Out) :: velo_example          (NFFT)
  Real*8,Intent(Out) :: velo_example_unfold   (NFFT)
  integer                              :: ibin,iibin
  real*8                               :: st!, SH_magnitude_UV, SH_total, Std_turb, W_air
  real*8                               :: Noise, Npts, Nmean, Nmax
  real*8, dimension(:),    allocatable :: Vb_l, diff_Vb_l
  real*8, dimension(:),    allocatable :: Vd, Sv
  real*8, dimension(:),    allocatable :: velo_st, turb_spectra
  real*8, dimension(:),    allocatable :: spectra_velo, noise_turb_spectra, snr_turb_spectra
  real*8, dimension(:),    allocatable :: spectra_velo_unfold, noise_turb_spectra_unfold
    !______________________________________________________________________
    !nbin = nd
    ! Define the velocity of the hydrometeors at the boundaries of the
    ! droplet bin sizes - we need this for the estimation of the radar
    ! Doppler spectrum
    allocate(Vb_l     (nd+1))     !!! Particle fall velocity at bin boundaries
    allocate(diff_Vb_l(nd  ))     !!! dU for each bin size
     
    if(MaxVal(fvel)==MinVal(fvel) ) then
      diff_Vb_l(1:nd) = 1.d0 
    else
      call interp1(diam, fvel(1:nd), size(diam), ddiam, Vb_l(1:nd+1), size(ddiam), .true.)
      diff_Vb_l(1:nd) = Vb_l(2:nd+1) - Vb_l(1:nd)   !!! this is needed for S(v) = N(D)s(D)dD/dU
      ! To deal with constant fall velocity with diameter (avoid diff_Vb_l = 0)
      if (MinVal(diff_Vb_l)<=1.d-10)then
        do ibin=1,nd
          if(diff_Vb_l(ibin)<=1.d-10)then
            if(ibin>1) diff_Vb_l(ibin) = diff_Vb_l(ibin-1) 
          endif
        enddo
      endif
    endif
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! In CR-SIM, skip the estimatting attenuation.!
    ! To compute Attenuation at each diameter bin, 
    ! use Ah,Av from compute_polalim_var or compute_rad_forward_var
    !  E.g. Sv_att(k,:) = Sv/(10**(Ah/10.0))
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!! STEP 5: Doppler spectra simulator      !!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     
    !!!! STEP 5.1: Memory allocation !!!!!!!!!!!!!!!!
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!! Turbulence variable !!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    allocate(Vd             (nd)  )     !!! Vertical drop velocity incl. vertical motion
    allocate(Sv             (nd)  )     !!! Doppler spectrum (Z/(ms^-1)) 
    allocate(spectra_velo(NFFT), noise_turb_spectra(NFFT),snr_turb_spectra(NFFT))
    allocate(spectra_velo_unfold(NFFT),noise_turb_spectra_unfold(NFFT))
    
    !!!!! END OF STEP 5.1 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!! START DOPPLER SPECTRA SIMULATION LOOP !!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!! we estimate the radar Doppler spectrum at the model resolution
      
    !!!!! STEP 5.2: Define radar volume spatial filters and noise !!!!!!!!!!!!!!
    Noise = NOISE_1km*(range_m/1000.d0)**2.d0

    !!!!! STEP 5.3: Estimate the Quiet-Air (no turbulence/shear) Doppler spectrum !!!!!!
    Sv          = (N_d*zhh_d/diff_Vb_l) !*1.0e6       !!! N(D)*sb(D)*dD/dV! 
    !! replace Sv with 0 when diff_Vb_l = 0
    if(MinVal(diff_Vb_l)<=1.0d-10)then
      do ibin=1,nd
        if(diff_Vb_l(ibin)<=1.0d-10) Sv(ibin)=0.d0
      enddo
    endif
    !!!! Here we add the gaseous and liquid attenuation !!!!!!!!!!
    !!Sv_att(k,:) = Sv/(10**(A_liquid_accum/10.0))/(10**(A_gases/10.0))
    
    !!!! In the particular model used, the vertical air motion is positive up and the particle fall velocity
    !!!! is positive down. In the simulated spectra we use ARM convection (toward the radar/ground is negative)
    Vd           = ww - (fvel(:)*dsin(elev*pi/180.d0)) !!!! the sign used here maybe model-dependent
    if((elev>=-90.d0) .and. (elev<0.d0)) Vd = Vd * (-1.d0) ! added for negative elevation by oue 20190911
    
    !!!! if fall velocity is constant
    if(dabs(MaxVal(Vd)-MinVal(Vd))<0.01d0) then !if all Vd are constant
      Sv(2)=Sum(Sv(:))/0.01d0; Sv(3:nd)=0.d0; Sv(1)=0.d0
      Vd(1)=MaxVal(Vd)+0.01d0; Vd(3)=MaxVal(Vd)-0.01d0
    else ! if Vd are constant for larger diameters
      ibin=1
      do while (ibin<nd)
        if(Vd(ibin)==Vd(ibin+1))then
          do iibin=ibin+1,nd
            if(Vd(ibin)/=Vd(iibin)) exit
            Sv(ibin)= Sv(ibin)+ Sv(iibin)
            Sv(iibin)=0.d0
          enddo
          ibin = iibin
        else
          ibin=ibin+1
        endif
      enddo
    end if

    !if(MaxVal(N_d)>0)then
      !do ibin=1,nd
        !write(*,*) ibin,Sv(ibin),Vd(ibin),fvel(ibin),diff_Vb_l(ibin),zhh_d(ibin),N_d(ibin)
      !enddo
    !endif
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!! STEP 5.4: Estimate Shear and subgrid turbulence spectral broadening !!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! In CR-SIM, spectral broadening has been calculated. This is an input value.    !!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
    !!!! Total spectral broadening estimation    
    st = sw_dyn
    if(st .le. 0.0d0) st = 0.0000001d0

    !!!!!! STEP 5.5 Perform turbulence convolution  !!!!!!!!!!!!!!!
    call Turbulence_Convolution(st, Vd, Sv(:), size(Vd), velo_st, turb_spectra)

    !!!!!! STEP 5.6: Estimate the radar Doppler spectra with Noise + broadening !!!
    call spectra_generator(NFFT, Noise,Nave, VNyquist,                   &
                           velo_st, turb_spectra, size(velo_st),             &
                           spectra_velo, noise_turb_spectra, snr_turb_spectra)

!!!! added for negative elevation by oue 20190911
    if((elev>=-90.d0) .and. (elev<0.d0)) then
        spectra_velo       = spectra_velo(NFFT:1:-1)
        noise_turb_spectra = noise_turb_spectra(NFFT:1:-1)
        snr_turb_spectra   = snr_turb_spectra(NFFT:1:-1)
    endif
!!!!
    
    spectra_example       (1:NFFT) = noise_turb_spectra
    velo_example          (1:NFFT) = spectra_velo

    !write(*,*) 'spectra_generator',10*dlog10(Sum(Sv*diff_Vb_l)),10*dlog10(Sum(turb_spectra*0.01))
     
    !!!!!! STEP 5.7: Estimate the noise floor of the radar Doppler spectrum !!! 
    call hsmethod(noise_turb_spectra, Nave ,NFFT, Npts, Nmean, Nmax)

    !!!!!! STEP 5.8: Unfold the radar Doppler spectra before moment estimation !!!
    call spectra_unfolding(spectra_velo, noise_turb_spectra, VNyquist, NFFT, &
                           spectra_velo_unfold, noise_turb_spectra_unfold    )

    spectra_example_unfold(1:NFFT) = noise_turb_spectra_unfold
    velo_example_unfold   (1:NFFT) = spectra_velo_unfold

    do ibin=1,NFFT
      if(spectra_example(ibin)<=0.d0) then
        spectra_example(ibin)=0.d0
      endif   
      if(spectra_example_unfold(ibin)<=0.d0) then
        spectra_example_unfold(ibin)=0.d0
      endif
      !if(velo_example(ibin)<-999.0) spectra_example(ibin)=VNyquist*(-1)
      !if(velo_example_unfold(ibin)<-999.0) spectra_example(ibin)=VNyquist*(-1)
    end do
    deallocate(velo_st)
    deallocate(turb_spectra)
      
    deallocate(Vb_l     )
    deallocate(diff_Vb_l)
     
    deallocate(spectra_velo, noise_turb_spectra,snr_turb_spectra)
    deallocate(spectra_velo_unfold,noise_turb_spectra_unfold)
   
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! END OF STEP 5 !!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     
  end subroutine compute_doppler_spectrum
  !   
  subroutine Turbulence_Convolution(st, VL, L, nsizes, vel, spectra)
  !
  ! This funcion adds turbulence to the spectra.
  ! The non-turbulent radar Doppler spectrum  is convoluted with a Gaussian function
  ! with standard deviation equal to the total kinematical broadening of the radar Doppler
  ! spectrum according to Gossard, 1994.
  ! 
  Use crsim_mod
  Use phys_param_mod, ONLY : pi
  Implicit None
  ! 
  !------------------------------------------------------
  interface
    subroutine conv2(a, n1, b, n2, c)
      integer, intent(in) :: n1, n2
      real*8, dimension(n1), intent(in) :: a
      real*8, dimension(n2), intent(in) :: b
      real*8, dimension(:), allocatable, intent(out) :: c
    end subroutine conv2
  end interface
  !------------------------------------------------------ 
  ! 
  integer, intent(in) :: nsizes
  real*8, intent(in) :: st
  real*8, dimension(nsizes), intent(in) :: VL, L
  real*8, dimension(:), allocatable, intent(out) :: vel, spectra
  real*8, parameter :: dV = 0.01d0              !!! resolution of 0.1 cm/sec (1 cm/s?)
  integer :: L_v, L_t, idx, idx_off
  real*8 :: vmin, vmax, rr
  real*8, dimension(:), allocatable :: Sv_hr, t, turb, C
    ! 
    ! turbulence broadening standard deviation st, typical range [0.1 - 0.4]
    ! the fall velocity of the liquid (m/sec) 
    ! the reflectivity of liquid (linear units)
    ! st is the kinematic spectral broadening m s-
    ! 
    !!!!! Interpolate S(v)dV to a high resolution, regular V-space
    vmin = minval(VL)-5.d0*st
    vmax = maxval(VL)+5.d0*st
    L_v  = int((vmax-vmin)/dV)+1
     
    allocate(Sv_hr  (L_v))
    allocate(vel    (L_v))
    allocate(spectra(L_v))
     
    vel = vmin + dV*(/(idx, idx=0, L_v-1, 1)/)
    call interp1(VL(nsizes:1:-1), L(nsizes:1:-1), nsizes, vel, Sv_hr, L_v, .false.)
    !write(*,*) 'Turb',10*dlog10(Sum(Sv_hr*0.01))
    !!! interpolate N(D)s(D)dD/dV to high resolution, regular bins 
    !!!! add turbulence to the spectra !!!!
    L_t = int(4.d0/dV)+1
    rr  = st/dV    !!! width in high resolution bins  
    
    allocate(t   (L_t))
    allocate(turb(L_t))
    
    t     = -2.d0/dV + (/(idx, idx=0, L_t-1, 1)/)
    turb  = 1.d0/(dsqrt(2.d0*pi)*rr)*dexp(-(t**2)/(rr**2))
    
    call conv2(Sv_hr, L_v, turb, L_t, C)
   
    idx_off = int(2.d0/dV)+1
    spectra = C(idx_off:idx_off+L_v-1)
    
  end subroutine Turbulence_Convolution
  !===============================================================================
  subroutine spectra_generator(NFFT, Noise, Nave, VNyquist, VL, L, nsize,   &
                               sp_velo, noise_turb_spectra, snr_turb_spectra)
  ! DESCRIPTION
  ! This function calculates the final form of the Doppler spectrum. 
  ! 
  implicit none
  
  integer, intent(in) :: NFFT, Nave, nsize
  real*8, intent(in) :: Noise, VNyquist
  real*8, dimension(nsize), intent(in) :: VL, L
  real*8, dimension(NFFT), intent(out) :: sp_velo, noise_turb_spectra, snr_turb_spectra
  
  real*8, parameter :: Pn = 1.d0
  real*8, dimension(3*NFFT)    :: spectra_velo, spectra_power
  real*8, dimension(NFFT)      :: sp_power, x_n
  real*8, dimension(Nave,NFFT) :: noise_n_turb_spectra
  real*8 :: dV, Signal, SNR, K
  integer :: i

    !! number of FFT points in the Doppler spectrum 
    !! Noise power (in mm^6/m^3 linear units)
    !! number of average spectra for noise variance reduction, typical range [1 40]
    !! NyquistVelocity in m/sec --
    !! the Doppler velocity of the distributed targets (m/sec)
    !! the reflectivity (mm^6/m^3/(ms^-1)) of the targets (linear units)
    
    !!!!! Create three times the velocity domain to treat aliasing... 
    
    dV = 2.d0*VNyquist/NFFT                         !!! Velocity resolution 
    spectra_velo  = -3.d0*VNyquist + dV*(/(i, i=0, 3*NFFT-1, 1)/)  !!! three times wide velocity domain !!!           
    call interp1(VL, L, nsize, spectra_velo, spectra_power, 3*NFFT, .false.) !!!interpolate N(D)s(D)dD/dV to spectral bins !!!
    
    !!!! Now fold the spectrum back to -VNyquist to +VNyquist-dv domain !!!!! 
    sp_power = spectra_power(       1:  NFFT) &
             + spectra_power(  NFFT+1:2*NFFT) &
             + spectra_power(2*NFFT+1:3*NFFT)
    sp_velo  = -VNyquist + dV*(/(i, i=0, NFFT-1, 1)/)
    
    !!! Add noise to the spectrum !!!!!!!!!! 
    Signal = sum(sp_power)*dV
    ! Noise = 10^(Noise/10)!  
    
    !!!! Add SNR and noise to the spectra  
    if (Signal <= 0.d0) then
      K = 0.d0
    else
      SNR = 10.d0*dlog10(Signal/Noise)
      if(SNR.lt.-20.d0) then
        K = 0.d0
      else
        K = 10.d0**(SNR/10.d0)/Signal
      end if
    end if
    
    !===============================================================!          
    !== temporary added by oue to not use K     
    !snr_turb_spectra = K*sp_power + 2.0*Pn/NFFT 
    snr_turb_spectra = sp_power + 2.d0*Noise/NFFT
    !===============================================================!
    
    !!!!!! if Nave >0  stack spectra together to reduce noise !!!!! 
    !!! Memory allocation !!!!
    do i=1, Nave
      call random_number(x_n)
      x_n = x_n + 1.0d-6
      noise_n_turb_spectra(i,:) = -dlog(x_n)*snr_turb_spectra
    end do
    noise_turb_spectra = sum(noise_n_turb_spectra,1)/dble(Nave)
    
  end subroutine spectra_generator
  !===============================================================================
  !!!!! Doppler spectra unfolding algorithm !!!!
  subroutine spectra_unfolding(spectra_velo, noise_turb_spectra, VNyquist, NFFT, &
                               spectra_velo_unfold, noise_turb_spectra_unfold    )
  ! DESCRIPTION
  ! This function  unfolds the radar Doppler spectra only if the fall velocity 
  ! of the hydrometeors is large enough in order to produce partial or complete
  ! aliasing of the all the velocities (frequencies) of the hydrometeors.
  ! The main output of the spectral_unfolding function is the spectral velocity array 
  ! (spectra_velo_unfold), and the unfolded radar Doppler spectrum (noise_turb_spectra_unfold)
  implicit none
  !!! Estimate mean Doppler velocity:
  integer, intent(in) :: NFFT
  real*8,    intent(in) :: VNyquist
  real*8, dimension(NFFT), intent(in ) :: spectra_velo, noise_turb_spectra
  real*8, dimension(NFFT), intent(out) :: spectra_velo_unfold, noise_turb_spectra_unfold
  real*8 :: V_m
  integer :: LU, LA
  
    V_m = sum(spectra_velo*noise_turb_spectra)/sum(noise_turb_spectra)
    
    LU = NFFT/4
    LA = NFFT/4*3
    
    if (V_m > 0.4d0*VNyquist) then
      !!!!! Large downward Doppler velocities
      spectra_velo_unfold(   1:LA  ) = spectra_velo(LU+1:NFFT)
      spectra_velo_unfold(LA+1:NFFT) = spectra_velo(   1:LU  ) + 2*VNyquist
      noise_turb_spectra_unfold(   1:LA  ) = noise_turb_spectra(LU+1:NFFT)
      noise_turb_spectra_unfold(LA+1:NFFT) = noise_turb_spectra(   1:LU  )
    else
      if (V_m < -0.4d0*VNyquist) then
        !!!!! Updraft ......
        spectra_velo_unfold(   1:LU  ) = spectra_velo(LA+1:NFFT) - 2*VNyquist
        spectra_velo_unfold(LU+1:NFFT) = spectra_velo(   1:LA  )
        !spectra_velo_unfold(LU:NFFT-1) = spectra(LU:NFFT-1) 
        noise_turb_spectra_unfold(   1:LU  ) = noise_turb_spectra(LA+1:NFFT)
        noise_turb_spectra_unfold(LU+1:NFFT) = noise_turb_spectra(   1:LA  )
        !noise_turb_spectra_unfold(LU+1:NFFT)=noise_turb_spectra(LU+1:NFFT)
      else
        spectra_velo_unfold       = spectra_velo
        noise_turb_spectra_unfold = noise_turb_spectra
      end if
    end if
    ! 
  end subroutine spectra_unfolding
  !   
  subroutine interp1(xi, yi, ni, xo, yo, no, opt_extra)
  implicit none
  !   
  integer, intent(in) :: ni, no
  real*8, dimension(ni), intent(in)  :: xi, yi
  real*8, dimension(no), intent(in)  :: xo
  real*8, dimension(no), intent(out) :: yo
  logical, intent(in) :: opt_extra
  integer :: i, j
   
    do j=1, no
      ! extrapolation 
      !  if(xo(j).lt.xi(1)) then
      if(xo(j).lt.MinVal(xi)) then
        if(opt_extra.eqv..true.) then
          yo(j) = (yi(2)-yi(1))/(xi(2)-xi(1))*(xo(j)-xi(1))+yi(1)
        else
          yo(j) = 0.d0
        end if
      else if(xo(j).ge.xi(ni)) then
        if(opt_extra.eqv..true.) then
          yo(j) = (yi(ni)-yi(ni-1))/(xi(ni)-xi(ni-1))*(xo(j)-xi(ni))+yi(ni)
        else
          yo(j) = 0.d0
        end if
      else
        ! interpolation
        do i=1, ni-1
          if(xo(j).ge.xi(i).and.xo(j).lt.xi(i+1)) then
            yo(j) = (yi(i+1)-yi(i))/(xi(i+1)-xi(i))*(xo(j)-xi(i))+yi(i)
            exit
          end if
        end do
      end if
    end do
     
  end subroutine interp1
     
     
  subroutine hsmethod(S_n, nave, nffts, N_points, N_mean, N_max)
  ! After sorting the spectral points in order of ascending 
  ! power, the algorithm goes through all values of n from 0 to numPts-1 
  ! and for each n computes the values of both sides of the above equation.
  ! For early values of n, the left side is generally less than the 
  ! right and then becomes greater than the right side as higher-valued 
  ! points that do not belong to the original distribution are included
  ! in the n points. The last crossing from less than to greater than
  ! is considered the division between the desired set of points and
  ! "outliers" (e.g. noise vs. signal! noise or signal vs. interference,
  ! etc.).
  !
  ! Hildebrand, P. H., and R. S. Sekhon, Objective determination of
  ! the noise level in Doppler spectra, J. Appl. Meteorol., 13, 808, 1974.
  ! 
  ! sum(S^2)/N - (sum(S)/N)^2 = (mean(S))^2 / navg
  ! navg*[sum(S^2)/N - (sum(S))^2 / N^2] = (sum(S))^2 / N^2
  ! navg*[N*sum(S^2) - sum(S) * sum(S)] = sum(S) * sum(S)
  implicit none
  integer, intent(in) :: nave, nffts
  real*8, dimension(nffts), intent(in) :: S_n
  real*8, intent(out) :: N_points, N_mean, N_max
  real*8, dimension(nffts) :: P, a1, a3
  real*8 :: Nthld, sumLi, sumSq, sumNs, maxNs
  integer :: i, n, numNs
     
    P = S_n
    call quicksort(P,1,nffts)
     
    Nthld = 1.0d-10  !is it really needed? (had to decrease it from 10**-6)
   
    n     = 0
    sumLi = 0.d0
    sumSq = 0.d0
    
    numNs = 0
    sumNs = 0.d0
    maxNs = 0.d0
   
    a1 = 0.d0
    a3 = 0.d0
     
    do i=1, nffts
      if (P(i).gt.Nthld) then
        sumLi = sumLi + P(i)
        sumSq = sumSq + P(i)*P(i)
        n = n + 1
        a3(i) = sumLi*sumLi
        a1(i) = dsqrt(dble(nave))*(n*sumSq - a3(i))
        !a1(i) = nave*(n*sumSq-a3(i))!
        if (n.gt.nffts/4) then
          if (a1(i).le.a3(i)) then
            sumNs = sumLi
            numNs = n
            maxNs = P(i)
          end if
        else
          sumNs = sumLi
          numNs = n
          maxNs = P(i)
        end if
      end if
    end do
     
    if(numNs.gt.0) then
      N_mean   = sumNs/numNs
      N_max    = maxNs
      N_points = dble(numNs)
    else
      N_mean   = 0.d0
      N_max    = 0.d0
      N_points = 0.d0
    end if
     
  end subroutine hsmethod
     
  subroutine conv2(a, n1, b, n2, c)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! conv : return the convoltion of vectos a and b
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  implicit none
  
  integer, intent(in) :: n1, n2
  real*8, dimension(n1), intent(in) :: a
  real*8, dimension(n2), intent(in) :: b
  real*8, dimension(:), allocatable, intent(out) :: c
  integer :: n, i, i1, i2
  
    n = n1+n2-1
    allocate(c(n))
   
    do i=1, n
      i1 = max(i-n2+1,1)
      i2 = min(i,n1)
     
      !-- M. Oue modified below so that the calculated Zhh is equal to Sum(spectra)
      !-- This could be like weighting average filtering.
      !  c(i) = sum(a(i1:i2)*b(i+1-i1:i+1-i2:-1))
      c(i) = sum(a(i1:i2)*b(i+1-i1:i+1-i2:-1)) / sum(b(i+1-i1:i+1-i2:-1))
      ! 
    end do
    ! 
  end subroutine conv2
  !   
  recursive subroutine quicksort(a, first, last)
  ! quicksort.f -*-f90-*-
  ! Author: t-nissie
  ! License: GPLv3
  ! Gist: https://gist.github.com/t-nissie/479f0f16966925fa29ea
  !
  implicit none
     
  real*8 ::  a(*), x, t
  integer :: first, last
  integer :: i, j
     
    x = a( (first+last) / 2 )
    i = first
    j = last
     
    do
      do while (a(i) < x)
        i=i+1
      end do
    do while (x < a(j))
      j=j-1
    end do
    if (i >= j) exit
      t = a(i);  a(i) = a(j);  a(j) = t
      i = i+1
      j = j-1
    end do
     
    if (first < i-1) call quicksort(a, first, i-1)
    if (j+1 < last)  call quicksort(a, j+1, last)
    
  end subroutine quicksort
  !
  subroutine initialize_spectra_var(str,Frequency,ZMIN)
  Use crsim_mod
  Use phys_param_mod, ONLY: m999, c, &
                          & TimeSampling,NOISE_1km,NFFT, &
                          & PRF_S,PRF_C,PRF_X,PRF_Ka,PRF_Ku,PRF_W
  Implicit None
  Type(spectra_var),Intent(InOut)    :: str
  Real*8,Intent(In)                  :: Frequency      ! radar frequency in GHz
  Real*8,Intent(In)                  :: ZMIN      ! minimum ditectable Z at 1 km
                                                  ! [dBZ] specified in PARAMETERS
    !
    str%TimeSampling = TimeSampling
     
    str%PRF = m999
    if(Frequency == 3.0d0) str%PRF = PRF_S
    if(Frequency == 5.5d0) str%PRF = PRF_C
    if(Frequency == 9.5d0) str%PRF = PRF_X
    if(Frequency == 35.0d0)str%PRF = PRF_Ka
    if(Frequency == 13.6d0)str%PRF = PRF_Ku
    if(Frequency == 94.0d0)str%PRF = PRF_W
    
    if (str%PRF < 0.d0) then
      write(*,*) 'Problem with PRF in initialize_spectra_var'
      stop
    endif
     
    str%NOISE_1km    = NOISE_1km
    if ( (ZMIN > m999) .and. (ZMIN < -m999) ) &
      str%NOISE_1km = 10.d0 ** (0.1d0 * ZMIN)
    str%NFFT         = NFFT
    str%Lambda       = c/(Frequency*1.d+9)
    str%VNyquist     = 0.25d0* str%Lambda*str%PRF
    str%Nave = nint(TimeSampling*str%PRF/NFFT)
    !
    if (str%NFFT <= 0) then
      write(*,*) 'Problem with NFFT in initialize_spectra_var'
      stop
    endif
    ! 
  return
  end subroutine initialize_spectra_var
  !
  !======================================================!
  !======================================================!
  !=== END Subroutines for spectra simulations ==========!
  !======================================================!
  !======================================================!
     
