!! ----------------------------------------------------------------------------
!! ----------------------------------------------------------------------------
!!  *PROGRAM* WriteOutNetcdf
!!  @version crsim 4.x
!!
!!  *SRC_FILE*
!!  crsim/src/WriteOutNetcdf.f90
!!
!!
!!
!!  *LAST CHANGES*
!!
!!  Sep 17 2015  -  A.T.   a new variable, azimuth, is added in the output
!!                         variables list 
!!  Sep 23 2015  -  A.T.   changed the order of the output vars, existing doppler variables 
!!                         are no more vertical but for given elevation and azimuth;
!!                         introduced new variables DV90 and SPh90 for vertical (elevation=90) pointing,
!!                         added u,v wind components 
!!  Oct 15 2015  -  A.T.   added dry air density to the output list
!!
!!  Jan 04 2016  -  A.T    added missing lines for writting the u and v variables into the netcdf files for different species
!!
!!  Jan 08 2016 -  A. T    the description attributes for velocity variables are slightly modfied to account for change in notation:
!!                         all velocity fields are defined now as positive upward (away from the radar)
!!  Jan 12 2016 -  A. T    added the "missing_value" attribute (mainly for the proper plotting when using ncview)
!!  Jan 13 2016 -  A. T    xlat,xlong added in the output netcdf file
!!  Jan 14 2016 -  A. T    changed name of the spect. width variable SPh to SWh (and SPh90 to SWh90) and this is now sp. width due to
!!                         hydrometeor contribution. Added to the output are:SWt (spectrum width due to turbulence) and SWtot (total spectrum width)  
!!  Jan 15 2016  -  A. T   Added to the output SWs (spectrum width due to wind chear in radar volume) and SWv (spectrum width due to cross-wind)
!! 
!!  Feb 15 2016  -  A.T.   a new variable, range, is added in the output variables list  
!!  Feb 15 2016  -  A.T.   the attribute "missing_value" replaced by "_FillValue"
!!  Feb 17 2016  -  A.T.   added new optional output variables : lidar observed and true backscatter, extinction coeff and lidar ratio.
!!  Feb 17 2016  -  A.T.   added new output variable diff_back_phase in the polarimatric netcdf file
!!  Apr 18 2016  -  A.T.   added new variables in the output : rad_freq,rad_beamwidth,rad_range_resolution
!!  Apr 19 2016  -  A.T.   added new variables in the output: rad_ixc,rad_iyc,rad_zc 
!!
!!  May  2 2016  -  A.T.   bug corrected when writing the variables with dimension 1
!!  May  2 2016  -  A.T.   the name of ceilo variable "lidar_ratio" changed to "ceilo_lidar_ratio"  
!!  May  2 2016  -  A.T.   added new optional output variables : mpl observed and true backscatter, extinction coeff and lidar ratio for cloud and ice 
!!  May  9 2016  -  A.T.   added new optional output variables : aerosol observed and true backscatter, extinction coeff and aerosol lidar ratio
!!  May  9 2016  -  A.T.   corrected units for ceilo,mpl and aero extinction (m-1 sr-1 to m-1) and lidar ratio (unitless to sr)
!!  Mar 23 2017  -  M.O.   added arscl products and MWR LWP
!!  JUL 21 2017    -M.O     Incorporated  RAMS 2-moment microphysics.(8 categories) 
!!  May 27 2018    -M.O    added new output for Doppler spectra simulation
!!  Dec    2019    -A.T    added new polarimetric variable in the output: cross-correlation coefficient RHOhv 
!!  Oct    2023    -M.O.   added new output for airborne radar simulation (Zsfc, Dopp_airborne)
!!  Jul    2024    -M.O.   added new output Avap: water vapor attenuation and Atot: 2-way total attenuation
!!
!!  *DESCRIPTION* 
!!
!!  This program  contains the subroutines needed for writing the output
!!  netcdf data files 
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
subroutine WriteOutNetcdf_mp(iht,OutFileName,rmout,mout,qq,env,conf,lout,mpl,aero,arscl,mwr,spectra,airborne,status) !arscl&mwr are added by oue 2017.03.23, added spectra by oue in May 2018
Use netcdf
Use crsim_mod
Use postprocess_mod !added by oue 2017.03.23
Implicit None

integer, intent(in)                    :: iht
character(len=*),Intent(in)            :: OutFileName
type(rmout_var),Intent(In)             :: rmout
type(mout_var),Intent(In)              :: mout
Real*8,Intent(In)                      :: qq(rmout%nx,rmout%ny,rmout%nz) ! water content [kg/m^3]
type(env_var), Intent(in)              :: env
type(conf_var), Intent(in)             :: conf
type(lout_var),Intent(in)              :: lout
type(mpl_var),Intent(in)               :: mpl
type(aero_var),Intent(in)              :: aero
type(arscl_var),Intent(in)             :: arscl !added by oue 2017.03.23
type(mwr_var),Intent(in)               :: mwr !added by oue 2017.03.25
type(spectra_var),Intent(in)           :: spectra !added by oue 2018.05.27
type(airborne_var),Intent(in)          :: airborne !added by oue Oct 2023
Integer,Intent(out)                    :: status
!
character(len=160)                     :: err_msg
!
integer  :: nxID,nyID,nzID,nhtID 
integer  :: ZhhID,ZvvID,ZvhID
integer  :: ZdrID,LDRhID,RHOhvID
integer  :: DVhID
integer  :: SWhID,SWh90ID
integer  :: SWtID,SWtotID
integer  :: SWsID,SWvID
integer  :: DoppID,Dopp90ID
integer  :: KdpID,AdpID,AhID,AvID
integer  :: elevID,azimID,rangeID,ZminID,AvapID,AtotID
integer  :: zID, tempID,wID,qhydroID,uID,vID,rho_dID
integer  :: xID,yID
integer  :: xlatID,xlongID
integer  :: ceilo_back_trueID,ceilo_extID,ceilo_back_obsID,lidar_ratioID
integer  :: ceilo_first_cloud_baseID
integer  :: mpl_wavelID,mpl_back_trueID,mpl_extID,mpl_back_obsID, mpl_lidar_ratioID, mpl_rayleigh_backID
integer  :: aero_back_trueID,aero_extID,aero_back_obsID,aero_lidar_ratioID
integer  :: diff_back_phaseID
integer  :: oneID, rad_freqID,rad_beamwidthID,rad_range_resolutionID
integer  :: ixcID,iycID,zcID
!
!---------------------------
! For ARSCL added by oue 2017.03.23
integer  :: nlayersID
integer  :: arscl_cloud_maskID,arscl_cloud_detection_flagID,arscl_cloud_base_heightID,arscl_cloud_top_heightID
!--------------------------
! For MWR LWP added by oue 2017.03.25
integer  :: true_lwpID, mwr_lwpID, mwr_samplesID
!--------------------------
! For Doppler spectra added by oue in May 2018
Integer  :: nfftID
integer  :: spectra_zhhID, spectra_zvhID, spectra_zvvID, velbinsID
!--------------------------
! For airborne added by oue in Oct 2023
integer  :: Zsfc_oceanID, Zsfc_land_coarseID, Zsfc_land_flatID, DV_airborneID
!--------------------------
!
integer  :: ncid
!
integer                                  :: ihc
Character(len=nf90_max_name)             :: att_string
Character(len=10)                        :: str1,str2
integer                                  :: j1,j2
integer,parameter                        :: one=1
integer,Allocatable,Dimension(:)         :: iwork 
real,Allocatable,Dimension(:)            :: work
!
!
status =  nf90_create(path =OutFileName, cmode = NF90_NETCDF4, ncid = ncid)
If (status /= 0) then ; err_msg='Error in nf90_create' ;  Goto 999 ; endif

! define dimensions
status = nf90_def_dim(ncid, "nx",rmout%nx,nxID)
status = nf90_def_dim(ncid, "ny",rmout%ny,nyID)
status = nf90_def_dim(ncid, "nz",rmout%nz,nzID)
status = nf90_def_dim(ncid, "nht",rmout%nht,nhtID)
status = nf90_def_dim(ncid, "one",one,oneID)
if (conf%spectraID==1) then ! For Doppler spectra by oue in May 2018
status = nf90_def_dim(ncid, "nfft",spectra%NFFT,nfftID)
endif
if ((conf%arsclID==1) .and. (iht==0)) then ! For ARSCL products by oue 2017.03.23
status = nf90_def_dim(ncid, "n_layers",arscl%n_layers,nlayersID)
endif
!
! define variables
!
status = nf90_def_var(ncid, "xlat",  nf90_float, (/nxID,nyID/),xlatID)
status = nf90_def_var(ncid, "xlong", nf90_float, (/nxID,nyID/),xlongID)
!
status = nf90_def_var(ncid, "Zhh", nf90_float, (/nxID,nyID,nzID/), ZhhID)
status = nf90_def_var(ncid, "Zvv", nf90_float, (/nxID,nyID,nzID/), ZvvID)
status = nf90_def_var(ncid, "Zvh", nf90_float, (/nxID,nyID,nzID/), ZvhID)
status = nf90_def_var(ncid, "Zdr", nf90_float, (/nxID,nyID,nzID/), ZdrID)
status = nf90_def_var(ncid, "LDRh", nf90_float, (/nxID,nyID,nzID/),LDRhID)
status = nf90_def_var(ncid, "RHOhv", nf90_float, (/nxID,nyID,nzID/),RHOhvID)
status = nf90_def_var(ncid, "DV", nf90_float, (/nxID,nyID,nzID/), DoppID)
status = nf90_def_var(ncid, "SWh", nf90_float, (/nxID,nyID,nzID/), SWhID)
status = nf90_def_var(ncid, "SWt", nf90_float, (/nxID,nyID,nzID/), SWtID)
status = nf90_def_var(ncid, "SWs", nf90_float, (/nxID,nyID,nzID/), SWsID)
status = nf90_def_var(ncid, "SWv", nf90_float, (/nxID,nyID,nzID/), SWvID)
status = nf90_def_var(ncid, "SWtot", nf90_float, (/nxID,nyID,nzID/), SWtotID)
status = nf90_def_var(ncid, "DV90", nf90_float, (/nxID,nyID,nzID/), Dopp90ID)
status = nf90_def_var(ncid, "SWh90", nf90_float, (/nxID,nyID,nzID/), SWh90ID)
status = nf90_def_var(ncid, "RWV", nf90_float, (/nxID,nyID,nzID/), DVhID)
!
status = nf90_def_var(ncid, "Kdp", nf90_float,(/nxID,nyID,nzID/),KdpID)
status = nf90_def_var(ncid, "Adp", nf90_float,(/nxID,nyID,nzID/),AdpID)
status = nf90_def_var(ncid, "Ah", nf90_float, (/nxID,nyID,nzID/),AhID)
status = nf90_def_var(ncid, "Av", nf90_float, (/nxID,nyID,nzID/),AvID)
status = nf90_def_var(ncid, "diff_back_phase",nf90_float,(/nxID,nyID,nzID/),diff_back_phaseID)
!
!---------------------------
! CEILO modified by oue for RAMS
IF ((conf%ceiloID==1) .and. (iht==0)) THEN
status = nf90_def_var(ncid, "ceilo_back_obs", nf90_float, (/nxID,nyID,nzID/),ceilo_back_obsID)
status = nf90_def_var(ncid, "ceilo_back_true",nf90_float, (/nxID,nyID,nzID/),ceilo_back_trueID)
status = nf90_def_var(ncid, "ceilo_ext",nf90_float,(/nxID,nyID,nzID/),ceilo_extID) 
status = nf90_def_var(ncid, "ceilo_lidar_ratio",nf90_float, (/nxID,nyID,nzID/),lidar_ratioID)
status = nf90_def_var(ncid, "ceilo_first_cloud_base",nf90_float, (/nxID,nyID/),ceilo_first_cloud_baseID)
ENDIF
!---------------------------
!---------------------------
! MPL
IF (conf%mplID>0) THEN 
IF (iht==0) THEN
status = nf90_def_var(ncid, "mpl_rayleigh_back",    nf90_float,(/nxID,nyID,nzID/),mpl_rayleigh_backID)
ENDIF
!
IF ( (iht==0) .or. (iht==1) .or. (iht==3).or.(iht==7) ) THEN
status = nf90_def_var(ncid, "mpl_wavel",  nf90_float,(/oneID/),mpl_wavelID)
status = nf90_def_var(ncid, "mpl_back_obs",     nf90_float, (/nxID,nyID,nzID/),mpl_back_obsID)
status = nf90_def_var(ncid, "mpl_back_true",    nf90_float, (/nxID,nyID,nzID/),mpl_back_trueID)
status = nf90_def_var(ncid, "mpl_ext",          nf90_float, (/nxID,nyID,nzID/),mpl_extID)
status = nf90_def_var(ncid, "mpl_lidar_ratio",  nf90_float, (/nxID,nyID,nzID/),mpl_lidar_ratioID)
ENDIF
! AEROSOL
IF(conf%aeroID==1) THEN
IF (iht==0) THEN
status = nf90_def_var(ncid, "aero_back_obs",     nf90_float,(/nxID,nyID,nzID/),aero_back_obsID)
status = nf90_def_var(ncid, "aero_back_true",    nf90_float,(/nxID,nyID,nzID/),aero_back_trueID)
status = nf90_def_var(ncid, "aero_ext",          nf90_float,(/nxID,nyID,nzID/),aero_extID)
status = nf90_def_var(ncid, "aero_lidar_ratio",  nf90_float,(/nxID,nyID,nzID/),aero_lidar_ratioID)
ENDIF
ENDIF
!---------------------------
ENDIF !if conf%mplID>0
!---------------------------
!---------------------------
! Doppler spectra by oue in May 2018
IF (conf%spectraID>0) THEN 
status = nf90_def_var(ncid, "velocity_bins",  nf90_float,(/nfftID/),velbinsID)
status = nf90_def_var(ncid, "spectra_zhh",    nf90_float,(/nxID,nyID,nzID,nfftID/),spectra_zhhID)
status = nf90_def_var(ncid, "spectra_zvh",    nf90_float,(/nxID,nyID,nzID,nfftID/),spectra_zvhID)
status = nf90_def_var(ncid, "spectra_zvv",    nf90_float,(/nxID,nyID,nzID,nfftID/),spectra_zvvID)
ENDIF
!---------------------------
! Airborne by oue in Oct 2023
IF ((conf%airborne==1) .and. (iht==0)) THEN 
status = nf90_def_var(ncid, "Zsfc_Ocean",        nf90_float,(/nxID,nyID/),Zsfc_oceanID)
status = nf90_def_var(ncid, "Zsfc_Land_Coarse",  nf90_float,(/nxID,nyID/),Zsfc_land_coarseID)
status = nf90_def_var(ncid, "Zsfc_Land_Flat",    nf90_float,(/nxID,nyID/),Zsfc_land_flatID)
status = nf90_def_var(ncid, "DV_airborne",       nf90_float,(/nxID,nyID,nzID/),DV_airborneID)
ENDIF
!---------------------------
! ARSCL added by oue 2017.03.23
if ((conf%arsclID==1) .and. (conf%mplID>0) .and. (conf%ceiloID==1)) then 
if (iht==0) then
status = nf90_def_var(ncid, "arscl_cloud_mask", nf90_short,(/nxID,nyID,nzID/),arscl_cloud_maskID)
status = nf90_def_var(ncid, "arscl_cloud_source_flag", nf90_short,(/nxID,nyID,nzID/),arscl_cloud_detection_flagID)
status = nf90_def_var(ncid, "arscl_cloud_layer_base_height", nf90_float,(/nxID,nyID,nlayersID/),arscl_cloud_base_heightID)
status = nf90_def_var(ncid, "arscl_cloud_layer_top_height", nf90_float,(/nxID,nyID,nlayersID/),arscl_cloud_top_heightID)
endif
endif
!---------------------------
! MWR LWP added by oue 2017.03.25
if ((conf%mwrID==1) .and. (iht==0)) then
status = nf90_def_var(ncid, "model_lwp", nf90_float,(/nxID,nyID/),true_lwpID)
status = nf90_def_var(ncid, "mwr_lwp", nf90_float,(/nxID,nyID/),mwr_lwpID)
status = nf90_def_var(ncid, "number_of_gridpoints_mwrlwp", nf90_int,(/nxID,nyID,nzID/),mwr_samplesID)
endif
!---------------------------
!--------------------------- 

!
if (iht==0) then
status = nf90_def_var(ncid, "Zmin",  nf90_float, (/nxID,nyID,nzID/),ZminID)
status = nf90_def_var(ncid, "temp",  nf90_float, (/nxID,nyID,nzID/),tempID)
status = nf90_def_var(ncid, "rho_d",  nf90_float, (/nxID,nyID,nzID/),rho_dID)
status = nf90_def_var(ncid, "u",     nf90_float, (/nxID,nyID,nzID/),uID)
status = nf90_def_var(ncid, "v",     nf90_float, (/nxID,nyID,nzID/),vID)
status = nf90_def_var(ncid, "w",     nf90_float, (/nxID,nyID,nzID/),wID)
status = nf90_def_var(ncid, "Avap", nf90_float, (/nxID,nyID,nzID/),AvapID)
If(dabs(conf%elev)==90) Then
status = nf90_def_var(ncid, "Atot", nf90_float, (/nxID,nyID,nzID/),AtotID)
endif
endif
!
status = nf90_def_var(ncid, "rad_freq",  nf90_float, (/oneID/),rad_freqID)
status = nf90_def_var(ncid, "rad_beamwidth", nf90_float, (/oneID/),rad_beamwidthID)
status = nf90_def_var(ncid, "rad_range_resolution",nf90_float,(/oneID/),rad_range_resolutionID)
status = nf90_def_var(ncid, "rad_ixc",nf90_int,(/oneID/),ixcID)
status = nf90_def_var(ncid, "rad_iyc",nf90_int,(/oneID/),iycID)
status = nf90_def_var(ncid, "rad_zc",nf90_float,(/oneID/),zcID)
!
status = nf90_def_var(ncid, "elev",  nf90_float, (/nxID,nyID,nzID/),elevID)
status = nf90_def_var(ncid, "azim",  nf90_float, (/nxID,nyID,nzID/),azimID)
status = nf90_def_var(ncid, "range",  nf90_float, (/nxID,nyID,nzID/),rangeID)
status = nf90_def_var(ncid, "height",nf90_float, (/nxID,nyID,nzID/),zID)
status = nf90_def_var(ncid, "wcont", nf90_float, (/nxID,nyID,nzID/),qhydroID)
!
status = nf90_def_var(ncid, "x_scene",  nf90_float, (/nxID/),xID)
status = nf90_def_var(ncid, "y_scene",  nf90_float, (/nyID/),yID)
!

!
! Attributes
!

    status = nf90_put_att(ncid,xlatID,'description','Latitude (North is positive)')
    status = nf90_put_att(ncid,xlatID,'units','deg')
    status = nf90_put_att(ncid,xlatID,'_FillValue',-999.e0)
!
    status = nf90_put_att(ncid,xlongID,'description','Longitude (West is positive)')
    status = nf90_put_att(ncid,xlongID,'units','deg')
    status = nf90_put_att(ncid,xlongID,'_FillValue',-999.e0)
!
    status = nf90_put_att(ncid,ZhhID,'description','Reflectivity at hh polarization')
    status = nf90_put_att(ncid,ZhhID,'units','dBZ')
    status = nf90_put_att(ncid,ZhhID,'_FillValue',-999.e0)
!
    status = nf90_put_att(ncid,ZvvID,'description','Reflectivity at vv polarization')
    status = nf90_put_att(ncid,ZvvID,'units','dBZ')
    status = nf90_put_att(ncid,ZvvID,'_FillValue',-999.e0)
!
    status = nf90_put_att(ncid,ZvhID,'description','Reflectivity at vh polarization')
    status = nf90_put_att(ncid,ZvhID,'units','dBZ')
    status = nf90_put_att(ncid,ZvhID,'_FillValue',-999.e0)
!
    status = nf90_put_att(ncid,ZdrID,'description','Differential reflectivity')
    status = nf90_put_att(ncid,ZdrID,'units','dB')
    status = nf90_put_att(ncid,ZdrID,'_FillValue',-999.e0)
!
    status = nf90_put_att(ncid,LDRhID,'description','Linear depolarization ratio')
    status = nf90_put_att(ncid,LDRhID,'units','dB')
    status = nf90_put_att(ncid,LDRhID,'_FillValue',-999.e0)
!
    status = nf90_put_att(ncid,RHOhvID,'description','Cross-correlation coefficient')
    status = nf90_put_att(ncid,RHOhvID,'units','-')
    status = nf90_put_att(ncid,RHOhvID,'_FillValue',-999.e0)
!
    status = nf90_put_att(ncid,DoppID,'description','Radial Doppler velocity, positive upward')
    status = nf90_put_att(ncid,DoppID,'units','m/s')
    status = nf90_put_att(ncid,DoppID,'_FillValue',-999.e0)
!
    status = nf90_put_att(ncid,SWhID,'description','Spectrum width due to hydrometeors')
    status = nf90_put_att(ncid,SWhID,'units','m/s')
    status = nf90_put_att(ncid,SWhID,'_FillValue',-999.e0)
!
    status = nf90_put_att(ncid,SWtID,'description','Spectrum width due to turbulence')
    status = nf90_put_att(ncid,SWtID,'units','m/s')
    status = nf90_put_att(ncid,SWtID,'_FillValue',-999.e0)
!
    status = nf90_put_att(ncid,SWsID,'description','Spectrum width due to wind shear')
    status = nf90_put_att(ncid,SWsID,'units','m/s')
    status = nf90_put_att(ncid,SWsID,'_FillValue',-999.e0)
!
    status = nf90_put_att(ncid,SWvID,'description','Spectrum width due to cross wind')
    status = nf90_put_att(ncid,SWvID,'units','m/s')
    status = nf90_put_att(ncid,SWvID,'_FillValue',-999.e0)
!
    status = nf90_put_att(ncid,SWtotID,'description','Total spectrum width')
    status = nf90_put_att(ncid,SWtotID,'units','m/s')
    status = nf90_put_att(ncid,SWtotID,'_FillValue',-999.e0)
!
    status = nf90_put_att(ncid,Dopp90ID,'description','Vertical Doppler velocity (el=90), positive upward')
    status = nf90_put_att(ncid,Dopp90ID,'units','m/s')
    status = nf90_put_att(ncid,Dopp90ID,'_FillValue',-999.e0)
!
    status = nf90_put_att(ncid,SWh90ID,'description','Spectrum width due to hydrometeors, elevation=90')
    status = nf90_put_att(ncid,SWh90ID,'units','m/s')
    status = nf90_put_att(ncid,SWh90ID,'_FillValue',-999.e0)
!
    status = nf90_put_att(ncid,DVhID,'description','Reflectivity weighted velocity (el=90deg, w=0m/s)')
    status = nf90_put_att(ncid,DVhID,'units','m/s')
    status = nf90_put_att(ncid,DVhID,'_FillValue',-999.e0)
!   
!
    status = nf90_put_att(ncid,KdpID,'description','Specific differencial phase')
    status = nf90_put_att(ncid,KdpID,'units','deg/km')
    status = nf90_put_att(ncid,KdpID,'_FillValue',-999.e0)
!
    status = nf90_put_att(ncid,AdpID,'description','Differencial attenuation')
    status = nf90_put_att(ncid,AdpID,'units','dB/km')
    status = nf90_put_att(ncid,AdpID,'_FillValue',-999.e0)
!
    status = nf90_put_att(ncid,AhID,'description','Specific horizontal attenuation')
    status = nf90_put_att(ncid,AhID,'units','dB/km')
    status = nf90_put_att(ncid,AhID,'_FillValue',-999.e0)
!
    status = nf90_put_att(ncid,AvID,'description','Specific vertical attenuation')
    status = nf90_put_att(ncid,AvID,'units','dB/km')
    status = nf90_put_att(ncid,AvID,'_FillValue',-999.e0)
!
    status = nf90_put_att(ncid,diff_back_phaseID,'description','Differential backscatter phase')
    status = nf90_put_att(ncid,diff_back_phaseID,'units','deg')
    status = nf90_put_att(ncid,diff_back_phaseID,'_FillValue',-999.e0)
!
!------------------------------------------
! CEILO
IF  ((conf%ceiloID==1) .and. (iht==0)) THEN
!
    status = nf90_put_att(ncid,ceilo_back_obsID,'description','Ceilo observed backscatter')
    status = nf90_put_att(ncid,ceilo_back_obsID,'units','[m sr]^-1')
    status = nf90_put_att(ncid,ceilo_back_obsID,'_FillValue',-999.e0)
!
    status = nf90_put_att(ncid,ceilo_back_trueID,'description','Ceilo true backscatter')
    status = nf90_put_att(ncid,ceilo_back_trueID,'units','[m sr]^-1')
    status = nf90_put_att(ncid,ceilo_back_trueID,'_FillValue',-999.e0)
!
    status = nf90_put_att(ncid,ceilo_extID,'description','Ceilo extinction coefficient')
    status = nf90_put_att(ncid,ceilo_extID,'units','[m]^-1')
    status = nf90_put_att(ncid,ceilo_extID,'_FillValue',-999.e0)
!
    status = nf90_put_att(ncid,lidar_ratioID,'description','Ceilo lidar ratio')
    status = nf90_put_att(ncid,lidar_ratioID,'units','sr')
    status = nf90_put_att(ncid,lidar_ratioID,'_FillValue',-999.e0)
!
    status = nf90_put_att(ncid,ceilo_first_cloud_baseID,'description','Height of the first cloud base')
    status = nf90_put_att(ncid,ceilo_first_cloud_baseID,'units','m')
    status = nf90_put_att(ncid,ceilo_first_cloud_baseID,'_FillValue',-999.e0)
ENDIF
!------------------------------------------
! MPL
IF (conf%mplID>0) THEN
IF ( (iht==0) .or. (iht==1) .or. (iht==3).or.(iht==7) ) THEN
!
    status = nf90_put_att(ncid,mpl_wavelID,'description','MPL wavelength')
    status = nf90_put_att(ncid,mpl_wavelID,'units','nm')
!
    if (iht==0) status = nf90_put_att(ncid,mpl_back_obsID,'description','MPL observed backscatter (cloud+ice+aerosol+atm.moleculs)')
    if ((iht==1) .or.(iht==3).or.(iht==7))  status = nf90_put_att(ncid,mpl_back_obsID,&
                               'description','MPL observed backscatter (hydrometeor only)')
    status = nf90_put_att(ncid,mpl_back_obsID,'units','[m sr]^-1')
    status = nf90_put_att(ncid,mpl_back_obsID,'_FillValue',-999.e0)
!
    if (iht==0) status = nf90_put_att(ncid,mpl_back_trueID,'description','MPL true backscatter (cloud+ice+aerosol+atm.moleculs)')
    if ((iht==1) .or.(iht==3).or.(iht==7))  status = nf90_put_att(ncid,mpl_back_trueID,&
         'description','MPL true backscatter (hydrometeor only)')
    status = nf90_put_att(ncid,mpl_back_trueID,'units','[m sr]^-1')
    status = nf90_put_att(ncid,mpl_back_trueID,'_FillValue',-999.e0)
!
    if (iht==0) status = nf90_put_att(ncid,mpl_rayleigh_backID,'description','MPL molecular backscatter)')
    status = nf90_put_att(ncid,mpl_rayleigh_backID,'units','[m sr]^-1')
    status = nf90_put_att(ncid,mpl_rayleigh_backID,'_FillValue',-999.e0)
!
    if (iht==0) status =  nf90_put_att(ncid,mpl_extID,'description','MPL extinction coefficient (cloud+ice+aerosol+atm.moleculs)')
    if ((iht==1) .or.(iht==3).or.(iht==7)) status = nf90_put_att(ncid,mpl_extID,&
         'description','MPL extinction coefficient (hydrometeor only)')
    status = nf90_put_att(ncid,mpl_extID,'units','[m]^-1')
    status = nf90_put_att(ncid,mpl_extID,'_FillValue',-999.e0)
!
    if (iht==0) status = nf90_put_att(ncid,mpl_lidar_ratioID,'description','MPL Lidar ratio (cloud+ice+aerosol+atm.moleculs)')
    if ((iht==1) .or.(iht==3).or.(iht==7)) status = nf90_put_att(ncid,mpl_lidar_ratioID,&
         'description','MPL Lidar ratio (hydrometeor only)')
    status = nf90_put_att(ncid,mpl_lidar_ratioID,'units','sr')
    status = nf90_put_att(ncid,mpl_lidar_ratioID,'_FillValue',-999.e0)
!
ENDIF
!
IF (conf%aeroID==1) THEN
IF (iht==0) THEN
!
    status = nf90_put_att(ncid,aero_back_obsID,'description','MPL aerosol observed backscatter')
    status = nf90_put_att(ncid,aero_back_obsID,'units','[m sr]^-1')
    status = nf90_put_att(ncid,aero_back_obsID,'_FillValue',-999.e0)
!   
    status = nf90_put_att(ncid,aero_back_trueID,'description','MPL aerosol true backscatter')
    status = nf90_put_att(ncid,aero_back_trueID,'units','[m sr]^-1')
    status = nf90_put_att(ncid,aero_back_trueID,'_FillValue',-999.e0)
!   
    status = nf90_put_att(ncid,aero_extID,'description','MPL aerosol extinction coefficient')
    status = nf90_put_att(ncid,aero_extID,'units','[m]^-1')
    status = nf90_put_att(ncid,aero_extID,'_FillValue',-999.e0)
!   
    status = nf90_put_att(ncid,aero_lidar_ratioID,'description','MPL aerosol Lidar ratio')
    status = nf90_put_att(ncid,aero_lidar_ratioID,'units','sr')
    status = nf90_put_att(ncid,aero_lidar_ratioID,'_FillValue',-999.e0)
!
ENDIF !  iht==0
ENDIF !  conf%aeroID==1
!
ENDIF !conf%mplID>0

!------------------------------------------
! Doppler spectra by oue in May 2018
IF (conf%spectraID>0) THEN 
    status = nf90_put_att(ncid,velbinsID,'description','Doppler spectrum velocity')
    status = nf90_put_att(ncid,velbinsID,'units','[m s]^-1')
    status = nf90_put_att(ncid,velbinsID,'_FillValue',-999.e0)

    status = nf90_put_att(ncid,spectra_zhhID,'description','copolar Doppler spectra for horizontal polarization')
    status = nf90_put_att(ncid,spectra_zhhID,'units','dB')
    status = nf90_put_att(ncid,spectra_zhhID,'_FillValue',-999.e0)

    status = nf90_put_att(ncid,spectra_zvhID,'description','crosspolar Doppler spectra for horizontal polarization')
    status = nf90_put_att(ncid,spectra_zvhID,'units','dB')
    status = nf90_put_att(ncid,spectra_zvhID,'_FillValue',-999.e0)

    status = nf90_put_att(ncid,spectra_zvvID,'description','copolar Doppler spectra for vertical polarization')
    status = nf90_put_att(ncid,spectra_zvvID,'units','dB')
    status = nf90_put_att(ncid,spectra_zvvID,'_FillValue',-999.e0)
ENDIF

!------------------------------------------
! Airborne by oue in Oct 2023
IF ((conf%airborne==1) .and. (iht==0)) THEN 
    status = nf90_put_att(ncid,Zsfc_oceanID,'description','Surface reflectivity over ocean')
    status = nf90_put_att(ncid,Zsfc_oceanID,'units','dBZ')
    status = nf90_put_att(ncid,Zsfc_oceanID,'_FillValue',-999.e0)

    status = nf90_put_att(ncid,Zsfc_land_coarseID,'description','Surface reflectivity over coarse land (i.e. wet snow)')
    status = nf90_put_att(ncid,Zsfc_land_coarseID,'units','dBZ')
    status = nf90_put_att(ncid,Zsfc_land_coarseID,'_FillValue',-999.e0)

    status = nf90_put_att(ncid,Zsfc_land_flatID,'description','Surface reflectivity over flat land (i.e. flat road)')
    status = nf90_put_att(ncid,Zsfc_land_flatID,'units','dBZ')
    status = nf90_put_att(ncid,Zsfc_land_flatID,'_FillValue',-999.e0)

    status = nf90_put_att(ncid,DV_airborneID,'description','Radial airmotion related to airborne moving')
    status = nf90_put_att(ncid,DV_airborneID,'units','[m s]^-1')
    status = nf90_put_att(ncid,DV_airborneID,'_FillValue',-999.e0)
ENDIF

!---------------------------
! ARSCL attributes added by oue 2017.03.23
if ((conf%arsclID==1) .and. (conf%mplID>0) .and. (conf%ceiloID==1)) then 
if (iht==0) then
    status = nf90_put_att(ncid,arscl_cloud_maskID,'description','Cloud mask from radar, mpl, ceilometer observations')
    status = nf90_put_att(ncid,arscl_cloud_maskID,'units','unitless')
    status = nf90_put_att(ncid,arscl_cloud_maskID,'flag_0_description','clear')
    status = nf90_put_att(ncid,arscl_cloud_maskID,'flag_1_description','cloudy')
!
    status = nf90_put_att(ncid,arscl_cloud_detection_flagID,'description','Instrument source flag for cloud detections')
    status = nf90_put_att(ncid,arscl_cloud_detection_flagID,'units','unitless')
    status = nf90_put_att(ncid,arscl_cloud_detection_flagID,'_FillValue',int2(-999))
    status = nf90_put_att(ncid,arscl_cloud_detection_flagID,'flag_values','1s, 2s, 3s, 4s')
    status = nf90_put_att(ncid,arscl_cloud_detection_flagID,'flag_1_description','Clear according to radar and MPL')
    status = nf90_put_att(ncid,arscl_cloud_detection_flagID,'flag_2_description','Cloud detected by radar and MPL')
    status = nf90_put_att(ncid,arscl_cloud_detection_flagID,'flag_3_description','Cloud detected by radar only')
    status = nf90_put_att(ncid,arscl_cloud_detection_flagID,'flag_4_description','Cloud detected by MPL only')
!
    status = nf90_put_att(ncid,arscl_cloud_base_heightID,'description','Base height of cloudy layers for up to 10 layers')
    status =  nf90_put_att(ncid,arscl_cloud_base_heightID,'comments','based on combined radar, ceilometer, and mpl observations')
    status = nf90_put_att(ncid,arscl_cloud_base_heightID,'units','m')
    status = nf90_put_att(ncid,arscl_cloud_base_heightID,'_FillValue',-999.e0)
!
    status = nf90_put_att(ncid,arscl_cloud_top_heightID,'description','Top height of cloudy layers for up to 10 layers')
    status = nf90_put_att(ncid,arscl_cloud_top_heightID,'comments','based on combined radar, ceilometer, and mpl observations')
    status = nf90_put_att(ncid,arscl_cloud_top_heightID,'units','m')
    status = nf90_put_att(ncid,arscl_cloud_top_heightID,'_FillValue',-999.e0)
endif
endif !conf%arsclID==1
!---------------------------
! MWR LWP attributes added by oue 2017.03.25
if ((conf%mwrID==1) .and. (iht==0)) then
    status = nf90_put_att(ncid,true_lwpID,'description','Model liquid water path (cloud+rain)')
    status = nf90_put_att(ncid,true_lwpID,'units','kg m^-2')

    status = nf90_put_att(ncid,mwr_lwpID,'description','Microwave radiometer liquid water path taking account of field of view')
    status = nf90_put_att(ncid,mwr_lwpID,'units','kg m^-2')
    status = nf90_put_att(ncid,mwr_lwpID,'field_of_view_deg',conf%mwr_view)
    status = nf90_put_att(ncid,mwr_lwpID,'altitude_of_mwr',conf%mwr_alt)

    status = nf90_put_att(ncid,mwr_samplesID,'description','Number of gridpoints within MWR field of view')
    status = nf90_put_att(ncid,mwr_samplesID,'units','unitless')
endif !conf%mwrID==1
!---------------------------
!---------------------------



!
if (iht==0) Then
!
    status = nf90_put_att(ncid,ZminID,'description','Radar sensitivity limitation with range')
    status = nf90_put_att(ncid,ZminID,'units','dBZ')
    status = nf90_put_att(ncid,ZminID,'_FillValue',-999.e0)
!
    status = nf90_put_att(ncid,tempID,'description','Temperature')
    status = nf90_put_att(ncid,tempID,'units','C')
!
    status = nf90_put_att(ncid,rho_dID,'description','Dry air density')
    status = nf90_put_att(ncid,rho_dID,'units','kg m^-3')
!
    status = nf90_put_att(ncid,uID,'description','U horizontal wind component')
    status = nf90_put_att(ncid,uID,'units','m/s')
!
    status = nf90_put_att(ncid,vID,'description','V horizontal wind component')
    status = nf90_put_att(ncid,vID,'units','m/s')
!
    status = nf90_put_att(ncid,wID,'description','Vertical air velocity')
    status = nf90_put_att(ncid,wID,'units','m/s')
!
    status = nf90_put_att(ncid,AvapID,'description','Specific water vapor attenuation')
    status = nf90_put_att(ncid,AvapID,'units','dB/km')
    status = nf90_put_att(ncid,AvapID,'_FillValue',-999.e0)
!
If(dabs(conf%elev)==90) Then
    status = nf90_put_att(ncid,AtotID,'description','Two-way total attenuation (Ah+Avap)')
    status = nf90_put_att(ncid,AtotID,'units','dB')
    status = nf90_put_att(ncid,AtotID,'_FillValue',-999.e0)
endif
Endif
!------------------------------------------
!------------------------------------------
!
    status = nf90_put_att(ncid,elevID,'description','Elevation')
    status = nf90_put_att(ncid,elevID,'units','degrees')
    status = nf90_put_att(ncid,elevID,'_FillValue',-999.e0)
!
    status = nf90_put_att(ncid,azimID,'description','Azimuth, E=0, W=180, N=90, S=270 degs')
    status = nf90_put_att(ncid,azimID,'units','degrees')
    status = nf90_put_att(ncid,azimID,'_FillValue',-999.e0)
!
    status = nf90_put_att(ncid,rangeID,'description','Radar range')
    status = nf90_put_att(ncid,rangeID,'units','m')
    status = nf90_put_att(ncid,rangeID,'_FillValue',-999.e0)
!
    status = nf90_put_att(ncid,zID,'description','Height')
    status = nf90_put_att(ncid,zID,'units','m')
!
    status = nf90_put_att(ncid,qhydroID,'description','Water content')
    status = nf90_put_att(ncid,qhydroID,'units','kg/m^3')
    status = nf90_put_att(ncid,qhydroID,'_FillValue',-999.e0)
!
    status = nf90_put_att(ncid,rad_freqID,'description','radar frequency')
    status = nf90_put_att(ncid,rad_freqID,'long_name','configuration parameter')
    status = nf90_put_att(ncid,rad_freqID,'units','GHz')
!
    status = nf90_put_att(ncid,rad_beamwidthID,'description','radar antenna beamwidth')
    status = nf90_put_att(ncid,rad_beamwidthID,'long_name','configuration parameter')
    status = nf90_put_att(ncid,rad_beamwidthID,'units','deg')
!
    status = nf90_put_att(ncid,rad_range_resolutionID,'description','radar range resolution') 
    status = nf90_put_att(ncid,rad_range_resolutionID,'long_name','configuration parameter')
    status = nf90_put_att(ncid,rad_range_resolutionID,'units','m')
!
    status = nf90_put_att(ncid,ixcID,'description','index of position in x direction of radar origin')
    status = nf90_put_att(ncid,ixcID,'long_name','configuration parameter')
    status = nf90_put_att(ncid,ixcID,'units','-')
!
    status = nf90_put_att(ncid,iycID,'description','index of position in y direction of radar origin')
    status = nf90_put_att(ncid,iycID,'long_name','configuration parameter')
    status = nf90_put_att(ncid,iycID,'units','-')
!
    status = nf90_put_att(ncid,zcID,'description','radar height')
    status = nf90_put_att(ncid,zcID,'long_name','configuration parameter')
    status = nf90_put_att(ncid,zcID,'units','m')
!
    status = nf90_put_att(ncid,xID,'description','Scene extent in E-W direction')
    status = nf90_put_att(ncid,xID,'units','m')
!
    status = nf90_put_att(ncid,yID,'description','Scene extent in N-S direction')
    status = nf90_put_att(ncid,yID,'units','m')
!
!---------------------------------------------------------------------------------------
! GLOBAL ATTRIBUTES

   att_string="forward scanning radar simulator output" 
   status = nf90_put_att(ncid, NF90_GLOBAL,'description',att_string) 
   !
   status = nf90_put_att(ncid, NF90_GLOBAL,'model_version','crsim_v4.0')
   !
   att_string=Trim(conf%WRFInputFile)
   status = nf90_put_att(ncid, NF90_GLOBAL,'WRF_input_file',att_string) 
   !
   write( str1, '(i5)' ) conf%MP_PHYSICS
   att_string=Trim(str1)
   status = nf90_put_att(ncid, NF90_GLOBAL,'MP_PHYSICS',att_string)
   !
   if (Min(conf%ix_start,conf%ix_end)<0) then
   j1=1 ; j2=rmout%nx
   else
   j1=conf%ix_start ; j2=conf%ix_end
   endif
   write( str1, '(i10)' ) j1
   write( str2, '(i10)' ) j2
   att_string=Trim(Adjustl(str1))//" - "//Trim(Adjustl(str2))
   status = nf90_put_att(ncid, NF90_GLOBAL,'x_indices_of_WRF_extracted_scene',att_string)
   !
   if (Min(conf%iy_start,conf%iy_end)<0) then
   j1=1 ; j2=rmout%ny
   else
   j1=conf%iy_start ; j2=conf%iy_end
   endif
   write( str1, '(i10)' ) j1
   write( str2, '(i10)' ) j2
   att_string=Trim(Adjustl(str1))//" - "//Trim(Adjustl(str2))
   status = nf90_put_att(ncid, NF90_GLOBAL,'y_indices_of_WRF_extracted_scene',att_string)
   !
   if (Min(conf%iz_start,conf%iz_end)<0) then
   j1=1 ; j2=rmout%nz
   else
   j1=conf%iz_start ; j2=conf%iz_end
   endif
   write( str1, '(i10)' ) j1
   write( str2, '(i10)' ) j2
   att_string=Trim(Adjustl(str1))//" - "//Trim(Adjustl(str2))
   status = nf90_put_att(ncid, NF90_GLOBAL,'z_indices_of_WRF_extracted_scene',att_string)
   !
   write( str1, '(i10)' ) conf%it
   att_string=Trim(Adjustl(str1))
   status = nf90_put_att(ncid, NF90_GLOBAL,'scene_extracted_at_time_step',att_string)
   !
   write( str1, '(f10.1)' ) conf%freq
   att_string=Trim(Adjustl(str1))//" GHz"
   status = nf90_put_att(ncid, NF90_GLOBAL,'radar_frequency',att_string)
   !
   write( str1, '(i10)' ) conf%ixc
   write( str2, '(i10)' ) conf%iyc
   att_string=Trim(Adjustl(str1))//" , "//Trim(Adjustl(str2))
   status = nf90_put_att(ncid,NF90_GLOBAL,'x_and_y_indices_of_radar_position',att_string)
   !
   write( str1, '(f10.1)' ) conf%zc
   att_string=Trim(Adjustl(str1))//" m"
   status = nf90_put_att(ncid,NF90_GLOBAL,'height_of_radar',att_string)
   !
   if (conf%elev<0.d0) then
   att_string="elevation of each scene pixel is relative to the radar origin"
   else
   write( str1, '(f5.1)' ) conf%elev
   att_string="elevation is set to "//Trim(Adjustl(str1))//" degrees"
   endif
   status = nf90_put_att(ncid,NF90_GLOBAL,'scanning_mode',att_string)
   !
   !------------------------------------------
   ! Doppler spectra by oue in May 2018
   if(conf%spectraID==1) then
   status = nf90_put_att(ncid, NF90_GLOBAL,'PRF_for_Doppler_spectra',spectra%PRF)
   status = nf90_put_att(ncid, NF90_GLOBAL,'TimeSampling_for_Doppler_spectra',spectra%TimeSampling)
   status = nf90_put_att(ncid, NF90_GLOBAL,'Noise1km_for_Doppler_spectra',spectra%NOISE_1km)
   status = nf90_put_att(ncid, NF90_GLOBAL,'Nyquist_velocity_for_Doppler_spectra',spectra%VNyquist)
   status = nf90_put_att(ncid, NF90_GLOBAL,'Nave_for_Doppler_spectra',spectra%Nave)
   end if
   !------------------------------------------
   ! Airborne by oue in Oct 2023
   if((conf%airborne==1) .and. (iht==0)) then
   status = nf90_put_att(ncid, NF90_GLOBAL,'Pulse_length_m',conf%pulse_len)
   status = nf90_put_att(ncid, NF90_GLOBAL,'Airborne_moving_speed_mps',conf%airborne_spd)
   status = nf90_put_att(ncid, NF90_GLOBAL,'Airborne_moving_azimuth_direction',conf%airborne_azdeg)
   status = nf90_put_att(ncid, NF90_GLOBAL,'Airborne_moving_elevation_direction',conf%airborne_eldeg)
   end if
   !------------------------------------------
   !
   status = nf90_put_att(ncid, NF90_GLOBAL,'created_by','Aleksandra Tatarevic, Mariko Oue')
   status = nf90_put_att(ncid,NF90_GLOBAL,'institute','Stony Brook University')
   status = nf90_put_att(ncid,NF90_GLOBAL,'websites',&
        'https://www.bnl.gov/CMAS/cr-sim.php; http://radarscience.weebly.com/radar-simulators.html')
!---------------------------------------------------------------------------------------

status = nf90_enddef(ncid)
If (status /= 0) then ; err_msg='Error in nf90_enddef'; Goto 999 ; endif

!------------------------------------------------------------
IF (iht==0) THEN
! write variables
!
status =  nf90_put_var(ncid,xlatID,real(env%xlat))
status =  nf90_put_var(ncid,xlongID,real(env%xlong))
!
status =  nf90_put_var(ncid,ZhhID,real(rmout%Zhh_tot))
status =  nf90_put_var(ncid,ZvvID,real(mout%Zvv_tot))
status =  nf90_put_var(ncid,ZvhID,real(mout%Zvh_tot))
status =  nf90_put_var(ncid,ZdrID,real(mout%Zdr_tot))
status =  nf90_put_var(ncid,LDRhID,real(mout%LDRh_tot))
status =  nf90_put_var(ncid,RHOhvID,real(mout%RHOhv_tot))
status =  nf90_put_var(ncid,DoppID,real(rmout%Dopp_tot))
status =  nf90_put_var(ncid,SWhID,real(rmout%dDVh_tot))
status =  nf90_put_var(ncid,SWtID,real(rmout%sw_t))
status =  nf90_put_var(ncid,SWsID,real(rmout%sw_s))
status =  nf90_put_var(ncid,SWvID,real(rmout%sw_v))
status =  nf90_put_var(ncid,SWtotID,real(rmout%SWt_tot))

status =  nf90_put_var(ncid,Dopp90ID,real(rmout%Dopp90_tot))
status =  nf90_put_var(ncid,SWh90ID,real(rmout%dDVh90_tot))
status =  nf90_put_var(ncid,DVhID,real(rmout%DVh_tot))
!
status =  nf90_put_var(ncid,KdpID,real(mout%Kdp_tot))
status =  nf90_put_var(ncid,AdpID,real(mout%Adp_tot))
status =  nf90_put_var(ncid,AhID,real(rmout%Ah_tot))
status =  nf90_put_var(ncid,AvID,real(mout%Av_tot))
status =  nf90_put_var(ncid,diff_back_phaseID,real(mout%diff_back_phase_tot))
!
!---------------------------------------
! CEILO modified by oue for RAMS
IF (conf%ceiloID==1)THEN 
status =  nf90_put_var(ncid,ceilo_back_obsID,real(lout%ceilo_back_obs_tot))
status =  nf90_put_var(ncid,ceilo_back_trueID,real(lout%ceilo_back_true_tot))
status =  nf90_put_var(ncid,ceilo_extID,real(lout%ceilo_ext_tot))
status =  nf90_put_var(ncid,lidar_ratioID,real(lout%lidar_ratio_tot))
status =  nf90_put_var(ncid,ceilo_first_cloud_baseID,real(lout%ceilo_first_cloud_base))
ENDIF
!---------------------------------------
! MPL
IF (conf%mplID>0) THEN
allocate(work(one)); work(1)=real(mpl%wavel)
status =  nf90_put_var(ncid,mpl_wavelID,work)
deallocate(work)
status =  nf90_put_var(ncid,mpl_back_obsID,real(mpl%back_obs_tot))
status =  nf90_put_var(ncid,mpl_back_trueID,real(mpl%back_true_tot))
status =  nf90_put_var(ncid,mpl_rayleigh_backID,real(mpl%rayleigh_back))
status =  nf90_put_var(ncid,mpl_extID,real(mpl%ext_tot))
status =  nf90_put_var(ncid,mpl_lidar_ratioID,real(mpl%lidar_ratio_tot))
!! aerosol
IF (conf%aeroID==1) THEN
status =  nf90_put_var(ncid,aero_back_obsID,real(aero%back_obs))
status =  nf90_put_var(ncid,aero_back_trueID,real(aero%back_true))
status =  nf90_put_var(ncid,aero_extID,real(aero%ext))
status =  nf90_put_var(ncid,aero_lidar_ratioID,real(aero%lidar_ratio))
ENDIF ! if (conf%aeroID==1)
ENDIF ! if (conf%mplID>0)
!---------------------------
!------------------------------------------
! Doppler spectra by oue in May 2018
if(conf%spectraID==1) then
status =  nf90_put_var(ncid,velbinsID,    real(spectra%vel_bins))
status =  nf90_put_var(ncid,spectra_zhhID,real(spectra%zhh_spectra_tot))
status =  nf90_put_var(ncid,spectra_zvhID,real(spectra%zvh_spectra_tot))
status =  nf90_put_var(ncid,spectra_zvvID,real(spectra%zvv_spectra_tot))
end if
!---------------------------
! Airborne added by oue Oct 2023
if((conf%airborne==1) .and. (iht==0)) then
status =  nf90_put_var(ncid,Zsfc_oceanID,      real(airborne%Zsfc_ocean))
status =  nf90_put_var(ncid,Zsfc_land_coarseID,real(airborne%Zsfc_land_coarse))
status =  nf90_put_var(ncid,Zsfc_land_flatID,  real(airborne%Zsfc_land_flat))
status =  nf90_put_var(ncid,DV_airborneID,     real(airborne%Dopp_airborne))
end if
!---------------------------
! ARSCL  added by oue 2017.03.23
if ((conf%arsclID==1) .and. (conf%mplID>0) .and. (conf%ceiloID==1)) then 
status =  nf90_put_var(ncid,arscl_cloud_maskID,arscl%cloud_mask)
status =  nf90_put_var(ncid,arscl_cloud_detection_flagID,arscl%cloud_detection_flag)
status =  nf90_put_var(ncid,arscl_cloud_base_heightID,real(arscl%cloud_layer_base))
status =  nf90_put_var(ncid,arscl_cloud_top_heightID,real(arscl%cloud_layer_top))
endif ! if (conf%arsclID==1)
!--------------------------- 
! MWR LWP added by oue 2017.03.25 
if (conf%mwrID==1) then
status =  nf90_put_var(ncid,true_lwpID,mwr%true_lwp)
status =  nf90_put_var(ncid,mwr_lwpID,mwr%mwr_lwp)
status =  nf90_put_var(ncid,mwr_samplesID,mwr%n_samples)
endif ! if (conf%mwrID==1)   
!---------------------------------------
!---------------------------------------

status =  nf90_put_var(ncid,ZminID, real(conf%ZMIN + 20.d0 *  dlog10(1.d-3*rmout%range)) ) 
status =  nf90_put_var(ncid,tempID,real(env%temp))
status =  nf90_put_var(ncid,rho_dID,real(env%rho_d))
status =  nf90_put_var(ncid,uID,real(env%u))
status =  nf90_put_var(ncid,vID,real(env%v))
status =  nf90_put_var(ncid,wID,real(env%w))
status =  nf90_put_var(ncid,AvapID,real(rmout%Avap))
If(dabs(conf%elev)==90) Then
status =  nf90_put_var(ncid,AtotID,real(rmout%Atot))
endif
!
status =  nf90_put_var(ncid,elevID,real(rmout%elev))
status =  nf90_put_var(ncid,azimID,real(rmout%azim))
status =  nf90_put_var(ncid,rangeID,real(rmout%range))
status =  nf90_put_var(ncid,zID,real(env%z))
status =  nf90_put_var(ncid,qhydroID,real(qq)) ! kg/m3
!
allocate(work(one),iwork(one)) 
work(1)=real(conf%freq)
status =  nf90_put_var(ncid,rad_freqID,work)
work(1)=real(conf%Theta1)
status =  nf90_put_var(ncid,rad_beamwidthID,work)
work(1)=real(conf%dr)
status =  nf90_put_var(ncid,rad_range_resolutionID,work)
iwork(1)=conf%ixc
status =  nf90_put_var(ncid,ixcID,iwork)
iwork(1)=conf%iyc
status =  nf90_put_var(ncid,iycID,iwork)
work(1)=real(conf%zc)
status =  nf90_put_var(ncid,zcID,work)
deallocate(work,iwork)
!
status =  nf90_put_var(ncid,xID,real(env%x))
status =  nf90_put_var(ncid,yID,real(env%y))
!------------------------------------------------------------
ELSE ! HYDROMETEORS
! write variables

status =  nf90_put_var(ncid,xlatID,real(env%xlat))
status =  nf90_put_var(ncid,xlongID,real(env%xlong))
!
status =  nf90_put_var(ncid,ZhhID,real(rmout%Zhh(:,:,:,iht)))
status =  nf90_put_var(ncid,ZvvID,real(mout%Zvv(:,:,:,iht)))
status =  nf90_put_var(ncid,ZvhID,real(mout%Zvh(:,:,:,iht)))
status =  nf90_put_var(ncid,ZdrID,real(mout%Zdr(:,:,:,iht)))
status =  nf90_put_var(ncid,LDRhID,real(mout%LDRh(:,:,:,iht)))
status =  nf90_put_var(ncid,RHOhvID,real(mout%RHOhv(:,:,:,iht)))
status =  nf90_put_var(ncid,DoppID,real(rmout%Dopp(:,:,:,iht)))
status =  nf90_put_var(ncid,SWhID,real(rmout%dDVh(:,:,:,iht)))
status =  nf90_put_var(ncid,SWtID,real(rmout%sw_t))
status =  nf90_put_var(ncid,SWsID,real(rmout%sw_s))
status =  nf90_put_var(ncid,SWvID,real(rmout%sw_v))
status =  nf90_put_var(ncid,SWtotID,real(rmout%SWt(:,:,:,iht)))
status =  nf90_put_var(ncid,Dopp90ID,real(rmout%Dopp90(:,:,:,iht)))
status =  nf90_put_var(ncid,SWh90ID,real(rmout%dDVh90(:,:,:,iht)))
status =  nf90_put_var(ncid,DVhID,real(rmout%DVh(:,:,:,iht)))
!
status =  nf90_put_var(ncid,KdpID,real(mout%Kdp(:,:,:,iht)))
status =  nf90_put_var(ncid,AdpID,real(mout%Adp(:,:,:,iht)))
status =  nf90_put_var(ncid,AhID,real(rmout%Ah(:,:,:,iht)))
status =  nf90_put_var(ncid,AvID,real(mout%Av(:,:,:,iht)))
status =  nf90_put_var(ncid,diff_back_phaseID,real(mout%diff_back_phase(:,:,:,iht)))
!
!---------------------------------------
! MPL
IF (conf%mplID>0) THEN
IF ((iht==1) .or. (iht==3).or.(iht==7)) THEN

if (iht==1) ihc=1
if (iht==3) ihc=2
if (iht==7) ihc=3
allocate(work(one)); work(1) =real(mpl%wavel)
status =  nf90_put_var(ncid,mpl_wavelID,work)
deallocate(work)
status =  nf90_put_var(ncid,mpl_back_obsID,real(mpl%back_obs(:,:,:,ihc)))
status =  nf90_put_var(ncid,mpl_back_trueID,real(mpl%back_true(:,:,:,ihc)))
status =  nf90_put_var(ncid,mpl_extID,real(mpl%ext(:,:,:,ihc)))
status =  nf90_put_var(ncid,mpl_lidar_ratioID,real(mpl%lidar_ratio(:,:,:,ihc)))

ENDIF
ENDIF
!---------------------------------------
!------------------------------------------
! Doppler spectra by oue in May 2018
if(conf%spectraID==1) then
status =  nf90_put_var(ncid,velbinsID,    real(spectra%vel_bins))
status =  nf90_put_var(ncid,spectra_zhhID,real(spectra%zhh_spectra(:,:,:,:,iht)))
status =  nf90_put_var(ncid,spectra_zvhID,real(spectra%zvh_spectra(:,:,:,:,iht)))
status =  nf90_put_var(ncid,spectra_zvvID,real(spectra%zvv_spectra(:,:,:,:,iht)))
end if
!---------------------------
!
!---------------------------------------
status =  nf90_put_var(ncid,elevID,real(rmout%elev))
status =  nf90_put_var(ncid,azimID,real(rmout%azim))
status =  nf90_put_var(ncid,rangeID,real(rmout%range))
status =  nf90_put_var(ncid,zID,real(env%z))
status =  nf90_put_var(ncid,qhydroID,real(qq))! kg/m3
!
allocate(work(one),iwork(one))
work(1)=real(conf%freq)
status =  nf90_put_var(ncid,rad_freqID,work)
work(1)=real(conf%Theta1)
status =  nf90_put_var(ncid,rad_beamwidthID,work)
work(1)=real(conf%dr)
status =  nf90_put_var(ncid,rad_range_resolutionID,work)
iwork(1)=conf%ixc
status =  nf90_put_var(ncid,ixcID,iwork)
iwork(1)=conf%iyc
status =  nf90_put_var(ncid,iycID,iwork)
work(1)=real(conf%zc)
status =  nf90_put_var(ncid,zcID,work)
deallocate(work,iwork)
!
!
status =  nf90_put_var(ncid,xID,real(env%x))
status =  nf90_put_var(ncid,yID,real(env%y))
!
ENDIF
!------------------------------------------------------------

999 If (status /= 0) then
write(*,*) 'Error in WriteOutNetcdf'
stop
endif
!
status = nf90_close(ncid)
return
end subroutine WriteOutNetcdf_mp
!
!
subroutine WriteOutNetcdf_rmp(iht,OutFileName,rmout,qq,env,conf,lout,mpl,aero,arscl,mwr,spectra,airborne,status) !arscl&mwr are added by oue 2017.03.23, added spectra by oue in May 2018
Use netcdf
Use crsim_mod
Use postprocess_mod !added by oue 2017.03.23
Implicit None

integer, intent(in)                    :: iht
character(len=*),Intent(in)            :: OutFileName
type(rmout_var),Intent(In)             :: rmout
Real*8,Intent(In)                      :: qq(rmout%nx,rmout%ny,rmout%nz) ! water content [kg/m^3]
type(env_var), Intent(in)              :: env
type(conf_var), Intent(in)             :: conf
type(lout_var),Intent(In)              :: lout
type(mpl_var),Intent(In)               :: mpl
type(aero_var),Intent(In)              :: aero
type(arscl_var),Intent(in)             :: arscl !added by oue 2017.03.23
type(mwr_var),Intent(in)               :: mwr !added by oue 2017.03.25
type(spectra_var),Intent(in)           :: spectra !added by oue 2018.05.27
type(airborne_var),Intent(in)          :: airborne !added by oue Oct 2023
Integer,Intent(out)                    :: status
!
character(len=160)                     :: err_msg
!
!integer  :: nx,ny,nz,nht
integer  :: nxID,nyID,nzID,nhtID
integer  :: ZhhID
integer  :: DVhID
integer  :: DoppID, Dopp90ID
integer  :: SWhID,SWh90ID
integer  :: SWtID,SWtotID
integer  :: SWsID,SWvID
integer  :: AhID
integer  :: elevID,azimID,rangeID,ZminID,AvapID,AtotID
integer  :: zID, tempID,uID,vID,wID,qhydroID,rho_dID
integer  :: xID,yID
integer  :: xlatID,xlongID
integer  :: ceilo_back_trueID,ceilo_extID,ceilo_back_obsID,lidar_ratioID
integer  :: mpl_wavelID,mpl_back_trueID,mpl_extID,mpl_back_obsID,mpl_lidar_ratioID,mpl_rayleigh_backID
integer  :: aero_back_trueID,aero_extID,aero_back_obsID,aero_lidar_ratioID
integer  :: ceilo_first_cloud_baseID
integer  :: oneID, rad_freqID,rad_beamwidthID,rad_range_resolutionID
integer  :: ixcID,iycID,zcID
!
!---------------------------
! For ARSCL added by oue 2017.03.23
integer  :: nlayersID
integer  :: arscl_cloud_maskID,arscl_cloud_detection_flagID,arscl_cloud_base_heightID,arscl_cloud_top_heightID
!--------------------------
! For MWR LWP added by oue 2017.03.25
integer  :: true_lwpID, mwr_lwpID, mwr_samplesID
!--------------------------
! For Doppler spectra added by oue in May 2018
Integer  :: nfftID
integer  :: spectra_zhhID, velbinsID
!--------------------------
! For airborne added by oue in Oct 2023    
integer  :: Zsfc_oceanID, Zsfc_land_coarseID, Zsfc_land_flatID, DV_airborneID
!--------------------------
!
integer  :: ncid
!
integer                                  :: ihc
Character(len=nf90_max_name)             :: att_string
Character(len=10)                        :: str1,str2
integer                                  :: j1,j2
integer,parameter                        :: one=1
real,Dimension(:),Allocatable            :: work
integer,Dimension(:),Allocatable         :: iwork

!
status =  nf90_create(path =OutFileName, cmode = NF90_NETCDF4, ncid = ncid)
If (status /= 0) then ; err_msg='Error in nf90_create' ;  Goto 999 ; endif

! define dimensions
status = nf90_def_dim(ncid, "nx",rmout%nx,nxID)
status = nf90_def_dim(ncid, "ny",rmout%ny,nyID)
status = nf90_def_dim(ncid, "nz",rmout%nz,nzID)
status = nf90_def_dim(ncid, "nht",rmout%nht,nhtID)
status = nf90_def_dim(ncid, "one",one,oneID)
if (conf%spectraID==1) then ! For Doppler spectra by oue in May 2018
status = nf90_def_dim(ncid, "nfft",spectra%NFFT,nfftID)
endif
if ((conf%arsclID==1) .and. (iht==0)) then ! For ARSCL products by oue 2017.03.23
status = nf90_def_dim(ncid, "n_layers",arscl%n_layers,nlayersID)
endif
!
! define variables
!
status = nf90_def_var(ncid, "xlat",  nf90_float, (/nxID,nyID/),xlatID)
status = nf90_def_var(ncid, "xlong", nf90_float, (/nxID,nyID/),xlongID)
!
status = nf90_def_var(ncid, "Zhh", nf90_float, (/nxID,nyID,nzID/), ZhhID)

status = nf90_def_var(ncid, "DV", nf90_float, (/nxID,nyID,nzID/), DoppID)
status = nf90_def_var(ncid, "SWh", nf90_float, (/nxID,nyID,nzID/), SWhID)
status = nf90_def_var(ncid, "SWt", nf90_float, (/nxID,nyID,nzID/), SWtID)
status = nf90_def_var(ncid, "SWs", nf90_float, (/nxID,nyID,nzID/), SWsID)
status = nf90_def_var(ncid, "SWv", nf90_float, (/nxID,nyID,nzID/), SWvID)
status = nf90_def_var(ncid, "SWtot", nf90_float, (/nxID,nyID,nzID/), SWtotID)
!
status = nf90_def_var(ncid, "DV90", nf90_float, (/nxID,nyID,nzID/), Dopp90ID)
status = nf90_def_var(ncid, "SWh90", nf90_float, (/nxID,nyID,nzID/), SWh90ID)
!
status = nf90_def_var(ncid, "RWV", nf90_float, (/nxID,nyID,nzID/), DVhID)
!
status = nf90_def_var(ncid, "Ah", nf90_float, (/nxID,nyID,nzID/),AhID)
!
!---------------------------
! CEILO
IF ((conf%ceiloID==1).and.(iht==0)) THEN
status = nf90_def_var(ncid, "ceilo_back_obs",     nf90_float, (/nxID,nyID,nzID/),ceilo_back_obsID)
status = nf90_def_var(ncid, "ceilo_back_true",    nf90_float, (/nxID,nyID,nzID/),ceilo_back_trueID)
status = nf90_def_var(ncid, "ceilo_ext",          nf90_float, (/nxID,nyID,nzID/),ceilo_extID)
status = nf90_def_var(ncid, "ceilo_lidar_ratio"   ,     nf90_float, (/nxID,nyID,nzID/),lidar_ratioID)
status = nf90_def_var(ncid, "ceilo_first_cloud_base",nf90_float, (/nxID,nyID/),ceilo_first_cloud_baseID)
ENDIF
!---------------------------
!---------------------------
! MPL
IF (conf%mplID>0) THEN
!
IF (iht==0) THEN
status = nf90_def_var(ncid, "mpl_rayleigh_back",nf90_float,(/nxID,nyID,nzID/),mpl_rayleigh_backID)
ENDIF
!
IF ( (iht==0) .or. (iht==1) .or. (iht==3).or.(iht==7) ) THEN
status = nf90_def_var(ncid, "mpl_wavel", nf90_float,(/oneID/),mpl_wavelID)
status = nf90_def_var(ncid, "mpl_back_obs",  nf90_float, (/nxID,nyID,nzID/),mpl_back_obsID)
status = nf90_def_var(ncid, "mpl_back_true", nf90_float, (/nxID,nyID,nzID/),mpl_back_trueID)
status = nf90_def_var(ncid, "mpl_ext", nf90_float, (/nxID,nyID,nzID/),mpl_extID)
status = nf90_def_var(ncid, "mpl_lidar_ratio", nf90_float, (/nxID,nyID,nzID/),mpl_lidar_ratioID)
ENDIF
! AEROSOL
IF(conf%aeroID==1) THEN
IF (iht==0) THEN
status = nf90_def_var(ncid, "aero_back_obs",nf90_float,(/nxID,nyID,nzID/),aero_back_obsID)
status = nf90_def_var(ncid, "aero_back_true",nf90_float,(/nxID,nyID,nzID/),aero_back_trueID)
status = nf90_def_var(ncid, "aero_ext",nf90_float,(/nxID,nyID,nzID/),aero_extID)
status = nf90_def_var(ncid, "aero_lidar_ratio",nf90_float,(/nxID,nyID,nzID/),aero_lidar_ratioID)
ENDIF
ENDIF
!
ENDIF !if conf%mplID>0
!---------------------------
!---------------------------
! Doppler spectra by oue in May 2018
IF (conf%spectraID>0) THEN
status = nf90_def_var(ncid, "velocity_bins",  nf90_float,(/nfftID/),velbinsID)
status = nf90_def_var(ncid, "spectra_zhh",    nf90_float,(/nxID,nyID,nzID,nfftID/),spectra_zhhID)
ENDIF
!---------------------------
! Airborne by oue in Oct 2023
IF ((conf%airborne==1) .and. (iht==0)) THEN
   status = nf90_def_var(ncid, "Zsfc_Ocean",        nf90_float,(/nxID,nyID/),Zsfc_oceanID)
   status = nf90_def_var(ncid, "Zsfc_Land_Coarse",  nf90_float,(/nxID,nyID/),Zsfc_land_coarseID)
   status = nf90_def_var(ncid, "Zsfc_Land_Flat",    nf90_float,(/nxID,nyID/),Zsfc_land_flatID)
   status = nf90_def_var(ncid, "DV_airborne",       nf90_float,(/nxID,nyID,nzID/),DV_airborneID)
ENDIF
!---------------------------
! ARSCL added by oue 2017.03.23
if ((conf%arsclID==1) .and. (conf%mplID>0) .and. (conf%ceiloID==1)) then 
if (iht==0) then
status = nf90_def_var(ncid, "arscl_cloud_mask", nf90_short,(/nxID,nyID,nzID/),arscl_cloud_maskID)
status = nf90_def_var(ncid, "arscl_cloud_source_flag", nf90_short,(/nxID,nyID,nzID/),arscl_cloud_detection_flagID)
status = nf90_def_var(ncid, "arscl_cloud_layer_base_height", nf90_float,(/nxID,nyID,nlayersID/),arscl_cloud_base_heightID)
status = nf90_def_var(ncid, "arscl_cloud_layer_top_height", nf90_float,(/nxID,nyID,nlayersID/),arscl_cloud_top_heightID)
endif
endif
!---------------------------
! MWR LWP added by oue 2017.03.25
if ((conf%mwrID==1) .and. (iht==0)) then
status = nf90_def_var(ncid, "model_lwp", nf90_float,(/nxID,nyID/),true_lwpID)
status = nf90_def_var(ncid, "mwr_lwp", nf90_float,(/nxID,nyID/),mwr_lwpID)
status = nf90_def_var(ncid, "number_of_gridpoints_mwrlwp", nf90_int,(/nxID,nyID,nzID/),mwr_samplesID)
endif
!---------------------------
!---------------------------
!
!---------------------------
if (iht==0) Then
status = nf90_def_var(ncid, "Zmin",  nf90_float, (/nxID,nyID,nzID/),ZminID)
status = nf90_def_var(ncid, "temp",  nf90_float, (/nxID,nyID,nzID/),tempID)
status = nf90_def_var(ncid, "rho_d",  nf90_float, (/nxID,nyID,nzID/),rho_dID)
status = nf90_def_var(ncid, "u",     nf90_float, (/nxID,nyID,nzID/),uID)
status = nf90_def_var(ncid, "v",     nf90_float, (/nxID,nyID,nzID/),vID)
status = nf90_def_var(ncid, "w",     nf90_float, (/nxID,nyID,nzID/),wID)
status = nf90_def_var(ncid, "Avap",  nf90_float, (/nxID,nyID,nzID/),AvapID)
If(dabs(conf%elev)==90) Then
status = nf90_def_var(ncid, "Atot",  nf90_float, (/nxID,nyID,nzID/),AtotID)
endif
endif 
!---------------------------
!
status = nf90_def_var(ncid, "rad_freq",  nf90_float, (/oneID/),rad_freqID)
status = nf90_def_var(ncid, "rad_beamwidth", nf90_float, (/oneID/),rad_beamwidthID)
status = nf90_def_var(ncid, "rad_range_resolution",nf90_float,(/oneID/),rad_range_resolutionID)
status = nf90_def_var(ncid, "rad_ixc",nf90_int,(/oneID/),ixcID)
status = nf90_def_var(ncid, "rad_iyc",nf90_int,(/oneID/),iycID)
status = nf90_def_var(ncid, "rad_zc",nf90_float,(/oneID/),zcID)
!
status = nf90_def_var(ncid, "elev",  nf90_float, (/nxID,nyID,nzID/),elevID)
status = nf90_def_var(ncid, "azim",  nf90_float, (/nxID,nyID,nzID/),azimID)
status = nf90_def_var(ncid, "range",  nf90_float, (/nxID,nyID,nzID/),rangeID)
status = nf90_def_var(ncid, "height",nf90_float, (/nxID,nyID,nzID/),zID)
status = nf90_def_var(ncid, "wcont", nf90_float, (/nxID,nyID,nzID/),qhydroID)
!
status = nf90_def_var(ncid, "x_scene",  nf90_float, (/nxID/),xID)
status = nf90_def_var(ncid, "y_scene",  nf90_float, (/nyID/),yID)
!
! Attributes
!
    status = nf90_put_att(ncid,xlatID,'description','Latitude (North is positive)')
    status = nf90_put_att(ncid,xlatID,'units','deg')
    status = nf90_put_att(ncid,xlatID,'_FillValue',-999.e0)
!
    status = nf90_put_att(ncid,xlongID,'description','Longitude (West is positive)')
    status = nf90_put_att(ncid,xlongID,'units','deg')
    status = nf90_put_att(ncid,xlongID,'_FillValue',-999.e0)
!
    status = nf90_put_att(ncid,ZhhID,'description','Reflectivity at hh polarization')
    status = nf90_put_att(ncid,ZhhID,'units','dBZ')
    status = nf90_put_att(ncid,ZhhID,'_FillValue',-999.e0)
!
    status = nf90_put_att(ncid,DoppID,'description','Radial Doppler velocity, positive upward')
    status = nf90_put_att(ncid,DoppID,'units','m/s')
    status = nf90_put_att(ncid,DoppID,'_FillValue',-999.e0)
!
    status = nf90_put_att(ncid,SWhID,'description','Spectrum width due to hydrometeors')
    status = nf90_put_att(ncid,SWhID,'units','m/s')
    status = nf90_put_att(ncid,SWhID,'_FillValue',-999.e0)
!
    status = nf90_put_att(ncid,SWtID,'description','Spectrum width due to turbulence')
    status = nf90_put_att(ncid,SWtID,'units','m/s')
    status = nf90_put_att(ncid,SWtID,'_FillValue',-999.e0)
!
    status = nf90_put_att(ncid,SWsID,'description','Spectrum width due to wind shear')
    status = nf90_put_att(ncid,SWsID,'units','m/s')
    status = nf90_put_att(ncid,SWsID,'_FillValue',-999.e0)
!
    status = nf90_put_att(ncid,SWvID,'description','Spectrum width due to cross wind')
    status = nf90_put_att(ncid,SWvID,'units','m/s')
    status = nf90_put_att(ncid,SWvID,'_FillValue',-999.e0)
!
    status = nf90_put_att(ncid,SWtotID,'description','Total spectrum width')
    status = nf90_put_att(ncid,SWtotID,'units','m/s')
    status = nf90_put_att(ncid,SWtotID,'_FillValue',-999.e0)
!   
    status = nf90_put_att(ncid,Dopp90ID,'description','Vertical Doppler velocity (el=90), positive upward')
    status = nf90_put_att(ncid,Dopp90ID,'units','m/s')
    status = nf90_put_att(ncid,Dopp90ID,'_FillValue',-999.e0)
!
    status = nf90_put_att(ncid,SWh90ID,'description','Spectrum width due to hydrometeors, elevation=90')
    status = nf90_put_att(ncid,SWh90ID,'units','m/s')
    status = nf90_put_att(ncid,SWh90ID,'_FillValue',-999.e0)
!
    status = nf90_put_att(ncid,DVhID,'description','Reflectivity weighted velocity (el=90deg, w=0m/s)')
    status = nf90_put_att(ncid,DVhID,'units','m/s')
    status = nf90_put_att(ncid,DVhID,'_FillValue',-999.e0)
!   
    status = nf90_put_att(ncid,AhID,'description','Specific horizontal attenuation')
    status = nf90_put_att(ncid,AhID,'units','dB/km')
    status = nf90_put_att(ncid,AhID,'_FillValue',-999.e0)
!
!------------------------------------------
! CEILO
IF ( (conf%ceiloID==1).and. (iht==0)) THEN
!
    status = nf90_put_att(ncid,ceilo_back_obsID,'description','Ceilo observed backscatter')
    status = nf90_put_att(ncid,ceilo_back_obsID,'units','[m sr]^-1')
    status = nf90_put_att(ncid,ceilo_back_obsID,'_FillValue',-999.e0)
!
    status = nf90_put_att(ncid,ceilo_back_trueID,'description','Ceilo true backscatter')
    status = nf90_put_att(ncid,ceilo_back_trueID,'units','[m sr]^-1')
    status = nf90_put_att(ncid,ceilo_back_trueID,'_FillValue',-999.e0)
!
    status = nf90_put_att(ncid,ceilo_extID,'description','Ceilo extinction coefficient')
    status = nf90_put_att(ncid,ceilo_extID,'units','[m]^-1')
    status = nf90_put_att(ncid,ceilo_extID,'_FillValue',-999.e0)
!
    status = nf90_put_att(ncid,lidar_ratioID,'description','Ceilo lidar ratio')
    status = nf90_put_att(ncid,lidar_ratioID,'units','sr')
    status = nf90_put_att(ncid,lidar_ratioID,'_FillValue',-999.e0)
!
    status = nf90_put_att(ncid,ceilo_first_cloud_baseID,'description','Height of the first cloud base')
    status = nf90_put_att(ncid,ceilo_first_cloud_baseID,'units','m')
    status = nf90_put_att(ncid,ceilo_first_cloud_baseID,'_FillValue',-999.e0)
!
ENDIF
!------------------------------------------
!------------------------------------------ 
! MPL
IF (conf%mplID>0) THEN
IF ( (iht==0) .or. (iht==1) .or. (iht==3).or.(iht==7) ) THEN
!
    status = nf90_put_att(ncid,mpl_wavelID,'description','MPL wavelength')
    status = nf90_put_att(ncid,mpl_wavelID,'units','nm')
!
    if (iht==0) status = nf90_put_att(ncid,mpl_back_obsID,'description','MPL observed backscatter (cloud+ice+aerosol+atm.moleculs)')
    if ((iht==1) .or.(iht==3).or.(iht==7))  status = nf90_put_att(ncid,mpl_back_obsID, &
                               'description','MPL observed backscatter (hydrometeor only)')
    status = nf90_put_att(ncid,mpl_back_obsID,'units','[m sr]^-1')
    status = nf90_put_att(ncid,mpl_back_obsID,'_FillValue',-999.e0)
!   
    if (iht==0) status = nf90_put_att(ncid,mpl_back_trueID,'description','MPL true backscatter (cloud+ice+aerosol+atm.moleculs)')
    if ((iht==1) .or.(iht==3).or.(iht==7))  status =nf90_put_att(ncid,mpl_back_trueID,&
         'description','MPL true backscatter (hydrometeor only)')
    status = nf90_put_att(ncid,mpl_back_trueID,'units','[m sr]^-1')
    status = nf90_put_att(ncid,mpl_back_trueID,'_FillValue',-999.e0)
!  
    if (iht==0) status = nf90_put_att(ncid,mpl_rayleigh_backID,'description','MPL molecular backscatter')
    status = nf90_put_att(ncid,mpl_rayleigh_backID,'units','[m sr]^-1')
    status = nf90_put_att(ncid,mpl_rayleigh_backID,'_FillValue',-999.e0)
!  
    if (iht==0) status =  nf90_put_att(ncid,mpl_extID,'description','MPL extinction coefficient (cloud+ice+aerosol+atm.moleculs)')
    if ((iht==1) .or.(iht==3).or.(iht==7)) status = nf90_put_att(ncid,mpl_extID,&
         'description','MPL extinction coefficient (hydrometeor only)')
    status = nf90_put_att(ncid,mpl_extID,'units','[m]^-1')
    status = nf90_put_att(ncid,mpl_extID,'_FillValue',-999.e0)
!   
    if (iht==0) status = nf90_put_att(ncid,mpl_lidar_ratioID,'description','MPL Lidar ratio (cloud+ice+aerosol+atm.moleculs)')
    if ((iht==1) .or.(iht==3).or.(iht==7)) status = nf90_put_att(ncid,mpl_lidar_ratioID,&
         'description','MPL Lidar ratio (hydrometeor only)')
    status = nf90_put_att(ncid,mpl_lidar_ratioID,'units','sr')
    status = nf90_put_att(ncid,mpl_lidar_ratioID,'_FillValue',-999.e0)
!
ENDIF
IF (conf%aeroID==1) THEN
IF (iht==0) THEN
!
    status = nf90_put_att(ncid,aero_back_obsID,'description','MPL aerosol observed backscatter')
    status = nf90_put_att(ncid,aero_back_obsID,'units','[m sr]^-1')
    status = nf90_put_att(ncid,aero_back_obsID,'_FillValue',-999.e0)
!   
    status = nf90_put_att(ncid,aero_back_trueID,'description','MPL aerosol true backscatter')
    status = nf90_put_att(ncid,aero_back_trueID,'units','[m sr]^-1')
    status = nf90_put_att(ncid,aero_back_trueID,'_FillValue',-999.e0)
!   
    status = nf90_put_att(ncid,aero_extID,'description','MPL aerosol extinction coefficient')
    status = nf90_put_att(ncid,aero_extID,'units','[m]^-1')
    status = nf90_put_att(ncid,aero_extID,'_FillValue',-999.e0)
!   
    status = nf90_put_att(ncid,aero_lidar_ratioID,'description','MPL aerosol Lidar ratio')
    status = nf90_put_att(ncid,aero_lidar_ratioID,'units','sr')
    status = nf90_put_att(ncid,aero_lidar_ratioID,'_FillValue',-999.e0)
!
ENDIF !  iht==0
ENDIF !  conf%aeroID==1
!
!
ENDIF !conf%mplID>0
!------------------------------------------
!------------------------------------------
! Doppler spectra by oue in May 2018
IF (conf%spectraID>0) THEN
    status = nf90_put_att(ncid,velbinsID,'description','Doppler spectrum velocity')
    status = nf90_put_att(ncid,velbinsID,'units','[m s]^-1')
    status = nf90_put_att(ncid,velbinsID,'_FillValue',-999.e0)

    status = nf90_put_att(ncid,spectra_zhhID,'description','copolar Doppler spectra for horizontal polarization')
    status = nf90_put_att(ncid,spectra_zhhID,'units','dB')
    status = nf90_put_att(ncid,spectra_zhhID,'_FillValue',-999.e0)
ENDIF
!---------------------------
! Airborne by oue in Oct 2023
IF ((conf%airborne==1) .and. (iht==0)) THEN
    status = nf90_put_att(ncid,Zsfc_oceanID,'description','Surface reflectivity over ocean')
    status = nf90_put_att(ncid,Zsfc_oceanID,'units','dBZ')
    status = nf90_put_att(ncid,Zsfc_oceanID,'_FillValue',-999.e0)

    status = nf90_put_att(ncid,Zsfc_land_coarseID,'description','Surface reflectivity over coarse land (i.e. wet snow)')
    status = nf90_put_att(ncid,Zsfc_land_coarseID,'units','dBZ')
    status = nf90_put_att(ncid,Zsfc_land_coarseID,'_FillValue',-999.e0)

    status = nf90_put_att(ncid,Zsfc_land_flatID,'description','Surface reflectivity over flat land (i.e. flat road)')
    status = nf90_put_att(ncid,Zsfc_land_flatID,'units','dBZ')
    status = nf90_put_att(ncid,Zsfc_land_flatID,'_FillValue',-999.e0)

    status = nf90_put_att(ncid,DV_airborneID,'description','Radial airmotion related to airborne moving')
    status = nf90_put_att(ncid,DV_airborneID,'units','[m s]^-1')
    status = nf90_put_att(ncid,DV_airborneID,'_FillValue',-999.e0)
ENDIF
!---------------------------
!---------------------------
! ARSCL atteributes added by oue 2017.03.23
if ((conf%arsclID==1) .and. (conf%mplID>0) .and. (conf%ceiloID==1)) then 
if (iht==0) then
    status = nf90_put_att(ncid,arscl_cloud_maskID,'description','Cloud mask from radar, mpl, ceilometer observations')
    status = nf90_put_att(ncid,arscl_cloud_maskID,'units','unitless')
    status = nf90_put_att(ncid,arscl_cloud_maskID,'flag_0_description','clear')
    status = nf90_put_att(ncid,arscl_cloud_maskID,'flag_1_description','cloudy')
!
    status = nf90_put_att(ncid,arscl_cloud_detection_flagID,'description','Instrument source flag for cloud detections')
    status = nf90_put_att(ncid,arscl_cloud_detection_flagID,'units','unitless')
    status = nf90_put_att(ncid,arscl_cloud_detection_flagID,'_FillValue',int2(-999))
    status = nf90_put_att(ncid,arscl_cloud_detection_flagID,'flag_values','1s, 2s, 3s, 4s')
    status = nf90_put_att(ncid,arscl_cloud_detection_flagID,'flag_1_description','Clear according to radar and MPL')
    status = nf90_put_att(ncid,arscl_cloud_detection_flagID,'flag_2_description','Cloud detected by radar and MPL')
    status = nf90_put_att(ncid,arscl_cloud_detection_flagID,'flag_3_description','Cloud detected by radar only')
    status = nf90_put_att(ncid,arscl_cloud_detection_flagID,'flag_4_description','Cloud detected by MPL only')
!
    status = nf90_put_att(ncid,arscl_cloud_base_heightID,'description','Base height of cloudy layers for up to 10 layers')
    status =  nf90_put_att(ncid,arscl_cloud_base_heightID,'comments','based on combined radar, ceilometer, and mpl observations')
    status = nf90_put_att(ncid,arscl_cloud_base_heightID,'units','m')
    status = nf90_put_att(ncid,arscl_cloud_base_heightID,'_FillValue',-999.e0)
!
    status = nf90_put_att(ncid,arscl_cloud_top_heightID,'description','Top height of cloudy layers for up to 10 layers')
    status = nf90_put_att(ncid,arscl_cloud_top_heightID,'comments','based on combined radar, ceilometer, and mpl observations')
    status = nf90_put_att(ncid,arscl_cloud_top_heightID,'units','m')
    status = nf90_put_att(ncid,arscl_cloud_top_heightID,'_FillValue',-999.e0)
endif
endif !conf%arsclID==1
!---------------------------
! MWR LWP attributes added by oue 2017.03.25
if ((conf%mwrID==1) .and. (iht==0)) then
    status = nf90_put_att(ncid,true_lwpID,'description','Model liquid water path (cloud+rain)')
    status = nf90_put_att(ncid,true_lwpID,'units','kg m^-2')

    status = nf90_put_att(ncid,mwr_lwpID,'description','Microwave radiometer liquid water path taking account of field of view')
    status = nf90_put_att(ncid,mwr_lwpID,'units','kg m^-2')
    status = nf90_put_att(ncid,mwr_lwpID,'field_of_view_deg',conf%mwr_view)
    status = nf90_put_att(ncid,mwr_lwpID,'field_of_mwr',conf%mwr_alt)

    status = nf90_put_att(ncid,mwr_samplesID,'description','Number of gridpoints within MWR field of view')
    status = nf90_put_att(ncid,mwr_samplesID,'units','unitless')
endif !conf%mwrID==1
!---------------------------
!---------------------------
!
if (iht==0) then
!
    status = nf90_put_att(ncid,ZminID,'description','Radar sensitivity limitation with range')
    status = nf90_put_att(ncid,ZminID,'units','dBZ')
    status = nf90_put_att(ncid,ZminID,'_FillValue',-999.e0)
!
    status = nf90_put_att(ncid,tempID,'description','Temperature')
    status = nf90_put_att(ncid,tempID,'units','C')
!
    status = nf90_put_att(ncid,rho_dID,'description','Dry air density')
    status = nf90_put_att(ncid,rho_dID,'units','kg m^-3')
!
    status = nf90_put_att(ncid,uID,'description','U horizontal wind component')
    status = nf90_put_att(ncid,uID,'units','m/s')
!
    status = nf90_put_att(ncid,vID,'description','V horizontal wind component')
    status = nf90_put_att(ncid,vID,'units','m/s')
!
    status = nf90_put_att(ncid,wID,'description','Vertical air velocity')
    status = nf90_put_att(ncid,wID,'units','m/s')
!
    status = nf90_put_att(ncid,AvapID,'description','Specific water vapor attenuation')
    status = nf90_put_att(ncid,AvapID,'units','dB/km')
    status = nf90_put_att(ncid,AvapID,'_FillValue',-999.e0)
!
If(dabs(conf%elev)==90) Then
    status = nf90_put_att(ncid,AtotID,'description','Two-way total attenuation (Ah+Avap)')
    status = nf90_put_att(ncid,AtotID,'units','dB')
    status = nf90_put_att(ncid,AtotID,'_FillValue',-999.e0)
endif
 endif
!---------------------------
!
    status = nf90_put_att(ncid,elevID,'description','Elevation')
    status = nf90_put_att(ncid,elevID,'units','degrees')
    status = nf90_put_att(ncid,elevID,'_FillValue',-999.e0)
!
    status = nf90_put_att(ncid,azimID,'description','Azimuth, E=0, W=180, N=90, S=270 degs')
    status = nf90_put_att(ncid,azimID,'units','degrees')
    status = nf90_put_att(ncid,azimID,'_FillValue',-999.e0)
!
    status = nf90_put_att(ncid,rangeID,'description','Radar range')
    status = nf90_put_att(ncid,rangeID,'units','m')
    status = nf90_put_att(ncid,rangeID,'_FillValue',-999.e0)
!
    status = nf90_put_att(ncid,zID,'description','Height')
    status = nf90_put_att(ncid,zID,'units','m')
!
    status = nf90_put_att(ncid,qhydroID,'description','Water content')
    status = nf90_put_att(ncid,qhydroID,'units','kg/m^3')
    status = nf90_put_att(ncid,qhydroID,'_FillValue',-999.e0)
!
    status = nf90_put_att(ncid,rad_freqID,'description','radar frequency')
    status = nf90_put_att(ncid,rad_freqID,'long_name','configuration parameter')
    status = nf90_put_att(ncid,rad_freqID,'units','GHz')
!
    status = nf90_put_att(ncid,rad_beamwidthID,'description','radar antenna beamwidth')
    status = nf90_put_att(ncid,rad_beamwidthID,'long_name','configuration parameter')
    status = nf90_put_att(ncid,rad_beamwidthID,'units','deg')
!
    status = nf90_put_att(ncid,rad_range_resolutionID,'description','radar range resolution')
    status = nf90_put_att(ncid,rad_range_resolutionID,'long_name','configuration parameter')
    status = nf90_put_att(ncid,rad_range_resolutionID,'units','m')
!
    status = nf90_put_att(ncid,ixcID,'description','index of position in x direction of radar origin')
    status = nf90_put_att(ncid,ixcID,'long_name','configuration parameter')
    status = nf90_put_att(ncid,ixcID,'units','-')
!
    status = nf90_put_att(ncid,iycID,'description','index of position in y direction of radar origin')
    status = nf90_put_att(ncid,iycID,'long_name','configuration parameter')
    status = nf90_put_att(ncid,iycID,'units','-')
!
    status = nf90_put_att(ncid,zcID,'description','radar height')
    status = nf90_put_att(ncid,zcID,'long_name','configuration parameter')
    status = nf90_put_att(ncid,zcID,'units','m')
!
    status = nf90_put_att(ncid,xID,'description','Scene extent in E-W direction')
    status = nf90_put_att(ncid,xID,'units','m')
!
    status = nf90_put_att(ncid,yID,'description','Scene extent in N-S direction')
    status = nf90_put_att(ncid,yID,'units','m')
!
!---------------------------------------------------------------------------------------
! GLOBAL ATTRIBUTES

   att_string="forward radar simulator output"
   status = nf90_put_att(ncid, NF90_GLOBAL,'description',att_string)
   !
   status = nf90_put_att(ncid, NF90_GLOBAL,'model_version','crsim_v4.0')
   !
   att_string=Trim(conf%WRFInputFile)
   status = nf90_put_att(ncid, NF90_GLOBAL,'WRF_input_file',att_string)
   !
   write( str1, '(i5)' ) conf%MP_PHYSICS
   att_string=Trim(str1)
   status = nf90_put_att(ncid, NF90_GLOBAL,'MP_PHYSICS',att_string)
   !
   if (Min(conf%ix_start,conf%ix_end)<0) then
   j1=1 ; j2=rmout%nx
   else
   j1=conf%ix_start ; j2=conf%ix_end
   endif
   write( str1, '(i10)' ) j1
   write( str2, '(i10)' ) j2
   att_string=Trim(Adjustl(str1))//" - "//Trim(Adjustl(str2))
   status = nf90_put_att(ncid,NF90_GLOBAL,'x_indices_of_WRF_extracted_scene',att_string)
   !
   if (Min(conf%iy_start,conf%iy_end)<0) then
   j1=1 ; j2=rmout%ny
   else
   j1=conf%iy_start ; j2=conf%iy_end
   endif
   write( str1, '(i10)' ) j1
   write( str2, '(i10)' ) j2
   att_string=Trim(Adjustl(str1))//" - "//Trim(Adjustl(str2))
   status = nf90_put_att(ncid,NF90_GLOBAL,'y_indices_of_WRF_extracted_scene',att_string)
   !
   if (Min(conf%iz_start,conf%iz_end)<0) then
   j1=1 ; j2=rmout%nz
   else
   j1=conf%iz_start ; j2=conf%iz_end
   endif
   write( str1, '(i10)' ) j1
   write( str2, '(i10)' ) j2
   att_string=Trim(Adjustl(str1))//" - "//Trim(Adjustl(str2))
   status = nf90_put_att(ncid,NF90_GLOBAL,'z_indices_of_WRF_extracted_scene',att_string)
   !
   write( str1, '(i10)' ) conf%it
   att_string=Trim(Adjustl(str1))
   status = nf90_put_att(ncid,NF90_GLOBAL,'scene_extracted_at_time_step',att_string)
   !
  !
   write( str1, '(f10.1)' ) conf%freq
   att_string=Trim(Adjustl(str1))//" GHz"
   status = nf90_put_att(ncid, NF90_GLOBAL,'radar_frequency',att_string)
   !
   write( str1, '(i10)' ) conf%ixc
   write( str2, '(i10)' ) conf%iyc
   att_string=Trim(Adjustl(str1))//" , "//Trim(Adjustl(str2))
   status =nf90_put_att(ncid,NF90_GLOBAL,'x_and_y_indices_of_radar_position',att_string)
   !
   write( str1, '(f10.1)' ) conf%zc
   att_string=Trim(Adjustl(str1))//" m"
   status = nf90_put_att(ncid,NF90_GLOBAL,'height_of_radar',att_string)
   !
   if (conf%elev<0.d0) then
   att_string="elevation of each scene pixel is relative to the radar origin"
   else
   write( str1, '(f5.1)' ) conf%elev
   att_string="elevation is set to "//Trim(Adjustl(str1))//" degrees"
   endif
   status = nf90_put_att(ncid,NF90_GLOBAL,'scanning_mode',att_string)
   !
   !------------------------------------------
   ! Doppler spectra by oue in May 2018
   if(conf%spectraID==1) then
   status = nf90_put_att(ncid, NF90_GLOBAL,'PRF_for_Doppler_spectra',spectra%PRF)
   status = nf90_put_att(ncid, NF90_GLOBAL,'TimeSampling_for_Doppler_spectra',spectra%TimeSampling)
   status = nf90_put_att(ncid, NF90_GLOBAL,'Noise1km_for_Doppler_spectra',spectra%NOISE_1km)
   status = nf90_put_att(ncid, NF90_GLOBAL,'Nyquist_velocity_for_Doppler_spectra',spectra%VNyquist)
   status = nf90_put_att(ncid, NF90_GLOBAL,'Nave_for_Doppler_spectra',spectra%Nave)
   end if
   !------------------------------------------
   ! Airborne by oue in Oct 2023
   if((conf%airborne==1) .and. (iht==0)) then
   status = nf90_put_att(ncid, NF90_GLOBAL,'Pulse_length_m',conf%pulse_len)
   status = nf90_put_att(ncid, NF90_GLOBAL,'Airborne_moving_speed_mps',conf%airborne_spd)
   status = nf90_put_att(ncid, NF90_GLOBAL,'Airborne_moving_azimuth_direction',conf%airborne_azdeg)
   status = nf90_put_att(ncid, NF90_GLOBAL,'Airborne_moving_elevation_direction',conf%airborne_eldeg)
   end if
   !------------------------------------------
   !
   status = nf90_put_att(ncid, NF90_GLOBAL,'created_by','Aleksandra Tatarevic')
   status = nf90_put_att(ncid,NF90_GLOBAL,'institute','http://wwww.clouds.mcgill.ca')
   status = nf90_put_att(ncid,NF90_GLOBAL,'websites',&
        'https://www.bnl.gov/CMAS/cr-sim.php; https://you.stonybrook.edu/radar/research/radar-simulators/')

!---------------------------------------------------------------------------------------

status = nf90_enddef(ncid)
If (status /= 0) then ; err_msg='Error in nf90_enddef'; Goto 999 ; endif

!------------------------------------------------------------
IF (iht==0) THEN
! write variables
status =  nf90_put_var(ncid,xlatID,real(env%xlat))
status =  nf90_put_var(ncid,xlongID,real(env%xlong))
!
status =  nf90_put_var(ncid,ZhhID,real(rmout%Zhh_tot))
status =  nf90_put_var(ncid,DoppID,real(rmout%Dopp_tot))
status =  nf90_put_var(ncid,SWhID,real(rmout%dDVh_tot))
status =  nf90_put_var(ncid,SWtID,real(rmout%sw_t))
status =  nf90_put_var(ncid,SWsID,real(rmout%sw_s))
status =  nf90_put_var(ncid,SWvID,real(rmout%sw_v))
status =  nf90_put_var(ncid,SWtotID,real(rmout%SWt_tot))
status =  nf90_put_var(ncid,Dopp90ID,real(rmout%Dopp90_tot))
status =  nf90_put_var(ncid,SWh90ID,real(rmout%dDVh90_tot))
status =  nf90_put_var(ncid,DVhID,real(rmout%DVh_tot))
status =  nf90_put_var(ncid,AhID,real(rmout%Ah_tot))
!
!---------------------------------------
! CEILO modified by oue for RAMS
IF (conf%ceiloID==1) THEN
status =  nf90_put_var(ncid,ceilo_back_obsID,real(lout%ceilo_back_obs_tot))
status =  nf90_put_var(ncid,ceilo_back_trueID,real(lout%ceilo_back_true_tot))
status =  nf90_put_var(ncid,ceilo_extID,real(lout%ceilo_ext_tot))
status =  nf90_put_var(ncid,lidar_ratioID,real(lout%lidar_ratio_tot))
status =  nf90_put_var(ncid,ceilo_first_cloud_baseID,real(lout%ceilo_first_cloud_base))
ENDIF
!---------------------------------------
! MPL
IF (conf%mplID>0) THEN 
allocate(work(one)); work(1)=real(mpl%wavel)
status =  nf90_put_var(ncid,mpl_wavelID,work)
deallocate(work)
status =  nf90_put_var(ncid,mpl_back_obsID,real(mpl%back_obs_tot))
status =  nf90_put_var(ncid,mpl_back_trueID,real(mpl%back_true_tot))
status =  nf90_put_var(ncid,mpl_rayleigh_backID,real(mpl%rayleigh_back))
status =  nf90_put_var(ncid,mpl_extID,real(mpl%ext_tot))
status =  nf90_put_var(ncid,mpl_lidar_ratioID,real(mpl%lidar_ratio_tot))
!! aerosol
IF (conf%aeroID==1) THEN
status =  nf90_put_var(ncid,aero_back_obsID,real(aero%back_obs))
status =  nf90_put_var(ncid,aero_back_trueID,real(aero%back_true))
status =  nf90_put_var(ncid,aero_extID,real(aero%ext))
status =  nf90_put_var(ncid,aero_lidar_ratioID,real(aero%lidar_ratio))
ENDIF ! if (conf%aeroID==1)
ENDIF ! if (conf%mplID>0)
!---------------------------------------
!------------------------------------------
! Doppler spectra by oue in May 2018
if(conf%spectraID==1) then
status =  nf90_put_var(ncid,velbinsID,    real(spectra%vel_bins))
status =  nf90_put_var(ncid,spectra_zhhID,real(spectra%zhh_spectra_tot))
end if
!---------------------------
! Airborne added by oue Oct 2023
if((conf%airborne==1) .and. (iht==0)) then
status =  nf90_put_var(ncid,Zsfc_oceanID,      real(airborne%Zsfc_ocean))
status =  nf90_put_var(ncid,Zsfc_land_coarseID,real(airborne%Zsfc_land_coarse))
status =  nf90_put_var(ncid,Zsfc_land_flatID,  real(airborne%Zsfc_land_flat))
status =  nf90_put_var(ncid,DV_airborneID,     real(airborne%Dopp_airborne))
end if
!---------------------------
! ARSCL  added by oue 2017.03.23
if ((conf%arsclID==1) .and. (conf%mplID>0) .and. (conf%ceiloID==1)) then 
status =  nf90_put_var(ncid,arscl_cloud_maskID,arscl%cloud_mask)
status =  nf90_put_var(ncid,arscl_cloud_detection_flagID,arscl%cloud_detection_flag)
status =  nf90_put_var(ncid,arscl_cloud_base_heightID,real(arscl%cloud_layer_base))
status =  nf90_put_var(ncid,arscl_cloud_top_heightID,real(arscl%cloud_layer_top))
endif ! if (conf%arsclID==1) 
!--------------------------- 
! MWR LWP added by oue 2017.03.25 
if (conf%mwrID==1) then
status =  nf90_put_var(ncid,true_lwpID,mwr%true_lwp)
status =  nf90_put_var(ncid,mwr_lwpID,mwr%mwr_lwp)
status =  nf90_put_var(ncid,mwr_samplesID,mwr%n_samples)
endif ! if (conf%mwrID==1)   
!---------------------------------------
!---------------------------------------
!
status =  nf90_put_var(ncid,ZminID, real(conf%ZMIN + 20.d0 *  dlog10(1.d-3 * rmout%range)))
status =  nf90_put_var(ncid,tempID,real(env%temp))
status =  nf90_put_var(ncid,rho_dID,real(env%rho_d))
status =  nf90_put_var(ncid,uID,real(env%u))
status =  nf90_put_var(ncid,vID,real(env%v))
status =  nf90_put_var(ncid,wID,real(env%w))
status =  nf90_put_var(ncid,AvapID,real(rmout%Avap))
If(dabs(conf%elev)==90) Then
status =  nf90_put_var(ncid,AtotID,real(rmout%Atot))
endif
!---------------------------------------
!
status =  nf90_put_var(ncid,elevID,real(rmout%elev))
status =  nf90_put_var(ncid,azimID,real(rmout%azim))
status =  nf90_put_var(ncid,rangeID,real(rmout%range))
status =  nf90_put_var(ncid,zID,real(env%z))
status =  nf90_put_var(ncid,qhydroID,real(qq))! kg/m3
!
allocate(work(one),iwork(one))
work(1)=real(conf%freq)
status =  nf90_put_var(ncid,rad_freqID,work)
work(1)=real(conf%Theta1)
status =  nf90_put_var(ncid,rad_beamwidthID,work)
work(1)=real(conf%dr)
status =  nf90_put_var(ncid,rad_range_resolutionID,work)
iwork(1)=conf%ixc
status =  nf90_put_var(ncid,ixcID,iwork)
iwork(1)=conf%iyc
status =  nf90_put_var(ncid,iycID,iwork)
work(1)=real(conf%zc)
status =  nf90_put_var(ncid,zcID,work)
deallocate(work,iwork)
!
status =  nf90_put_var(ncid,xID,real(env%x))
status =  nf90_put_var(ncid,yID,real(env%y))
!------------------------------------------------------------
ELSE ! HYDROMETEORS
! write variables
status =  nf90_put_var(ncid,xlatID,real(env%xlat))
status =  nf90_put_var(ncid,xlongID,real(env%xlong))
!
status =  nf90_put_var(ncid,ZhhID,real(rmout%Zhh(:,:,:,iht)))
status =  nf90_put_var(ncid,DoppID,real(rmout%Dopp(:,:,:,iht)))
status =  nf90_put_var(ncid,SWhID,real(rmout%dDVh(:,:,:,iht)))
status =  nf90_put_var(ncid,SWtID,real(rmout%sw_t))
status =  nf90_put_var(ncid,SWsID,real(rmout%sw_s))
status =  nf90_put_var(ncid,SWvID,real(rmout%sw_v))
status =  nf90_put_var(ncid,SWtotID,real(rmout%SWt(:,:,:,iht)))
status =  nf90_put_var(ncid,Dopp90ID,real(rmout%Dopp90(:,:,:,iht)))
status =  nf90_put_var(ncid,SWh90ID,real(rmout%dDVh90(:,:,:,iht)))
status =  nf90_put_var(ncid,DVhID,real(rmout%DVh(:,:,:,iht)))
status =  nf90_put_var(ncid,AhID,real(rmout%Ah(:,:,:,iht)))
!---------------------------------------
status =  nf90_put_var(ncid,elevID,real(rmout%elev))
status =  nf90_put_var(ncid,azimID,real(rmout%azim))
status =  nf90_put_var(ncid,rangeID,real(rmout%range))
status =  nf90_put_var(ncid,zID,real(env%z))
status =  nf90_put_var(ncid,qhydroID,real(qq)) ! kg/m3
!
!---------------------------------------
! MPL
IF (conf%mplID>0) THEN
IF ((iht==1) .or. (iht==3).or.(iht==7)) THEN
!
if (iht==1) ihc=1
if (iht==3) ihc=2
if (iht==7) ihc=3
!
allocate(work(one)); work(1)=real(mpl%wavel)
status =  nf90_put_var(ncid,mpl_wavelID,work)
deallocate(work)
status =  nf90_put_var(ncid,mpl_back_obsID,real(mpl%back_obs(:,:,:,ihc)))
status =  nf90_put_var(ncid,mpl_back_trueID,real(mpl%back_true(:,:,:,ihc)))
status =  nf90_put_var(ncid,mpl_extID,real(mpl%ext(:,:,:,ihc)))
status =  nf90_put_var(ncid,mpl_lidar_ratioID,real(mpl%lidar_ratio(:,:,:,ihc)))
!
ENDIF
ENDIF
!---------------------------------------
!------------------------------------------
! Doppler spectra by oue in May 2018
if(conf%spectraID==1) then
status =  nf90_put_var(ncid,velbinsID,    real(spectra%vel_bins))
status =  nf90_put_var(ncid,spectra_zhhID,real(spectra%zhh_spectra(:,:,:,:,iht)))
end if
!---------------------------
!
allocate(work(one),iwork(one))
work(1)=real(conf%freq)
status =  nf90_put_var(ncid,rad_freqID,work)
work(1)=real(conf%Theta1)
status =  nf90_put_var(ncid,rad_beamwidthID,work)
work(1)=real(conf%dr)
status =  nf90_put_var(ncid,rad_range_resolutionID,work)
iwork(1)=conf%ixc
status =  nf90_put_var(ncid,ixcID,iwork)
iwork(1)=conf%iyc
status =  nf90_put_var(ncid,iycID,iwork)
work(1)=real(conf%zc)
status =  nf90_put_var(ncid,zcID,work)
deallocate(work,iwork)
!
status =  nf90_put_var(ncid,xID,real(env%x))
status =  nf90_put_var(ncid,yID,real(env%y))
!
ENDIF
!------------------------------------------------------------

999 If (status /= 0) then
write(*,*) 'Error in WriteOutNetcdf'
stop
endif
!
status = nf90_close(ncid)
return
end subroutine WriteOutNetcdf_rmp




