!! ----------------------------------------------------------------------------
!! ----------------------------------------------------------------------------
!!  *PROGRAM* ReadConfParameters
!!  @version crsim 3.x
!!
!!  *SRC_FILE*
!!  crsim/crsim/ReadConfParameters.f90
!!
!!
!!
!!  *LAST CHANGES*
!!
!! Jan 12 2016             - A. T.   included new parameters in Configuration file, the minimum values of the input mixing ratios
!!                                   for species. Mixing ration values smaller or equal than specified thresholds are set to 0.
!! Jan 13 2016             - A. T.   included new conf parameters: radar beamwidth and range resolution
!! Feb 15 2016             - A. T.   included new conf parameter : coefficient ZMIN in equation DBZmin [dBZ]=ZMIN [dBZ]+ 20 log10 (radar range [km])
!! Feb 17 2016             - A. T.   included new conf parameter : ceiloID (=1 introduce ceilometer vertical measurements)
!! Apr 29 2016             - A. T.   included new conf parameter : mplID (=1 or 2  introduce MPL vertical measurements at 353 and 532 nm respectivelly,
!!                                    <=0 do not do it). 
!! May 06 2016             - A. T.   introduced new MPL parameters: aeroID, aero_tau and aero_lidar_ratio
!! Mar 23 2017             - M. O.   introduced ARSCL and MWR LWP parameters: arsclID/mwrID, mwr_view (MWR firld of view)
!! JUL 17 2017             - M. O.   Added MP_PHYSICS=30 for ICON (6 categories inclusing hail)
!! JUL 21 2017             - M. O.   Added MP_PHYSICS=40 for RAMS (8 categories inclusing hail)
!! OCT 11 2017             - K. YU.  Added OpenMP thread number settging. It has more priority to the environment variable OMP_NUM_THREADS
!! Oct 30 2017             - D. W.   incoporated P3 (if str%MP_PHYSICS==50, str%nht=6) (added by DW)
!! Apr    2018             - M.O     Incorporated SAM warm bin microphysics (MP_PHYSICS=70)
!! MAY 27 2018             - M.O     Incorporated Doppler spectra generator for bulk moment microphysics (spectraID)
!!                                   Added an option to use hail instead of graupel for MP_PHYSICS=10 (WRF Morrison 2 moment)!!
!! Aug    2019             - A.T     Added new paramater airborne if ixc <= -999 or iyc<=-999
!!                                   for down looking radar at aeroplane flying  along x axis at y=yc and scanning perpendicularly from -90  to +90
!!                                   or aeroplane flying  along y axis at x=xc and scanning perpendicularly from -90  to +90
!!                                   with zero elevation angle in vertical
!! Oct    2023             -M.O.     Incorporated airborne radar simulator
!!
!!  *DESCRIPTION* 
!!
!!  This program  contains a  subroutine needed for reading information 
!!  given in the User Parameter File "PARAMETERS" 
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
!!  You should have received a copy of the GNU General Public License
!!  along with this program.  If not, see <http://www.gnu.org/licenses/>.
!!
!!----------------------------------------------------------------------------------------------------------------- 
!!-----------------------------------------------------------------------------------------------------------------
!!
subroutine ReadConfParameters(InpFile,str)
Use crsim_mod
Implicit None
!
Character(len=*),Intent(in)                              :: InpFile
Type(conf_var),Intent(out)                               :: str
Integer                                                  :: i
!
integer,parameter                                        :: id=1
character(len=220)                                       :: inputline='#'
character(len=940)                                       :: PWD
!
integer                                                  :: work
!
!
open(unit=id,file=InpFile)
!
inputline='#'
do while (inputline(1:1).eq.'#')
read(id,'(a)') inputline
enddo
read(inputline,*) str%OMPThreadNum ! Number of OpenMP threads
!
inputline='#'
do while (inputline(1:1).eq.'#')
read(id,'(a)') inputline
enddo
read(inputline,*) str%ix_start,str%ix_end
!
inputline='#'
do while (inputline(1:1).eq.'#')
read(id,'(a)') inputline
enddo
read(inputline,*) str%iy_start,str%iy_end
!
inputline='#'
do while (inputline(1:1).eq.'#')
read(id,'(a)') inputline
enddo
read(inputline,*) str%iz_start, str%iz_end
!
inputline='#'
do while (inputline(1:1).eq.'#')
read(id,'(a)') inputline
enddo
read(inputline,*) str%it
!
inputline='#'
do while (inputline(1:1).eq.'#')
read(id,'(a)') inputline
enddo
read(inputline,*) work
!
str%MP_PHYSICS=work
str%snow_spherical=1
str%MP10HailOption=0

if (work==901) then
   str%MP_PHYSICS=9
   str%snow_spherical=0
endif
if (work==101) then
   str%MP_PHYSICS=10
   str%MP10HailOption=1
endif

str%nht=5
if (str%MP_PHYSICS==9) str%nht=6
if (str%MP_PHYSICS==30) str%nht=6 !added by oue 2017/07/17 for ICON
if (str%MP_PHYSICS==40) str%nht=8 !added by oue 2017/07/21 for RAMS
if (str%MP_PHYSICS==50) str%nht=6 !added by DW  2017/10/30 for P3
if (str%MP_PHYSICS==70) str%nht=2 !added by oue 2018 April, for SAM warm bin
!
inputline='#'
do while (inputline(1:1).eq.'#')
read(id,'(a)') inputline
enddo
read(inputline,*) str%ixc,str%iyc

! AUG2019
str%airborne = 0 ! by default radar is looking up 
!check if airborne
if ( (str%ixc <= -999) .or. (str%iyc <= -999) ) then
  ! airborne flying along x axis at y=yc and scanning perpendicularly from -90 to +90   
  ! or airborne flying along x axis at x=yc and scanning perpendicularly from
  ! -90 to +90  
  str%airborne = 1
  !write(*,*) 'the option is not implemented yet'
  !stop
endif
!
inputline='#'
do while (inputline(1:1).eq.'#')
read(id,'(a)') inputline
enddo
read(inputline,*) str%zc
!
str%hydro_luts(:)="NULL" ! init 
do i=1,str%nht
inputline='#'
do while (inputline(1:1).eq.'#')
read(id,'(a)') inputline
enddo
read(inputline,*) str%hydro_luts(i)
enddo
!
str%thr_mix_ratio(:)=0.d0 ! init
do i=1,str%nht
inputline='#'
do while (inputline(1:1).eq.'#')
read(id,'(a)') inputline
enddo
read(inputline,*) str%thr_mix_ratio(i)
enddo
!
str%horientID(:)=3 ! init 
str%sigma(:)=40.d0 ; str%sigma(1:3)=10.d0 ! init
do i=1,str%nht
inputline='#'
do while (inputline(1:1).eq.'#')
read(id,'(a)') inputline
enddo
read(inputline,*) str%horientID(i),str%sigma(i)
enddo
!
where (str%horientID<0) str%horientID=3
where (str%sigma(1:3)<0) str%sigma(1:3)=10.d0
where (str%sigma(4:str%nht)<0) str%sigma(4:str%nht)=40.d0 
!
inputline='#'
do while (inputline(1:1).eq.'#')
read(id,'(a)') inputline
enddo
read(inputline,*) str%freq
!
inputline='#'
do while (inputline(1:1).eq.'#')
read(id,'(a)') inputline
enddo
read(inputline,*) str%elev
!
inputline='#'
do while (inputline(1:1).eq.'#')
read(id,'(a)') inputline
enddo
read(inputline,*) str%radID
!
inputline='#'
do while (inputline(1:1).eq.'#')
read(id,'(a)') inputline
enddo
read(inputline,*) str%Theta1
!
inputline='#'
do while (inputline(1:1).eq.'#')
read(id,'(a)') inputline
enddo
read(inputline,*) str%dr
!
inputline='#'
do while (inputline(1:1).eq.'#')
read(id,'(a)') inputline
enddo
read(inputline,*) str%ZMIN
!
inputline='#'
do while (inputline(1:1).eq.'#')
read(id,'(a)') inputline
enddo
read(inputline,*) str%ceiloID
!
inputline='#'
do while (inputline(1:1).eq.'#')
read(id,'(a)') inputline
enddo
read(inputline,*) str%mplID
!
inputline='#'
do while (inputline(1:1).eq.'#')
read(id,'(a)') inputline
enddo
read(inputline,*) str%aeroID
!
inputline='#'
do while (inputline(1:1).eq.'#')
read(id,'(a)') inputline
enddo
read(inputline,*) str%aero_tau
!
inputline='#'
do while (inputline(1:1).eq.'#')
read(id,'(a)') inputline
enddo
read(inputline,*) str%aero_lidar_ratio
!
!- Doppler spectra ID (spectraID): ==1: simulate Doppler spectra; /=1: not simulate Doppler spectra
inputline='#'
do while (inputline(1:1).eq.'#')
read(id,'(a)') inputline
enddo
read(inputline,*) str%spectraID
!
!- airborne added by oue Oct 2023
!- Airborne simulation ID (airborne): ==1: simulate airborne radar variables (surface backscatter and airplane radial velocity; /=1: not simulate airborne radar variables
inputline='#'
do while (inputline(1:1).eq.'#')
read(id,'(a)') inputline
enddo
read(inputline,*) str%airborne
!- Airborne simulation: airplane speed (m/s) and az/el directions (degrees from north/ degrees from horizon)
str%airborne_spd=0.0
str%airborne_azdeg=0.0
str%airborne_eldeg=0.0
inputline='#'
do while (inputline(1:1).eq.'#')
read(id,'(a)') inputline
enddo
read(inputline,*) str%airborne_spd,str%airborne_azdeg,str%airborne_eldeg
!- Airborne simulation: pulse length (m). default is 100 m
str%pulse_len = 100
inputline='#'
do while (inputline(1:1).eq.'#')
read(id,'(a)') inputline
enddo
read(inputline,*) str%pulse_len
!
!- ARSCL parameter added by oue 2017.03.23
inputline='#'
do while (inputline(1:1).eq.'#')
read(id,'(a)') inputline
enddo
read(inputline,*) str%arsclID
!
!- Microwave radiometer liquid water path  added by oue 2017.03.24
inputline='#'
do while (inputline(1:1).eq.'#')
read(id,'(a)') inputline
enddo
read(inputline,*) str%mwrID
!
if (str%mwrID == 1) then
inputline='#'
do while (inputline(1:1).eq.'#')
read(id,'(a)') inputline
enddo
read(inputline,*) str%mwr_view
!
inputline='#'
do while (inputline(1:1).eq.'#')
read(id,'(a)') inputline
enddo
read(inputline,*) str%mwr_alt
endif
!
!
!
inputline='#'


close(id)
!
str%sigma_r=0.35d0*str%dr
str%sigma_theta=str%Theta1/(4.d0*dsqrt(dlog(2.d0)))
!
str%mpl_wavel=0.d0
if (str%mplID==1) str%mpl_wavel=353.d0
if (str%mplID==2) str%mpl_wavel=532.d0
!
! define important paths
!
! INIT
str%path_back_scatt_libs="NULL"
str%file_ice_parameter="NULL"
do i=1,4
str%WRF_MP_PHYSICS_20_InputFiles(i)="NULL"
enddo
!
! LUTS
call getenv('PWD',pwd)
str%path_back_scatt_libs=Trim(PWD)//"/../aux/LLUTS3/"
write(*,*) 'LUTS=',Trim(str%path_back_scatt_libs)
!
! file_ice_parameter
if (str%MP_PHYSICS==50) then !added by DW  2017/10/30 for P3
str%file_ice_parameter=Trim(str%path_back_scatt_libs)//"../P3/p3_ice_psd_para_lookup_table.dat"
endif
! 
! WRF_MP_PHYSICS_20_InputFiles
if (str%MP_PHYSICS==20) then
str%WRF_MP_PHYSICS_20_InputFiles(1)= Trim(str%path_back_scatt_libs)//"../SBM_20/bulkradii.asc_s_0_03_0_9" ! diams
str%WRF_MP_PHYSICS_20_InputFiles(2)= Trim(str%path_back_scatt_libs)//"../SBM_20/bulkdens.asc_s_0_03_0_9" ! densities
str%WRF_MP_PHYSICS_20_InputFiles(3)=Trim(str%path_back_scatt_libs)//"../SBM_20/masses.asc"  ! mass
str%WRF_MP_PHYSICS_20_InputFiles(4)=Trim(str%path_back_scatt_libs)//"../SBM_20/termvels.asc"  ! fall velocities
endif

!==========================================================
!- Massages added by oue April 2018
     if (str%MP_PHYSICS==8) then
write(*,*) 'Selected microphysics: 8 (WRF Thompson bulk moment)'
else if (str%MP_PHYSICS==9) then
write(*,*) 'Selected microphysics: 9 (WRF Milbranndt-Yau)'
if(str%snow_spherical==0) write(*,*) ' Non spherical snow'
else if (str%MP_PHYSICS==10) then
write(*,*) 'Selected microphysics: 10 (WRF Morrison 2 moment)'
if(str%MP10HailOption==1) write(*,*) ' Use hail instead of graupel'
else if (str%MP_PHYSICS==20) then
write(*,*) 'Selected microphysics: 20 (WRF Fan bin microphysics)'
else if (str%MP_PHYSICS==30) then
write(*,*) 'Selected microphysics: 30 (ICON 2 moment)'
else if (str%MP_PHYSICS==40) then
write(*,*) 'Selected microphysics: 40 (RAMS 2 moment)'
else if (str%MP_PHYSICS==50) then
write(*,*) 'Selected microphysics: 50 (WRF P3)'
else if (str%MP_PHYSICS==70) then
write(*,*) 'Selected microphysics: 70 (SAM warm bin)'
else if (str%MP_PHYSICS==75) then
write(*,*) 'Selected microphysics: 75 (SAM Morrison 2 moment)'
else if (str%MP_PHYSICS==80) then
write(*,*) 'Selected microphysics: 80 (CM1 Morrison 2 moment)'
else
write(*,*) 'Unknown microphysics:',str%MP_PHYSICS
end if
if(str%radID /=1)     write(*,*) 'Include polarimetry'
if(str%spectraID ==1) write(*,*) 'Include Doppler spectra'
if(str%airborne ==1)  then ! airborne added by oue Oct 2023
   write(*,*) 'Include airborne radar: the option is temporary inplemented (evaluation stage)'
   write(*,*) 'Airplane speed: ',str%airborne_spd,'m/s'
   write(*,*) 'Airplane azimuth direction: ',str%airborne_azdeg,'deg'
   write(*,*) 'Airplane elevation direction: ',str%airborne_eldeg,'deg'
   write(*,*) 'Airplane radar pulse length: ',str%pulse_len,'m'
end if
!==========================================================

!

return
end subroutine ReadConfParameters

