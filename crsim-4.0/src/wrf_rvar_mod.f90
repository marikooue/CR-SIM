!! ----------------------------------------------------------------------------
!! ----------------------------------------------------------------------------
!!  *PROGRAM* wrf_rvar_mod
!!  @version crsim 4.* 
!!
!!  *SRC_FILE*
!!  crsim/src/wrf_rvar_mod.f90
!!
!!
!!  *LAST CHANGES*

!!  *DESCRIPTION* 
!!
!!  This module contains a number of subroutines for allocating, deallocating and
!!  nullifying (or setting to the missed value) data types associated to the
!!  input WRF data fields.
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
!!  GNU  General Public License for more details.
!!
!!  You should have received a copy of the GNU General Public License
!!  along with this program.  If not, see <http://www.gnu.org/licenses/>.
!!
!!----------------------------------------------------------------------------------------------------------------- 
!!-----------------------------------------------------------------------------------------------------------------
!!
module wrf_rvar_mod
Implicit None
!
Type wrf_rvar
!
integer                  :: nt
integer                  :: nx,nxp1
integer                  :: ny,nyp1
integer                  :: nz,nzp1
!
real*4,Dimension(:),Allocatable            :: dx    ! [1/m] -> [m] (inverse) x-grid length
real*4,Dimension(:),Allocatable            :: dy    ! [1/m] -> [m] (inverse) y-grid length
real*4,Dimension(:,:,:,:),Allocatable      :: u     ! [m/s] u horizontal wind componenet 
real*4,Dimension(:,:,:,:),Allocatable      :: v     ! [m/s] v horizontal wind componenet 
real*4,Dimension(:,:,:,:),Allocatable      :: w     ! [m/s] z-wind component
real*4,Dimension(:,:,:,:),Allocatable      :: pb    ! [Pa] base state pressure
real*4,Dimension(:,:,:,:),Allocatable      :: p     ! [Pa] perturbation pressure 
real*4,Dimension(:,:,:,:),Allocatable      :: rho_d ! [m^3/kg]->[kg/m^3] (inverse) dry air density
real*4,Dimension(:,:,:,:),Allocatable      :: theta ! [ K] perturbation potential temperature
real*4,Dimension(:,:,:,:),Allocatable      :: phb   ! [m^2/s^2] base state geopotential
real*4,Dimension(:,:,:,:),Allocatable      :: ph    ! [m^2/s^2] perturbation geopotential
!
real*4,Dimension(:,:,:,:),Allocatable      :: press       ! [mb] pressure
real*4,Dimension(:,:,:,:),Allocatable      :: temp        ! [C] temperature
real*4,Dimension(:,:,:,:),Allocatable      :: geop_height ! [m] geopotential height
real*4,Dimension(:,:),Allocatable      :: topo ! [m] topography height
real*8,Dimension(:),Allocatable        :: hgtm ! [m] height (used for RAMS and SAM)
real*4,Dimension(:,:,:,:),Allocatable      :: qvapor      ! [kg/kg] Water vapor mixing ratio 
!
real*4,Dimension(:,:,:),Allocatable        :: xlat  ! [deg] latitude; south is negative
real*4,Dimension(:,:,:),Allocatable        :: xlong ! [deg] longitude; west is negative
!
real*4,Dimension(:,:,:,:),Allocatable      :: tke   ! [m^2/s^2]  turbulence kinetic energy

real*4,Dimension(:,:,:),Allocatable      :: sfctemp ! [K] Sea surface temperature
real*4,Dimension(:,:,:),Allocatable      :: sfcuu  ! [m/s] surface u
real*4,Dimension(:,:,:),Allocatable      :: sfcvv  ! [m/s] surface v

!
End Type wrf_rvar

Type wrf_rvar_mp09
!
integer                  :: nt
integer                  :: nx
integer                  :: ny
integer                  :: nz
!
real*4,Dimension(:,:,:,:),Allocatable      :: qcloud ! [kg/kg] Cloud water mixing ratio
real*4,Dimension(:,:,:,:),Allocatable      :: qrain  ! [kg/kg] Rain water mixing ratio
real*4,Dimension(:,:,:,:),Allocatable      :: qice   ! [kg/kg] Ice water mixing ratio
real*4,Dimension(:,:,:,:),Allocatable      :: qsnow  ! [kg/kg] Snow water mixing ratio
real*4,Dimension(:,:,:,:),Allocatable      :: qgraup ! [kg/kg] Graupel water mixing ratio
real*4,Dimension(:,:,:,:),Allocatable      :: qhail  ! [kg/kg] Hail water mixing ratio

real*4,Dimension(:,:,:,:),Allocatable      :: qncloud ! [1/kg] Cloud Number concentration
real*4,Dimension(:,:,:,:),Allocatable      :: qnrain  ! [1/kg] Rain Number concentration
real*4,Dimension(:,:,:,:),Allocatable      :: qnice   ! [1/kg] Ice  Number concentration
real*4,Dimension(:,:,:,:),Allocatable      :: qnsnow  ! [1/kg] Snow  Number concentration
real*4,Dimension(:,:,:,:),Allocatable      :: qngraup ! [1/kg] Graupel  Number concentration
real*4,Dimension(:,:,:,:),Allocatable      :: qnhail  ! [1/kg] Hail  Number concentration
!
End Type wrf_rvar_mp09
!!
Type wrf_rvar_mp10
!
integer                  :: nt
integer                  :: nx
integer                  :: ny
integer                  :: nz
!
real*4,Dimension(:,:,:,:),Allocatable      :: qcloud ! [kg/kg] Cloud water mixing ratio
real*4,Dimension(:,:,:,:),Allocatable      :: qrain  ! [kg/kg] Rain water mixing ratio
real*4,Dimension(:,:,:,:),Allocatable      :: qice   ! [kg/kg] Ice water mixing ratio
real*4,Dimension(:,:,:,:),Allocatable      :: qsnow  ! [kg/kg] Snow water mixing ratio
real*4,Dimension(:,:,:,:),Allocatable      :: qgraup ! [kg/kg] Graupel water mixing ratio

real*4,Dimension(:,:,:,:),Allocatable      :: qncloud ! [1/kg] Cloud Number concentration
real*4,Dimension(:,:,:,:),Allocatable      :: qnrain  ! [1/kg] Rain Number concentration
real*4,Dimension(:,:,:,:),Allocatable      :: qnice   ! [1/kg] Ice  Number concentration
real*4,Dimension(:,:,:,:),Allocatable      :: qnsnow  ! [1/kg] Snow  Number concentration
real*4,Dimension(:,:,:,:),Allocatable      :: qngraup ! [1/kg] Graupel  Number concentration
!
End Type wrf_rvar_mp10
!!
!!! Added by oue 2016/09/19--
Type wrf_rvar_mp08
!
integer                  :: nt
integer                  :: nx
integer                  :: ny
integer                  :: nz
!
real*4,Dimension(:,:,:,:),Allocatable      :: qcloud ! [kg/kg] Cloud water mixing ratio
real*4,Dimension(:,:,:,:),Allocatable      :: qrain  ! [kg/kg] Rain water mixing ratio
real*4,Dimension(:,:,:,:),Allocatable      :: qice   ! [kg/kg] Ice water mixing ratio
real*4,Dimension(:,:,:,:),Allocatable      :: qsnow  ! [kg/kg] Snow water mixing ratio
real*4,Dimension(:,:,:,:),Allocatable      :: qgraup ! [kg/kg] Graupel water mixing ratio

real*4,Dimension(:,:,:,:),Allocatable      :: qncloud ! [1/kg] Cloud Number concentration
real*4,Dimension(:,:,:,:),Allocatable      :: qnrain  ! [1/kg] Rain Number concentration
real*4,Dimension(:,:,:,:),Allocatable      :: qnice   ! [1/kg] Ice  Number concentration
real*4,Dimension(:,:,:,:),Allocatable      :: qnsnow  ! [1/kg] Snow  Number concentration
real*4,Dimension(:,:,:,:),Allocatable      :: qngraup ! [1/kg] Graupel  Number concentration
!
End Type wrf_rvar_mp08
!!!-- added by oue
!!
!!! Added by oue 2017/07/17--for ICON
Type wrf_rvar_mp30
!
integer                  :: nt
integer                  :: nx
integer                  :: ny
integer                  :: nz
!
real*4,Dimension(:,:,:,:),Allocatable      :: qcloud ! [kg/kg] Cloud water mixing ratio
real*4,Dimension(:,:,:,:),Allocatable      :: qrain  ! [kg/kg] Rain water mixing ratio
real*4,Dimension(:,:,:,:),Allocatable      :: qice   ! [kg/kg] Ice water mixing ratio
real*4,Dimension(:,:,:,:),Allocatable      :: qsnow  ! [kg/kg] Snow water mixing ratio
real*4,Dimension(:,:,:,:),Allocatable      :: qgraup ! [kg/kg] Graupel water mixing ratio
real*4,Dimension(:,:,:,:),Allocatable      :: qhail  ! [kg/kg] Hail water mixing ratio

real*4,Dimension(:,:,:,:),Allocatable      :: qncloud ! [1/kg] Cloud Number concentration
real*4,Dimension(:,:,:,:),Allocatable      :: qnrain  ! [1/kg] Rain Number concentration
real*4,Dimension(:,:,:,:),Allocatable      :: qnice   ! [1/kg] Ice  Number concentration
real*4,Dimension(:,:,:,:),Allocatable      :: qnsnow  ! [1/kg] Snow  Number concentration
real*4,Dimension(:,:,:,:),Allocatable      :: qngraup ! [1/kg] Graupel  Number concentration
real*4,Dimension(:,:,:,:),Allocatable      :: qnhail  ! [1/kg] Hail Number concentration
!
End Type wrf_rvar_mp30
!!!-- added by oue
!!
!!! Added by oue 2017/07/17--for RAMS
Type wrf_rvar_mp40
!
integer                  :: nt
integer                  :: nx
integer                  :: ny
integer                  :: nz
!
real*4,Dimension(:,:,:,:),Allocatable      :: qcloud ! [kg/kg] Cloud water mixing ratio
real*4,Dimension(:,:,:,:),Allocatable      :: qrain  ! [kg/kg] Rain water mixing ratio
real*4,Dimension(:,:,:,:),Allocatable      :: qice   ! [kg/kg] Ice water mixing ratio
real*4,Dimension(:,:,:,:),Allocatable      :: qsnow  ! [kg/kg] Snow water mixing ratio
real*4,Dimension(:,:,:,:),Allocatable      :: qgraup ! [kg/kg] Graupel water mixing ratio
real*4,Dimension(:,:,:,:),Allocatable      :: qhail  ! [kg/kg] Hail water mixing ratio
real*4,Dimension(:,:,:,:),Allocatable      :: qdrzl  ! [1/kg] Drizzle  Number concentration 
real*4,Dimension(:,:,:,:),Allocatable      :: qaggr  ! [1/kg] Aggregate Number concentration 

real*4,Dimension(:,:,:,:),Allocatable      :: qncloud ! [1/kg] Cloud Number concentration
real*4,Dimension(:,:,:,:),Allocatable      :: qnrain  ! [1/kg] Rain Number concentration
real*4,Dimension(:,:,:,:),Allocatable      :: qnice   ! [1/kg] Ice  Number concentration
real*4,Dimension(:,:,:,:),Allocatable      :: qnsnow  ! [1/kg] Snow  Number concentration
real*4,Dimension(:,:,:,:),Allocatable      :: qngraup ! [1/kg] Graupel  Number concentration
real*4,Dimension(:,:,:,:),Allocatable      :: qnhail  ! [1/kg] Hail Number concentration
real*4,Dimension(:,:,:,:),Allocatable      :: qndrzl  ! [1/kg] Drizzle  Number concentration
real*4,Dimension(:,:,:,:),Allocatable      :: qnaggr  ! [1/kg] Aggregate Number concentration
!
End Type wrf_rvar_mp40
!!!-- added by oue
!!
! Type wrf_rvar_mp50 is added by DW
Type wrf_rvar_mp50
!
integer                  :: nt
integer                  :: nx
integer                  :: ny
integer                  :: nz
!
real*4,Dimension(:,:,:,:),Allocatable      :: qcloud  ! [kg/kg] Cloud water mixing ratio
real*4,Dimension(:,:,:,:),Allocatable      :: qrain   ! [kg/kg] Rain water mixing ratio
real*4,Dimension(:,:,:,:),Allocatable      :: qice    ! [kg/kg] Ice water mixing ratio
real*4,Dimension(:,:,:,:),Allocatable      :: qncloud ! [1/kg]  Cloud Number concentration
real*4,Dimension(:,:,:,:),Allocatable      :: qnrain  ! [1/kg]  Rain Number concentration
real*4,Dimension(:,:,:,:),Allocatable      :: qnice   ! [1/kg]  Ice  Number concentration
real*4,Dimension(:,:,:,:),Allocatable      :: qir     ! [kg/kg] Rime ice mass mixing ratio
real*4,Dimension(:,:,:,:),Allocatable      :: qib     ! [m3/kg] Rime ice volume mixing ratio
!
End Type wrf_rvar_mp50
!!
Type wrf_rvar_mp20
!
integer                  :: nt
integer                  :: nx
integer                  :: ny
integer                  :: nz
integer                  :: nbins
!
! ff1,ff5,ff6 are bin number concentrations multiplied by bin mass and divided by
! density of dry air [kg/kg] = [1/m^3] * [kg]  / [kg/m^3 ]
real*4,Dimension(:,:,:,:,:),Allocatable      :: ff1 ! [kg/kg] cloud/rain bin number concentration
real*4,Dimension(:,:,:,:,:),Allocatable      :: ff5 ! [kg/kg] ice/snow bin number concentration
real*4,Dimension(:,:,:,:,:),Allocatable      :: ff6 ! [kg/kg] graupel  bin number concentration
!
End Type wrf_rvar_mp20
!!
Type wrf_rvar_mp70
!
integer                  :: nt
integer                  :: nx
integer                  :: ny
integer                  :: nz
integer                  :: nbins
!
real*4,Dimension(:,:,:,:,:),Allocatable    :: fm1 ! [g/kg] cloud/rain bin number concentration
real*4,Dimension(:,:,:,:,:),Allocatable    :: fn1 ! [1/cm3] ice/snow bin number concentration
real*4,Dimension(:,:,:,:),Allocatable      :: qcloud  ! [g/kg] Cloud water mixing ratio
real*4,Dimension(:,:,:,:),Allocatable      :: qrain   ! [g/kg] Rain water mixing ratio
real*4,Dimension(:,:,:,:),Allocatable      :: qncloud ! [1/cm3]  Cloud Number concentration
real*4,Dimension(:,:,:,:),Allocatable      :: qnrain  ! [1/cm3]  Rain Number concentration
real*4,Dimension(:,:,:,:),Allocatable      :: naero   ! [1/cm3] Aerosol Number concentration
!
End Type wrf_rvar_mp70
!!
!!
Type wrf_rvar_mp75
!
integer                  :: nt
integer                  :: nx
integer                  :: ny
integer                  :: nz
!
real*4,Dimension(:,:,:,:),Allocatable      :: qcloud ! [g/kg] Cloud water mixing ratio
real*4,Dimension(:,:,:,:),Allocatable      :: qrain  ! [g/kg] Rain water mixing ratio
real*4,Dimension(:,:,:,:),Allocatable      :: qice   ! [g/kg] Ice water mixing ratio
real*4,Dimension(:,:,:,:),Allocatable      :: qsnow  ! [g/kg] Snow water mixing ratio
real*4,Dimension(:,:,:,:),Allocatable      :: qgraup ! [g/kg] Graupel water mixing ratio
real*4,Dimension(:,:,:,:),Allocatable      :: qnoprep ! [g/kg] Non-precipitating Condensate (Water+Ice) 
real*4,Dimension(:,:,:,:),Allocatable      :: qprep ! [g/kg] Precipitating Water (Rain+Snow)

real*4,Dimension(:,:,:,:),Allocatable      :: qncloud ! [1/cm3] Cloud Number concentration
real*4,Dimension(:,:,:,:),Allocatable      :: qnrain  ! [1/cm3] Rain Number concentration
real*4,Dimension(:,:,:,:),Allocatable      :: qnice   ! [1/cm3] Ice  Number concentration
real*4,Dimension(:,:,:,:),Allocatable      :: qnsnow  ! [1/cm3] Snow  Number concentration
real*4,Dimension(:,:,:,:),Allocatable      :: qngraup ! [1/cm3] Graupel  Number concentration
real*4,Dimension(:,:,:,:),Allocatable      :: naero   ! [1/cm3] Aerosol Number concentration 
!
End Type wrf_rvar_mp75
!!
!!
!!!-- Added by oue for CM1 Apr 2020
Type wrf_rvar_mp80
!
integer                  :: nt
integer                  :: nx
integer                  :: ny
integer                  :: nz
!
real*4,Dimension(:,:,:,:),Allocatable      :: qcloud ! [kg/kg] Cloud water mixing ratio
real*4,Dimension(:,:,:,:),Allocatable      :: qrain  ! [kg/kg] Rain water mixing ratio
real*4,Dimension(:,:,:,:),Allocatable      :: qice   ! [kg/kg] Ice water mixing ratio
real*4,Dimension(:,:,:,:),Allocatable      :: qsnow  ! [kg/kg] Snow water mixing ratio
real*4,Dimension(:,:,:,:),Allocatable      :: qgraup ! [kg/kg] Graupel water mixing ratio

real*4,Dimension(:,:,:,:),Allocatable      :: qncloud ! [1/kg] Cloud Number concentration
real*4,Dimension(:,:,:,:),Allocatable      :: qnrain  ! [1/kg] Rain Number concentration
real*4,Dimension(:,:,:,:),Allocatable      :: qnice   ! [1/kg] Ice  Number concentration
real*4,Dimension(:,:,:,:),Allocatable      :: qnsnow  ! [1/kg] Snow  Number concentration
real*4,Dimension(:,:,:,:),Allocatable      :: qngraup ! [1/kg] Graupel  Number concentration
!
End Type wrf_rvar_mp80
!!!- added by oue
!!



Contains

subroutine allocate_wrf_rvar(str)
Implicit None
Type(wrf_rvar),Intent(InOut)     :: str
Integer                         :: nt,nx,ny,nz
Integer                         :: nxp1,nyp1,nzp1
!
nt=str%nt
nx=str%nx ; nxp1=str%nxp1
ny=str%ny ; nyp1=str%nyp1
nz=str%nz ; nzp1=str%nzp1
!
Allocate(str%dx(nt))
Allocate(str%dy(nt))
Allocate(str%u(nxp1,ny,nz,nt))
Allocate(str%v(nx,nyp1,nz,nt))
Allocate(str%w(nx,ny,nzp1,nt))
Allocate(str%pb(nx,ny,nz,nt))
Allocate(str%p(nx,ny,nz,nt))
Allocate(str%rho_d(nx,ny,nz,nt))
Allocate(str%theta(nx,ny,nz,nt))
Allocate(str%phb(nx,ny,nzp1,nt))
Allocate(str%ph(nx,ny,nzp1,nt))
!
Allocate(str%press(nx,ny,nz,nt))
Allocate(str%temp(nx,ny,nz,nt))
Allocate(str%geop_height(nx,ny,nzp1,nt))
Allocate(str%topo(nx,ny))
Allocate(str%hgtm(nz))
!
Allocate(str%qvapor(nx,ny,nz,nt))
!
Allocate(str%xlat(nx,ny,nt))
Allocate(str%xlong(nx,ny,nt))
Allocate(str%tke(nx,ny,nz,nt))
!
! added for surface backscatter Oct 2023
Allocate(str%sfctemp(nx,ny,nt))
Allocate(str%sfcuu(nx,ny,nt))
Allocate(str%sfcvv(nx,ny,nt))
!
return
end subroutine allocate_wrf_rvar
!!
!!
subroutine allocate_wrf_rvar_mp09(str)
Implicit None
Type(wrf_rvar_mp09),Intent(InOut)     :: str
Integer                              :: nt,nx,ny,nz
!
nt=str%nt
nx=str%nx
ny=str%ny
nz=str%nz
!
Allocate(str%qcloud(nx,ny,nz,nt))
Allocate(str%qrain(nx,ny,nz,nt))
Allocate(str%qice(nx,ny,nz,nt))
Allocate(str%qsnow(nx,ny,nz,nt))
Allocate(str%qgraup(nx,ny,nz,nt))
Allocate(str%qhail(nx,ny,nz,nt))
!
Allocate(str%qncloud(nx,ny,nz,nt))
Allocate(str%qnrain(nx,ny,nz,nt))
Allocate(str%qnice(nx,ny,nz,nt))
Allocate(str%qnsnow(nx,ny,nz,nt))
Allocate(str%qngraup(nx,ny,nz,nt))
Allocate(str%qnhail(nx,ny,nz,nt))
!
return
end subroutine allocate_wrf_rvar_mp09
!!
subroutine allocate_wrf_rvar_mp10(str)
Implicit None
Type(wrf_rvar_mp10),Intent(InOut)     :: str
Integer                              :: nt,nx,ny,nz
!
nt=str%nt
nx=str%nx 
ny=str%ny
nz=str%nz
!
Allocate(str%qcloud(nx,ny,nz,nt))
Allocate(str%qrain(nx,ny,nz,nt))
Allocate(str%qice(nx,ny,nz,nt))
Allocate(str%qsnow(nx,ny,nz,nt))
Allocate(str%qgraup(nx,ny,nz,nt))
!
Allocate(str%qncloud(nx,ny,nz,nt))
Allocate(str%qnrain(nx,ny,nz,nt))
Allocate(str%qnice(nx,ny,nz,nt))
Allocate(str%qnsnow(nx,ny,nz,nt))
Allocate(str%qngraup(nx,ny,nz,nt))
!
return
end subroutine allocate_wrf_rvar_mp10
!!
!!! Added by oue 2016/09/19---
subroutine allocate_wrf_rvar_mp08(str)
Implicit None
Type(wrf_rvar_mp08),Intent(InOut)     :: str
Integer                              :: nt,nx,ny,nz
!
nt=str%nt
nx=str%nx 
ny=str%ny
nz=str%nz
!
Allocate(str%qcloud(nx,ny,nz,nt))
Allocate(str%qrain(nx,ny,nz,nt))
Allocate(str%qice(nx,ny,nz,nt))
Allocate(str%qsnow(nx,ny,nz,nt))
Allocate(str%qgraup(nx,ny,nz,nt))
!
Allocate(str%qncloud(nx,ny,nz,nt))
Allocate(str%qnrain(nx,ny,nz,nt))
Allocate(str%qnice(nx,ny,nz,nt))
Allocate(str%qnsnow(nx,ny,nz,nt))
Allocate(str%qngraup(nx,ny,nz,nt))
!
return
end subroutine allocate_wrf_rvar_mp08
!!!-- added by oue
!
!!! Added by oue 2017/07/17---for ICON
subroutine allocate_wrf_rvar_mp30(str)
Implicit None
Type(wrf_rvar_mp30),Intent(InOut)     :: str
Integer                              :: nt,nx,ny,nz
!
nt=str%nt
nx=str%nx 
ny=str%ny
nz=str%nz
!
Allocate(str%qcloud(nx,ny,nz,nt))
Allocate(str%qrain(nx,ny,nz,nt))
Allocate(str%qice(nx,ny,nz,nt))
Allocate(str%qsnow(nx,ny,nz,nt))
Allocate(str%qgraup(nx,ny,nz,nt))
Allocate(str%qhail(nx,ny,nz,nt))
!
Allocate(str%qncloud(nx,ny,nz,nt))
Allocate(str%qnrain(nx,ny,nz,nt))
Allocate(str%qnice(nx,ny,nz,nt))
Allocate(str%qnsnow(nx,ny,nz,nt))
Allocate(str%qngraup(nx,ny,nz,nt))
Allocate(str%qnhail(nx,ny,nz,nt))
!
return
end subroutine allocate_wrf_rvar_mp30
!!!-- added by oue
!
!!! Added by oue 2017/07/21---for RAMS
subroutine allocate_wrf_rvar_mp40(str)
Implicit None
Type(wrf_rvar_mp40),Intent(InOut)     :: str
Integer                              :: nt,nx,ny,nz
!
nt=str%nt
nx=str%nx 
ny=str%ny
nz=str%nz
!
Allocate(str%qcloud(nx,ny,nz,nt))
Allocate(str%qdrzl(nx,ny,nz,nt))
Allocate(str%qrain(nx,ny,nz,nt))
Allocate(str%qice(nx,ny,nz,nt))
Allocate(str%qsnow(nx,ny,nz,nt))
Allocate(str%qgraup(nx,ny,nz,nt))
Allocate(str%qhail(nx,ny,nz,nt))
Allocate(str%qaggr(nx,ny,nz,nt))
!
Allocate(str%qncloud(nx,ny,nz,nt))
Allocate(str%qndrzl(nx,ny,nz,nt))
Allocate(str%qnrain(nx,ny,nz,nt))
Allocate(str%qnice(nx,ny,nz,nt))
Allocate(str%qnsnow(nx,ny,nz,nt))
Allocate(str%qnaggr(nx,ny,nz,nt))
Allocate(str%qngraup(nx,ny,nz,nt))
Allocate(str%qnhail(nx,ny,nz,nt))
!
return
end subroutine allocate_wrf_rvar_mp40
!!!-- added by oue
!
!subroutine allocate_wrf_rvar_mp50(str) is added by DW --for P3
!
subroutine allocate_wrf_rvar_mp50(str)
!
Implicit None
Type(wrf_rvar_mp50),Intent(InOut)     :: str
Integer                              :: nt,nx,ny,nz
!
nt=str%nt
nx=str%nx
ny=str%ny
nz=str%nz
!
Allocate(str%qcloud(nx,ny,nz,nt))
Allocate(str%qrain(nx,ny,nz,nt))
Allocate(str%qice(nx,ny,nz,nt))
!
Allocate(str%qncloud(nx,ny,nz,nt))
Allocate(str%qnrain(nx,ny,nz,nt))
Allocate(str%qnice(nx,ny,nz,nt))
Allocate(str%qir(nx,ny,nz,nt))
Allocate(str%qib(nx,ny,nz,nt))
!
return
end subroutine allocate_wrf_rvar_mp50
!
subroutine allocate_wrf_rvar_mp20(str)
Implicit None
Type(wrf_rvar_mp20),Intent(InOut)     :: str
Integer                              :: nt,nx,ny,nz,nbins
!
nt=str%nt
nx=str%nx
ny=str%ny
nz=str%nz
nbins=str%nbins
!
Allocate(str%ff1(nx,ny,nz,nt,nbins))
Allocate(str%ff5(nx,ny,nz,nt,nbins))
Allocate(str%ff6(nx,ny,nz,nt,nbins))
!
return
end subroutine allocate_wrf_rvar_mp20
!!
!
!!! Added by oue---for SAM
subroutine allocate_wrf_rvar_mp70(str)
Implicit None
Type(wrf_rvar_mp70),Intent(InOut)     :: str
Integer                              :: nt,nx,ny,nz,nbins
!
nt=str%nt
nx=str%nx 
ny=str%ny
nz=str%nz
nbins=str%nbins
!
Allocate(str%fm1(nx,ny,nz,nt,nbins))
Allocate(str%fn1(nx,ny,nz,nt,nbins))
Allocate(str%qcloud (nx,ny,nz,nt))
Allocate(str%qrain  (nx,ny,nz,nt))
Allocate(str%qncloud(nx,ny,nz,nt))
Allocate(str%qnrain (nx,ny,nz,nt))
Allocate(str%naero  (nx,ny,nz,nt))
!
return
end subroutine allocate_wrf_rvar_mp70
!!!-- added by oue
subroutine allocate_wrf_rvar_mp75(str)
Implicit None
Type(wrf_rvar_mp75),Intent(InOut)     :: str
Integer                              :: nt,nx,ny,nz
!
nt=str%nt
nx=str%nx 
ny=str%ny
nz=str%nz
!
Allocate(str%qnoprep(nx,ny,nz,nt))
Allocate(str%qprep(nx,ny,nz,nt))
Allocate(str%qcloud(nx,ny,nz,nt))
Allocate(str%qrain(nx,ny,nz,nt))
Allocate(str%qice(nx,ny,nz,nt))
Allocate(str%qsnow(nx,ny,nz,nt))
Allocate(str%qgraup(nx,ny,nz,nt))
!
Allocate(str%qncloud(nx,ny,nz,nt))
Allocate(str%qnrain(nx,ny,nz,nt))
Allocate(str%qnice(nx,ny,nz,nt))
Allocate(str%qnsnow(nx,ny,nz,nt))
Allocate(str%qngraup(nx,ny,nz,nt))
Allocate(str%naero  (nx,ny,nz,nt))
!
return
end subroutine allocate_wrf_rvar_mp75

!
!!
!!! Added by our for CM1 Apr 2020
subroutine allocate_wrf_rvar_mp80(str)
Implicit None
Type(wrf_rvar_mp80),Intent(InOut)     :: str
Integer                              :: nt,nx,ny,nz
!
nt=str%nt
nx=str%nx 
ny=str%ny
nz=str%nz
!
Allocate(str%qcloud(nx,ny,nz,nt))
Allocate(str%qrain(nx,ny,nz,nt))
Allocate(str%qice(nx,ny,nz,nt))
Allocate(str%qsnow(nx,ny,nz,nt))
Allocate(str%qgraup(nx,ny,nz,nt))
!
Allocate(str%qncloud(nx,ny,nz,nt))
Allocate(str%qnrain(nx,ny,nz,nt))
Allocate(str%qnice(nx,ny,nz,nt))
Allocate(str%qnsnow(nx,ny,nz,nt))
Allocate(str%qngraup(nx,ny,nz,nt))
!
return
end subroutine allocate_wrf_rvar_mp80
!!!-added by oue

!!
subroutine initialize_wrf_rvar(str)
Implicit None
Type(wrf_rvar),Intent(InOut)     :: str
ReaL*4,parameter                :: m999=-999.e0
!
str%dx=m999
str%dy=m999
str%u=m999
str%v=m999
str%w=m999
str%pb=m999
str%p=m999
str%rho_d=m999
str%theta=m999
str%phb=m999
str%ph=m999
!
str%press=m999
str%temp=m999
str%geop_height=m999
str%topo=0.e0
str%hgtm=0.e0
!
str%qvapor=m999
!
str%xlat=m999
str%xlong=m999
str%tke=m999
!
! added for surface backscatter Oct. 2023
str%sfctemp=m999
str%sfcuu=m999
str%sfcvv=m999
!
return
end subroutine initialize_wrf_rvar
!!
!!
subroutine initialize_wrf_rvar_mp09(str)
Implicit None
Type(wrf_rvar_mp09),Intent(InOut)     :: str
ReaL*4,parameter                     :: m999=-999.e0
!
str%qcloud=m999
str%qrain=m999
str%qice=m999
str%qsnow=m999
str%qgraup=m999
str%qhail=m999
!
str%qncloud=m999
str%qnrain=m999
str%qnice=m999
str%qnsnow=m999
str%qngraup=m999
str%qnhail=m999
!
return
end subroutine initialize_wrf_rvar_mp09
!!
subroutine initialize_wrf_rvar_mp10(str)
Implicit None
Type(wrf_rvar_mp10),Intent(InOut)     :: str
ReaL*4,parameter                     :: m999=-999.e0     
!
str%qcloud=m999
str%qrain=m999
str%qice=m999
str%qsnow=m999
str%qgraup=m999
!
str%qncloud=m999
str%qnrain=m999
str%qnice=m999
str%qnsnow=m999
str%qngraup=m999
!
return
end subroutine initialize_wrf_rvar_mp10
!!
!!! Added by oue 2016/09/19 --
subroutine initialize_wrf_rvar_mp08(str)
Implicit None
Type(wrf_rvar_mp08),Intent(InOut)     :: str
ReaL*4,parameter                     :: m999=-999.e0     
!
str%qcloud=m999
str%qrain=m999
str%qice=m999
str%qsnow=m999
str%qgraup=m999
!
str%qncloud=m999
str%qnrain=m999
str%qnice=m999
str%qnsnow=m999
str%qngraup=m999
!
return
end subroutine initialize_wrf_rvar_mp08
!!!-- added by oue
!
!!! Added by oue 2017/07/17 --for ICON
subroutine initialize_wrf_rvar_mp30(str)
Implicit None
Type(wrf_rvar_mp30),Intent(InOut)     :: str
ReaL*4,parameter                     :: m999=-999.e0     
!
str%qcloud=m999
str%qrain=m999
str%qice=m999
str%qsnow=m999
str%qgraup=m999
str%qhail=m999
!
str%qncloud=m999
str%qnrain=m999
str%qnice=m999
str%qnsnow=m999
str%qngraup=m999
str%qnhail=m999
!
return
end subroutine initialize_wrf_rvar_mp30
!!!-- added by oue
!
!!! Added by oue 2017/07/21 --for RAMS
subroutine initialize_wrf_rvar_mp40(str)
Implicit None
Type(wrf_rvar_mp40),Intent(InOut)     :: str
ReaL*4,parameter                     :: m999=-999.e0     
!
str%qcloud=m999
str%qrain=m999
str%qice=m999
str%qsnow=m999
str%qgraup=m999
str%qhail=m999
str%qdrzl=m999
str%qaggr=m999
!
str%qncloud=m999
str%qnrain=m999
str%qnice=m999
str%qnsnow=m999
str%qngraup=m999
str%qnhail=m999
str%qndrzl=m999
str%qnaggr=m999
!
return
end subroutine initialize_wrf_rvar_mp40
!!!-- added by oue
!
!subroutine initialize_wrf_rvar_mp50(str) is added by DW --for P3
!
subroutine initialize_wrf_rvar_mp50(str)
!
Implicit None
Type(wrf_rvar_mp50),Intent(InOut)     :: str
ReaL*4,parameter                     :: m999=-999.e0
!
str%qcloud=m999
str%qrain=m999
str%qice=m999
!
str%qncloud=m999
str%qnrain=m999
str%qnice=m999
str%qir=m999
str%qib=m999
!
return
end subroutine initialize_wrf_rvar_mp50
!------------------------

subroutine initialize_wrf_rvar_mp20(str)
Implicit None
Type(wrf_rvar_mp20),Intent(InOut)     :: str
ReaL*4,parameter                     :: m999=-999.e0
!
str%ff1=m999
str%ff5=m999
str%ff6=m999
!
return
end subroutine initialize_wrf_rvar_mp20
!!
subroutine initialize_wrf_rvar_mp70(str)
Implicit None
Type(wrf_rvar_mp70),Intent(InOut)     :: str
ReaL*4,parameter                     :: m999=-999.e0
!
str%fm1=m999
str%fn1=m999
str%qcloud=m999
str%qrain=m999
str%qncloud=m999
str%qnrain=m999
str%naero=m999
!
return
end subroutine initialize_wrf_rvar_mp70
!!
subroutine initialize_wrf_rvar_mp75(str)
Implicit None
Type(wrf_rvar_mp75),Intent(InOut)     :: str
ReaL*4,parameter                     :: m999=-999.e0
!
str%qcloud=m999
str%qrain=m999
str%qice=m999
str%qsnow=m999
str%qgraup=m999
str%qnoprep=m999
str%qprep=m999
!  
str%qncloud=m999
str%qnrain=m999
str%qnice=m999
str%qnsnow=m999
str%qngraup=m999
str%naero=m999
!
return
end subroutine initialize_wrf_rvar_mp75
!!
!!! Added by oue for CM1 Apr 2020
subroutine initialize_wrf_rvar_mp80(str)
Implicit None
Type(wrf_rvar_mp80),Intent(InOut)     :: str
ReaL*4,parameter                     :: m999=-999.e0     
!
str%qcloud=m999
str%qrain=m999
str%qice=m999
str%qsnow=m999
str%qgraup=m999
!
str%qncloud=m999
str%qnrain=m999
str%qnice=m999
str%qnsnow=m999
str%qngraup=m999
!
return
end subroutine initialize_wrf_rvar_mp80
!!

!!
subroutine deallocate_wrf_rvar(str)
Implicit None
Type(wrf_rvar),Intent(InOut)     :: str
!
Deallocate(str%dx)
Deallocate(str%dy)
Deallocate(str%u)
Deallocate(str%v)
Deallocate(str%w)
Deallocate(str%pb)
Deallocate(str%p)
Deallocate(str%rho_d)
Deallocate(str%theta)
Deallocate(str%phb)
Deallocate(str%ph)
!
Deallocate(str%press)
Deallocate(str%temp)
Deallocate(str%geop_height)
Deallocate(str%topo)
Deallocate(str%hgtm)
!
Deallocate(str%qvapor)
!
Deallocate(str%xlat)
Deallocate(str%xlong)
Deallocate(str%tke)
!
!added for surface backscatter Oct. 2023
Deallocate(str%sfctemp)
Deallocate(str%sfcuu)
Deallocate(str%sfcvv)
!
str%nt=0
str%nx=0 ; str%nxp1=0
str%ny=0 ; str%nyp1=0
str%nz=0 ; str%nzp1=0
!
return
end subroutine deallocate_wrf_rvar
!!
!!
subroutine deallocate_wrf_rvar_mp09(str)
Implicit None
Type(wrf_rvar_mp09),Intent(InOut)     :: str
!
Deallocate(str%qcloud)
Deallocate(str%qrain)
Deallocate(str%qice)
Deallocate(str%qsnow)
Deallocate(str%qgraup)
Deallocate(str%qhail)
!
Deallocate(str%qncloud)
Deallocate(str%qnrain)
Deallocate(str%qnice)
Deallocate(str%qnsnow)
Deallocate(str%qngraup)
Deallocate(str%qnhail)
!
str%nt=0
str%nx=0
str%ny=0
str%nz=0
!
return
end subroutine deallocate_wrf_rvar_mp09
!!
subroutine deallocate_wrf_rvar_mp10(str)
Implicit None
Type(wrf_rvar_mp10),Intent(InOut)     :: str
!
Deallocate(str%qcloud)
Deallocate(str%qrain)
Deallocate(str%qice)
Deallocate(str%qsnow)
Deallocate(str%qgraup)
!
Deallocate(str%qncloud)
Deallocate(str%qnrain)
Deallocate(str%qnice)
Deallocate(str%qnsnow)
Deallocate(str%qngraup)
!
str%nt=0
str%nx=0 
str%ny=0 
str%nz=0 
!
return
end subroutine deallocate_wrf_rvar_mp10
!
!!Added by oue 2016/9/19 --
subroutine deallocate_wrf_rvar_mp08(str)
Implicit None
Type(wrf_rvar_mp08),Intent(InOut)     :: str
!
Deallocate(str%qcloud)
Deallocate(str%qrain)
Deallocate(str%qice)
Deallocate(str%qsnow)
Deallocate(str%qgraup)
!
Deallocate(str%qncloud)
Deallocate(str%qnrain)
Deallocate(str%qnice)
Deallocate(str%qnsnow)
Deallocate(str%qngraup)
!
str%nt=0
str%nx=0 
str%ny=0 
str%nz=0 
!
return
end subroutine deallocate_wrf_rvar_mp08
!!--added by oue
!
!!Added by oue 2017/7/17 --for ICON
subroutine deallocate_wrf_rvar_mp30(str)
Implicit None
Type(wrf_rvar_mp30),Intent(InOut)     :: str
!
Deallocate(str%qcloud)
Deallocate(str%qrain)
Deallocate(str%qice)
Deallocate(str%qsnow)
Deallocate(str%qgraup)
Deallocate(str%qhail)
!
Deallocate(str%qncloud)
Deallocate(str%qnrain)
Deallocate(str%qnice)
Deallocate(str%qnsnow)
Deallocate(str%qngraup)
Deallocate(str%qnhail)
!
str%nt=0
str%nx=0 
str%ny=0 
str%nz=0 
!
return
end subroutine deallocate_wrf_rvar_mp30
!!--added by oue
!
!!Added by oue 2017/7/21 --for RAMS
subroutine deallocate_wrf_rvar_mp40(str)
Implicit None
Type(wrf_rvar_mp40),Intent(InOut)     :: str
!
Deallocate(str%qcloud)
Deallocate(str%qrain)
Deallocate(str%qice)
Deallocate(str%qsnow)
Deallocate(str%qgraup)
Deallocate(str%qhail)
Deallocate(str%qdrzl)
Deallocate(str%qaggr)
!
Deallocate(str%qncloud)
Deallocate(str%qnrain)
Deallocate(str%qnice)
Deallocate(str%qnsnow)
Deallocate(str%qngraup)
Deallocate(str%qnhail)
Deallocate(str%qndrzl)
Deallocate(str%qnaggr)
!
str%nt=0
str%nx=0 
str%ny=0 
str%nz=0 
!
return
end subroutine deallocate_wrf_rvar_mp40
!!--added by oue
!
!subroutine deallocate_wrf_rvar_mp50(str) is added by DW --for P3
subroutine deallocate_wrf_rvar_mp50(str)
!
Implicit None
Type(wrf_rvar_mp50),Intent(InOut)     :: str
!
Deallocate(str%qcloud)
Deallocate(str%qrain)
Deallocate(str%qice)
!
Deallocate(str%qncloud)
Deallocate(str%qnrain)
Deallocate(str%qnice)
Deallocate(str%qir)
Deallocate(str%qib)
!
str%nt=0
str%nx=0
str%ny=0
str%nz=0
!
return
end subroutine deallocate_wrf_rvar_mp50
!--------------------------

subroutine deallocate_wrf_rvar_mp20(str)
Implicit None
Type(wrf_rvar_mp20),Intent(InOut)     :: str
!
Deallocate(str%ff1)
Deallocate(str%ff5)
Deallocate(str%ff6)
!
str%nt=0
str%nx=0
str%ny=0
str%nz=0
str%nbins=0
!
return
end subroutine deallocate_wrf_rvar_mp20
!!

subroutine deallocate_wrf_rvar_mp70(str)
Implicit None
Type(wrf_rvar_mp70),Intent(InOut)     :: str
!
Deallocate(str%fm1)
Deallocate(str%fn1)
Deallocate(str%qcloud)
Deallocate(str%qrain)
Deallocate(str%qncloud)
Deallocate(str%qnrain)
Deallocate(str%naero)

str%nt=0
str%nx=0
str%ny=0
str%nz=0
str%nbins=0
!
return
end subroutine deallocate_wrf_rvar_mp70
!!
subroutine deallocate_wrf_rvar_mp75(str)
Implicit None
Type(wrf_rvar_mp75),Intent(InOut)     :: str
!
Deallocate(str%qcloud)
Deallocate(str%qrain)
Deallocate(str%qice)
Deallocate(str%qsnow)
Deallocate(str%qnoprep)
Deallocate(str%qprep)
!
Deallocate(str%qncloud)
Deallocate(str%qnrain)
Deallocate(str%qnice)
Deallocate(str%qnsnow)
Deallocate(str%qngraup)
Deallocate(str%naero)
!
str%nt=0
str%nx=0 
str%ny=0 
str%nz=0 
!
return
end subroutine deallocate_wrf_rvar_mp75
!
!!! Added by oue for CM1 Apr 2020
!!
subroutine deallocate_wrf_rvar_mp80(str)
Implicit None
Type(wrf_rvar_mp80),Intent(InOut)     :: str
!
Deallocate(str%qcloud)
Deallocate(str%qrain)
Deallocate(str%qice)
Deallocate(str%qsnow)
Deallocate(str%qgraup)
!
Deallocate(str%qncloud)
Deallocate(str%qnrain)
Deallocate(str%qnice)
Deallocate(str%qnsnow)
Deallocate(str%qngraup)
!
str%nt=0
str%nx=0 
str%ny=0 
str%nz=0 
!
return
end subroutine deallocate_wrf_rvar_mp80
!
!!
end module wrf_rvar_mod
!!
