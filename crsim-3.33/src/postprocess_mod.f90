!! ----------------------------------------------------------------------------
!! ----------------------------------------------------------------------------
!!  *PROGRAM* arscl_mod
!!  @version crsim 3.x.
!!
!!  *SRC_FILE*
!!  crsim/src/arscl_mod.f90
!!
!!
!!  *LAST CHANGES*
!!
!!  MAY 23 2017  :: M.O.    Implemented ARSCL simulation and MWR LWP
!!
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
module postprocess_mod
Implicit None


!
Type arscl_var ! output arscl variables 
!
Integer                               :: nx,ny,nz,n_layers
!
Integer*2, Dimension(:,:,:),Allocatable :: cloud_detection_flag ! cloud detection flag
Integer*2, Dimension(:,:,:),Allocatable :: cloud_mask ! cloud detection flag
Real*8, Dimension(:,:,:),Allocatable  :: cloud_layer_base ! cloud layer base height [m]
Real*8, Dimension(:,:,:),Allocatable  :: cloud_layer_top  ! cloud layer top height [m]
End type arscl_var
!!
!!
!------------
Type mwr_var ! output mwr variables 
Integer  :: nx,ny,nz
Real*8, Dimension(:,:),Allocatable :: true_lwp ! Model (original) liquid water path [kg m-3]
Real*8, Dimension(:,:),Allocatable :: mwr_lwp ! Microwave radiometer liquid water path [kg m-3]
Integer, Dimension(:,:,:),Allocatable :: n_samples ! Number of samples within MWR field of view
End type mwr_var
!-----------
!!                                                                          
!!!

Contains
!!
subroutine allocate_arscl_var(str)
Implicit None
Type(arscl_var),Intent(InOut)   :: str
Integer                         :: nx,ny,nz,n_layers
!
nx=str%nx
ny=str%ny
nz=str%nz
n_layers = str%n_layers
!
Allocate(str%cloud_detection_flag(nx,ny,nz))
Allocate(str%cloud_mask(nx,ny,nz))
Allocate(str%cloud_layer_base(nx,ny,n_layers))
Allocate(str%cloud_layer_top(nx,ny,n_layers))
!
return
end subroutine allocate_arscl_var
!
!------
subroutine allocate_mwr_var(str)
Implicit None
Type(mwr_var),Intent(InOut)   :: str
Integer                       :: nx,ny,nz
!
nx=str%nx
ny=str%ny
nz=str%nz
!
Allocate(str%true_lwp(nx,ny))
Allocate(str%mwr_lwp(nx,ny))
Allocate(str%n_samples(nx,ny,nz))
!
return
end subroutine allocate_mwr_var
!-----


subroutine m999_arscl_var(str)
Implicit None
Type(arscl_var),Intent(InOut)   :: str
Real*8,Parameter               :: m999=0.d0
!
str%cloud_detection_flag=-999
str%cloud_mask=-999
str%cloud_layer_base=m999
str%cloud_layer_top=m999
!
return
end subroutine m999_arscl_var
!
!------
subroutine m999_mwr_var(str)
Implicit None
Type(mwr_var),Intent(InOut)   :: str
Real*8,Parameter              :: m999=0.d0
!
str%true_lwp=m999
str%mwr_lwp=m999
str%n_samples=0
!
return
end subroutine m999_mwr_var
!-------

subroutine deallocate_arscl_var(str)
Implicit None
Type(arscl_var),Intent(InOut)   :: str
!
Deallocate(str%cloud_detection_flag)
Deallocate(str%cloud_mask)
Deallocate(str%cloud_layer_base)
Deallocate(str%cloud_layer_top)
!
str%nx=0; str%ny=0; str%nz=0 
!
return
end subroutine deallocate_arscl_var

!---------
subroutine deallocate_mwr_var(str)
Implicit None
Type(mwr_var),Intent(InOut)   :: str
!
Deallocate(str%true_lwp)
Deallocate(str%mwr_lwp)
Deallocate(str%n_samples)
!
str%nx=0; str%ny=0; str%nx=0
!
return
end subroutine deallocate_mwr_var
!--------
!!

end module postprocess_mod

