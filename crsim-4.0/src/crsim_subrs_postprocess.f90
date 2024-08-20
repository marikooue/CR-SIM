 !! ----------------------------------------------------------------------------
 !! ----------------------------------------------------------------------------
 !!  *PROGRAM* crsim_subrs_postprocess 
 !!  @version crsim 4.x
 !!
 !!  *SRC_FILE*
 !!  crsim/src/crsim_subrs_arscl.f90
 !!
 !!
 !!
 !!  *LAST CHANGES*
 !!
 !! Mar  23  2017   -M.O     Implemented Post processing: ARSCL simulation and MWR LWP
 !! Aug 20   2019   -A.T     Corrected precission issues (real to double precission)
 !! Aug 21   2019   -A.T     Introduced indentation.
 !!
 !!  *DESCRIPTION* 
 !!
 !!  This program  contains the subroutines needed for postp processing: 
 !!  computation of ARSCL products and MWR LWP
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
 !!  You should have received a copy of the GNU  General Public License
 !!  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 !!
 !!
 !!----------------------------------------------------------------------------------------------------------------- 
 !!-----------------------------------------------------------------------------------------------------------------
 !!
 subroutine estimate_arscl_cloudmask(nz,zz,n_layers,&
                           radar_zhh_true,radar_zhh_min,radar_zhh_atten,&
                           mpl_back_hydro, mpl_back_obs,mpl_back_aerosol,mpl_back_rayleigh,&
                           first_cloud_base,&
                           cloud_detection_flag,cloud_mask,cloud_layer_base,cloud_layer_top)
 ! Lidar assumed to be at the first grid level
 Implicit None
 !!
 integer, Intent(in)   :: nz
 Real*8, Intent(In)    :: zz(nz) ! height [m]
 integer, Intent(in)   :: n_layers ! maximum nimber of cloud layers
 Real*8, Intent(In)    :: radar_zhh_true(nz) ! true radar reflectivity (vertical pointing) [dBZ] 
 Real*8, Intent(In)    :: radar_zhh_min(nz)  ! radar sesitivity [dBZ] 
 Real*8, Intent(In)    :: radar_zhh_atten(nz) ! radar specific attenuation [dB km-1] 
 Real*8, Intent(In)    :: mpl_back_hydro(nz)    ! lidar hydrometeor backscatter (no attenuation)[m-1 sr-1] 
 Real*8, Intent(In)    :: mpl_back_obs(nz)      ! lidar observed backscatter (back_obs+aerosol+rayleigh) [m-1 sr-1] 
 Real*8, Intent(In)    :: mpl_back_aerosol(nz)  ! lidar aerosol backscatter [m-1 sr-1] 
 Real*8, Intent(In)    :: mpl_back_rayleigh(nz) ! lidar molecular backscatter [m-1 sr-1] 
 Real*8, Intent(In)    :: first_cloud_base      ! ceilometer first cloud base height [m] 
 Integer*2, Intent(OuT)  :: cloud_detection_flag(nz) ! cloud detection flag
 Integer*2, Intent(OuT)  :: cloud_mask(nz) ! cloud detection flag
 Real*8, Intent(OuT)   :: cloud_layer_base(n_layers) ! cloud layer base height [m]
 Real*8, Intent(OuT)   :: cloud_layer_top(n_layers) ! cloud layer top height [m]
 !
 Real*8,Parameter        :: m999=-999.d0
 Integer                 :: iz
 Integer                 :: k,cloud_layer_flag
 real*8,dimension(:),Allocatable   ::  radar_zhh_obs  ! dBZ
 real*8,dimension(:),Allocatable   ::  mpl_back_hydro_obs  ! [m-1 sr-1]  
 Real*8                  :: radar_zhh_atten2,atten_tot !2-way attenuation [dB]
 Real*8                  :: dz !dz [m]
 !
   Allocate(radar_zhh_obs(nz),mpl_back_hydro_obs(nz))
   !=======
   cloud_layer_flag = 0
   k=1
   cloud_layer_base(:)=m999
   cloud_layer_top(:)=m999
   !==MPL observed hydrometeor backscatter and
   !==radar observed reflectivity attenuated by precipitation
   != Two way precipitation attenuation at the first grid box
   if (radar_zhh_atten(1) > 0.d0) then
     atten_tot = radar_zhh_atten(1) * (zz(2)-zz(1))*0.001d0 * 2.d0 ! dB
   else
     atten_tot = 0.d0
   endif
   !
   do iz=1,nz
     mpl_back_hydro_obs(iz) = mpl_back_obs(iz)
     if (mpl_back_hydro(iz) <= 0.d0) then
      mpl_back_hydro_obs(iz) = 0.d0
     endif
     if (mpl_back_obs(iz) <= (mpl_back_rayleigh(iz)+mpl_back_aerosol(iz))) then
        mpl_back_hydro_obs(iz) = 0.d0
     endif
     if (iz >= 2) then
       dz = zz(iz)-zz(iz-1)
       if (radar_zhh_atten(iz) > 0.d0) then
         != Two way radar reflectivity attenuation at a grid box
         radar_zhh_atten2 = radar_zhh_atten(iz) * dz*0.001d0 * 2.d0 ! dB
         != Two way total radar reflectivity attenuation 
         atten_tot =  atten_tot + radar_zhh_atten2 ! dB
       endif
     endif
     != Observed radar reflectivity attenuated by precipitation
     radar_zhh_obs(iz) = radar_zhh_true(iz) - atten_tot
     if (radar_zhh_obs(iz) < radar_zhh_min(iz)) then
       radar_zhh_obs(iz)=m999
     endif
     != Cloud detection flag 1:Clear, 2:KAZR&MPL, 3: KAZR, 4:MPL
     cloud_detection_flag(iz) = int2(1)
     != Cloud mask 0: Clear grid box, 1:Cloudy grid box
     cloud_mask(iz) = int2(0)
     if (zz(iz) < first_cloud_base) then
       continue
     else
       if (mpl_back_hydro_obs(iz) >0.d0) then
         cloud_detection_flag(iz) = int2(4)
         if (radar_zhh_obs(iz) > m999) then
           cloud_detection_flag(iz) = int2(2)
         endif
       elseif (radar_zhh_obs(iz) > m999) then
         cloud_detection_flag(iz) = int2(3)
       endif
       if (cloud_detection_flag(iz) > int2(1)) then
         cloud_mask(iz) = int2(1)
       endif
     endif
     != Cloud layers
     if (k<=n_layers) then
       if (cloud_mask(iz)==int2(1)) then
         if (cloud_layer_flag == 0) then
           cloud_layer_base(k) = zz(iz) !cloud layer base heights           
         endif
         cloud_layer_flag = 1
       endif
       if (cloud_mask(iz)==int2(0)) then
         if (cloud_layer_flag == 1) then
           cloud_layer_top(k) = zz(iz)   !cloud layer top heights 
            k=k+1
         endif
         cloud_layer_flag = 0
       endif
     endif
   end do
   !=======

   Deallocate(radar_zhh_obs,mpl_back_hydro_obs)
   !
 return
 end subroutine estimate_arscl_cloudmask
 !
 !-----------------------------------
 ! Estimate microwave radiometer (MWR) liquid water path (LWP)
 ! taking account of field of view
 subroutine estimate_mwr_lwp(nx,ny,nz,zz,dx,dy,Theta,mwr_alt,&
                             lwc,true_lwp,mwr_lwp,n_samples)
 Use phys_param_mod, ONLY: pi
 Implicit None
 !
 integer, Intent(in)   :: nx,ny,nz
 Real*8, Intent(In)    :: zz(nx,ny,nz) ! height [m]
 Real*8, Intent(In)    :: dx,dy
 Real*8, Intent(In)    :: Theta ! Field of view of MWR [deg]
 Real*8, Intent(In)    :: mwr_alt ! height of MWR [m]
 Real*8, Intent(In)    :: lwc(nx,ny,nz) ! Liquid water content (rain+cloud) [kg m-3] 
 Real*8, Intent(Out)    :: mwr_lwp(nx,ny) ! MWR LWP [kg m-2]
 Real*8, Intent(Out)    :: true_lwp(nx,ny) ! original LWP [kg m-2]
 Integer, Intent(Out)   :: n_samples(nx,ny,nz) ! Number of samples within field of view
 
 Integer:: ix,iy,iz
 Integer:: iix,iiy,ix1,ix2,iy1,iy2
 Real*8 :: Theta_rad,dthetax,dthetay ! radian
 Real*8 :: dz
 Real*8 :: alt, hwhm, dis
 Real*8 :: wfac,wx,wy ! weight according to a Gaussian distribution
 Real*8 :: Wlwc,Wtot
 Real*8,Dimension(:,:,:),Allocatable :: mlwc ! Liquid water content (rain+cloud) [kg m-3]
   !
   Allocate(mlwc(nx,ny,nz))
   mlwc = lwc
   where (lwc <= 0.d0) mlwc=0.d0
  
   Theta_rad = Theta * pi /180.d0
   wfac=pi*pi/(-2.d0*dlog(2.d0))
   !
   do ix=1,nx
     do iy=1,ny
       true_lwp(ix,iy) = 0.d0
       do iz=1,nz
         n_samples(ix,iy,iz)=0
           if (iz==1) then
             dz = zz(ix,iy,2)-zz(ix,iy,1)
           else
             dz =  zz(ix,iy,iz)-zz(ix,iy,iz-1)
           endif
           ! original LWP
           true_lwp(ix,iy) = true_lwp(ix,iy) + mlwc(ix,iy,iz)*dz
           ! height from MWR
           alt = zz(ix,iy,iz) - mwr_alt 
           ! half width at half maximum (HWHM), which is a function of height
           hwhm = alt * dtan(Theta_rad/2.d0); ! in m
           ! extract horizontal domain covered by the MWR field of view
           ix1 = ix - nint(hwhm / dx) -1
           if (ix1 < 1) ix1=0
           ix2 = ix + nint(hwhm / dx) +1
           if (ix2 > nx) ix2=nx
           iy1 = iy - nint(hwhm / dy) -1
           if (iy1 < 1) iy1=0
           iy2 = iy + nint(hwhm / dy) +1
           if (iy2 > ny) iy2=ny
           Wlwc = 0.d0
           Wtot = 0.d0
           do iix=ix1,ix2
             do iiy=iy1,iy2
               dis = dsqrt((dble(ix-iix)*dx)**2.d0+(dble(iy-iiy)*dy)**2.d0)
               ! if (dis > hwhm) then
               if (datan(dis/alt) > (0.5d0*Theta_rad)) then
                 continue
               else
                 ! Counting number of samples within field of view 
                 n_samples(ix,iy,iz)=n_samples(ix,iy,iz)+1
                 ! angle difference in x and y
                 dthetax = datan(dabs(dble(iix-ix))*dx/alt)
                 dthetay = datan(dabs(dble(iiy-iy))*dy/alt)
                 ! weight depending on angle
                 wx = dexp(wfac * dthetax*dthetax/Theta_rad/Theta_rad)
                 wy = dexp(wfac * dthetay*dthetay/Theta_rad/Theta_rad)
                 Wlwc = Wlwc + (wx+wy)*mlwc(iix,iiy,iz)
                 Wtot = Wtot + (wx+wy)
               endif
             enddo
           enddo
           mwr_lwp(ix,iy) = mwr_lwp(ix,iy) + (Wlwc/Wtot)*dz
         enddo !iz
       enddo
     enddo
   Deallocate(mlwc)
   
 return
 end subroutine estimate_mwr_lwp
!!-----------------------------------                                                  
