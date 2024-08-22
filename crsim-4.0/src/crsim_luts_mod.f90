  !!
  !!  *PROGRAM* crsim_luts_mod
  !!  @version crsim 4.x
  !!
  !!  *SRC_FILE*
  !!  crsim/src/crsim_luts_mod.f90
  !!
  !!
  !!  *LAST CHANGES*
  !!
  !! JUL 17 2017            - M.O     Added graupel density (900) 
  !! Oct 30 2017            - D.W.    Incorporated P3 microphysics (added by DW)
  !! Mar 17 2018            - M.O     Added graupel density (500)
  !! Dec 09 2018            - A.T     Added snow_ar0.60 new densities (10 and 50) 
  !!                                  and n_den_l_snow modified from 5 to 7. Those
  !!                                  changes are directed by modifications in 
  !!                                  routines for Thompson microphysics made by M.O.
  !! Dec 10 2018            - A.T.    Introduced indentation.
  !!
  !!
  !! 
  !!  *DESCRIPTION* 
  !!
  !!  This module contains stored information regarding the LUTs
  !!
  !!
  !!  Copyright (C) 2024 Stony Brook University - Brookhaven National Laboratory
  !!  Radar Science Group (//radarscience.weebly.com//)
  !!  Contact email address: mariko.oue@stonybrook.edu pavlos.kollias@stonybrook.edu
  !!
  !!  This file is free software: you can redistribute it and/or modify
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
  module crsim_luts_mod
  Implicit None
  
  integer, parameter          ::  n_lfreq=6
  integer, parameter          ::  n_lelev=91
  
  integer, parameter          ::  n_ltemp_cld=23
  integer, parameter          ::  n_ltemp_rain=11
  integer, parameter          ::  n_ltemp_ice=1
  integer, parameter          ::  n_ltemp_snow=1
  integer, parameter          ::  n_ltemp_graup=1
  integer, parameter          ::  n_ltemp_smallice=1 ! added by DW 2017/10/30 for P3
  integer, parameter          ::  n_ltemp_unrice=1   ! added by DW 
  integer, parameter          ::  n_ltemp_graupP3=1  ! added by DW 
  integer, parameter          ::  n_ltemp_parice=1   ! added by DW 
  
  integer,parameter           :: n_lden_cld=1
  integer,parameter           :: n_lden_rain=1
  integer,parameter           :: n_lden_ice=6
  integer,parameter           :: n_lden_snow=7 !5 !7  !  5 modified by AT following modifications specifeed by MO
  integer,parameter           :: n_lden_graup=3 !2 !modified by oue 2017/07/17, 2018/03/13
  integer,parameter           :: n_lden_smallice=2  ! added by DW 2017/10/30 for P3
  integer,parameter           :: n_lden_unrice=26   ! added by DW 
  integer,parameter           :: n_lden_graupP3=29  ! added by DW 
  integer,parameter           :: n_lden_parice=29   ! added by DW 
  
  Real*8                      :: lfreq(n_lfreq)
  Real*8                      :: lelev(n_lelev)
  
  Real*8                      :: ltemp_cld(n_ltemp_cld)
  Real*8                      :: ltemp_rain(n_ltemp_rain)
  Real*8                      :: ltemp_ice(n_ltemp_ice)
  Real*8                      :: ltemp_snow(n_ltemp_snow)
  Real*8                      :: ltemp_graup(n_ltemp_graup)
  Real*8                      :: ltemp_smallice(n_ltemp_smallice) ! added by DW  2017/10/30 for P3
  Real*8                      :: ltemp_unrice(n_ltemp_unrice)     ! added by DW 
  Real*8                      :: ltemp_graupP3(n_ltemp_graupP3)   ! added by DW 
  Real*8                      :: ltemp_parice(n_ltemp_parice)     ! added by DW 
  
  Real*8                      :: lden_cld(n_lden_cld)
  Real*8                      :: lden_rain(n_lden_rain)
  Real*8                      :: lden_ice(n_lden_ice)
  Real*8                      :: lden_snow(n_lden_snow)
  Real*8                      :: lden_graup(n_lden_graup)
  Real*8                      :: lden_smallice(n_lden_smallice)  ! added by DW 2017/10/30 for P3
  Real*8                      :: lden_unrice(n_lden_unrice)      ! added by DW 
  Real*8                      :: lden_graupP3(n_lden_graupP3)    ! added by DW 
  Real*8                      :: lden_parice(n_lden_parice)      ! added by DW 
  

  data lfreq/3.d0, 5.5d0, 9.5d0, 13.6d0, 35.d0, 94.d0/
  
  data lelev/ 0.d0, 1.d0, 2.d0, 3.d0, 4.d0, 5.d0, 6.d0, 7.d0, 8.d0 ,9.d0,&
              10.d0, 11.d0, 12.d0, 13.d0, 14.d0, 15.d0, 16.d0, 17.d0, 18.d0 ,19.d0,&
              20.d0, 21.d0, 22.d0, 23.d0, 24.d0, 25.d0, 26.d0, 27.d0, 28.d0 ,29.d0,&
              30.d0, 31.d0, 32.d0, 33.d0, 34.d0, 35.d0, 36.d0, 37.d0, 38.d0 ,39.d0,&
              40.d0, 41.d0, 42.d0, 43.d0, 44.d0, 45.d0, 46.d0, 47.d0, 48.d0 ,49.d0,&
              50.d0, 51.d0, 52.d0, 53.d0, 54.d0, 55.d0, 56.d0, 57.d0, 58.d0 ,59.d0,&
              60.d0, 61.d0, 62.d0, 63.d0, 64.d0, 65.d0, 66.d0, 67.d0, 68.d0 ,69.d0,&
              70.d0, 71.d0, 72.d0, 73.d0, 74.d0, 75.d0, 76.d0, 77.d0, 78.d0 ,79.d0,&
              80.d0, 81.d0, 82.d0, 83.d0, 84.d0, 85.d0, 86.d0, 87.d0, 88.d0 ,89.d0, 90.d0/
  
  data ltemp_cld/244.16d0, 246.16d0, 248.16d0, 252.16d0,&
                 254.16d0, 256.16d0, 258.16d0, 260.16d0, 262.16d0, 264.16d0, 266.16d0,&
                 268.16d0, 270.16d0, 272.16d0, 274.16d0, 276.16d0, 278.16d0, 280.16d0, 282.16d0,&
                 284.16d0, 286.16d0, 288.16d0, 292.16d0/
  
  data ltemp_rain/272.16d0, 274.16d0, 276.16d0, 278.16d0,&
                  280.16d0, 282.16d0, 284.16d0, 286.16d0, 288.16d0, 290.16d0, 292.16d0/
  data ltemp_ice/243.16d0/
  data ltemp_snow/243.16d0/
  data ltemp_graup/243.16d0/
  data ltemp_smallice/243.16d0/  ! added by DW 2017/10/30 for P3
  data ltemp_unrice/243.16d0/    ! added by DW 
  data ltemp_graupP3/243.16d0/   ! added by DW 
  data ltemp_parice/243.16d0/    ! added by DW 
  
  data lden_cld/1000.d0/
  data lden_rain/1000.d0/
  data lden_ice/400.d0, 500.d0, 600.d0, 700.d0, 800.d0, 900.d0/
  data lden_snow/10.d0, 50.d0, 100.d0, 200.d0, 300.d0, 400.d0, 500.d0/ ! modified by AT following additions by MO
  !data lden_snow/100.d0, 200.d0, 300.d0, 400.d0, 500.d0/ 
 
  data lden_graup/400.d0, 500.d0, 900.d0/  !modified by oue 2017/07/17, 2018/03/13
  data lden_smallice/1.d0, 900.d0/   ! added by DW 2017/10/30 for P3
  data lden_unrice/1.d0, 5.d0, 10.d0, 20.d0, 30.d0, 40.d0, 50.d0, 60.d0, 70.d0, 80.d0,&
                   90.d0, 100.d0, 150.d0, 200.d0, 250.d0, 300.d0, 350.d0, 400.d0,&
                   450.d0, 500.d0, 550.d0, 600.d0, 650.d0, 700.d0, 750.d0, 800.d0/  ! added by DW 
  data lden_graupP3/1.d0, 2.d0, 3.d0, 4.d0, 5.d0, 6.d0, 7.d0, 8.d0, 9.d0, 10.d0, 20.d0,&
                   30.d0, 40.d0, 50.d0, 60.d0, 70.d0, 80.d0, 90.d0, 100.d0, 150.d0, 200.d0,&
                   250.d0, 300.d0, 350.d0, 400.d0, 450.d0, 500.d0, 550.d0, 600.d0/ ! added by DW 
  data lden_parice/1.d0, 2.d0, 3.d0, 4.d0, 5.d0, 6.d0, 7.d0, 8.d0, 9.d0, 10.d0, 20.d0,&
                   30.d0, 40.d0, 50.d0, 60.d0, 70.d0, 80.d0, 90.d0, 100.d0, 150.d0, 200.d0,&
                   250.d0, 300.d0, 350.d0, 400.d0, 450.d0, 500.d0, 550.d0, 600.d0/  ! added by DW 
  
  end module crsim_luts_mod
