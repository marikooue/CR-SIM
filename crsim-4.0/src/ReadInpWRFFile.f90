  !! ----------------------------------------------------------------------------
  !! ----------------------------------------------------------------------------
  !!  *PROGRAM* ReadInpWRFFile
  !!  @version crsim 4.x
  !!
  !!  *SRC_FILE*
  !!  crsim/src/ReadInpWRFFile.f90
  !!
  !!
  !!  *LAST CHANGES*
  !!  
  !!   Sep 15 2015  : A.T   u,v horizontal wind components are included in the
  !!                        input.They are needed for computation of radial
  !!                        Doppler velocity. Consequently,
  !!                        dimension nxp1 and nyp1 are read. It is assumed that
  !!                        nxp1=nx+1 and nyp1=ny+1.
  !!   Jan 13 2015  : A.T   xlog,xlat, and tke are included in the required WRF input. 
  !!                        tke is needed for the spectrum width computation.
  !!
  !!   Sep 19 2016  : M.O   Incorporateed Thompson microphysics.
  !!   Jul 17 2017  : M.O   Incorporateed ICON 2-moment microphysics.
  !!
  !!   Oct 30 2016  : D.W.  Incorporate P3 microphysics, adding ReadInpWRF_var_p3 (added by DW)
  !!   Jul 09 2019  : M.O.  Modified ReadInpWRF_MP_PHYSICS_50 to read QNCLOUD
  !!   Aug 20 2019  : A.T.  Corrected real-to-double conversion
  !!   Aug 26 2019  : A.T.  Introduced new module wrf_rvar_mod to read WRF
  !!                        variables in real format which are after reading
  !!                        transformed to double precission variables via wrf_var_mod. 
  !!   Aug 26 2019  : A.T.  Introduced indentation.
  !!   Mar    2024  : M.O.  Increased the number of bins for SBM (MP20)
  !!
  !!  *DESCRIPTION* 
  !!
  !!  This program  contains the subroutines needed for reading the input WRF  
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
  !!  You should have received a copy of the GNU General Public License
  !!  along with this program.  If not, see <http://www.gnu.org/licenses/>.
  !!
  !!----------------------------------------------------------------------------------------------------------------- 
  !!-----------------------------------------------------------------------------------------------------------------
  !!
  subroutine ReadInpWRF_dim(InpFile,str,status)
  Use netcdf
  Use wrf_var_mod
  Implicit None
  !
  Character(len=*),Intent(in)                              :: InpFile
  Type(wrf_var),Intent(InOut)                              :: str
  Integer, Intent(out)                                     :: status
  !
  Integer                                                  :: ncid,iDim,length
  Integer                                                  :: nDims,nVars
  Character(len=nf90_max_name)                             :: name, err_msg
  Integer, Allocatable, Dimension(:)                       :: dim_lengths
  Character(len=nf90_max_name), Allocatable , Dimension(:) :: dim_names
  Integer                                                  :: aa
    !
    !!
    !! Init
    str%nt=0
    str%nx=0
    str%ny=0
    str%nz=0
    str%nzp1=0
    !
    str%nxp1=0
    str%nyp1=0
    !
    status=  nf90_open(Trim(InpFile), NF90_NOWRITE, ncid)
    if (status.ne.0) then ; err_msg='Error in nf90_open' ; goto 999 ; endif
    !
    status= nf90_inquire(ncid, nDims, nVars)
    If (status/=0) Then ; err_msg = 'Error in nf90_inquire' ; goto 999 ; endif
    !
    Allocate(dim_names(1:nDims))
    Allocate(dim_lengths(1:nDims))
    !
    ! Get info on the dimensions
    ! Loop over the number of dimensions
    !   
    Do iDim = 1,nDims
      !
      status=nf90_inquire_dimension(ncid, idim, name=name, len=length)
      !If (status/=0) Then ; err_msg = 'Error in nf90_inquire_dimension' ; Goto 999 ; End If
      !
      dim_names(iDim)   = name
      dim_lengths(iDim) = length
      !
      Select Case (name)
      Case ('Time')             ;  str%nt    = length
      Case ('west_east')        ;  str%nx    = length
      Case ('west_east_stag')   ;  str%nxp1  = length
      Case ('south_north')      ;  str%ny    = length
      Case ('south_north_stag') ;  str%nyp1  = length
      Case ('bottom_top')       ;  str%nz    = length
      Case ('bottom_top_stag')  ;  str%nzp1  = length
      !
      Case ('DateStrLen')         ;  aa              = length
      Case ('soil_layers_stag')   ;  aa              = length
      !
      !Case default ; status = 1 ; err_msg = 'Unrecognised dimension: ' ; Goto 999
      End Select
      !
    EndDo   ! iDim
    !
    Deallocate(dim_names,dim_lengths)
    !
    status=  nf90_close(ncid)
    if (status.ne.0) then ; err_msg='Error in nf90_close' ; goto 999 ; endif
    !
999 If (status.ne.0) Then
      write(*,*) err_msg
      Call Exit(1)
    Endif
    !
    If (str%nxp1/=str%nx+1) then
      write(*,*) 'Problem in input dimensions nx and nxp1'
      Call Exit(1)
    EndIf
    !
    If (str%nyp1/=str%ny+1) then
      write(*,*) 'Problem in input dimensions ny and nyp1'
      Call Exit(1)
    EndIf
    !
  return
  end subroutine ReadInpWRF_dim
  !
  !!
  subroutine ReadInpWRF_var(InpFile,str,status)
  Use netcdf
  Use wrf_var_mod
  Use wrf_rvar_mod
  Implicit None
  !
  Character(len=*),Intent(In)     :: InpFile
  Type(wrf_var),Intent(InOut)     :: str
  Integer,Intent(out)             :: status
  !!
  Type(wrf_rvar)                                           :: strr ! input real*4 vars
  integer                                                  :: nDims, nVars
  integer                                                  :: ncid, iVar
  Character(len=nf90_max_name), Allocatable , Dimension(:) :: var_names
  Character(len=nf90_max_name)                             :: name, err_msg
  Real*8 :: dx,dy !added by oue
  Real*4, Allocatable , Dimension(:,:,:) :: landmask
    !
    !
    ! define strr dimensions
    strr%nt=str%nt
    strr%nx=str%nx
    strr%ny=str%ny
    strr%nz=str%nz
    strr%nzp1=str%nzp1
    strr%nxp1=str%nxp1
    strr%nyp1=str%nyp1
    !
    ! allocate strr vars
    call allocate_wrf_rvar(strr)
    call initialize_wrf_rvar(strr)
    !
    Allocate(landmask(strr%nx,strr%ny,strr%nt))
    !
    status=  nf90_open(Trim(InpFile), NF90_NOWRITE, ncid)
    If (status/=0) Then
      err_msg = 'Error in my_nf90_open'
      Goto 999
    End If
    !
    ! Get info on nDims and nVars
    status= nf90_inquire(ncid, nDims, nVars)
    If (status/=0) Then
      err_msg = 'Error in my_nf90_inquire'
      Goto 999
    End If
    !
    !
    ! Create arrays that contain the variables
    Allocate(var_names(1:nVars))
    !
    Do  iVar=1,nVars
      !
      status= nf90_inquire_variable(ncid, iVar, name=name)
      If (status/=0) Then
        err_msg = 'Error in my_nf90_inquire_variable: '//name
        Goto 999
      Endif
      !
      Select Case (Trim(Adjustl(name)))
      !
      Case ('RDX')
      status=nf90_get_var(ncid,iVar,strr%dx)
      str%dx=dble(strr%dx)
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: RDX'
        Goto 999
      End If
      !
      Case ('RDY')
      status=nf90_get_var(ncid,iVar,strr%dy)
      str%dy=dble(strr%dy)
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: RDY'
        Goto 999
      End If
      !
      Case ('U')
      status=nf90_get_var(ncid,iVar,strr%u)
      str%u=dble(strr%u)
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: U'
        Goto 999
      End If
      !
      Case ('V')
      status=nf90_get_var(ncid,iVar,strr%v)
      str%v=dble(strr%v)
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: V'
        Goto 999
      End If
      !
      Case ('W')
      status=nf90_get_var(ncid,iVar,strr%w)
      str%w=dble(strr%w)
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: W'
        Goto 999
      End If
      !
      Case ('PB')
      status=nf90_get_var(ncid,iVar,strr%pb)
      str%pb=dble(strr%pb)
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: PB'
        Goto 999
      End If
      !
      Case ('P')
      status=nf90_get_var(ncid,iVar,strr%p)
      str%p=dble(strr%p)
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: P'
        Goto 999
      End If
      !
      Case ('ALT')
      status=nf90_get_var(ncid,iVar,strr%rho_d)
      str%rho_d=dble(strr%rho_d)
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: ALT'
        Goto 999
      End If
      !
      Case ('T')
      status=nf90_get_var(ncid,iVar,strr%theta)
      str%theta=dble(strr%theta)
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: T'
        Goto 999
      End If
      !
      Case ('PHB')
      status=nf90_get_var(ncid,iVar,strr%phb)
      str%phb=dble(strr%phb)
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: PHB'
        Goto 999
      End If
      !
      Case ('PH')
      status=nf90_get_var(ncid,iVar,strr%ph)
      str%ph=dble(strr%ph)
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: PH'
        Goto 999
      End If
      !
      Case ('QVAPOR')
      status=nf90_get_var(ncid,iVar,strr%QVAPOR)
      str%QVAPOR=dble(strr%QVAPOR)
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: QVAPOR'
        Goto 999
      End If
      !
      Case ('XLAT')
      status=nf90_get_var(ncid,iVar,strr%xlat)
      str%xlat=dble(strr%xlat)
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: xlat'
        Goto 999
      End If
      !
      Case ('XLONG')
      status=nf90_get_var(ncid,iVar,strr%xlong)
      str%xlong=dble(strr%xlong)
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: xlong'
        Goto 999
      End If
      !
      Case ('TKE')
      status=nf90_get_var(ncid,iVar,strr%tke)
      str%tke=dble(strr%tke)
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: tke'
        Goto 999
      End If
      !
      Case ('TKE_PBL')    ! In P3, the TKE_PBL is the same as the TKE in other microphysics  (added by DW)
      status=nf90_get_var(ncid,iVar,strr%tke)         !(added by DW)
      str%tke=dble(strr%tke)
      If (status/=0) Then                            !(added by DW)
        err_msg = 'Error in my_nf90_get_var: tke_pbl'  !(added by DW)
        Goto 999                                       !(added by DW)
      End If                                         !(added by DW)
      !
      Case ('SSTSK')
      status=nf90_get_var(ncid,iVar,strr%sfctemp)
      str%sfctemp=dble(strr%sfctemp)
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: sfctemp'
        Goto 999
      End If
      !
      Case ('U10')
      status=nf90_get_var(ncid,iVar,strr%sfcuu)
      str%sfcuu=dble(strr%sfcuu)
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: sfuu (U10)'
        Goto 999
      End If
      !
      Case ('V10')
      status=nf90_get_var(ncid,iVar,strr%sfcvv)
      str%sfcvv=dble(strr%sfcvv)
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: sfcvv (V10)'
        Goto 999
      End If
      !
      Case ('XLAND')
      status=nf90_get_var(ncid,iVar,landmask)
      str%topo=dble(landmask(:,:,1))
      deallocate(landmask)
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: XLAND'
        Goto 999
      End If
      !
      End Select
      !
    End Do
    !
    !------------------------------------
    !if there is not RDX oy RDY, get from global attributes added by oue for ICON
    if (str%dx(1) <= 0.d0) then
      status=nf90_get_att(ncid,NF90_GLOBAL,"DX",dx)
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_att: DX'
        Goto 999
      endif
      str%dx=1.d0/dx
    endif
    if (str%dy(1) <= 0.d0) then
      status=nf90_get_att(ncid,NF90_GLOBAL,"DY",dy)
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_att: DY'
        Goto 999
      endif
      str%dy=1.d0/dy
    endif
    !------------------------------------
    ! 
    Deallocate(var_names)
    !
999 If (status.Ne.0) Then
      write(*,*) err_msg
      Call Exit(1)
    Else
      status=nf90_close(ncid)
    Endif
    !
    call deallocate_wrf_rvar(strr)
    !
  return
  end subroutine ReadInpWRF_var
  !!
  subroutine ReadInpWRF_MP_PHYSICS_09(InpFile,str,status)
  Use netcdf
  Use typeSizes
  Use wrf_rvar_mod
  Use wrf_var_mod
  Implicit None
  !
  Character(len=*),Intent(in)                              :: InpFile
  Type(wrf_var_mp09),Intent(InOut)                         :: str
  Integer,Intent(out)                                      :: status
  !
  Type(wrf_rvar_mp09)                                      :: strr
  integer                                                  :: nDims, nVars
  integer                                                  :: ncid, iVar
  Character(len=nf90_max_name), Allocatable , Dimension(:) :: var_names
  Character(len=nf90_max_name)                             :: name, err_msg
    !
    !
    ! define strr dimensions
    strr%nt=str%nt
    strr%nx=str%nx
    strr%ny=str%ny
    strr%nz=str%nz
    !
    ! allocate strr vars
    call allocate_wrf_rvar_mp09(strr)
    call initialize_wrf_rvar_mp09(strr)
    !
    !Open the file
    status=  nf90_open(Trim(InpFile), NF90_NOWRITE, ncid)
    If (status/=0) Then
      err_msg = 'Error in my_nf90_open'
      Goto 999
    End If
    !
    ! Get info on nDims and nVars
    status= nf90_inquire(ncid, nDims, nVars)
    If (status/=0) Then
      err_msg = 'Error in my_nf90_inquire'
      Goto 999
    End If
    !
    !
    ! Create arrays that contain the variables
    Allocate(var_names(1:nVars))
    !
    Do  iVar=1,nVars
      !
      status= nf90_inquire_variable(ncid, iVar, name=name)
      If (status/=0) Then
        err_msg = 'Error in my_nf90_inquire_variable: '//name
        Goto 999
      Endif
      !
      Select Case (Trim(Adjustl(name)))
      !
      Case ('QCLOUD')
      status=nf90_get_var(ncid,iVar,strr%QCLOUD)
      str%QCLOUD = dble(strr%QCLOUD)
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: QCLOUD'
        Goto 999
      End If
      !
      Case ('QRAIN')
      status=nf90_get_var(ncid,iVar,strr%QRAIN)
      str%QRAIN  = dble(strr%QRAIN)
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: QRAIN'
        Goto 999
      End If
      !
      Case ('QICE')
      status=nf90_get_var(ncid,iVar,strr%QICE)
      str%QICE  = dble(strr%QICE)
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: QICE'
        Goto 999
      End If
      !
      Case ('QSNOW')
      status=nf90_get_var(ncid,iVar,strr%QSNOW)
      str%QSNOW  = dble(strr%QSNOW)
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: QSNOW'
        Goto 999
      End If
      !
      Case ('QGRAUP')
      status=nf90_get_var(ncid,iVar,strr%QGRAUP)
      str%QGRAUP  = dble(strr%QGRAUP)
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: QGRAUP'
        Goto 999
      End If
      !
      Case ('QHAIL')
      status=nf90_get_var(ncid,iVar,strr%QHAIL)
      str%QHAIL  = dble(strr%QHAIL)
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: QHAIL'
        Goto 999
      End If
      !
      Case ('QNCLOUD')
      status=nf90_get_var(ncid,iVar,strr%QNCLOUD)
      str%QNCLOUD  = dble(strr%QNCLOUD)
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: QNCLOUD'
        Goto 999
      End If
      !
      Case ('QNRAIN')
      status=nf90_get_var(ncid,iVar,strr%QNRAIN)
      str%QNRAIN  = dble(strr%QNRAIN)
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: QNRAIN'
        Goto 999
      End If
      !
      Case ('QNICE')
      status=nf90_get_var(ncid,iVar,strr%QNICE)
      str%QNICE  = dble(strr%QNICE)
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: QNICE'
        Goto 999
      End If
      !
      Case ('QNSNOW')
      status=nf90_get_var(ncid,iVar,strr%QNSNOW)
      str%QNSNOW  = dble(strr%QNSNOW)
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: QNSNOW'
        Goto 999
      End If
      !
      Case ('QNGRAUPEL')
      status=nf90_get_var(ncid,iVar,strr%QNGRAUP)
      str%QNGRAUP  = dble(strr%QNGRAUP)
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: QNGRAUPEL'
        Goto 999
      End If
      !
      Case ('QNHAIL')
      status=nf90_get_var(ncid,iVar,strr%QNHAIL)
      str%QNHAIL  = dble(strr%QNHAIL)
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: QNHAIL'
        Goto 999
      End If
      !
      End Select
      !
    End Do
    !
    Deallocate(var_names)
    
999 If (status.Ne.0) Then
      write(*,*)  err_msg
      Call Exit(1)
    Else
      status=nf90_close(ncid)
    Endif
    !
    !
    go to 101
    write(*,*) '-----------------'
    write(*,*) '---mix ratio---'
    write(*,*) 'cloud',minval(str%qcloud),maxval(str%qcloud)
    write(*,*) 'rain',minval(str%qrain),maxval(str%qrain)
    write(*,*) 'ice',minval(str%qice),maxval(str%qice)
    write(*,*) 'snow',minval(str%qsnow),maxval(str%qsnow)
    write(*,*) 'graupel',minval(str%qgraup),maxval(str%qgraup)
    write(*,*) 'hail',minval(str%qhail),maxval(str%qhail)
    
    write(*,*) '---concentr---'
    write(*,*) 'cloud',minval(str%qncloud),maxval(str%qncloud)
    write(*,*) 'rain',minval(str%qnrain),maxval(str%qnrain)
    write(*,*) 'ice',minval(str%qnice),maxval(str%qnice)
    write(*,*) 'snow',minval(str%qnsnow),maxval(str%qnsnow)
    write(*,*) 'graupel',minval(str%qngraup)!,maxval(str%qngraup)
    write(*,*) 'hail',minval(str%qnhail),maxval(str%qnhail)
    ! 
    !write(*,*) '-----------------'
    !pause
    101 continue
    !
    call deallocate_wrf_rvar_mp09(strr)   
    ! 
  return
  end subroutine ReadInpWRF_MP_PHYSICS_09
  !!
  !!
  subroutine ReadInpWRF_MP_PHYSICS_10(InpFile,str,status)
  Use netcdf
  Use typeSizes
  Use wrf_rvar_mod
  Use wrf_var_mod
  Implicit None
  !
  Character(len=*),Intent(in)                              :: InpFile
  Type(wrf_var_mp10),Intent(InOut)                         :: str
  Integer,Intent(out)                                      :: status
  !
  Type(wrf_rvar_mp10)                                      :: strr  
  integer                                                  :: nDims, nVars
  integer                                                  :: ncid, iVar
  Character(len=nf90_max_name), Allocatable , Dimension(:) :: var_names
  Character(len=nf90_max_name)                             :: name, err_msg
    !
    !- Massages added by oue April 2018
    write(*,*) 'Selected microphysics: 10 (Morrison 2 moment)'
    write(*,*) 'Required hydrometeor parameters: QCLOUD, QRAIN, QICE, QSNOW, QGAUP '
    write(*,*) 'Required hydrometeor parameters: QNCLOUD(if any), QNRAIN, QNICE, QNSNOW, QNGAUPEL'
    write(*,*) 'Reading hydrometeor mixing ratios and number concentrations...'
    !
    ! define strr dimensions
    strr%nt=str%nt
    strr%nx=str%nx
    strr%ny=str%ny
    strr%nz=str%nz
    !
    ! allocate strr vars
    call allocate_wrf_rvar_mp10(strr)
    call initialize_wrf_rvar_mp10(strr)
    !
    !Open the file
    status=  nf90_open(Trim(InpFile), NF90_NOWRITE, ncid)
    If (status/=0) Then
      err_msg = 'Error in my_nf90_open'
      Goto 999
    End If
    !
    ! Get info on nDims and nVars
    status= nf90_inquire(ncid, nDims, nVars)
    If (status/=0) Then
      err_msg = 'Error in my_nf90_inquire'
      Goto 999
    End If
    !
    !
    ! Create arrays that contain the variables
    Allocate(var_names(1:nVars))
    !
    Do  iVar=1,nVars
      !
      status= nf90_inquire_variable(ncid, iVar, name=name)
      If (status/=0) Then
        err_msg = 'Error in my_nf90_inquire_variable: '//name
        Goto 999
      Endif
      !
      Select Case (Trim(Adjustl(name)))
      !
      Case ('QCLOUD')
      status=nf90_get_var(ncid,iVar,strr%QCLOUD) ; write(*,*) 'QCLOUD is found.'
      str%QCLOUD = dble(strr%QCLOUD)
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: QCLOUD'
        Goto 999
      End If
      !
      Case ('QRAIN')
      status=nf90_get_var(ncid,iVar,strr%QRAIN) ; write(*,*) 'QRAIN is found.'
      str%QRAIN = dble(strr%QRAIN)
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: QRAIN'
        Goto 999
      End If
      !
      Case ('QICE')
      status=nf90_get_var(ncid,iVar,strr%QICE) ; write(*,*) 'QICE is found.'
      str%QICE = dble(strr%QICE)
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: QICE'
        Goto 999
      End If
      !
      Case ('QSNOW')
      status=nf90_get_var(ncid,iVar,strr%QSNOW) ; write(*,*) 'QSNOW is found.'
      str%QSNOW = dble(strr%QSNOW)
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: QSNOW'
        Goto 999
      End If
      !
      Case ('QGRAUP')
      status=nf90_get_var(ncid,iVar,strr%QGRAUP) ; write(*,*) 'QGRAUP is found.'
      str%QGRAUP = dble(strr%QGRAUP)
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: QGRAUP'
        Goto 999
      End If
      !
      Case ('QNCLOUD')
      status=nf90_get_var(ncid,iVar,strr%QNCLOUD) ; write(*,*) 'QNCLOUD is found.'
      str%QNCLOUD = dble(strr%QNCLOUD)
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: QNCLOUD'
        Goto 999
      End If
      !
      Case ('QNRAIN')
      status=nf90_get_var(ncid,iVar,strr%QNRAIN) ; write(*,*) 'QNRAIN is found.'
      str%QNRAIN = dble(strr%QNRAIN)
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: QNRAIN'
        Goto 999
      End If
      !
      Case ('QNICE')
      status=nf90_get_var(ncid,iVar,strr%QNICE) ; write(*,*) 'QNICE is found.'
      str%QNICE = dble(strr%QNICE)
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: QNICE'
        Goto 999
      End If
      !
      Case ('QNSNOW')
      status=nf90_get_var(ncid,iVar,strr%QNSNOW) ; write(*,*) 'QNSNOW is found.'
      str%QNSNOW = dble(strr%QNSNOW)
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: QNSNOW'
        Goto 999
      End If
      !
      Case ('QNGRAUPEL')
      status=nf90_get_var(ncid,iVar,strr%QNGRAUP) ; write(*,*) 'QNGRAUPEL is found.'
      str%QNGRAUP = dble(strr%QNGRAUP)
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: QNGRAUPEL'
        Goto 999
      End If
      !
      End Select
      !
    End Do
    !
    Deallocate(var_names)
    
999 If (status.Ne.0) Then
      write(*,*)  err_msg
      Call Exit(1)
    Else
      status=nf90_close(ncid)
    Endif
    !
    call deallocate_wrf_rvar_mp10(strr)
    !
  return
  end subroutine ReadInpWRF_MP_PHYSICS_10
  !
  !-------------------------------------------------------------------------------------------------------------------
  !!! Added by oue, 2016/09/19
  subroutine ReadInpWRF_MP_PHYSICS_08(InpFile,str,status)
  Use netcdf
  Use typeSizes
  Use wrf_var_mod
  Use wrf_rvar_mod
  Implicit None
  !
  Character(len=*),Intent(in)                              :: InpFile
  Type(wrf_var_mp08),Intent(InOut)                         :: str
  Integer,Intent(out)                                      :: status
  !
  Type(wrf_rvar_mp08)                                      :: strr
  integer                                                  :: nDims, nVars
  integer                                                  :: ncid, iVar
  Character(len=nf90_max_name), Allocatable , Dimension(:) :: var_names
  Character(len=nf90_max_name)                             :: name, err_msg
    !
    !
    ! define strr dimensions
    strr%nt=str%nt
    strr%nx=str%nx
    strr%ny=str%ny
    strr%nz=str%nz
    !
    ! allocate strr vars
    call allocate_wrf_rvar_mp08(strr)
    call initialize_wrf_rvar_mp08(strr)
    !
    !Open the file
    status=  nf90_open(Trim(InpFile), NF90_NOWRITE, ncid)
    If (status/=0) Then
      err_msg = 'Error in my_nf90_open'
      Goto 999
    End If
    !
    ! Get info on nDims and nVars
    status= nf90_inquire(ncid, nDims, nVars)
    If (status/=0) Then
      err_msg = 'Error in my_nf90_inquire'
      Goto 999
    End If
    !
    !
    ! Create arrays that contain the variables
    Allocate(var_names(1:nVars))
    !
    Do  iVar=1,nVars
      !
      status= nf90_inquire_variable(ncid, iVar, name=name)
      If (status/=0) Then
        err_msg = 'Error in my_nf90_inquire_variable: '//name
        Goto 999
      Endif
      !
      Select Case (Trim(Adjustl(name)))
      !
      Case ('QCLOUD')
      status=nf90_get_var(ncid,iVar,strr%QCLOUD)
      str%QCLOUD = dble(strr%QCLOUD)
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: QCLOUD'
        Goto 999
      End If
      !
      Case ('QRAIN')
      status=nf90_get_var(ncid,iVar,strr%QRAIN)
      str%QRAIN = dble(strr%QRAIN)
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: QRAIN'
        Goto 999
      End If
      !
      Case ('QICE')
      status=nf90_get_var(ncid,iVar,strr%QICE)
      str%QICE = dble(strr%QICE)
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: QICE'
        Goto 999
      End If
      !
      Case ('QSNOW')
      status=nf90_get_var(ncid,iVar,strr%QSNOW)
      str%QSNOW = dble(strr%QSNOW)
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: QSNOW'
        Goto 999
      End If
      !
      Case ('QGRAUP')
      status=nf90_get_var(ncid,iVar,strr%QGRAUP)
      str%QGRAUP = dble(strr%QGRAUP)
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: QGRAUP'
        Goto 999
      End If
      !
      Case ('QNCLOUD')
      status=nf90_get_var(ncid,iVar,strr%QNCLOUD)
      str%QNCLOUD = dble(strr%QNCLOUD)
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: QNCLOUD'
        Goto 999
      End If
      !
      Case ('QNRAIN')
      status=nf90_get_var(ncid,iVar,strr%QNRAIN)
      str%QNRAIN = dble(strr%QNRAIN)
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: QNRAIN'
        Goto 999
      End If
      !
      Case ('QNICE')
      status=nf90_get_var(ncid,iVar,strr%QNICE)
      str%QNICE = dble(strr%QNICE)
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: QNICE'
        Goto 999
      End If
      !
      Case ('QNSNOW')
      status=nf90_get_var(ncid,iVar,strr%QNSNOW)
      str%QNSNOW = dble(strr%QNSNOW)
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: QNSNOW'
        Goto 999
      End If
      !
      Case ('QNGRAUPEL')
      status=nf90_get_var(ncid,iVar,strr%QNGRAUP)
      str%QNGRAUP = dble(strr%QNGRAUP)
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: QNGRAUPEL'
        Goto 999
      End If
      !
      End Select
      !
    End Do
    !
    Deallocate(var_names)
    
999 If (status.Ne.0) Then
      write(*,*)  err_msg
      Call Exit(1)
    Else
      status=nf90_close(ncid)
    Endif
    !
    call deallocate_wrf_rvar_mp08(strr)
    !
  return
  end subroutine ReadInpWRF_MP_PHYSICS_08
  !!-- added by oue
  !----------------------------------------------------------------
  !-------------------------------------------------------------------------------------------------------------------
  !!! Added by oue, 2017/07/17 for ICON
  subroutine ReadInpWRF_MP_PHYSICS_30(InpFile,str,status)
  Use netcdf
  Use typeSizes
  Use wrf_rvar_mod
  Use wrf_var_mod
  Implicit None
  !
  Character(len=*),Intent(in)                              :: InpFile
  Type(wrf_var_mp30),Intent(InOut)                         :: str
  Integer,Intent(out)                                      :: status
  !
  Type(wrf_rvar_mp30)                                      :: strr
  integer                                                  :: nDims, nVars
  integer                                                  :: ncid, iVar
  Character(len=nf90_max_name), Allocatable , Dimension(:) :: var_names
  Character(len=nf90_max_name)                             :: name, err_msg
    !
    !
    ! define strr dimensions
    strr%nt=str%nt
    strr%nx=str%nx
    strr%ny=str%ny
    strr%nz=str%nz
    !
    ! allocate strr vars
    call allocate_wrf_rvar_mp30(strr)
    call initialize_wrf_rvar_mp30(strr)
    !
    !Open the file
    status=  nf90_open(Trim(InpFile), NF90_NOWRITE, ncid)
    If (status/=0) Then
      err_msg = 'Error in my_nf90_open'
      Goto 999
    End If
    !
    ! Get info on nDims and nVars
    status= nf90_inquire(ncid, nDims, nVars)
    If (status/=0) Then
      err_msg = 'Error in my_nf90_inquire'
      Goto 999
    End If
    !
    !
    ! Create arrays that contain the variables
    Allocate(var_names(1:nVars))
    !
    Do  iVar=1,nVars
      !
      status= nf90_inquire_variable(ncid, iVar, name=name)
      If (status/=0) Then
        err_msg = 'Error in my_nf90_inquire_variable: '//name
        Goto 999
      Endif
      !
      Select Case (Trim(Adjustl(name)))
      !
      Case ('QCLOUD')
      status=nf90_get_var(ncid,iVar,strr%QCLOUD)
      str%QCLOUD = dble(strr%QCLOUD)
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: QCLOUD'
        Goto 999
      End If
      !
      Case ('QRAIN')
      status=nf90_get_var(ncid,iVar,strr%QRAIN)
      str%QRAIN = dble(strr%QRAIN)
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: QRAIN'
        Goto 999
      End If
      !
      Case ('QICE')
      status=nf90_get_var(ncid,iVar,strr%QICE)
      str%QICE = dble(strr%QICE)
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: QICE'
        Goto 999
      End If
      !
      Case ('QSNOW')
      status=nf90_get_var(ncid,iVar,strr%QSNOW)
      str%QSNOW = dble(strr%QSNOW)
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: QSNOW'
        Goto 999
      End If
      !
      Case ('QGRAUP')
      status=nf90_get_var(ncid,iVar,strr%QGRAUP)
      str%QGRAUP = dble(strr%QGRAUP)
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: QGRAUP'
        Goto 999
      End If
      !
      Case ('QHAIL')
      status=nf90_get_var(ncid,iVar,strr%QHAIL)
      str%QHAIL = dble(strr%QHAIL)
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: QHAIL'
        Goto 999
      End If
      !
      Case ('QNCLOUD')
      status=nf90_get_var(ncid,iVar,strr%QNCLOUD)
      str%QNCLOUD = dble(strr%QNCLOUD)
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: QNCLOUD'
        Goto 999
      End If
      !
      Case ('QNRAIN')
      status=nf90_get_var(ncid,iVar,strr%QNRAIN)
      str%QNRAIN = dble(strr%QNRAIN)
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: QNRAIN'
        Goto 999
      End If
      !
      Case ('QNICE')
      status=nf90_get_var(ncid,iVar,strr%QNICE)
      str%QNICE = dble(strr%QNICE)
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: QNICE'
        Goto 999
      End If
      !
      Case ('QNSNOW')
      status=nf90_get_var(ncid,iVar,strr%QNSNOW)
      str%QNSNOW = dble(strr%QNSNOW)
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: QNSNOW'
        Goto 999
      End If
      !
      Case ('QNGRAUPEL')
      status=nf90_get_var(ncid,iVar,strr%QNGRAUP)
      str%QNGRAUP = dble(strr%QNGRAUP)
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: QNGRAUPEL'
        Goto 999
      End If
      !
      Case ('QNGRAUP')
      status=nf90_get_var(ncid,iVar,strr%QNGRAUP)
      str%QNGRAUP = dble(strr%QNGRAUP)
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: QNGRAUPEL'
        Goto 999
      End If
      !
      Case ('QNHAIL')
      status=nf90_get_var(ncid,iVar,strr%QNHAIL)
      str%QNHAIL = dble(strr%QNHAIL)
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: QNHAIL'
        Goto 999
      End If
      !
      End Select
      !
    End Do
    !
    Deallocate(var_names)
    
999 If (status.Ne.0) Then
      write(*,*)  err_msg
      Call Exit(1)
    Else
      status=nf90_close(ncid)
    Endif
    !
    call deallocate_wrf_rvar_mp30(strr)
    !
  return
  end subroutine ReadInpWRF_MP_PHYSICS_30
  !!-- added by oue
  !----------------------------------------------------------------
  !----------------------------------------------------------------
  ! The whole subroutine ReadInpWRF_MP_PHYSICS_50(InpFile,str,status) is added by DW 2017/10/30 for P3
  
  subroutine ReadInpWRF_MP_PHYSICS_50(InpFile,str,status)
  Use netcdf
  Use typeSizes
  Use wrf_rvar_mod
  Use wrf_var_mod
  Implicit None
  !
  Character(len=*),Intent(in)                              :: InpFile
  Type(wrf_var_mp50),Intent(InOut)                         :: str
  Integer,Intent(out)                                      :: status
  !
  Type(wrf_rvar_mp50)                                       :: strr
  integer                                                  :: nDims, nVars
  integer                                                  :: ncid, iVar
  Character(len=nf90_max_name), Allocatable , Dimension(:) :: var_names
  Character(len=nf90_max_name)                             :: name, err_msg
    !
    ! define strr dimensions
    strr%nt=str%nt
    strr%nx=str%nx
    strr%ny=str%ny
    strr%nz=str%nz
    !
    ! allocate strr vars
    call allocate_wrf_rvar_mp50(strr)
    call initialize_wrf_rvar_mp50(strr)
    !
    !
    !Open the file
    status=  nf90_open(Trim(InpFile), NF90_NOWRITE, ncid)
    If (status/=0) Then
      err_msg = 'Error in my_nf90_open'
      Goto 999
    End If
    !
    ! Get info on nDims and nVars
    status= nf90_inquire(ncid, nDims, nVars)
    If (status/=0) Then
      err_msg = 'Error in my_nf90_inquire'
      Goto 999
    End If
    !
    !
    ! Create arrays that contain the variables
    Allocate(var_names(1:nVars))
    !
    Do  iVar=1,nVars
      !
      status= nf90_inquire_variable(ncid, iVar, name=name)
      If (status/=0) Then
        err_msg = 'Error in my_nf90_inquire_variable: '//name
        Goto 999
      Endif
      !
      Select Case (Trim(Adjustl(name)))
      !
      !-- For an output for MC3E May 20 case
      Case ('QICE1') ! read QICE1 instead of QICE
      status=nf90_get_var(ncid,iVar,strr%QICE)
      str%QICE = dble(strr%QICE)
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: QICE1'
        Goto 999
      End If
      !
      Case ('QNGRAUPEL1') ! Read QNGRAUPEL1 instead of QIR
      status=nf90_get_var(ncid,iVar,strr%QIR)
      str%QIR = dble(strr%QIR)
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: QNGRAUPEL1'
        Goto 999
      End If
      !
      Case ('QVGRAUPEL1') ! Read QVGRAUPEL1 instead of QIB
      status=nf90_get_var(ncid,iVar,strr%QIB)
      str%QIB = dble(strr%QIB)
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: QVGRAUPEL1'
        Goto 999
      End If
      !
      Case ('QNICE1') ! Read QNICE1 instead of QNICE
      status=nf90_get_var(ncid,iVar,strr%QNICE)
      str%QNICE = dble(strr%QNICE)
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: QNICE1'
        Goto 999
      End If
      !----------------------------------------
      !
      Case ('QCLOUD')
      status=nf90_get_var(ncid,iVar,strr%QCLOUD)
      str%QCLOUD = dble(strr%QCLOUD)
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: QCLOUD'
        Goto 999
      End If
      !
      Case ('QRAIN')
      status=nf90_get_var(ncid,iVar,strr%QRAIN)
      str%QRAIN = dble(strr%QRAIN)
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: QRAIN'
        Goto 999
      End If
      !
      Case ('QICE')
      status=nf90_get_var(ncid,iVar,strr%QICE)
      str%QICE = dble(strr%QICE)
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: QICE'
        Goto 999
      End If
      !
      ! read out QNCLOUD if any (available for P3 with unfixed cloud number concentration option), WRFv4.1, 20190709
      Case ('QNCLOUD')
      status=nf90_get_var(ncid,iVar,strr%QNCLOUD)
      str%QNCLOUD = dble(strr%QNCLOUD)
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: QNCLOUD'
        Goto 999
      End If
      !
      Case ('QNRAIN')
      status=nf90_get_var(ncid,iVar,strr%QNRAIN)
      str%QNRAIN = dble(strr%QNRAIN)
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: QNRAIN'
        Goto 999
      End If
      !
      Case ('QNICE')
      status=nf90_get_var(ncid,iVar,strr%QNICE)
      str%QNICE = dble(strr%QNICE)
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: QNICE'
        Goto 999
      End If
      !
      Case ('QIR')
      status=nf90_get_var(ncid,iVar,strr%QIR)
      str%QIR = dble(strr%QIR)
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: QIR'
        Goto 999
      End If
      !
      Case ('QIB')
      status=nf90_get_var(ncid,iVar,strr%QIB)
      str%QIB = dble(strr%QIB)
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: QIB'
        Goto 999
      End If
      !
      End Select
      !
    End Do
    !
    Deallocate(var_names)

999 If (status.Ne.0) Then
      write(*,*)  err_msg
      Call Exit(1)
    Else
      status=nf90_close(ncid)
    Endif
    !
    call deallocate_wrf_rvar_mp50(strr)
    !
  return
  end subroutine ReadInpWRF_MP_PHYSICS_50
  !
  !-----------------------------------------------------------
  !
  !
  subroutine ReadInpWRF_MP_PHYSICS_20(InpFile,str,status)
  Use netcdf
  Use typeSizes
  Use wrf_rvar_mod
  Use wrf_var_mod
  Implicit None
  !
  Character(len=*),Intent(in)                              :: InpFile
  Type(wrf_var_mp20),Intent(InOut)                         :: str
  Integer,Intent(out)                                      :: status
  !
  Type(wrf_rvar_mp20)                                      :: strr
  integer                                                  :: nDims, nVars,length
  integer                                                  :: ncid, iVar
  Character(len=nf90_max_name), Allocatable , Dimension(:) :: var_names
  Character(len=nf90_max_name)                             :: name, err_msg
    !
    !
    ! define strr dimensions
    strr%nt=str%nt
    strr%nx=str%nx
    strr%ny=str%ny
    strr%nz=str%nz
    strr%nbins=str%nbins

    !Open the file
    status=  nf90_open(Trim(InpFile), NF90_NOWRITE, ncid)
    If (status/=0) Then
      err_msg = 'Error in my_nf90_open'
      Goto 999
    End If
    !
    ! Get info on nDims and nVars
    status= nf90_inquire(ncid, nDims, nVars)
    If (status/=0) Then
      err_msg = 'Error in my_nf90_inquire'
      Goto 999
    End If
    !
    ! allocate strr vars
    call allocate_wrf_rvar_mp20(strr)
    call initialize_wrf_rvar_mp20(strr)
    !
     !
    ! Create arrays that contain the variables
    Allocate(var_names(1:nVars))
    !
    Do  iVar=1,nVars
      !
      status= nf90_inquire_variable(ncid, iVar, name=name)
      If (status/=0) Then
        err_msg = 'Error in my_nf90_inquire_variable: '//name
        Goto 999
      Endif
      !
      Select Case (Trim(Adjustl(name)))
         
      !
      Case ('ff1i01')
         write(*,*) name
      status=nf90_get_var(ncid,iVar,strr%ff1(:,:,:,:,1))
      str%ff1(:,:,:,:,1) = dble( strr%ff1(:,:,:,:,1) )
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: ff1i01'
        Goto 999
      End If
      !
      Case ('ff1i02')
      status=nf90_get_var(ncid,iVar,strr%ff1(:,:,:,:,2))
      str%ff1(:,:,:,:,2) = dble( strr%ff1(:,:,:,:,2) )
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: ff1i02'
        Goto 999
      End If
      !
      Case ('ff1i03')
      status=nf90_get_var(ncid,iVar,strr%ff1(:,:,:,:,3))
      str%ff1(:,:,:,:,3) = dble( strr%ff1(:,:,:,:,3) )
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: ff1i03'
        Goto 999
      End If
      !
      Case ('ff1i04')
      status=nf90_get_var(ncid,iVar,strr%ff1(:,:,:,:,4))
      str%ff1(:,:,:,:,4) = dble( strr%ff1(:,:,:,:,4) )
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: ff1i04'
        Goto 999
      End If
      !
      Case ('ff1i05')
      status=nf90_get_var(ncid,iVar,strr%ff1(:,:,:,:,5))
      str%ff1(:,:,:,:,5) = dble( strr%ff1(:,:,:,:,5) )
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: ff1i05'
        Goto 999
      End If
      !
      Case ('ff1i06')
      status=nf90_get_var(ncid,iVar,strr%ff1(:,:,:,:,6))
      str%ff1(:,:,:,:,6) = dble( strr%ff1(:,:,:,:,6) )
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: ff1i06'
        Goto 999
      End If
      !
      Case ('ff1i07')
      status=nf90_get_var(ncid,iVar,strr%ff1(:,:,:,:,7))
      str%ff1(:,:,:,:,7) = dble( strr%ff1(:,:,:,:,7) )
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: ff1i07'
        Goto 999
      End If
      !
      Case ('ff1i08')
      status=nf90_get_var(ncid,iVar,strr%ff1(:,:,:,:,8))
      str%ff1(:,:,:,:,8) = dble( strr%ff1(:,:,:,:,8) )
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: ff1i08'
        Goto 999
      End If
      !
      Case ('ff1i09')
      status=nf90_get_var(ncid,iVar,strr%ff1(:,:,:,:,9))
      str%ff1(:,:,:,:,9) = dble( strr%ff1(:,:,:,:,9) )
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: ff1i09'
        Goto 999
      End If
      !
      Case ('ff1i10')
      status=nf90_get_var(ncid,iVar,strr%ff1(:,:,:,:,10))
      str%ff1(:,:,:,:,10) = dble( strr%ff1(:,:,:,:,10) )
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: ff1i10'
        Goto 999
      End If
      !
      Case ('ff1i11')
      status=nf90_get_var(ncid,iVar,strr%ff1(:,:,:,:,11))
      str%ff1(:,:,:,:,11) = dble( strr%ff1(:,:,:,:,11) )
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: ff1i11'
        Goto 999
      End If
      !
      Case ('ff1i12')
      status=nf90_get_var(ncid,iVar,strr%ff1(:,:,:,:,12))
      str%ff1(:,:,:,:,12) = dble( strr%ff1(:,:,:,:,12) )
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: ff1i12'
        Goto 999
      End If
      !
      Case ('ff1i13')
      status=nf90_get_var(ncid,iVar,strr%ff1(:,:,:,:,13))
      str%ff1(:,:,:,:,13) = dble( strr%ff1(:,:,:,:,13) )
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: ff1i13'
        Goto 999
      End If
      !
      Case ('ff1i14')
      status=nf90_get_var(ncid,iVar,strr%ff1(:,:,:,:,14))
      str%ff1(:,:,:,:,14) = dble( strr%ff1(:,:,:,:,14) )
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: ff1i14'
        Goto 999
      End If
      !
      Case ('ff1i15')
      status=nf90_get_var(ncid,iVar,strr%ff1(:,:,:,:,15))
      str%ff1(:,:,:,:,15) = dble( strr%ff1(:,:,:,:,15) )
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: ff1i15'
        Goto 999
      End If
      !
      Case ('ff1i16')
      status=nf90_get_var(ncid,iVar,strr%ff1(:,:,:,:,16))
      str%ff1(:,:,:,:,16) = dble( strr%ff1(:,:,:,:,16) )
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: ff1i16'
        Goto 999
      End If
      !
      Case ('ff1i17')
      status=nf90_get_var(ncid,iVar,strr%ff1(:,:,:,:,17))
      str%ff1(:,:,:,:,17) = dble( strr%ff1(:,:,:,:,17) )
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: ff1i17'
        Goto 999
      End If
      !
      Case ('ff1i18')
      status=nf90_get_var(ncid,iVar,strr%ff1(:,:,:,:,18))
      str%ff1(:,:,:,:,18) = dble( strr%ff1(:,:,:,:,18) )
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: ff1i18'
        Goto 999
      End If
      !
      Case ('ff1i19')
      status=nf90_get_var(ncid,iVar,strr%ff1(:,:,:,:,19))
      str%ff1(:,:,:,:,19) = dble( strr%ff1(:,:,:,:,19) )
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: ff1i19'
        Goto 999
      End If
      !
      Case ('ff1i20')
      status=nf90_get_var(ncid,iVar,strr%ff1(:,:,:,:,20))
      str%ff1(:,:,:,:,20) = dble( strr%ff1(:,:,:,:,20) )
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: ff1i20'
        Goto 999
      End If
      !
      Case ('ff1i21')
      status=nf90_get_var(ncid,iVar,strr%ff1(:,:,:,:,21))
      str%ff1(:,:,:,:,21) = dble( strr%ff1(:,:,:,:,21) )
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: ff1i21'
        Goto 999
      End If
      !
      Case ('ff1i22')
      status=nf90_get_var(ncid,iVar,strr%ff1(:,:,:,:,22))
      str%ff1(:,:,:,:,22) = dble( strr%ff1(:,:,:,:,22) )
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: ff1i22'
        Goto 999
      End If
      !
      Case ('ff1i23')
      status=nf90_get_var(ncid,iVar,strr%ff1(:,:,:,:,23))
      str%ff1(:,:,:,:,23) = dble( strr%ff1(:,:,:,:,23) )
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: ff1i23'
        Goto 999
      End If
      !
      Case ('ff1i24')
      status=nf90_get_var(ncid,iVar,strr%ff1(:,:,:,:,24))
      str%ff1(:,:,:,:,24) = dble( strr%ff1(:,:,:,:,24) )
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: ff1i24'
        Goto 999
      End If
      !
      Case ('ff1i25')
      status=nf90_get_var(ncid,iVar,strr%ff1(:,:,:,:,25))
      str%ff1(:,:,:,:,25) = dble( strr%ff1(:,:,:,:,25) )
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: ff1i25'
        Goto 999
      End If
      !
      Case ('ff1i26')
      status=nf90_get_var(ncid,iVar,strr%ff1(:,:,:,:,26))
      str%ff1(:,:,:,:,26) = dble( strr%ff1(:,:,:,:,26) )
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: ff1i26'
        Goto 999
      End If
      !
      Case ('ff1i27')
      status=nf90_get_var(ncid,iVar,strr%ff1(:,:,:,:,27))
      str%ff1(:,:,:,:,27) = dble( strr%ff1(:,:,:,:,27) )
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: ff1i27'
        Goto 999
      End If
      !
      Case ('ff1i28')
      status=nf90_get_var(ncid,iVar,strr%ff1(:,:,:,:,28))
      str%ff1(:,:,:,:,28) = dble( strr%ff1(:,:,:,:,28) )
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: ff1i28'
        Goto 999
      End If
      !
      Case ('ff1i29')
      status=nf90_get_var(ncid,iVar,strr%ff1(:,:,:,:,29))
      str%ff1(:,:,:,:,29) = dble( strr%ff1(:,:,:,:,29) )
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: ff1i29'
        Goto 999
      End If
      !
      Case ('ff1i30')
      status=nf90_get_var(ncid,iVar,strr%ff1(:,:,:,:,30))
      str%ff1(:,:,:,:,30) = dble( strr%ff1(:,:,:,:,30) )
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: ff1i30'
        Goto 999
      End If
      !
      Case ('ff1i31')
      status=nf90_get_var(ncid,iVar,strr%ff1(:,:,:,:,31))
      str%ff1(:,:,:,:,31) = dble( strr%ff1(:,:,:,:,31) )
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: ff1i31'
        Goto 999
      End If
      !
      Case ('ff1i32')
      status=nf90_get_var(ncid,iVar,strr%ff1(:,:,:,:,32))
      str%ff1(:,:,:,:,32) = dble( strr%ff1(:,:,:,:,32) )
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: ff1i32'
        Goto 999
      End If
      !
      Case ('ff1i33')
      status=nf90_get_var(ncid,iVar,strr%ff1(:,:,:,:,33))
      str%ff1(:,:,:,:,33) = dble( strr%ff1(:,:,:,:,33) )
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: ff1i33'
        Goto 999
      End If
      !
      Case ('ff1i34')
      status=nf90_get_var(ncid,iVar,strr%ff1(:,:,:,:,34))
      str%ff1(:,:,:,:,34) = dble( strr%ff1(:,:,:,:,34) )
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: ff1i34'
        Goto 999
      End If
      !
      Case ('ff1i35')
      status=nf90_get_var(ncid,iVar,strr%ff1(:,:,:,:,35))
      str%ff1(:,:,:,:,35) = dble( strr%ff1(:,:,:,:,35) )
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: ff1i35'
        Goto 999
      End If
      !
      Case ('ff1i36')
      status=nf90_get_var(ncid,iVar,strr%ff1(:,:,:,:,36))
      str%ff1(:,:,:,:,36) = dble( strr%ff1(:,:,:,:,36) )
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: ff1i36'
        Goto 999
      End If
      !
      Case ('ff1i37')
      status=nf90_get_var(ncid,iVar,strr%ff1(:,:,:,:,37))
      str%ff1(:,:,:,:,37) = dble( strr%ff1(:,:,:,:,37) )
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: ff1i37'
        Goto 999
      End If
      !
      Case ('ff1i38')
      status=nf90_get_var(ncid,iVar,strr%ff1(:,:,:,:,38))
      str%ff1(:,:,:,:,38) = dble( strr%ff1(:,:,:,:,38) )
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: ff1i38'
        Goto 999
      End If
      !
      Case ('ff1i39')
      status=nf90_get_var(ncid,iVar,strr%ff1(:,:,:,:,39))
      str%ff1(:,:,:,:,39) = dble( strr%ff1(:,:,:,:,39) )
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: ff1i39'
        Goto 999
      End If
      !
      Case ('ff1i40')
      status=nf90_get_var(ncid,iVar,strr%ff1(:,:,:,:,40))
      str%ff1(:,:,:,:,40) = dble( strr%ff1(:,:,:,:,40) )
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: ff1i40'
        Goto 999
      End If
      !
      Case ('ff1i41')
      status=nf90_get_var(ncid,iVar,strr%ff1(:,:,:,:,41))
      str%ff1(:,:,:,:,41) = dble( strr%ff1(:,:,:,:,41) )
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: ff1i41'
        Goto 999
      End If
      !
      Case ('ff1i42')
      status=nf90_get_var(ncid,iVar,strr%ff1(:,:,:,:,42))
      str%ff1(:,:,:,:,42) = dble( strr%ff1(:,:,:,:,42) )
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: ff1i40'
        Goto 999
      End If
      !
      Case ('ff1i43')
      status=nf90_get_var(ncid,iVar,strr%ff1(:,:,:,:,43))
      str%ff1(:,:,:,:,43) = dble( strr%ff1(:,:,:,:,43) )
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: ff1i43'
        Goto 999
      End If
      !
      !-------------------------------------
      !
      Case ('ff5i01')
      status=nf90_get_var(ncid,iVar,strr%ff5(:,:,:,:,1))
      str%ff5(:,:,:,:,1) = dble( strr%ff5(:,:,:,:,1) )
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: ff5i01'
        Goto 999
      End If
      !
      Case ('ff5i02')
      status=nf90_get_var(ncid,iVar,strr%ff5(:,:,:,:,2))
      str%ff5(:,:,:,:,2) = dble( strr%ff5(:,:,:,:,2) )
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: ff5i02'
        Goto 999
      End If
      !
      Case ('ff5i03')
      status=nf90_get_var(ncid,iVar,strr%ff5(:,:,:,:,3))
      str%ff5(:,:,:,:,3) = dble( strr%ff5(:,:,:,:,3) )
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: ff5i03'
        Goto 999
      End If
      !
      Case ('ff5i04')
      status=nf90_get_var(ncid,iVar,strr%ff5(:,:,:,:,4))
      str%ff5(:,:,:,:,4) = dble( strr%ff5(:,:,:,:,4) )
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: ff5i04'
        Goto 999
      End If
      !
      Case ('ff5i05')
      status=nf90_get_var(ncid,iVar,strr%ff5(:,:,:,:,5))
      str%ff5(:,:,:,:,5) = dble( strr%ff5(:,:,:,:,5) )
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: ff5i05'
        Goto 999
      End If
      !
      Case ('ff5i06')
      status=nf90_get_var(ncid,iVar,strr%ff5(:,:,:,:,6))
      str%ff5(:,:,:,:,6) = dble( strr%ff5(:,:,:,:,6) )
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: ff5i06'
        Goto 999
      End If
      !
      Case ('ff5i07')
      status=nf90_get_var(ncid,iVar,strr%ff5(:,:,:,:,7))
      str%ff5(:,:,:,:,7) = dble( strr%ff5(:,:,:,:,7) )
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: ff5i07'
        Goto 999
      End If
      !
      Case ('ff5i08')
      status=nf90_get_var(ncid,iVar,strr%ff5(:,:,:,:,8))
      str%ff5(:,:,:,:,8) = dble( strr%ff5(:,:,:,:,8) )
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: ff5i08'
        Goto 999
      End If
      !
      Case ('ff5i09')
      status=nf90_get_var(ncid,iVar,strr%ff5(:,:,:,:,9))
      str%ff5(:,:,:,:,9) = dble( strr%ff5(:,:,:,:,9) )
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: ff5i09'
        Goto 999
      End If
      !
      Case ('ff5i10')
      status=nf90_get_var(ncid,iVar,strr%ff5(:,:,:,:,10))
      str%ff5(:,:,:,:,10) = dble( strr%ff5(:,:,:,:,10) )
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: ff5i10'
        Goto 999
      End If
      !
      Case ('ff5i11')
      status=nf90_get_var(ncid,iVar,strr%ff5(:,:,:,:,11))
      str%ff5(:,:,:,:,11) = dble( strr%ff5(:,:,:,:,11) )
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: ff5i11'
        Goto 999
      End If
      !
      Case ('ff5i12')
      status=nf90_get_var(ncid,iVar,strr%ff5(:,:,:,:,12))
      str%ff5(:,:,:,:,12) = dble( strr%ff5(:,:,:,:,12) )
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: ff5i12'
        Goto 999
      End If
      !
      Case ('ff5i13')
      status=nf90_get_var(ncid,iVar,strr%ff5(:,:,:,:,13))
      str%ff5(:,:,:,:,13) = dble( strr%ff5(:,:,:,:,13) )
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: ff5i13'
        Goto 999
      End If
      !
      Case ('ff5i14')
      status=nf90_get_var(ncid,iVar,strr%ff5(:,:,:,:,14))
      str%ff5(:,:,:,:,14) = dble( strr%ff5(:,:,:,:,14) )
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: ff5i14'
        Goto 999
      End If
      !
      Case ('ff5i15')
      status=nf90_get_var(ncid,iVar,strr%ff5(:,:,:,:,15))
      str%ff5(:,:,:,:,15) = dble( strr%ff5(:,:,:,:,15) )
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: ff5i15'
        Goto 999
      End If
      !
      Case ('ff5i16')
      status=nf90_get_var(ncid,iVar,strr%ff5(:,:,:,:,16))
      str%ff5(:,:,:,:,16) = dble( strr%ff5(:,:,:,:,16) )
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: ff5i16'
        Goto 999
      End If
      !
      !
      Case ('ff5i17')
      status=nf90_get_var(ncid,iVar,strr%ff5(:,:,:,:,17))
      str%ff5(:,:,:,:,17) = dble( strr%ff5(:,:,:,:,17) )
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: ff5i17'
        Goto 999
      End If
      !
      Case ('ff5i18')
      status=nf90_get_var(ncid,iVar,strr%ff5(:,:,:,:,18))
      str%ff5(:,:,:,:,18) = dble( strr%ff5(:,:,:,:,18) )
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: ff5i18'
        Goto 999
      End If
      !
      Case ('ff5i19')
      status=nf90_get_var(ncid,iVar,strr%ff5(:,:,:,:,19))
      str%ff5(:,:,:,:,19) = dble( strr%ff5(:,:,:,:,19) )
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: ff5i19'
        Goto 999
      End If
      !
      Case ('ff5i20')
      status=nf90_get_var(ncid,iVar,strr%ff5(:,:,:,:,20))
      str%ff5(:,:,:,:,20) = dble( strr%ff5(:,:,:,:,20) )
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: ff5i20'
        Goto 999
      End If
      !
      Case ('ff5i21')
      status=nf90_get_var(ncid,iVar,strr%ff5(:,:,:,:,21))
      str%ff5(:,:,:,:,21) = dble( strr%ff5(:,:,:,:,21) )
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: ff5i21'
        Goto 999
      End If
      !
      Case ('ff5i22')
      status=nf90_get_var(ncid,iVar,strr%ff5(:,:,:,:,22))
      str%ff5(:,:,:,:,22) = dble( strr%ff5(:,:,:,:,22) )
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: ff5i22'
        Goto 999
      End If
      !
      Case ('ff5i23')
      status=nf90_get_var(ncid,iVar,strr%ff5(:,:,:,:,23))
      str%ff5(:,:,:,:,23) = dble( strr%ff5(:,:,:,:,23) )
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: ff5i23'
        Goto 999
      End If
      !
      Case ('ff5i24')
      status=nf90_get_var(ncid,iVar,strr%ff5(:,:,:,:,24))
      str%ff5(:,:,:,:,24) = dble( strr%ff5(:,:,:,:,24) )
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: ff5i24'
        Goto 999
      End If
      !
      Case ('ff5i25')
      status=nf90_get_var(ncid,iVar,strr%ff5(:,:,:,:,25))
      str%ff5(:,:,:,:,25) = dble( strr%ff5(:,:,:,:,25) )
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: ff5i25'
        Goto 999
      End If
      !
      Case ('ff5i26')
      status=nf90_get_var(ncid,iVar,strr%ff5(:,:,:,:,26))
      str%ff5(:,:,:,:,26) = dble( strr%ff5(:,:,:,:,26) )
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: ff5i26'
        Goto 999
      End If
      !
      Case ('ff5i27')
      status=nf90_get_var(ncid,iVar,strr%ff5(:,:,:,:,27))
      str%ff5(:,:,:,:,27) = dble( strr%ff5(:,:,:,:,27) )
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: ff5i27'
        Goto 999
      End If
      !
      Case ('ff5i28')
      status=nf90_get_var(ncid,iVar,strr%ff5(:,:,:,:,28))
      str%ff5(:,:,:,:,28) = dble( strr%ff5(:,:,:,:,28) )
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: ff5i28'
        Goto 999
      End If
      !
      Case ('ff5i29')
      status=nf90_get_var(ncid,iVar,strr%ff5(:,:,:,:,29))
      str%ff5(:,:,:,:,29) = dble( strr%ff5(:,:,:,:,29) )
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: ff5i29'
        Goto 999
      End If
      !
      Case ('ff5i30')
      status=nf90_get_var(ncid,iVar,strr%ff5(:,:,:,:,30))
      str%ff5(:,:,:,:,30) = dble( strr%ff5(:,:,:,:,30) )
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: ff5i30'
        Goto 999
      End If
      !
      Case ('ff5i31')
      status=nf90_get_var(ncid,iVar,strr%ff5(:,:,:,:,31))
      str%ff5(:,:,:,:,31) = dble( strr%ff5(:,:,:,:,31) )
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: ff5i31'
        Goto 999
      End If
      !
      Case ('ff5i32')
      status=nf90_get_var(ncid,iVar,strr%ff5(:,:,:,:,32))
      str%ff5(:,:,:,:,32) = dble( strr%ff5(:,:,:,:,32) )
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: ff5i32'
        Goto 999
      End If
      !
      Case ('ff5i33')
      status=nf90_get_var(ncid,iVar,strr%ff5(:,:,:,:,33))
      str%ff5(:,:,:,:,33) = dble( strr%ff5(:,:,:,:,33) )
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: ff5i33'
        Goto 999
      End If
      !
       Case ('ff5i34')
      status=nf90_get_var(ncid,iVar,strr%ff5(:,:,:,:,34))
      str%ff5(:,:,:,:,34) = dble( strr%ff5(:,:,:,:,34) )
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: ff5i34'
        Goto 999
      End If
      !
      Case ('ff5i35')
      status=nf90_get_var(ncid,iVar,strr%ff5(:,:,:,:,35))
      str%ff5(:,:,:,:,35) = dble( strr%ff5(:,:,:,:,35) )
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: ff5i35'
        Goto 999
      End If
      !
      Case ('ff5i36')
      status=nf90_get_var(ncid,iVar,strr%ff5(:,:,:,:,36))
      str%ff5(:,:,:,:,36) = dble( strr%ff5(:,:,:,:,36) )
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: ff5i36'
        Goto 999
      End If
      !
      Case ('ff5i37')
      status=nf90_get_var(ncid,iVar,strr%ff5(:,:,:,:,37))
      str%ff5(:,:,:,:,37) = dble( strr%ff5(:,:,:,:,37) )
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: ff5i37'
        Goto 999
      End If
      !
      Case ('ff5i38')
      status=nf90_get_var(ncid,iVar,strr%ff5(:,:,:,:,38))
      str%ff5(:,:,:,:,38) = dble( strr%ff5(:,:,:,:,38) )
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: ff5i38'
        Goto 999
      End If
      !
      Case ('ff5i39')
      status=nf90_get_var(ncid,iVar,strr%ff5(:,:,:,:,39))
      str%ff5(:,:,:,:,39) = dble( strr%ff5(:,:,:,:,39) )
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: ff5i39'
        Goto 999
      End If
      !
      Case ('ff5i40')
      status=nf90_get_var(ncid,iVar,strr%ff5(:,:,:,:,40))
      str%ff5(:,:,:,:,40) = dble( strr%ff5(:,:,:,:,40) )
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: ff5i40'
        Goto 999
      End If
      !
      Case ('ff5i41')
      status=nf90_get_var(ncid,iVar,strr%ff5(:,:,:,:,41))
      str%ff5(:,:,:,:,41) = dble( strr%ff5(:,:,:,:,41) )
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: ff5i41'
        Goto 999
      End If
      !
      Case ('ff5i42')
      status=nf90_get_var(ncid,iVar,strr%ff5(:,:,:,:,42))
      str%ff5(:,:,:,:,42) = dble( strr%ff5(:,:,:,:,42) )
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: ff5i40'
        Goto 999
      End If
      !
      Case ('ff5i43')
      status=nf90_get_var(ncid,iVar,strr%ff5(:,:,:,:,43))
      str%ff5(:,:,:,:,43) = dble( strr%ff5(:,:,:,:,43) )
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: ff5i43'
        Goto 999
      End If
      !
      !-------------------------------------
      !
      Case ('ff6i01')
      status=nf90_get_var(ncid,iVar,strr%ff6(:,:,:,:,1))
      str%ff6(:,:,:,:,1) = dble( strr%ff6(:,:,:,:,1) )
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: ff6i01'
        Goto 999
      End If
      !
      Case ('ff6i02')
      status=nf90_get_var(ncid,iVar,strr%ff6(:,:,:,:,2))
      str%ff6(:,:,:,:,2) = dble( strr%ff6(:,:,:,:,2) )
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: ff6i02'
        Goto 999
      End If
      !
      Case ('ff6i03')
      status=nf90_get_var(ncid,iVar,strr%ff6(:,:,:,:,3))
      str%ff6(:,:,:,:,3) = dble( strr%ff6(:,:,:,:,3) )
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: ff6i03'
        Goto 999
      End If
      !
      Case ('ff6i04')
      status=nf90_get_var(ncid,iVar,strr%ff6(:,:,:,:,4))
      str%ff6(:,:,:,:,4) = dble( strr%ff6(:,:,:,:,4) )
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: ff6i04'
        Goto 999
      End If
      !
      Case ('ff6i05')
      status=nf90_get_var(ncid,iVar,strr%ff6(:,:,:,:,5))
      str%ff6(:,:,:,:,5) = dble( strr%ff6(:,:,:,:,5) )
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: ff6i05'
        Goto 999
      End If
      !
      Case ('ff6i06')
      status=nf90_get_var(ncid,iVar,strr%ff6(:,:,:,:,6))
      str%ff6(:,:,:,:,6) = dble( strr%ff6(:,:,:,:,6) )
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: ff6i06'
        Goto 999
      End If
      !
      Case ('ff6i07')
      status=nf90_get_var(ncid,iVar,strr%ff6(:,:,:,:,7))
      str%ff6(:,:,:,:,7) = dble( strr%ff6(:,:,:,:,7) )
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: ff6i07'
        Goto 999
      End If
      !
      Case ('ff6i08')
      status=nf90_get_var(ncid,iVar,strr%ff6(:,:,:,:,8))
      str%ff6(:,:,:,:,8) = dble( strr%ff6(:,:,:,:,8) )
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: ff6i08'
        Goto 999
      End If
      !
      Case ('ff6i09')
      status=nf90_get_var(ncid,iVar,strr%ff6(:,:,:,:,9))
      str%ff6(:,:,:,:,9) = dble( strr%ff6(:,:,:,:,9) )
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: ff6i09'
        Goto 999
      End If
      !
      Case ('ff6i10')
      status=nf90_get_var(ncid,iVar,strr%ff6(:,:,:,:,10))
      str%ff6(:,:,:,:,10) = dble( strr%ff6(:,:,:,:,10) )
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: ff6i10'
        Goto 999
      End If
      !
      Case ('ff6i11')
      status=nf90_get_var(ncid,iVar,strr%ff6(:,:,:,:,11))
      str%ff6(:,:,:,:,11) = dble( strr%ff6(:,:,:,:,11) )
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: ff6i11'
        Goto 999
      End If
      !
      Case ('ff6i12')
      status=nf90_get_var(ncid,iVar,strr%ff6(:,:,:,:,12))
      str%ff6(:,:,:,:,12) = dble( strr%ff6(:,:,:,:,12) )
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: ff6i12'
        Goto 999
      End If
      !
      Case ('ff6i13')
      status=nf90_get_var(ncid,iVar,strr%ff6(:,:,:,:,13))
      str%ff6(:,:,:,:,13) = dble( strr%ff6(:,:,:,:,13) )
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: ff6i13'
        Goto 999
      End If
      !
      Case ('ff6i14')
      status=nf90_get_var(ncid,iVar,strr%ff6(:,:,:,:,14))
      str%ff6(:,:,:,:,14) = dble( strr%ff6(:,:,:,:,14) )
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: ff6i14'
        Goto 999
      End If
      !
      Case ('ff6i15')
      status=nf90_get_var(ncid,iVar,strr%ff6(:,:,:,:,15))
      str%ff6(:,:,:,:,15) = dble( strr%ff6(:,:,:,:,15) )
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: ff6i15'
        Goto 999
      End If
      !
      Case ('ff6i16')
      status=nf90_get_var(ncid,iVar,strr%ff6(:,:,:,:,16))
      str%ff6(:,:,:,:,16) = dble( strr%ff6(:,:,:,:,16) )
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: ff6i16'
        Goto 999
      End If
      !
      Case ('ff6i17')
      status=nf90_get_var(ncid,iVar,strr%ff6(:,:,:,:,17))
      str%ff6(:,:,:,:,17) = dble( strr%ff6(:,:,:,:,17) )
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: ff6i17'
        Goto 999
      End If
      !
      Case ('ff6i18')
      status=nf90_get_var(ncid,iVar,strr%ff6(:,:,:,:,18))
      str%ff6(:,:,:,:,18) = dble( strr%ff6(:,:,:,:,18) )
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: ff6i18'
        Goto 999
      End If
      !
      Case ('ff6i19')
      status=nf90_get_var(ncid,iVar,strr%ff6(:,:,:,:,19))
      str%ff6(:,:,:,:,19) = dble( strr%ff6(:,:,:,:,19) )
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: ff6i19'
        Goto 999
      End If
      !
      Case ('ff6i20')
      status=nf90_get_var(ncid,iVar,strr%ff6(:,:,:,:,20))
      str%ff6(:,:,:,:,20) = dble( strr%ff6(:,:,:,:,20) )
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: ff6i20'
        Goto 999
      End If
      !
      Case ('ff6i21')
      status=nf90_get_var(ncid,iVar,strr%ff6(:,:,:,:,21))
      str%ff6(:,:,:,:,21) = dble( strr%ff6(:,:,:,:,21) )
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: ff6i21'
        Goto 999
      End If
      !
      Case ('ff6i22')
      status=nf90_get_var(ncid,iVar,strr%ff6(:,:,:,:,22))
      str%ff6(:,:,:,:,22) = dble( strr%ff6(:,:,:,:,22) )
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: ff6i22'
        Goto 999
      End If
      !
      Case ('ff6i23')
      status=nf90_get_var(ncid,iVar,strr%ff6(:,:,:,:,23))
      str%ff6(:,:,:,:,23) = dble( strr%ff6(:,:,:,:,23) )
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: ff6i23'
        Goto 999
      End If
      !
      Case ('ff6i24')
      status=nf90_get_var(ncid,iVar,strr%ff6(:,:,:,:,24))
      str%ff6(:,:,:,:,24) = dble( strr%ff6(:,:,:,:,24) )
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: ff6i24'
        Goto 999
      End If
      !
      Case ('ff6i25')
      status=nf90_get_var(ncid,iVar,strr%ff6(:,:,:,:,25))
      str%ff6(:,:,:,:,25) = dble( strr%ff6(:,:,:,:,25) )
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: ff6i25'
        Goto 999
      End If
      !
      Case ('ff6i26')
      status=nf90_get_var(ncid,iVar,strr%ff6(:,:,:,:,26))
      str%ff6(:,:,:,:,26) = dble( strr%ff6(:,:,:,:,26) )
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: ff6i26'
        Goto 999
      End If
      !
      Case ('ff6i27')
      status=nf90_get_var(ncid,iVar,strr%ff6(:,:,:,:,27))
      str%ff6(:,:,:,:,27) = dble( strr%ff6(:,:,:,:,27) )
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: ff6i27'
        Goto 999
      End If
      !
      Case ('ff6i28')
      status=nf90_get_var(ncid,iVar,strr%ff6(:,:,:,:,28))
      str%ff6(:,:,:,:,28) = dble( strr%ff6(:,:,:,:,28) )
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: ff6i28'
        Goto 999
      End If
      !
      Case ('ff6i29')
      status=nf90_get_var(ncid,iVar,strr%ff6(:,:,:,:,29))
      str%ff6(:,:,:,:,29) = dble( strr%ff6(:,:,:,:,29) )
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: ff6i29'
        Goto 999
      End If
      !
      Case ('ff6i30')
      status=nf90_get_var(ncid,iVar,strr%ff6(:,:,:,:,30))
      str%ff6(:,:,:,:,30) = dble( strr%ff6(:,:,:,:,30) )
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: ff6i30'
        Goto 999
      End If
      !
      Case ('ff6i31')
      status=nf90_get_var(ncid,iVar,strr%ff6(:,:,:,:,31))
      str%ff6(:,:,:,:,31) = dble( strr%ff6(:,:,:,:,31) )
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: ff6i31'
        Goto 999
      End If
      !
      Case ('ff6i32')
      status=nf90_get_var(ncid,iVar,strr%ff6(:,:,:,:,32))
      str%ff6(:,:,:,:,32) = dble( strr%ff6(:,:,:,:,32) )
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: ff6i32'
        Goto 999
      End If
      !
      Case ('ff6i33')
      status=nf90_get_var(ncid,iVar,strr%ff6(:,:,:,:,33))
      str%ff6(:,:,:,:,33) = dble( strr%ff6(:,:,:,:,33) )
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: ff6i33'
        Goto 999
      End If
      !
       Case ('ff6i34')
      status=nf90_get_var(ncid,iVar,strr%ff6(:,:,:,:,34))
      str%ff6(:,:,:,:,34) = dble( strr%ff6(:,:,:,:,34) )
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: ff6i34'
        Goto 999
      End If
      !
      Case ('ff6i35')
      status=nf90_get_var(ncid,iVar,strr%ff6(:,:,:,:,35))
      str%ff6(:,:,:,:,35) = dble( strr%ff6(:,:,:,:,35) )
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: ff6i35'
        Goto 999
      End If
      !
      Case ('ff6i36')
      status=nf90_get_var(ncid,iVar,strr%ff6(:,:,:,:,36))
      str%ff6(:,:,:,:,36) = dble( strr%ff6(:,:,:,:,36) )
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: ff6i36'
        Goto 999
      End If
      !
      Case ('ff6i37')
      status=nf90_get_var(ncid,iVar,strr%ff6(:,:,:,:,37))
      str%ff6(:,:,:,:,37) = dble( strr%ff6(:,:,:,:,37) )
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: ff6i37'
        Goto 999
      End If
      !
      Case ('ff6i38')
      status=nf90_get_var(ncid,iVar,strr%ff6(:,:,:,:,38))
      str%ff6(:,:,:,:,38) = dble( strr%ff6(:,:,:,:,38) )
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: ff6i38'
        Goto 999
      End If
      !
      Case ('ff6i39')
      status=nf90_get_var(ncid,iVar,strr%ff6(:,:,:,:,39))
      str%ff6(:,:,:,:,39) = dble( strr%ff6(:,:,:,:,39) )
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: ff6i39'
        Goto 999
      End If
      !
      Case ('ff6i40')
      status=nf90_get_var(ncid,iVar,strr%ff6(:,:,:,:,40))
      str%ff6(:,:,:,:,40) = dble( strr%ff6(:,:,:,:,40) )
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: ff6i40'
        Goto 999
      End If
      !
      Case ('ff6i41')
      status=nf90_get_var(ncid,iVar,strr%ff6(:,:,:,:,41))
      str%ff6(:,:,:,:,41) = dble( strr%ff6(:,:,:,:,41) )
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: ff6i41'
        Goto 999
      End If
      !
      Case ('ff6i42')
      status=nf90_get_var(ncid,iVar,strr%ff6(:,:,:,:,42))
      str%ff6(:,:,:,:,42) = dble( strr%ff6(:,:,:,:,42) )
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: ff6i40'
        Goto 999
      End If
      !
      Case ('ff6i43')
      status=nf90_get_var(ncid,iVar,strr%ff6(:,:,:,:,43))
      str%ff6(:,:,:,:,43) = dble( strr%ff6(:,:,:,:,43) )
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: ff6i43'
        Goto 999
      End If
      !     !-------------------------------------
      
      End Select
      !
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var:'
        Goto 999
      End If
      !
    End Do
    !
    Deallocate(var_names)
999 If (status.Ne.0) Then
      write(*,*) err_msg
      Call Exit(1)
    Else
      status=nf90_close(ncid)
    Endif
    !
    call deallocate_wrf_rvar_mp20(strr)
    !
  return
  end subroutine ReadInpWRF_MP_PHYSICS_20
  !
