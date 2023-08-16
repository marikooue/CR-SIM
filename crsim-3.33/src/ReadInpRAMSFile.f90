  !! ----------------------------------------------------------------------------
  !! ----------------------------------------------------------------------------
  !!  *PROGRAM* ReadInpWRFFile
  !!  @version crsim 3.x
  !!
  !!  *SRC_FILE*
  !!  crsim/src/ReadInpRAMSFile.f90
  !!
  !!
  !!  *LAST CHANGES*
  !!  
  !!   Jul 21 2017  : M.O   Read RAMS Input and calculate env variables
  !!
  !!   Apr    2018  : M.O   Read SAM warm bin (liquid: cloud/rain)
  !!   Jun 17 2018  : M.O   Read SAM Morrison 2-moment microphysics (MP_PHYSICS=75)
  !!   Aug 20 2019  : A.T   Corrected conversion issues double-to-real and similar
  !!   Aug 26 2019  : A.T.  Introduced new module wrf_rvar_mod to read WRF
  !!                        variables in real format which are after reading
  !!                        transformed to double precission variables via wrf_var_mod. 
  !!   Aug 26 2019  : A.T.  Introduced indentation.
  !!   Apr    2020  : M.O.  Read CM1 with Morrison 2-moment microphysics (MP_PHYSICS=80)
  !!
  !!
  !!  *DESCRIPTION* 
  !!
  !!  This program  contains the subroutines needed for reading the input WRF  
  !!  netcdf data files 
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
  subroutine ReadInpRAMS_dim(InpFile,str,status)
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
      Case ('phony_dim_0')        ;  str%nx    = length
      Case ('phony_dim_1')       ;  str%nz    = length
      
      Case ('phony_dim_2')       ;  aa              = length
      Case ('phony_dim_3')       ;  aa              = length
      Case ('phony_dim_4')       ;  aa              = length
      !
      !Case default ; status = 1 ; err_msg = 'Unrecognised dimension: ' ; Goto 999
      End Select
      !
    enddo   ! iDim
    !
    str%nt    = 1
    str%ny    = str%nx 
    str%nxp1  = str%nx
    str%nyp1  = str%ny
    str%nzp1  = str%nz
    
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
    If (str%nxp1/=str%nx) then
      write(*,*) 'Problem in input dimensions nx and nxp1'
      write(*,*) 'nx=nxp1 for RAMS input'
      Call Exit(1)
    EndIf
    !
    If (str%nyp1/=str%ny) then
      write(*,*) 'Problem in input dimensions ny and nyp1'
      write(*,*) 'nx=ny, ny=nyp1 for RAMS input'
      Call Exit(1)
    EndIf
    !
    !
  return
  end subroutine ReadInpRAMS_dim
  !
  !!
  subroutine ReadInpRAMS_var(InpFile,InpText,str,status)
  Use netcdf
  Use wrf_rvar_mod
  Use wrf_var_mod
  Use phys_param_mod, Only: Rearth
  Implicit None
  !
  Character(len=*),Intent(In)     :: InpFile
  Character(len=*),Intent(In)     :: InpText
  Type(wrf_var),Intent(InOut)     :: str
  Integer,Intent(out)             :: status
  !!
  Type(wrf_rvar)                  :: strr
  integer                                                  :: nDims, nVars
  integer                                                  :: ncid, iVar
  Character(len=nf90_max_name), Allocatable , Dimension(:) :: var_names
  Character(len=nf90_max_name)                             :: name, err_msg
  Real*8           ::  dlon,dlat,dlon_km,dlat_km,dx,dy
  Real*8           ::  pi
  Integer          ::  ixc,iyc
  LOGICAL          :: file_exists
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
    call allocate_wrf_rvar(strr)
    call initialize_wrf_rvar(strr)
    
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
      str%dx = dble(strr%dx)
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: RDX'
        Goto 999
      End If
      !
      Case ('RDY')
      status=nf90_get_var(ncid,iVar,strr%dy)
      str%dy = dble(strr%dy)
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: RDY'
        Goto 999
      End If
      !
      Case ('UC')
      status=nf90_get_var(ncid,iVar,strr%u)
      str%u = dble(strr%u)
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: UC'
        Goto 999
      End If
      !
      Case ('VC')
      status=nf90_get_var(ncid,iVar,strr%v)
      str%v = dble(strr%v)
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: VC'
        Goto 999
      End If
      !
      Case ('WC')
      status=nf90_get_var(ncid,iVar,strr%w)
      str%w = dble(strr%w)
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: WC'
        Goto 999
      End If
      !
      Case ('PI')
      status=nf90_get_var(ncid,iVar,strr%pb)
      str%pb = dble(strr%pb)
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: PI'
        Goto 999
      End If
      !
      Case ('PC')
      status=nf90_get_var(ncid,iVar,strr%p)
      str%p = dble(strr%p)
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: PC'
        Goto 999
      End If
      !
      !Case ('DN0') 
      !status=nf90_get_var(ncid,iVar,strr%rho_d)
      !str%rho_d = dble(strr%rho_d)
      !If (status/=0) Then
        !err_msg = 'Error in my_nf90_get_var: ALT'
        !Goto 999
      !End If
      !
      Case ('THETA')
      status=nf90_get_var(ncid,iVar,strr%theta)
      str%theta = dble(strr%theta)
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: T'
        Goto 999
      End If
      !
      Case ('TOPT')
      status=nf90_get_var(ncid,iVar,strr%topo)
      str%topo = dble(strr%topo)
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: TOPT'
        Goto 999
      End If
      !
      Case ('RV')
      status=nf90_get_var(ncid,iVar,strr%QVAPOR)
      str%QVAPOR = dble(strr%QVAPOR)
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: QVAPOR'
        Goto 999
      End If
      !
      Case ('GLAT')
      status=nf90_get_var(ncid,iVar,strr%xlat)
      str%xlat = dble(strr%xlat)
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: xlat'
        Goto 999
      End If
      !
      Case ('GLON')
      status=nf90_get_var(ncid,iVar,strr%xlong)
      str%xlong = dble(strr%xlong)
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: xlong'
        Goto 999
      End If
      !
      Case ('TKE')
      status=nf90_get_var(ncid,iVar,strr%tke)
      str%tke = dble(strr%tke)
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: tke'
        Goto 999
      End If
      !
      End Select
      !
    End Do
    !
    !------------------------------------
    !if there is not RDX or RDY, calculate from GLON/GLAT
    INQUIRE(FILE=Trim(InpText), EXIST=file_exists) 
    If (file_exists) Then
      write(*,*) 'Reading dx, z from ', Trim(InpText)
      call ReadInpRAMS_text(InpText,str%nz,dx,str%hgtm)
      str%dx=1.d0/dx
      str%dy=1.d0/dx !modified by oue Apr 2020; RAMS uses same dx=dy.
      write(*,*) 'RAMS dx,dy',dx
    endif
    
    if ((str%dx(1) <=0.d0) .or. (str%dy(1) <=0.d0)) then
      write(*,*) 'Calculating dx, dy'
      ixc=nint(str%nx/2.d0)
      iyc=nint(str%ny/2.d0)
      dlon=str%xlong(ixc+1,iyc,1)-str%xlong(ixc,iyc,1)
      dlat=str%xlat(ixc,iyc+1,1)-str%xlat(ixc,iyc,1)
      pi=dacos(-1.d0)
      dlat_km=0.008993d0
      dlon_km=360.d0 / (2.d0 * pi * Rearth * dcos(str%xlat(ixc,iyc,1) /180.d0 * pi))
      dx=dlon/dlon_km * 1000.d0 !m
      dy=dlat/dlat_km * 1000.d0 !m
      str%dx=1.d0/dx
      str%dy=1.d0/dy
      write(*,*) 'dx,dy',dx,dy
    endif
    
    !------------------------------------
    
    Deallocate(var_names)
    
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
  end subroutine ReadInpRAMS_var
  !!
  !----------------------------------------------------------------
  !-------------------------------------------------------------------------------------------------------------------
  !!! Added by oue, 2017/07/21 for RAMS
  subroutine ReadInpRAMS_MP_PHYSICS_40(InpFile,str,status)
  Use netcdf
  Use typeSizes
  Use wrf_var_mod
  Use wrf_rvar_mod
  Implicit None
  !
  Character(len=*),Intent(in)                              :: InpFile
  Type(wrf_var_mp40),Intent(InOut)                         :: str
  Integer,Intent(out)                                      :: status
  !
  Type(wrf_rvar_mp40)                                      :: strr
  integer                                                  :: nDims, nVars
  integer                                                  :: ncid, iVar
  Character(len=nf90_max_name), Allocatable , Dimension(:) :: var_names
  Character(len=nf90_max_name)                             :: name, err_msg
    
    ! define strr dimensions
    strr%nt=str%nt
    strr%nx=str%nx
    strr%ny=str%ny
    strr%nz=str%nz
    !
    call allocate_wrf_rvar_mp40(strr)
    call initialize_wrf_rvar_mp40(strr)    
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
      Case ('RCP')
      status=nf90_get_var(ncid,iVar,strr%QCLOUD)
      str%QCLOUD = dble(strr%QCLOUD)
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: QCLOUD'
        Goto 999
      End If
      !
      Case ('RDP')
      status=nf90_get_var(ncid,iVar,strr%QDRZL)
      str%QDRZL = dble(strr%QDRZL)
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: QDRIZZLE'
        Goto 999
      End If
      !
      Case ('RRP')
      status=nf90_get_var(ncid,iVar,strr%QRAIN)
      str%QRAIN = dble(strr%QRAIN)
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: QRAIN'
        Goto 999
      End If
      !
      Case ('RPP')
      status=nf90_get_var(ncid,iVar,strr%QICE)
      str%QICE = dble(strr%QICE)
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: QICE'
        Goto 999
      End If
      !
      Case ('RSP')
      status=nf90_get_var(ncid,iVar,strr%QSNOW)
      str%QSNOW = dble(strr%QSNOW)
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: QSNOW'
        Goto 999
      End If
      !
      Case ('RAP')
      status=nf90_get_var(ncid,iVar,strr%QAGGR)
      str%QAGGR = dble(strr%QAGGR)
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: QAGGREGATE'
        Goto 999
      End If
      !
      Case ('RGP')
      status=nf90_get_var(ncid,iVar,strr%QGRAUP)
      str%QGRAUP = dble(strr%QGRAUP)
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: QGRAUP'
        Goto 999
      End If
      !
      Case ('RHP')
      status=nf90_get_var(ncid,iVar,strr%QHAIL)
      str%QHAIL = dble(strr%QHAIL)
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: QHAIL'
        Goto 999
      End If
      !
      Case ('CCP')
      status=nf90_get_var(ncid,iVar,strr%QNCLOUD)
      str%QNCLOUD = dble(strr%QNCLOUD)
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: QNCLOUD'
        Goto 999
      End If
      !
      Case ('CDP')
      status=nf90_get_var(ncid,iVar,strr%QNDRZL)
      str%QNDRZL = dble(strr%QNDRZL)
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: QNDRIZZLE'
        Goto 999
      End If
      !
      Case ('CRP')
      status=nf90_get_var(ncid,iVar,strr%QNRAIN)
      str%QNRAIN = dble(strr%QNRAIN)
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: QNRAIN'
        Goto 999
      End If
      !
      Case ('CPP')
      status=nf90_get_var(ncid,iVar,strr%QNICE)
      str%QNICE = dble(strr%QNICE)
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: QNICE'
        Goto 999
      End If
      !
      Case ('CSP')
      status=nf90_get_var(ncid,iVar,strr%QNSNOW)
      str%QNSNOW = dble(strr%QNSNOW)
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: QNSNOW'
        Goto 999
      End If
      !
      Case ('CAP')
      status=nf90_get_var(ncid,iVar,strr%QNAGGR)
      str%QNAGGR = dble(strr%QNAGGR)
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: QNAGGREGATE'
        Goto 999
      End If
      !
      Case ('CGP')
      status=nf90_get_var(ncid,iVar,strr%QNGRAUP)
      str%QNGRAUP = dble(strr%QNGRAUP)
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: QNGRAUPEL'
        Goto 999
      End If
      !
      Case ('CHP')
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
    call deallocate_wrf_rvar_mp40(strr)
    !
  return
  end subroutine ReadInpRAMS_MP_PHYSICS_40
  !!-- added by oue
  !----------------------------------------------------------------
  !!
  subroutine get_env_vars_rams(conf,wrf,env)
  Use wrf_var_mod
  Use crsim_mod
  Use phys_param_mod, ONLY: Rd, eps,p0,cp,grav
  Implicit None
  !
  Type(conf_var),Intent(In)                  :: conf
  Type(wrf_var),Intent(InOut)                :: wrf
  Type(env_var),Intent(InOut)                :: env
  !
  Integer                                    :: ix1,ix2, iy1,iy2,iz1,iz2,it
  Real*8,Dimension(:,:,:,:), Allocatable     :: rho_d
  Real*8,Dimension(:),Allocatable            :: xtrack,ytrack
  Real*8,Dimension(:,:,:),Allocatable        :: gheight,ww
  Real*8,Dimension(:,:,:),Allocatable        :: uu,vv
  Real*8,Dimension(:,:,:),Allocatable        :: Ku,Kv,Kw
  Real*8                                     :: hres,kap
  !!
  Real*8             :: dp,dz
  
  Integer            :: ix,iy,iz,izmin,izmax,idz,iz_inv
  Real*8,Dimension(:),Allocatable            :: hgt
  
    !-------------------------------------------------------------
    !
    ix1=1 ; ix2=env%nx
    iy1=1 ; iy2=env%ny
    iz1=1 ; iz2=env%nz
    it = 1
    !           
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
    kap=cp/Rd
    do ix=1,wrf%nx
      do iy=1,wrf%ny
        do iz=1,wrf%nz
          wrf%press(ix,iy,iz,it)=p0 *((wrf%pb(ix,iy,iz,it)+wrf%p(ix,iy,iz,it))/cp) ** kap   ! Pa   
        enddo
      enddo
    enddo
    wrf%temp=wrf%theta*(wrf%press/100000.d0)**(287.d0/1004.d0) ! K  
    wrf%geop_height=(wrf%phb+wrf%ph)/9.81d0 ! m  

    izmin=1 ; izmax=wrf%nz ; idz=1
    if(wrf%press(1,1,1,1) < wrf%press(1,1,wrf%nz,1)) then
      izmin=wrf%nz ; izmax=1 ; idz=-1
    endif

    if(MaxVal(wrf%hgtm(:))<=0.d0) then
      !-Calculate z using hydrostatic equilibrium
      Allocate(rho_d(wrf%nx,wrf%ny,wrf%nz,wrf%nt))
      if (MaxVal(wrf%rho_d)<-888.d0) then
        rho_d=(wrf%press*wrf%qvapor)/(eps+wrf%qvapor) ! e water vapor pressure in Pa  
        rho_d=(wrf%press-rho_d)/(Rd * wrf%temp)   !  kg/m^3   
      else
        rho_d=wrf%rho_d
      endif
      do ix=1,wrf%nx
        do iy=1,wrf%ny
          dz=wrf%topo(ix,iy)
            do iz=izmin,izmax,idz
              if(iz==1) then
                dp=wrf%press(ix,iy,iz,1)-wrf%press(ix,iy,iz+1,1)
              else
                dp=wrf%press(ix,iy,iz-1,1)-wrf%press(ix,iy,iz,1)
              endif
              if (wrf%press(1,1,1,1) < wrf%press(1,1,wrf%nz,1)) dp=dp*(-1.d0)
              dz=dz+dp/grav/rho_d(ix,iy,iz,1)
              wrf%geop_height(ix,iy,iz,1)=dz
            enddo
          enddo
        enddo
        Deallocate(rho_d)
      else
      !- z from RAMS text file
      Allocate(hgt(wrf%nz))
      do iz=izmin,izmax,idz
        iz_inv=iz
        if(wrf%press(1,1,1,1) < wrf%press(1,1,wrf%nz,1)) iz_inv=wrf%nz-iz+1
        if(wrf%hgtm(iz_inv)<0.d0) wrf%hgtm(iz_inv)=0.d0
        hgt(iz)=wrf%hgtm(iz_inv)
      enddo
      do ix=1,wrf%nx
        do iy=1,wrf%ny
          wrf%geop_height(ix,iy,:,1)=wrf%topo(ix,iy)+hgt(:)
        enddo
      enddo
      Deallocate(hgt)
    endif
    !  
    !------------------------------------------------------------------------------------  
    !-------------------------------------------------------------------------------       
    !
    !saska
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
    hres=1.d0/wrf%dy(it)
    env%dy=hres
    do iy=1,wrf%ny
      ytrack(iy)=dble(iy-1)*hres
    enddo
    !!----------------
    do iz=1,wrf%nz
      gheight(:,:,iz)=wrf%geop_height(:,:,iz,it)
      ww(:,:,iz)=wrf%w(:,:,iz,it)
      if(iz<wrf%nz) then
        Kw(:,:,iz)= (wrf%w(:,:,iz+1,it)-wrf%w(:,:,iz,it)) / (wrf%geop_height(:,:,iz+1,it)-wrf%geop_height(:,:,iz,it))
      else
        Kw(:,:,iz)= (wrf%w(:,:,iz,it)-wrf%w(:,:,iz-1,it)) / (wrf%geop_height(:,:,iz,it)-wrf%geop_height(:,:,iz-1,it))
      endif
      if (wrf%geop_height(1,1,1,1) > wrf%geop_height(1,1,wrf%nz,1)) Kw(:,:,iz)=Kw(:,:,iz)*(-1.d0)
    enddo
    !
    !---------------- 
    do ix=1,wrf%nx
      uu(ix,:,:)=wrf%u(ix,:,:,it)
      if(ix<wrf%nx)then
        Ku(ix,:,:)=wrf%dx(it) * (wrf%u(ix+1,:,:,it)-wrf%u(ix,:,:,it))
      else
        Ku(ix,:,:)=wrf%dx(it) * (wrf%u(ix,:,:,it)-wrf%u(ix-1,:,:,it))
      endif
    enddo
    !
    do iy=1,wrf%ny
      vv(:,iy,:)=wrf%v(:,iy,:,it)
      if(iy<wrf%ny)then
        Kv(:,iy,:)=wrf%dy(it) * (wrf%v(:,iy+1,:,it)-wrf%v(:,iy,:,it))
      else
        Kv(:,iy,:)=wrf%dy(it) * (wrf%v(:,iy,:,it)-wrf%v(:,iy-1,:,it))
      endif
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
    
    !---------------------------------------------------------------------------------
    Deallocate(xtrack,ytrack,gheight,ww,uu,vv,Kw,Ku,Kv)
    !---------------------------------------------------------------------------------
    !--------------------------------------------------------------------------------- 
    !
    write(*,*) 'Info: Getting press,temp,rho_d,qvapor,tke'
    env%press(1:env%nx,1:env%ny,1:env%nz)=wrf%press(ix1:ix2,iy1:iy2,iz1:iz2,it) ! Pa 
    env%temp(1:env%nx,1:env%ny,1:env%nz)=wrf%temp(ix1:ix2,iy1:iy2,iz1:iz2,it) ! K 
    env%qvapor(1:env%nx,1:env%ny,1:env%nz)=wrf%qvapor(ix1:ix2,iy1:iy2,iz1:iz2,it) ! kg/kg 
    env%tke(1:env%nx,1:env%ny,1:env%nz)=wrf%tke(ix1:ix2,iy1:iy2,iz1:iz2,it) ! m^2/s^2    
    if (MaxVal(wrf%rho_d)<-888.d0) then
      ! if ALT not present in the input file                                           
      env%rho_d=(env%press*env%qvapor)/(eps+env%qvapor) ! e water vapor pressure in Pa
      env%rho_d=(env%press-env%rho_d)/(Rd * env%temp)   !  kg/m^3    
    else
      !if ALT present in the input file                          
      env%rho_d(1:env%nx,1:env%ny,1:env%nz)=wrf%rho_d(ix1:ix2,iy1:iy2,iz1:iz2,it)
      ! m^3/kg 
      !env%rho_d=1.d0/env%rho_d !  kg/m^3  
    endif
    !
    ! convert T from K to C 
    env%temp=env%temp-273.15d0 ! C 
    ! convert press frpm Pa to mb   
    env%press=env%press*1.d-2  ! mb 
    
    write(*,*) 'press [mb]',MinVal(env%press),MaxVal(env%press)
    write(*,*) 'temp [C]',MinVal(env%temp),MaxVal(env%temp)
    write(*,*) 'rho_d [kg/m^3]',MinVal(env%rho_d),MaxVal(env%rho_d)
    write(*,*) 'tke [m^2/s^2]',MinVal(env%tke),MaxVal(env%tke)
    ! 
  return
  end subroutine get_env_vars_rams
  !!
  subroutine ReadInpRAMS_text(InpFile,nz,dx,hgtm)
  Use wrf_var_mod
  Implicit None
  !
  Character(len=*),Intent(in)    :: InpFile
  Integer,Intent(in)             :: nz
  Real*8,Intent(out)             :: dx
  Real*8,Intent(out)             :: hgtm(nz)
  !
  integer                        :: nz1,iz,n
  integer,parameter              :: nrow=10272 !should be changed depending on InpFile 
  integer,parameter              :: id=1
  character(len=220)             :: inputline='#'
    !
    dx=-999.d0
    hgtm=0.d0
    
    open(unit=id,file=InpFile)
    !
    inputline='#'
    A1: do n=1,nrow
      read(id,'(a)') inputline
      if (inputline(1:9).eq.'__deltaxn') then
        read(id,'(a)') inputline
        read(id,'(a)') inputline
        read(id,'(a)') inputline
        read(id,'(a)') inputline
        read(inputline,*) dx
      endif
    
      if (inputline(1:7).eq.'__ztn01') then
        read(id,'(a)') inputline
        read(inputline,*) nz1
        write(*,*) 'NZ=',nz1, nz
        do iz=1,nz
          read(id,'(a)') inputline
          read(inputline,*) hgtm(iz)
        enddo
        exit A1
      endif
    
    enddo A1
    
    close(id)
    
  return
  end subroutine ReadInpRAMS_text
  !---------------------------------------------------------------------
  subroutine ReadInpSAM_dim(InpFile,str,status)
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
      Case ('x')        ;  str%nx    = length
      Case ('y')        ;  str%ny    = length
      Case ('z')       ;  str%nz    = length
      Case ('time')       ;  str%nt    = length
      End Select
      !
    enddo   ! iDim
    !
    !str%nt    = 1
    str%nxp1  = str%nx
    str%nyp1  = str%ny
    str%nzp1  = str%nz
    
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
    If (str%nxp1/=str%nx) then
       write(*,*) 'Problem in input dimensions nx and nxp1'
       write(*,*) 'nx=nxp1 for SAM input'
       Call Exit(1)
    EndIf
    !
    If (str%nyp1/=str%ny) then
      write(*,*) 'Problem in input dimensions ny and nyp1'
      write(*,*) 'nx=ny, ny=nyp1 for SAM input'
      Call Exit(1)
    EndIf
    !
  return
  end subroutine ReadInpSAM_dim
  !
  !---------------------------------------------------------------------
  !subroutine ReadInpSAM_var(InpFile,str,status)
  subroutine ReadInpSAM_var(InpFile,InpProfile,str,status)
  Use netcdf
  Use wrf_var_mod
  Use wrf_rvar_mod
  Implicit None
  !
  Character(len=*),Intent(in)                              :: InpFile
  Character(len=*),Intent(in)                              :: InpProfile
  Type(wrf_var),Intent(InOut)                              :: str
  Integer, Intent(out)                                     :: status
  !
  Type(wrf_rvar)                                           :: strr
  Integer                                                  :: ncid,iDim,length
  Integer                                                  :: nDims,nVars,iVar
  Character(len=nf90_max_name)                             :: name, err_msg
  Character(len=nf90_max_name), Allocatable , Dimension(:) :: var_names
  
  Integer :: nzp,ltime
  Real*4, Allocatable , Dimension(:,:) :: rho_p
  !Integer, Allocatable, Dimension(:)                       :: dim_lengths
  !Character(len=nf90_max_name), Allocatable , Dimension(:) :: dim_names
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
    call allocate_wrf_rvar(strr)
    call initialize_wrf_rvar(strr)
    ! 
    status=  nf90_open(Trim(InpFile), NF90_NOWRITE, ncid)
    if (status.ne.0) then ; err_msg='Error in nf90_open' ; goto 999 ; endif
    !
    status= nf90_inquire(ncid, nDims, nVars)
    If (status/=0) Then ; err_msg = 'Error in nf90_inquire' ; goto 999 ; endif
    !
    !
    strr%p=0.e0
    str%p=0.d0
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
      !
      Case ('U')
      status=nf90_get_var(ncid,iVar,strr%u)
      str%u = dble (strr%u)
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: U'
        Goto 999
      End If
      !
      Case ('V')
      status=nf90_get_var(ncid,iVar,strr%v)
      str%v = dble (strr%v)
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: V'
        Goto 999
      End If
      !
      Case ('W')
      status=nf90_get_var(ncid,iVar,strr%w)
      str%w = dble (strr%w)
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: W'
        Goto 999
      End If
      !
      Case ('p')
      status=nf90_get_var(ncid,iVar,strr%pb(1,1,:,1))
      str%pb(1,1,:,1) = dble (strr%pb(1,1,:,1))
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: P'
        Goto 999
      End If
      !
      Case ('PP')
      status=nf90_get_var(ncid,iVar,strr%p)
      str%p = dble (strr%p)
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: PP'
        Goto 999
      End If
      !
      Case ('TABS')
      status=nf90_get_var(ncid,iVar,strr%temp)
      str%temp = dble (strr%temp)
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: TABS'
        Goto 999
      End If
      !
      Case ('z')
      status=nf90_get_var(ncid,iVar,strr%hgtm)
      str%hgtm = dble (strr%hgtm)
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: z'
        Goto 999
      End If
      !
      Case ('QV')
      status=nf90_get_var(ncid,iVar,strr%QVAPOR)
      str%QVAPOR = dble (strr%QVAPOR)
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: QVAPOR'
        Goto 999
      End If
      !
      Case ('x')
      status=nf90_get_var(ncid,iVar,strr%xlong(:,1,1))
      str%xlong(:,1,1) = dble (strr%xlong(:,1,1))
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: x'
        Goto 999
      End If
      !
      Case ('y')
      status=nf90_get_var(ncid,iVar,strr%xlat(1,:,1))
      str%xlat(1,:,1) = dble (strr%xlat(1,:,1)) !- fixed by oue Apr 2020
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: y'
        Goto 999
      End If
      !
      !
      End Select
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
    
    !------------------------------------
    !get dx/dy
    str%dx=1./(str%xlong(2,1,1)-str%xlong(1,1,1)) !- fixed by oue Apr 2020
    str%dy=1./(str%xlat(1,2,1)-str%xlat(1,1,1)) !- fixed by oue Apr 2020
    write(*,*) 'SAM dx,dy',str%dx,str%dy
    !------------------------------------
    
    !------------------------------------
    ! Read profile data
    !------------------------------------
    status=  nf90_open(Trim(InpProfile), NF90_NOWRITE, ncid)
    if (status.ne.0) then ; err_msg='No profiling data' ; goto 998 ; endif
    !!
    status= nf90_inquire(ncid, nDims, nVars)
    If (status/=0) Then ; err_msg = 'Error in nf90_inquire' ; goto 998 ; endif
    !!
    !Allocate(dim_names(1:nDims))
    !Allocate(dim_lengths(1:nDims))
    !
    ! Get info on the dimensions
    ! Loop over the number of dimensions
    !   
    nzp = str%nz ; ltime = 3060
    !
    Do iDim = 1,nDims
      !!
      status=nf90_inquire_dimension(ncid, idim, name=name, len=length)
      If (status/=0) Then ; err_msg = 'Failure to read profiling data' ; Goto 998 ; End If
      !!
      Select Case (name)
      Case ('z')       ;  nzp    = length
      Case ('time')    ;  ltime    = length
      End Select
      !!
    enddo   ! iDim
    !!
    if (str%nz /= nzp) then ; err_msg='Mismatch z dimension; Failure to read profiling data ' ; goto 998 ; endif!
    !
    Allocate(rho_p(nzp,ltime))
    !
    Allocate(var_names(1:nVars))
    !!
    Do  iVar=1,nVars
      !!
      status= nf90_inquire_variable(ncid, iVar, name=name)
      If (status/=0) Then
        err_msg = 'Error in my_nf90_inquire_variable: Failure to read profiling data'//name
        Goto 998
      Endif
      !!
      Select Case (Trim(Adjustl(name)))
      !!
      Case ('RHO')
      status=nf90_get_var(ncid,iVar,rho_p) 
      write(*,*) 'Read RHO from first profile data, assuming that RHO is constant with time.' 
      If (status/=0) Then
        err_msg = 'Cannot read my_nf90_get_var: RHO; RHO will be calculated.'
        Goto 998
      End If
      !!
      End Select
      !!
    EndDo
    !
998 If (status.Ne.0) Then
      write(*,*) err_msg
      status=nf90_close(ncid)
      !Call Exit(1)
    Else
      status=nf90_close(ncid)
      str%rho_d(1,1,1:str%nz,1) = dble(rho_p(1:str%nz,1)) 
    Endif
    !
    call deallocate_wrf_rvar(strr)
    !
  return
  end subroutine ReadInpSAM_var
  !
  !-------------------------------------------------------------------------------------------------------------------
  !!! Added by oue, for SAM warm bin
  subroutine ReadInpSAM_MP_PHYSICS_70(InpFile,nkr,str,status)
  Use netcdf
  Use typeSizes
  Use wrf_rvar_mod
  Use wrf_var_mod
  Implicit None
  !
  Character(len=*),Intent(in)                              :: InpFile
  Integer,Intent(In)                                       :: nkr ! number of bins
  Type(wrf_var_mp70),Intent(InOut)                         :: str
  Integer,Intent(out)                                      :: status
  !
  Type(wrf_rvar_mp70)                                      :: strr
  integer                                                  :: nDims, nVars
  integer                                                  :: ncid, iVar
  Character(len=nf90_max_name), Allocatable , Dimension(:) :: var_names
  Character(len=nf90_max_name)                             :: name, err_msg
  Integer                     :: im,in
    !
    !- Massages 
    write(*,*) 'Selected microphysics: 70 (SAM liquid bin)'
    write(*,*) 'Required hydrometeor parameters: Bin data for mass, M0*'
    write(*,*) 'Required hydrometeor parameters: Bin data for number, N0*'
    write(*,*) 'Required hydrometeor parameters: QC, QR, NC, NR, NA'
    write(*,*) 'Reading hydrometeor mixing ratios and number concentrations...'
    im=0 ; in=0 !counting the number of bins
    !
    strr%nt = str%nt
    strr%nx = str%nx
    strr%ny = str%ny
    strr%nz = str%nz
    strr%nbins = str%nbins
    !
    call allocate_wrf_rvar_mp70(strr)
    call initialize_wrf_rvar_mp70(strr)
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
      Case ('M01')
      status=nf90_get_var(ncid,iVar,strr%fm1(:,:,:,:,1)) ; im=im+1 !write(*,*) 'M01 is found.'
      str%fm1(:,:,:,:,1) = dble(strr%fm1(:,:,:,:,1))
      If (status/=0) Then ; err_msg = 'Error in my_nf90_get_var: MXX' ; Goto 999 ; End If
      !
      Case ('M02')
      status=nf90_get_var(ncid,iVar,strr%fm1(:,:,:,:,2)) ; im=im+1 ! write(*,*) 'M02 is found.'
      str%fm1(:,:,:,:,2) = dble(strr%fm1(:,:,:,:,2))
      If (status/=0) Then ; err_msg = 'Error in my_nf90_get_var: MXX' ; Goto 999 ; End If
      !
      Case ('M03')
      status=nf90_get_var(ncid,iVar,strr%fm1(:,:,:,:,3)) ; im=im+1 ! write(*,*) 'M03 is found.'
      str%fm1(:,:,:,:,3) = dble(strr%fm1(:,:,:,:,3))
      If (status/=0) Then ; err_msg = 'Error in my_nf90_get_var: MXX' ; Goto 999 ; End If
      !
      Case ('M04')
      status=nf90_get_var(ncid,iVar,strr%fm1(:,:,:,:,4)) ;  im=im+1 !write(*,*) 'M04 is found.'
      str%fm1(:,:,:,:,4) = dble(strr%fm1(:,:,:,:,4))
      If (status/=0) Then ; err_msg = 'Error in my_nf90_get_var: MXX' ; Goto 999 ; End If
      !
      Case ('M05')
      status=nf90_get_var(ncid,iVar,strr%fm1(:,:,:,:,5)) ; im=im+1 ! write(*,*) 'M05 is found.'
      str%fm1(:,:,:,:,5) = dble(strr%fm1(:,:,:,:,5))
      If (status/=0) Then ; err_msg = 'Error in my_nf90_get_var: MXX' ; Goto 999 ; End If
      !
      Case ('M06')
      status=nf90_get_var(ncid,iVar,strr%fm1(:,:,:,:,6)) ; im=im+1 ! write(*,*) 'M06 is found.'
      str%fm1(:,:,:,:,6) = dble(strr%fm1(:,:,:,:,6))
      If (status/=0) Then ; err_msg = 'Error in my_nf90_get_var: MXX' ; Goto 999 ; End If
      !
      Case ('M07')
      status=nf90_get_var(ncid,iVar,strr%fm1(:,:,:,:,7)) ; im=im+1 ! write(*,*) 'M07 is found.'
      str%fm1(:,:,:,:,7) = dble(strr%fm1(:,:,:,:,7))
      If (status/=0) Then ; err_msg = 'Error in my_nf90_get_var: MXX' ; Goto 999 ; End If
      !
      Case ('M08')
      status=nf90_get_var(ncid,iVar,strr%fm1(:,:,:,:,8)) ; im=im+1 ! write(*,*) 'M08 is found.'
      str%fm1(:,:,:,:,8) = dble(strr%fm1(:,:,:,:,8))
      If (status/=0) Then ; err_msg = 'Error in my_nf90_get_var: MXX' ; Goto 999 ; End If
      !
      Case ('M09')
      status=nf90_get_var(ncid,iVar,strr%fm1(:,:,:,:,9)) ; im=im+1 ! write(*,*) 'M09 is found.'
      str%fm1(:,:,:,:,9) = dble(strr%fm1(:,:,:,:,9))
      If (status/=0) Then ; err_msg = 'Error in my_nf90_get_var: MXX' ; Goto 999 ; End If
      !
      Case ('M10')
      status=nf90_get_var(ncid,iVar,strr%fm1(:,:,:,:,10)) ; im=im+1 ! write(*,*) 'M10 is found.'
      str%fm1(:,:,:,:,10) = dble(strr%fm1(:,:,:,:,10))
      If (status/=0) Then ; err_msg = 'Error in my_nf90_get_var: MXX' ; Goto 999 ; End If
      !
      Case ('M11')
      status=nf90_get_var(ncid,iVar,strr%fm1(:,:,:,:,11)) ; im=im+1 ! write(*,*) 'M11 is found.'
      str%fm1(:,:,:,:,11) = dble(strr%fm1(:,:,:,:,11))
      If (status/=0) Then ; err_msg = 'Error in my_nf90_get_var: MXX' ; Goto 999 ; End If
      !
      Case ('M12')
      status=nf90_get_var(ncid,iVar,strr%fm1(:,:,:,:,12)) ; im=im+1 ! write(*,*) 'M12 is found.'
      str%fm1(:,:,:,:,12) = dble(strr%fm1(:,:,:,:,12))
      If (status/=0) Then ; err_msg = 'Error in my_nf90_get_var: MXX' ; Goto 999 ; End If
      !
      Case ('M13')
      status=nf90_get_var(ncid,iVar,strr%fm1(:,:,:,:,13)) ; im=im+1 ! write(*,*) 'M13 is found.'
      str%fm1(:,:,:,:,13) = dble(strr%fm1(:,:,:,:,13))
      If (status/=0) Then ; err_msg = 'Error in my_nf90_get_var: MXX' ; Goto 999 ; End If
      !
      Case ('M14')
      status=nf90_get_var(ncid,iVar,strr%fm1(:,:,:,:,14)) ; im=im+1 ! write(*,*) 'M14 is found.'
      str%fm1(:,:,:,:,14) = dble(strr%fm1(:,:,:,:,14))
      If (status/=0) Then ; err_msg = 'Error in my_nf90_get_var: MXX' ; Goto 999 ; End If
      !
      Case ('M15')
      status=nf90_get_var(ncid,iVar,strr%fm1(:,:,:,:,15)) ; im=im+1 ! write(*,*) 'M15 is found.'
      str%fm1(:,:,:,:,15) = dble(strr%fm1(:,:,:,:,15))
      If (status/=0) Then ; err_msg = 'Error in my_nf90_get_var: MXX' ; Goto 999 ; End If
      !
      Case ('M16')
      status=nf90_get_var(ncid,iVar,strr%fm1(:,:,:,:,16)) ; im=im+1 ! write(*,*) 'M16 is found.'
      str%fm1(:,:,:,:,16) = dble(strr%fm1(:,:,:,:,16))
      If (status/=0) Then ; err_msg = 'Error in my_nf90_get_var: MXX' ; Goto 999 ; End If
      !
      Case ('M17')
      status=nf90_get_var(ncid,iVar,strr%fm1(:,:,:,:,17)) ; im=im+1 ! write(*,*) 'M17 is found.'
      str%fm1(:,:,:,:,17) = dble(strr%fm1(:,:,:,:,17))
      If (status/=0) Then ; err_msg = 'Error in my_nf90_get_var: MXX' ; Goto 999 ; End If
      !
      Case ('M18')
      status=nf90_get_var(ncid,iVar,strr%fm1(:,:,:,:,18)) ; im=im+1 ! write(*,*) 'M18 is found.'
      str%fm1(:,:,:,:,18) = dble(strr%fm1(:,:,:,:,18))
      If (status/=0) Then ; err_msg = 'Error in my_nf90_get_var: MXX' ; Goto 999 ; End If
      !
      Case ('M19')
      status=nf90_get_var(ncid,iVar,strr%fm1(:,:,:,:,19)) ; im=im+1 ! write(*,*) 'M19 is found.'
      str%fm1(:,:,:,:,19) = dble(strr%fm1(:,:,:,:,19))
      If (status/=0) Then ; err_msg = 'Error in my_nf90_get_var: MXX' ; Goto 999 ; End If
      !
      Case ('M20')
      status=nf90_get_var(ncid,iVar,strr%fm1(:,:,:,:,20)) ; im=im+1 ! write(*,*) 'M20 is found.'
      str%fm1(:,:,:,:,20) = dble(strr%fm1(:,:,:,:,20))
      If (status/=0) Then ; err_msg = 'Error in my_nf90_get_var: MXX' ; Goto 999 ; End If
      !
      Case ('M21')
      status=nf90_get_var(ncid,iVar,strr%fm1(:,:,:,:,21)) ; im=im+1 ! write(*,*) 'M21 is found.'
      str%fm1(:,:,:,:,21) = dble(strr%fm1(:,:,:,:,21))
      If (status/=0) Then ; err_msg = 'Error in my_nf90_get_var: MXX' ; Goto 999 ; End If
      !
      Case ('M22')
      status=nf90_get_var(ncid,iVar,strr%fm1(:,:,:,:,22)) ; im=im+1 ! write(*,*) 'M22 is found.'
      str%fm1(:,:,:,:,22) = dble(strr%fm1(:,:,:,:,22))
      If (status/=0) Then ; err_msg = 'Error in my_nf90_get_var: MXX' ; Goto 999 ; End If
      !
      Case ('M23')
      status=nf90_get_var(ncid,iVar,strr%fm1(:,:,:,:,23)) ; im=im+1 ! write(*,*) 'M23 is found.'
      str%fm1(:,:,:,:,23) = dble(strr%fm1(:,:,:,:,23))
      If (status/=0) Then ; err_msg = 'Error in my_nf90_get_var: MXX' ; Goto 999 ; End If
      !
      Case ('M24')
      status=nf90_get_var(ncid,iVar,strr%fm1(:,:,:,:,24)) ; im=im+1 ! write(*,*) 'M24 is found.'
      str%fm1(:,:,:,:,24) = dble(strr%fm1(:,:,:,:,24))
      If (status/=0) Then ; err_msg = 'Error in my_nf90_get_var: MXX' ; Goto 999 ; End If
      !
      Case ('M25')
      status=nf90_get_var(ncid,iVar,strr%fm1(:,:,:,:,25)) ; im=im+1 ! write(*,*) 'M25 is found.'
      str%fm1(:,:,:,:,25) = dble(strr%fm1(:,:,:,:,25))
      If (status/=0) Then ; err_msg = 'Error in my_nf90_get_var: MXX' ; Goto 999 ; End If
      !
      Case ('M26')
      status=nf90_get_var(ncid,iVar,strr%fm1(:,:,:,:,26)) ; im=im+1 ! write(*,*) 'M26 is found.'
      str%fm1(:,:,:,:,26) = dble(strr%fm1(:,:,:,:,26))
      If (status/=0) Then ; err_msg = 'Error in my_nf90_get_var: MXX' ; Goto 999 ; End If
      !
      Case ('M27')
      status=nf90_get_var(ncid,iVar,strr%fm1(:,:,:,:,27)) ; im=im+1 ! write(*,*) 'M27 is found.'
      str%fm1(:,:,:,:,27) = dble(strr%fm1(:,:,:,:,27))
      If (status/=0) Then ; err_msg = 'Error in my_nf90_get_var: MXX' ; Goto 999 ; End If
      !
      Case ('M28')
      status=nf90_get_var(ncid,iVar,strr%fm1(:,:,:,:,28)) ; im=im+1 ! write(*,*) 'M28 is found.'
      str%fm1(:,:,:,:,28) = dble(strr%fm1(:,:,:,:,28))
      If (status/=0) Then ; err_msg = 'Error in my_nf90_get_var: MXX' ; Goto 999 ; End If
      !
      Case ('M29')
      status=nf90_get_var(ncid,iVar,strr%fm1(:,:,:,:,29)) ; im=im+1 ! write(*,*) 'M29 is found.'
      str%fm1(:,:,:,:,29) = dble(strr%fm1(:,:,:,:,29))
      If (status/=0) Then ; err_msg = 'Error in my_nf90_get_var: MXX' ; Goto 999 ; End If
      !
      Case ('M30')
      status=nf90_get_var(ncid,iVar,strr%fm1(:,:,:,:,30)) ; im=im+1 ! write(*,*) 'M30 is found.'
      str%fm1(:,:,:,:,30) = dble(strr%fm1(:,:,:,:,30))
      If (status/=0) Then ; err_msg = 'Error in my_nf90_get_var: MXX' ; Goto 999 ; End If
      !
      Case ('M31')
      status=nf90_get_var(ncid,iVar,strr%fm1(:,:,:,:,31)) ; im=im+1 ! write(*,*) 'M31 is found.'
      str%fm1(:,:,:,:,31) = dble(strr%fm1(:,:,:,:,31))
      If (status/=0) Then ; err_msg = 'Error in my_nf90_get_var: MXX' ; Goto 999 ; End If
      !
      Case ('M32')
      status=nf90_get_var(ncid,iVar,strr%fm1(:,:,:,:,32)) ; im=im+1 ! write(*,*) 'M32 is found.'
      str%fm1(:,:,:,:,32) = dble(strr%fm1(:,:,:,:,32))
      If (status/=0) Then ; err_msg = 'Error in my_nf90_get_var: MXX' ; Goto 999 ; End If
      !
      Case ('M33')
      status=nf90_get_var(ncid,iVar,strr%fm1(:,:,:,:,33)) ; im=im+1 ! write(*,*) 'M33 is found.'
      str%fm1(:,:,:,:,33) = dble(strr%fm1(:,:,:,:,33))
      If (status/=0) Then ; err_msg = 'Error in my_nf90_get_var: MXX' ; Goto 999 ; End If
      !
      Case ('M34')
      status=nf90_get_var(ncid,iVar,strr%fm1(:,:,:,:,34)) ; im=im+1 ! write(*,*) 'M34 is found.'
      str%fm1(:,:,:,:,34) = dble(strr%fm1(:,:,:,:,34))
      If (status/=0) Then ; err_msg = 'Error in my_nf90_get_var: MXX' ; Goto 999 ; End If
      !
      Case ('M35')
      status=nf90_get_var(ncid,iVar,strr%fm1(:,:,:,:,35)) ; im=im+1 ! write(*,*) 'M35 is found.'
      str%fm1(:,:,:,:,35) = dble(strr%fm1(:,:,:,:,35))
      If (status/=0) Then ; err_msg = 'Error in my_nf90_get_var: MXX' ; Goto 999 ; End If
      !
      Case ('M36')
      status=nf90_get_var(ncid,iVar,strr%fm1(:,:,:,:,36)) ; im=im+1 ! write(*,*) 'M36 is found.'
      str%fm1(:,:,:,:,36) = dble(strr%fm1(:,:,:,:,36))
      If (status/=0) Then ; err_msg = 'Error in my_nf90_get_var: MXX' ; Goto 999 ; End If
      !
      Case ('N01')
      status=nf90_get_var(ncid,iVar,strr%fn1(:,:,:,:,1)) ; in=in+1 ! write(*,*) 'N01 is found.'
      str%fn1(:,:,:,:,1) = dble(strr%fn1(:,:,:,:,1))
      If (status/=0) Then ; err_msg = 'Error in my_nf90_get_var: NXX' ; Goto 999 ; End If
      !
      Case ('N02')
      status=nf90_get_var(ncid,iVar,strr%fn1(:,:,:,:,2)) ; in=in+1 ! write(*,*) 'N02 is found.'
      str%fn1(:,:,:,:,2) = dble(strr%fn1(:,:,:,:,2))
      If (status/=0) Then ; err_msg = 'Error in my_nf90_get_var: NXX' ; Goto 999 ; End If
      !
      Case ('N03')
      status=nf90_get_var(ncid,iVar,strr%fn1(:,:,:,:,3)) ; in=in+1 ! write(*,*) 'N03 is found.'
      str%fn1(:,:,:,:,3) = dble(strr%fn1(:,:,:,:,3))
      If (status/=0) Then ; err_msg = 'Error in my_nf90_get_var: NXX' ; Goto 999 ; End If
      !
      Case ('N04')
      status=nf90_get_var(ncid,iVar,strr%fn1(:,:,:,:,4)) ; in=in+1 ! write(*,*) 'N04 is found.'
      str%fn1(:,:,:,:,4) = dble(strr%fn1(:,:,:,:,4))
      If (status/=0) Then ; err_msg = 'Error in my_nf90_get_var: NXX' ; Goto 999 ; End If
      !
      Case ('N05')
      status=nf90_get_var(ncid,iVar,strr%fn1(:,:,:,:,5)) ; in=in+1 ! write(*,*) 'N05 is found.'
      str%fn1(:,:,:,:,5) = dble(strr%fn1(:,:,:,:,5))
      If (status/=0) Then ; err_msg = 'Error in my_nf90_get_var: NXX' ; Goto 999 ; End If
      !
      Case ('N06')
      status=nf90_get_var(ncid,iVar,strr%fn1(:,:,:,:,6)) ; in=in+1 ! write(*,*) 'N06 is found.'
      str%fn1(:,:,:,:,6) = dble(strr%fn1(:,:,:,:,6))
      If (status/=0) Then ; err_msg = 'Error in my_nf90_get_var: NXX' ; Goto 999 ; End If
      !
      Case ('N07')
      status=nf90_get_var(ncid,iVar,strr%fn1(:,:,:,:,7)) ; in=in+1 ! write(*,*) 'N07 is found.'
      str%fn1(:,:,:,:,7) = dble(strr%fn1(:,:,:,:,7))
      If (status/=0) Then ; err_msg = 'Error in my_nf90_get_var: NXX' ; Goto 999 ; End If
      !
      Case ('N08')
      status=nf90_get_var(ncid,iVar,strr%fn1(:,:,:,:,8)) ; in=in+1 ! write(*,*) 'N08 is found.'
      str%fn1(:,:,:,:,8) = dble(strr%fn1(:,:,:,:,8))
      If (status/=0) Then ; err_msg = 'Error in my_nf90_get_var: NXX' ; Goto 999 ; End If
      !
      Case ('N09')
      status=nf90_get_var(ncid,iVar,strr%fn1(:,:,:,:,9)) ; in=in+1 ! write(*,*) 'N09 is found.'
      str%fn1(:,:,:,:,9) = dble(strr%fn1(:,:,:,:,9))
      If (status/=0) Then ; err_msg = 'Error in my_nf90_get_var: NXX' ; Goto 999 ; End If
      !
      Case ('N10')
      status=nf90_get_var(ncid,iVar,strr%fn1(:,:,:,:,10)) ; in=in+1 ! write(*,*) 'N10 is found.'
      str%fn1(:,:,:,:,10) = dble(strr%fn1(:,:,:,:,10))
      If (status/=0) Then ; err_msg = 'Error in my_nf90_get_var: NXX' ; Goto 999 ; End If
      !
      Case ('N11')
      status=nf90_get_var(ncid,iVar,strr%fn1(:,:,:,:,11)) ; in=in+1 ! write(*,*) 'N11 is found.'
      str%fn1(:,:,:,:,11) = dble(strr%fn1(:,:,:,:,11))
      If (status/=0) Then ; err_msg = 'Error in my_nf90_get_var: NXX' ; Goto 999 ; End If
      !
      Case ('N12')
      status=nf90_get_var(ncid,iVar,strr%fn1(:,:,:,:,12)) ; in=in+1 ! write(*,*) 'N12 is found.'
      str%fn1(:,:,:,:,12) = dble(strr%fn1(:,:,:,:,12))
      If (status/=0) Then ; err_msg = 'Error in my_nf90_get_var: NXX' ; Goto 999 ; End If
      !
      Case ('N13')
      status=nf90_get_var(ncid,iVar,strr%fn1(:,:,:,:,13)) ; in=in+1 ! write(*,*) 'N13 is found.'
      str%fn1(:,:,:,:,13) = dble(strr%fn1(:,:,:,:,13))
      If (status/=0) Then ; err_msg = 'Error in my_nf90_get_var: NXX' ; Goto 999 ; End If
      !
      Case ('N14')
      status=nf90_get_var(ncid,iVar,strr%fn1(:,:,:,:,14)) ; in=in+1 ! write(*,*) 'N14 is found.'
      str%fn1(:,:,:,:,14) = dble(strr%fn1(:,:,:,:,14))
      If (status/=0) Then ; err_msg = 'Error in my_nf90_get_var: NXX' ; Goto 999 ; End If
      !
      Case ('N15')
      status=nf90_get_var(ncid,iVar,strr%fn1(:,:,:,:,15)) ; in=in+1 ! write(*,*) 'N15 is found.'
      str%fn1(:,:,:,:,15) = dble(strr%fn1(:,:,:,:,15))
      If (status/=0) Then ; err_msg = 'Error in my_nf90_get_var: NXX' ; Goto 999 ; End If
      !
      Case ('N16')
      status=nf90_get_var(ncid,iVar,strr%fn1(:,:,:,:,16)) ; in=in+1 ! write(*,*) 'N16 is found.'
      str%fn1(:,:,:,:,16) = dble(strr%fn1(:,:,:,:,16))
      If (status/=0) Then ; err_msg = 'Error in my_nf90_get_var: NXX' ; Goto 999 ; End If
      !
      Case ('N17')
      status=nf90_get_var(ncid,iVar,strr%fn1(:,:,:,:,17)) ; in=in+1 ! write(*,*) 'N17 is found.'
      str%fn1(:,:,:,:,17) = dble(strr%fn1(:,:,:,:,17))
      If (status/=0) Then ; err_msg = 'Error in my_nf90_get_var: NXX' ; Goto 999 ; End If
      !
      Case ('N18')
      status=nf90_get_var(ncid,iVar,strr%fn1(:,:,:,:,18)) ; in=in+1 ! write(*,*) 'N18 is found.'
      str%fn1(:,:,:,:,18) = dble(strr%fn1(:,:,:,:,18))
      If (status/=0) Then ; err_msg = 'Error in my_nf90_get_var: NXX' ; Goto 999 ; End If
      !
      Case ('N19')
      status=nf90_get_var(ncid,iVar,strr%fn1(:,:,:,:,19)) ; in=in+1 ! write(*,*) 'N19 is found.'
      str%fn1(:,:,:,:,19) = dble(strr%fn1(:,:,:,:,19))
      If (status/=0) Then ; err_msg = 'Error in my_nf90_get_var: NXX' ; Goto 999 ; End If
      !
      Case ('N20')
      status=nf90_get_var(ncid,iVar,strr%fn1(:,:,:,:,20)) ; in=in+1 ! write(*,*) 'N20 is found.'
      str%fn1(:,:,:,:,20) = dble(strr%fn1(:,:,:,:,20))
      If (status/=0) Then ; err_msg = 'Errr in my_nf90_get_var: NXX' ; Goto 999 ; End If
      !
      Case ('N21')
      status=nf90_get_var(ncid,iVar,strr%fn1(:,:,:,:,21)) ; in=in+1 ! write(*,*) 'N21 is found.'
      str%fn1(:,:,:,:,21) = dble(strr%fn1(:,:,:,:,21))
      If (status/=0) Then ; err_msg = 'Error in my_nf90_get_var: NXX' ; Goto 999 ; End If
      !
      Case ('N22')
      status=nf90_get_var(ncid,iVar,strr%fn1(:,:,:,:,22)) ; in=in+1 ! write(*,*) 'N22 is found.'
      str%fn1(:,:,:,:,22) = dble(strr%fn1(:,:,:,:,22))
      If (status/=0) Then ; err_msg = 'Error in my_nf90_get_var: NXX' ; Goto 999 ; End If
      !
      Case ('N23')
      status=nf90_get_var(ncid,iVar,strr%fn1(:,:,:,:,23)) ; in=in+1 ! write(*,*) 'N23 is found.'
      str%fn1(:,:,:,:,23) = dble(strr%fn1(:,:,:,:,23))
      If (status/=0) Then ; err_msg = 'Error in my_nf90_get_var: NXX' ; Goto 999 ; End If
      !
      Case ('N24')
      status=nf90_get_var(ncid,iVar,strr%fn1(:,:,:,:,24)) ; in=in+1 ! write(*,*) 'N24 is found.'
      str%fn1(:,:,:,:,24) = dble(strr%fn1(:,:,:,:,24))
      If (status/=0) Then ; err_msg = 'Error in my_nf90_get_var: NXX' ; Goto 999 ; End If
      !
      Case ('N25')
      status=nf90_get_var(ncid,iVar,strr%fn1(:,:,:,:,25)) ; in=in+1 ! write(*,*) 'N25 is found.'
      str%fn1(:,:,:,:,25) = dble(strr%fn1(:,:,:,:,25))
      If (status/=0) Then ; err_msg = 'Error in my_nf90_get_var: NXX' ; Goto 999 ; End If
      !
      Case ('N26')
      status=nf90_get_var(ncid,iVar,strr%fn1(:,:,:,:,26)) ; in=in+1 ! write(*,*) 'N26 is found.'
      str%fn1(:,:,:,:,26) = dble(strr%fn1(:,:,:,:,26))
      If (status/=0) Then ; err_msg = 'Error in my_nf90_get_var: NXX' ; Goto 999 ; End If
      !
      Case ('N27')
      status=nf90_get_var(ncid,iVar,strr%fn1(:,:,:,:,27)) ; in=in+1 ! write(*,*) 'N27 is found.'
      str%fn1(:,:,:,:,27) = dble(strr%fn1(:,:,:,:,27))
      If (status/=0) Then ; err_msg = 'Error in my_nf90_get_var: NXX' ; Goto 999 ; End If
      !
      Case ('N28')
      status=nf90_get_var(ncid,iVar,strr%fn1(:,:,:,:,28)) ; in=in+1 ! write(*,*) 'N28 is found.'
      str%fn1(:,:,:,:,28) = dble(strr%fn1(:,:,:,:,28))
      If (status/=0) Then ; err_msg = 'Error in my_nf90_get_var: NXX' ; Goto 999 ; End If
      !
      Case ('N29')
      status=nf90_get_var(ncid,iVar,strr%fn1(:,:,:,:,29)) ; in=in+1 ! write(*,*) 'N29 is found.'
      str%fn1(:,:,:,:,29) = dble(strr%fn1(:,:,:,:,29))
      If (status/=0) Then ; err_msg = 'Error in my_nf90_get_var: NXX' ; Goto 999 ; End If
      !
      Case ('N30')
      status=nf90_get_var(ncid,iVar,strr%fn1(:,:,:,:,30)) ; in=in+1 ! write(*,*) 'N30 is found.'
      str%fn1(:,:,:,:,30) = dble(strr%fn1(:,:,:,:,30))
      If (status/=0) Then ; err_msg = 'Error in my_nf90_get_var: NXX' ; Goto 999 ; End If
      !
      Case ('N31')
      status=nf90_get_var(ncid,iVar,strr%fn1(:,:,:,:,31)) ; in=in+1 ! write(*,*) 'N31 is found.'
      str%fn1(:,:,:,:,31) = dble(strr%fn1(:,:,:,:,31))
      If (status/=0) Then ; err_msg = 'Error in my_nf90_get_var: NXX' ; Goto 999 ; End If
      !
      Case ('N32')
      status=nf90_get_var(ncid,iVar,strr%fn1(:,:,:,:,32)) ; in=in+1 ! write(*,*) 'N32 is found.'
      str%fn1(:,:,:,:,32) = dble(strr%fn1(:,:,:,:,32))
      If (status/=0) Then ; err_msg = 'Error in my_nf90_get_var: NXX' ; Goto 999 ; End If
      !
      Case ('N33')
      status=nf90_get_var(ncid,iVar,strr%fn1(:,:,:,:,33)) ; in=in+1 ! write(*,*) 'N33 is found.'
      str%fn1(:,:,:,:,33) = dble(strr%fn1(:,:,:,:,33))
      If (status/=0) Then ; err_msg = 'Error in my_nf90_get_var: NXX' ; Goto 999 ; End If
      !
      Case ('N34')
      status=nf90_get_var(ncid,iVar,strr%fn1(:,:,:,:,34)) ; in=in+1 ! write(*,*) 'N34 is found.'
      str%fn1(:,:,:,:,34) = dble(strr%fn1(:,:,:,:,34))
      If (status/=0) Then ; err_msg = 'Error in my_nf90_get_var: NXX' ; Goto 999 ; End If
      !
      Case ('N35')
      status=nf90_get_var(ncid,iVar,strr%fn1(:,:,:,:,35)) ; in=in+1 ! write(*,*) 'N35 is found.'
      str%fn1(:,:,:,:,35) = dble(strr%fn1(:,:,:,:,35))
      If (status/=0) Then ; err_msg = 'Error in my_nf90_get_var: NXX' ; Goto 999 ; End If
      !
      Case ('N36')
      status=nf90_get_var(ncid,iVar,strr%fn1(:,:,:,:,36)) ; in=in+1 ! write(*,*) 'N36 is found.'
      str%fn1(:,:,:,:,36) = dble(strr%fn1(:,:,:,:,36))
      If (status/=0) Then ; err_msg = 'Error in my_nf90_get_var: NXX' ; Goto 999 ; End If
      !
      Case ('NA')
      status=nf90_get_var(ncid,iVar,strr%NAERO) ; write(*,*) 'NA is found.'
      str%NAERO = dble( strr%NAERO)
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: QAERO'
        Goto 999
      End If
      !
      Case ('QR')
      status=nf90_get_var(ncid,iVar,strr%QRAIN) ; write(*,*) 'QR is found.'
      str%QRAIN = dble( strr%QRAIN)
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: QRAIN'
        Goto 999
      End If
      !
      Case ('QC')
      status=nf90_get_var(ncid,iVar,strr%QCLOUD) ; write(*,*) 'QC is found.'
      str%QCLOUD = dble( strr%QCLOUD)
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: QCLOUD'
        Goto 999
      End If
      !
      Case ('NC')
      status=nf90_get_var(ncid,iVar,strr%QNCLOUD) ; write(*,*) 'NC is found.'
      str%QNCLOUD = dble( strr%QNCLOUD)
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: QNCLOUD'
        Goto 999
      End If
      !!
      Case ('NR')
      status=nf90_get_var(ncid,iVar,strr%QNRAIN) ; write(*,*) 'NR is found.'
      str%QNRAIN = dble( strr%QNRAIN)
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: QNRAIN'
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
    write(*,*) 'Number of mass bins read=',im
    write(*,*) 'Number of number bins read=',in
    !
    call deallocate_wrf_rvar_mp70(strr)
    !
  return
  end subroutine ReadInpSAM_MP_PHYSICS_70
  !!-- added by oue
  !----------------------------------------------------------------
  !
  !!! Added by oue, for SAM morrison microphysics
  subroutine ReadInpSAM_MP_PHYSICS_75(InpFile,str,status)
  Use netcdf
  Use typeSizes
  Use wrf_rvar_mod
  Use wrf_var_mod
  Implicit None
  !
  Character(len=*),Intent(in)                              :: InpFile
  Type(wrf_var_mp75),Intent(InOut)                         :: str
  Integer,Intent(out)                                      :: status
  !
  Type(wrf_rvar_mp75)                                      :: strr
  integer                                                  :: nDims, nVars
  integer                                                  :: ncid, iVar
  Character(len=nf90_max_name), Allocatable , Dimension(:) :: var_names
  Character(len=nf90_max_name)                             :: name, err_msg
  Integer                     :: im,in
    !
    !- Massages 
    write(*,*) 'Selected microphysics: 75 (SAM Morrison microphysics)'
    write(*,*) 'Required hydrometeor parameters: QN, QP, QR, QI, QS, NC, NR, NI, NS, NG'
    write(*,*) 'Reading hydrometeor mixing ratios and number concentrations...'
    im=0 ; in=0 !counting the number of bins
    !
    strr%nt=str%nt
    strr%nx=str%nx
    strr%ny=str%ny
    strr%nz=str%nz
    !
    call allocate_wrf_rvar_mp75(strr)
    call initialize_wrf_rvar_mp75(strr)
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
      !!
      
      Case ('NA')
      status=nf90_get_var(ncid,iVar,strr%NAERO) ; write(*,*) 'NA is found.'
      str%NAERO = dble( strr%NAERO)
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: QAERO'
        Goto 999
      End If
      !
      Case ('QR')
      status=nf90_get_var(ncid,iVar,strr%QRAIN) ; write(*,*) 'QR is found.'
      str%QRAIN = dble( strr%QRAIN)
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: QRAIN'
        Goto 999
      End If
      !
      Case ('QI')
      status=nf90_get_var(ncid,iVar,strr%QICE) ; write(*,*) 'QI is found.'
      str%QICE = dble( strr%QICE)
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: QICE'
        Goto 999
      End If
      !
      Case ('QS')
      status=nf90_get_var(ncid,iVar,strr%QSNOW) ; write(*,*) 'QS is found.'
      str%QSNOW = dble( strr%QSNOW)
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: QSNOW'
        Goto 999
      End If
      !
      Case ('QN')
      status=nf90_get_var(ncid,iVar,strr%QNOPREP) ; write(*,*) 'QN is found.'
      str%QNOPREP = dble( strr%QNOPREP)
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: QNOPREP (Non-precipitating Condensate (Water+Ice))'
        Goto 999
      End If
      !
      Case ('QP')
      status=nf90_get_var(ncid,iVar,strr%QPREP) ; write(*,*) 'QP is found.'
      str%QPREP = dble( strr%QPREP)
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: QPREP (Precipitating Water (Rain+Snow))'
        Goto 999
      End If
      !
      Case ('NC')
      status=nf90_get_var(ncid,iVar,strr%QNCLOUD) ; write(*,*) 'NC is found.'
      str%QNCLOUD = dble( strr%QNCLOUD)
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: QNCLOUD'
        Goto 999
      End If
      !!
      Case ('NR')
      status=nf90_get_var(ncid,iVar,strr%QNRAIN) ; write(*,*) 'NR is found.'
      str%QNRAIN = dble( strr%QNRAIN)
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: QNRAIN'
        Goto 999
      End If
      !
      Case ('NI')
      status=nf90_get_var(ncid,iVar,strr%QNICE) ; write(*,*) 'NI is found.'
      str%QNICE = dble( strr%QNICE)
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: QNICE'
        Goto 999
      End If
      !
      Case ('NS')
      status=nf90_get_var(ncid,iVar,strr%QNSNOW) ; write(*,*) 'NS is found.'
      str%QNSNOW = dble( strr%QNSNOW)
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: QNSNOW'
        Goto 999
      End If
      !
      Case ('NG')
      status=nf90_get_var(ncid,iVar,strr%QNGRAUP) ; write(*,*) 'NG is found.'
      str%QNGRAUP = dble( strr%QNGRAUP)
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: QNGRAUP'
        Goto 999
      End If
      !
      End Select
      !
    End Do
    !
    Deallocate(var_names)
    !
    !- Get QC and QG
    str%QCLOUD = str%QNOPREP - str%QICE
    str%QGRAUP = str%QPREP - str%QRAIN - str%QSNOW
    
999 If (status.Ne.0) Then
      write(*,*)  err_msg
      Call Exit(1)
    Else
      status=nf90_close(ncid)
    Endif
    !
    call deallocate_wrf_rvar_mp75(strr)
    !
  return
  end subroutine ReadInpSAM_MP_PHYSICS_75
  !!-- added by oue
  !------------------------------------------------------------------------------------
  subroutine get_env_vars_sam(conf,wrf,env)
  Use wrf_var_mod
  Use crsim_mod
  Use phys_param_mod, ONLY: Rd, eps,p0,cp,grav
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
  
  Integer            :: ix,iy,iz,izmin,izmax,idz
    !-------------------------------------------------------------
    !
    ix1=1 ; ix2=env%nx
    iy1=1 ; iy2=env%ny
    iz1=1 ; iz2=env%nz
    !it = 1
    !           
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
    !If (wrf%nt>1) it=conf%it
    !
    ! ------------------------------------------------------------------------------  
    ! ------------------------------------------------------------------------------ 
    ! get meteor. fields from wrf vars    
    !    
    !kap=cp/Rd
    
    izmin=1 ; izmax=wrf%nz ; idz=1
    if(wrf%pb(1,1,1,1) < wrf%pb(1,1,wrf%nz,1)) then
      izmin=wrf%nz ; izmax=1 ; idz=-1
    endif
    !  
    !------------------------------------------------------------------------------------  
    !-------------------------------------------------------------------------------       
    !
    it=conf%it
    ! x, y , z        
    Allocate(xtrack(wrf%nx),ytrack(wrf%ny),gheight(wrf%nx,wrf%ny,wrf%nz),&
            Kw(wrf%nx,wrf%ny,wrf%nz),Ku(wrf%nx,wrf%ny,wrf%nz),Kv(wrf%nx,wrf%ny,wrf%nz),&
            uu(wrf%nx,wrf%ny,wrf%nz),vv(wrf%nx,wrf%ny,wrf%nz),ww(wrf%nx,wrf%ny,wrf%nz))
    !                
    !!----------------
    do iz=1,wrf%nz
      ww(:,:,iz)=wrf%w(:,:,iz,it)
      if(iz<wrf%nz) then
        Kw(:,:,iz)= (wrf%w(:,:,iz+1,it)-wrf%w(:,:,iz,it)) / (wrf%hgtm(iz+1)-wrf%hgtm(iz))
      else
        Kw(:,:,iz)= (wrf%w(:,:,iz,it)-wrf%w(:,:,iz-1,it)) / (wrf%hgtm(iz)-wrf%hgtm(iz-1))
      endif
      if(wrf%hgtm(1) > wrf%hgtm(wrf%nz)) Kw(:,:,iz)=Kw(:,:,iz)*(-1.d0) 
    enddo
    !
    !---------------- 
    do ix=1,wrf%nx
      uu(ix,:,:)=wrf%u(ix,:,:,it)
      if(ix<wrf%nx)then
        Ku(ix,:,:)=wrf%dx(it) * (wrf%u(ix+1,:,:,it)-wrf%u(ix,:,:,it))
      else
        Ku(ix,:,:)=wrf%dx(it) * (wrf%u(ix,:,:,it)-wrf%u(ix-1,:,:,it))
      endif
    enddo
    !
    do iy=1,wrf%ny
      vv(:,iy,:)=wrf%v(:,iy,:,it)
      if(iy<wrf%ny)then
        Kv(:,iy,:)=wrf%dy(it) * (wrf%v(:,iy+1,:,it)-wrf%v(:,iy,:,it))
      else
        Kv(:,iy,:)=wrf%dy(it) * (wrf%v(:,iy,:,it)-wrf%v(:,iy-1,:,it))
      endif
    enddo
    !
    !-----------------
    !
    env%x(1:env%nx)=wrf%xlong(ix1:ix2,1,1) - wrf%xlong(1,1,1) !correct x so that x(1)=0 modified by oue Apr 2020
    env%y(1:env%ny)=wrf%xlat(1,iy1:iy2,1) - wrf%xlat(1,1,1) !correct y so that y(1)=0 modified by oue Apr 2020
    do ix=1,env%nx
      do iy=1,env%ny
        env%z(ix,iy,1:env%nz)=wrf%hgtm(iz1:iz2) ! m  
        env%press(ix,iy,1:env%nz)=wrf%pb(1,1,iz1:iz2,1) * 100.d0 ! + wrf%p(ix,iy,iz,1) ! Pa 
      enddo
    enddo
    !
    env%u(1:env%nx,1:env%ny,1:env%nz)=wrf%u(ix1:ix2,iy1:iy2,iz1:iz2,conf%it)
    env%v(1:env%nx,1:env%ny,1:env%nz)=wrf%v(ix1:ix2,iy1:iy2,iz1:iz2,conf%it)
    env%w(1:env%nx,1:env%ny,1:env%nz)=wrf%w(ix1:ix2,iy1:iy2,iz1:iz2,conf%it)
    !
    !env%xlat(1:env%nx,1:env%ny)=wrf%xlat(ix1:ix2,iy1:iy2,it) ! deg   
    !env%xlong(1:env%nx,1:env%ny)=wrf%xlong(ix1:ix2,iy1:iy2,it) ! deg  
    !
    env%Kw(1:env%nx,1:env%ny,1:env%nz)=Kw(ix1:ix2,iy1:iy2,iz1:iz2)
    env%Ku(1:env%nx,1:env%ny,1:env%nz)=Ku(ix1:ix2,iy1:iy2,iz1:iz2)
    env%Kv(1:env%nx,1:env%ny,1:env%nz)=Kv(ix1:ix2,iy1:iy2,iz1:iz2)
    !
    ! added by oue Apr 2020
    env%dx = 1./wrf%dx(it)
    env%dy = 1./wrf%dy(it)
    !--
    
    !---------------------------------------------------------------------------------
    Deallocate(ww,uu,vv,Kw,Ku,Kv)
    !---------------------------------------------------------------------------------
    !--------------------------------------------------------------------------------- 
    !
    write(*,*) 'Info: Getting press,temp,rho_d,qvapor'
    !env%press(1:env%nx,1:env%ny,1:env%nz)=wrf%press(ix1:ix2,iy1:iy2,iz1:iz2,it)) ! Pa 
    env%temp(1:env%nx,1:env%ny,1:env%nz)=wrf%temp(ix1:ix2,iy1:iy2,iz1:iz2,it)   ! K 
    env%qvapor(1:env%nx,1:env%ny,1:env%nz)=wrf%qvapor(ix1:ix2,iy1:iy2,iz1:iz2,it)*1.d-3 ! g/kg -> kg/kg 
    env%tke(1:env%nx,1:env%ny,1:env%nz)=wrf%tke(ix1:ix2,iy1:iy2,iz1:iz2,it) ! m^2/s^2    
    
    if (MaxVal(wrf%rho_d(1,1,:,1)) > 0.d0) then
      do ix=1,env%nx
        do iy=1,env%ny
          env%rho_d(ix,iy,1:env%nz)=wrf%rho_d(1,1,iz1:iz2,1) ! kg/m^3  
        enddo
      enddo
    else
      write(*,*) 'Computing RHO from Input data'
      env%rho_d=(env%press*env%qvapor)/(eps+env%qvapor) ! e water vapor pressure in Pa
      env%rho_d=(env%press-env%rho_d)/(Rd * env%temp)   !  kg/m^3
    endif    
    !
    ! convert T from K to C 
    env%temp=env%temp-273.15d0 ! C 
    ! convert press frpm Pa to mb   
    env%press=env%press*1.d-2  ! mb 
    
    write(*,*) 'press [mb]',MinVal(env%press),MaxVal(env%press)
    write(*,*) 'temp [C]',MinVal(env%temp),MaxVal(env%temp)
    write(*,*) 'rho_d [kg/m^3]',MinVal(env%rho_d),MaxVal(env%rho_d)
    !write(*,*) 'tke [m^2/s^2]',MinVal(env%tke),MaxVal(env%tke)
    ! 
  return
  end subroutine get_env_vars_sam
  !!
  !----------------------------------------------------------------
  !----------------------------------------------------------------
  ! Read CM1 output added by oue Apr 2020
  !---------------------------------------------------------------------
  subroutine ReadInpCM1_dim(InpFile,str,status)
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
      Case ('ni')        ;  str%nx    = length
      Case ('nj')        ;  str%ny    = length
      Case ('nk')       ;  str%nz    = length
      Case ('time')       ;  str%nt    = length
      End Select
      !
    enddo   ! iDim
    !
    !str%nt    = 1
    str%nxp1  = str%nx
    str%nyp1  = str%ny
    str%nzp1  = str%nz
    
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
    If (str%nxp1/=str%nx) then
       write(*,*) 'Problem in input dimensions nx and nxp1'
       write(*,*) 'nx=nxp1 for CM1 input'
       Call Exit(1)
    EndIf
    !
    If (str%nyp1/=str%ny) then
      write(*,*) 'Problem in input dimensions ny and nyp1'
      write(*,*) 'nx=ny, ny=nyp1 for CM1 input'
      Call Exit(1)
    EndIf
    !
  return
  end subroutine ReadInpCM1_dim
  !
  !---------------------------------------------------------------------
  subroutine ReadInpCM1_var(InpFile,str,status)
  Use netcdf
  Use wrf_var_mod
  Use wrf_rvar_mod
  Implicit None
  !
  Character(len=*),Intent(in)                              :: InpFile
  Type(wrf_var),Intent(InOut)                              :: str
  Integer, Intent(out)                                     :: status
  !
  Type(wrf_rvar)                                           :: strr
  Integer                                                  :: ncid
  Integer                                                  :: nDims,nVars,iVar
  Character(len=nf90_max_name)                             :: name, err_msg
  Character(len=nf90_max_name), Allocatable , Dimension(:) :: var_names
  
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
    call allocate_wrf_rvar(strr)
    call initialize_wrf_rvar(strr)
    ! 
    status=  nf90_open(Trim(InpFile), NF90_NOWRITE, ncid)
    if (status.ne.0) then ; err_msg='Error in nf90_open' ; goto 999 ; endif
    !
    status= nf90_inquire(ncid, nDims, nVars)
    If (status/=0) Then ; err_msg = 'Error in nf90_inquire' ; goto 999 ; endif
    !
    !
    strr%p=0.e0
    str%p=0.d0
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
      !
      Case ('uinterp')
      status=nf90_get_var(ncid,iVar,strr%u)
      str%u = dble (strr%u)
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: uinterp'
        Goto 999
      End If
      !
      Case ('vinterp')
      status=nf90_get_var(ncid,iVar,strr%v)
      str%v = dble (strr%v)
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: vinterp'
        Goto 999
      End If
      !
      Case ('winterp')
      status=nf90_get_var(ncid,iVar,strr%w)
      str%w = dble (strr%w)
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: winterp'
        Goto 999
      End If
      !
      Case ('prs')
      status=nf90_get_var(ncid,iVar,strr%press)
      str%press = dble (strr%press)
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: prs'
        Goto 999
      End If
      !
      !Case ('PP')
      !status=nf90_get_var(ncid,iVar,strr%p)
      !str%p = dble (strr%p)
      !If (status/=0) Then
      !  err_msg = 'Error in my_nf90_get_var: PP'
      !  Goto 999
      !End If
      !
      Case ('th')
      status=nf90_get_var(ncid,iVar,strr%theta)
      str%theta = dble (strr%theta)
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: th (potential temperature)'
        Goto 999
      End If
      !
      Case ('z')
      status=nf90_get_var(ncid,iVar,strr%hgtm)
      str%hgtm = dble (strr%hgtm) * 1.d3
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: z'
        Goto 999
      End If
      !
      Case ('qv')
      status=nf90_get_var(ncid,iVar,strr%QVAPOR)
      str%QVAPOR = dble (strr%QVAPOR)
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: QVAPOR'
        Goto 999
      End If
      !
      Case ('xh')
      status=nf90_get_var(ncid,iVar,strr%xlong(:,1,1))
      str%xlong(:,1,1) = dble (strr%xlong(:,1,1)) * 1.d3
      str%xlong(:,1,1) = str%xlong(:,1,1) 
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: x'
        Goto 999
      End If
      !
      Case ('yh')
      status=nf90_get_var(ncid,iVar,strr%xlat(1,:,1))
      str%xlat(1,:,1) = dble (strr%xlat(1,:,1)) * 1.d3
      str%xlat(1,:,1) = str%xlat(1,:,1) 
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: y'
        Goto 999
      End If
      !
      !
      End Select
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
    
    !------------------------------------
    !get dx/dy
    str%dx=1./(str%xlong(2,1,1)-str%xlong(1,1,1))
    str%dy=1./(str%xlat(1,2,1)-str%xlat(1,1,1))
    write(*,*) 'CM1 dx,dy',1./str%dx(1),1./str%dy(1)
    !------------------------------------
    !
    !
    call deallocate_wrf_rvar(strr)
    !
  return
  end subroutine ReadInpCM1_var
  !
  !----------------------------------------------------------------
  !
  !!! Added by oue, for CM1 morrison microphysics
  subroutine ReadInpCM1_MP_PHYSICS_80(InpFile,str,status)
  Use netcdf
  Use typeSizes
  Use wrf_rvar_mod
  Use wrf_var_mod
  Implicit None
  !
  Character(len=*),Intent(in)                              :: InpFile
  Type(wrf_var_mp80),Intent(InOut)                         :: str
  Integer,Intent(out)                                      :: status
  !
  Type(wrf_rvar_mp80)                                      :: strr
  integer                                                  :: nDims, nVars
  integer                                                  :: ncid, iVar
  Character(len=nf90_max_name), Allocatable , Dimension(:) :: var_names
  Character(len=nf90_max_name)                             :: name, err_msg
  Integer                     :: im,in
    !
    !- Massages 
    write(*,*) 'Selected microphysics: 80 (CM1 Morrison microphysics)'
    write(*,*) 'Required hydrometeor parameters: QC, QR, QI, QS, QG, NR, NI, NS, NG'
    write(*,*) 'Reading hydrometeor mixing ratios and number concentrations...'
    im=0 ; in=0 !counting the number of bins
    !
    strr%nt=str%nt
    strr%nx=str%nx
    strr%ny=str%ny
    strr%nz=str%nz
    !
    call allocate_wrf_rvar_mp80(strr)
    call initialize_wrf_rvar_mp80(strr)
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
      !!
      
      !
      Case ('qc')
      status=nf90_get_var(ncid,iVar,strr%QCLOUD) ; write(*,*) 'qc is found.'
      str%QCLOUD = dble( strr%QCLOUD)
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: QCLOUD'
        Goto 999
      End If
      !
      Case ('qr')
      status=nf90_get_var(ncid,iVar,strr%QRAIN) ; write(*,*) 'qr is found.'
      str%QRAIN = dble( strr%QRAIN)
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: QRAIN'
        Goto 999
      End If
      !
      Case ('qi')
      status=nf90_get_var(ncid,iVar,strr%QICE) ; write(*,*) 'qi is found.'
      str%QICE = dble( strr%QICE)
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: QICE'
        Goto 999
      End If
      !
      Case ('qs')
      status=nf90_get_var(ncid,iVar,strr%QSNOW) ; write(*,*) 'qs is found.'
      str%QSNOW = dble( strr%QSNOW)
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: QSNOW'
        Goto 999
      End If
      !
      Case ('qg')
      status=nf90_get_var(ncid,iVar,strr%QGRAUP) ; write(*,*) 'qg is found.'
      str%QGRAUP = dble( strr%QGRAUP)
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: QGRAUPEL'
        Goto 999
      End If
      !
      Case ('ncc')
      status=nf90_get_var(ncid,iVar,strr%QNCLOUD) ; write(*,*) 'ncc is found.'
      str%QNCLOUD = dble( strr%QNCLOUD)
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: QNCLOUD'
        Goto 999
      End If
      !!
      Case ('ncr')
      status=nf90_get_var(ncid,iVar,strr%QNRAIN) ; write(*,*) 'ncr is found.'
      str%QNRAIN = dble( strr%QNRAIN)
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: QNRAIN'
        Goto 999
      End If
      !
      Case ('nci')
      status=nf90_get_var(ncid,iVar,strr%QNICE) ; write(*,*) 'nci is found.'
      str%QNICE = dble( strr%QNICE)
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: QNICE'
        Goto 999
      End If
      !
      Case ('ncs')
      status=nf90_get_var(ncid,iVar,strr%QNSNOW) ; write(*,*) 'ncs is found.'
      str%QNSNOW = dble( strr%QNSNOW)
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: QNSNOW'
        Goto 999
      End If
      !
      Case ('ncg')
      status=nf90_get_var(ncid,iVar,strr%QNGRAUP) ; write(*,*) 'ncg is found.'
      str%QNGRAUP = dble( strr%QNGRAUP)
      If (status/=0) Then
        err_msg = 'Error in my_nf90_get_var: QNGRAUP'
        Goto 999
      End If
      !
      End Select
      !
    End Do
    !
    Deallocate(var_names)
    !
    
999 If (status.Ne.0) Then
      write(*,*)  err_msg
      Call Exit(1)
    Else
      status=nf90_close(ncid)
    Endif
    !
    call deallocate_wrf_rvar_mp80(strr)
    !
  return
  end subroutine ReadInpCM1_MP_PHYSICS_80
  !!-- added by oue Apr 2020
  !------------------------------------------------------------------------------------
  subroutine get_env_vars_cm1(conf,wrf,env)
  Use wrf_var_mod
  Use crsim_mod
  Use phys_param_mod, ONLY: Rd, eps,p0,cp,grav
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
  
  Integer            :: ix,iy,iz,izmin,izmax,idz
    !-------------------------------------------------------------
    !
    ix1=1 ; ix2=env%nx
    iy1=1 ; iy2=env%ny
    iz1=1 ; iz2=env%nz
    !it = 1
    !           
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
    !If (wrf%nt>1) it=conf%it
    !
    ! ------------------------------------------------------------------------------  
    ! ------------------------------------------------------------------------------ 
    ! get meteor. fields from wrf vars    
    !    
    !kap=cp/Rd
    
    izmin=1 ; izmax=wrf%nz ; idz=1
    if(wrf%press(1,1,1,1) < wrf%press(1,1,wrf%nz,1)) then
      izmin=wrf%nz ; izmax=1 ; idz=-1
    endif
    !  
    !------------------------------------------------------------------------------------
    !-- temperature from theta
    wrf%temp=wrf%theta*(wrf%press/100000.d0)**(287.d0/1004.d0) ! K 
    !-------------------------------------------------------------------------------       
    !
    it=conf%it
    ! x, y , z        
    Allocate(xtrack(wrf%nx),ytrack(wrf%ny),gheight(wrf%nx,wrf%ny,wrf%nz),&
            Kw(wrf%nx,wrf%ny,wrf%nz),Ku(wrf%nx,wrf%ny,wrf%nz),Kv(wrf%nx,wrf%ny,wrf%nz),&
            uu(wrf%nx,wrf%ny,wrf%nz),vv(wrf%nx,wrf%ny,wrf%nz),ww(wrf%nx,wrf%ny,wrf%nz))
    !                
    !!----------------
    do iz=1,wrf%nz
      ww(:,:,iz)=wrf%w(:,:,iz,it)
      if(iz<wrf%nz) then
        Kw(:,:,iz)= (wrf%w(:,:,iz+1,it)-wrf%w(:,:,iz,it)) / (wrf%hgtm(iz+1)-wrf%hgtm(iz))
      else
        Kw(:,:,iz)= (wrf%w(:,:,iz,it)-wrf%w(:,:,iz-1,it)) / (wrf%hgtm(iz)-wrf%hgtm(iz-1))
      endif
      if(wrf%hgtm(1) > wrf%hgtm(wrf%nz)) Kw(:,:,iz)=Kw(:,:,iz)*(-1.d0) 
    enddo
    !
    !---------------- 
    do ix=1,wrf%nx
      uu(ix,:,:)=wrf%u(ix,:,:,it)
      if(ix<wrf%nx)then
        Ku(ix,:,:)=wrf%dx(it) * (wrf%u(ix+1,:,:,it)-wrf%u(ix,:,:,it))
      else
        Ku(ix,:,:)=wrf%dx(it) * (wrf%u(ix,:,:,it)-wrf%u(ix-1,:,:,it))
      endif
    enddo
    !
    do iy=1,wrf%ny
      vv(:,iy,:)=wrf%v(:,iy,:,it)
      if(iy<wrf%ny)then
        Kv(:,iy,:)=wrf%dy(it) * (wrf%v(:,iy+1,:,it)-wrf%v(:,iy,:,it))
      else
        Kv(:,iy,:)=wrf%dy(it) * (wrf%v(:,iy,:,it)-wrf%v(:,iy-1,:,it))
      endif
    enddo
    !
    !-----------------
    !
    env%x(1:env%nx)=wrf%xlong(ix1:ix2,1,1) - wrf%xlong(1,1,1) !correct x so that x(1)=0
    env%y(1:env%ny)=wrf%xlat(1,iy1:iy2,1) - wrf%xlat(1,1,1) !correct y so that y(1)=0
    do ix=1,env%nx
      do iy=1,env%ny
        env%z(ix,iy,1:env%nz)=wrf%hgtm(iz1:iz2) ! m  
      enddo
    enddo
    !
    env%u(1:env%nx,1:env%ny,1:env%nz)=wrf%u(ix1:ix2,iy1:iy2,iz1:iz2,conf%it)
    env%v(1:env%nx,1:env%ny,1:env%nz)=wrf%v(ix1:ix2,iy1:iy2,iz1:iz2,conf%it)
    env%w(1:env%nx,1:env%ny,1:env%nz)=wrf%w(ix1:ix2,iy1:iy2,iz1:iz2,conf%it)
    !
    !env%xlat(1:env%nx,1:env%ny)=wrf%xlat(ix1:ix2,iy1:iy2,it) ! deg   
    !env%xlong(1:env%nx,1:env%ny)=wrf%xlong(ix1:ix2,iy1:iy2,it) ! deg  
    !
    env%Kw(1:env%nx,1:env%ny,1:env%nz)=Kw(ix1:ix2,iy1:iy2,iz1:iz2)
    env%Ku(1:env%nx,1:env%ny,1:env%nz)=Ku(ix1:ix2,iy1:iy2,iz1:iz2)
    env%Kv(1:env%nx,1:env%ny,1:env%nz)=Kv(ix1:ix2,iy1:iy2,iz1:iz2)
    !
    env%dx = 1./wrf%dx(it)
    env%dy = 1./wrf%dy(it)
    !
    
    !---------------------------------------------------------------------------------
    Deallocate(ww,uu,vv,Kw,Ku,Kv)
    !---------------------------------------------------------------------------------
    !--------------------------------------------------------------------------------- 
    !
    write(*,*) 'Info: Getting press,temp,rho_d,qvapor'
    env%press(1:env%nx,1:env%ny,1:env%nz)=wrf%press(ix1:ix2,iy1:iy2,iz1:iz2,it) ! Pa 
    env%temp(1:env%nx,1:env%ny,1:env%nz)=wrf%temp(ix1:ix2,iy1:iy2,iz1:iz2,it)   ! K 
    env%qvapor(1:env%nx,1:env%ny,1:env%nz)=wrf%qvapor(ix1:ix2,iy1:iy2,iz1:iz2,it) ! kg/kg 
    env%tke(1:env%nx,1:env%ny,1:env%nz)=wrf%tke(ix1:ix2,iy1:iy2,iz1:iz2,it) ! m^2/s^2    
    
    write(*,*) 'Computing RHO from Input data'
    env%rho_d=(env%press*env%qvapor)/(eps+env%qvapor) ! e water vapor pressure in Pa
    env%rho_d=(env%press-env%rho_d)/(Rd * env%temp)   !  kg/m^3
    !
    ! convert T from K to C 
    env%temp=env%temp-273.15d0 ! C 
    ! convert press frpm Pa to mb   
    env%press=env%press*1.d-2  ! mb 
    
    write(*,*) 'press [mb]',MinVal(env%press),MaxVal(env%press)
    write(*,*) 'temp [C]',MinVal(env%temp),MaxVal(env%temp)
    write(*,*) 'rho_d [kg/m^3]',MinVal(env%rho_d),MaxVal(env%rho_d)
    !write(*,*) 'tke [m^2/s^2]',MinVal(env%tke),MaxVal(env%tke)
    ! 
  return
  end subroutine get_env_vars_cm1
  !!
  !----------------------------------------------------------------
