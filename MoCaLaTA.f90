!================================  Modules  ===================================

module kinds
  integer, parameter:: i4b     = selected_int_kind(9)
  integer, parameter:: sp      = kind(1.0)
  integer, parameter:: dp      = kind(1.0d0)
  integer, parameter:: qp      = 2 * dp
  integer, parameter:: posKind = kind(1d0)   !Double is probably necessary
end module kinds

!------------------------------------------------------------------------------

module PhysicalConstants
  use kinds
  implicit none
  !Astrophysical parameters
  real(dp)::           H_0                        !Hubble constant
  real(dp)::           Omega_M                    !Matter density parameter
  real(dp)::           Omega_L                    !Dark energy density parameter
  real(dp),parameter:: kpc      = 3.08567758131d21!Kiloparsec; cm
  real(dp),parameter:: pc       = 3.08567758131d18!Parsec; cm
  real(dp),parameter:: c        = 2.99792458d10   !Speed of light; cm/s
  real(dp),parameter:: M_sun    = 1.98892d33      !Solar mass; g

  !Atomic parameters
  real(dp),parameter:: m_H      = 1.673534d-24    !Hydrogen mass; g
  real(dp),parameter:: m_e      = 9.1093897d-28   !Electron mass; g
  real(dp),parameter:: e        = 4.803206d-10    !Electron charge; esu
  real(dp),parameter:: h_Pl     = 6.6260755d-27   !Planck's constant; erg s
  real(dp),parameter:: k_B      = 1.380658d-16    !Boltzmann's constant; erg/K
  real(dp),parameter:: eV       = 1.60217653d-12  !Electron volt; erg

  !Lya parameters
  real(dp),parameter:: f_12     = 0.4162          !Oscillator strength
  real(dp),parameter:: nu_0     = 2.46606d15      !Lya line-center frequency; Hz
  real(dp),parameter:: lambda_0 = c / nu_0 * 1d8  !Line center wavelength; Angstrom
  real(dp),parameter:: E_Lya    = h_Pl * nu_0     !Lya energy; erg
  real(dp),parameter:: Dnu_L    = 9.936d7         !Lya natural line width; Hz
  real(dp),parameter:: eta      = .71             !Slab-to-cube converter
  real(dp),parameter:: at_lim   = 2000.           !Value of a * tau_0 beyond which to apply Neufeld solution

  !Dust parameters:
  real(dp)::           albedo                     !Pscat / (Pabs + Pscat) on dust
  real(dp)::           g                          !Dust scattering asymmetry parameter
  real(dp)::           ta2EBV                     !tau_a-to-E(B-V) conversion parameter
  real(dp)::           f_ion                      !Fraction of ionized hydrogen that contributes to dust density
  real(dp),parameter:: zeta     = 0.525           !Neufeld f_esc fitting parameter

  !Numerical constants:
  real(dp),parameter:: pi       = acos(-1d0)                                      !3.14159265358979!pi
  real(dp),parameter:: pib2     = pi / 2                                          !1.57079632679490!pi/2
  real(dp),parameter:: tpi      = 2 * pi                                          !6.28318530717959!2pi
  real(dp),parameter:: sqpi     = sqrt(pi)                                        !1.77245385090552!sqrt(pi)
  real(dp),parameter:: sd2FW    = 2 * sqrt(2 * log(2d0))                          !2.354820046
end module PhysicalConstants

!------------------------------------------------------------------------------

module CellStructure
  use kinds
  implicit none
  type:: Cell
    real(sp)::            n_HI, Dnu_D, U_bulk(3)
    real(kind=posKind)::  C_pos(3)
    type(Cell), pointer:: Child(:,:,:)
#   ifdef AMR
    integer(kind=1)::     Level                 !Max. 127 levels of refinement
    logical(kind=1)::     Refined
#   endif
#   ifdef dust
    real(sp)::            n_d
#   endif
#   ifdef multi
    integer::             clnum
    integer(kind=1)::     phase !0 = "void", ie vacuum,
                                !1 = "cloud", 
                                !2 = "ICM"
                                !3 = "shell"
#   endif
  end type Cell
end module CellStructure

!------------------------------------------------------------------------------

module DataArrays
  use kinds
  use CellStructure
  real(sp),pointer::        n_HIString(:),Dnu_DString(:), &
                            U_xString(:),U_yString(:),U_zString(:), &
                            InitLya(:,:),InitFUV(:,:)
# ifdef AMR
  integer(kind=1),pointer:: LevelString(:)
# endif
# ifdef multi
  integer(kind=1),pointer:: PhaseString(:)
  integer,pointer::         CloudString(:)
# endif
# ifdef dust
  real(sp),pointer::        n_dString(:),n_HIIString(:)
# endif
  type(Cell),target::       BaseGrid            !Mother grid with resolution ni,nj,nk:
  real(dp),allocatable::    MoreClouds(:,:),CloudVelDisp(:,:)
  real(sp),allocatable::    Clouds(:,:)
end module DataArrays

!------------------------------------------------------------------------------

module GlobalData
  use kinds
  integer::             ni,nj,nk                !Base grid resolution
  integer::             wrtph,wrtfl
  integer::             BaseRes(3)
  real(dp)::            dx0,dy0,dz0,dxTiny,dxWeeny!Base grid cell size / negligible distance
  real(dp)::            z_0                     !Half thickness of slab/cube/sphere
  real(dp),allocatable,&
         dimension(:):: dx,dxH,dy,dyH,dz,dzH    !Cell size / Half cell size
  integer::             N_cells,i_cell          !Total # of cells; counter
  integer,parameter::   n_phTypes = 2           !For now: Lya + FUV = 2 types of photons
  integer::             n_phtot,n_ph(n_phTypes) !Total # of photon packets to be launched/...per type
  integer::             n_packets(n_phTypes)    !Packets launched so far per photon type
  real(dp)::            n_eff(n_phTypes)        !Effective number of photons launched so far
  real(dp)::            n_FUV,NBaper,n_absFUV   !FUV photons that are sufficiently far from the line center not to scatter (significantly)
  logical::             trueFUV                 !Flag if outside NBaper
  real(dp)::            n_abs(n_phTypes)        !Number of absorbed photons
  integer::             emdimL,emdimF           !
  integer::             n_phases                !# of different phases in multiphase medium, excluding "void"
  real(dp),allocatable::FFcounter(:)            !Counter for cell phases
  integer::             inlun,outlun,cloudlun,windlun,N0lun!Logical unit numbers
  integer::             L_max                   !Maximum level of refinement
  real(dp)::            D_box,R_box             !Phys. size/half size of Box
  real(dp)::            r_gal,axes(3)           !Radius of model galaxy, ellipsoid axes in terms of r_gal
  real(dp)::            r_inner,r_outer         !Inner and outer radius of shell
  real(dp)::            H1,H2,exp1,exp2,H_lum(3)!Scale height and exp. factor of model galaxy along and perp. to plane (dens. propto exp((r/H1,2)^exp1,2))
  real(dp)::            CloudCorr               !Emission cloud correlation factor ([0,1], where 0 is "no emission from clouds" and 1 is "no emission from ICM")
  real(dp)::            OA,cosOAb2sq,n_jet(3)   !Wind openening angle and direction
  real(dp)::            X_QSO(3)                !Background QSO position
  character(len=50)::   X_QSOc
  real(dp)::            Q_Lya,Q_FUV,Q_FUV_AA,Q_tot !Total # of photons (line; cont. in BW; cont. per Angstrom, total in BW)
  real(dp)::            L_Lya,L_FUV_Hz,L_FUV,L_tot
  real(dp)::            nu1,nu2
  integer::             N_cl
  real(dp)::            FF                      !Cloud filling factor
  real(dp)::            rf,r_cl,rmin,rmax,beta
  real(dp)::            n_HI_cl, T_cl, Dnu_D_cl, n_d_cl, Z_cl
  real(dp)::            n_HI_ICM,T_ICM,Dnu_D_ICM,n_d_ICM,Z_ICM
  real(dp)::            n_HI_sh, T_sh, Dnu_D_sh, n_d_sh, Z_sh
  real(dp)::            V_in,V_out,V_inf,sigV_cl,V_z,Z0,EBV
  logical::             intxc,maxxc,recoil,isoem,empos(5),followGrad,writeN0,nullifyx,ArchitectureOfAggression,fullboxRT
  integer::             N_los                   !# of sightlines for calculating average quantities
  character(len=200)::  X_init
  real(dp)::            xc_max,nhat_0(3),X_i0(3)
  real(dp)::            EW_int,Lya2tot          !Initrinsic equivalent width; line-to-[line+cont.] ratio (depends on BW)
  real(dp)::            z,d_L                   !Redshift and luminosity distance of system
  real(dp)::            t1,t2                   !Time
  character(len=500)::  FullPathString
  integer,allocatable:: FullPathArray(:,:)
  character(len=200)::  DataDir,Subdir,CellData !Data directories
  character(len=200)::  Cloudfile,Windfile      !File for writing cloud positions
  character(len=200)::  parfile                 !Input parameter file
  character(len=50)::   model                   !
  character(len=20)::   DustType,x_injType,x_critType,Vprof
  real(dp)::            userXsecD,x_inj,sigma_V
  real(dp)::            atA(10),xcA(20),fescxc(10,20)!Tables for core escape probabilities
  real(dp)::            dat,dxc                 !Steps in above tables
  real(dp)::            Clumpiness              !Clumpiness factor
end module GlobalData

!------------------------------------------------------------------------------

module IJK
  integer:: i,j,k
end module IJK

!------------------------------------------------------------------------------

module ObservationalParameters
  use kinds
  implicit none
  real(dp), &
  allocatable, &
  dimension(:,:,:)::   CCDxp,CCDyp,CCDzp,CCDxm,CCDym,CCDzm
  logical::            WhichCCDs(6)
  real(dp), &
  allocatable, &
  dimension(:)::       specxp,specyp,speczp,specxm,specym,speczm
  character(len=200):: specin,specout,outlabel,ext

  real(dp)::           BW(2),dlam,DDlam         !All in Angstrom
  real(dp)::           w_esc,w_ph               !Photon weights
  real(dp),parameter:: tau_max = 30.            !tau threshold before abandoning RT from point of scattering to box edge
  integer::            pix,SpecRes2D,SpecRes1D  !CCD pixels; spectral resolution
  real(dp)::           D_obs,R_obs,o1,o2        !Area covered by CCD
  real(dp)::           obspar(4)                !(<NHI>,<EBV>,NHI_QSO,EBV_QSO)
end module ObservationalParameters

!------------------------------------------------------------------------------

module CurrentCellParameters
  use kinds
  real(dp):: C_pos(3)                           !Position of cell center
  real(dp):: n_HI                               !Density of neutral hydrogen
  real(dp):: n_d                                !"Dust density", i.e. rescaled hydrogen density
  real(dp):: Dnu_D                              !Doppler width of line
  real(dp):: U_bulk(3)                          !Gas bulk velocity in Doppler widths
# ifdef multi
  integer::  clnum                              !Cloud number in array "Clouds"
# endif
  real(dp):: tau_0                              !Optical depth of cell from center to edge for a x = 0 photon
  real(dp):: x_crit                             !Critical x-value for determining atom velocity
  real(dp):: xp,xm,yp,ym,zp,zm,Faces(6)
end module CurrentCellParameters

!------------------------------------------------------------------------------

module DerivedParameters
  use kinds
  real(dp)::  a,sigma_x,sigma_d,Faces_eff(6),R_cube,z_frac,x_adj,logatr
  real(dp)::  at,tau_a,at_eff                     !From CENTER to face
  integer:: Adjac
  real(dp)::  g_rec
  real(dp)::  x_cw
end module DerivedParameters

!------------------------------------------------------------------------------

module AtomParameters
  use kinds
  real(dp)::              u_II,u_1,u_2            !Atom velocities paral./perp. to incident photon
  real(dp),dimension(3):: U,uhat_1,uhat_2         !Atom vel. and unit vectors of perp. velocities in lab system
end module AtomParameters

!------------------------------------------------------------------------------

module iArrayCopy

  interface ArrayCopy
    module procedure ArrayCopy_r, ArrayCopy_d, ArrayCopy_i
  end interface

  contains
    subroutine ArrayCopy_r(src,dest,n_copied,n_not_copied)
    use kinds
    real(sp), dimension(:), intent(in)::  src
    real(sp), dimension(:), intent(out):: dest
    integer(i4b), intent(out)::           n_copied, n_not_copied
    n_copied         = min(size(src),size(dest))
    n_not_copied     = size(src) - n_copied
    dest(1:n_copied) = src(1:n_copied)
    end subroutine ArrayCopy_r

    subroutine ArrayCopy_d(src,dest,n_copied,n_not_copied)
    use kinds
    real(dp), dimension(:), intent(in)::  src
    real(dp), dimension(:), intent(out):: dest
    integer(i4b), intent(out)::           n_copied, n_not_copied
    n_copied         = min(size(src),size(dest))
    n_not_copied     = size(src) - n_copied
    dest(1:n_copied) = src(1:n_copied)
    end subroutine ArrayCopy_d

    subroutine ArrayCopy_i(src,dest,n_copied,n_not_copied)
    use kinds
    integer(i4b), dimension(:), intent(in)::  src
    integer(i4b), dimension(:), intent(out):: dest
    integer(i4b), intent(out)::               n_copied, n_not_copied
    n_copied         = min(size(src),size(dest))
    n_not_copied     = size(src) - n_copied
    dest(1:n_copied) = src(1:n_copied)
    end subroutine ArrayCopy_i
end module iArrayCopy

!------------------------------------------------------------------------------

module iAssignParameters
  interface
    subroutine AssignParameters(CurrentCell)
      use kinds
      use CellStructure
      implicit none
      type(Cell), target::  CurrentCell
    end subroutine AssignParameters
  end interface
end module iAssignParameters

!------------------------------------------------------------------------------

module iv2u
  interface v2u
    function v2u_1s2s(v,Dnu_D)
      use kinds
      implicit none
      real(sp),intent(in):: v,Dnu_D
      real(sp)::            v2u_1s2s
    end function v2u_1s2s
    function v2u_1d2d(v,Dnu_D)
      use kinds
      implicit none
      real(dp),intent(in):: v,Dnu_D
      real(dp)::            v2u_1d2d
    end function v2u_1d2d
    function v2u_1sv2sv(v,Dnu_D)
      use kinds
      implicit none
      real(sp),intent(in):: v(:),Dnu_D(:)
      real(sp)::            v2u_1sv2sv(size(v))
    end function v2u_1sv2sv
    function v2u_1dv2dv(v,Dnu_D)
      use kinds
      implicit none
      real(dp),intent(in):: v(:),Dnu_D(:)
      real(dp)::            v2u_1dv2dv(size(v))
    end function v2u_1dv2dv
    function v2u_1s2sv(v,Dnu_D)
      use kinds
      implicit none
      real(sp),intent(in):: v,Dnu_D(:)
      real(sp)::            v2u_1s2sv(size(Dnu_D))
    end function v2u_1s2sv
    function v2u_1d2dv(v,Dnu_D)
      use kinds
      implicit none
      real(dp),intent(in):: v,Dnu_D(:)
      real(dp)::            v2u_1d2dv(size(Dnu_D))
    end function v2u_1d2dv
    function v2u_1sv2s(v,Dnu_D)
      use kinds
      implicit none
      real(sp),intent(in):: v(:),Dnu_D
      real(sp)::            v2u_1sv2s(size(v))
    end function v2u_1sv2s
    function v2u_1dv2d(v,Dnu_D)
      use kinds
      implicit none
      real(dp),intent(in):: v(:),Dnu_D
      real(dp)::            v2u_1dv2d(size(v))
    end function v2u_1dv2d
    function v2u_1s2d(v,Dnu_D)
      use kinds
      implicit none
      real(sp),intent(in):: v
      real(dp),intent(in):: Dnu_D
      real(dp)::            v2u_1s2d
    end function v2u_1s2d
    function v2u_1d2s(v,Dnu_D)
      use kinds
      implicit none
      real(dp),intent(in):: v
      real(sp),intent(in):: Dnu_D
      real(dp)::            v2u_1d2s
    end function v2u_1d2s
    function v2u_1sv2d(v,Dnu_D)
      use kinds
      implicit none
      real(sp),intent(in):: v(:)
      real(dp),intent(in):: Dnu_D
      real(dp)::            v2u_1sv2d(size(v))
    end function v2u_1sv2d
    function v2u_1dv2s(v,Dnu_D)
      use kinds
      implicit none
      real(dp),intent(in):: v(:)
      real(sp),intent(in):: Dnu_D
      real(dp)::            v2u_1dv2s(size(v))
    end function v2u_1dv2s
    function v2u_1s2dv(v,Dnu_D)
      use kinds
      implicit none
      real(sp),intent(in):: v
      real(dp),intent(in):: Dnu_D(:)
      real(dp)::            v2u_1s2dv(size(Dnu_D))
    end function v2u_1s2dv
    function v2u_1d2sv(v,Dnu_D)
      use kinds
      implicit none
      real(dp),intent(in):: v
      real(sp),intent(in):: Dnu_D(:)
      real(dp)::            v2u_1d2sv(size(Dnu_D))
    end function v2u_1d2sv
    function v2u_1sv2dv(v,Dnu_D)
      use kinds
      implicit none
      real(sp),intent(in):: v(:)
      real(dp),intent(in):: Dnu_D(:)
      real(dp)::            v2u_1sv2dv(size(v))
    end function v2u_1sv2dv
    function v2u_1dv2sv(v,Dnu_D)
      use kinds
      implicit none
      real(dp),intent(in):: v(:)
      real(sp),intent(in):: Dnu_D(:)
      real(dp)::            v2u_1dv2sv(size(v))
    end function v2u_1dv2sv
  end interface
end module iv2u

!------------------------------------------------------------------------------

module igasdev
  interface gasdev
    subroutine gasdev_s(harvest)
    use kinds
    real(sp), intent(out):: harvest
    end subroutine gasdev_s

    subroutine gasdev_sv(harvest)
    use kinds
    real(sp), dimension(:), intent(out):: harvest
    end subroutine gasdev_sv

    subroutine gasdev_d(harvest)
    use kinds
    real(dp), intent(out):: harvest
    end subroutine gasdev_d

    subroutine gasdev_dv(harvest)
    use kinds
    real(dp), dimension(:), intent(out):: harvest
    end subroutine gasdev_dv
  end interface
end module igasdev

!------------------------------------------------------------------------------

module iLocateHostCell
  interface
    subroutine LocateHostCell(X_pos,CurrentCell,i,j,k,escape)
      use kinds
      use CellStructure
      implicit none
      real(dp), intent(in)::   X_pos(3)
      type(Cell), pointer::  CurrentCell        !intent(OUT)
      integer, intent(out):: i,j,k
      logical,intent(out)::  escape
    end subroutine LocateHostCell
  end interface
end module iLocateHostCell

!------------------------------------------------------------------------------

module iMonteCarlo
  interface MonteCarlo
    subroutine MonteCarlo_s(harvest)
    use kinds
    real(sp), intent(out):: harvest
    end subroutine MonteCarlo_s

    subroutine MonteCarlo_sv(harvest)
    use kinds
    real(sp), dimension(:), intent(out):: harvest
    end subroutine MonteCarlo_sv

    subroutine MonteCarlo_d(harvest)
    use kinds
    real(dp), intent(out):: harvest
    end subroutine MonteCarlo_d

    subroutine MonteCarlo_dv(harvest)
    use kinds
    real(dp), dimension(:), intent(out):: harvest
    end subroutine MonteCarlo_dv
  end interface
end module iMonteCarlo

!------------------------------------------------------------------------------

module iBuildFlorent
  interface
    recursive subroutine BuildFlorent(CurrentCell,CurrentLevel)
      use CellStructure
      implicit none
      type(Cell), target::  CurrentCell
      integer, intent(in):: CurrentLevel
    end subroutine BuildFlorent
  end interface
end module iBuildFlorent

!------------------------------------------------------------------------------

module iBuildNestedUniverse
  interface
    recursive subroutine BuildNestedUniverse(CurrentCell,CurrentLevel)
      use CellStructure
      implicit none
      type(Cell), target::  CurrentCell
      integer, intent(in):: CurrentLevel
    end subroutine BuildNestedUniverse
  end interface
end module iBuildNestedUniverse

!------------------------------------------------------------------------------

module iBuildMultiphase
  interface
    recursive subroutine BuildMultiphase(CurrentCell,CurrentLevel)
      use CellStructure
      implicit none
      type(Cell), target::  CurrentCell
      integer, intent(in):: CurrentLevel
    end subroutine BuildMultiphase
  end interface
end module iBuildMultiphase

!------------------------------------------------------------------------------

module iBuildSemireal
  interface
    recursive subroutine BuildSemireal(CurrentCell,CurrentLevel)
      use CellStructure
      implicit none
      type(Cell), target::  CurrentCell
      integer, intent(in):: CurrentLevel
    end subroutine BuildSemireal
  end interface
end module iBuildSemireal

!------------------------------------------------------------------------------

module iBuildShell
  interface
    recursive subroutine BuildShell(CurrentCell,CurrentLevel)
      use CellStructure
      implicit none
      type(Cell), target::  CurrentCell
      integer, intent(in):: CurrentLevel
    end subroutine BuildShell
  end interface
end module iBuildShell

!------------------------------------------------------------------------------

module iBuildSlab
  interface
    recursive subroutine BuildSlab(CurrentCell,CurrentLevel)
      use CellStructure
      implicit none
      type(Cell), target::  CurrentCell
      integer, intent(in):: CurrentLevel
    end subroutine BuildSlab
  end interface
end module iBuildSlab

!------------------------------------------------------------------------------

module iDigDeeper
  interface
    recursive subroutine DigDeeper(X_pos,CurrentCell)
      use kinds
      use CellStructure
      implicit none
      real(dp), intent(IN)::   X_pos(3)
      type(Cell), pointer::  CurrentCell        !intent(INOUT)
    end subroutine DigDeeper
  end interface
end module iDigDeeper

!------------------------------------------------------------------------------

module iDoppler
  interface Doppler
    function Doppler_s(T)
      use kinds
      use PhysicalConstants
      implicit none
      real(sp),intent(in):: T
      real(sp)::            Doppler_s
    end function Doppler_s
    function Doppler_sv(T)
      use kinds
      use PhysicalConstants
      implicit none
      real(sp),intent(in):: T(:)
      real(sp)::            Doppler_sv(size(T))
    end function Doppler_sv
    function Doppler_d(T)
      use kinds
      use PhysicalConstants
      implicit none
      real(dp),intent(in):: T
      real(dp)::            Doppler_d
    end function Doppler_d
    function Doppler_dv(T)
      use kinds
      use PhysicalConstants
      implicit none
      real(dp),intent(in):: T(:)
      real(dp)::            Doppler_dv(size(T))
    end function Doppler_dv
  end interface
end module iDoppler

!------------------------------------------------------------------------------

module iCalcFF
  interface
    recursive subroutine CalcFF(CurrentCell)
      use kinds
      use CellStructure
      use GlobalData
      implicit none
      type(Cell), target:: CurrentCell
    end subroutine CalcFF
  end interface
end module iCalcFF

!------------------------------------------------------------------------------

module iTemperature
  interface Temperature
    function Temperature_s(Dnu_D)
      use kinds
      use PhysicalConstants
      implicit none
      real(sp),intent(in):: Dnu_D
      real(sp)::            Temperature_s
    end function Temperature_s
    function Temperature_sv(Dnu_D)
      use kinds
      use PhysicalConstants
      implicit none
      real(sp),intent(in):: Dnu_D(:)
      real(sp)::            Temperature_sv(size(Dnu_D))
    end function Temperature_sv
    function Temperature_d(Dnu_D)
      use kinds
      use PhysicalConstants
      implicit none
      real(dp),intent(in):: Dnu_D
      real(dp)::            Temperature_d
    end function Temperature_d
    function Temperature_dv(Dnu_D)
      use kinds
      use PhysicalConstants
      implicit none
      real(dp),intent(in):: Dnu_D(:)
      real(dp)::            Temperature_dv(size(Dnu_D))
    end function Temperature_dv
  end interface
end module iTemperature

!------------------------------------------------------------------------------

module iGetCellData
  interface
    recursive subroutine GetCellData(HostCell,DirVec)!,oldCell,DirVec)
      use CellStructure
      implicit none
      type(Cell), pointer::  HostCell!,OldCell   !intent(in)
      integer, intent(IN)::  DirVec(3)
    end subroutine GetCellData
  end interface
end module iGetCellData

!------------------------------------------------------------------------------

module iCellPath
  interface
    subroutine CellPath(CurrentCell)
      use GlobalData
      use CellStructure
      implicit none
      type(Cell), pointer::            CurrentCell!intent(in)
    end subroutine CellPath
  end interface
end module iCellPath

!------------------------------------------------------------------------------

module iEnterNeighbor
  interface
    recursive subroutine EnterNeighbor(CurrentCell,X_cut,DirVec,i,j,k,escape)
      use kinds
      use CellStructure
      use GlobalData
      use iLocateHostCell
      implicit none
      type(Cell), pointer::   CurrentCell     !intent(INOUT)
      real(dp),intent(IN)::     X_cut(3)
      integer, intent(IN)::   DirVec(3)
      integer,intent(INOUT):: i,j,k
      logical,intent(OUT)::   escape
    end subroutine EnterNeighbor
  end interface
end module iEnterNeighbor

!------------------------------------------------------------------------------

module iEveryoneShouldHaveACenter
  interface
    recursive subroutine EveryoneShouldHaveACenter(CurrentCell,i,j,k,ParentPos)
      use kinds
      use CellStructure
      use GlobalData
      implicit none
      type(Cell), target::           CurrentCell
      integer,intent(in)::            i,j,k
      real(kind=posKind),intent(in):: ParentPos(3)
    end subroutine EveryoneShouldHaveACenter
  end interface
end module iEveryoneShouldHaveACenter

!------------------------------------------------------------------------------

module iPrintCell
  interface
    subroutine PrintCell(CurrentCell)
      use CellStructure
      implicit none
      type(Cell), pointer::  CurrentCell              !intent(in)
    end subroutine PrintCell
  end interface
end module iPrintCell

!------------------------------------------------------------------------------

module iScatPhoton
  interface
    recursive subroutine ScatPhoton(x,nhat_i,nhat_f,HostCell,absorb)
      use kinds
      use AtomParameters
      use DerivedParameters
      use CellStructure
      use CurrentCellParameters
      use PhysicalConstants
      implicit none
      real(dp), intent(inout):: x
      real(dp), intent(in)::    nhat_i(3)
      real(dp), intent(out)::   nhat_f(3)
      logical, intent(out)::  absorb
      type(cell), pointer::   HostCell          !intent IN!
    end subroutine ScatPhoton
  end interface
end module iScatPhoton

!------------------------------------------------------------------------------

module iExposeCCDs
  interface
    subroutine ExposeCCDs(x,X_pos,nhat_i,ScatCell,FreeEsc)
      use kinds
      use AtomParameters
      use DerivedParameters
      use CellStructure!check
      use CurrentCellParameters
      use PhysicalConstants
      implicit none
      real(dp), intent(in)::  x
      real(dp), intent(in)::  X_pos(3),nhat_i(3)
      logical, intent(in):: FreeEsc
      type(cell), pointer:: ScatCell, StartCell        !intent IN!
    end subroutine ExposeCCDs
  end interface
end module iExposeCCDs

!------------------------------------------------------------------------------

module iPropagate
  interface
    subroutine Propagate(x,X_start,DirVec,nhat_i,tau_esc,PropCell,FreeEsc)
      use kinds
      use CellStructure
      use PhysicalConstants
      use iPrintCell
      implicit none
      real(dp), intent(in)::  tau_esc,x
      real(dp), intent(in)::  X_start(3),nhat_i(3)
      integer,intent(in)::  DirVec(3)
      logical, intent(in):: FreeEsc
      type(cell), pointer:: PropCell                  !intent IN!
    end subroutine Propagate
  end interface
end module iPropagate

!------------------------------------------------------------------------------

module iLorentzTransform
  interface
    subroutine LorentzTransform(x,nhat_i,oldCell,newCell)
      use kinds
      use CellStructure
      implicit none
      real(dp), intent(INOUT):: x
      real(dp), intent(IN)::    nhat_i(3)
      type(Cell), pointer::   oldCell, newCell
    end subroutine LorentzTransform
  end interface
end module iLorentzTransform

!------------------------------------------------------------------------------

module iGoToCell
  interface
    subroutine GoToCell(CurrentCell)
      use GlobalData
      use CellStructure
      implicit none
      type(Cell), pointer:: CurrentCell       !intent(IN)
    end subroutine GoToCell
  end interface
end module iGoToCell

!------------------------------------------------------------------------------

!==============================================================================

program MoCaLaTA

!-----------------  Declarations  ---------------------------------------------

!use ran_state, only: ran_seed !NumRec

use kinds
use PhysicalConstants
use DataArrays
use GlobalData
use IJK
use ObservationalParameters
use CurrentCellParameters
use DerivedParameters
use AtomParameters
use iLocateHostCell
use iDigDeeper
use iGetCellData
use iCellPath
use iEnterNeighbor
use iScatPhoton
use iPrintCell
use iLorentzTransform
use iExposeCCDs
use iGoToCell
use iDoppler
use iTemperature

implicit none
integer::               n                       !Photon counter
integer,allocatable::   N_0(:)                  !# of clouds met for a given photon
real(dp),dimension(3):: X_pos,nhat_i,nhat_f,X_i,X_f,X_if,X_esc
real(dp)::              xc,ss,xx,x,x_old ,tau,r,dist,ddd(6)
real(dp)::              sigma,XsecD,norm,Voigt,percentile,CoreSkip
integer(kind=8)::       scats,s1,s2
logical::               LeaveCell,Scatter,escape,absorb
integer::               photon                  !Photon type (1 = Lya, 2 = FUV)
real(dp)::              R1,R2,R3,ran
real(dp)::              oldD,oldU(3)
integer::               DirVec(3)
integer::               oldCloud,newCloud
type(Cell), pointer::   HostCell => null(), oldCell => null()

!-----------------  Initialization  -------------------------------------------

call ReadInput
call AllocateCCDs

select case (trim(model))
  case ('data');       call ReadData
  case ('slab');       call MakeSlab
  case ('shell');      call MakeShell
  case ('multiphase'); call ReadMultiphase
  case ('semireal');   call MakeSemireal
  case ('florent');    call MakeFlorent
  case default;        stop 'Argh unknown model in Initialization!'
end select

D_obs = D_obs * kpc
if (D_obs .gt. D_box) D_obs = D_box
R_obs = D_obs / 2d0
o1    = R_box - R_obs                           !"Left" coordinate of obs. box
o2    = R_box + R_obs                           !"Right" coordinate of obs. box

call CalcDX
call ConstructUniverse
!call CalcClumpiness
call CalcHImass((/1.04_dp,2.25_dp,-0.34_dp/), 16.7_dp)

if (N_los .gt. 0) then
  if (abs(axes(1)-axes(2)).lt.1e-4 .and. abs(axes(1)-axes(3)).lt.1e-4) then
    call CalcAverageQuantities('all')
  else
    call CalcAverageQuantities('xyz')
  endif
endif
if (trim(X_QSOc).ne.'nil') call CalcAverageQuantities('qso')

if (writeN0) then
  allocate(N_0(n_phtot))
  N_0 = 0
  call getlun(N0lun)
  open(N0lun,file=trim(DataDir)//trim(Subdir)//'N0'//trim(outlabel)//'.dat',action='write',form='formatted',status='replace')
endif

! deallocate data arrays
if (associated(n_HIString))  deallocate(n_HIString)
if (associated(Dnu_DString)) deallocate(Dnu_DString)
if (associated(U_xString))   deallocate(U_xString)
if (associated(U_yString))   deallocate(U_yString)
if (associated(U_zString))   deallocate(U_zString)
# ifdef AMR
if (associated(LevelString)) deallocate(LevelString)
# endif
# ifdef multi
if (associated(PhaseString)) deallocate(PhaseString)
if (associated(CloudString)) deallocate(CloudString)
# endif
# ifdef dust
if (associated(n_dString))   deallocate(n_dString)
# endif

ArchitectureOfAggression = .true.
fullboxRT                = .true.
NBaper    = 5.                                  !Aperture/2 for Lya photons. Outside [lambda_0-NBaper, lambda_0+NBaper], FUV photons don't scatter (significantly)
n_eff     = 0d0
n_FUV     = 0d0
n_packets = 0
n_abs     = 0d0
n_absFUV  = 0d0
sigma_d   = 0d0                                 !Stays zero if dust is not included

!-----------------  Photon creation  -------------------------------------------

call cpu_time(t1)                               !Initial time
do n = 1,n_phtot                                !Loop over every photon
  if (writeN0) N_0(n) = 0
  oldCloud = -666
  dist = 0
  call which(photon)
  call InitPos(n,photon,X_pos)
  X_esc = X_pos
  call LocateHostCell(X_pos,HostCell,i,j,k,escape)
  DirVec = (/7,9,13/)
  call GetCellData(HostCell,DirVec)             !Physical parameters of cell

  call InitDir(nhat_i)                          !Emit isotropically
  call InitTau(tau)                             !Optical depth that photon will reach
  call InitFreq(Dnu_D,nhat_i,x,photon,U_bulk)   !Assign initial frequency
  if (len_trim(specin).gt.0) call sample(inlun,x,nhat_i,U_bulk,Dnu_D,X_pos,dist,photon)

  sigma_x = sigma(Dnu_D,x)                      !Atomic cross section
# ifdef dust
  sigma_d = XsecD(Dnu_D,x)                      !Dust cross section
# endif
  LeaveCell = .false.

  call ExposeCCDs(x,X_pos,nhat_i,HostCell,.true.)

!-----------------  Begin journey  --------------------------------------------

  do                                            !Loop until photon escapes Box
    r     = min(2d0*D_box,tau / (n_HI*sigma_x + n_d*sigma_d))  !Distance traveled
    X_i   = X_pos
    X_pos = X_pos + r*nhat_i                    !Update position
    X_f   = X_pos

    call CheckIntersect(X_i,X_f,nhat_i,DirVec,Scatter,LeaveCell)

    if (LeaveCell) then                         !Photon escapes cell
      oldCell => HostCell
      call EnterNeighbor(HostCell,X_f,DirVec,i,j,k,escape)
#    ifdef multi
     if (writeN0) then
       newCloud = HostCell%clnum
       if (HostCell%phase.eq.1 .and. newCloud.ne.oldCloud) then
         N_0(n)   = N_0(n) + 1
         oldCloud = newCloud
       endif
     endif
#    endif
      X_if = X_f - X_i
      r    = norm(X_if)                         !Distance to edge of cell
      dist = dist + r                           !Update total distance for current photon
      if (escape) then                          !Photon escapes Universe
        if (len_trim(specout).gt.0) call sample(outlun,x,nhat_i,U_bulk,Dnu_D,X_esc,dist,photon)
        exit
      endif
!     if (r/kpc .lt. 1d-15) Scatter = .true.    !Check, and possibly make r/dxTiny
      tau   = tau - r*(n_HI*sigma_x+n_d*sigma_d)!Remaining tau
      X_pos = X_f                               !Update position

      call GetCellData(HostCell,DirVec)         !Physical parameters of cell ijk
      call LorentzTransform(x,nhat_i,oldCell,HostCell)
    endif

    if (Scatter) then                           !Photon stays in cell and is scattered
      dist = dist + r                           !Update total distance for current photon
      x_old = x
      if (.not. ArchitectureOfAggression) then
        if (r .gt. .1*minval(ddd)) then
          ddd   = (/xp-X_pos(1), X_pos(1)-xm, & !Distances from photon to
                    yp-X_pos(2), X_pos(2)-ym, & !edges of host cell
                    zp-X_pos(3), X_pos(3)-zm/)  !
          tau_0 = minval(ddd) * n_HI * sigma(Dnu_D,0d0) !From photon to nearest cell face
          at    = a * tau_0                     !From photon to nearest cell face
          if (intxc) x_crit = CoreSkip(at)      !Threshold for core-skipping acceleration scheme
          if (maxxc) x_crit = min(x_crit,xc_max)
        endif
      endif

      call ScatPhoton(x,nhat_i,nhat_f,HostCell,absorb) !Assign new direction and frequency
      if (absorb) then
        n_abs(photon) = n_abs(photon) + w_ph
        if (trueFUV) n_absFUV = n_absFUV + w_ph
        exit
      endif
      X_esc = X_pos
      call ExposeCCDs(x_old,X_pos,nhat_i,HostCell,.false.)

      nhat_i = nhat_f

      call InitTau(tau)                       !...that photon will reach
    endif

    sigma_x = sigma(Dnu_D,x)                    !Re-calc cross section
#   ifdef dust
    sigma_d = XsecD(Dnu_D,x)                    !
#   endif
  enddo

  if (writeN0)                                             &
  write(N0lun,'(i4,4f8.2)') N_0(n),                        &!# of clouds for this photon
                    sum(N_0(1:n))/real(n),                 &!Mean so far
                    percentile(n,real(N_0(1:n),dp),50d0),  &!Median so far
                    percentile(n,real(N_0(1:n),dp),15.9d0),&!16th perc.
                    percentile(n,real(N_0(1:n),dp),84.1d0)  !84th perc.

  if (mod(n,wrtfl).eq.0 .or. n.eq.n_phtot) call ReadOutCCD(n)
enddo

end program MoCaLaTA

!==============================  Subroutines  =================================

subroutine AllocateCCDs

!Allocate arrays for 3D CCDs and 1D spectra, in 6 directions.

use kinds
use PhysicalConstants
use ObservationalParameters

implicit none
real(dp)::  lambda
integer:: i

if (WhichCCDs(1)) then
  allocate(CCDxp(pix,pix,SpecRes2D))
  allocate(specxp(SpecRes1D))
  CCDxp  = 0
  specxp = 0
endif

if (WhichCCDs(2)) then
  allocate(CCDxm(pix,pix,SpecRes2D))
  allocate(specxm(SpecRes1D))
  CCDxm  = 0
  specxm = 0
endif

if (WhichCCDs(3)) then
  allocate(CCDyp(pix,pix,SpecRes2D))
  allocate(specyp(SpecRes1D))
  CCDyp  = 0
  specyp = 0
endif

if (WhichCCDs(4)) then
  allocate(CCDym(pix,pix,SpecRes2D))
  allocate(specym(SpecRes1D))
  CCDym  = 0
  specym = 0
endif

if (WhichCCDs(5)) then
  allocate(CCDzp(pix,pix,SpecRes2D))
  allocate(speczp(SpecRes1D))
  CCDzp  = 0
  speczp = 0
endif

if (WhichCCDs(6)) then
  allocate(CCDzm(pix,pix,SpecRes2D))
  allocate(speczm(SpecRes1D))
  CCDzm  = 0
  speczm = 0
endif

DDlam = BW(2) - BW(1)
dlam  = DDlam / SpecRes1D

end subroutine AllocateCCDs

!------------------------------------------------------------------------------

subroutine ConstructUniverse

use kinds
use GlobalData
use DataArrays
use PhysicalConstants
use ObservationalParameters
use iBuildNestedUniverse
use iBuildFlorent
use iBuildMultiphase
use iBuildSemireal
use iBuildSlab
use iBuildShell
use iEveryoneShouldHaveACenter
use iAssignParameters

implicit none
integer::            i,j,k
real(kind=posKind):: DummyPos(3)
integer::            mm,nn,kk,lun
character(len=10)::  ani,f,ch
character(len=50)::  pre,epi

interface
  recursive subroutine SeekAndDestroy(CurrentCell)
    use DataArrays
    use CellStructure
    use iv2u
    use PhysicalConstants
    use GlobalData
    implicit none
    type(Cell), target:: CurrentCell
  end subroutine SeekAndDestroy
end interface

write(*,*) 'Constructing Universe'

write(*,*) ' - Initializing base grid...'
call InitializeBaseGrid

write(*,*) ' - Building fully threaded structure...'
select case (trim(model))
  case ('data')
    i_cell = 0
    do i = 1,ni
      do j = 1,nj
        do k = 1,nk
          call BuildNestedUniverse(BaseGrid%Child(i,j,k),0)
        enddo
      enddo
    enddo

  case ('slab')
    i_cell = 0
    do i = 1,ni
      do j = 1,nj
        do k = 1,nk
          call BuildSlab(BaseGrid%Child(i,j,k),0)
        enddo
      enddo
    enddo

  case ('shell')
    i_cell = 0
    do i = 1,ni
      do j = 1,nj
        do k = 1,nk
          call BuildShell(BaseGrid%Child(i,j,k),0)
        enddo
      enddo
    enddo

  case ('multiphase')
    i_cell = 0
    do i = 1,ni
      do j = 1,nj
        do k = 1,nk
          call BuildMultiphase(BaseGrid%Child(i,j,k),0)
        enddo
      enddo
    enddo

  case ('semireal')
    if (len_trim(Cloudfile) .gt. 0) then
      call getlun(cloudlun)
      open(cloudlun,file=trim(DataDir)//trim(Subdir)//trim(Cloudfile), status='replace',action='write')
    endif
    if (len_trim(Windfile) .gt. 0) then
      call getlun(windlun)
      open(windlun,file=trim(DataDir)//trim(Subdir)//trim(Windfile), status='replace',action='write')
    endif
    nn  = ni + 1
    call Int2Char(nn,ani,.false.)
    pre = 'The more clouds, the slower'
    epi = '|'
    mm  = len_trim(pre)
    kk  = len_trim(epi)
    call Int2Char(mm + 2*(nn+kk), ch,.false.)
    f   = '(a' // trim(ch)  // ')'
    write(*,f,advance='no') trim(pre)//repeat(' ',nn)//trim(epi)//repeat(achar(8),nn+kk)
    i_cell = 0
    do i = 1,ni
      write(*,'(a)',advance='no') '.'; flush(6)
      do j = 1,nj
        do k = 1,nk
          call BuildSemireal(BaseGrid%Child(i,j,k),0)
        enddo
      enddo
    enddo
    write(*,'(a)',advance='yes') '.'
    ! To write cloud center positions, uncomment this block
    ! open (99,file=trim(DataDir)//trim(Subdir)//'cloudcenters.dat', status='replace',action='write')
    ! do cl = 1,N_cl
    !   write(99,*) Clouds(cl,:) / kpc
    ! enddo
    ! close(99)

  case ('florent')
    i_cell = 0
    do i = 1,ni
      do j = 1,nj
        do k = 1,nk
          call BuildFlorent(BaseGrid%Child(i,j,k),0)
        enddo
      enddo
    enddo

  case default
    stop 'Argh unknown model in ConstructUniverse!'
end select

# ifdef AMR
write(*,*) ' - Assigning cell centers...'
DummyPos = (/6,6,6/)
do i = 1,ni
  do j = 1,nj
    do k = 1,nk
      call EveryoneShouldHaveACenter(BaseGrid%Child(i,j,k),i,j,k,DummyPos)
    enddo
  enddo
enddo
# endif

if (trim(model) .eq. 'multiphase') then 
  do i = 1,ni
    do j = 1,nj
      do k = 1,nk
        call AssignParameters(BaseGrid%Child(i,j,k))
      enddo
    enddo
  enddo
endif

! !FIRST Build, THEN Everyone, THEN Seek...
! i_cell = 0
! do i = 1,ni
!   do j = 1,nj
!     do k = 1,nk
!       call SeekAndDestroy(BaseGrid%Child(i,j,k))
!     enddo
!   enddo
! enddo
! stop

! call getlun(lun)
! open(lun,file='fescxc.dat',status='old',action='read')
! read(lun,*) fescxc
! close(lun)

fescxc(1:10, 1) = (/1.0000000, 1.0000000, 1.0000000, 1.0000000, 1.0000000, 1.0000000, 1.0000000, 1.0000000, 1.0000000, 1.0000000/)
fescxc(1:10, 2) = (/1.0000000, 1.0000000, 1.0000000, 1.0000000, 1.0000000, 1.0000000, 1.0000000, 1.0000000, 1.0000000, 1.0000000/)
fescxc(1:10, 3) = (/1.0000000, 1.0000000, 1.0000000, 1.0000000, 1.0000000, 1.0000000, 1.0000000, 1.0000000, 1.0000000, 1.0000000/)
fescxc(1:10, 4) = (/1.0000000, 1.0000000, 1.0000000, 1.0000000, 1.0000000, 1.0000000, 1.0000000, 1.0000000, 1.0000000, 1.0000000/)
fescxc(1:10, 5) = (/0.9434176, 0.9749770, 0.9923456, 1.0000000, 1.0000000, 1.0000000, 1.0000000, 1.0000000, 1.0000000, 1.0000000/)
fescxc(1:10, 6) = (/0.5217663, 0.6939473, 0.9034173, 0.9538222, 0.9749770, 0.9990770, 1.0000000, 1.0000000, 1.0000000, 1.0000000/)
fescxc(1:10, 7) = (/0.1135932, 0.2421210, 0.5217668, 0.7250708, 0.9331264, 0.9749770, 1.0000000, 1.0000000, 1.0000000, 1.0000000/)
fescxc(1:10, 8) = (/0.0304603, 0.0963607, 0.2701889, 0.4885344, 0.7829334, 0.9434176, 0.9802983, 0.9941230, 1.0000000, 1.0000000/)
fescxc(1:10, 9) = (/0.0152628, 0.0544752, 0.1561327, 0.3115964, 0.6287170, 0.8641639, 0.9544670, 0.9815567, 0.9991945, 1.0000000/)
fescxc(1:10,10) = (/0.0109832, 0.0383502, 0.0942698, 0.2169690, 0.5104452, 0.7915681, 0.9187273, 0.9740932, 0.9991945, 1.0000000/)
fescxc(1:10,11) = (/0.0090154, 0.0275971, 0.0685851, 0.1649342, 0.4099015, 0.7093384, 0.8865787, 0.9642301, 0.9966555, 1.0000000/)
fescxc(1:10,12) = (/0.0068533, 0.0209786, 0.0532935, 0.1353842, 0.3255715, 0.6497490, 0.8512126, 0.9496225, 0.9941230, 1.0000000/)
fescxc(1:10,13) = (/0.0055641, 0.0168463, 0.0418672, 0.1099162, 0.2701887, 0.5634046, 0.8007910, 0.9328597, 0.9865640, 1.0000000/)
fescxc(1:10,14) = (/0.0047200, 0.0138281, 0.0314793, 0.0854086, 0.2217811, 0.4624638, 0.7509746, 0.9025099, 0.9740932, 1.0000000/)
fescxc(1:10,15) = (/0.0038743, 0.0117303, 0.0255573, 0.0716615, 0.1881362, 0.3923066, 0.6735696, 0.8798374, 0.9666865, 1.0000000/)
fescxc(1:10,16) = (/0.0033965, 0.0097349, 0.0226529, 0.0601271, 0.1527447, 0.3364626, 0.6201661, 0.8490497, 0.9496225, 1.0000000/)
fescxc(1:10,17) = (/0.0029130, 0.0078174, 0.0181908, 0.0521369, 0.1253787, 0.2917505, 0.5752749, 0.8110432, 0.9352361, 0.9991945/)
fescxc(1:10,18) = (/0.0026683, 0.0070053, 0.0161237, 0.0452085, 0.1051981, 0.2529801, 0.5130416, 0.7688471, 0.9140642, 0.9890773/)
fescxc(1:10,19) = (/0.0024441, 0.0064168, 0.0146070, 0.0409589, 0.0912184, 0.2169690, 0.4563781, 0.7307031, 0.8933716, 0.9840572/)
fescxc(1:10,20) = (/0.0022635, 0.0058777, 0.0138284, 0.0387731, 0.0782337, 0.1840541, 0.4012590, 0.6770058, 0.8820788, 0.9765747/)

do i = 1,10
  atA = log10(4.7017128) + (i-1)/2d0            ! [.67, 1.17, ..., 5.17]
enddo
do i = 1,20
  xcA = (i-1) * .5d0                            ! [0, .5, ..., 9.5]
enddo
dat = atA(2) - atA(1)                           ! = 0.5
dxc = xcA(2) - xcA(1)                           ! = 0.5

end subroutine ConstructUniverse

!------------------------------------------------------------------------------

subroutine CalcDX

use CellStructure
use GlobalData
use DataArrays

implicit none
integer:: L

# ifdef AMR
L_max = maxval(LevelString)                     !Maximum level of refinement
# else
L_max = 0
# endif

allocate(dx(0:L_max))
allocate(dy(0:L_max))
allocate(dz(0:L_max))
allocate(dxH(0:L_max))
allocate(dyH(0:L_max))
allocate(dzH(0:L_max))

dx0 = D_box / ni
dy0 = D_box / nj
dz0 = D_box / nk

do L = 0,L_max
  dx(L) = dx0 / 2d0**L
  dy(L) = dy0 / 2d0**L
  dz(L) = dz0 / 2d0**L
enddo

dxH = dx / 2d0
dyH = dy / 2d0
dzH = dz / 2d0

dxTiny  = 0.1d0 * dx(L_max)
dxWeeny = 1d-6  * dx(L_max)

allocate(FullPathArray(3,0:L_max-1))

end subroutine CalcDX

!------------------------------------------------------------------------------

subroutine CalcAverageQuantities(dir)

use kinds
use CellStructure
use CurrentCellParameters
use PhysicalConstants
use GlobalData
use DataArrays
use ObservationalParameters
use iLocateHostCell
use iEnterNeighbor
use iGetCellData
use iLorentzTransform
use iCalcFF
use iPrintCell
use iMonteCarlo

implicit none
character(len=3), intent(in):: dir
integer::               i,j,k,n,DirVec(3),seed(20),oldCloud,newCloud,ncycles,ic
integer::               ccc,stalls,N_loseff
integer(1)::            onefix
real(dp),dimension(3):: nhat,X_pos,X_i,X_f
real(dp)::              x,r,norm,sigma,XsecD,percentile,FF_actual,X_if(3),eps,sigd,sigx,R1
type(Cell), pointer::   HostCell => null(), oldCell => null()
logical::               esc,scat,leave
logical,allocatable::   CF(:)
real(dp),allocatable::  tau_a(:),clLOS(:),NN_HI(:),t0(:)
real(dp):: ari_t0,ari_N_HII,ari_tau_a,ari_EBV,ari_f_c,ari_N_0,ari_CF
real(dp):: med_t0,med_N_HII,med_tau_a,med_EBV,med_f_c,lo_f_c,hi_f_c,med_N_0,lo_N_HII,hi_N_HII

if (dir .eq. 'qso') then
  write(*,*) 'Re-setting N_los to 3.'
  N_los = 3                                       !# of sightlines
  onefix = 0
  if (X_QSO(1) .lt. dxWeeny) then
    X_QSO(1) = dxWeeny
    nhat     = (/1,0,0/)
    onefix   = onefix + 1
  endif
  if (X_QSO(2) .lt. dxWeeny) then
    X_QSO(2) = dxWeeny
    nhat     = (/0,1,0/)
    onefix   = onefix + 1
  endif
  if (X_QSO(3) .lt. dxWeeny) then
    X_QSO(3) = dxWeeny
    nhat     = (/0,0,1/)
    onefix   = onefix + 1
  endif
  if (D_box-X_QSO(1) .lt. dxWeeny) then
    X_QSO(1) = D_box - dxWeeny
    nhat     = (/-1,0,0/)
    onefix   = onefix + 1
  endif
  if (D_box-X_QSO(2) .lt. dxWeeny) then
    X_QSO(2) = D_box - dxWeeny
    nhat     = (/0,-1,0/)
    onefix   = onefix + 1
  endif
  if (D_box-X_QSO(3) .lt. dxWeeny) then
    X_QSO(2) = D_box - dxWeeny
    nhat     = (/0,0,-1/)
    onefix   = onefix + 1
  endif
  if (onefix .ne. 1_4) then
    write(*,*) "Something's wrong with the QSO position:"
    write(*,*) X_QSO / kpc
    write(*,*) 'Also, it was fixed', onefix, 'times.'
    stop
  endif
endif

allocate(t0(N_los))   
allocate(NN_HI(N_los))         !Array for collecting N_HI
allocate(tau_a(N_los))         !Array for collecting tau_a
# ifdef multi
allocate(clLOS(N_los))         !Array for collecting # clouds
allocate(CF(N_los))            !Array for >= 1-cloud LOSs

if (.not.allocated(FFcounter)) allocate(FFcounter(0:n_phases))
if (dir.ne.'qso' .and. n_phases.ge.0) then
  FFcounter = 0d0
  do i = 1,ni
    do j = 1,nj
      do k = 1,nk
        call CalcFF(BaseGrid%Child(i,j,k))
      enddo
    enddo
  enddo 
  write(*,*) 'Filling factor of clouds:       ', real(FFcounter(1) / sum(FFcounter(1:n_phases)))
  write(*,*) 'Non-void volume fraction of box:', &
           real(sum(FFcounter(1:n_phases)) / sum(FFcounter))
endif
# endif

ncycles = merge(3,1,dir .eq. 'xyz')

do ic = 1,ncycles
  write(*,*) '--------------------- Average quantities for ' // dir(ic:ic) // '-direction ---------------------'
  t0    = 0
  NN_HI = 0      
  tau_a = 0      
# ifdef multi
  clLOS = 0d0    
  CF    = .false.
# endif

  do n = 1,N_los
    ccc      = 0
    oldCloud = -666
    if (dir .eq. 'all') then
    ! call InitPos(0,1,X_pos) ! <-- use this for starting sightlines the same
                              !     locations as photons, but make sure to reset
                              !     n_eff(photon) afterwards (look in InitPos...)
      call CentralEmission(X_pos)
      call IsoDir(nhat)
    elseif (dir .eq. 'xyz') then
      call NearCentralEmission(ic,X_pos)
      call MonteCarlo(R1)
      nhat     = (/0,0,0/)
      nhat(ic) = merge(-1,1,R1.gt..5)
    elseif (dir .eq. 'qso') then
      X_pos = X_QSO
    else
      stop 'What... what direction...?'
    endif

    call LocateHostCell(X_pos,HostCell,i,j,k,esc)
    DirVec = (/7,9,13/)
    call GetCellData(HostCell,DirVec)!,oldCell,DirVec)
    x = 0

    do
      ccc   = ccc + 1
      X_i   = X_pos
      X_pos = X_pos + D_box*nhat                  !Somewhere outside box
      X_f   = X_pos

      call CheckIntersect(X_i,X_f,nhat,DirVec,scat,leave)

      oldCell => HostCell
      call EnterNeighbor(HostCell,X_f,DirVec,i,j,k,esc)

      X_if     = X_f - X_i
      r        = norm(X_if)
      NN_HI(n) = NN_HI(n) + r*n_HI
      t0(n)    = t0(n)    + r*n_HI * sigma(Dnu_D,0d0)
      tau_a(n) = tau_a(n) + r*n_d    * XsecD(Dnu_D,x) * (1-albedo)
#     ifdef multi
      newCloud = HostCell%clnum
      if (HostCell%phase.eq.1 .and. newCloud.ne.oldCloud) then
        clLOS(n) = clLOS(n) + 1d0
        CF(n)    = .true.
        oldCloud = newCloud
      endif
    ! The following may over-count a cloud if LOS grazes its outer cells:
    ! if (HostCell%phase.eq.1 .and. oldCell%phase.eq.2) then !Going from ICM to cloud
    !   clLOS(n) = clLOS(n) + 1d0
    !   CF(n)    = .true.
    ! endif
#     endif
      if (esc) exit

      X_pos = X_f
      call GetCellData(HostCell,DirVec)!,oldCell,DirVec)
      call LorentzTransform(x,nhat,oldCell,HostCell)
    enddo
! print*, ccc
  enddo

  ari_t0    = sum(t0)    / N_los
  ari_N_HII = sum(NN_HI) / N_los
  ari_tau_a = sum(tau_a) / N_los
  ari_EBV   = ta2EBV * ari_tau_a
#   ifdef multi
  ari_f_c   = sum(clLOS) / N_los
  ari_CF    = count(CF) / real(N_los)
#   endif

  med_t0    = percentile(N_los,t0,50d0)
  med_N_HII = percentile(N_los,NN_HI,50d0)
  lo_N_HII = percentile(N_los,NN_HI,15.8655d0)
  hi_N_HII = percentile(N_los,NN_HI,84.1345d0)
  med_tau_a = percentile(N_los,tau_a,50d0)
  med_EBV   = ta2EBV * med_tau_a
#   ifdef multi
  med_f_c   = percentile(N_los,clLOS,50d0)
  lo_f_c    = percentile(N_los,clLOS,15.8655d0)
  hi_f_c    = percentile(N_los,clLOS,84.1345d0)
#   endif

  write(*,*) '           Ari. mean        Median        16%               84%'
  write(*,*) '<tau_0> ', real((/ari_t0,    med_t0/))
  write(*,*) '<N_HI>  ', real((/ari_N_HII, med_N_HII, lo_N_HII,hi_N_HII/))
  write(*,*) '<tau_a> ', real((/ari_tau_a, med_tau_a/))
  write(*,*) 'E(B-V)  ', real((/ari_EBV  , med_EBV/))
#   ifdef multi
  eps   = 0
#   ifdef dust
  Dnu_D = Dnu_DString(1)
  sigd  = XsecD(Dnu_D,x_inj)
  sigx  = sigma(Dnu_D,x_inj)
  eps   = (1-albedo) * n_dString(1)*sigd / &
                     (n_HIString(1)*sigx + n_dString(1)*sigd)
#   endif
  write(*,*) 'f_c     ', real((/ari_f_c  , med_f_c, lo_f_c,hi_f_c/))
  write(*,*) 'C.F.    ', real(ari_CF)
#   endif
  write(*,*)
enddo

if (dir .eq. 'qso') then
  obspar(3) = ari_N_HII
  obspar(4) = ari_EBV
else
  obspar(1) = ari_N_HII
  obspar(2) = ari_EBV
endif

deallocate(t0)
deallocate(NN_HI)
deallocate(tau_a)
# ifdef multi
  deallocate(clLOS)
  deallocate(CF)
# endif

end subroutine CalcAverageQuantities

!------------------------------------------------------------------------------

subroutine CalcClumpiness

use kinds
use GlobalData
use DataArrays
use PhysicalConstants

implicit none
real(dp),allocatable::  Clumpy(:)
real(dp)::              r1,r2,Mean,MeanOfSquares,MeanSquared,Swn,Swn2
real(dp)::              V,V_tot,M_HI,M_d,Zmet
integer::               i,j,k

r1     = r_inner
r2     = r_outer

M_HI   = 0
M_d    = 0
V_tot  = 0
Swn    = 0
Swn2   = 0

do i = 1,ni
  do j = 1,nj
    do k = 1,nk
      call SampleClumpsInCC(BaseGrid%Child(i,j,k))
    enddo
  enddo
enddo

MeanOfSquares = Swn2 / V_tot
Mean          = Swn  / V_tot
MeanSquared   = Mean**2

Clumpiness    = MeanOfSquares / MeanSquared
write(*,*) 'Clumpiness factor C = <n^2> / <n>^2:', real(Clumpiness,sp)
write(*,*) 'Total HI mass in galaxy [M_sun]:    ', real(M_HI / M_sun)
write(*,*) 'Total dust mass in galaxy [M_sun]:  ', real(M_d  / M_sun), &
  '(assuming dust-to-metal ratio 1/3)'

contains
  recursive subroutine SampleClumpsInCC(CurrentCell)
  use CellStructure
  implicit none
  type(Cell),target::  CurrentCell
  integer::            ii,jj,kk,L
  real(dp)::           X(3),XwrtC(3),d,n_HI,n_d,V,M_HIi,M_mi,norm
  real(dp),parameter:: Z_sun = .02              !Solar metallicity
  real(dp),parameter:: psi   = .76              !Hydrogen mass fraction
  real(dp),parameter:: D2M   = .3333            !Dust-to-metal ratio (D2M=1/3 => dust-to-nondustmetals = 1/2)
! real(dp),parameter:: x_HI  = 1.               !Dust-to-metal ratio (D2M=1/3 => dust-to-nondustmetals = 1/2)
  real(dp)::           x_HI                     !Dust-to-metal ratio (D2M=1/3 => dust-to-nondustmetals = 1/2)
# ifdef AMR
# ifdef multi
  if (CurrentCell%Refined) then
    do ii = 1,2
      do jj = 1,2
        do kk = 1,2
          call SampleClumpsInCC(CurrentCell%Child(ii,jj,kk))
        enddo
      enddo
    enddo
  else
    X     = CurrentCell%C_pos
    XwrtC = X - R_box
    d     = norm(XwrtC)
    if (d.ge.r_inner .and. d.le.r_outer) then
      L     = CurrentCell%Level
      V     = dx(L) * dy(L) * dz(L)             !Cell volume in cm3
      n_HI  = CurrentCell%n_HI                  !Cell HI density
      M_HIi = V * n_HI * m_H                    !Total HI mass in current cell
      M_HI  = M_HI + M_HIi                      !Accumulated galaxy HI mass so far
#     ifdef dust
      Zmet  = CurrentCell%n_d/n_HI * Z0 * Z_sun !Cell metallicity
    ! M_d   = M_d + Zmet*(V*n_HI*m_H/psi) / (1-Zmet) * D2M
      x_HI  = merge(1d0,1d-8,CurrentCell%phase.eq.1)
      M_mi  = Zmet * M_HIi/(x_HI*psi)           !Total metal mass in current cell
      M_d   = M_d + M_mi*D2M                    !Accumulated galaxy dust mass so far
#     endif
      Swn   = Swn   + V * n_HI
      Swn2  = Swn2  + V * n_HI**2
      V_tot = V_tot + V
    endif
  endif
# endif
# endif
  end subroutine SampleClumpsInCC

end subroutine CalcClumpiness

!------------------------------------------------------------------------------

subroutine CalcHImass(X_cen,r_vir)

use kinds
use GlobalData
use DataArrays
use PhysicalConstants

implicit none
real(dp),intent(in):: X_cen(3),r_vir            !Gal center and vir. rad. in kpc
real(dp),allocatable::  Clumpy(:)
real(dp)::              r1,r2,Mean,MeanOfSquares,MeanSquared,Swn,Swn2
real(dp)::              V,V_tot,M_HI,M_d,Zmet
integer::               i,j,k

M_HI   = 0
M_d    = 0
V_tot  = 0
Swn    = 0
Swn2   = 0

do i = 1,ni
  do j = 1,nj
    do k = 1,nk
      call SampleClumps(BaseGrid%Child(i,j,k),X_cen,r_vir)
    enddo
  enddo
enddo

MeanOfSquares = Swn2 / V_tot
Mean          = Swn  / V_tot
MeanSquared   = Mean**2

Clumpiness    = MeanOfSquares / MeanSquared
write(*,*) 'Clumpiness factor C = <n^2> / <n>^2:', real(Clumpiness,sp)
write(*,*) 'Total HI mass in galaxy [M_sun]:    ', real(M_HI / M_sun)
write(*,*) 'Total dust mass in galaxy [M_sun]:  ', real(M_d  / M_sun), &
  '(assuming dust-to-metal ratio 1/3)'

contains
  recursive subroutine SampleClumps(CurrentCell,X_cen,r_vir)
  use CellStructure
  implicit none
  type(Cell),target::   CurrentCell
  integer::             ii,jj,kk,L
  real(dp),intent(in):: X_cen(3),r_vir          !Gal center and vir. rad. in kpc
  real(dp)::            X(3),XwrtC(3),d,n_HI,n_d,V,M_HIi,M_mi,norm
  real(dp),parameter::  Z_sun = .02            !Solar metallicity
  real(dp),parameter::  psi   = .76            !Hydrogen mass fraction
  real(dp),parameter::  D2M   = .3333          !Dust-to-metal ratio (D2M=1/3 => dust-to-nondustmetals = 1/2)
! real(dp),parameter::  x_HI  = 1.             !Dust-to-metal ratio (D2M=1/3 => dust-to-nondustmetals = 1/2)
  real(dp)::            x_HI                   !Dust-to-metal ratio (D2M=1/3 => dust-to-nondustmetals = 1/2)
# ifdef AMR
  if (CurrentCell%Refined) then
    do ii = 1,2
      do jj = 1,2
        do kk = 1,2
          call SampleClumps(CurrentCell%Child(ii,jj,kk),X_cen,r_vir)
        enddo
      enddo
    enddo
  else
    X     = CurrentCell%C_pos
    XwrtC = X - R_box - X_cen*kpc
    d     = norm(XwrtC)
    if (d/kpc .lt. r_vir) then
      L     = CurrentCell%Level
      V     = dx(L) * dy(L) * dz(L)             !Cell volume in cm3
      n_HI  = CurrentCell%n_HI                  !Cell HI density
      if (L.gt.4) write(1,*) (X-R_box- X_cen*kpc)/kpc, L, n_HI, V
      M_HIi = V * n_HI * m_H                    !Total HI mass in current cell
      M_HI  = M_HI + M_HIi                      !Accumulated galaxy HI mass so far
! #     ifdef dust
!       Zmet  = CurrentCell%n_d/n_HI * Z0 * Z_sun !Cell metallicity
!     ! M_d   = M_d + Zmet*(V*n_HI*m_H/psi) / (1-Zmet) * D2M
!       x_HI  = merge(1d0,1d-8,CurrentCell%phase.eq.1)
!       M_mi  = Zmet * M_HIi/(x_HI*psi)           !Total metal mass in current cell
!       M_d   = M_d + M_mi*D2M                    !Accumulated galaxy dust mass so far
! #     endif
      Swn   = Swn   + V * n_HI
      Swn2  = Swn2  + V * n_HI**2
      V_tot = V_tot + V
    endif
  endif
# endif
  end subroutine SampleClumps

end subroutine CalcHImass

!------------------------------------------------------------------------------

subroutine CellPath(CurrentCell)

use CellStructure
use GlobalData

implicit none
type(Cell), pointer:: CurrentCell                !intent(IN)
integer::             l,e,le,i,zero
character(LEN=5)::    idx(3)
type(Cell), pointer:: TempCell => null()

zero           = iachar('0')
le             = len(idx(1))
FullPathString = ''
FullPathArray  = 0
TempCell       => CurrentCell

! do l = 0,0-1!CurrentCell%Level,0,-1
!   FullPathArray(:,l) = TempCell%indexico
!
!   do i = 1,3
!     idx(i) = ''
!     do e = 1,le
!       idx(i) = achar(zero + mod(TempCell%indexico(i),10**e)/10**(e-1)) // &
!                trim(idx(i))
!       if (TempCell%indexico(i) .lt. 10**e) exit
!     enddo
!   enddo
!
!   FullPathString = '%Child(' // trim(idx(1)) // ',' &
!                              // trim(idx(2)) // ',' &
!                              // trim(idx(3)) // ')' &
!                              // trim(FullPathString)
!   TempCell => TempCell%Parent
! enddo

FullPathString = 'BaseGrid' // trim(FullPathString)

FullPathString = 'Fix CellPath subroutine'

end subroutine CellPath

!------------------------------------------------------------------------------

subroutine CentralEmission(X_pos)

!Return position negligibly far from the center of the box. If a photon is
!positioned exactly at the face/corner of a cell, it may end up acquiring cell
!boundaries from one cell, but have its position in another, resulting in the
!wrath of Beelzebub.

use kinds
use GlobalData
use iMonteCarlo

implicit none
real(dp), intent(out):: X_pos(3)
real(dp)::              n(3),rsq

do
  call MonteCarlo(n)
  n   = 2*n - 1
  rsq = sum(n**2)
  if (rsq .le. 1) exit
enddo

n = n / sqrt(rsq)

X_pos = R_box + n*dxWeeny

end subroutine CentralEmission

!------------------------------------------------------------------------------

subroutine ChangeBasis(nhat_i)

!Change basis from expressing velocity of scattering atom relative to
!direction of incident photon to expressing it in Cell frame (= lab frame if
!U_bulk = 0).
!Return also unit vectors perpendicular to u_II

use kinds
use AtomParameters

implicit none
real(dp)::              nhat_i(3)
real(dp),dimension(3):: n_ixy,nhat_ixy,UU_II,UU_1,UU_2
real(dp)::              norm

!!n_ixy    = nhat_i                               !Projection of nhat_i onto
!!n_ixy(3) = 0.                                   !x-y plane
!!
!!nhat_ixy = n_ixy / norm(n_ixy)                  !Normalize to unit vector
!!
!!uhat_1(1) = -nhat_ixy(2)                        !Rotate nhat_ixy by pi/2 in
!!uhat_1(2) = nhat_ixy(1)                         !x-y plane to get first unit
!!uhat_1(3) = 0.                                  !vector perpendicular to u_II
!!
!!call CrossProduct(nhat_i,uhat_1,uhat_2)         !Take cross product of nhat_i
!!                                                !and uhat_1 to get last unit
!!                                                !vector
!!
!!UU_II = u_II * nhat_i                           !Velocity components parallel
!!UU_1  = u_1  * uhat_1                           !and perpendicular to nhat_i in
!!UU_2  = u_2  * uhat_2                           !reference frame of bulk motion
!!
!!U = UU_II + UU_1 + UU_2                         !Final velocity in Cell frame

if (abs(nhat_i(3)-1d0) .gt. 1d-12) then
  n_ixy    = (/nhat_i(1),nhat_i(2),0d0/)        !Projection of nhat_i onto x-y plane
  nhat_ixy = n_ixy / norm(n_ixy)                !Normalize to unit vector
  uhat_1   = (/-nhat_ixy(2),nhat_ixy(1),0d0/)   !Rotate nhat_ixy by pi/2 in x-y plane to get first unit vector perpendicular to u_II
  call CrossProduct(nhat_i,uhat_1,uhat_2)       !nhat_i X uhat_1 = final unit vector
else                                            !However, if nhat_i = (0,0,1),
  uhat_1 = (/1,0,0/)                            ! then simply use x and y axes
  uhat_2 = (/0,1,0/)
endif

UU_II = u_II * nhat_i                           !Velocity components parallel
UU_1  = u_1  * uhat_1                           !and perpendicular to nhat_i in
UU_2  = u_2  * uhat_2                           !reference frame of bulk motion

U     = UU_II + UU_1 + UU_2                     !Final velocity in Cell frame

end subroutine ChangeBasis

!------------------------------------------------------------------------------

subroutine CheckIntersect(X_i,X_f,nhat_i,DirVec,Scatter,LeaveCell)

!Determines the point of intersection between the path of the photon and the
!face of the cell. If the photon escapes the cell, it is put at the face and a
!flag is set for updating parameters. If it stays in the the cell, a flag is
!set for scattering the photon.

use kinds
use CurrentCellParameters
use PhysicalConstants
use IJK
use GlobalData

implicit none
real(dp),intent(IN)::    X_i(3),nhat_i(3)
real(dp),intent(INOUT):: X_f(3)
integer, intent(OUT):: DirVec(3)
logical, intent(OUT):: Scatter,LeaveCell
real(dp)::               X_cut(3),norm

Scatter   = .true.                              !Flag values if photon does not
LeaveCell = .false.                             !leave cell

!Check if photon crosses yz-plane in positive direction
if (X_f(1) .ge. xp) then
   X_cut(1) = xp                                           !Intersection
   X_cut(2) = X_i(2) + (xp - X_i(1)) * nhat_i(2)/nhat_i(1) !with
   X_cut(3) = X_i(3) + (xp - X_i(1)) * nhat_i(3)/nhat_i(1) !plane

   if (X_cut(2) .lt. yp .and. &                 !Check if it happened
       X_cut(2) .ge. ym .and. &                 !through face of cell
       X_cut(3) .lt. zp .and. &
       X_cut(3) .ge. zm)        then

      X_f       = X_cut                         !Place photon at face
      DirVec    = (/1,0,0/)                     !Continue in neighboring cell
      LeaveCell = .true.                        !Set flag for updating params
      Scatter   = .false.                       !Unset flag for scattering
   endif
endif

!Check if photon crosses yz-plane in negative direction
if (X_f(1) .lt. xm) then
   X_cut(1) = xm                                           !Intersection
   X_cut(2) = X_i(2) + (xm - X_i(1)) * nhat_i(2)/nhat_i(1) !with
   X_cut(3) = X_i(3) + (xm - X_i(1)) * nhat_i(3)/nhat_i(1) !plane

   if (X_cut(2) .lt. yp .and. &                 !If so, check if it happened
       X_cut(2) .ge. ym .and. &                 !through face of cell
       X_cut(3) .lt. zp .and. &
       X_cut(3) .ge. zm)         then

      X_f       = X_cut                         !Place photon at face
      DirVec    = (/-1,0,0/)                    !Continue in neighboring cell
      LeaveCell = .true.                        !Set flag for updating params
      Scatter   = .false.                       !Unset flag for scattering
   endif
endif

!Check if photon crosses xz-plane in positive direction
if (X_f(2) .ge. yp) then
   X_cut(1) = X_i(1) + (yp - X_i(2)) * nhat_i(1)/nhat_i(2) !Intersection
   X_cut(2) = yp                                           !with
   X_cut(3) = X_i(3) + (yp - X_i(2)) * nhat_i(3)/nhat_i(2) !plane

   if (X_cut(1) .lt. xp .and. &                 !If so, check if it happened
       X_cut(1) .ge. xm .and. &                 !through face of cell
       X_cut(3) .lt. zp .and. &
       X_cut(3) .ge. zm)        then

      X_f       = X_cut                         !Place photon at face
      DirVec    = (/0,1,0/)                     !Continue in neighboring cell
      LeaveCell = .true.                        !Set flag for updating params
      Scatter   = .false.                       !Unset flag for scattering
   endif
endif

!Check if photon crosses xz-plane in negative direction
if (X_f(2) .lt. ym) then
   X_cut(1) = X_i(1) + (ym - X_i(2)) * nhat_i(1)/nhat_i(2) !Intersection
   X_cut(2) = ym                                           !with
   X_cut(3) = X_i(3) + (ym - X_i(2)) * nhat_i(3)/nhat_i(2) !plane

   if (X_cut(1) .lt. xp .and. &                 !If so, check if it happened
       X_cut(1) .ge. xm .and. &                 !through face of cell
       X_cut(3) .lt. zp .and. &
       X_cut(3) .ge. zm)        then

      X_f       = X_cut                         !Place photon at face
      DirVec    = (/0,-1,0/)                    !Continue in neighboring cell
      LeaveCell = .true.                        !Set flag for updating params
      Scatter   = .false.                       !Unset flag for scattering
   endif
endif

!Check if photon crosses xy-plane in positive direction
if (X_f(3) .ge. zp) then
   X_cut(1) = X_i(1) + (zp - X_i(3)) * nhat_i(1)/nhat_i(3) !Intersection
   X_cut(2) = X_i(2) + (zp - X_i(3)) * nhat_i(2)/nhat_i(3) !with
   X_cut(3) = zp                                           !plane

   if (X_cut(1) .lt. xp .and. &                 !If so, check if it happened
       X_cut(1) .ge. xm .and. &                 !through face of cell
       X_cut(2) .lt. yp .and. &
       X_cut(2) .ge. ym)        then

      X_f       = X_cut                         !Place photon at face
      DirVec    = (/0,0,1/)                     !Continue in neighboring cell
      LeaveCell = .true.                        !Set flag for updating params
      Scatter   = .false.                       !Unset flag for scattering
   endif
endif

!Check if photon crosses xy-plane in negative direction
if (X_f(3) .lt. zm) then
   X_cut(1) = X_i(1) + (zm - X_i(3)) * nhat_i(1)/nhat_i(3) !Intersection
   X_cut(2) = X_i(2) + (zm - X_i(3)) * nhat_i(2)/nhat_i(3) !with
   X_cut(3) = zm                                           !plane

   if (X_cut(1) .lt. xp .and. &                 !If so, check if it happened
       X_cut(1) .ge. xm .and. &                 !through face of cell
       X_cut(2) .lt. yp .and. &
       X_cut(2) .ge. ym)        then

      X_f       = X_cut                         !Place photon at face
      DirVec    = (/0,0,-1/)                    !Continue in neighboring cell
      LeaveCell = .true.                        !Set flag for updating params
      Scatter   = .false.                       !Unset flag for scattering
   endif
endif

end subroutine CheckIntersect

!------------------------------------------------------------------------------

subroutine CrossProduct(A,B,C)

use kinds

!Calculates as C the cross product of two vectors A and B

implicit none
real(dp),dimension(3),intent(in)::  A,B
real(dp),dimension(3),intent(out):: C

C(1) = A(2)*B(3) - A(3)*B(2)
C(2) = A(3)*B(1) - A(1)*B(3)
C(3) = A(1)*B(2) - A(2)*B(1)

end subroutine CrossProduct

!------------------------------------------------------------------------------

subroutine Dipole(nhat_i,nhat_1,nhat_2,nhat_f)

!Randomly selects a unit vector for the outgoing photon according to a dipole
!distribution. nhat_i is a unit vector in the direction of the incident photon,
!and nhat_1 and nhat_2 are perpendicular unit vectors.

use kinds
use GlobalData
use iMonteCarlo

implicit none
real(dp),dimension(3):: nhat_i,nhat_1,nhat_2,nhat_f
real(dp)::              rsq,R1,R2(2),fact,mu,n_1,n_2,frevert
real(dp),parameter::    ot = 1./3.,tt = 2./3.

call MonteCarlo(R1)                             !Generate cosine of polar
fact = 2.- 4.*R1 + sqrt(5.+16.*(-1+R1)*R1)      !scattering angle (in atom
mu   = (-1 + fact**tt) / (fact**ot)             !frame)
if (mu .lt. -.9999999d0) mu = -.9999999d0       !Avoid floating exception
if (mu .gt.  .9999999d0) mu =  .9999999d0       !in frevert

do
   call MonteCarlo(R2)                          !Make random values along
   R2 = 2*R2 - 1                                !nhat_1 and nhat_2
   rsq = sum(R2**2)                             !Use r only if |r|<1 in order
   if (rsq .le. 1.) exit                        !not to favorize corners of
enddo                                           !square

frevert = sqrt((1. - mu**2) / rsq)              !Make x,y-components touch
n_1     = frevert * R2(1)                       !unit sphere
n_2     = frevert * R2(2)

nhat_f  = n_1*nhat_1 + n_2*nhat_2 + mu*nhat_i

end subroutine Dipole

!-------------------------------------------------------------------------------

subroutine EffectiveCell(X_pos)

!Find boundaries and "radius" of effective cell, i.e. the cell to be used for
!the Neufeld solution. The photon is located at the center of this cell.

use kinds
use CurrentCellParameters
use DerivedParameters

implicit none
real(dp), intent(in)::  X_pos(3)
real(dp)                d_face(6)

d_face(1) = xp - X_pos(1)                       !
d_face(2) = X_pos(1) - xm                       !
d_face(3) = yp - X_pos(2)                       !Distances from photon to
d_face(4) = X_pos(2) - ym                       !edges of host cell
d_face(5) = zp - X_pos(3)                       !
d_face(6) = X_pos(3) - zm                       !

Adjac     = sum(minloc(d_face))                 !Index of adjacent face; 'sum' to convert array to scalar
R_cube    = .99 * d_face(Adjac)                 !Distance from center to edge
x_adj     = Faces(Adjac)                        !x-, y-, OR z-value of adjacent face

! Faces_eff(1) = X_pos(1) + R_cube                !xp_eff
! Faces_eff(2) = X_pos(1) - R_cube                !xm_eff
! Faces_eff(3) = X_pos(2) + R_cube                !yp_eff
! Faces_eff(4) = X_pos(2) - R_cube                !ym_eff
! Faces_eff(5) = X_pos(3) + R_cube                !zp_eff
! Faces_eff(6) = X_pos(3) - R_cube                !zm_eff

! Faces_eff(Adjac) = Faces(Adjac)     !Make sure Adjac face coincide _exactly_ with face of host cell

end subroutine EffectiveCell

!------------------------------------------------------------------------------

subroutine EnterNeighbor(CurrentCell,X_cut,DirVec,i,j,k,escape)

!Locate neighboring cell.

use kinds
use CellStructure
use GlobalData
use PhysicalConstants
use DataArrays
use iLocateHostCell
use iPrintCell

implicit none
type(Cell), pointer::   CurrentCell             !intent(INOUT)
real(dp),intent(in)::   X_cut(3)
integer, intent(in)::   DirVec(3)
integer,intent(inout):: i,j,k
logical,intent(out)::   escape
real(dp)::              peek(3)
integer::               NeighbInd(3)

! NeighbInd = CurrentCell%indexico + DirVec

! if (ALL(NeighbInd.ge.(/1,1,1/) .and. NeighbInd.le.(/2,2,2/))) then
!   CurrentCell => CurrentCell%Parent%Child(NeighbInd(1), &
!                                           NeighbInd(2), &
!                                           NeighbInd(3))
! else
    peek = X_cut + DirVec*dxTiny
    call LocateHostCell(peek,CurrentCell,i,j,k,escape)
! endif

end subroutine EnterNeighbor

!------------------------------------------------------------------------------

subroutine ExposeCCDs(x,X_pos,nhat_i,ScatCell,FreeEsc)

use kinds
use AtomParameters
use CellStructure
use CurrentCellParameters
use PhysicalConstants
use ObservationalParameters
use iPropagate
use iPrintCell

implicit none
real(dp),intent(in):: x
real(dp),intent(in):: X_pos(3),nhat_i(3)
logical,intent(in)::  FreeEsc
real(dp)::            X_start(3),xx,yy,zz
type(cell),pointer::  ScatCell                 !intent IN!
type(cell),pointer::  StartCell
real(dp)::            x_cell,tau_esc,sigma,XsecD
logical::             inxtube,inytube,inztube
integer::             DirVec(3)

xx      = X_pos(1)
yy      = X_pos(2)
zz      = X_pos(3)
inxtube = yy.gt.o1 .and. yy.lt.o2 .and. &
          zz.gt.o1 .and. zz.lt.o2
inytube = xx.gt.o1 .and. xx.lt.o2 .and. &
          zz.gt.o1 .and. zz.lt.o2
inztube = xx.gt.o1 .and. xx.lt.o2 .and. &
          yy.gt.o1 .and. yy.lt.o2

if (WhichCCDs(1) .and. inxtube) then
  if (FreeEsc) then
    x_cell = x                !Photon already has Doppler shift from emitting atom
  else
    x_cell  = x - u_II + U(1) !x of photon going in xp-dir IN CELL REFERENCE!!!
  endif
  tau_esc = (xp - xx) * (n_HI*sigma(Dnu_D,x_cell) + n_d*XsecD(Dnu_D,x_cell))
  if (tau_esc .lt. tau_max) then
    StartCell => ScatCell
    X_start   = (/xp,yy,zz/)
    DirVec    = (/1,0,0/)
    call Propagate(x_cell,X_start,DirVec,nhat_i,tau_esc,StartCell,FreeEsc)
  endif
endif

if (WhichCCDs(2) .and. inxtube) then
  if (FreeEsc) then
    x_cell = x
  else
    x_cell  = x - u_II - U(1)
  endif
  tau_esc = (xx - xm) * (n_HI*sigma(Dnu_D,x_cell) + n_d*XsecD(Dnu_D,x_cell))
  if (tau_esc .lt. tau_max) then
    StartCell => ScatCell
    X_start   = (/xm,yy,zz/)
    DirVec    = (/-1,0,0/)
    call Propagate(x_cell,X_start,DirVec,nhat_i,tau_esc,StartCell,FreeEsc)
  endif
endif

if (WhichCCDs(3) .and. inytube) then
  if (FreeEsc) then
    x_cell = x
  else
    x_cell  = x - u_II + U(2)
  endif
  tau_esc = (yp - yy) * (n_HI*sigma(Dnu_D,x_cell) + n_d*XsecD(Dnu_D,x_cell))
  if (tau_esc .lt. tau_max) then
    StartCell => ScatCell
    X_start   = (/xx,yp,zz/)
    DirVec    = (/0,1,0/)
    call Propagate(x_cell,X_start,DirVec,nhat_i,tau_esc,StartCell,FreeEsc)
  endif
endif

if (WhichCCDs(4) .and. inytube) then
  if (FreeEsc) then
    x_cell = x
  else
    x_cell  = x - u_II - U(2)
  endif
  tau_esc = (yy - ym) * (n_HI*sigma(Dnu_D,x_cell) + n_d*XsecD(Dnu_D,x_cell))
  if (tau_esc .lt. tau_max) then
    StartCell => ScatCell
    X_start   = (/xx,ym,zz/)
    DirVec    = (/0,-1,0/)
    call Propagate(x_cell,X_start,DirVec,nhat_i,tau_esc,StartCell,FreeEsc)
  endif
endif

if (WhichCCDs(5) .and. inztube) then
  if (FreeEsc) then
    x_cell = x
  else
    x_cell  = x - u_II + U(3)
  endif
  tau_esc = (zp - zz) * (n_HI*sigma(Dnu_D,x_cell) + n_d*XsecD(Dnu_D,x_cell))
  if (tau_esc .lt. tau_max) then
    StartCell => ScatCell
    X_start   = (/xx,yy,zp/)
    DirVec    = (/0,0,1/)
    call Propagate(x_cell,X_start,DirVec,nhat_i,tau_esc,StartCell,FreeEsc)
  endif
endif

if (WhichCCDs(6) .and. inztube) then
  if (FreeEsc) then
    x_cell = x
  else
    x_cell  = x - u_II - U(3)
  endif
  tau_esc = (zz - zm) * (n_HI*sigma(Dnu_D,x_cell) + n_d*XsecD(Dnu_D,x_cell))
  if (tau_esc .lt. tau_max) then
    StartCell => ScatCell
    X_start   = (/xx,yy,zm/)
    DirVec    = (/0,0,-1/)
    call Propagate(x_cell,X_start,DirVec,nhat_i,tau_esc,StartCell,FreeEsc)
  endif
endif

end subroutine ExposeCCDs

!------------------------------------------------------------------------------

subroutine gasdev_s(harvest)

!Return a random number from a Gaussian distribution with a standard deviation of 1.

use kinds
use iMonteCarlo
implicit none
real(sp), intent(out):: harvest
real(sp)::              rsq,v(2)
real(sp), save::        G
logical, save::         G_old=.false.

if (G_old) then
  harvest = G
  G_old   = .false.
else
  do
    call MonteCarlo(v)
    v   = 2.0_sp * v - 1.0_sp
    rsq = sum(v**2)
    if (rsq .gt. 0.0 .and. rsq .lt. 1.0) exit
  enddo
  rsq     = sqrt(-2.0_sp*log(rsq)/rsq)
  harvest = v(1) * rsq
  G       = v(2) * rsq
  G_old   = .true.
endif

end subroutine gasdev_s
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

subroutine gasdev_sv(harvest)
use kinds
use iArrayCopy
use iMonteCarlo
implicit none
real(sp), dimension(:), intent(out)::       harvest
real(sp), dimension(size(harvest))::        rsq,v1,v2
real(sp), allocatable, dimension(:), save:: G
integer(i4b)::                              n,ng,nn,m
integer(i4b), save::                        last_allocated=0
logical, save::                             G_old=.false.
logical, dimension(size(harvest))::         mask

n = size(harvest)
if (n .ne. last_allocated) then
  if (last_allocated .ne. 0) deallocate(G)
  allocate(G(n))
  last_allocated = n
  G_old          = .false.
endif
if (G_old) then
  harvest = G
  G_old   = .false.
else
  ng = 1
  do
    if (ng .gt. n) exit
    call MonteCarlo(v1(ng:n))
    call MonteCarlo(v2(ng:n))
    v1(ng:n)   = 2.0_sp * v1(ng:n) - 1.0_sp
    v2(ng:n)   = 2.0_sp * v2(ng:n) - 1.0_sp
    rsq(ng:n)  = v1(ng:n)**2 + v2(ng:n)**2
    mask(ng:n) = (rsq(ng:n).gt.0.0 .and. rsq(ng:n).lt.1.0)
    call ArrayCopy(pack(v1(ng:n),mask(ng:n)),v1(ng:),nn,m)
    v2(ng:ng+nn-1)  = pack(v2(ng:n),mask(ng:n))
    rsq(ng:ng+nn-1) = pack(rsq(ng:n),mask(ng:n))
    ng              = ng + nn
  enddo
  rsq     = sqrt(-2.0_sp * log(rsq)/rsq)
  harvest = v1 * rsq
  G       = v2 * rsq
  G_old   = .true.
endif

end subroutine gasdev_sv
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

subroutine gasdev_d(harvest)
use kinds
use iMonteCarlo
implicit none
real(dp), intent(out):: harvest
real(dp)::              rsq,v(2)
real(dp), save::        G
logical, save::         G_old=.false.

if (G_old) then
  harvest      = G
  G_old = .false.
else
  do
    call MonteCarlo(v)
    v   = 2.0_dp * v  - 1.0_dp
    rsq = sum(v**2)
    if (rsq .gt. 0.0 .and. rsq .lt. 1.0) exit
  enddo
  rsq     = sqrt(-2.0_dp*log(rsq)/rsq)
  harvest = v(1) * rsq
  G       = v(2) * rsq
  G_old   = .true.
endif

end subroutine gasdev_d
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

subroutine gasdev_dv(harvest)

use kinds
use iArrayCopy
use iMonteCarlo
implicit none
real(dp), dimension(:), intent(out)::       harvest
real(dp), dimension(size(harvest))::        rsq,v1,v2
real(dp), allocatable, dimension(:), save:: G
integer(i4b)::                              n,ng,nn,m
integer(i4b), save::                        last_allocated=0
logical, save::                             G_old=.false.
logical, dimension(size(harvest))::         mask

n = size(harvest)
if (n .ne. last_allocated) then
  if (last_allocated .ne. 0) deallocate(G)
  allocate(G(n))
  last_allocated = n
  G_old = .false.
endif
if (G_old) then
  harvest      = G
  G_old = .false.
else
  ng = 1
  do
    if (ng .gt. n) exit
    call MonteCarlo(v1(ng:n))
    call MonteCarlo(v2(ng:n))
    v1(ng:n)   = 2.0_dp * v1(ng:n) - 1.0_dp
    v2(ng:n)   = 2.0_dp * v2(ng:n) - 1.0_dp
    rsq(ng:n)  = v1(ng:n)**2 + v2(ng:n)**2
    mask(ng:n) = (rsq(ng:n).gt.0.0 .and. rsq(ng:n).lt.1.0)
    call ArrayCopy(pack(v1(ng:n),mask(ng:n)),v1(ng:),nn,m)
    v2(ng:ng+nn-1)  = pack(v2(ng:n),mask(ng:n))
    rsq(ng:ng+nn-1) = pack(rsq(ng:n),mask(ng:n))
    ng              = ng + nn
  enddo
  rsq     = sqrt(-2.0_dp * log(rsq)/rsq)
  harvest = v1 * rsq
  G       = v2 * rsq
  G_old   = .true.
endif

end subroutine gasdev_dv

!------------------------------------------------------------------------------

subroutine getlun(u)

!Return an available logical unit number, starting from 10, since I often use
!single-digit numbers during debugging without openening the files first.

implicit none
integer,intent(out):: u
integer::             i
logical::             aaben

do i = 10,101
  u = i
  inquire(u, opened=aaben)
  if (.not.aaben) exit
enddo

if (u .eq. 101) stop "No more luns"

end subroutine getlun

!------------------------------------------------------------------------------

subroutine GoToCell(CurrentCell)

!Returns a pointer to the cell specified by 'FullPathArray'. Intended for
!debugging.

use CellStructure
use GlobalData
use DataArrays

implicit none
type(Cell), pointer:: CurrentCell       !intent(IN)
integer::             ii,jj,kk,L

CurrentCell => BaseGrid

do L = 0,L_max-1
  if (FullPathArray(1,L) .eq. 0) exit
  if (FullPathArray(1,L) .lt. 1 .or. FullPathArray(1,L) .gt. ni) STOP "Rahsgu!"

  ii = FullPathArray(1,L)
  jj = FullPathArray(2,L)
  kk = FullPathArray(3,L)

  CurrentCell => CurrentCell%Child(ii,jj,kk)
enddo

end subroutine GoToCell

!------------------------------------------------------------------------------

subroutine GetCellData(HostCell,DirVec)

use kinds
use CurrentCellParameters
use DataArrays
use GlobalData
use DerivedParameters
use PhysicalConstants

implicit none
type(cell), pointer:: HostCell
integer, intent(IN):: DirVec(3)
real(dp)::              loga,XsecD,sigma,sigma_0,CoreSkip
integer::             L

n_HI    = dble(HostCell%n_HI)
Dnu_D   = dble(HostCell%Dnu_D)
U_bulk  = dble(HostCell%U_bulk)

# ifdef AMR
L = HostCell%Level
# else
L = 0
# endif

# ifdef dust
n_d = dble(HostCell%n_d)
# else
n_d = 0d0
# endif

# ifdef multi
  clnum = HostCell%clnum
# endif

if (recoil) g_rec = (nu_0/c)**2*h_Pl/(m_H*Dnu_D)!Recoil factor (Field 1959, Adams 1971), approximated by nu -> nu_0
sigma_0 = sigma(Dnu_D,0d0)                      !Line center hydrogen cross section
tau_0   = dx(L) * n_HI * sigma_0                !Across WHOLE cell
a       = Dnu_L / (2 * Dnu_D)                   !Damping parameter
at      = .5d0 * a * tau_0                      !From CENTER to face of cell
loga    = log10(a)
x_cw    = 1.59 - 0.60*loga - 0.03*loga**2       !Core/wing transition

if (intxc) x_crit = CoreSkip(at)                !Threshold for core-skipping acceleration scheme
if (maxxc) x_crit = min(x_crit,xc_max)

# ifdef dust
tau_a = dxH(L) * n_d * XsecD(Dnu_D,0d0) * (1. - albedo)!Dust abs. from CENTER to face of cell

if (tau_a .gt. 0d0) then
  logatr = log10(at/tau_a)                      !For recovering core-history
else
  logatr = 2 * atA(10)
endif

if (logatr .lt. atA(1))  logatr = atA(1)  + 1e-6
# endif

if (DirVec(1) .eq. 1) then                      !This block ensures that
  xm = xp                                       !adjacent cell faces have
  xp = HostCell%C_pos(1) + dxH(L)               !exactly the same coordinates.
  yp = HostCell%C_pos(2) + dyH(L)               !If not, then photons may get
  ym = HostCell%C_pos(2) - dyH(L)               !stuck between cells.
  zp = HostCell%C_pos(3) + dzH(L)
  zm = HostCell%C_pos(3) - dzH(L)
elseif (DirVec(1) .eq. -1) then
  xp = xm
  xm = HostCell%C_pos(1) - dxH(L)
  yp = HostCell%C_pos(2) + dyH(L)
  ym = HostCell%C_pos(2) - dyH(L)
  zp = HostCell%C_pos(3) + dzH(L)
  zm = HostCell%C_pos(3) - dzH(L)
elseif (DirVec(2) .eq. 1)  then
  xp = HostCell%C_pos(1) + dxH(L)
  xm = HostCell%C_pos(1) - dxH(L)
  ym = yp
  yp = HostCell%C_pos(2) + dyH(L)
  zp = HostCell%C_pos(3) + dzH(L)
  zm = HostCell%C_pos(3) - dzH(L)
elseif (DirVec(2) .eq. -1) then
  xp = HostCell%C_pos(1) + dxH(L)
  xm = HostCell%C_pos(1) - dxH(L)
  yp = ym
  ym = HostCell%C_pos(2) - dyH(L)
  zp = HostCell%C_pos(3) + dzH(L)
  zm = HostCell%C_pos(3) - dzH(L)
elseif (DirVec(3) .eq. 1)  then
  xp = HostCell%C_pos(1) + dxH(L)
  xm = HostCell%C_pos(1) - dxH(L)
  yp = HostCell%C_pos(2) + dyH(L)
  ym = HostCell%C_pos(2) - dyH(L)
  zm = zp
  zp = HostCell%C_pos(3) + dzH(L)
elseif (DirVec(3) .eq. -1) then
  xp = HostCell%C_pos(1) + dxH(L)
  xm = HostCell%C_pos(1) - dxH(L)
  yp = HostCell%C_pos(2) + dyH(L)
  ym = HostCell%C_pos(2) - dyH(L)
  zp = zm
  zm = HostCell%C_pos(3) - dzH(L)
elseif (DirVec(1) .eq. 7) then                  !If photon does not come from
  xp = HostCell%C_pos(1) + dxH(L)               !neighboring cell (i.e. at
  xm = HostCell%C_pos(1) - dxH(L)               !injection).
  yp = HostCell%C_pos(2) + dyH(L)
  ym = HostCell%C_pos(2) - dyH(L)
  zp = HostCell%C_pos(3) + dzH(L)
  zm = HostCell%C_pos(3) - dzH(L)
else
  stop "WTF?!"
endif

Faces  = (/xp,xm,yp,ym,zp,zm/)

end subroutine GetCellData

!------------------------------------------------------------------------------

subroutine HenyeyGreenstein(nhat_i,g,nhat_f)

use kinds
use iMonteCarlo

!Randomly selects a unit vector nhat_f for the outgoing photon according to the
!anisotropic distribution given by Henyey & Greenstein (1941).
!nhat_i is a unit vector in the direction of the incident photon, and nhat_f is
!that of the outgoing photon.
!g = <cos(theta)> = <mu> is the average of the deflection angle from nhat_i
!(the asymmetry parameter).
!
!TODO: For |g| >~ 0.995, corresponding to ~5 degrees, the algorithm becomes
!extremely slow, so here g is approximated as 1, i.e. nhat_f = nhat_i. I can't
!think of any realistic physical conditions where this will be a problem, but
!for the sake of beauty, as well as for various test cases, this should be
!fixed. Since the "cap" suspended by the ~5-degree-cone is very close to being
!flat, it will probably suffice to draw a random point in a circular plane of
!radius = sin(theta) = sin(acos(mu)).

implicit none
real(dp),intent(in)::  nhat_i(3),g
real(dp),intent(out):: nhat_f(3)
real(dp)::             nhat(3)
real(dp)::             Pmax,P,R1,R2(2),mu,n1,n2,frevert,rsq

if (abs(g) .lt. .995) then
  Pmax = .5 * (1 - g**2) / (1 + g**2 - 2*g)**1.5

  do
    call MonteCarlo(R1)
    mu = (2*R1 - 1)!* Pmax
    P  = .5 * (1 - g**2) / (1 + g**2 - 2*g*mu)**1.5

    call MonteCarlo(R1)
    if (R1 .le. P/Pmax) then                    !Accept mu if randomly chosen
      do                                        !point under comparison function
        call MonteCarlo(R2)                     !lies also under P(mu)
        R2  = 2*R2 - 1
        rsq = sum(R2**2)
        if (rsq .le. 1d0) exit                  !Accept nhat if inside unit sphere
      enddo

      frevert = sqrt((1d0 - mu**2) / rsq)       !Make x,y-components touch
      n1      = frevert * R2(1)                 !unit sphere
      n2      = frevert * R2(2)

      nhat    = (/n1,n2,mu/)                    !Unit vector of outgoing photon
                                                !in ref.frame where nhat_i lies
                                                !along the z-axis

      call RotateVector(nhat_i,nhat,nhat_f)     !Transform to global coordinates
      exit
    endif
  enddo
else
  nhat_f = sign(nhat_i,g)
endif

end subroutine HenyeyGreenstein

!------------------------------------------------------------------------------

subroutine holdon(t)

!Pause execution for t seconds. NB: Argument is single precision!

use kinds

implicit none
real(sp),intent(in):: t
real(sp)::            t1,t2

call cpu_time(t1)

do
  call cpu_time(t2)
  if (t2-t1 .gt. t) exit
enddo

end subroutine holdon

!------------------------------------------------------------------------------

subroutine InitializeBaseGrid

use CellStructure
use GlobalData
use DataArrays

implicit none
integer:: i,j,k

allocate(BaseGrid%Child(ni,nj,nk))

do i = 1, ni
  do j = 1, nj
    do k = 1, nk
#     ifdef AMR
      BaseGrid%Child(i,j,k)%Refined = .false.
      BaseGrid%Child(i,j,k)%Level   = 0
#     endif
      BaseGrid%Child(i,j,k)%C_pos   = ((/i,j,k/) - .5d0) * (/dx0,dy0,dz0/)
      nullify(BaseGrid%Child(i,j,k)%Child)
#     ifdef multi
      BaseGrid%Child(i,j,k)%phase   = 0
#     endif
    enddo
  enddo
enddo

end subroutine InitializeBaseGrid

!------------------------------------------------------------------------------

subroutine InitDir(nhat_i)

use kinds
use GlobalData

implicit none
real(dp),intent(out):: nhat_i(3)

if (isoem) then
  call IsoDir(nhat_i)
else
  nhat_i = nhat_0
endif

end subroutine InitDir

!------------------------------------------------------------------------------

subroutine InitFreq(Dnu_D,nhat_i,x,photon,U_bulk)

!Initial frequency. In the case of photons produced in the ISM (Lya), the
!frequency is according to random thermal motion of atom AND random natural
!broadening (i.e. a Voigt profile). In the case of continuum photons from stars,
!it is according to the given spectrum (at the moment only a flat spectrum in
!nu).

use kinds
use AtomParameters
use PhysicalConstants
use DerivedParameters
use ObservationalParameters
use GlobalData
use igasdev
use iv2u
use iMonteCarlo

implicit none
real(dp),intent(in)::  Dnu_D,nhat_i(3),U_bulk(3)
integer,intent(in)::   photon
real(dp),intent(out):: x
real(dp)::             R,R3(3),nu,v,lam

trueFUV = .false.

if (photon .eq. 1) then                         !Lya photon
  if (trim(x_injType) .eq. 'proper') then
    call gasdev(R3)
    U    = R3 / sqrt(2d0)                       !=> f(x) = exp(-x^2) / sqrt(pi)
    u_II = dot_product(U,nhat_i)
    call MonteCarlo(R)
    x = a*tan(pi*R + pi/a) + u_II               !=> Voigt
  elseif (trim(x_injType) .eq. 'fixtgauss') then
    call gasdev(R)
    v = R * sigma_V ! DON'T divide by sqrt(2)
    x = v2u(v,Dnu_D)
  elseif (trim(x_injType) .eq. 'x_inj') then
    x = x_inj
  else
    write(*,*) 'x_injType error in "InitFreq":', x_injType
    stop
  endif
elseif (photon .eq. 2) then                         !FUV photon
  call MonteCarlo(R)
  nu  = nu1 + R * (nu2-nu1)
  x   = (nu - nu_0) / Dnu_D
  lam = c/nu / 1d8
  if (photon.eq.2 .and. abs(lam-lambda_0).gt.NBaper) then
    n_FUV   = n_FUV + 1d0
    trueFUV = .true.
  endif
else
  stop 'Unknown photon'
endif

end subroutine InitFreq

!------------------------------------------------------------------------------

subroutine InitPos(n,photon,X_pos)

!Find initial position of a photon

use kinds
use GlobalData
use ObservationalParameters
use DataArrays
use iMonteCarlo
use PhysicalConstants
use iLocateHostCell
use igasdev

implicit none
integer, intent(in)::   n
integer,intent(in)::    photon
real(dp), intent(out):: X_pos(3)
integer::               nn,idum,jdum,kdum,EmPhase
real(dp)::              XwrtC(3),R3(3),R2(2),R1,rsq,P,left,right,norm,d
real(dp),parameter::    shellcorr = 0.0     !Relative prob. of emission in cone
real(dp),parameter::    ICMcorr   = 1.      !Cloud correlation factor (if pos. is NOT a cloud, a random number must exeed this number, OR ELSE...)
logical::               escdum,okay
type(Cell), pointer::   TempCell

if (empos(1)) then                              !X_init depends on cell luminosities
  nn = n_packets(photon)
  if (photon .eq. 1) then
    X_pos = InitLya(nn,1:3)
    if (emdimL .eq. 3) then
      w_ph = 1d0
    else
      w_ph = InitLya(nn,4)
    endif
  elseif (photon .eq. 2) then
    X_pos = InitFUV(nn,1:3)
    if (emdimF .eq. 3) then
      w_ph = 1d0
    else
      w_ph = InitFUV(nn,4)
    endif
  endif
elseif (empos(2)) then                          !Photons are emitted centrally
  call CentralEmission(X_pos)
  w_ph = 1d0
elseif (empos(3)) then                          !Photons are emitted homogeneously IN SPHERE
  do
    do
      call MonteCarlo(R3)
      R3  = 2*R3 - 1
      rsq = sum(R3**2)
      if (rsq.le.1) exit
    enddo
    X_pos = R_box*R3 + R_box

#   ifdef multi
    call LocateHostCell(X_pos,TempCell,idum,jdum,kdum,escdum)
    if (TempCell%phase .eq. 0) then             !Void
      cycle
    elseif (TempCell%phase .eq. 1) then         !Cloud
      exit
    elseif (TempCell%phase .eq. 2) then         !ICM
      call MonteCarlo(R1)
      if (R1 .le. ICMcorr) exit
    elseif (TempCell%phase .eq. 3) then         !Shell/jet
      call MonteCarlo(R1)
      if (R1 .le. shellcorr) exit
    else
      stop 'Not ready for phase > 3.'
    endif
    cycle
#   else
    exit ! If not multi-phase, empos(3) just returns random point in full sphere/box
#   endif
  enddo
  w_ph = 1d0
elseif (empos(4)) then                          !Photons are emitted according
  do                                            !to user-given location
    do
      call MonteCarlo(R3)
      R3  = 2*R3 - 1
      rsq = sum(R3**2)
      if (rsq .le. 1) exit
    enddo
    R3    = R3 / sqrt(rsq)
    X_pos = X_i0 + R3*dxWeeny
    if (all(X_pos.gt.0d0) .and. all(X_pos.lt.D_box)) exit
  enddo
  w_ph = 1d0
elseif (empos(5)) then                          !Exponentially decreasing emission probability
  call MonteCarlo(R1)
  if (CloudCorr .gt. -.5) EmPhase = merge(1,2,R1 .lt. CloudCorr)!Photon is emitted from cloud with P = CloudCorr

  do
    do
      call gasdev(R3)
      X_pos = R3 * H_lum                          !Position wrt. center
      if (sum((X_pos/(axes*r_gal))**2) .lt. 1) exit
    enddo
    X_pos = X_pos + R_box

#   ifdef multi
      call LocateHostCell(X_pos,TempCell,idum,jdum,kdum,escdum)
      if (CloudCorr.lt.-.5 .or. TempCell%phase .eq. EmPhase) exit
#   else
      exit ! If not multi-phase, empos(5) just returns random point inside r_gal
#   endif
  enddo

  w_ph  = 1d0
else
  stop "From where should the photons be emitted?"
endif

n_eff(photon) = n_eff(photon) + w_ph

if (mod(n,wrtph) .eq. 0) write(*,*) 'Launched photon #',n,'of',n_phtot
flush(6)

end subroutine InitPos

!------------------------------------------------------------------------------

subroutine InitTau(tau)

!Randomly selects the optical depth that a photon will reach, weighted with an
!exponential probability density.

use kinds
use iMonteCarlo

implicit none
real(dp):: R,tau

call MonteCarlo(R)
tau = -log(R)

end subroutine InitTau

!------------------------------------------------------------------------------

subroutine Int2Char(i,a,LZ)

implicit none
integer, intent(IN)::           i
character(LEN=*), intent(OUT):: a
logical, intent(IN)::           LZ              !Leading zeros?
integer, parameter::            zero = iachar('0')
integer::                       n, N_digits, length
character(LEN=len(a))::         b

a        = ''
b        = ''
length   = len(a)
N_digits = 0

do n = 1,length
  N_digits = N_digits + 1
  a = achar(zero + mod(i,10**n)/10**(n-1)) // trim(a)
  if (i .lt. 10**n) exit
enddo

if (LZ) then
  b = repeat('0',len(a) - N_digits)
  a = trim(b) // trim(a)
endif

end subroutine Int2Char

!------------------------------------------------------------------------------

subroutine IsoDir(nhat_i)

!Randomly selects a unit vector with isotropical distribution

use kinds
use GlobalData
use iMonteCarlo

implicit none
real(dp),intent(out):: nhat_i(3)
real(dp)::             R3(3),rsq

do
  call MonteCarlo(R3)                           !Make random x-, y- and z-
  R3  = 2*R3 - 1                                !values between -1 and 1
  rsq = sum(R3**2)
  if (rsq .le. 1.) exit                         !Use r only if |r|<1 in order
enddo                                           !not to favor corners of box

nhat_i = R3 / sqrt(rsq)                         !Normalize to unit vector

end subroutine IsoDir

!------------------------------------------------------------------------------

subroutine IsoUp(nhat_i)

!Randomly selects a unit vector with semi-isotropical, upward distribution,
!i.e. nhat_i(3) is always >= 0

use kinds
use GlobalData
use iMonteCarlo

implicit none
real(dp),intent(out):: nhat_i(3)
real(dp)::             R3(3),rsq

do
  call MonteCarlo(R3)                           !Make random x-, y- and z-
  R3  = 2*R3 - 1                                !values between -1 and 1
  rsq = sum(R3**2)
  if (rsq .le. 1.) exit                         !Use r only if |r|<1 in order
enddo                                           !not to favor corners of box

R3(3) = abs(R3(3))

nhat_i = R3 / sqrt(rsq)                         !Normalize to unit vector

end subroutine IsoUp

!------------------------------------------------------------------------------

subroutine InverseNeufeld(x)

!Calculate frequency of photon escaping a sufficiently thick cube as a function
!of initial frequency and value of a * tau_0.

use kinds
use DerivedParameters
use GlobalData
use PhysicalConstants
use iMonteCarlo

implicit none
real(dp), intent(inout):: x
real(dp)::                R,base

call MonteCarlo(R)

base = sqrt(54./pi**3)  *  eta*at_eff  *  log(tan(pi*R/2)) + x**3

if (base .lt. 0d0) then
  base = -base
  x    = -base**.33333333
else
  x    =  base**.33333333
endif

end subroutine InverseNeufeld

!------------------------------------------------------------------------------

subroutine LocateHostCell(X_pos,CurrentCell,i,j,k,escape)

use kinds
use CellStructure
use GlobalData
use PhysicalConstants
use DataArrays
use iDigDeeper

implicit none
real(dp), intent(in):: X_pos(3)
type(Cell), pointer::  CurrentCell              !intent(out)
integer, intent(out):: i,j,k
logical,intent(out)::  escape
real(dp)::             d,norm,XwrtC(3)
integer,parameter::    indilow(3) = (/1,1,1/)
integer::              indicur(3)

i = ceiling(X_pos(1) / dx0)                     !
j = ceiling(X_pos(2) / dy0)                     !Indices of mother cell
k = ceiling(X_pos(3) / dz0)                     !
indicur = (/i,j,k/)

if (fullboxRT) then
  if (all(indicur .ge. indilow) .and. &
      all(indicur .le. BaseRes)) then
    escape = .false.                            !Photon is inside box
  else
    escape = .true.                             !Photon is outside box
  endif                                         !so this is final statement
else
  if (trim(model) .ne. 'slab') then
    XwrtC  = X_pos - R_box
    d      = norm(XwrtC)
    escape = (d .ge. R_box)
  else
    d      = abs(X_pos(3)-R_box)
    escape = (d .gt. z_0)
  endif
endif

if (.not. escape) then
  CurrentCell => BaseGrid%Child(i,j,k)
# ifdef AMR
  if (CurrentCell%Refined) call DigDeeper(X_pos,CurrentCell)
# endif
endif

end subroutine LocateHostCell

!------------------------------------------------------------------------------

subroutine LocationOfEscape(xy_esc)

use kinds
use GlobalData
use DerivedParameters
use igasdev

implicit none
real(dp), intent(out)::  xy_esc(2)
real(dp)::               R(2)

do
  call gasdev(R)
  R = R * .55d0                           !sigma_SB ~ 0.55
  if (all(abs(R) .lt. 1d0)) exit
enddo

xy_esc = R * R_cube

end subroutine LocationOfEscape

!------------------------------------------------------------------------------

subroutine LorentzTransform(x,nhat_i,oldCell,newCell)

!Lorentz transform x according to temperature of velocity of new cell

use kinds
use CellStructure

implicit none
real(dp), intent(INOUT):: x
real(dp), intent(IN)::    nhat_i(3)
type(Cell), pointer::   oldCell, newCell
real(dp)::                oldD,oldU(3),newD,newU(3),u1n,u2n

oldD = dble(oldCell%Dnu_D)
oldU = dble(oldCell%U_bulk)
newD = dble(newCell%Dnu_D)
newU = dble(newCell%U_bulk)

u1n  = dot_product(oldU,nhat_i)
u2n  = dot_product(newU,nhat_i)

x    = (x + u1n) * oldD/newD - u2n

end subroutine LorentzTransform

!------------------------------------------------------------------------------

subroutine MakeFlorent

use kinds
use PhysicalConstants
use GlobalData
use DataArrays
use ObservationalParameters
use CurrentCellParameters
use iDoppler

implicit none
integer:: i,j,k,n,p,q
real(dp)::  R,R3(3),nhat(3),dd

Dnu_DString   = Doppler((/1d5, T_cl, T_ICM/))        !Doppler width in couds; Hz
n_HIString(0) = 1e-37                           !Emptiness outside shell
n_HIString(1) = max(n_HI_cl,1e-37)              !Neutral hydrogen density in clouds
n_HIString(2) = max(n_HI_ICM,1e-37)             !Neutral hydrogen density in ICM
# ifdef dust
n_dString(0)  = 1e-37                           !Emptiness outside shell
n_dString(1)  = max(n_HI_cl  * Z_cl/Z0, 1e-37)  !"Dust density" in clouds
n_dString(2)  = max(n_HI_ICM * Z_ICM/Z0,1e-37)  !"Dust density" in ICM
# endif

end subroutine MakeFlorent

!------------------------------------------------------------------------------

subroutine MakeSemireal

use kinds
use PhysicalConstants
use GlobalData
use DataArrays
use ObservationalParameters
use iMonteCarlo
use iDoppler
use igasdev

implicit none
integer:: i,j,k,n,q,cl
real(dp)::  R,R2(2),R3(3),rsq,P,X_cl(3),d_cl,X(3),nhat(3),dd
real(dp)::  tt1,tt2

allocate(Clouds(N_cl,3))
allocate(CloudVelDisp(N_cl,3))

Dnu_DString   = Doppler((/T_sh, T_cl, T_ICM, T_sh/)) !Doppler width in couds; Hz
n_HIString(0) = 1e-37                           !Emptiness outside shell
n_HIString(1) = max(n_HI_cl,1d-37)              !Cloud HI density
n_HIString(2) = max(n_HI_ICM,1d-37)             !ICM HI density, in center
n_HIString(3) = max(n_HI_sh,1d-37)              !Shell HI density
# ifdef dust
n_dString(0)  = 1e-37                           !Emptiness outside shell
n_dString(1)  = max(n_HI_cl  * Z_cl  / Z0,1d-37)!Cloud "dust density"
n_dString(2)  = max(n_HI_ICM * Z_ICM / Z0,1d-37)!ICM "dust density", in center
n_dString(3)  = max(n_HI_sh  * Z_sh  / Z0,1d-37)!Shell "dust density"
# endif

if (N_cl .gt. 0) then
  cl = 1
  do                                            !Put cloud centers in "Clouds"
    do                                          !density gradient
      call MonteCarlo(R3)
      R3  = 2*R3 - 1
      rsq = sum(R3**2)
      if (rsq .le. 1) exit
    enddo
    R3 = R3 / sqrt(rsq)                         !Normalized random unit vector

    do                                          !
      call MonteCarlo(R2)                       !Normalized gradient weighted
      P = exp(-(R2(1)/(H1/R_box))**exp1)                !distance in [0,1] (truncated!)
      if (R2(2) .le. P) exit                    !
    enddo                                       !

    X_cl = R_box + (R3 * R2(1)*r_gal)
    d_cl = sqrt(sum((X_cl-R_box)**2))           !Distance from box center

    if (d_cl .gt. r_gal) cycle                  !Make sure they are in galaxy
  ! Uncomment folowing line  if you don't want the possibility of a central cloud:
  ! if (d_cl .le. r_cl+D_box/ni*sqrt(3.)/2) cycle!...and don't overlap box center
    Clouds(cl,:) = X_cl
    call gasdev(R3)
    CloudVelDisp(cl,:) = R3 * sigV_cl           !Cloud vel. disp. std. dev. in km/s
    cl = cl + 1
    if (cl .gt. N_cl) exit
  enddo
endif

end subroutine MakeSemireal

!------------------------------------------------------------------------------

subroutine MakeShell

use kinds
use PhysicalConstants
use GlobalData
use DataArrays
use ObservationalParameters
use CurrentCellParameters

implicit none
integer:: i,j,k,n,p,q
real(dp)::  R,R3(3),nhat(3),dd

Dnu_DString  = Dnu_D_cl                         !All cells should have same T, in order for x of escaping photons to make sense
n_HIString   = (/1d-37, n_HI_cl/)
# ifdef dust
n_dString    = (/0d0, n_HI_cl*Z_cl/Z0/)
# endif

end subroutine MakeShell

!------------------------------------------------------------------------------

subroutine MakeSlab

use kinds
use PhysicalConstants
use GlobalData
use DataArrays
use ObservationalParameters
use CurrentCellParameters
use iv2u

implicit none
integer:: i,j,k,n,p,q
real(dp)::  R,R3(3),nhat(3),dd
real(dp)::  U_z

allocate(U_zString(0:1))

Dnu_DString  = Dnu_D_cl                         !All cells should have same T, in order for x of escaping photons to make sense
n_HIString  = (/1d-37, n_HI_cl/)
U_z         = v2u(V_z,Dnu_DString(1))
U_zString   = (/0d0, U_z/) ! (/U_z, U_z/)
# ifdef dust
n_dString   = (/0d0, n_HI_cl*Z_cl/Z0/)
# endif

end subroutine MakeSlab

!------------------------------------------------------------------------------

subroutine MonteCarlo_s(harvest)

!Return a random number in [0,1] in (every element of) harvest. Harvest can be
!a single or double precision scalar or 1D-vector.
!We have found no differences except for statistical noise in the final output
!using "only" Fortrans intrinsic subroutine "random_number", but if one wishes,
!it is straightforward to replace it by a higher quality random number
!generator, e.g. Numerical Recipes' "ran1". Note, however, that in this case
!the "true random initialization" in the subroutine "ReadInput" should also be
!changed, if it is to be used.

use kinds
implicit none
real(sp),intent(out):: harvest

call random_number(harvest)

end subroutine MonteCarlo_s
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

subroutine MonteCarlo_sv(harvest)
use kinds
implicit none
real(sp),dimension(:),intent(out):: harvest

call random_number(harvest)

end subroutine MonteCarlo_sv
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

subroutine MonteCarlo_d(harvest)
use kinds
implicit none
real(dp),intent(out):: harvest

call random_number(harvest)

end subroutine MonteCarlo_d
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

subroutine MonteCarlo_dv(harvest)
use kinds
implicit none
real(dp),dimension(:),intent(out):: harvest

call random_number(harvest)

end subroutine MonteCarlo_dv

!------------------------------------------------------------------------------

subroutine NearCentralEmission(dir,X_pos)

!Return position negligibly far from the center of the box in the direction
!given by keyword 'dir' (1=x, 2=y, 3=z), but within a circle of radius
!'HowClose' in the other two directions. This is useful for calculating
!approximate average quantities along an axis, where a little more statistics
!than just AT axis is needed.

use kinds
use GlobalData
use iMonteCarlo

implicit none
integer, intent(in)::   dir
real(dp), intent(out):: X_pos(3)
real(dp)::              R1,R2(2),rsq,HowClose

HowClose = 0.01 * R_box

do
  call MonteCarlo(R2)
  R2   = 2*R2 - 1
  rsq = sum(R2**2)
  if (rsq .le. 1) exit
enddo
R2 = R2 / sqrt(rsq)

call MonteCarlo(R1)
R1 = 2*R1 - 1

select case (dir)
  case (1)
    X_pos(1) = R1    * dxWeeny
    X_pos(2) = R2(1) * HowClose
    X_pos(3) = R2(2) * HowClose
  case (2)
    X_pos(1) = R2(1) * HowClose
    X_pos(2) = R1    * dxWeeny
    X_pos(3) = R2(2) * HowClose
  case (3)
    X_pos(1) = R2(1) * HowClose
    X_pos(2) = R2(2) * HowClose
    X_pos(3) = R1    * dxWeeny
  case default; stop 'Rogntudju!'
end select

X_pos = X_pos + R_box

end subroutine NearCentralEmission

!------------------------------------------------------------------------------

subroutine Neufeld(x,X_pos,nhat_f,absorb)!,dist,scats)

!For a sufficiently thick host cell, this subroutine constructs the largest
!possible cell inside the host cell in which the photon is in the center and
!calculates analytically the characteristics of the escaping photon

use kinds
use CurrentCellParameters
use DerivedParameters
use GlobalData
use PhysicalConstants
use iMonteCarlo

implicit none
real(dp), intent(inout):: x
!real(dp), intent(inout):: dist
!integer(kind=8), intent(inout):: scats
real(dp), intent(inout):: X_pos(3), nhat_f(3)
logical,intent(out)::   absorb
real(dp)::                arg,at3ta
real(dp)::                R,ExitVec(3),xy_esc(2)
real(dp)::                P_esc
integer::               face

at3ta = (at_eff)**.333333 * (tau_a*z_frac) * eta**1.333333
arg   = sqrt(3.) / (zeta*pi**.41666667) * (at3ta)**.55
P_esc = 1d0 / cosh(arg)

call MonteCarlo(R)

if (R .gt. P_esc) then ! R>P => no esc.
  absorb = .true.
else
  absorb = .false.
  call Thomson(ExitVec)
  call InverseNeufeld(x)
  call LocationOfEscape(xy_esc)

  call MonteCarlo(R)
  face = ceiling(6 * R)

  if     (face .eq. 1) then                   !Photon escapes in xp-direction
    nhat_f    =  ExitVec

    X_pos(1)  =  X_pos(1) + R_cube!x_adj
    X_pos(2)  =  X_pos(2) + xy_esc(1)
    X_pos(3)  =  X_pos(3) + xy_esc(2)
  elseif (face .eq. 2) then                   !Photon escapes in xm-direction
    nhat_f(1) = -ExitVec(1)
    nhat_f(2) =  ExitVec(2)
    nhat_f(3) =  ExitVec(3)

    X_pos(1)  =  X_pos(1) - R_cube!x_adj
    X_pos(2)  =  X_pos(2) + xy_esc(1)
    X_pos(3)  =  X_pos(3) + xy_esc(2)
  elseif (face .eq. 3) then                   !Photon escapes in yp-direction
    nhat_f(1) =  ExitVec(2)
    nhat_f(2) =  ExitVec(1)
    nhat_f(3) =  ExitVec(3)

    X_pos(1)  =  X_pos(1) + xy_esc(1)
    X_pos(2)  =  X_pos(2) + R_cube!x_adj
    X_pos(3)  =  X_pos(3) + xy_esc(2)
  elseif (face .eq. 4) then                   !Photon escapes in ym-direction
    nhat_f(1) =  ExitVec(2)
    nhat_f(2) = -ExitVec(1)
    nhat_f(3) =  ExitVec(3)

    X_pos(1)  =  X_pos(1) + xy_esc(1)
    X_pos(2)  =  X_pos(2) - R_cube!x_adj
    X_pos(3)  =  X_pos(3) + xy_esc(2)
  elseif (face .eq. 5) then                   !Photon escapes in zp-direction
    nhat_f(1) =  ExitVec(2)
    nhat_f(2) =  ExitVec(3)
    nhat_f(3) =  ExitVec(1)

    X_pos(1)  =  X_pos(1) + xy_esc(1)
    X_pos(2)  =  X_pos(2) + xy_esc(2)
    X_pos(3)  =  X_pos(3) + R_cube!x_adj
  elseif (face .eq. 6) then                   !Photon escapes in zm-direction
    nhat_f(1) =  ExitVec(2)
    nhat_f(2) =  ExitVec(3)
    nhat_f(3) = -ExitVec(1)

    X_pos(1)  =  X_pos(1) + xy_esc(1)
    X_pos(2)  =  X_pos(2) + xy_esc(2)
    X_pos(3)  =  X_pos(3) - R_cube!x_adj
  else
    stop 'Komnuon!'
  endif
endif

end subroutine Neufeld

!------------------------------------------------------------------------------

subroutine PrintCell(CurrentCell)

!Write cell data to stdout (intended for debugging)

use CellStructure
use GlobalData
use PhysicalConstants
use iCellPath

implicit none
type(Cell), pointer:: CurrentCell              !intent(in)
real(dp)::            norm

call CellPath(CurrentCell)
write(*,*) "Data for '" // trim(FullPathString) // "'"

# ifdef multi
write(*,*) 'Phase:   ', CurrentCell%phase
write(*,*) 'clnum:   ', CurrentCell%clnum
# endif

# ifdef AMR
write(*,*) 'Level:   ', CurrentCell%Level
# endif

write(*,*) 'n_HI:    ', real(CurrentCell%n_HI)

# ifdef dust
write(*,*) '"n_d":   ', real(CurrentCell%n_d)
# endif

write(*,*) 'T:       ', real(m_H/2d0/k_B * (CurrentCell%Dnu_D*c/nu_0)**2)
write(*,*) 'V_bulk:  ', real(CurrentCell%U_bulk*CurrentCell%Dnu_D*c/nu_0/1e5)
write(*,*) 'C_pos:   ', real(CurrentCell%C_pos / kpc)
write(*,*) 'd:       ', real(norm((CurrentCell%C_pos-R_box) / kpc))

write(*,*)

end subroutine PrintCell

!------------------------------------------------------------------------------

subroutine Propagate(x,X_start,DirVec,nhat_i,tau_esc,PropCell,FreeEsc)

use kinds
use CellStructure
use PhysicalConstants
use GlobalData
use ObservationalParameters
use iPrintCell
use iEnterNeighbor
use iLorentzTransform

implicit none
real(dp)::              tau_esc,x,mu,DDD(3)
real(dp), intent(inout)::  X_start(3)
real(dp), intent(in)::  nhat_i(3)
integer,intent(in)::  DirVec(3)
logical,intent(in)::  FreeEsc
type(cell), pointer:: PropCell                  !intent IN!
type(cell), pointer:: PrevCell
real(dp)::              x_cell,n_HI,n_d,Dnu_D,sigma_x,sigma,sigma_d,XsecD
real(dp)::              lambda
integer::             i,j,k,L,m,n,s1,s2
integer::             ax1,ax2
logical::             skip,escape
integer::aaa

skip = .false.

do
  PrevCell => PropCell
  call EnterNeighbor(PropCell,X_start,DirVec,i,j,k,escape)

  if (escape) exit
  DDD = dble(DirVec)
  call LorentzTransform(x,DDD,PrevCell,PropCell)

  n_HI    = PropCell%n_HI
# ifdef dust
  n_d     = PropCell%n_d
# endif
  Dnu_D   = PropCell%Dnu_D
  sigma_x = sigma(Dnu_D,x)

# ifdef AMR
  L       = PropCell%Level
# else
  L       = 0
# endif

# ifdef dust
  sigma_d = XsecD(Dnu_D,x)
  tau_esc = tau_esc + dx(L) * (n_HI*sigma_x + n_d*sigma_d)
# else
  tau_esc = tau_esc + dx(L) * (n_HI*sigma_x)
# endif

  X_start = X_start + DirVec*dx(L)   !Advance position

  if (tau_esc .gt. tau_max) then
    skip = .true.
    exit
  endif
enddo

if (.not. skip) then
  x      = x + dot_product(real(DirVec),PropCell%U_bulk) !x is now in global ref. frame, but with cell's Dnu_D
  lambda = c / (x * PropCell%Dnu_D + nu_0) * 1d8         !lambda in AA in global ref. frame
  s1     = ceiling((lambda - BW(1))/DDlam * SpecRes1D)
  s2     = ceiling((lambda - BW(1))/DDlam * SpecRes2D)
  if (s1.lt.1)         s1 = 1
  if (s1.gt.SpecRes1D) s1 = SpecRes1D
  if (s2.lt.1)         s2 = 1
  if (s2.gt.SpecRes2D) s2 = SpecRes2D

  !The following if statement makes sure that when galaxy is viewed from a given
  !direction, the two spatial axes always have their values increasing from the
  !lower left corner to the right (for the "x" axis) and up (for the "y" axis):
  if (DirVec(1).eq. 1) then; ax1 = 2; ax2 = 3; mu =  nhat_i(1) ; elseif &
     (DirVec(1).eq.-1) then; ax1 = 3; ax2 = 2; mu = -nhat_i(1) ; elseif &
     (DirVec(2).eq. 1) then; ax1 = 3; ax2 = 1; mu =  nhat_i(2) ; elseif &
     (DirVec(2).eq.-1) then; ax1 = 1; ax2 = 3; mu = -nhat_i(2) ; elseif &
     (DirVec(3).eq. 1) then; ax1 = 1; ax2 = 2; mu =  nhat_i(3) ; elseif &
     (DirVec(3).eq.-1) then; ax1 = 2; ax2 = 1; mu = -nhat_i(3) ; else; stop 'Bad DirVec'; endif
  m = ceiling((X_start(ax1)-(R_box-R_obs))/D_obs * pix)
  n = ceiling((X_start(ax2)-(R_box-R_obs))/D_obs * pix)

  if (FreeEsc) then
    w_esc = (2)* .5d0             * exp(-tau_esc) * w_ph
  else
    w_esc = (2)* .375d0*(1+mu**2) * exp(-tau_esc) * w_ph
  endif

  if (DirVec(1).eq. 1) then; specxp(s1) = specxp(s1) + w_esc ; CCDxp(m,n,s2) = CCDxp(m,n,s2) + w_esc ; elseif &
     (DirVec(1).eq.-1) then; specxm(s1) = specxm(s1) + w_esc ; CCDxm(m,n,s2) = CCDxm(m,n,s2) + w_esc ; elseif &
     (DirVec(2).eq. 1) then; specyp(s1) = specyp(s1) + w_esc ; CCDyp(m,n,s2) = CCDyp(m,n,s2) + w_esc ; elseif &
     (DirVec(2).eq.-1) then; specym(s1) = specym(s1) + w_esc ; CCDym(m,n,s2) = CCDym(m,n,s2) + w_esc ; elseif &
     (DirVec(3).eq. 1) then; speczp(s1) = speczp(s1) + w_esc ; CCDzp(m,n,s2) = CCDzp(m,n,s2) + w_esc ; elseif &
     (DirVec(3).eq.-1) then; speczm(s1) = speczm(s1) + w_esc ; CCDzm(m,n,s2) = CCDzm(m,n,s2) + w_esc
  else
    stop 'Bad DirVec'
  endif
endif

end subroutine Propagate

!------------------------------------------------------------------------------

subroutine Quicksort(n,arr)

use kinds

implicit none
integer,intent(in)::   n
real(dp),intent(inout):: arr(n)
integer,parameter::    m=7, nstack=50
integer::              i,ir,j,jstack,k,l,istack(nstack)
real(dp)::               a,temp
logical::              flag

!This subroutine is a f90-ed version of the f77 subroutine given in Numerical
!Recipes in Fortran.
!It sorts an array arr(1:n) into ascending numerical order using the Quicksort
!algorithm. n is input; arr is replaced on output by its sorted rearrangement.
!Parameters: m is the size of subarrays sorted by straight insertion and nstack
!is the required auxiliary storage.

jstack = 0
l      = 1
ir     = n

do
  if (ir-l .lt. m) then     !Insertion sort when subarray small enough.
    do j = l+1,ir
      a = arr(j)

      flag = .true.
      do i = j-1,l,-1
        if (arr(i) .le. a) then
          flag = .false.
          exit
        endif
        arr(i+1) = arr(i)
      enddo

      if (flag) i = l-1      !This flag shouldn't be necessary, if "i" always
      arr(i+1) = a           !has its final value when loop is exited. Does it?
    enddo

    if (jstack .eq. 0) exit

    ir     = istack(jstack)  !Pop stack and begin a new round of partitioning.
    l      = istack(jstack-1)
    jstack = jstack - 2
  else
    k        = (l+ir)/2      !Choose median of left, center, and right elements
                             !as partitioning element a. Also rearrange so that
                             !a(l) .le. a(l+1) .le. a(ir).
    temp     = arr(k)
    arr(k)   = arr(l+1)
    arr(l+1) = temp

    if (arr(l) .gt. arr(ir)) then
      temp    = arr(l)
      arr(l)  = arr(ir)
      arr(ir) = temp
    endif

    if (arr(l+1) .gt. arr(ir)) then
      temp     = arr(l+1)
      arr(l+1) = arr(ir)
      arr(ir)  = temp
    endif

    if (arr(l) .gt. arr(l+1)) then
      temp     = arr(l)
      arr(l)   = arr(l+1)
      arr(l+1) = temp
    endif

    i = l + 1                !Initialize pointers for partitioning.
    j = ir
    a = arr(l+1)             !Partitioning element.

    do                       !Beginning of innermost loop.
      i = i+1                !Scan up to find element > a.
      if (arr(i) .lt. a) cycle

      do
        j = j - 1            !Scan down to find element < a.
        if (arr(j) .le. a) exit
      enddo

      if (j .lt. i) exit     !Pointers crossed. Exit with partitioning complete.

      temp   = arr(i)        !Exchange elements.
      arr(i) = arr(j)
      arr(j) = temp
    enddo                    !End of innermost loop.

    arr(l+1) = arr(j)        !Insert partitioning element.
    arr(j)   = a
    jstack   = jstack + 2
                             !Push pointers to larger subarray on stack,
                             !process smaller subarray immediately.
    if (jstack .gt. nstack) stop 'nstack too small in sort'

    if (ir-i+1 .ge. j-l) then
      istack(jstack)   = ir
      istack(jstack-1) = i
      ir               = j - 1
    else
      istack(jstack)   = j - 1
      istack(jstack-1) = l
      l                = i
    endif
  endif
enddo

end subroutine Quicksort

!------------------------------------------------------------------------------

subroutine RandomUparal(x)

!Randomly selects the velocity of the scattering atom parallel to the
!direction of the incident photon, weighted with the convolution of a
!Maxwellian and a Lorentzian probability density

use kinds
use PhysicalConstants
use AtomParameters
use GlobalData
use DerivedParameters
use iMonteCarlo

implicit none
real(dp)::  x
real(dp)::  x_orig,u_0,p,theta_0,theta
real(dp)::  R1,R2,R3
real(dp)::  e_u02

x_orig = x                                      !Save original x
x      = abs(x)                                 !Use a positive x

if (x .ge. 0. .and. x .lt. 0.2) then            !
   u_0 = 0.                                     !
elseif (x .ge. 0.2 .and. x .lt. x_cw + 2.3) then!Determine transition
   u_0 = x - .01*a**.1666667*exp(1.2*x)         !velocity u_0 between the
else                                            !two comparison distributions
   u_0 = 4.5                                    !
endif                                           !

e_u02   = exp(-u_0**2)
theta_0 = atan((u_0-x) / a)

p = (theta_0 + pib2) / ((1 - e_u02)*theta_0 + (1 + e_u02)*pib2)

do
  call MonteCarlo(R1)

  if (R1 .le. p) then
    call MonteCarlo(R2)
    theta = -pib2 + R2*(theta_0 + pib2)
    u_II  = a*TAN(theta) + x
    call MonteCarlo(R3)
    if (R3 .le. exp(-u_II**2)) exit
  else
    call MonteCarlo(R2)
    theta = theta_0 + R2*(pib2-theta_0)
    u_II  = a*tan(theta) + x
    call MonteCarlo(R3)
    if (R3 .le. exp(-u_II**2)/e_u02) exit
  endif
enddo

u_II = sign(1d0,x_orig) * u_II                  !Revert to original sign
x    = x_orig

end subroutine RandomUparal

!------------------------------------------------------------------------------

subroutine RandomUperp(x,absorb)

!Randomly selects the velocities of the scattering atom perpendicular to the
!direction of the incident photon, weighted with a Gaussian probability
!density

use kinds
use PhysicalConstants
use GlobalData
use AtomParameters
use CurrentCellParameters
use iMonteCarlo

implicit none
real(dp),intent(in)::   x
logical,intent(out):: absorb
real(dp)::              R2(2)

if (abs(x) .ge. x_crit) then
  call MonteCarlo(R2)
  u_1 = sqrt(-log(R2(1))) * cos(tpi*R2(2))
  u_2 = sqrt(-log(R2(1))) * sin(tpi*R2(2))
else
# ifdef dust
  call RecoverCoreHistory(absorb)
# else
  absorb = .false.
# endif
  if (.not. absorb) then
    call MonteCarlo(R2)
    u_1 = sqrt(x_crit**2 - log(R2(1))) * cos(tpi*R2(2))
    u_2 = sqrt(x_crit**2 - log(R2(1))) * sin(tpi*R2(2))
  endif
endif

! BU: if (abs(x) .ge. x_crit) then
! BU:   call MonteCarlo(R2)
! BU:   u_1 = sqrt(-log(R2(1))) * cos(tpi*R2(2))
! BU:   u_2 = sqrt(-log(R2(1))) * sin(tpi*R2(2))
! BU: else
! BU: # ifdef dust
! BU:   call RecoverCoreHistory(absorb)
! BU:   if (.not. absorb) then
! BU: # endif
! BU:     call MonteCarlo(R2)
! BU:     u_1 = sqrt(x_crit**2 - log(R2(1))) * cos(tpi*R2(2))
! BU:     u_2 = sqrt(x_crit**2 - log(R2(1))) * sin(tpi*R2(2))
! BU: # ifdef dust
! BU:   endif
! BU: # endif
! BU: endif

end subroutine RandomUperp

!------------------------------------------------------------------------------

subroutine ReadData

!Read data from file into the following variables and arrays:
!-------------------------------
!Rec# | Variable    | Explanation
!-------------------------------
!  1  | N_cells     | Total number of cells (=ni*nj*nk for non-AMR).
!     | D_box       | Side length of computational box in kpc.
!     | ni,nj,nk    | Resolution of base gris in x-, y-, and z-direction.
!     |             | Note that at the moment MoCaLaTA is not fully prepared
!     |             | for these parameters not having the same value.
!     | n_ph        | Array containing the number of photon packets of the
!     |             | various types of photons. At the moment, the number of
!     |             | types is two --- Lya and FUV photons.
!     | emdimL      | Size of the second dimension of the InitLya array;
!     |             | emdimL = 3 if all photons are weighted equally, otherwise
!     |             | emdimL = 4.
!     | emdimF      | Same for InitFUV.
!  2  | InitLya     | Rank 2 array for the initial positions of every Lya
!     |             | photon. Each photon has three positions, plus possibly
!     |             | an associated weight.
!  3  | InitFUV     | Same for FUV photons.
!  4  | n_HIString  | Density of neutral hydrogen in cm-3.
!  5  | Dnu_DString | Gas temperature in Kelvin.
!  6  | U_xString   | Gas bulk velocity in x-direction in km/s.
!  7  | U_yString   |        -"-           y       -"-
!  8  | U_zString   |        -"-           z       -"-
! (9) | LevelString | Cell refinement levels (optional; in case of AMR).
!(10) | n_dString   | Gas metallicity in terms of Solar metallicity (optional;
!     |             | in case of including dust).
!(11) | n_HIIString | Density of ionized hydrogen in cm-3 (optional; in case of
!     |             | including dust in ionized gas).

use kinds
use PhysicalConstants
use ObservationalParameters
use GlobalData
use DataArrays
use iDoppler
use iv2u

implicit none
character(len=300):: infile
character(len=2)::   cn
integer::            i,l,n,n_str,n_recL,n_recF,ifracL,ifracF
real(dp)::           BWtemp(2)

infile = trim(DataDir) // trim(subdir) // trim(Celldata)
call getlun(l)
open(l,file=trim(infile),form='unformatted',status='old',action='read')

read(l) N_cells,D_box,ni,nj,nk,L_tot,BWtemp,n_ph,emdimL,emdimF,n_recL,n_recF!Remember: n_ph = n_ph(2)
if (any((/emdimL,emdimF/) .ne. 3)) stop 'Not ready for photon weighting scheme. Sorry.'
if (n_phtot .ne. sum(n_ph)) then
  write(*,*) "Error: Keyword 'n_phtot' has been changed in '", trim(parfile), "' since"
  write(*,*) "it was used with 'InitialEmission', from", sum(n_ph), 'to', n_phtot, '.'
  write(*,*) "Don't do that, okay? You're messing things up."
  stop
endif
if (any(abs(BWtemp-BW) .gt. 1e-6)) then
  write(*,*) "Error: Keyword 'BW' has been changed in '", trim(parfile), "' since"
  write(*,*) "it was used with 'InitialEmission', from ", real(BWtemp), 'to'
  write(*,*) real(BW)
  write(*,*) "Don't do that, okay? You're messing things up."
  stop
endif

D_box   = D_box * kpc
R_box   = D_box / 2d0
BaseRes = (/ni,nj,nk/)
n_phtot = sum(n_ph)
Lya2tot = real(n_ph(1)) / n_phtot               !Probability for an emitted photon to be a continuum photon

n_str = 7                                       !Number of "strings"/records
allocate(InitLya(n_ph(1),emdimL))
allocate(InitFUV(n_ph(2),emdimF))
allocate(n_HIString(N_cells))
allocate(Dnu_DString(N_cells))
allocate(U_xString(N_cells))
allocate(U_yString(N_cells))
allocate(U_zString(N_cells))
# ifdef AMR
  allocate(LevelString(N_cells))
  n_str = n_str + 1
# endif
# ifdef dust
  allocate(n_dString(N_cells))
  allocate(n_HIIString(N_cells))
  n_str = n_str + 2
# endif

call Int2Char(2*n_str+15,cn,.false.)
write(*,'(a'//cn//')',advance='no') &
  'Loading data' // repeat(' ',n_str) // '|' // repeat(achar(8),n_str+1)

ifracL = n_ph(1) / n_recL
ifracF = n_ph(2) / n_recF

do i = 1,n_recL
  if (i .lt. n_recL) then
    read(l) (InitLya(n,1:emdimL), n = (i-1)*ifracL+1, i*ifracL)
  else
    read(l) (InitLya(n,1:emdimL), n = (i-1)*ifracL+1, n_ph(1))
  endif
enddo
write(*,'(a)',advance='no') '.'

do i = 1,n_recF
  if (i .lt. n_recF) then
    read(l) (InitFUV(n,1:emdimF), n = (i-1)*ifracF+1, i*ifracF)
  else
    read(l) (InitFUV(n,1:emdimF), n = (i-1)*ifracF+1, n_ph(2))
  endif
enddo
write(*,'(a)',advance='no') '.'

# ifdef AMR
read(l) LevelString; write(*,'(a)',advance='no') '.'
# endif
read(l) n_HIString;  write(*,'(a)',advance='no') '.'
read(l) Dnu_DString; write(*,'(a)',advance='no') '.'
read(l) U_xString;   write(*,'(a)',advance='no') '.'
read(l) U_yString;   write(*,'(a)',advance='no') '.'
read(l) U_zString;   write(*,'(a)',advance='no') '.'
# ifdef dust
read(l) n_dString;   write(*,'(a)',advance='no') '.'
read(l) n_HIIString; write(*,'(a)',advance='no') '.'
# endif
write(*,'(a)',advance='yes') '.'
close(l)

# ifdef AMR
  L_max = maxval(LevelString)                   !Maximum level of refinement
# else
  L_max = 0
# endif

where(n_HIString .lt. 1e-37) n_HIString = 1e-37 !Make sure there are no zeros

!           ------------   Convert parameters:   -----------

!Express positions in cm rather than kpc:
InitLya(:,1:3) = InitLya(:,1:3) * kpc
InitFUV(:,1:3) = InitFUV(:,1:3) * kpc

!Convert temperature in K to "width of the associated Doppler broadened Lya
!line in frequency", i.e. Dnu_DString now contains Dnu_D:
Dnu_DString = Doppler(Dnu_DString)

!Convert velocity to "velocity in terms of 'width of the associated Doppler
!broadened Lya line in velocity'", i.e. U_[xyz]String now contains u_[xyz] =
!v_[xyz] / v_th:
U_xString = v2u(U_xString,Dnu_DString)
U_yString = v2u(U_yString,Dnu_DString)
U_zString = v2u(U_zString,Dnu_DString)

!Convert metallicity to "dust density", i.e. a rescaled hydrogen density such
!that multiplying this number by the dust cross section per  hydrogen nucleus
!(in the subroutine "XsecD") gives the optical depth of dust per unit distance:
# ifdef dust
n_dString = (n_HIString + f_ion*n_HIIString) * n_dString/Z0
deallocate(n_HIIString)
# endif

end subroutine ReadData

!------------------------------------------------------------------------------

subroutine ReadMultiphase

use kinds
use PhysicalConstants
use GlobalData
use DataArrays
use ObservationalParameters
use iMonteCarlo
use iDoppler
use igasdev

implicit none
integer::            i,j,k,l,n,cl
real(dp)::           R,R3(3),X_cl(3),d_cl,X(3),nhat(3),dd
real(dp)::           tt1,tt2
character(len=300):: infile

Dnu_DString   = Doppler((/T_ICM, T_cl, T_ICM/)) !Doppler width in couds; Hz
n_HIString(0) = 1e-37                           !Emptiness outside shell
n_HIString(1) = max(n_HI_cl,1d-37)              !Neutral hydrogen density in clouds
n_HIString(2) = max(n_HI_ICM,1d-37)             !Neutral hydrogen density in ICM
# ifdef dust
n_dString(0)  = 1e-37                           !Emptiness outside shell
n_dString(1)  = max(n_HI_cl  * Z_cl  / Z0,1d-37)!"Dust density" in clouds
n_dString(2)  = max(n_HI_ICM * Z_ICM / Z0,1d-37)!"Dust density" in ICM
# endif

infile = trim(DataDir) // trim(subdir) // trim(Celldata)
call getlun(l)
open(l,file=trim(infile),form='unformatted',status='old',action='read')

read(l) N_cells

allocate(Clouds(N_cl,4))
allocate(CloudVelDisp(N_cl,3))
# ifdef AMR
allocate(LevelString(N_cells))
read(l) LevelString
L_max = maxval(LevelString)                     !Maximum level of refinement
# endif
# ifdef multi
allocate(PhaseString(N_cells))
allocate(CloudString(N_cells))
read(l) PhaseString
read(l) CloudString
# endif
read(l) Clouds                                  !Cloud x,y,z,r in kpc

close(l)

Clouds = Clouds * kpc

do n = 1,N_cl
  call gasdev(R3)
  CloudVelDisp(n,:) = R3 * sigV_cl              !Assign x,y,z vel disp to clouds
enddo

end subroutine ReadMultiphase

!------------------------------------------------------------------------------

subroutine ReadInput

use kinds
use PhysicalConstants
use GlobalData
use ObservationalParameters
use CurrentCellParameters
use DerivedParameters
use iDoppler
use DataArrays

implicit none
character(len=200)::          dummy,n_init
real(dp)::                    sigma,sigma_0,CoreSkip
real(dp)::                    NN_HI,NN_HI_sh,norm,n0(3),SFR,LumDist,zero,temp
logical::                     alright,TrueRan
integer::                     narg,l,n_slab,c1,c2
character(len=12),parameter:: numsi = '-.0123456789'!Possible characters in a number
character(len=2), parameter:: delim = ' ,'          !Delimiter between initial direction vector indicees
character(len=200)::          command
integer::                     seedsize
integer,allocatable::         seed(:),time(:)

x_critType = ''
x_injType  = ''
FF         = 1d0
at         = 0d0
tau_a      = -1
Z_cl       = -1
EBV        = -1
X_QSOc     = 'nil'
sigV_cl    = 0.
writeN0    = .false.
nullifyx   = .false.

call get_command_argument(1,parfile)
narg = command_argument_count()
if (narg .eq. 0) then
  call get_command(command)
  write(*,*)
  write(*,*) 'The proper syntax is:'
  write(*,*) '   ' // trim(command) // ' params.in'
  write(*,*) "where 'params.in' is the name of your input parameter file. Note that there is"
  write(*,*) "is no '<' before the argument. This is the new, sweet 2003 standard, allowing"
  write(*,*) "user-given input from standard input at a later point in the execution."
  write(*,*)
  stop
endif
call getlun(l)
open (l,file=trim(parfile),form='formatted',status='old',action='read')

read(l,*) dummy
read(l,*) model

  select case (trim(model))
    case ('data')
    ! I/O
      read(l,*) dummy
      read(l,*) DataDir      !Mother directory
      read(l,*) Subdir       !Subdirectory
      read(l,*) dummy        !Name of data file
      read(l,*) CellData     !Name of data file
      read(l,*) wrtph        !Write to std. output for every wrtph phot
      read(l,*) wrtfl        !Write file for every wrtfl photon
      read(l,*) outlabel     !Label for main output files
      read(l,*) specin       !Write input spectrum character
      read(l,*) specout      !-"- for output spectrum in Doppler widths

    ! Observations
      read(l,*) dummy
      read(l,*) BW           !Wavelength interval bounds in Angstrom
      read(l,*) pix          !# of pixels/side in CCD
      read(l,*) SpecRes1D    !Resolution of 1D spectrum; bins
      read(l,*) SpecRes2D    !Resolution of 2D spectrum; bins
      read(l,*) WhichCCDs    !(xp,xm,yp,ym,zp,zm)

    ! Simulation
      read(l,*) dummy
      read(l,*) D_obs        !Side length of area covered by CCD (centered on box center)
      read(l,*) dummy        !Side length of box, outside which lum. is set to zero (centered on box center)
      read(l,*) n_phtot      !Total # of photons; Lya + cont.
      read(l,*) dummy
      read(l,*) dummy
      read(l,*) x_critType   !'intrinsic', 'global', or user-given value (as string) with prefix 'm' for maximum
      read(l,*) x_injType    !Initial frequency: 'proper', '...
      read(l,*) X_init       !Initial position: 'central', 'homo', 'lumdep', or '(x y z)'
      read(l,*) n_init       !Initial direction: 'iso' or '(nx ny nz)'
      read(l,*) TrueRan      !Use "true" random numbers, given by current date and time
      read(l,*) N_los        !# of sightlines for calculating average quantities

    ! Gas
      read(l,*) dummy
      read(l,*) recoil       !Include atom recoil

    ! Dust
      read(l,*) dummy
      read(l,*) DustType     !'SMC', 'LMC' or user-given cross section
      read(l,*) f_ion        !Fraction of ionized hydrogen that contributes to dust density
      read(l,*) albedo       !Dust albedo
      read(l,*) g            !Dust scattering asymmetry parameter

    ! Cosmology
      read(l,*) dummy
      read(l,*) z            !Redshift of snapshot
      read(l,*) H_0          !Hubble constant, km/s/Mpc
      read(l,*) Omega_M      !Matter density parameter
      read(l,*) Omega_L      !Dark energy density parameter

    case ('slab')
    ! I/O
      read(l,*) dummy
      read(l,*) DataDir      !Mother directory
      read(l,*) Subdir       !Subdirectory
      read(l,*) wrtph        !Write to std. output for every wrtph phot
      read(l,*) wrtfl        !Write file for every wrtfl photon
      read(l,*) outlabel     !Label for main output files
      read(l,*) specin       !-"- for input spectrum in Doppler widths
      read(l,*) specout      !-"- for output spectrum in Doppler widths

    ! Observations
      read(l,*) dummy
      read(l,*) BW           !Wavelength interval bounds in Angstrom
      read(l,*) pix          !# of pixels/side in CCD
      read(l,*) SpecRes1D    !Resolution of 1D spectrum; bins
      read(l,*) SpecRes2D    !Resolution of 2D spectrum; bins
      read(l,*) WhichCCDs    !(xp,xm,yp,ym,zp,zm)

    ! Simulation
      read(l,*) dummy
      read(l,*) ni           !Base grid resolution, x direction
      read(l,*) nj           !                      y direction
      read(l,*) nk           !                      z direction
      read(l,*) n_slab       !Total slab thickness in # of cells
      read(l,*) D_box        !Box size in kpc
      read(l,*) D_obs        !Side length of area covered by CCD (centered on box center)
      read(l,*) n_phtot      !Total # of emitted photon packets; Lya + cont.
      read(l,*) SFR          !Star formation rate in Msun/yr (sets normalization)
      read(l,*) EW_int       !Intrinsic EW in Angstrom
      read(l,*) x_critType   !'intrinsic', 'global', or user-given value (as string) with prefix 'm' for maximum
      read(l,*) x_injType    !Initial frequency: 'proper', '...
      read(l,*) X_init       !Initial position: 'central', 'homo', 'lumdep', or '(x y z)'
      read(l,*) n_init       !Initial direction: 'iso' or '(nx ny nz)'
      read(l,*) TrueRan      !Use "true" random numbers, given by current date and time
      read(l,*) N_los        !# of sightlines for calculating average quantities

    ! Gas                 !Give only compatible values. To let system decide, give a value -1.
      read(l,*) dummy
      read(l,*) at           !a * tau_0 from center to face of slab
      read(l,*) T_cl         !Gas temperature
      read(l,*) n_HI_cl      !Neutral hydrogen density
      read(l,*) NN_HI        !Neutral hydrogen column density from center to face of slab
      read(l,*) V_z          !Slab velocity in z direction in km/s
      read(l,*) recoil       !Include recoil

    ! Dust
      read(l,*) dummy
      read(l,*) tau_a        !Dust absorption optical depth from center to face of slab
      read(l,*) Z_cl         !Metallicity in terms of Solar
      read(l,*) EBV          !Color excess from center to face of slab
      read(l,*) DustType     !'SMC', 'LMC', or user-given cross section
      read(l,*) albedo       !Dust albedo
      read(l,*) g            !Dust scattering asymmetry parameter

    ! Cosmology
      read(l,*) dummy
      read(l,*) z            !Redshift of snapshot
      read(l,*) H_0          !Hubble constant, km/s/Mpc
      read(l,*) Omega_M      !Matter density parameter
      read(l,*) Omega_L      !Dark energy density parameter

      n_phases = 1

      z_0 = .5d0 * n_slab/real(nk) * D_box*kpc

    case ('shell')
    ! I/O
      read(l,*) dummy
      read(l,*) DataDir      !Mother directory
      read(l,*) Subdir       !Subdirectory
      read(l,*) wrtph        !Write to std. output for every wrtph phot
      read(l,*) wrtfl        !Write file for every wrtfl photon
      read(l,*) outlabel     !Label for main output files
      read(l,*) specin          !-"- for input spectrum in Doppler widths
      read(l,*) specout         !-"- for output spectrum in Doppler widths

    ! Observations
      read(l,*) dummy
      read(l,*) BW           !Wavelength interval bounds in Angstrom
      read(l,*) pix          !# of pixels/side in CCD
      read(l,*) SpecRes1D    !Resolution of 1D spectrum; bins
      read(l,*) SpecRes2D    !Resolution of 2D spectrum; bins
      read(l,*) WhichCCDs    !(xp,xm,yp,ym,zp,zm)

    ! Simulation
      read(l,*) dummy
      read(l,*) ni           !Base grid resolution, x direction
      read(l,*) nj           !                      y direction
      read(l,*) nk           !                      z direction
      read(l,*) D_box        !Box size in kpc
      read(l,*) D_obs        !Side length of area covered by CCD (centered on box center)
      read(l,*) r_inner      !Radius of inner shell surface in kpc
      read(l,*) r_outer      !Radius of outer shell surface in kpc
      read(l,*) n_phtot      !Total # of emitted photon packets; Lya + cont.
      read(l,*) SFR          !Star formation rate in Msun/yr (sets normalization)
      read(l,*) EW_int       !Intrinsic EW in Angstrom
      read(l,*) x_critType   !'intrinsic', 'global', or user-given value (as string) with prefix 'm' for maximum
      read(l,*) x_injType    !Initial frequency: 'proper', '...
      read(l,*) X_init       !Initial position: 'central', 'homo', 'lumdep', or '(x y z)'
      read(l,*) n_init       !Initial direction: 'iso' or '(nx ny nz)'
      read(l,*) TrueRan      !Use "true" random numbers, given by current date and time
      read(l,*) N_los        !# of sightlines for calculating average quantities

    ! Gas                     !Give only compatible values. To let system decide, give a value -1.
      read(l,*) dummy
      read(l,*) at           !a * tau_0 from center to outer face of shell
      read(l,*) T_cl         !Gas temperature
      read(l,*) n_HI_cl      !Neutral hydrogen density
      read(l,*) sigV_cl      !
      read(l,*) NN_HI        !Neutral hydrogen column density from center to outer face of shell
      read(l,*) V_out        !Outflow velocity in km/s (<0 for collapse)
      read(l,*) Vprof        !Velocity profile
      read(l,*) recoil       !Include atom recoil

    ! Dust
      read(l,*) dummy
      read(l,*) tau_a        !Cloud metallicity in terms of Solar
      read(l,*) Z_cl         !Metallicity in terms of Solar
      read(l,*) EBV          !Color excess
      read(l,*) DustType     !'SMC', 'LMC' or user-given cross section
      read(l,*) albedo       !Dust albedo
      read(l,*) g            !Dust scattering asymmetry parameter

    ! Cosmology
      read(l,*) dummy
      read(l,*) z            !Redshift of snapshot
      read(l,*) H_0          !Hubble constant, km/s/Mpc
      read(l,*) Omega_M      !Matter density parameter
      read(l,*) Omega_L      !Dark energy density parameter

      n_phases = 1

      r_inner = r_inner * kpc
      r_outer = r_outer * kpc
      z_0     = r_outer - r_inner

    case ('multiphase')
    ! I/O
      read(l,*) dummy
      read(l,*) DataDir    !Mother directory
      read(l,*) Subdir     !Subdirectory with with input file
      read(l,*) dummy      !Particle data
      read(l,*) dummy      !Particle data, formatted for visualization; '' to omit
      read(l,*) CellData   !Cell data; levels, phases, and associated clouds
      read(l,*) dummy      !Cell data, formatted for visualization; '' to omit
      read(l,*) wrtph      !Write to std. output for every wrtph photon
      read(l,*) wrtfl      !Write file for every wrtfl photon
      read(l,*) outlabel   !Label for main output files
      read(l,*) specin     !-"- for input spectrum in Doppler widths; '' to omit
      read(l,*) specout    !-"- for output spectrum in Doppler widths; '' to omit

    ! Observations
      read(l,*) dummy
      read(l,*) BW         !Wavelength interval bounds in Angstrom
      read(l,*) pix        !# of pixels/side in CCD
      read(l,*) SpecRes1D  !Resolution of 1D spectrum; bins
      read(l,*) SpecRes2D  !Resolution of 2D spectrum; bins
      read(l,*) WhichCCDs  !(xp,xm,yp,ym,zp,zm)

    ! Simulation
      read(l,*) dummy
      read(l,*) ni         !Base grid resolution, x direction
      read(l,*) nj         !                      y direction
      read(l,*) nk         !                      z direction
      read(l,*) D_box      !Box size in kpc
      read(l,*) D_obs      !Side length of area covered by CCD (centered on box center)
      read(l,*) n_phtot    !Total # of photons; Lya + cont.
      read(l,*) SFR        !Star formation rate in Msun/yr (sets normalization)
      read(l,*) EW_int     !Intrinsic EW in Angstrom
      read(l,*) x_critType !'intrinsic', 'global', '<value>', or 'max<value>', where <value> is a real number
      read(l,*) x_injType  !'proper', 'fwhm<value/kms>', or 'x_inj<value>'
      read(l,*) X_init     !Initial position: 'lumdep', 'central', 'homo', or '(x y z)', with x,y,z given in terms of D_box
      read(l,*) n_init     !Initial direction: 'iso' or '(nx ny nz)'
      read(l,*) CloudCorr  !Emission cloud correlation factor ([0,1], where 0 is "no emission from clouds" and 1 is "no emission from ICM")
      read(l,*) TrueRan    !Use "true" random number (given by current date and time)?
      read(l,*) N_los        !# of sightlines for calculating average quantities
      read(l,*) writeN0    !Write N0 (# of cloud scattering)?
      read(l,*) dummy      !Dots per clouds in particle representation

    ! Galaxy
      read(l,*) dummy
      read(l,*) r_gal      !Maximum radius of galaxy
      read(l,*) axes       !Ellipsoid radii relations (will be normalized to maxval)
      read(l,*) H1         !Disk scale length in plane
      read(l,*) exp1       !Density gradient exponent in plane
      read(l,*) H2         !Disk scale height perpendicular to plane
      read(l,*) exp2       !Density gradient exponent perp. to plane
      read(l,*) N_cl       !# of clouds
      read(l,*) rmin       !Minimum cloud radius in kpc
      read(l,*) rmax       !Maximum cloud radius in kpc
      read(l,*) beta       !Cloud radius distribution power law slope
      read(l,*) followGrad !Should cloud density follow galactic density gradient?
      read(l,*) T_cl       !Cloud gas temperature
      read(l,*) n_HI_cl    !Cloud neutral hydrogen density
      read(l,*) Z_cl       !Cloud metallicity in terms of Solar
      read(l,*) sigV_cl    !Velocity dispersion standard deviation of clouds in km/s
      read(l,*) T_ICM      !ICM gas temperature
      read(l,*) n_HI_ICM   !ICM neutral hydrogen density in the center
      read(l,*) Z_ICM      !ICM metallicity in terms of Solar
      read(l,*) V_in       !Cloud and ICM infall velocity in km/s (>0 for outflow)

    ! Wind
      read(l,*) dummy      
      read(l,*) r_inner    !Radius of inner shell in kpc (or cells if preceded by a minus)
      read(l,*) r_outer    !Radius of outer shell     -- " --
      read(l,*) NN_HI_sh   !HI column density of shell
      read(l,*) Z_sh       !Shell metallicity in terms of Solar
      read(l,*) T_sh       !Shell gas temperature
      read(l,*) V_out      !Outflow velocity in km/s (>0 for collapse)
      read(l,*) Vprof        !Velocity profile
      read(l,*) OA         !(Full) opening angle in degrees
      read(l,*) n_jet      !Direction of one of the jets (the other is opposite)

    ! Gas and dust    
      read(l,*) dummy      
      read(l,*) recoil     !Include atom recoil?
      read(l,*) DustType   !'SMC', 'LMC', or fixed cross section per hydrogen nucleus, e.g. '1.03e-21'
      read(l,*) albedo     !Dust albedo
      read(l,*) g          !Dust scattering asymmetry parameter

    ! Background QSO
      read(l,*) dummy
      read(l,*) X_QSOc     !QSO position in kpc (box center is [0,0,0]) => calulate N_HI from X_QSO through box

    ! Cosmology
      read(l,*) dummy
      read(l,*) z          !Redshift of snapshot
      read(l,*) H_0        !Hubble constant, km/s/Mpc
      read(l,*) Omega_M    !Matter density parameter
      read(l,*) Omega_L    !Dark energy density parameter

      H1       = H1      * kpc
      H2       = H2      * kpc
      r_inner  = r_inner * kpc
      r_outer  = r_outer * kpc
      r_gal    = r_gal   * kpc                  !Maximum galaxy radius of ICM and clouds
      axes     = axes    / maxval(axes)         !Make one axis equal unity
      n_phases = 2


    case ('semireal')
    ! I/O
      read(l,*) dummy
      read(l,*) DataDir      !Mother directory
      read(l,*) Subdir       !Subdirectory
      read(l,*) wrtph        !Write to std. output for every wrtph phot
      read(l,*) wrtfl        !Write file for every wrtfl photon
      read(l,*) outlabel     !Label for main output files
      read(l,*) specin       !Write input spectrum character
      read(l,*) specout      !-"- for output spectrum in Doppler widths
      read(l,*) Cloudfile    !File for writing cloud cells; empty string to not write
      read(l,*) Windfile     !File for writing wind cells; empty string to not write

    ! Observations
      read(l,*) dummy
      read(l,*) BW            !Wavelength interval bounds in Angstrom
      read(l,*) pix           !# of pixels/side in CCD
      read(l,*) SpecRes1D     !Resolution of 1D spectrum; bins
      read(l,*) SpecRes2D     !Resolution of 2D spectrum; bins
      read(l,*) WhichCCDs     !(xp,xm,yp,ym,zp,zm)

    ! Simulation
      read(l,*) dummy
      read(l,*) ni           !Base grid resolution, x direction
      read(l,*) nj           !                      y direction
      read(l,*) nk           !                      z direction
      read(l,*) D_box        !Box size in kpc
      read(l,*) D_obs        !Side length of area covered by CCD (centered on box center)
      read(l,*) n_phtot      !Total # of emitted photon packets; Lya + cont.
      read(l,*) SFR          !Star formation rate in Msun/yr (sets normalization)
      read(l,*) EW_int       !Intrinsic EW in Angstrom
      read(l,*) x_critType   !'intrinsic', 'global', or user-given value (as string) with prefix 'm' for maximum
      read(l,*) x_injType    !Initial frequency: 'proper', '...
      read(l,*) X_init       !Initial position: 'central', 'homo', 'lumdep', or '(x y z)'
      read(l,*) n_init       !Initial direction: 'iso' or '(nx ny nz)'
      read(l,*) TrueRan      !Use "true" random numbers, given by current date and time
      read(l,*) N_los        !# of sightlines for calculating average quantities

    ! Galaxy
      read(l,*) dummy
      read(l,*) r_gal
      read(l,*) H1
      read(l,*) exp1
      read(l,*) H2
      read(l,*) exp2
      read(l,*) N_cl       !# of clouds
      read(l,*) rf         !Cloud radius in terms of box radius
      read(l,*) followGrad
      read(l,*) T_cl       !Cloud gas temperature
      read(l,*) n_HI_cl    !Cloud neutral hydrogen density
      read(l,*) Z_cl       !Cloud metallicity in terms of Solar
      read(l,*) sigV_cl    !
      read(l,*) T_ICM      !ICM gas temperature
      read(l,*) n_HI_ICM   !ICM neutral hydrogen density
      read(l,*) Z_ICM      !ICM metallicity in terms of Solar
      read(l,*) V_in       !Infall velocity for clouds and ICM

    ! Wind
      read(l,*) dummy
      read(l,*) r_inner    !Radius of inner shell surface in kpc
      read(l,*) r_outer    !Radius of outer shell surface in kpc
      read(l,*) NN_HI_sh   !HI Column density of shell
      read(l,*) Z_sh       !Shell metallicity in terms of Solar
      read(l,*) T_sh       !Shell gas temperature
      read(l,*) V_out      !Outflow velocity in km/s (>0 for collapse)
      read(l,*) Vprof        !Velocity profile
      read(l,*) OA         !(Full) opening angle in degrees
      read(l,*) n_jet      !Direction of one of the jets (not necessarily unit vector)

    ! Gas and dust
      read(l,*) dummy
      read(l,*) recoil     !Include atom recoil
      read(l,*) DustType   !'SMC', 'LMC' or user-given cross section
      read(l,*) albedo     !Dust albedo
      read(l,*) g          !Dust scattering asymmetry parameter

    ! Background QSO
      read(l,*) dummy
      read(l,*) X_QSOc     !QSO position in kpc (box center is [0,0,0]) => calulate N_HI from X_QSO through box

    ! Cosmology
      read(l,*) dummy
      read(l,*) z          !Redshift of snapshot
      read(l,*) H_0        !Hubble constant, km/s/Mpc
      read(l,*) Omega_M    !Matter density parameter
      read(l,*) Omega_L    !Dark energy density parameter

      n_phases = 3

      if (r_inner .lt. 0) then
        r_inner = (-r_inner -1) * D_box/ni
        r_outer = (-r_outer)    * D_box/ni
      endif
      H1        = H1    * kpc
      H2        = H2    * kpc
      r_gal     = r_gal * kpc                   !Galaxy "radius" of ICM and clouds
      r_cl      = rf * D_box/2 * kpc            !Cloud radius in kpc
      r_inner   = r_inner * kpc
      r_outer   = r_outer * kpc
      z_0       = r_outer - r_inner             !Shell thickness in kpc
      n_jet     = n_jet / norm(n_jet)           !Make unit vector
      cosOAb2sq = cos(OA/2 / 180 * pi)**2       !cos^2(OA/2) in rad
      n_HI_sh   = NN_HI_sh / z_0                !HI density in shell
      write(*,*) 'n_HI in shell:', real(n_HI_sh)

    case ('florent')
    ! I/O
      read(l,*) dummy
      read(l,*) DataDir      !Mother directory
      read(l,*) Subdir       !Subdirectory
      read(l,*) wrtph        !Write to std. output for every wrtph phot
      read(l,*) wrtfl        !Write file for every wrtfl photon
      read(l,*) outlabel     !Label for main output files
      read(l,*) specin          !-"- for input spectrum in Doppler widths
      read(l,*) specout         !-"- for output spectrum in Doppler widths

    ! Observations
      read(l,*) dummy
      read(l,*) BW           !Wavelength interval bounds in Angstrom
      read(l,*) pix          !# of pixels/side in CCD
      read(l,*) SpecRes1D    !Resolution of 1D spectrum; bins
      read(l,*) SpecRes2D    !Resolution of 2D spectrum; bins
      read(l,*) WhichCCDs    !(xp,xm,yp,ym,zp,zm)

    ! Simulation
      read(l,*) dummy
      read(l,*) ni           !Base grid resolution, x direction
      read(l,*) nj           !                      y direction
      read(l,*) nk           !                      z direction
      read(l,*) D_box        !Box size in kpc
      read(l,*) D_obs        !Side length of area covered by CCD (centered on box center)
      read(l,*) r_inner      !Radius of inner shell surface in kpc
      read(l,*) r_outer      !Radius of outer shell surface in kpc
      read(l,*) n_phtot      !Total # of emitted photon packets; Lya + cont.
      read(l,*) SFR          !Star formation rate in Msun/yr (sets normalization)
      read(l,*) EW_int       !Intrinsic EW in Angstrom
      read(l,*) x_critType   !'intrinsic', 'global', or user-given value (as string) with prefix 'm' for maximum
      read(l,*) x_injType    !Initial frequency: 'proper', '...
      read(l,*) X_init       !Initial position: 'central', 'homo', 'lumdep', or '(x y z)'
      read(l,*) n_init       !Initial direction: 'iso' or '(nx ny nz)'
      read(l,*) TrueRan      !Use "true" random numbers, given by current date and time
      read(l,*) N_los        !# of sightlines for calculating average quantities

    ! Gas                 !Give only compatible values. To let system decide, give a value -1.
      read(l,*) dummy
      read(l,*) FF           !Cloud filling factor
      read(l,*) T_cl         !Cloud gas temperature
      read(l,*) n_HI_cl      !Cloud neutral hydrogen density
      read(l,*) sigV_cl      !
      read(l,*) T_ICM        !ICM gas temperature
      read(l,*) n_HI_ICM     !ICM neutral hydrogen density
      read(l,*) V_out        !Outflow velocity in km/s (<0 for collapse)
      read(l,*) Vprof        !Velocity profile
      read(l,*) recoil       !Include atom recoil

    ! Dust
      read(l,*) dummy
      read(l,*) tau_a        !Dust absorption optical depth from center to face
      read(l,*) EBV          !Color excess form center to face
      read(l,*) Z_cl         !Cloud metallicity in terms of Solar
      read(l,*) Z_ICM        !ICM metallicity in terms of Solar
      read(l,*) DustType     !'SMC', 'LMC' or user-given cross section
      read(l,*) albedo       !Dust albedo
      read(l,*) g            !Dust scattering asymmetry parameter

    ! Cosmology
      read(l,*) dummy
      read(l,*) z            !Redshift of snapshot
      read(l,*) H_0          !Hubble constant, km/s/Mpc
      read(l,*) Omega_M      !Matter density parameter
      read(l,*) Omega_L      !Dark energy density parameter

      n_phases = 2

      if (r_inner .lt. 0) then
        r_inner = (-r_inner -1) * D_box/ni
        r_outer = (-r_outer)    * D_box/ni
      endif
      r_inner = r_inner * kpc
      r_outer = r_outer * kpc
      z_0     = r_outer - r_inner
      NN_HI   = (FF*n_HI_cl + (1-FF)*n_HI_ICM) * z_0

    case default
      stop 'Argh unknown model in ReadInput!'
  end select

close(l)

! Convert to system units
H_0   = H_0 * 1d2                               !Convert from km/s/MpC to cm/s/kpc
d_L   = LumDist(z)                              !Luminosity distance in kpc
D_box = D_box * kpc
R_box = D_box / 2d0

! Check input consistency
if (trim(model).eq.'slab' .or. trim(model).eq.'shell') then

  alright = .false.
  if (at.gt.0) then
    if (any((/n_HI_cl,NN_HI/).gt.0)) stop "If 'at' is given, please set n_HI_cl and NN_HI to -1. T_cl is optional."
    if (T_cl .lt. 0) then
      T_cl = 10d0
      write(*,*) 'T_cl is set to 10 K.'
    endif
    Dnu_D_cl = Doppler(T_cl)
    sigma_0  = sigma(Dnu_D_cl,0d0)
    a        = Dnu_L / (2 * Dnu_D_cl)
    tau_0    = at / a                            !From center to face
    NN_HI    = tau_0 / sigma_0                   !From center to face
    n_HI_cl  = NN_HI / z_0                       !In all slab cells
    alright  = .true.
  else
    if (n_HI_cl.gt.0) then
      if (T_cl.lt.0) stop "If n_HI_cl is given, please give also T_cl"
      if (NN_HI .gt. 0) stop "Please use either n_HI_cl or NN_HI, not both"
      alright = .true.
    elseif (NN_HI.gt.0) then
      if (T_cl.lt.0) stop "If NN_HI is given, please give also T_cl"
      if (n_HI_cl .gt. 0) stop "Please use either n_HI_cl or NN_HI, not both"
      n_HI_cl    = NN_HI / z_0
      alright = .true.
    endif
    Dnu_D_cl   = Doppler(T_cl)
    sigma_0 = sigma(Dnu_D_cl,0d0)
    a       = Dnu_L / (2 * Dnu_D_cl)
    NN_HI   = n_HI_cl * z_0
    tau_0   = NN_HI * sigma_0
    at      = a * tau_0
  endif
  if (.not.alright) then
    write(*,*) "Please give either of the sets (at), (at,T_cl), (n_HI_cl,T_cl), or (NN_HI,T_cl),"
    write(*,*) "and set the remaining parameters to -1."
    stop
  endif

  ! Make sure BW is in Angstrom
  if (BW(1).lt.0 .or. BW(1).gt.lambda_0) then     !BW is probably given in x-units
    temp  = BW(1)
    BW(1) = BW(2)
    BW(2) = temp
    BW    = c / (BW*Dnu_D_cl + nu_0) * 1e8
  endif
endif

if (T_ICM .gt. 0) then
  if (abs(T_ICM-T_cl).gt.1e-6 .and. len_trim(specout).gt.0) then
    write(*,*) 'For the sampling of x to make sense, all temperatures must be equal, INCLUDING'
    write(*,*) "'void' temperature, which is set not in the input file, but in the subroutine"
    write(*,*) "'Make[***]', where [***] corresponds to the appropriate model."
    call holdon(.2)
    write(*,*) 'Continuing, but writing x = 0.'
    call holdon(.1)
    nullifyx = .true.
  endif
endif

! Determine dust cross section
if (trim(DustType) .eq. 'SMC') then
  sigma_d = .395083d-21
  Z0      = 10**(-.6)                           !Reference metallicity in terms of Solar
  ta2EBV  = 0.0640968
elseif (trim(DustType) .eq. 'LMC') then
  sigma_d = .723314d-21
  Z0      = 10**(-.3)
  ta2EBV  = 0.0576759
else
  read(DustType,*) userXsecD
  sigma_d = userXsecD
  Z0      = 1d0
  ta2EBV  = .11                                  !According to ver06
endif

zero = -1e-6
if (trim(model).eq.'slab'          .or. &
    trim(model).eq.'shell'         .or. &
    trim(model).eq.'florent') then
  if (tau_a.gt.zero .and. any((/Z_cl, EBV/) .gt.zero) .or. &
      Z_cl .gt.zero .and. any((/tau_a,EBV/) .gt.zero) .or. &
      EBV  .gt.zero .and. any((/tau_a,Z_cl/).gt.zero)) stop 'Please only provide one of tau_a, Z_cl, and EBV'

  if (EBV .gt. zero) then
    tau_a = EBV / ta2EBV
    Z_cl  = Z0 * tau_a / ((1-albedo) * NN_HI * sigma_d)
  endif
  if (tau_a .gt. zero) then
    Z_cl  = Z0 * tau_a / ((1-albedo) * NN_HI * sigma_d)
    EBV   = ta2EBV * tau_a
  endif
  if (Z_cl .gt. zero) then
    tau_a = (1-albedo) * NN_HI * Z_cl/Z0 * sigma_d
    EBV   = ta2EBV * tau_a
  endif

# ifdef dust
# else
  if (tau_a .gt. 1e-10) then
    write(*,'("WARNING: MoCaLaTA has been compiled without dust option, but tau_a = ",f5.2)') tau_a
    call holdon(2.)
  endif
# endif
endif

! Set 'intxc' and calculate global x_crit, if necessary
if (trim(x_critType) .eq. 'intrinsic') then
  intxc  = .true.
  maxxc  = .false.
  write(*,*) 'Using intrinsic x_crit function'
elseif (trim(x_critType) .eq. 'global') then    !Calculate x_crit on basis of a*tau_0 from center to face
  if (abs(FF-1d0) .gt. 1d-7) stop "x_critType 'global' only compatible with filling factor = 1."
  intxc  = .false.
  maxxc  = .false.
  Dnu_D_cl   = Doppler(T_cl)
  sigma_0 = sigma(Dnu_D_cl,0d0)
  a       = Dnu_L / (2 * Dnu_D_cl)
  NN_HI   = n_HI_cl * z_0
  tau_0   = NN_HI * sigma_0
  at      = a * tau_0
  x_crit = CoreSkip(at)                         !...in all cells
  write(*,*) 'x_crit set globally to tau_0-dependent value', real(x_crit)
elseif (index(x_critType,'max') .eq. 1) then
  maxxc  = .true.
  intxc  = .true.
  read(x_critType(4:),*) xc_max
  write(*,*) 'Using intrinsic x_crit function, but with a max of', real(xc_max)
elseif (.not.trim(x_critType).eq.'') then
  intxc = .false.
  maxxc = .false.
  read(x_critType,*) x_crit
  write(*,*) 'x_crit set globally to userdefined value', real(x_crit)
else
  write(*,*) "Key word 'x_critType' has no value in " // trim(parfile)
  stop
endif

! Determine initial frequency distribution
if (trim(x_injType) .eq. 'proper') then
  write(*,*) 'Using Voigt function appropriate for cell.'
elseif (index(x_injType,'fwhm') .eq. 1) then
  read(x_injType(5:),*) sigma_V
  sigma_V   = sigma_V / sd2FW
  x_injType = 'fixtgauss'
elseif (index(x_injType,'sdev') .eq. 1) then
  read(x_injType(5:),*) sigma_V
  x_injType = 'fixtgauss'
elseif (index(x_injType,'x') .eq. 1) then
  read(x_injType(2:),*) x_inj
  x_injType = 'x_inj'
else
  write(*,*) 'Unknown x_injType', x_injType
  stop
endif

! Manipulate random number generator
if (TrueRan) then
  call random_seed(size=seedsize)               !Number of elements in seed
  allocate(time(8))
  allocate(seed(seedsize))
  call random_seed(get=seed)                    !Current seed
  call date_and_time(values=time)
  if (sum(time).eq.0 .or. time(1).eq.-huge(time)) stop "System doesn't know date and time."
  where (time .le. 0) time = 1
  call random_seed(put=seed+product(time))      !Add to first 8 indices
  deallocate(time)
  deallocate(seed)
endif

! Direction of initial emission

if (trim(n_init) .eq. 'iso') then
  isoem = .true.
elseif (index(n_init,'(') .eq. 1) then
  isoem = .false.
  c1 = scan(n_init,numsi)
  c2 = scan(n_init(c1:),delim) + c1-1
  read(n_init(c1:c2-1),*) n0(1)

  c1 = scan(n_init(c2:),numsi) + c2-1
  c2 = scan(n_init(c1:),delim) + c1-1
  read(n_init(c1:c2-1),*) n0(2)

  c1 = scan(n_init(c2:),numsi) + c2-1
  c2 = scan(n_init(c1:),delim//')') + c1-1
  read(n_init(c1:c2-1),*) n0(3)

  nhat_0 = n0 / norm(n0)
else
  write(*,*) "Don't understand n_init = '", trim(n_init), "'."
  stop
endif

! Initial position

if (CloudCorr-1d0 .gt. 1d-6) stop 'CloudCorr should be [0, 1] (or equal to -1 to ignore clouds.)'

if (trim(X_init) .eq. 'lumdep') then
  if (CloudCorr.gt.1d-6) stop "CloudCorr > 0, but X_init = 'lumdep'"
  empos  = (/.true., .false., .false., .false., .false./)
elseif (trim(X_init) .eq. 'central') then
  if (CloudCorr.gt.1d-6) stop "CloudCorr > 0, but X_init = 'central'"
  empos  = (/.false., .true., .false., .false., .false./)
elseif (trim(X_init) .eq. 'homo') then
  empos  = (/.false., .false., .true., .false., .false./)
elseif (index(X_init,'(') .eq. 1) then
  if (CloudCorr.gt.1d-6) then
    write(*,*) "CloudCorr > 0, but X_init = "//trim(X_init)
    stop
  endif
  empos  = (/.false., .false., .false., .true., .false./)
  c1 = scan(X_init,numsi)
  c2 = scan(X_init(c1:),delim) + c1-1
  read(X_init(c1:c2-1),*) X_i0(1)

  c1 = scan(X_init(c2:),numsi) + c2-1
  c2 = scan(X_init(c1:),delim) + c1-1
  read(X_init(c1:c2-1),*) X_i0(2)

  c1 = scan(X_init(c2:),numsi) + c2-1
  c2 = scan(X_init(c1:),delim//')') + c1-1
  read(X_init(c1:c2-1),*) X_i0(3)

  X_i0 = D_box * X_i0
elseif (index(X_init,'grad(') .eq. 1) then
  empos  = (/.false., .false., .false., .false., .true./)
  c1 = scan(X_init,numsi)
  c2 = scan(X_init(c1:),delim) + c1-1
  read(X_init(c1:c2-1),*) H_lum(1)

  c1 = scan(X_init(c2:),numsi) + c2-1
  c2 = scan(X_init(c1:),delim) + c1-1
  read(X_init(c1:c2-1),*) H_lum(2)

  c1 = scan(X_init(c2:),numsi) + c2-1
  c2 = scan(X_init(c1:),delim//')') + c1-1
  read(X_init(c1:c2-1),*) H_lum(3)

  H_lum = H_lum * kpc
elseif (index(X_init,'grad') .eq. 1) then
  empos  = (/.false., .false., .false., .false., .true./)
  read(X_init(5:),*) H_lum(1)                   !Put value in first index
  H_lum = H_lum(1) * kpc                        !Copy to the other indicies and convert to kpc
else
  write(*,*) "Don't understand X_init = '", trim(X_init), "'."
  stop
endif

! Lya/FUV photon packets

nu1 = c / BW(2) * 1d8                           !Hz
nu2 = c / BW(1) * 1d8                           !Hz

if (trim(model) .ne. 'data') then
  if (EW_int .gt. 0d0) then
    L_Lya    = SFR * 1.1d42                     !erg/s; Assume Salpeter IMF
    L_FUV_Hz = SFR / 1.4d-28 / (EW_int/76.178)  !erg/s/Hz
  else
    L_Lya    = 0d0                              !erg/s; Assume Salpeter IMF
    L_FUV_Hz = SFR / 1.4d-28                    !erg/s/Hz
  endif

  BaseRes  = (/ni,nj,nk/)

  Q_Lya    = L_Lya / E_Lya                      !ph/s
  Q_FUV    = L_FUV_Hz/h_Pl * (log(nu2)-log(nu1))!ph/s
  Q_tot    = Q_Lya + Q_FUV                      !ph/s
  n_ph(1)  = nint(Q_Lya/Q_tot * n_phtot)        !Total # of Lya photon packets
  n_ph(2)  = n_phtot - n_ph(1)                  !Total # of FUV photon packets
  Lya2tot  = Q_Lya / Q_tot                      !Probability for an emitted photon to be a continuum photon
  L_FUV    = L_FUV_Hz * (nu2-nu1)
  L_tot    = L_Lya + L_FUV
endif

! Allocate arrays

if (trim(model) .ne. 'data') then
  allocate(Dnu_DString(0:n_phases))
  allocate(n_HIString(0:n_phases))
# ifdef dust
  allocate(n_dString(0:n_phases))
# endif
endif

DataDir = trim(DataDir) // '/'
Subdir  = trim(Subdir)  // '/'
ext     = trim(outlabel) // '.bin'
if (len_trim(specin) .gt.0) then
  call getlun(inlun)
  open (inlun,file=trim(DataDir)//trim(Subdir)//trim(specin),status='replace')
endif
if (len_trim(specout).gt.0) then
  call getlun(outlun)
  open (outlun,file=trim(DataDir)//trim(Subdir)//trim(specout),status='replace')
endif

if (trim(X_QSOc).ne.'nil') then
  if (index(X_QSOc,'(') .eq. 1) then
    c1 = scan(X_QSOc,numsi)
    c2 = scan(X_QSOc(c1:),delim) + c1-1
    read(X_QSOc(c1:c2-1),*) X_QSO(1)

    c1 = scan(X_QSOc(c2:),numsi) + c2-1
    c2 = scan(X_QSOc(c1:),delim) + c1-1
    read(X_QSOc(c1:c2-1),*) X_QSO(2)

    c1 = scan(X_QSOc(c2:),numsi) + c2-1
    c2 = scan(X_QSOc(c1:),delim//')') + c1-1
    read(X_QSOc(c1:c2-1),*) X_QSO(3)

    X_QSO = X_QSO*kpc + D_box/2                 !Convert origo from [0] to [R_box] in kpc
  else
    write(*,*) "If QSO sightline is desired, encap position in parentheses:"
    write(*,*) "'"//trim(X_QSOc)//"' => '("//trim(X_QSOc)//")'"
    write(*,*) "For no QSO, set keyword to 'nil'"
    stop
  endif
endif

! More stuff

if (V_out .lt. 0) write(*,*) 'WARNING: V_out < 0 gives a funny profile'

end subroutine ReadInput

!------------------------------------------------------------------------------

subroutine ReadOutCCD(n)

!use ran_state, only: ran_seed !NumRec
use kinds
use ObservationalParameters
use GlobalData
use PhysicalConstants

implicit none
integer::            i,n_CCDs,l
integer,intent(in):: n
character(len=7)::   frame
real(dp)::           Total,Sxp,Sxm,Syp,Sym,Szp,Szm
real(dp)::           n_sofar
real(dp)::           f_esc(n_phTypes),f_abs(n_phTypes),f_escFUV,b

write(*,*)
write(*,*) 'Reading out CCD...'

Total    = 0d0
n_CCDs   = 0
n_sofar  = sum(n_eff)
f_abs    = n_abs / n_eff!n_sofar
f_esc    = 1 - f_abs
f_escFUV = 1 - n_absFUV/n_FUV
b        = f_esc(1) / f_escFUV

write(*,*) '             Lya              FUV'
write(*,*) 'n_eff  ',  real(n_eff)
write(*,*) 'f_abs  ',  real(f_abs)
write(*,*) 'f_esc  ',  real(f_esc(1)), real(f_escFUV), '~=~>  EW boost ~',real(b)!f_esc(1)/f_esc(2))
write(*,*)

call getlun(l)

if (WhichCCDs(1)) then
  open (l, file=trim(DataDir)//trim(Subdir)// 'xp' //trim(ext), form='unformatted',status='replace')
  write(l) n_sofar
  write(l) pix,SpecRes1D,SpecRes2D,n,L_tot,BW,z,d_L,R_obs/kpc
  write(l) obspar
  write(l) specxp
  write(l) CCDxp
  close(l)
  Sxp    = sum(specxp)
  Total  = Total + Sxp
  n_CCDs = n_CCDs + 1
  write(*,*) 'xp:        ', real(Sxp/n_sofar)
endif

if (WhichCCDs(2)) then
  open (l, file=trim(DataDir)//trim(Subdir)// 'xm' //trim(ext), form='unformatted',status='replace')
  write(l) n_sofar
  write(l) pix,SpecRes1D,SpecRes2D,n,L_tot,BW,z,d_L,R_obs/kpc
  write(l) obspar
  write(l) specxm
  write(l) CCDxm
  close(l)
  Sxm    = sum(specxm)
  Total  = Total + Sxm
  n_CCDs = n_CCDs + 1
  write(*,*) 'xm:        ', real(Sxm/n_sofar)
endif

if (WhichCCDs(3)) then
  open (l, file=trim(DataDir)//trim(Subdir)// 'yp' //trim(ext), form='unformatted',status='replace')
  write(l) n_sofar
  write(l) pix,SpecRes1D,SpecRes2D,n,L_tot,BW,z,d_L,R_obs/kpc
  write(l) obspar
  write(l) specyp
  write(l) CCDyp
  close(l)
  Syp    = sum(specyp)
  Total  = Total + Syp
  n_CCDs = n_CCDs + 1
  write(*,*) 'yp:        ', real(Syp/n_sofar)
endif

if (WhichCCDs(4)) then
  open (l, file=trim(DataDir)//trim(Subdir)// 'ym' //trim(ext), form='unformatted',status='replace')
  write(l) n_sofar
  write(l) pix,SpecRes1D,SpecRes2D,n,L_tot,BW,z,d_L,R_obs/kpc
  write(l) obspar
  write(l) specym
  write(l) CCDym
  close(l)
  Sym    = sum(specym)
  Total  = Total + Sym
  n_CCDs = n_CCDs + 1
  write(*,*) 'ym:        ', real(Sym/n_sofar)
endif

if (WhichCCDs(5)) then
  open (l, file=trim(DataDir)//trim(Subdir)// 'zp' //trim(ext), form='unformatted',status='replace')
  write(l) n_sofar
  write(l) pix,SpecRes1D,SpecRes2D,n,L_tot,BW,z,d_L,R_obs/kpc
  write(l) obspar
  write(l) speczp
  write(l) CCDzp
  close(l)
  Szp    = sum(speczp)
  Total  = Total + Szp
  n_CCDs = n_CCDs + 1
  write(*,*) 'zp:        ', real(Szp/n_sofar)
endif

if (WhichCCDs(6)) then
  open (l, file=trim(DataDir)//trim(Subdir)// 'zm' //trim(ext), form='unformatted',status='replace')
  write(l) n_sofar
  write(l) pix,SpecRes1D,SpecRes2D,n,L_tot,BW,z,d_L,R_obs/kpc
  write(l) obspar
  write(l) speczm
  write(l) CCDzm
  close(l)
  Szm    = sum(speczm)
  Total  = Total + Szm
  n_CCDs = n_CCDs + 1
  write(*,*) 'zm:        ', real(Szm/n_sofar)
endif

call cpu_time(t2)

write(*,*) 'Average:   ', real(Total / n_CCDs / n_sofar)
write(*,*) 'f_esc,iso: ', real(f_esc)
write(*,*) 'Secs/ph.:  ', real((t2 - t1) / n_sofar)
write(*,*) 'n_ph,eff:  ', real(n_sofar)

if (len_trim(specin) .gt.0) flush(inlun)
if (len_trim(specout).gt.0) flush(outlun)

end subroutine ReadOutCCD

!-------------------------------------------------------------------------------

subroutine RecoverCoreHistory(absorb)

!Calculate probability that photon would have been absorbed and number of
!scatterings that photon would have experienced, had x_crit-acceleration
!scheme not been invoked, making a bilinear interpolation of tabled values
!coming from Monte Carlo simulations.

use kinds
use CurrentCellParameters
use GlobalData
use DerivedParameters
use iMonteCarlo

implicit none
logical, intent(out):: absorb
integer::              j,k
real(dp)::               y1,y2,y3,y4,t,u,f_esc,R

if (logatr .gt. atA(10)) then
  absorb = .false.
else
  j  = ceiling((logatr - atA(1)) / dat)
  k  = ceiling((x_crit - xcA(1)) / dxc)

  y1 = fescxc(j,k)
  y2 = fescxc(j+1,k)
  y3 = fescxc(j+1,k+1)
  y4 = fescxc(j,k+1)

  t  = (logatr - atA(j)) / (atA(j+1) - atA(j))
  u  = (x_crit - xcA(k)) / (xcA(k+1) - xcA(k))

  f_esc = (1 - t)*(1 - u)*y1 + t*(1 - u)*y2 + t*u*y3 + (1 - t)*u*y4

  call MonteCarlo(R)
  if (R .gt. f_esc) then
    absorb = .true.
    write(*,*) x_crit, 10.**logatr
  else
    absorb = .false.
  endif
endif

end subroutine RecoverCoreHistory

!------------------------------------------------------------------------------

subroutine RotateVector(nhat_i,nhat,nhat_f)

!Rotate unit vector of outgoing photon so that it is given according to
!incident photon.
!nhat_i is direction of incoming photon, nhat is direction of deflection from
!from nhat_i, i.e. direction of outgoing photon is nhat_i lies along the global
!z-axis. nhat_f is true direction of outgoing photon in global coordinates.

use kinds

implicit none
real(dp),intent(in)::  nhat_i(3),nhat(3)
real(dp),intent(out):: nhat_f(3)
real(dp)::             sq,r,cost,sint,cosp,sinp
real(dp)::             xi,yi,zi,xf,yf,zf

xi = nhat_i(1)
yi = nhat_i(2)
zi = nhat_i(3)

cost = nhat(3)                                  !cos(theta)
sint = sqrt(1d0 - nhat(3)**2)                   !sin(theta)
r    = sqrt(nhat(1)**2 + nhat(2)**2)
cosp = nhat(1) / r                              !cos(phi)
sinp = nhat(2) / r                              !sin(phi)

if (1d0 - abs(zi) .ge. 1d-4) then
  sq   = sqrt(1d0 - zi**2)
  xf = sint/sq * (xi*zi*cosp - yi*sinp) + xi*cost
  yf = sint/sq * (yi*zi*cosp + xi*sinp) + yi*cost
  zf = -sint*cosp*sq + zi*cost
else                                            !If nhat_i is close to z-axis
  xf = sint*cosp
  yf = sint*sinp
  zf = zi*cost
endif

nhat_f = (/xf,yf,zf/)

end subroutine RotateVector

!------------------------------------------------------------------------------

subroutine sample(lun,x,nhat,U_bulk,Dnu_D,X_pos,dist,photon)

use kinds
use PhysicalConstants
use ObservationalParameters
use GlobalData

implicit none
integer,intent(in)::  lun,photon
real(dp),intent(in):: x,nhat(3),U_bulk(3),Dnu_D,X_pos(3),dist
real(dp)::            x_lab,lambda_lab

1 format (es12.3,2x,f10.2,f12.4,4x,3f10.5,2x,3f12.3,2x,f12.3,2x,i1)

x_lab      = x + dot_product(nhat,U_bulk)
lambda_lab = c / (x_lab * Dnu_D + nu_0) * 1e8
if (nullifyx) x_lab = 0.                        !To avoid confusion, if x doesn't make sense, set to zero

write(lun,1) w_ph,       & !1
             x_lab,      & !2
             lambda_lab, & !3
             nhat,       & !4:6
             X_pos/kpc,  & !7:9
             dist/kpc,   & !10
             photon        !11

!header: '#     w_ph         x_lab   lambda_lab       nhat1     nhat2    nhat3          X_pos1      X_pos2      X_pos3         dist   photon'

end subroutine sample

!------------------------------------------------------------------------------

subroutine ScatPhoton(x,nhat_i,nhat_f,HostCell,absorb)

!Calculate frequency and direction of outgoing photon, on basis of atomic
!velocity

use kinds
use AtomParameters
use DerivedParameters
use CellStructure
use CurrentCellParameters
use PhysicalConstants
use GlobalData
use iMonteCarlo

implicit none
real(dp), intent(inout):: x
real(dp), intent(in)::    nhat_i(3)
real(dp), intent(out)::   nhat_f(3)
logical, intent(out)::    absorb
type(cell), pointer::     HostCell
real(dp)::                R,R1,R2,R3,P_H,k_HI,k_d

absorb = .false.

k_HI = n_HI * sigma_x                           !H scattering coefficient
k_d  = n_d  * sigma_d                           !Dust interaction coefficient
P_H  = k_HI / (k_HI + k_d)                      !Probability of being scattered
                                                !by a hydrogen atom
!write(1,*) x, 1d0/(1d0-P_H)

call MonteCarlo(R)

if (R .le. P_H) then
  call RandomUparal(x)                          !Velocity components parallel
  call RandomUperp(x,absorb)                    !and perpendicular to incident photon in terms of Doppler velocity X sqrt2 IN REF. FRAME OF U_BULK

  if (.not.absorb) then
    call ChangeBasis(nhat_i)                    !Express total velocity U in cell frame.
                                                !Calc. also unit vectors perpendicular to u_II
!   if (abs(x) .lt. x_cw) then !Revert this
      call IsoDir(nhat_f)                       !Scatter core photons isotropically
!   else
!     call Dipole(nhat_i,uhat_1,uhat_2,nhat_f)  !Scatter wing photons with dipole distribution
!   endif

    x = x - u_II + dot_product(nhat_f,U)        !Lorentz transform frequency due to atom velocity
    if (recoil) x = x + g_rec*(dot_product(nhat_i,nhat_f) - 1d0) !Add recoil
  endif
else
  call MonteCarlo(R)
  if (R .le. albedo) then                       !Photon is scattered
    absorb = .false.
    call HenyeyGreenstein(nhat_i,g,nhat_f)      !Scatter anisotropically, but coherently
  else                                          !Photon is absorbed
    absorb = .true.
  endif
endif

end subroutine ScatPhoton

!------------------------------------------------------------------------------

subroutine Thomson(ExitVec)

!Calculate direction vector of photon escaping a cube *in the xp-direction*
!according to Thomson distribution: P(mu) = 6/7 (mu + 2mu^2)

use kinds
use GlobalData
use iMonteCarlo

implicit none
real(dp), intent(out):: ExitVec(3)
real(dp)::              R2(2),rsq,frevert,n_1,n_2

do
  call MonteCarlo(R2)

  if (R2(2) .lt. .333333333 * (R2(1) + 2*R2(1)**2)) then
    ExitVec(1) = R2(1)

    do
      call MonteCarlo(R2)
      R2  = 2*R2 - 1
      rsq = sum(R2**2)
      if (rsq .le. 1) exit
    enddo

    frevert    = sqrt((1 - ExitVec(1)**2) / rsq) !Make x,y-components touch
    ExitVec(2) = frevert * R2(1)                  !unit sphere
    ExitVec(3) = frevert * R2(2)

    exit
  endif
enddo

end subroutine Thomson

!------------------------------------------------------------------------------

subroutine weight(photon)

! use kinds
! use GlobalData
! use ObservationalParameters
! use iMonteCarlo
! 
! implicit none
! integer,intent(in):: photon
! integer::            n
! real(dp)::           R
! 
! n = n_packets(photon)
! 
! if     (emdimL.eq.4 .and. photon.eq.1 .and. trim(model).eq.'data') then
!   w_ph = InitLya(n,4)
! elseif (emdimF.eq.4 .and. photon.eq.2 .and. trim(model).eq.'data') then
!   w_ph = InitFUV(n,4)
! else
!   w_ph = 1d0
! endif
! 
! n_eff(photon) = n_eff(photon) + w_ph

end subroutine weight

!------------------------------------------------------------------------------

subroutine which(photon)

!Determine which species of photon is emitted. The probability of emitting a
!given photon is proportional to the total number of photon packets of that
!species, until all those photons have been launched, after which only the
!remaining photons may be emitted. In this way, at all times during execution
!of the code, the spectrum will be representative of the final spectrum, except
!for smaller number statistics.

use kinds
use GlobalData
use ObservationalParameters
use iMonteCarlo

implicit none
integer,intent(out):: photon
real(dp)::            R

if (all(n_packets .lt. n_ph)) then
  call MonteCarlo(R)
  if (R .lt. Lya2tot) then
    photon = 1
  else
    photon = 2
  endif
elseif (n_packets(1).lt.n_ph(1) .and. n_packets(2).ge.n_ph(2)) then
  photon = 1
elseif (n_packets(1).ge.n_ph(1) .and. n_packets(2).lt.n_ph(2)) then
  photon = 2
endif

n_packets(photon) = n_packets(photon) + 1

end subroutine which

!------------------------------------------------------------------------------

recursive subroutine AssignParameters(CurrentCell)

use DataArrays
use CellStructure
use iv2u
use PhysicalConstants
use GlobalData

implicit none
type(Cell), target:: CurrentCell
real(dp)::           X(3),d,nhat(3),v(3),u(3),norm,XwrtC(3),Outflow,gradfac
integer::            ii,jj,kk,p,clnum

# ifdef AMR
# ifdef multi
if (CurrentCell%Refined) then
  do ii = 1,2
    do jj = 1,2
      do kk = 1,2
        call AssignParameters(CurrentCell%Child(ii,jj,kk))
      enddo
    enddo
  enddo
else
  p = CurrentCell%phase
  if (p .eq. 1) then                            !Cloud
    clnum = CurrentCell%clnum
    X     = Clouds(clnum,1:3)
    XwrtC = X - R_box
    d     = norm(XwrtC)
    nhat  = XwrtC / d
    v     = Outflow(d) * nhat
    v     = v + CloudVelDisp(clnum,:)
  else                                          !ICM
    X     = CurrentCell%C_pos
    XwrtC = X - R_box
    d     = norm(XwrtC)
    nhat  = XwrtC / d
    v     = Outflow(d) * nhat
  endif
  u = v2u(v,Dnu_DString(p))

  gradfac = exp(-(d/H1)**exp1)

  CurrentCell%n_HI   = n_HIString(p) * gradfac
  CurrentCell%Dnu_D  = Dnu_DString(p)
  CurrentCell%U_bulk = u
  CurrentCell%n_d    = n_dString(p)  * gradfac
endif
# endif
# endif

end subroutine AssignParameters

!------------------------------------------------------------------------------

recursive subroutine BuildFlorent(CurrentCell,CurrentLevel)

use kinds
use PhysicalConstants
use CellStructure
use GlobalData
use DataArrays
use iv2u
use iMonteCarlo
use igasdev

implicit none
type(Cell), target::     CurrentCell
integer, intent(in)::    CurrentLevel
integer::                ii,jj,kk,cl,p
!real(sp)::                 Dnu_D,n_HI,n_d
real(dp)::                 zz,X(3),norm,nhat(3),d,R,v(3),u(3),XwrtC(3),R3(3)

X     = CurrentCell%C_pos
XwrtC = X - R_box
d     = norm(XwrtC)
if (d.ge.r_inner .and. d.le.r_outer) then
  call MonteCarlo(R)
  p = merge(1,2,R.lt.FF)                        !1 if cell is in shell and
else                                            !  contains a cloud; 2 otherwise
  p = 0                                         !0 if outside shell (void)
endif

nhat = (X-R_box) / d                            !Unit vector from box center to cell center
v    = V_out * nhat                             !Cell velocity; km/s
if (trim(Vprof) .eq. 'linear') v = v * d/R_box
call gasdev(R3)
v = v + R3*sigV_cl
u = v2u(v,Dnu_DString(p))                       !Cell velocity; v_th

# ifdef AMR
  if (LevelString(i_cell) .eq. CurrentLevel) then
    CurrentCell%Level  = LevelString(p)
    CurrentCell%n_HI   = n_HIString(p)
    CurrentCell%Dnu_D  = Dnu_DString(p)
    CurrentCell%U_bulk = real(u)
#   ifdef dust
    CurrentCell%n_d    = n_dString(p)
#   endif
#   ifdef multi
    CurrentCell%phase  = p
#   endif
  elseif (LevelString(i_cell) .gt. CurrentLevel) then
    CurrentCell%Refined = .true.
    CurrentCell%Level   = CurrentLevel
    allocate(CurrentCell%Child(2,2,2))

    i_cell = i_cell - 1
    do ii = 1,2
      do jj = 1,2
        do kk = 1,2
          CurrentCell%Child(ii,jj,kk)%Refined = .false.
          nullify(CurrentCell%Child(ii,jj,kk)%Child)
          call BuildFlorent(CurrentCell%Child(ii,jj,kk), CurrentLevel+1)
        enddo
      enddo
    enddo
  else
    write(*,*) 'Aaargh! Error in levels', i_cell, LevelString(i_cell), CurrentLevel
    stop
  endif
# else
  ! if (p.eq.1) write(32,*) CurrentCell%C_pos/kpc
    CurrentCell%n_HI   = n_HIString(p)
    CurrentCell%Dnu_D  = Dnu_DString(p)
    CurrentCell%U_bulk = real(u)
#   ifdef dust
    CurrentCell%n_d    = n_dString(p)
#   endif
#   ifdef multi
    CurrentCell%phase  = p
#   endif
# endif

end subroutine BuildFlorent

!------------------------------------------------------------------------------

recursive subroutine BuildMultiphase(CurrentCell,CurrentLevel)

use kinds
use PhysicalConstants
use CellStructure
use GlobalData
use DataArrays
use iv2u

implicit none
type(Cell), target::     CurrentCell
integer, intent(in)::    CurrentLevel
integer::                ii,jj,kk,cl,p
!real(sp)::                 Dnu_D,n_HI,n_d
real(dp)::                 X(3),X_cl(3),nhat(3),d_cl,v(3),u(3),d

i_cell = i_cell + 1

# ifdef AMR
# ifdef multi
  if (LevelString(i_cell) .eq. CurrentLevel) then
    CurrentCell%Level = LevelString(i_cell)
    CurrentCell%phase = PhaseString(i_cell)
    CurrentCell%clnum = CloudString(i_cell)
  elseif (LevelString(i_cell) .gt. CurrentLevel) then
    CurrentCell%Refined = .true.
    CurrentCell%Level   = CurrentLevel
    allocate(CurrentCell%Child(2,2,2))

    i_cell = i_cell - 1
    do ii = 1,2
      do jj = 1,2
        do kk = 1,2
          CurrentCell%Child(ii,jj,kk)%Refined = .false.
          nullify(CurrentCell%Child(ii,jj,kk)%Child)
          call BuildMultiphase(CurrentCell%Child(ii,jj,kk), CurrentLevel+1)
        enddo
      enddo
    enddo
  else
    write(*,*) 'Aaargh! Error in levels', i_cell, LevelString(i_cell), CurrentLevel
    stop
  endif
# endif
# endif

end subroutine BuildMultiphase

!------------------------------------------------------------------------------

recursive subroutine BuildSemireal(CurrentCell,CurrentLevel)

use kinds
use PhysicalConstants
use CellStructure
use GlobalData
use DataArrays
use iv2u
use igasdev

implicit none
type(Cell), target::     CurrentCell
integer, intent(in)::    CurrentLevel
integer::                ii,jj,kk,cl,p,ijk_c(3),ijk_cl(3)
!real(sp)::                 Dnu_D,n_HI,n_d
real(dp)::               zz,norm,d_cl,d,gradfac
real(dp),dimension(3)::  R3,X,X_cl,nhat,v,dv,u,XwrtC(3)
real(dp)::               left,right
logical::                GoodEnough

X     = CurrentCell%C_pos
XwrtC = X - R_box
d     = norm(XwrtC)
p     = merge(2,0,d.le.r_gal)                   !Default phase: "ICM" if inside galaxy, "void" if outside
dv    = 0.

if (p.eq.2) then
  do cl = 1,N_cl
    X_cl = Clouds(cl,:)                         !Position of cloud cl's center
    d_cl = sqrt(sum((X-X_cl)**2))               !Distance of cell to cloud center
    !--------------------------------------------------------------------------
    !  Make sure cell is a cloud if hosting a cloud center, even if d_cl<r_cl.!
    !                 !!!! THIS ONLY WORKS FOR NON-AMR !!!!!                  !
                           ijk_c      = ceiling(X/dx0)                        !
                           ijk_cl     = ceiling(X_cl/dx0)                     !
                           GoodEnough = all(ijk_c .eq. ijk_cl)                !
      ! For AMR, use instead loop over Clouds-array with LocateHostCell.        !
    !--------------------------------------------------------------------------
    if (d_cl.le.r_cl .or. GoodEnough) then      !If inside cloud
      p  = 1                                    !Cloud phase
      dv = CloudVelDisp(cl,:)
      if (len_trim(Cloudfile) .gt. 0) write(cloudlun,*) X/kpc
      exit
    endif
  enddo
endif

if (d.ge.r_inner .and. d.le.r_outer) then
  left  = cosOAb2sq * dot_product((X-R_box)/kpc,(X-R_box)/kpc)
  right = dot_product((X-R_box)/kpc,n_jet)**2
  if (left .le. right) then
    p = 3                  !Set shell phase if inside cone
    if (len_trim(Windfile) .gt. 0) write(windlun,*) X/kpc
  endif
endif
! zz = 3*D_box/8 + X(2)/4 !Uncomment to cut upper half of shell, inclined 14 deg.
! if (X(3) .gt. zz) p = 0 !

nhat = XwrtC / norm(XwrtC)                      !Unit vector from box center to cell center
if (p .eq. 3) then
  v = V_out*nhat                                !
  if (trim(Vprof) .eq. 'linear') v = v * d/R_box
! if (.not.uniV) v = v * (d-r_inner)/(R_box-r_inner)
else                                            !Cell velocity; km/s
  v = V_in * nhat                               !
  if (trim(Vprof) .eq. 'linear') v = v * d/R_box
  v = v + dv
endif
u = v2u(v,Dnu_DString(p))                       !Cell velocity; v_th

if (p .eq. 2) then
  gradfac = exp(-(d/H1)**exp1)
else
  gradfac = 1d0
endif

# ifdef AMR
  if (LevelString(i_cell) .eq. CurrentLevel) then
    CurrentCell%Level  = LevelString(p)
    CurrentCell%n_HI   = n_HIString(p) * gradfac
    CurrentCell%Dnu_D  = Dnu_DString(p)
    CurrentCell%U_bulk = real(u)
#   ifdef dust
    CurrentCell%n_d    = n_dString(p) * gradfac
#   endif
#   ifdef multi
    CurrentCell%phase  = p
#   endif
  elseif (LevelString(i_cell) .gt. CurrentLevel) then
    CurrentCell%Refined = .true.
    CurrentCell%Level   = CurrentLevel
    allocate(CurrentCell%Child(2,2,2))

    i_cell = i_cell - 1
    do ii = 1,2
      do jj = 1,2
        do kk = 1,2
          CurrentCell%Child(ii,jj,kk)%Refined = .false.
          nullify(CurrentCell%Child(ii,jj,kk)%Child)
          call BuildSemireal(CurrentCell%Child(ii,jj,kk), CurrentLevel+1)
        enddo
      enddo
    enddo
  else
    write(*,*) 'Aaargh! Error in levels', i_cell, LevelString(i_cell), CurrentLevel
    stop
  endif
# else
    CurrentCell%n_HI   = n_HIString(p) * gradfac
    CurrentCell%Dnu_D  = Dnu_DString(p)
    CurrentCell%U_bulk = real(u)
#   ifdef dust
    CurrentCell%n_d    = n_dString(p) * gradfac
#   endif
#   ifdef multi
    CurrentCell%phase  = p
#   endif
# endif

end subroutine BuildSemireal

!------------------------------------------------------------------------------

recursive subroutine BuildNestedUniverse(CurrentCell,CurrentLevel)

use CellStructure
use GlobalData
use DataArrays

implicit none
type(Cell), target::     CurrentCell
integer, intent(IN)::    CurrentLevel
integer::                ii,jj,kk

i_cell = i_cell + 1

# ifdef AMR
  if (LevelString(i_cell) .eq. CurrentLevel) then
    CurrentCell%Level  = LevelString(i_cell)
    CurrentCell%n_HI   = n_HIString(i_cell)
    CurrentCell%Dnu_D  = Dnu_DString(i_cell)
    CurrentCell%U_bulk = (/U_xString(i_cell), &
                           U_yString(i_cell), &
                           U_zString(i_cell)/)
#   ifdef dust
    CurrentCell%n_d    = n_dString(i_cell)
#   endif
  elseif (LevelString(i_cell) .gt. CurrentLevel) then
    CurrentCell%Refined = .true.
    CurrentCell%Level   = CurrentLevel
    allocate(CurrentCell%Child(2,2,2))

    i_cell = i_cell - 1
    do ii = 1,2
      do jj = 1,2
        do kk = 1,2
          CurrentCell%Child(ii,jj,kk)%Refined = .false.
          nullify(CurrentCell%Child(ii,jj,kk)%Child)
          call BuildNestedUniverse(CurrentCell%Child(ii,jj,kk), CurrentLevel+1)
        enddo
      enddo
    enddo
  else
    write(*,*) 'Aaargh! Error in levels', i_cell, LevelString(i_cell), CurrentLevel
    stop
  endif
# else
    CurrentCell%n_HI   = n_HIString(i_cell)
    CurrentCell%Dnu_D  = Dnu_DString(i_cell)
    CurrentCell%U_bulk = (/U_xString(i_cell), &
                           U_yString(i_cell), &
                           U_zString(i_cell)/)
#   ifdef dust
    CurrentCell%n_d    = n_dString(i_cell)
#   endif
# endif

end subroutine BuildNestedUniverse

!------------------------------------------------------------------------------

recursive subroutine BuildShell(CurrentCell,CurrentLevel)

use kinds
use PhysicalConstants
use CellStructure
use GlobalData
use DataArrays
use iv2u

implicit none
type(Cell), target::     CurrentCell
integer, intent(in)::    CurrentLevel
integer::                ii,jj,kk,cl,p
!real(sp)::                 Dnu_D,n_HI,n_d
real(dp)::                 zz,X(3),norm,nhat(3),d,v(3),u(3),XwrtC(3)

X = CurrentCell%C_pos
XwrtC = X - R_box
d = norm(XwrtC)
p = merge(1,0,d.ge.r_inner .and. d.le.r_outer)!1 if in shell, 0 if outside

nhat = (X-R_box) / d                            !Unit vector from box center to cell center
v    = V_out * nhat                             !Cell velocity; km/s
if (trim(Vprof) .eq. 'linear') v = v * d/R_box
u = v2u(v,Dnu_DString(p))                     !Cell velocity; v_th

# ifdef AMR
  if (LevelString(i_cell) .eq. CurrentLevel) then
    CurrentCell%Level  = LevelString(p)
    CurrentCell%n_HI   = n_HIString(p)
    CurrentCell%Dnu_D  = Dnu_DString(p)
    CurrentCell%U_bulk = real(u)
#   ifdef dust
    CurrentCell%n_d    = n_dString(p)
#   endif
#   ifdef multi
    CurrentCell%phase  = p
#   endif
  elseif (LevelString(i_cell) .gt. CurrentLevel) then
    CurrentCell%Refined = .true.
    CurrentCell%Level   = CurrentLevel
    allocate(CurrentCell%Child(2,2,2))

    i_cell = i_cell - 1
    do ii = 1,2
      do jj = 1,2
        do kk = 1,2
          CurrentCell%Child(ii,jj,kk)%Refined = .false.
          nullify(CurrentCell%Child(ii,jj,kk)%Child)
          call BuildShell(CurrentCell%Child(ii,jj,kk), CurrentLevel+1)
        enddo
      enddo
    enddo
  else
    write(*,*) 'Aaargh! Error in levels', i_cell, LevelString(i_cell), CurrentLevel
    stop
  endif
# else
!   if (p.eq.1 .and.CurrentCell%C_pos(1).gt.R_box) write(32,*) d/kpc,CurrentCell%C_pos/kpc
    CurrentCell%n_HI   = n_HIString(p)
    CurrentCell%Dnu_D  = Dnu_DString(p)
    CurrentCell%U_bulk = real(u)
#   ifdef dust
    CurrentCell%n_d    = n_dString(p)
#   endif
#   ifdef multi
    CurrentCell%phase  = p
#   endif
# endif

end subroutine BuildShell

!------------------------------------------------------------------------------

recursive subroutine BuildSlab(CurrentCell,CurrentLevel)

use kinds
use PhysicalConstants
use CellStructure
use GlobalData
use DataArrays

implicit none
type(Cell), target::     CurrentCell
integer, intent(in)::    CurrentLevel
integer::                ii,jj,kk,cl,p
!real(sp)::                 Dnu_D,n_HI,n_d
real(dp)::                 zz,X_cl(3),nhat(3),d_cl,v(3),u(3)

p = merge(1,0,abs(CurrentCell%C_pos(3)-R_box) .lt. z_0)!1 if in slab, 0 if outside

# ifdef AMR
  if (LevelString(i_cell) .eq. CurrentLevel) then
    CurrentCell%Level  = LevelString(p)
    CurrentCell%n_HI   = n_HIString(p)
    CurrentCell%Dnu_D  = Dnu_DString(p)
    CurrentCell%U_bulk = (/0., 0., U_zString(p)/)
#   ifdef dust
    CurrentCell%n_d    = n_dString(p)
#   endif
#   ifdef multi
    CurrentCell%phase  = p
#   endif
  elseif (LevelString(i_cell) .gt. CurrentLevel) then
    CurrentCell%Refined = .true.
    CurrentCell%Level   = CurrentLevel
    allocate(CurrentCell%Child(2,2,2))

    i_cell = i_cell - 1
    do ii = 1,2
      do jj = 1,2
        do kk = 1,2
          CurrentCell%Child(ii,jj,kk)%Refined = .false.
          nullify(CurrentCell%Child(ii,jj,kk)%Child)
          call BuildSlab(CurrentCell%Child(ii,jj,kk), CurrentLevel+1)
        enddo
      enddo
    enddo
  else
    write(*,*) 'Aaargh! Error in levels', i_cell, LevelString(i_cell), CurrentLevel
    stop
  endif
# else
    CurrentCell%n_HI   = n_HIString(p)
    CurrentCell%Dnu_D  = Dnu_DString(p)
    CurrentCell%U_bulk = (/0., 0., U_zString(p)/)
#   ifdef dust
    CurrentCell%n_d    = n_dString(p)
#   endif
#   ifdef multi
    CurrentCell%phase  = p
#   endif
# endif

end subroutine BuildSlab

!------------------------------------------------------------------------------

recursive subroutine CalcFF(CurrentCell)

!Count number of cells of a given phase, weighted by the fractional volume
!they occupy (filling factor).
!
!NOTE: This subroutine only works for a multiphase medium.

use CellStructure
use GlobalData
use PhysicalConstants

implicit none
type(Cell), target:: CurrentCell
integer:: ii,jj,kk,p,L

# ifdef multi
#   ifdef AMR
      if (CurrentCell%Refined) then
        do ii = 1,2
          do jj = 1,2
            do kk = 1,2
              call CalcFF(CurrentCell%Child(ii,jj,kk))
            enddo
          enddo
        enddo
      else
        p = CurrentCell%phase
        L = CurrentCell%Level
        FFcounter(p) = FFcounter(p) + (1d0/2d0**(3d0*L))
      ! if (p .eq. 1) write(1,*) (CurrentCell%C_pos-R_box) / kpc
      endif
#   else
      p = CurrentCell%phase
      FFcounter(p) = FFcounter(p) + 1d0
#   endif
# endif

end subroutine CalcFF

!------------------------------------------------------------------------------

recursive subroutine DigDeeper(X_pos,CurrentCell)

!Check if cell is refined. If not, this is the host cell. If so, go one level
!deeper.

use kinds
use CellStructure
use GlobalData

implicit none
real(dp), intent(IN)::  X_pos(3)
type(Cell), pointer:: CurrentCell               !intent(INOUT)
integer::             ii,jj,kk

# ifdef AMR
ii = merge(1,2,X_pos(1) .lt. CurrentCell%C_pos(1))
jj = merge(1,2,X_pos(2) .lt. CurrentCell%C_pos(2))
kk = merge(1,2,X_pos(3) .lt. CurrentCell%C_pos(3))

CurrentCell => CurrentCell%Child(ii,jj,kk)

if (CurrentCell%Refined) call DigDeeper(X_pos,CurrentCell)
# endif

end subroutine DigDeeper

!------------------------------------------------------------------------------

recursive subroutine EveryoneShouldHaveACenter(CurrentCell,i,j,k,ParentPos)

!Assign center to all cells, even if refined.

use kinds
use CellStructure
use GlobalData
use PhysicalConstants

implicit none
type(Cell),target::             CurrentCell
integer,intent(in)::            i,j,k
real(kind=posKind),intent(in):: ParentPos(3)
integer::                       indexico(3),ii,jj,kk
real(dp)::                      R

indexico = (/i,j,k/)

# ifdef AMR
  if (CurrentCell%Level .eq. 0) then
    CurrentCell%C_pos = (indexico - .5d0) * dx(CurrentCell%Level)
  else
    CurrentCell%C_pos = ParentPos &
                      + 2*(indexico - 1.5d0) & ! -1 or +1
                      * dxH(CurrentCell%Level)
  endif

  if (CurrentCell%Refined) then
    do ii = 1,2
      do jj = 1,2
        do kk = 1,2
          call EveryoneShouldHaveACenter(CurrentCell%Child(ii,jj,kk),ii,jj,kk,CurrentCell%C_pos)
        enddo
      enddo
    enddo
  endif

# else
    CurrentCell%C_pos = (indexico - .5d0) * dx(0)
# endif

end subroutine EveryoneShouldHaveACenter

!------------------------------------------------------------------------------

recursive subroutine SeekAndDestroy(CurrentCell)

use DataArrays
use CellStructure
use iv2u
use PhysicalConstants
use GlobalData
use iTemperature
use iPrintCell

implicit none
type(Cell), target:: CurrentCell
type(Cell), pointer:: WallaCell
real(dp)::            X(3),d,nhat(3),v(3),u(3),norm,n_HI,T
integer::            ii, jj, kk

i_cell = i_cell + 1

# ifdef AMR
if (CurrentCell%Refined) then
  do ii = 1,2
    do jj = 1,2
      do kk = 1,2
        call SeekAndDestroy(CurrentCell%Child(ii,jj,kk))
      enddo
    enddo
  enddo
else
# endif
  i_cell = i_cell - 1
  X      = CurrentCell%C_pos - R_box
  d      = norm(X)
  nhat   = X / d
  u      = CurrentCell%U_bulk
  v      = CurrentCell%U_bulk*CurrentCell%Dnu_D*c / (nu_0*1e5)
! write(CurrentCell%phase+10,'(6es12.4)') d/kpc, norm(v), X/kpc, CurrentCell%n_HI
! write(CurrentCell%phase+100,*) X/kpc
# ifdef AMR
endif
# endif

end subroutine SeekAndDestroy

!------------------------------------------------------------------------------

!==============================  Functions  ===================================

function CoreSkip(at)

!Critical x-value for core-skipping acceleration scheme

use kinds

implicit none
real(dp):: CoreSkip,at

if (at .lt. 1.) then
  CoreSkip = 0.
elseif (at.ge.1  .and. at.lt.60) then
  CoreSkip = .02 * exp((0.6*log(at))**1.2)
else
  CoreSkip = .02 * exp((1.4*log(at))**0.6)
endif

end function CoreSkip
!------------------------------------------------------------------------------

function Doppler_s(T)

!Doppler frequency width

use kinds
use PhysicalConstants
implicit none
real(sp),intent(in):: T
real(sp)::            Doppler_s
Doppler_s = sqrt(2.*k_B*T / m_H) * nu_0/c
end function Doppler_s
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

function Doppler_sv(T)
use kinds
use PhysicalConstants
implicit none
real(sp),intent(in):: T(:)
real(sp)::            Doppler_sv(size(T))
Doppler_sv = sqrt(2.*k_B*T / m_H) * nu_0/c
end function Doppler_sv
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

function Doppler_d(T)
use kinds
use PhysicalConstants
implicit none
real(dp),intent(in):: T
real(dp)::            Doppler_d
Doppler_d = sqrt(2.*k_B*T / m_H) * nu_0/c
end function Doppler_d
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

function Doppler_dv(T)
use kinds
use PhysicalConstants
implicit none
real(dp),intent(in):: T(:)
real(dp)::            Doppler_dv(size(T))
Doppler_dv = sqrt(2.*k_B*T / m_H) * nu_0/c
end function Doppler_dv
!------------------------------------------------------------------------------

function H_z(z)

!Hubble parameter as a function of redshift DIVIDED BY present Hubble constant,
!.e. DIMENSIONLESS

use kinds
use PhysicalConstants

implicit none
real(dp):: H_z,z

H_z = (1d0 + z) * sqrt(1d0 + Omega_M*z + Omega_L*(1d0/(1d0 + z)**2 - 1d0))

end function
!------------------------------------------------------------------------------

function LumDist(z)

!Calulates luminosity distance in kpc as a function of redshift

use kinds
use PhysicalConstants

implicit none
real(dp),intent(in)::  z
integer::              n,j
real(dp)::             LumDist,dz,zz,f_odd,f_even,f_0,f_n
real(dp)::             H_z,I                    !H_z is in fact H(z)/H_0

n = 10000                                       !# of intervals; must be even
dz = z / n                                      !Interval width

f_odd = 0d0
do j=1,n-1,2
  zz    = j*dz
  f_odd = f_odd + 1 / H_z(zz)                   !Sum of terms with index odd
enddo

f_even = 0d0
do j=2,n-2,2
  zz    = j*dz
  f_even = f_even + 1 / H_z(zz)                 !Sum of terms with index even
enddo

f_0 = 1 / H_z(0d0)                              !Endpoint values
f_n = 1 / H_z(z)                                !

I   = dz/3 * (f_0 + 4*f_odd + 2*f_even + f_n)   !Evaluate integral using Simpson's rule

LumDist = c/H_0 * (1+z) * I                     !Luminosity distance in kpc

end function LumDist
!------------------------------------------------------------------------------

function percentile(n,arr,p)

!Returns the p'th percentile of the numbers in arr. For p = 50, the function
!returns the median.
!WARNING: No interpolation between elements is made => imprecise for small n.
!
!TODO:
! 1) Make p an array, so that several percentiles can be returned
!    simultaneously, thus minimizing # of sortings.
! 2) Make interpolation between arr-values encompassing the percentile:
!    percentile = lo + r*(hi-lo).
! 3) Make p optional, with default p = 50. Requires explicit interface.
!    This has been attempted and worked with ifort, g95, and f95, but with
!    gfortran the code simply stalls after entering function, but before
!    reaching first statement.

use kinds

implicit none
integer,intent(in)::  n                         !# of elements in arr
real(dp),intent(in):: arr(n)                    !Array of numbers to be sorted
real(dp),intent(in):: p                         !Percentile [0,100]
real(dp)::            percentile,temparr(n),pp
integer::             h
!optional::            p                         !If non-present, p = 50.
!NOTE: If 'optional' is implemented, make explicit interface

! if (present(p)) then
    pp = .01 * p
! else
!   pp = 0.5d0
! endif

h = nint(pp * n)

if (h .lt. 1) h = 1
if (h .gt. n) h = n

temparr = arr
call Quicksort(n,temparr)

percentile = temparr(h)

end function
!------------------------------------------------------------------------------

function norm(vec)

use kinds

implicit none
real(dp):: vec(3),norm

norm = sqrt(sum(vec**2))

end function
!------------------------------------------------------------------------------

function Outflow(r)

!Outflow speed (in km/s) as a function of distance (in cm) from the galactic
!center.

use kinds
use GlobalData
use PhysicalConstants

implicit none
real(dp):: Outflow,r,r_breakout,alpha,oma,A,term1,term2

!r_breakout and alpha should be given in in-file, and the rest of this block
!should be calc'd in ReadInput.
r_breakout = r_gal/10.!1d0 * kpc                          !Distance from center where wind starts
r_breakout = max(r_breakout,dxWeeny)            !Avoid infinities
alpha      = 1.5                                ![1,15,1.95] (Steidel, 2010, ApJ, 717, 289)
if (alpha .le. 1d0) stop 'Not ready for alpha =< 1.'
oma   = 1 - alpha
V_inf = V_out / sqrt(1-(r_gal/r_breakout)**oma) !V_out is at r_gal; V_inf is at infty
A     = -oma*V_inf**2 / (2*r_breakout**oma)     !"A constant that sets V_inf"
term1 = -2 * A / oma

if (r .ge. r_breakout) then
  select case (trim(Vprof))
    case ('momentum')
      term2   = r_breakout**oma - r**oma
      Outflow = sign(sqrt(term1*term2),V_inf)
    case ('linear') !Linearly from 0 to V_out as r goes from r_breakout to r_gal
      Outflow = V_out * (r-r_breakout) / (r_gal-r_breakout)
    case ('constant')
      Outflow = V_out
    case default
      write(*,*) "Unknown wind profile '"//trim(Vprof)//"'."
      stop
  end select
else
  Outflow = 0
endif

end function
!------------------------------------------------------------------------------

function Temperature_s(Dnu_D)

!Inverse function: Doppler -> Temperature

use kinds
use PhysicalConstants
implicit none
real(sp),intent(in):: Dnu_D
real(sp)::            Temperature_s
Temperature_s = m_H/2/k_B * (Dnu_D * c/nu_0)**2
end function Temperature_s
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

function Temperature_sv(Dnu_D)
use kinds
use PhysicalConstants
implicit none
real(sp),intent(in):: Dnu_D(:)
real(sp)::            Temperature_sv(size(Dnu_D))
Temperature_sv = m_H/2/k_B * (Dnu_D * c/nu_0)**2
end function Temperature_sv
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

function Temperature_d(Dnu_D)
use kinds
use PhysicalConstants
implicit none
real(dp),intent(in):: Dnu_D
real(dp)::            Temperature_d
Temperature_d = m_H/2/k_B * (Dnu_D * c/nu_0)**2
end function Temperature_d
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

function Temperature_dv(Dnu_D)
use kinds
use PhysicalConstants
implicit none
real(dp),intent(in):: Dnu_D(:)
real(dp)::            Temperature_dv(size(Dnu_D))
Temperature_dv = m_H/2/k_B * (Dnu_D * c/nu_0)**2
end function Temperature_dv
!------------------------------------------------------------------------------

function sigma(Dnu_D,x)

!Lyman alpha cross section

use kinds
use PhysicalConstants

implicit none
real(dp):: x,H,Voigt,Dnu_D,sigma

H     = Voigt(Dnu_D,x)
sigma = f_12 * sqpi * e**2 / (m_e * c * Dnu_D) * H

end function sigma
!------------------------------------------------------------------------------

function Voigt(Dnu_D,x)

!Voigt function, analytical approximation (Tasitsiomi, 2006)

use kinds
use PhysicalConstants

implicit none
real(dp):: Dnu_D,a,x,xt,z,q,Voigt

a  = Dnu_L / (2.*Dnu_D)
xt = x**2
z = (xt - .855) / (xt + 3.42)

if (z .le. 0.) then
  q = 0.
else
  q = z * (1. + 21./xt) * a / (pi * (xt+1.)) &
    * (.1117 + z*(4.421 + z*(-9.207 + 5.674*z)))
endif

Voigt = sqpi * q + exp(-xt)

end function Voigt
!------------------------------------------------------------------------------

function v2u_1s2s(v,Dnu_D)

!Generic. Elemental for 1D arrays.
!
!Converts a velocity in km/s to the dimensionless quantity "velocity in terms
!of velocity Doppler width".
!
!The function is called by the generic name "v2u", and the arguments may be
!scalar or 1D vectors, of either single or double precision, in any combination.
!
!The result will have kind = max(kind1,kind2), and rank = max(rank1,rank2).

use kinds
use PhysicalConstants
implicit none
real(sp),intent(in):: v,Dnu_D
real(sp)::            v2u_1s2s

v2u_1s2s = v*1e5 * nu_0 / (Dnu_D * c)

end function v2u_1s2s
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

function v2u_1d2d(v,Dnu_D)
use kinds
use PhysicalConstants
implicit none
real(dp),intent(in):: v,Dnu_D
real(dp)::            v2u_1d2d

v2u_1d2d = v*1d5 * nu_0 / (Dnu_D * c)

end function v2u_1d2d
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

function v2u_1sv2sv(v,Dnu_D)
use kinds
use PhysicalConstants
implicit none
real(sp),intent(in):: v(:),Dnu_D(:)
real(sp)::            v2u_1sv2sv(size(v))

if (size(v) .ne. size(Dnu_D)) stop "Size mismatch in v2u_1sv2sv"
v2u_1sv2sv = v*1e5 * nu_0 / (Dnu_D * c)

end function v2u_1sv2sv
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

function v2u_1dv2dv(v,Dnu_D)
use kinds
use PhysicalConstants
implicit none
real(dp),intent(in):: v(:),Dnu_D(:)
real(dp)::            v2u_1dv2dv(size(v))

if (size(v) .ne. size(Dnu_D)) stop "Size mismatch in v2u_1dv2dv"
v2u_1dv2dv = v*1d5 * nu_0 / (Dnu_D * c)

end function v2u_1dv2dv
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

function v2u_1s2sv(v,Dnu_D)
use kinds
use PhysicalConstants
implicit none
real(sp),intent(in):: v,Dnu_D(:)
real(sp)::            v2u_1s2sv(size(Dnu_D))

v2u_1s2sv = v*1e5 * nu_0 / (Dnu_D * c)

end function v2u_1s2sv
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

function v2u_1d2dv(v,Dnu_D)
use kinds
use PhysicalConstants
implicit none
real(dp),intent(in):: v,Dnu_D(:)
real(dp)::            v2u_1d2dv(size(Dnu_D))

v2u_1d2dv = v*1d5 * nu_0 / (Dnu_D * c)

end function v2u_1d2dv
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

function v2u_1sv2s(v,Dnu_D)
use kinds
use PhysicalConstants
implicit none
real(sp),intent(in):: v(:),Dnu_D
real(sp)::            v2u_1sv2s(size(v))

v2u_1sv2s = v*1e5 * nu_0 / (Dnu_D * c)

end function v2u_1sv2s
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

function v2u_1dv2d(v,Dnu_D)
use kinds
use PhysicalConstants
implicit none
real(dp),intent(in):: v(:),Dnu_D
real(dp)::            v2u_1dv2d(size(v))

v2u_1dv2d = v*1d5 * nu_0 / (Dnu_D * c)

end function v2u_1dv2d
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

function v2u_1s2d(v,Dnu_D)
use kinds
use PhysicalConstants
implicit none
real(sp),intent(in):: v
real(dp),intent(in):: Dnu_D
real(dp)::            v2u_1s2d

v2u_1s2d = v*1d5 * nu_0 / (Dnu_D * c)

end function v2u_1s2d
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

function v2u_1d2s(v,Dnu_D)
use kinds
use PhysicalConstants
implicit none
real(dp),intent(in):: v
real(sp),intent(in):: Dnu_D
real(dp)::            v2u_1d2s

v2u_1d2s = v*1d5 * nu_0 / (Dnu_D * c)

end function v2u_1d2s
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

function v2u_1sv2d(v,Dnu_D)
use kinds
use PhysicalConstants
implicit none
real(sp),intent(in):: v(:)
real(dp),intent(in):: Dnu_D
real(dp)::            v2u_1sv2d(size(v))

v2u_1sv2d = v*1d5 * nu_0 / (Dnu_D * c)

end function v2u_1sv2d
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

function v2u_1dv2s(v,Dnu_D)
use kinds
use PhysicalConstants
implicit none
real(dp),intent(in):: v(:)
real(sp),intent(in):: Dnu_D
real(dp)::            v2u_1dv2s(size(v))

v2u_1dv2s = v*1d5 * nu_0 / (Dnu_D * c)

end function v2u_1dv2s
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

function v2u_1s2dv(v,Dnu_D)
use kinds
use PhysicalConstants
implicit none
real(sp),intent(in):: v
real(dp),intent(in):: Dnu_D(:)
real(dp)::            v2u_1s2dv(size(Dnu_D))

v2u_1s2dv = v*1d5 * nu_0 / (Dnu_D * c)

end function v2u_1s2dv
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

function v2u_1d2sv(v,Dnu_D)
use kinds
use PhysicalConstants
implicit none
real(dp),intent(in):: v
real(sp),intent(in):: Dnu_D(:)
real(dp)::            v2u_1d2sv(size(Dnu_D))

v2u_1d2sv = v*1d5 * nu_0 / (Dnu_D * c)

end function v2u_1d2sv
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

function v2u_1sv2dv(v,Dnu_D)
use kinds
use PhysicalConstants
implicit none
real(sp),intent(in):: v(:)
real(dp),intent(in):: Dnu_D(:)
real(dp)::            v2u_1sv2dv(size(v))

if (size(v) .ne. size(Dnu_D)) stop "Size mismatch in v2u_1sv2dv"
v2u_1sv2dv = v*1d5 * nu_0 / (Dnu_D * c)

end function v2u_1sv2dv
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

function v2u_1dv2sv(v,Dnu_D)
use kinds
use PhysicalConstants
implicit none
real(dp),intent(in):: v(:)
real(sp),intent(in):: Dnu_D(:)
real(dp)::            v2u_1dv2sv(size(v))

if (size(v) .ne. size(Dnu_D)) stop "Size mismatch in v2u_1dv2sv"
v2u_1dv2sv = v*1d5 * nu_0 / (Dnu_D * c)

end function v2u_1dv2sv
!------------------------------------------------------------------------------

function XsecD(Dnu_D,x)

!Gnedin dust cross section

use kinds
use PhysicalConstants
use GlobalData

implicit none
real(dp):: XsecD, Dnu_D, x, ss
real(dp), parameter, dimension(7):: &
  lambdaSMC = (/0.042, 0.08, 0.22, 9.7, 18., 25., 0.067/),    &
       aSMC = (/185., 27., 0.005, 0.010, 0.012, 0.030, 10./), &
       bSMC = (/90., 15.5, -1.95, -1.95, -1.8, 0., 1.9/),     &
       pSMC = (/2., 4., 2., 2., 2., 2., 4./),                 &
       qSMC = (/2., 4., 2., 2., 2., 2., 15./),                &
  lambdaLMC = (/0.046, 0.08, 0.22, 9.7, 18., 25., 0.067/),    &
       aLMC = (/90., 19., 0.023, 0.005, 0.006, 0.02, 10./),   &
       bLMC = (/90., 21., -1.95, -1.95, -1.8, 0., 1.9/),      &
       pLMC = (/2., 4.5, 2., 2., 2., 2., 4./),                &
       qLMC = (/2., 4.5, 2., 2., 2., 2., 15./)

real(dp), parameter::    sigma_0SMC = 1d-22
real(dp), parameter::    sigma_0LMC = 3d-22
real(dp), dimension(7):: ll,aa,bb,pp,qq,xx,fff

! if (DustType .eq. 'SMC') then
!   ll = lambdaSMC
!   aa = aSMC
!   bb = bSMC
!   pp = pSMC
!   qq = qSMC
!   ss = sigma_0SMC
! else
!   ll = lambdaLMC
!   aa = aLMC
!   bb = bLMC
!   pp = pLMC
!   qq = qLMC
!   ss = sigma_0LMC
! endif

! xx    = (c / (x*Dnu_D + nu_0) * 1d4) / ll       !1d4 => wavelength is in micron
! fff   = aa / (xx**pp + xx**(-qq) + bb)
! XsecD = ss * sum(fff)

if (trim(DustType) .eq. 'SMC') then
  XsecD = (.395083 + 1.7164865d-16 * Dnu_D * x) * 1d-21
elseif (trim(DustType) .eq. 'LMC') then
  XsecD = (.723314 + 4.2156414d-16 * Dnu_D * x) * 1d-21
else
  XsecD = userXsecD
endif

end function XsecD
!------------------------------------------------------------------------------
