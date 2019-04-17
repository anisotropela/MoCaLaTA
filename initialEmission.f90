! From "RawData", calculate initial positions of photons, saving as CellData.
! 
!===========  Modules  ========================================================

module DataArrays

implicit none

integer(kind=1), dimension(:), pointer:: LevelString
real*4,    dimension(:), pointer:: n_HIString, TString, rhoString, n_HIIString, &
                                   n_HString, x_HIString, MetalString, &
                                   XString, YString, ZString,         &
                                   V_xString, V_yString, V_zString
real*4,    dimension(:), pointer:: logLya,logFUV
real*8,    dimension(:), pointer:: L_LyaString,L_FUVString
real*4,    dimension(:), pointer:: WeightStringL,WeightStringF

! real*4,    dimension(:), pointer::                  &
!                                  ! HeliumString,    &
!                                    CarbonString,    &
!                                    NitrogenString,  &
!                                    OxygenString,    &
!                                    MagnesiumString, &
!                                    SiliconString,   &
!                                    SulfurString,    &
!                                  ! CalciumString,   &
!                                    IronString

real*4,    dimension(:), pointer:: Dnu_DString
real*4,    dimension(:), pointer:: U_xString, U_yString, U_zString

end module DataArrays

!------------------------------------------------------------------------------

module PhysicalConstants
real*8,parameter:: pi      = 3.14159265359,   & !pi
                   sqpi    = 1.7724538509055, & !SQRT(pi)
                   kpc     = 3.08567758131d21,& !Kiloparsec; cm
                   Dnu_L   = 9.936d7,         & !Lya natural line width; Hz
                   nu_0    = 2.46607d15,      & !Lya line-center frequency; Hz
                   f_12    = 0.4162,          & !Oscillator strength
                   e       = 4.803206d-10,    & !Electron charge; esu
                   m_e     = 9.1093897d-28,   & !Electron mass; g
                   c       = 2.99792458d10,   & !Speed of light; cm/s
                   k_B     = 1.380658d-16,    & !Boltzmann's constant; erg/K
                   m_H     = 1.673534d-24,    & !Hydrogen mass; g
                   m_au    = 1.660538782d-24, & !Atomic mass unit
                   M_sun   = 1.989d33,        & !Mass of Sun; g
                   h_Pl    = 6.6260755d-27,   & !Planck's constant; erg s
                   E_Lya   = h_Pl * nu_0,     & !Lya energy; erg
                   H_0     = 71.d5,           & !Hubble constant; cm/s/Mpc
                   Omega_M = 0.3,             & !Matter density parameter
                   Omega_L = 0.7                !Dark energy density parameter

real*8,parameter:: n_He_H_sun = 8.511d-02,&     ! Number
                   n_C_H_sun  = 2.455d-04,&     ! densities
                   n_N_H_sun  = 6.026d-05,&     ! relative
                   n_O_H_sun  = 4.571d-04,&     ! to
                   n_Mg_H_sun = 3.388d-05,&     ! hydrogen;
                   n_Si_H_sun = 3.236d-05,&     ! Solar
                   n_S_H_sun  = 1.380d-05,&     ! values
                   n_Ca_H_sun = 2.042d-06,&     ! (i.e. = 1 for H)
                   n_Fe_H_sun = 2.818d-05       !

real*8,parameter::&! Z_He =                  &    !
                   Z_C  = 3.03e-3,         &    !
                   Z_N  = 9.15e-4,         &    ! Mass
                   Z_O  = 8.3d-3,          &    ! fractions;
                   Z_Mg = 5.15d-4,         &    ! Solar
                   Z_Si = 6.53e-4,         &    ! values
                   Z_S  = 4.55e-4,         &    !
                   Z_Ca = 5.99e-5,         &    !
                   Z_Fe = 1.2d-3                !

real*8,parameter:: m_He =  4.002602 * m_au,&    ! Atomic
                   m_C  = 12.0107   * m_au,&    ! masses
                   m_N  = 14.0067   * m_au,&    ! in
                   m_O  = 15.9994   * m_au,&    ! grams
                   m_Mg = 24.3050   * m_au,&    !
                   m_Si = 28.0855   * m_au,&    !
                   m_S  = 32.065    * m_au,&    !
                   m_Ca = 40.078    * m_au,&    !
                   m_Fe = 55.845    * m_au      !
end module PhysicalConstants

!------------------------------------------------------------------------------

module GlobalData

implicit none
integer::                           ni,nj,nk,N_cells,L_max,n_ph(2),n_phtot
real*8::                            D_box,R_box,D_obs,R_emit,D_emit
real*8::                            L_Lya       !Total Lya lum. in box
real*8::                            L_FUV_Hz    !Total FUV lum. in box per Hz
real*8::                            L_FUV       !Total FUV lum. in box
real*8::                            L_tot       !Total lum. in box (Lya + FUV)
real*8::                            Lya2tot     !Fraction of photons which are Lya
real*8::                            Q_Lya       !Total # of Lya photons (not packets), although possibly weighted
real*8::                            Q_FUV       !Total # of FUV photons (not packets), although possibly weighted
real*8::                            Q_tot       !Total # of all photons (not packets), although possibly weighted
real*8::                            wL_Lya,wL_FUV_Hz  !Total weighted luminosities
real*8::                            dx0,dy0,dz0,dxTiny
real*8,allocatable,dimension(:)::   dx,dxHalf,dy,dyHalf,dz,dzHalf
integer::                           n_m,m,Nc_beg
integer,allocatable::               Milestones(:),Retard(:)
real*4,allocatable,dimension(:,:):: LyaEmissionString,FUVEmissionString
character(200)::                    DataDir,SubDir,RawData,CellData,filein,fileout
character(len=3)::                  DustType
logical::                           writepos
real*8::                            f,rvir
real*8::                            bbL,bbF,BW(2)
integer::                           emdimL,emdimF

end module GlobalData

!------------------------------------------------------------------------------

program initialEmission

use DataArrays
use GlobalData

implicit none

writepos = .false.

call ReadInput
call ReadData
call CalcDX
call FindInitPos
!call SumMetals
call WriteData

end program initialEmission

! ==========  Subroutines  ====================================================

subroutine ReadInput

use GlobalData

implicit none
character(len=200):: dummy,X_init

  read*, dummy
  read*, dummy

! I/O
  read*, dummy
  read*, DataDir      !Mother directory
  read*, Subdir       !Subdirectory
  read*, RawData      !Name of data file
  read*, CellData     !Name of data file
  read*, dummy        !Print to std. output for every wrtph phot
  read*, dummy        !Write file for every wrtfl photon
  read*, dummy        !Label for main output files
  read*, dummy        !Write input spectrum character
  read*, dummy        !-"- for output spectrum in Doppler widths

! Observations
  read*, dummy
  read*, BW           !Wavelength interval bounds in Angstrom
  read*, dummy        !# of pixels/side in CCD
  read*, dummy        !Resolution of 1D spectrum; bins
  read*, dummy        !Resolution of 2D spectrum; bins
  read*, dummy        !(xp,xm,yp,ym,zp,zm)

! Simulation
  read*, dummy                          
  read*, dummy        !Side length of area covered by CCD
  read*, D_emit       !Side length of box, outside which lum. is set to zero
  read*, n_phtot
  read*, bbL
  read*, bbF
  read*, dummy        !'variable', 'global', or user-given value (as str
  read*, dummy        !Initial frequency: 'proper', '...
  read*, X_init       !Initial position: 'central', 'homo', 'lumdep', or
  read*, dummy        !Initial direction: 'iso' or '(nx ny nz)'
  read*, dummy        !Use "true" random numbers, given by current date 
  read*, dummy        !Fix physical parameters Interactively after Unive

if (trim(X_init) .ne. 'lumdep') then
  print*, "WARNING: Keyword X_init is set to '", trim(X_init), "', rather than 'lumdep'."
  call holdon(3.)
  print*, 'Continuing...'
  call holdon(1.)
endif

R_emit = D_emit / 2d0

filein  = trim(DataDir)//'/'// trim(subdir)//'/'// trim(RawData)
fileout = trim(DataDir)//'/'// trim(subdir)//'/'// trim(CellData)

end subroutine ReadInput

!------------------------------------------------------------------------------

subroutine ReadData

use DataArrays
use GlobalData

implicit none
real*8::  RRR,d,xx,yy,zz,pos(3),nhat(3),V(3)
integer:: i,mmm

open(14,file=trim(filein),action='read',status='old',form='unformatted')
read(14) N_cells, D_box, ni,nj,nk               !Tot. #cells, phys size in kpc,
                                                !and resolution
print*, N_cells, D_box, ni,nj,nk
R_box = D_box / 2
allocate(LevelString(N_cells))
allocate(XString(N_cells))
allocate(YString(N_cells))
allocate(ZString(N_cells))
allocate(L_LyaString(N_cells))
allocate(L_FUVString(N_cells))
allocate(logLya(N_cells))
allocate(logFUV(N_cells))
allocate(rhoString(N_cells))
allocate(n_HIString(N_cells))
allocate(n_HIIString(N_cells))
allocate(n_HString(N_cells))
allocate(x_HIString(N_cells))
allocate(TString(N_cells))
allocate(V_xString(N_cells))
allocate(V_yString(N_cells))
allocate(V_zString(N_cells))
!allocate(HeliumString(N_cells))
!allocate(CarbonString(N_cells))
!allocate(NitrogenString(N_cells))
!allocate(OxygenString(N_cells))
!allocate(MagnesiumString(N_cells))
!allocate(SiliconString(N_cells))
!allocate(SulfurString(N_cells))
!allocate(CalciumString(N_cells))
!allocate(IronString(N_cells))
allocate(MetalString(N_cells))

read(14) (LevelString(i), i=1,N_cells)! 0 - 7
read(14) (XString(i),     i=1,N_cells)! [0,D_box] in kpc
read(14) (YString(i),     i=1,N_cells)! [0,D_box] in kpc
read(14) (ZString(i),     i=1,N_cells)! [0,D_box] in kpc

read(14) (logLya(i),      i=1,N_cells)
read(14) (logFUV(i),      i=1,N_cells)
do i=1,N_cells
  if (logLya(i) .le. 0d0) logLya(i)  = -35.
  if (logFUV(i) .le. 0d0) logFUV(i)  = -35.
  L_LyaString(i) = 10d0**logLya(i)
  L_FUVString(i) = 10d0**logFUV(i)
enddo

 do i = 1,N_cells
   xx = XString(i)-R_box
   yy = YString(i)-R_box
   zz = ZString(i)-R_box
   if (any(abs((/xx,yy,zz/)) .gt. R_emit)) then
     L_LyaString(i) = 0d0
     L_FUVString(i) = 0d0
   endif
 enddo
 L_Lya    = sum(L_LyaString)
 L_FUV_Hz = sum(L_FUVString)
 write(*,'(a24,f5.0,a7,es16.9)') 'L_Lya,tot inside D_emit (', D_emit, ' kpc): ', L_Lya
 write(*,'(a24,f5.0,a7,es16.9)') 'L_FUV_Hz,tot inside D_emit (', D_emit, ' kpc): ', L_FUV_Hz

read(14) (rhoString(i),       i=1,N_cells)
read(14) (n_HString(i),       i=1,N_cells)
read(14) (x_HIString(i),      i=1,N_cells)
read(14) (TString(i),         i=1,N_cells)
read(14) (V_xString(i),       i=1,N_cells)! km/s
read(14) (V_yString(i),       i=1,N_cells)! km/s
read(14) (V_zString(i),       i=1,N_cells)! km/s
!read(14)(HeliumString(i),    i=1,N_cells)
!read(14) (CarbonString(i),    i=1,N_cells)
!read(14) (NitrogenString(i),  i=1,N_cells)
!read(14) (OxygenString(i),    i=1,N_cells)
!read(14) (MagnesiumString(i), i=1,N_cells)
!read(14) (SiliconString(i),   i=1,N_cells)
!read(14) (SulfurString(i),    i=1,N_cells)
!read(14)(CalciumString(i),   i=1,N_cells)
!read(14) (IronString(i),      i=1,N_cells)
read(14) (MetalString(i), i=1,N_cells)
close(14)

!do i = 1,N_cells
!  xx = XString(i)-R_box
!  yy = YString(i)-R_box
!  zz = ZString(i)-R_box
!  d  = sqrt(sum((/xx,yy,zz/)**2))
!  V  = (/V_xString(i), V_zString(i), V_zString(i)/)
!  write(1,*) d,sqrt(sum(V**2))
!enddo

do i = 1,N_cells
  n_HIString(i)  = n_HString(i) * x_HIString(i)
  n_HIIString(i) = n_HString(i) * (1. - x_HIString(i))
enddo

L_max = maxval(LevelString)

end subroutine ReadData

!------------------------------------------------------------------------------

subroutine WriteData

use DataArrays
use GlobalData
use PhysicalConstants

implicit none
integer:: i,n,L,n_recL,n_recF,ifracL,ifracF,mmm
real*8::  R, xx,yy,zz,dd,vx,vy,vz,vv,recmaxMb

recmaxMb = 10.                                  !Maximum record size in Mb
n_recL   = size(LyaEmissionString) / nint(1e6 * recmaxMb / kind(LyaEmissionString)) + 1
n_recF   = size(FUVEmissionString) / nint(1e6 * recmaxMb / kind(FUVEmissionString)) + 1
ifracL   = n_ph(1) / n_recL
ifracF   = n_ph(2) / n_recF

print*, 'n_recL', n_recL
print*, 'n_recF', n_recF
print*, 'ifracL', ifracL
print*, 'ifracF', ifracF

open(14,file=trim(fileout),action='write',status='replace',form='unformatted')

write(14) N_cells, D_box, ni,nj,nk, L_tot, BW, n_ph, emdimL, emdimF, n_recL, n_recF

do i = 1,n_recL
  if (i .lt. n_recL) then
    write(14) (LyaEmissionString(n,1:emdimL), n = (i-1)*ifracL+1, i*ifracL)
  else
    write(14) (LyaEmissionString(n,1:emdimL), n = (i-1)*ifracL+1, n_ph(1))
  endif
enddo

do i = 1,n_recF
  if (i .lt. n_recF) then
    write(14) (FUVEmissionString(n,1:emdimF), n = (i-1)*ifracF+1, i*ifracF)
  else
    write(14) (FUVEmissionString(n,1:emdimF), n = (i-1)*ifracF+1, n_ph(2))
  endif
enddo

write(14) (LevelString(i), i = 1,N_cells)          ! 0 - 7
write(14) (n_HIString(i),  i = 1,N_cells)
write(14) (TString(i),     i = 1,N_cells)
write(14) (V_xString(i),   i = 1,N_cells) !cgs
write(14) (V_yString(i),   i = 1,N_cells) !cgs
write(14) (V_zString(i),   i = 1,N_cells) !cgs
write(14) (MetalString(i), i = 1,N_cells)
write(14) (n_HIIString(i), i = 1,N_cells)
close(14)

! mmm = sum(maxloc(n_HString))
! do i = mmm-20,mmm+20
!   write(1,*) LevelString(i),n_HIString(i),TString(i)
! enddo
! stop

! do i = 1,N_cells
!   call random_number(R)
!   if (r.lt..005) then
!     xx = XString(i) - R_box
!     yy = YString(i) - R_box
!     zz = ZString(i) - R_box
!     dd = sqrt(sum((/xx**2, yy**2, zz**2/)))
! 
!     vx = V_xString(i)
!     vy = V_yString(i)
!     vz = V_zString(i)
!     vv = sqrt(sum((/vx**2, vy**2, vz**2/)))
! 
!     write(1,*) dd, vv
!   endif
! enddo

!rint*,'nmaxloc',maxloc(n_HIString)
!rint*,'Tmaxloc',maxloc(TString)

print*,'nmax',maxval(n_HIString)
print*,'Tmax',maxval(TString)
print*, L_max

end subroutine WriteData

!------------------------------------------------------------------------------

subroutine CalcDX

use DataArrays
use GlobalData

implicit none
integer:: L

allocate(dx(0:L_max))
allocate(dy(0:L_max))
allocate(dz(0:L_max))
allocate(dxHalf(0:L_max))
allocate(dyHalf(0:L_max))
allocate(dzHalf(0:L_max))

dx0 = D_box / ni
dy0 = D_box / nj
dz0 = D_box / nk

do L = 0,L_max
  dx(L) = dx0 / 2d0**L
  dy(L) = dy0 / 2d0**L
  dz(L) = dz0 / 2d0**L
enddo

dxHalf = dx / 2d0
dyHalf = dy / 2d0
dzHalf = dz / 2d0

dxTiny = 0.1d0 * dx(L_max)
  
end subroutine CalcDX

!------------------------------------------------------------------------------

subroutine FindInitPos

use DataArrays
use GlobalData
use PhysicalConstants

implicit none
integer:: i,l,n,level
real*8::  L_cum,F_cum,R,RRR(3),R1,R2,R3,X_pos(3)
real*8:: t1,t2,few
real*8:: logL_max,logF_max,logL_Lmax,logF_Fmax,wL,wF,nu1,nu2

emdimL = merge(3,4, bbL.gt.1e10) 
emdimF = merge(3,4, bbF.gt.1e10) 

if (emdimL .eq. 4) allocate(WeightStringL(N_cells))
if (emdimF .eq. 4) allocate(WeightStringF(N_cells))

nu1 = c / BW(2) * 1d8                      !Hz
nu2 = c / BW(1) * 1d8                      !Hz

L_Lya    = sum(L_LyaString)
L_FUV_Hz = sum(L_FUVString)
L_cum    = 0d0
F_cum    = 0d0

print*, 'L_Lya   ',L_Lya
print*, 'L_FUV_Hz',L_FUV_Hz
!print*, 'Cumulating luminosity...'
logL_max = log10(maxval(L_LyaString))
logF_max = log10(maxval(L_FUVString))

if (emdimL .eq. 4) then
  print*, 'Weighting Lya photons...'
  print*, '      -----------------------------------------------'
  print*, '     |                                               |'
  print*, '     | CHECK HOW n_ph DEPENDS ON DIFFERENT WEIGHTING |'
  print*, '     |                                               |'
  print*, '      -----------------------------------------------'
  stop
  do i = 1,N_cells
    logL_Lmax        = log10(L_LyaString(i)) - logL_max
    wL               = 10**(logL_Lmax / bbL)
    WeightStringL(i) = wL
    L_LyaString(i)   = L_LyaString(i) / wL
  enddo
endif

if (emdimF .eq. 4) then
  print*, 'Weighting FUV photons...'
  do i = 1,N_cells
    logF_Fmax        = log10(L_FUVString(i)) - logF_max
    wF               = 10**(logF_Fmax / bbF)
    WeightStringF(i) = wF
    L_FUVString(i)   = L_FUVString(i) / wF
  enddo
endif

if (emdimL .eq. 4) wL_Lya    = sum(L_LyaString)  !Artificial total Lya luminosity
if (emdimF .eq. 4) wL_FUV_Hz = sum(L_FUVString)  !Artificial total FUV luminosity

Q_Lya   = L_Lya / E_Lya                         !ph/s
Q_FUV   = L_FUV_Hz/h_Pl * (log(nu2)-log(nu1))   !ph/s
Q_tot   = Q_Lya + Q_FUV                         !ph/s
n_ph(1) = nint(Q_Lya/Q_tot * n_phtot)           !Total # of Lya photon packets
n_ph(2) = n_phtot - n_ph(1)                     !Total # of FUV photon packets
Lya2tot = Q_Lya / Q_tot                         !Probability for an emitted photon packet to be Lya
L_FUV   = L_FUV_Hz * (nu2-nu1)
L_tot   = L_Lya + L_FUV
  print*, 'n_ph    ', n_ph
  print*, 'Q_Lya   ', Q_Lya
  print*, 'Q_FUV   ', Q_FUV
  print*, 'Q_FUV ~ ', L_FUV_Hz/E_Lya * (nu2-nu1)
  print*, 'Q_tot   ', Q_tot
  print*, 'L_FUV   ', L_FUV
  print*, 'L_FUV ~ ', Q_FUV * E_Lya
! stop

allocate(LyaEmissionString(n_ph(1),emdimL))
allocate(FUVEmissionString(n_ph(2),emdimF))

n_m = 10000 !Optimally, I guess this should rather be a function of n_ph.
allocate(Milestones(n_m))
allocate(Retard(n_m))

!     --------------------  Find Lya positions  --------------------
  Milestones(1) = 1
  m             = 2
  do i = 1,N_cells
    L_cum          = L_cum + L_LyaString(i) 
    L_LyaString(i) = L_cum / L_Lya
    if (L_LyaString(i).gt.(m-1)/real(n_m) .and. m.le.n_m) then
      Milestones(m) = i
  !   if (Milestones(m) - Milestones(m-1) .eq. 1) Milestones(m) = i - 1
      m = m + 1
    endif
  enddo

  Retard = 0
  do m = 2,n_m
    if (Milestones(m) - Milestones(m-1) .eq. 1) Retard(m) = 1
  enddo
  do m = 1,n_m
    if (Retard(m) .eq. 1) Milestones(m) = Milestones(m-1)
  enddo

  call cpu_time(t1)
  l = 1
  do
    call random_number(R)
    m      = ceiling(R * n_m)
    Nc_beg = Milestones(m)

    do n = Nc_beg,N_cells
      i = n
      if (L_LyaString(i) .gt. R) exit
    enddo
    
    X_pos = (/XString(i), YString(i), ZString(i)/)

    level = LevelString(i)
    call random_number(RRR)
    X_pos = X_pos + (RRR-.5)*(/dx(level),dy(level),dz(level)/)
  ! if (any(abs(X_pos-R_box) .gt. R_emit)) cycle

    LyaEmissionString(l,1:3) = real(X_pos(1:3))
    if (writepos) then
      call random_number(R)
      if (R.lt..05) write(1,*) X_pos - R_box
    endif
    if (emdimL .eq. 4) LyaEmissionString(l,4) = real(WeightStringL(i))

    if (n_ph(1).gt.0) then
      if (mod(l,n_ph(1)/10) .eq. 0) print*, 'Found position for Lya photon',l,'of',n_ph(1)
    endif
    l = l + 1
    if (l .gt. n_ph(1)) exit
  enddo
  call cpu_time(t2)
  print*,'s/ph.:',(t2-t1)/n_ph(1)

!     --------------------  Find FUV positions  --------------------
  Milestones(1) = 1
  m             = 2
  do i = 1,N_cells
    F_cum          = F_cum + L_FUVString(i) 
    L_FUVString(i) = F_cum / L_FUV_Hz
    if (L_FUVString(i).gt.(m-1)/real(n_m) .and. m.le.n_m) then
      Milestones(m) = i
  !   if (Milestones(m) - Milestones(m-1) .eq. 1) Milestones(m) = i - 1
      m = m + 1
    endif
  enddo

  Retard = 0
  do m = 2,n_m
    if (Milestones(m) - Milestones(m-1) .eq. 1) Retard(m) = 1
  enddo
  do m = 1,n_m
    if (Retard(m) .eq. 1) Milestones(m) = Milestones(m-1)
  enddo

  call cpu_time(t1)
  l = 1
  do
    call random_number(R)
    m      = ceiling(R * n_m)
    Nc_beg = Milestones(m)

    do n = Nc_beg,N_cells
      i = n
      if (L_FUVString(i) .gt. R) exit
    enddo
    
    X_pos = (/XString(i), YString(i), ZString(i)/)

    call random_number(RRR)
    X_pos = X_pos + (RRR-.5) * dx(LevelString(i))
  ! if (any(abs(X_pos-R_box) .gt. R_emit)) cycle

    FUVEmissionString(l,1:3) = real(X_pos(1:3))
    if (writepos) then
      call random_number(R)
      if (R.lt..05) write(2,*) X_pos - R_box
    endif
    if (emdimF .eq. 4) FUVEmissionString(l,4) = real(WeightStringF(i))

    if (n_ph(2).gt.0) then
      if (mod(l,n_ph(2)/10) .eq. 0) print*, 'Found position for FUV photon',l,'of',n_ph(2)
    endif
    l = l + 1
    if (l .gt. n_ph(2)) exit
  enddo
  call cpu_time(t2)
  print*,'s/ph.:',(t2-t1)/n_ph(2)

end subroutine FindInitPos

!------------------------------------------------------------------------------

! subroutine SumMetals
! 
! use DataArrays
! use GlobalData
! use PhysicalConstants
! 
! implicit none
! integer:: j,L
! real*8::  Z,Z_sol,R
! 
! print*, 'Summing metals...'
! 
! allocate(MetalString(N_cells))
! 
! Z_sol =  Z_C + Z_N + Z_O + Z_Mg + Z_Si + Z_S + Z_Fe !+ Z_Ca
! 
! do j = 1,N_cells
!   Z = CarbonString(j)    * Z_C   &
!     + NitrogenString(j)  * Z_N   &
!     + OxygenString(j)    * Z_O   &
!     + MagnesiumString(j) * Z_Mg  &
!     + SiliconString(j)   * Z_Si  &
!     + SulfurString(j)    * Z_S   &
! !   + CalciumString(j)   * Z_Ca  &
!     + IronString(j)      * Z_Fe   
! 
!   MetalString(j) = Z / Z_sol
! ! call random_number(R)
! ! if (R.lt..005) write(1,*) Z
! enddo
! 
! end subroutine SumMetals
! 
!------------------------------------------------------------------------------

subroutine holdon(t)

!Pause execution for t seconds. Use only for short intervals, as the processor
!will keep spinning.
!NB: Argument is single precision!

implicit none
real*4,intent(in):: t
real*4::            t1,t2

call cpu_time(t1)

do
  call cpu_time(t2)
  if (t2-t1 .gt. t) exit
enddo

end subroutine holdon

!------------------------------------------------------------------------------
