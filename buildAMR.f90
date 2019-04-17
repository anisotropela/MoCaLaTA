
module localDefinitions

implicit none
integer, parameter :: SingleKind = kind(1.0)
integer, parameter :: RealKind = kind(1.0d0)
real(kind=RealKind), parameter :: pc      = 3.08568025e18
real(kind=RealKind), parameter :: kpc     = 1.e3*pc
real(kind=RealKind), parameter :: Mpc     = 1.e6*pc
real(kind=RealKind), parameter :: mp      = 1.6726231e-24
real(kind=RealKind), parameter :: mn      = 1.67492728e-24
real(kind=RealKind), parameter :: mh      = mp
!real(kind=RealKind), parameter :: mhe     = 2.*(mp+mn)
real(kind=RealKind), parameter :: psi     = 0.76
real(kind=RealKind), parameter :: G       = 4.3d-6
real(kind=RealKind), parameter :: M_sun   = 1.99d33
real(kind=RealKind), parameter :: Lya2SFR = 9.08d-43, FUV2SFR = 1.4d-28
real(kind=RealKind), parameter :: Lya2FUV = Lya2SFR/FUV2SFR, logLya2FUV = log10(Lya2FUV)

real(kind=RealKind):: m_SPH, L_tot

integer :: icosmic, ncosmic
real:: R

type :: zoneType
   real*4::  rho, tgas, HI, V(3)
 ! real*4::  HeI, HeII
   real*4::  Hydrogen, x
 ! real*4::  Helium
 ! real*4::  Carbon
 ! real*4::  Nitrogen
 ! real*4::  Oxygen
 ! real*4::  Magnesium
 ! real*4::  Silicon
 ! real*4::  Sulfur
 ! real*4::  Calcium
 ! real*4::  Iron
   real*4::  Z
   real*4 :: logL_Lya_HII,logL_Lya_cool
   logical(1) :: refined
   integer(1) :: level ! maximum 127 levels of cell refinement
   type(zoneType), pointer :: parent
   type(zoneType), dimension(:,:,:), pointer :: cell
end type zoneType

type(zoneType), target :: baseGrid

                                !cellArrayLevel is 4 byte, but what is WRITTEN
                                !is 1 byte: int(cellArrayLevel(i),1)
integer*4, dimension(:), pointer :: cellArrayLevel,cellArrayLevelORIG
real*4, dimension(:), pointer :: cellArrayXpos,  cellArrayYpos, cellArrayZpos, &
                                 cellArrayXvel,  cellArrayYvel, cellArrayZvel, &
                                 cellArrayHydrogen,    cellArrayx,             &
                                 cellArrayHI,                                  &
                               !                 cellArrayHeI,  cellArrayHeII, &
                               ! cellArrayHelium,                              &
                               ! cellArrayCarbon,                              &
                               ! cellArrayNitrogen,                            &
                               ! cellArrayOxygen,                              &
                               ! cellArrayMagnesium,                           &
                               ! cellArraySilicon,                             &
                               ! cellArraySulfur,                              &
                               ! cellArrayCalcium,                             &
                               ! cellArrayIron,                                &
                                 cellArrayZ,                                   &
                                 DammitAlex,                                &
                                 cellArrayTemp, cellArrayDensity
real*4, dimension(:), pointer :: cellArraylogL_Lya_HII,cellArraylogL_Lya_cool

type :: readLevelType
   integer :: ncell
   real*4, dimension(:,:), pointer :: pos, vel, abun, LLY
   real*4, dimension(:), pointer :: lT, lnH, lx
end type readLevelType

real*8:: R_box,dx0,sumHII,sumcool,VV,dd,V_tot,V_HI,V_HII,X_Gal(3),r_0

end module localDefinitions

!*****************************************************************************

program BuildAMR
use localDefinitions

implicit none
integer :: cc, i, j, k, icell, nx, ny, nz, i0, j0, k0, lmax, nlevels, level, itmp, nmetals
real(kind=RealKind) :: x0, y0, z0, xa, xb, ya, yb, za, zb, tmp1, tmp2, &
     xnew, ynew, znew, xpos, ypos, zpos, LLY
type(readLevelType), dimension(:), pointer :: readLevel
real:: loT
character(200):: DataDir,SubDir,RawData,filein,fileout,dummy

interface
  recursive subroutine placeCellProject(ParentCell,level,x0,y0,z0, &
                                        logtgas,lognh,logxneu,     &
                                      ! pCPHelium,                 &
                                      ! pCPCarbon,                 &
                                      ! pCPNitrogen,               &
                                      ! pCPOxygen,                 &
                                      ! pCPMagnesium,              &
                                      ! pCPSilicon,                &
                                      ! pCPSulfur,                 &
                                      ! pCPCalcium,                &
                                      ! pCPIron,                   &
                                        pCPZ,                      &
                                        vx,vy,vz,logHII,logcool)
    use localDefinitions
    implicit none
    integer, intent(in) ::             level
    real(kind=RealKind), intent(in) :: x0, y0, z0
    real*4, intent(in) ::              logtgas, lognh, logxneu, vx,vy,vz,logHII,logcool
  ! real*4, intent(in) ::              pCPHelium
  ! real*4, intent(in) ::              pCPCarbon
  ! real*4, intent(in) ::              pCPNitrogen
  ! real*4, intent(in) ::              pCPOxygen
  ! real*4, intent(in) ::              pCPMagnesium
  ! real*4, intent(in) ::              pCPSilicon
  ! real*4, intent(in) ::              pCPSulfur
  ! real*4, intent(in) ::              pCPCalcium
  ! real*4, intent(in) ::              pCPIron
    real*4, intent(in) ::              pCPZ
    type(zoneType), target ::          ParentCell
  end subroutine placeCellProject
end interface

interface
  recursive subroutine countCells(currentCell)
  use localDefinitions
  implicit none
  type(zoneType), target :: currentCell
  end subroutine countCells
end interface

read*, dummy        !'# Value...'
read*, dummy        !Model
read*, dummy        !I/O
read*, DataDir      !Mother directory
read*, Subdir       !Subdirectory
read*, RawData      !Name of data file

filein  = trim(DataDir)//'/' // trim(Subdir)// '/Jespers.bin'
fileout = trim(DataDir)//'/' // trim(Subdir)//'/' // trim(RawData)

! Read Jesper's dataset in binary format

!print*, "in: '"// trim(filein) //"'"
open(14,file=trim(filein), status='old',form='unformatted')

read(14) nlevels
print*, 'nlevels =', nlevels

allocate(readLevel(nlevels))

nmetals = 1 ! used to be 9

do level = 1, nlevels
  read(14) readLevel(level)%ncell
  write(*,*) 'level =', level, readLevel(level)%ncell

  allocate(readLevel(level)%pos(readLevel(level)%ncell,3))
  allocate(readLevel(level)%lT(readLevel(level)%ncell))
  allocate(readLevel(level)%lnH(readLevel(level)%ncell))
  allocate(readLevel(level)%lx(readLevel(level)%ncell))
  allocate(readLevel(level)%abun(readLevel(level)%ncell,nmetals))
  allocate(readLevel(level)%vel(readLevel(level)%ncell,3))
  allocate(readLevel(level)%LLY(readLevel(level)%ncell,2))

  read(14) (readLevel(level)%pos(icell,1),icell=1,readLevel(level)%ncell)
  read(14) (readLevel(level)%pos(icell,2),icell=1,readLevel(level)%ncell)
  read(14) (readLevel(level)%pos(icell,3),icell=1,readLevel(level)%ncell)
  read(14) (readLevel(level)%lT(icell),icell=1,readLevel(level)%ncell)
  read(14) (readLevel(level)%lnH(icell),icell=1,readLevel(level)%ncell)
  read(14) (readLevel(level)%lx(icell),icell=1,readLevel(level)%ncell)
! read(14) (readLevel(level)%abun(icell,1),icell=1,readLevel(level)%ncell) !He
! read(14) (readLevel(level)%abun(icell,2),icell=1,readLevel(level)%ncell) !C
! read(14) (readLevel(level)%abun(icell,3),icell=1,readLevel(level)%ncell) !N
! read(14) (readLevel(level)%abun(icell,4),icell=1,readLevel(level)%ncell) !O
! read(14) (readLevel(level)%abun(icell,5),icell=1,readLevel(level)%ncell) !Mg
! read(14) (readLevel(level)%abun(icell,6),icell=1,readLevel(level)%ncell) !Si
! read(14) (readLevel(level)%abun(icell,7),icell=1,readLevel(level)%ncell) !S
! read(14) (readLevel(level)%abun(icell,8),icell=1,readLevel(level)%ncell) !Ca
! read(14) (readLevel(level)%abun(icell,9),icell=1,readLevel(level)%ncell) !Fe
  read(14) (readLevel(level)%abun(icell,1),icell=1,readLevel(level)%ncell) !Z
  read(14) (readLevel(level)%vel(icell,1),icell=1,readLevel(level)%ncell)
  read(14) (readLevel(level)%vel(icell,2),icell=1,readLevel(level)%ncell)
  read(14) (readLevel(level)%vel(icell,3),icell=1,readLevel(level)%ncell)
  read(14) (readLevel(level)%LLY(icell,1:2),icell=1,readLevel(level)%ncell)
! read(14) (readLevel(level)%LLY(icell,2),icell=1,readLevel(level)%ncell)
enddo

close(14)

! Set up a fully nested 3D grid

i    = 0
itmp = 0
do while (itmp.lt.readLevel(1)%ncell)
  i    = i + 1
  itmp = i**3
enddo
if (itmp.ne.readLevel(1)%ncell) then
  write(*,*) 'base grid needs to be of size n^3', itmp, readLevel(1)%ncell
  stop
endif

nx = i
ny = i
nz = i
write(*,*) 'grid:', nx,'^3'

level = 1
xa =   1d10
xb = - 1d10
ya =   1d10
yb = - 1d10
za =   1d10
zb = - 1d10
do icell = 1, readLevel(level)%ncell
  xa = min(xa,dble(readLevel(level)%pos(icell,1)))
  xb = max(xb,dble(readLevel(level)%pos(icell,1)))
  ya = min(ya,dble(readLevel(level)%pos(icell,2)))
  yb = max(yb,dble(readLevel(level)%pos(icell,2)))
  za = min(za,dble(readLevel(level)%pos(icell,3)))
  zb = max(zb,dble(readLevel(level)%pos(icell,3)))
enddo

tmp1 = 0.5*(xa+xb)
tmp2 = 0.5*(xb-xa)*float(nx)/float(nx-1)
xa   = tmp1 - tmp2
xb   = tmp1 + tmp2

tmp1 = 0.5*(ya+yb)
tmp2 = 0.5*(yb-ya)*float(ny)/float(ny-1)
ya   = tmp1 - tmp2
yb   = tmp1 + tmp2

tmp1 = 0.5*(za+zb)
tmp2 = 0.5*(zb-za)*float(nz)/float(nz-1)
za   = tmp1 - tmp2
zb   = tmp1 + tmp2

write(*,*) 'edges of computational grid in kpc:'
write(*,*) xa, xb                                      !
write(*,*) ya, yb                                      !  [-R_box, +R_box]
write(*,*) za, zb                                      !

R_box = xb
dx0   = 2d0 * R_box / nx

do level = 1, nlevels
  do icell = 1, readLevel(level)%ncell
    readLevel(level)%pos(icell,1) = (readLevel(level)%pos(icell,1)-xa)/(xb-xa) !
    readLevel(level)%pos(icell,2) = (readLevel(level)%pos(icell,2)-ya)/(yb-ya) !  0.000  -  1.000
    readLevel(level)%pos(icell,3) = (readLevel(level)%pos(icell,3)-za)/(zb-za) ! (act.  0.00390640 -  0.996094)
  enddo
enddo

print*, 'setting up a fully nested 3D grid'

allocate(baseGrid%cell(nx,ny,nz))

do i = 1, nx
   do j = 1, ny
      do k = 1, nz
         baseGrid%cell(i,j,k)%refined   = .false.
         baseGrid%cell(i,j,k)%tgas      = 0e0
         baseGrid%cell(i,j,k)%rho       = 0e0
         baseGrid%cell(i,j,k)%Hydrogen  = 0e0
         baseGrid%cell(i,j,k)%x         = 0e0
         baseGrid%cell(i,j,k)%HI        = 0e0
       ! baseGrid%cell(i,j,k)%HeI       = 0e0
       ! baseGrid%cell(i,j,k)%HeII      = 0e0
       ! baseGrid%cell(i,j,k)%Helium    = 0e0
       ! baseGrid%cell(i,j,k)%Carbon    = 0e0
       ! baseGrid%cell(i,j,k)%Nitrogen  = 0e0
       ! baseGrid%cell(i,j,k)%Oxygen    = 0e0
       ! baseGrid%cell(i,j,k)%Magnesium = 0e0
       ! baseGrid%cell(i,j,k)%Silicon   = 0e0
       ! baseGrid%cell(i,j,k)%Sulfur    = 0e0
       ! baseGrid%cell(i,j,k)%Calcium   = 0e0
       ! baseGrid%cell(i,j,k)%Iron      = 0e0
         baseGrid%cell(i,j,k)%Z         = 0e0
         baseGrid%cell(i,j,k)%logL_Lya_HII = -35.
         baseGrid%cell(i,j,k)%logL_Lya_cool= -35.
         baseGrid%cell(i,j,k)%V         = (/0e0,0e0,0e0/)
         baseGrid%cell(i,j,k)%level     = 0
         baseGrid%cell(i,j,k)%parent    => baseGrid
         nullify(baseGrid%cell(i,j,k)%cell)
      enddo
   enddo
enddo

do level = 1, nlevels                           !For each level
   do icell = 1, readLevel(level)%ncell         !and each cell on this level

      x0 = readLevel(level)%pos(icell,1)        !Cell center coordinates if box
      y0 = readLevel(level)%pos(icell,2)        !were [0.0 - 1.0]
      z0 = readLevel(level)%pos(icell,3)        !

      i0 = int(x0*nx) + 1                       !
      j0 = int(y0*ny) + 1                       ! 1-128 (well, baseres)
      k0 = int(z0*nz) + 1                       !

      xnew = x0*float(nx) - float(i0-1)         !Fractional position in base
      ynew = y0*float(ny) - float(j0-1)         !cell (=> base cell has [.5,.5,.5]
      znew = z0*float(nz) - float(k0-1)         !L=1 cell has e.g. [.25,.75,.75]

      call placeCellProject(baseGrid%cell(i0,j0,k0),       &!
                            level,xnew,ynew,znew,          &!
                            readLevel(level)%lT(icell),    &!
                            readLevel(level)%lnH(icell),   &!
                            readLevel(level)%lx(icell),    &!
                          ! readLevel(level)%abun(icell,1),&! !He
                          ! readLevel(level)%abun(icell,2),&! !C
                          ! readLevel(level)%abun(icell,3),&! !N
                          ! readLevel(level)%abun(icell,4),&! !O
                          ! readLevel(level)%abun(icell,5),&! !Mg
                          ! readLevel(level)%abun(icell,6),&! !Si
                          ! readLevel(level)%abun(icell,7),&! !S
                          ! readLevel(level)%abun(icell,8),&! !Ca
                          ! readLevel(level)%abun(icell,9),&! !Fe
                            readLevel(level)%abun(icell,1),&! !Z
                            readLevel(level)%vel(icell,1), &!
                            readLevel(level)%vel(icell,2), &!
                            readLevel(level)%vel(icell,3), &!
                            readLevel(level)%LLY(icell,1), &!
                            readLevel(level)%LLY(icell,2))  !
   enddo
enddo

do level = 1, nlevels
   deallocate(readLevel(level)%pos)
   deallocate(readLevel(level)%lT)
   deallocate(readLevel(level)%lnH)
   deallocate(readLevel(level)%LLY)
   deallocate(readLevel(level)%lx)
   deallocate(readLevel(level)%vel)
enddo
deallocate(readLevel)

! Store Alex's hierarchy with space-filling curve in binary format

icosmic = 0
do i = 1, nx
   do j = 1, ny
      do k = 1, nz
         call countCells(baseGrid%cell(i,j,k))
      enddo
   enddo
enddo
ncosmic = icosmic
write(*,*) 'total number of cells =', ncosmic

  allocate (cellArrayLevel(ncosmic))
  allocate (cellArrayLevelORIG(ncosmic))
  allocate (cellArrayXpos(ncosmic))
  allocate (DammitAlex(ncosmic))
  allocate (cellArrayYpos(ncosmic))
  allocate (cellArrayZpos(ncosmic))
  allocate (cellArrayHI(ncosmic))
! allocate (cellArrayHeI(ncosmic))
! allocate (cellArrayHeII(ncosmic))
  allocate (cellArrayHydrogen(ncosmic))
  allocate (cellArrayx(ncosmic))
! allocate (cellArrayHelium(ncosmic))
! allocate (cellArrayCarbon(ncosmic))
! allocate (cellArrayNitrogen(ncosmic))
! allocate (cellArrayOxygen(ncosmic))
! allocate (cellArrayMagnesium(ncosmic))
! allocate (cellArraySilicon(ncosmic))
! allocate (cellArraySulfur(ncosmic))
! allocate (cellArrayCalcium(ncosmic))
! allocate (cellArrayIron(ncosmic))
  allocate (cellArrayZ(ncosmic))
  allocate (cellArrayTemp(ncosmic))
  allocate (cellArrayDensity(ncosmic))
  allocate (cellArrayXvel(ncosmic))
  allocate (cellArrayYvel(ncosmic))
  allocate (cellArrayZvel(ncosmic))
  allocate (cellArraylogL_Lya_HII(ncosmic))
  allocate (cellArraylogL_Lya_cool(ncosmic))

icosmic = 0
do i = 1, nx
  do j = 1, ny
    do k = 1, nz
      call writeCell(baseGrid%cell(i,j,k),0)
    enddo
  enddo
enddo

icosmic = 0
do i = 1, nx
   xpos = (DBLE(i)-0.5d0) / DBLE(nx)
   do j = 1, ny
      ypos = (DBLE(j)-0.5d0) / DBLE(ny)
      do k = 1, nz
         zpos = (DBLE(k)-0.5d0) / DBLE(nz)
         call computeCellCoordinates(xpos,ypos,zpos,0,1d0/DBLE(nx))
      enddo
   enddo
enddo

itmp = 0
do icosmic = 1, ncosmic
   cellArrayXpos(icosmic) = cellArrayXpos(icosmic)*(xb-xa) + xa
! if (cellArrayXpos(icosmic) .gt. 1d30) stop "+ inf"
! if (cellArrayXpos(icosmic) .lt. -1d30) stop "- inf"
   cellArrayYpos(icosmic) = cellArrayYpos(icosmic)*(yb-ya) + ya
   cellArrayZpos(icosmic) = cellArrayZpos(icosmic)*(zb-za) + za
   itmp = max(itmp,cellArrayLevel(icosmic))
enddo
print*, 'base grid level =', cellArrayLevel(1)
print*, 'max level of refinement =', itmp ! 0 - itmp

loT = 20.
cc  = 0
do i = 1,ncosmic
  if (cellArrayTemp(i) .lt. loT) then
    cellArrayTemp(i) = loT
    cc = cc + 1
  endif
enddo
print*, cc, 'cells were upgraded to T = ', loT

!================================================================

! -----------------------  write CellDataPreLum.bin

sumHII  = 0d0
sumcool = 0d0
do i = 1,ncosmic
  sumHII  = sumHII  + 10d0**dble(cellArraylogL_Lya_HII(i))
  sumcool = sumcool + 10d0**dble(cellArraylogL_Lya_cool(i))
enddo
write(*,'(a,es10.3,a)') 'L_Lya_cool,tot: ',sumcool          , ' erg/s'
write(*,'(a,es10.3,a)') 'L_Lya_HII,tot:  ',sumHII           , ' erg/s'
write(*,'(a,es10.3,a)') 'L_FUV,tot:      ',Lya2FUV * sumHII , ' erg/s/Hz'

do i = 1,ncosmic
  !---------------------------------------------------------------------!
  !              IMPORTANT! Pay attention to the following              !
  !                                                                     !
  cellArraylogL_Lya_cool(i) = log10(10d0**cellArraylogL_Lya_HII(i) &    !
                                  + 10d0**cellArraylogL_Lya_cool(i))    !
  !  Now cellArraylogL_Lya_cool contains both HII and cooling           !
  !                                                                     !
  cellArraylogL_Lya_HII(i) = cellArraylogL_Lya_HII(i) + logLya2FUV      !
  !  Now cellArraylogL_Lya_HII contains FUV                             !
  !                                                                     !
  !---------------------------------------------------------------------!

  cellArrayXpos(i) = cellArrayXpos(i) + R_box
  cellArrayYpos(i) = cellArrayYpos(i) + R_box
  cellArrayZpos(i) = cellArrayZpos(i) + R_box
  !Position are now in [0,D_box] rather than [-R_box,R_box]
enddo

print*, "Writing '", trim(fileout), "'"
open(14,file=trim(fileout),status='replace',form='unformatted')
  write(14) ncosmic, xb-xa, nx,ny,nz                      !Tot. #cells, phys size in kpc, resolution
  write(14) (int(cellArrayLevel(i),1),  i=1,ncosmic)              !L_min = 0
  write(14) (cellArrayXpos(i),          i=1,ncosmic) !in kpc; 0 - D_box
  write(14) (cellArrayYpos(i),          i=1,ncosmic) !in kpc
  write(14) (cellArrayZpos(i),          i=1,ncosmic) !in kpc
  !-------------------------------------------------------------------------!
  !                                                                         !
  !  ATTENTION: This is log(Lya,tot) in erg/s, i.e. both HII and cooling:   !
  write(14) (cellArraylogL_Lya_cool(i), i=1,ncosmic)                        !
  !                                                                         !
  !  ATTENTION: This is log(FUV) in erg/s/Hz                                !
  write(14) (cellArraylogL_Lya_HII(i),  i=1,ncosmic)                        !
  !                                                                         !
  !-------------------------------------------------------------------------!
  write(14) (cellArrayDensity(i),       i=1,ncosmic)
! write(14) (cellArrayHI(i),            i=1,ncosmic)
  write(14) (cellArrayHydrogen(i),      i=1,ncosmic)
  write(14) (cellArrayx(i),             i=1,ncosmic)
  write(14) (cellArrayTemp(i),          i=1,ncosmic)
  write(14) (cellArrayXvel(i),          i=1,ncosmic) !in km/s
  write(14) (cellArrayYvel(i),          i=1,ncosmic) !in km/s
  write(14) (cellArrayZvel(i),          i=1,ncosmic) !in km/s
! write(14) (cellArrayHelium(i),        i=1,ncosmic)
! write(14) (cellArrayCarbon(i),        i=1,ncosmic)
! write(14) (cellArrayNitrogen(i),      i=1,ncosmic)
! write(14) (cellArrayOxygen(i),        i=1,ncosmic)
! write(14) (cellArrayMagnesium(i),     i=1,ncosmic)
! write(14) (cellArraySilicon(i),       i=1,ncosmic)
! write(14) (cellArraySulfur(i),        i=1,ncosmic)
! write(14) (cellArrayCalcium(i),       i=1,ncosmic)
! write(14) (cellArrayIron(i),          i=1,ncosmic)
  write(14) (cellArrayZ(i),             i=1,ncosmic)
! write(14) (cellArrayHeI(i),           i=1,ncosmic)
! write(14) (cellArrayHeII(i),          i=1,ncosmic)
close(14)

V_tot = 0d0
V_HI  = 0d0
V_HII = 0d0

do i = 1,ncosmic
  dd = sqrt((cellArrayXpos(i)-R_box)**2 &
          + (cellArrayYpos(i)-R_box)**2 &
          + (cellArrayZpos(i)-R_box)**2)
  if (dd .gt. 0.9*R_box) cycle
  VV = (dx0/2d0**cellArrayLevel(i) * kpc)**3
  V_tot = V_tot + VV
  V_HI  = V_HI  + VV * cellArrayx(i)
  V_HII = V_HII + VV * (1d0 - cellArrayx(i))

  call random_number(R)
  if (R.lt.0.01 .and. cellArrayLevel(i).gt.2) write(1,*) &
  cellArraylogL_Lya_cool(i), cellArraylogL_Lya_HII(i)
! write(10+cellArrayLevel(i), *) cellArrayXpos(i),&
!                                cellArrayYpos(i),&
!                                cellArrayZpos(i)

! write(10+cellArrayLevel(i), *) cellArrayTemp(i),cellArrayDensity(i),cellArrayx(i),cellArrayHydrogen(i)
enddo

print*, 'V_tot', V_tot
print*, 'V_HI ', V_HI 
print*, 'V_HII', V_HII

end program BuildAMR

!*****************************************************************************

recursive subroutine placeCellProject(ParentCell,level,x0,y0,z0, &
                                      logtgas,lognh,logxneu,     &
                                    ! pCPHelium,                 &
                                    ! pCPCarbon,                 &
                                    ! pCPNitrogen,               &
                                    ! pCPOxygen,                 &
                                    ! pCPMagnesium,              &
                                    ! pCPSilicon,                &
                                    ! pCPSulfur,                 &
                                    ! pCPCalcium,                &
                                    ! pCPIron,                   &
                                      pCPZ,                      &
                                      vx,vy,vz,logHII,logcool)

use localDefinitions

implicit none
integer, intent(in) ::             level
integer ::                         inew, jnew, knew, i, j, k
real(kind=RealKind), intent(in) :: x0, y0, z0
real*4, intent(in) ::              logtgas, lognh, logxneu, vx,vy,vz,logHII,logcool
! real*4, intent(in) ::              pCPHelium
! real*4, intent(in) ::              pCPCarbon
! real*4, intent(in) ::              pCPNitrogen
! real*4, intent(in) ::              pCPOxygen
! real*4, intent(in) ::              pCPMagnesium
! real*4, intent(in) ::              pCPSilicon
! real*4, intent(in) ::              pCPSulfur
! real*4, intent(in) ::              pCPCalcium
! real*4, intent(in) ::              pCPIron
  real*4, intent(in) ::              pCPZ
real(kind=RealKind) ::             xnew, ynew, znew, nh, xneu!, nhe
type(zoneType), target ::          ParentCell
real*8::                             dx
integer::                          L

if (level.gt.1) then
   if (.not.ParentCell%refined) then
      allocate(ParentCell%cell(2,2,2))
      ParentCell%refined = .true.
      do i = 1, 2
         do j = 1, 2
            do k = 1, 2
               ParentCell%cell(i,j,k)%refined   = .false.
               ParentCell%cell(i,j,k)%level     = ParentCell%level + 1
               ParentCell%cell(i,j,k)%tgas      = ParentCell%tgas
               ParentCell%cell(i,j,k)%rho       = ParentCell%rho
               ParentCell%cell(i,j,k)%Hydrogen  = ParentCell%Hydrogen
               ParentCell%cell(i,j,k)%x         = ParentCell%x  
               ParentCell%cell(i,j,k)%HI        = ParentCell%HI
             ! ParentCell%cell(i,j,k)%HeI       = ParentCell%HeI
             ! ParentCell%cell(i,j,k)%HeII      = ParentCell%HeII
             ! ParentCell%cell(i,j,k)%Helium    = ParentCell%Helium   
             ! ParentCell%cell(i,j,k)%Carbon    = ParentCell%Carbon   
             ! ParentCell%cell(i,j,k)%Nitrogen  = ParentCell%Nitrogen 
             ! ParentCell%cell(i,j,k)%Oxygen    = ParentCell%Oxygen   
             ! ParentCell%cell(i,j,k)%Magnesium = ParentCell%Magnesium
             ! ParentCell%cell(i,j,k)%Silicon   = ParentCell%Silicon  
             ! ParentCell%cell(i,j,k)%Sulfur    = ParentCell%Sulfur   
             ! ParentCell%cell(i,j,k)%Calcium   = ParentCell%Calcium  
             ! ParentCell%cell(i,j,k)%Iron      = ParentCell%Iron     
               ParentCell%cell(i,j,k)%Z         = ParentCell%Z
               ParentCell%cell(i,j,k)%logL_Lya_HII = ParentCell%logL_Lya_HII
               ParentCell%cell(i,j,k)%logL_Lya_cool= ParentCell%logL_Lya_cool
               ParentCell%cell(i,j,k)%V         = ParentCell%V
               ParentCell%cell(i,j,k)%parent    => ParentCell
               nullify(ParentCell%cell(i,j,k)%cell)
            enddo
         enddo
      enddo
   endif

   if (x0.lt.0.5) then
      inew = 1
      xnew = 2.*x0
   else
      inew = 2
      xnew = 2.*x0 - 1.
   endif
   if (y0.lt.0.5) then
      jnew = 1
      ynew = 2.*y0
   else
      jnew = 2
      ynew = 2.*y0 - 1.
   endif
   if (z0.lt.0.5) then
      knew = 1
      znew = 2.*z0
   else
      knew = 2
      znew = 2.*z0 - 1.
   endif

   call placeCellProject(ParentCell%cell(inew,jnew,knew),level-1, &
        xnew,ynew,znew,logtgas,lognh,logxneu,                     &
                                     ! pCPHelium,                 &
                                     ! pCPCarbon,                 &
                                     ! pCPNitrogen,               &
                                     ! pCPOxygen,                 &
                                     ! pCPMagnesium,              &
                                     ! pCPSilicon,                &
                                     ! pCPSulfur,                 &
                                     ! pCPCalcium,                &
                                     ! pCPIron,                   &
                                       pCPZ,                      &
                                       vx,vy,vz,logHII,logcool)
else
   ParentCell%tgas      = 10.**logtgas
   nh                   = 10.**lognh
   xneu                 = 10.**logxneu
   ParentCell%rho       = nh * mh/psi
   ParentCell%Hydrogen  = nh
   ParentCell%x         = xneu
   ParentCell%HI        = nh * xneu
 ! nhe                  = (1.-psi) * ParentCell%rho / mhe
 ! ParentCell%HeI       = nhe * 1.
 ! ParentCell%HeII      = nhe * 0.
 ! ParentCell%Helium    = pCPHelium
 ! ParentCell%Carbon    = pCPCarbon
 ! ParentCell%Nitrogen  = pCPNitrogen
 ! ParentCell%Oxygen    = pCPOxygen
 ! ParentCell%Magnesium = pCPMagnesium
 ! ParentCell%Silicon   = pCPSilicon
 ! ParentCell%Sulfur    = pCPSulfur
 ! ParentCell%Calcium   = pCPCalcium
 ! ParentCell%Iron      = pCPIron
   ParentCell%Z         = pCPZ
   ParentCell%V(1)      = vx
   ParentCell%V(2)      = vy
   ParentCell%V(3)      = vz
   L                    = ParentCell%Level
   ParentCell%logL_Lya_HII = logHII
   ParentCell%logL_Lya_cool= logcool
endif

end subroutine placeCellProject

!-----------------------------------------------------------------------------

recursive subroutine countCells(currentCell)

use localDefinitions

implicit none
type(zoneType), target :: currentCell
integer :: i, j, k

if (currentCell%refined) then
   do i = 1, 2
      do j = 1, 2
         do k = 1, 2
            call countCells(currentCell%cell(i,j,k))
         enddo
      enddo
   enddo
else
   icosmic = icosmic + 1
endif

end subroutine countCells

!-----------------------------------------------------------------------------

recursive subroutine writeCell(currentCell,level)

use localDefinitions

implicit none
type(zoneType) :: currentCell
integer, intent(in) :: level
integer :: i, j, k

if (currentCell%refined) then
   do i = 1, 2
      do j = 1, 2
         do k = 1, 2
            call writeCell(currentCell%cell(i,j,k),level+1)
         enddo
      enddo
   enddo
else
   icosmic = icosmic + 1
   cellArrayLevel(icosmic)     = level
   cellArrayHI(icosmic)        = real(currentCell%HI)
 ! cellArrayHeI(icosmic)       = real(currentCell%HeI)
 ! cellArrayHeII(icosmic)      = real(currentCell%HeII)
 ! cellArrayHelium(icosmic)    = real(currentCell%Helium)
 ! cellArrayCarbon(icosmic)    = real(currentCell%Carbon)
 ! cellArrayNitrogen(icosmic)  = real(currentCell%Nitrogen)
 ! cellArrayOxygen(icosmic)    = real(currentCell%Oxygen)
 ! cellArrayMagnesium(icosmic) = real(currentCell%Magnesium)
 ! cellArraySilicon(icosmic)   = real(currentCell%Silicon)
 ! cellArraySulfur(icosmic)    = real(currentCell%Sulfur)
 ! cellArrayCalcium(icosmic)   = real(currentCell%Calcium)
 ! cellArrayIron(icosmic)      = real(currentCell%Iron)
   cellArrayZ(icosmic)         = real(currentCell%Z)
   cellArrayTemp(icosmic)      = real(currentCell%tgas)
   cellArrayDensity(icosmic)   = real(currentCell%rho)
   cellArrayHydrogen(icosmic)  = real(currentCell%Hydrogen)
   cellArrayx(icosmic)         = real(currentCell%x)
   cellArrayXvel(icosmic)      = real(currentCell%V(1))
   cellArrayYvel(icosmic)      = real(currentCell%V(2))
   cellArrayZvel(icosmic)      = real(currentCell%V(3))
   cellArraylogL_Lya_HII(icosmic) = real(currentCell%logL_Lya_HII)
   cellArraylogL_Lya_cool(icosmic)= real(currentCell%logL_Lya_cool)
endif

end subroutine writeCell

!-----------------------------------------------------------------------------

recursive subroutine computeCellCoordinates(x0,y0,z0,level,cellSize)

use localDefinitions

implicit none
integer, intent(in) :: level
real(kind=RealKind), intent(in) :: x0, y0, z0, cellSize
integer :: i, j, k
real(kind=RealKind) :: xnew, ynew, znew

icosmic = icosmic + 1

if (icosmic .eq. 10151) write(*,*) x0,y0,z0,level,cellSize

if (cellArrayLevel(icosmic).eq.level) then
   cellArrayXpos(icosmic) = x0                  !
   cellArrayYpos(icosmic) = y0                  ! Center of mother cells
   cellArrayZpos(icosmic) = z0                  !
else
   if (cellArrayLevel(icosmic).gt.level) then
      icosmic = icosmic - 1
      do i = 1, 2
         if (i.eq.1) then
            xnew = x0 - 0.25*cellSize
         else
            xnew = x0 + 0.25*cellSize
         endif
         do j = 1, 2
            if (j.eq.1) then
               ynew = y0 - 0.25*cellSize
            else
               ynew = y0 + 0.25*cellSize
            endif
            do k = 1, 2
               if (k.eq.1) then
                  znew = z0 - 0.25*cellSize
               else
                  znew = z0 + 0.25*cellSize
               endif
               call computeCellCoordinates(xnew,ynew,znew,level+1,cellSize/2d0)
            enddo
         enddo
      enddo
   else
      write(*,*) 'error in levels', icosmic, cellArrayLevel(icosmic), level
      stop
   endif
endif

end subroutine computeCellCoordinates

!-----------------------------------------------------------------------------

function dx(L)

use localDefinitions

real*8::  dx
integer:: L

dx  = dx0 / 2d0**L

end function

!-----------------------------------------------------------------------------
