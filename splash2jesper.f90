
! Documentation for the code splash2jesper.f90
! --------------------------------------------
!
! AUTHOR
! ------
!   Peter Laursen, January 23, 2019.
!
! PURPOSE
! -------
!   This code converts a batch of non-AMR SPLASH output files in ASCII format
!   to a format readable by the code BuildAMR.f90. The batch consists of 9
!   files containing the following data:
!     1   HI density
!     2   Temperature
!     3-5 Bulk velocity (1 files for each axis)
!     6   log(Lya luminosity) from stars
!     7   log(Lya luminosity) from cooling radiation
!     8   HII density
!     9   Metallicity
!  At some point I'll change this code to allow for not including file 7 and/or
!  files 8-9 (which are used for dust), but for now they must exist; if the
!  data doesn't exist, just create some fake files full of zeros.
!
! USAGE
! -----
!   To compile the code, use
!     > gfortran -O3 splash2jesper.f90 -o splash2jesper.x
!   The syntax is then
!     > ./splash2jesper.x splashness.in
!   where `splashness.in` is an input file of the format outlined in the
!   subroutine `readinput`.
!
! FORMAT
! ------
!   Note that the format of the SPLASH files are (after a line containing the
!   resolution) an (nj*nk) X ni array written with
!     do k=1,nk
!       do j=1,nj
!         write(*,*) (dat(i,j,k),i=1,ni)
!       enddo
!   enddo
!   In contrast, BuildAMR.f90 assumes a 1D array where the k index varies
!   fastest. Yes, I know...
!
! PERFORMANCE
! -----------
!   Reading a 150**3 grid takes ~2 sec. Writing is fast (because binary).
!
!   TODO: All the loops flattening the 3D SPLASH data to the 1D arrays needed
!   by BuildAMR should be doable using something like
!       1Darray = pack(cube, .true.)
!   perhaps together with transpose, but I haven't figured out how to do it in
!   the correct order, and the loops are pretty fast anyway.
!------------------------------------------------------------------------------

module kinds
  integer, parameter:: sp      = kind(1.0)
  integer, parameter:: dp      = kind(1.0d0)
end module kinds
!------------------------------------------------------------------------------

module globals
  use kinds
  character(200):: splashout,file1,nHIfile,nHIIfile,Tfile,vxfile,vyfile,vzfile,starLyafile,coolLyafile,Zfile,jesperfile
  character(20)::  f
  integer::        ni,nj,nk,icell,ncell
  integer::        ni1=-1,nj1=-1,nk1=-1
  real(sp), allocatable, dimension(:,:,:):: cube
  real(sp), allocatable, dimension(:,:) ::  pos, vel, abun, LLY
  real(sp), allocatable, dimension(:) ::    lT, lnH, lx ! <- total H num dens & neutral frac
  real(sp), allocatable, dimension(:) ::    nHI, nHIflat
  real(dp):: Lbox
end module globals
!------------------------------------------------------------------------------

program splash2jesper
  use kinds
  use globals
  implicit none

  f = 'formatted'
  call readinput

  call readcube(nHIfile);     call mkpos; call mknHI
  call readcube(Tfile);       call mkT
  call readcube(vxfile);      call mkv(1)
  call readcube(vyfile);      call mkv(2)
  call readcube(vzfile);      call mkv(3)
  call readcube(starLyafile); call mkLya(1)
  call readcube(coolLyafile); call mkLya(2)
  call readcube(nHIIfile);    call mknHandxHI
  call readcube(Zfile);       call mkZ(1)

  call writeJesper

end program splash2jesper
!------------------------------------------------------------------------------

subroutine readinput
  ! Read input file of the format indicated by the comments at the read
  ! statements 
  use globals
  implicit none
  character(200):: header

  call get_command_argument(1,splashout) !Use arg rather than redirect (i.e. no "<")
  open(1,file=trim(splashout),form='formatted',status='old',action='read')

  read(1,*) header
  read(1,*) Lbox         !Box size / kpc
  read(1,*) jesperfile   !Output filename, readable by BuildAMR.f90
  read(1,*) nHIfile      !Filename for HI number density / cm-3
  read(1,*) Tfile        !Filename for temperature / K
  read(1,*) vxfile       !Filename for bulk velocity in x-direction / km/s
  read(1,*) vyfile       !Filename for bulk velocity in y-direction / km/s
  read(1,*) vzfile       !Filename for bulk velocity in z-direction / km/s
  read(1,*) starLyafile  !Filename for log(Lya luminosity) from stars / erg/s (or just s-1)
  read(1,*) coolLyafile  !Filename for log(Lya luminosity) from cooling / erg/s (or just s-1)
  read(1,*) nHIIfile     !Filename for HII number density / cm-3
  read(1,*) Zfile        !Filename for metallicity / Zsun

  close(1)
end subroutine readinput
!------------------------------------------------------------------------------

subroutine readcube(cubefile)
  !Read SPLASH output file.
  !f should be either 'formatted' (for text files) or 'unformatted' (for binary)
  use globals
  implicit none
  character(200),intent(in):: cubefile
  character(200):: line
  integer:: i,j,k

  write(*,*) 'Loading ' // trim(cubefile) // '...'
  open(1,file=trim(cubefile),form=trim(f),status='old',action='read')

  ! Read header and array size
  if (f .eq. 'unformatted') then
    read(1) ni,nj,nk
  else
    do i=1,25
      read(1,*) line
      if (index(trim(line), '#') .ne. 1) then
        backspace(1)
        read(1,*) ni,nj,nk
        exit
      endif
      if (i.eq.25) stop "Doesn't seem like a SPLASH file..."
    enddo
  endif

  ! Save grid size and file name of first file for later checks
  if (ni1 .eq. -1) then
    ni1   = ni
    nj1   = nj
    nk1   = nk
    file1 = cubefile
    ncell = ni*nj*nk
    allocate(cube(ni,nj,nk))
  endif

  ! Check if current grid matches first grid
  if ((ni.ne.ni1) .or. (nj.ne.nj1) .or. (nk.ne.nk1)) then
    write(*,*) 'GRID ERROR:'
    write(*,*) '  Grid size of ' // trim(cubefile) // " doesn't match that of " // trim(file1)
    stop
  endif

  ! Read cube
  if (f == 'unformatted') then
    read(1) cube
  else
    do k=1,nk
      do j=1,nj
        read(1,*) (cube(i,j,k), i=1,ni)
      enddo
    enddo
  endif

  print*, trim(cubefile), cube(ni/2,nj/2,nk/2)

  close(1)
end subroutine readcube
!------------------------------------------------------------------------------

subroutine mkpos
  ! Construct an ncell X 3 array of cell center positions from the box size and
  ! resolution
  use globals
  implicit none
  real(dp):: dx,dy,dz
  integer:: i,j,k

  ! Cell size (in each direction... maybe some day we'll implement ni != nj != nk)
  dx = Lbox / ni
  dy = Lbox / nj
  dz = Lbox / nk

  allocate(pos(ncell,3))

  ! Calculate cell center positions
  icell = 0
  do k = 1,nk
    do j = 1,nj
      do i = 1,ni
        icell      = icell + 1
        pos(icell,1) = (i-.5) * dx - Lbox/2
        pos(icell,2) = (j-.5) * dy - Lbox/2
        pos(icell,3) = (k-.5) * dz - Lbox/2
      enddo
    enddo
  enddo
end subroutine mkpos
!------------------------------------------------------------------------------

subroutine mknHI
  !Flatten nHI cube to 1D array.
  use globals
  implicit none
  integer:: i,j,k

  where(cube .lt. 1e-33) cube = 1e-33 !Make sure there are no zeros
  allocate(nHI(ncell))
  allocate(nHIflat(ncell))

  icell = 0
  do i = 1,ni
    do j = 1,nj
      do k = 1,nk
        icell      = icell + 1
        nHI(icell) = cube(i,j,k)
      enddo
    enddo
  enddo

end subroutine mknHI
!------------------------------------------------------------------------------

subroutine mkT
  !Flatten T cube to 1D array.
  use globals
  implicit none
  integer:: i,j,k

  where(cube .lt. 1e-33) cube = 1e-33           !Make sure there are no zeros
  cube = log10(cube)                            !Since we want logT
  allocate(lT(ncell))

  icell = 0
  do i = 1,ni
    do j = 1,nj
      do k = 1,nk
        icell     = icell + 1
        lT(icell) = cube(i,j,k)
      enddo
    enddo
  enddo
end subroutine mkT
!------------------------------------------------------------------------------

subroutine mkv(axis)
  !Flatten velocity to 1D array.
  use globals
  implicit none
  integer, intent(in):: axis
  integer:: i,j,k

  if (.not. allocated(vel)) allocate(vel(ncell,3))

  icell = 0
  do i = 1,ni
    do j = 1,nj
      do k = 1,nk
        icell           = icell + 1
        vel(icell,axis) = cube(i,j,k)
      enddo
    enddo
  enddo
end subroutine mkv
!------------------------------------------------------------------------------

subroutine mkLya(process)
  !Flatten L_Lya to 1D array.
  use globals
  implicit none
  integer, intent(in):: process                 !1:stars or 2:cooling
  integer:: i,j,k

! where(cube .lt. 1e-33) cube = 1e-33           !Make sure there are no zeros
! cube = log10(cube)                            !Since we want logL_Lya
  if (.not. allocated(LLY)) allocate(LLY(ncell,2))

  icell = 0
  do i = 1,ni
    do j = 1,nj
      do k = 1,nk
        icell              = icell + 1
        LLY(icell,process) = cube(i,j,k)
      enddo
    enddo
  enddo
end subroutine mkLya
!------------------------------------------------------------------------------

subroutine mknHandxHI
  !Create HII and xHI arrays from H cube and (previously made) HI array
  use globals
  implicit none
  integer::  i,j,k
  real(dp):: nHII,nH,xHI

  where(cube .lt. 1e-33) cube = 1e-33 !Make sure there are no zeros
  allocate(lnH(ncell))
  allocate(lx(ncell))

  icell = 0
  do i = 1,ni
    do j = 1,nj
      do k = 1,nk
        icell      = icell + 1
        nHII       = cube(i,j,k)       !HII number density
        nH         = nHI(icell) + nHII !Total H number density
        xHI        = nHI(icell) / nH   !Neutral fraction
        lnh(icell) = log10(nH)         !Logarithms for Jesper
        lx(icell)  = log10(xHI)        !     -"-
      enddo
    enddo
  enddo
end subroutine mknHandxHI
!------------------------------------------------------------------------------

subroutine mkZ(element)
  !Flatten Z cube to 1D array.
  use globals
  implicit none
  integer, intent(in):: element !atm just 1, but may be extended to 1-7
  integer:: i,j,k

  where(cube .lt. 1e-33) cube = 1e-33               !Make sure there are no zeros
  if (.not. allocated(abun)) allocate(abun(ncell,1))!Use only one metallicity, but write 7 times

  icell = 0
  do i = 1,ni
    do j = 1,nj
      do k = 1,nk
        icell               = icell + 1
        abun(icell,element) = cube(i,j,k)
      enddo
    enddo
  enddo
end subroutine mkZ
!------------------------------------------------------------------------------

subroutine writeJesper
  use globals
  implicit none

  open(14,file=trim(jesperfile), status='replace',form='unformatted')

  write(*,*) 'Writing...'
  write(14) 1 !Only 1 level of refinement, since this is non-AMR
  write(14) ncell

  write(14) (pos(icell,1),  icell=1,ncell)
  write(14) (pos(icell,2),  icell=1,ncell)
  write(14) (pos(icell,3),  icell=1,ncell)
  write(14) (lT(icell),     icell=1,ncell)
  write(14) (lnH(icell),    icell=1,ncell)
  write(14) (lx(icell),     icell=1,ncell)
! write(14) (abun(icell,1), icell=1,ncell)
  write(14) (abun(icell,1), icell=1,ncell)
  write(14) (abun(icell,1), icell=1,ncell)
  write(14) (abun(icell,1), icell=1,ncell)
  write(14) (abun(icell,1), icell=1,ncell)
  write(14) (abun(icell,1), icell=1,ncell)
  write(14) (abun(icell,1), icell=1,ncell)
! write(14) (abun(icell,1), icell=1,ncell)
  write(14) (abun(icell,1), icell=1,ncell)
  write(14) (vel(icell,1),  icell=1,ncell)
  write(14) (vel(icell,2),  icell=1,ncell)
  write(14) (vel(icell,3),  icell=1,ncell)
  write(14) (LLY(icell,1:2),icell=1,ncell)
! write(14) (LLY(icell,2),  icell=1,ncell)

  close(14)
end subroutine writeJesper
!------------------------------------------------------------------------------
