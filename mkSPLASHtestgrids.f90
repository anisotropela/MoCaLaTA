
! A snippet to create test grids that look like SPLASH data
! ----------------------------------------------------------
!
! Note that the header won't have as many lines as a real SPLASH data file,
! but splash2jesperf90 can still read it.
!
!------------------------------------------------------------------------------

program mkSPLASHtestgrids
  implicit none
  integer, parameter::    sp = kind(1.0)
  integer, parameter::    dp = kind(1.0d0)
  integer::               i,j,k
  character(len=20)::     cni,nifmt
  integer, parameter::    nvar = 9
  character(len=200)::    filebase,testdir
  logical::               ex
  character(len=20)::     vars(nvar)
  integer::               v,ni,nj,nk,icell,ncell,lun
  real(sp), allocatable:: cube(:,:,:,:)
  real(dp), parameter::   kpc = 3.0856776d21!kpc/cm
  real(dp)::              Lbox,Rbox,dx,dy,dz,dcell
  real(dp)::              xcell(3)
  real(dp)::              NcolHI,nHI,nHII,T,Vout,starLya,coolLya,Z

  ! Set physical values
  ni      = 10                  !Resolution
  nj      = ni
  nk      = ni
  Lbox    = 20.                 !Box size in kpc
  Rbox    = Lbox / 2
  NcolHI  = 1e19                !HI column density from center to edge
  nHI     = NcolHI / (Rbox*kpc) !HI density (global)
  nHII    = 0.                  !HII density (global)
  T       = 1e4                 !Temperature (global)
  Vout    = 0.                  !Outflow velocity (global)
  starLya = 1e33                !Stellar Lya luminosity (central)
  coolLya = 0.                  !Cooling radiation (global)
  Z       = 0.                  !Metallicity (global)

  ! Filenames (filebase + vars(v) + '.dat'
  testdir  = 'testdir/'
  inquire(file=trim(testdir), exist=ex)
  if (.not. ex) call system('mkdir'//trim(testdir))
  filebase = trim(testdir) // 'SPLASHtestgrid_'
  vars(1)  = 'nHI'
  vars(2)  = 'nHII'
  vars(3)  = 'T'
  vars(4)  = 'vx'
  vars(5)  = 'vy'
  vars(6)  = 'vz'
  vars(7)  = 'starLya'
  vars(8)  = 'coolLya'
  vars(9)  = 'Z'

  !Initialize cubes
  allocate(cube(nvar,ni,nj,nk))
  cube = 0 !The following loop only assigns values inside sphere, so set to 0 outside

  ! Step through grids, assigning values
  dx = Lbox / ni
  dy = Lbox / nj
  dz = Lbox / nk

  do i=1,ni
    do j=1,nj
      do k=1,nk
        xcell = (/(i-.5)*dx, (j-.5)*dy, (k-.5)*dz/) - Rbox !Cell coordinates
        dcell = sqrt(sum(xcell**2))                        !Distance to center

        !If inside sphere, assign values
        if (dcell .lt. Rbox) then
          cube(1,i,j,k) = nHI
          cube(2,i,j,k) = nHII
          cube(3,i,j,k) = T
          cube(4,i,j,k) = xcell(1)/dcell * Vout
          cube(5,i,j,k) = xcell(2)/dcell * Vout
          cube(6,i,j,k) = xcell(3)/dcell * Vout
          cube(8,i,j,k) = coolLya
          cube(9,i,j,k) = Z
        endif

        !Emit from central 8 cells
        if (all(abs(xcell) .lt. (/dx,dy,dz/))) then
          cube(7,i,j,k) = starLya
        endif
      enddo
    enddo
  enddo

  ! Write files
  write(cni,'(i4.4)') ni
  nifmt = '(' // trim(cni) // 'es15.6)'
  do v=1,9
    lun = v + 10 ! Because Fortran...
    open(lun,file=trim(filebase)//trim(vars(v))//'.dat',form='formatted',status='replace',action='write')
    write(lun,'(a)') '# ' // trim(vars(v)) // ' test grid made to look like SPLASH data'
    write(lun,'(a)') '#'
    write(lun,'(a)') '# written in the form:'
    write(lun,'(a)') '#   do k=1,nz' !nx,ny,nz because that's what the SPLASH files look like
    write(lun,'(a)') '#      do j=1,ny'
    write(lun,'(a)') '#         write(*,*) (dat(i,j,k),i=1,nx)'
    write(lun,'(a)') '#      enddo'
    write(lun,'(a)') '#   enddo'
    write(lun,'(a)') '#'
    write(lun,'(a)') '# grid dimensions:'
    write(lun,'(a)') '#  nx   ny   nz'
    write(lun,'(3i5)') ni,nj,nk
    do k=1,nk
      do j=1,nj
        write(lun,nifmt) (cube(v,i,j,k), i=1,ni)
      enddo
    enddo
  enddo
end program mkSPLASHtestgrids
