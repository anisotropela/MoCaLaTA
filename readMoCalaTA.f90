! Test snippet for reading an IFU cube produced by MoCaLaTA
!
! Usage:
!   > gfortran -O3 readMoCaLaTA.f90 -o readMoCaLaTA.x
!   > ./readMoCaLaTA.x filename # (e.g. testdir/xp.bin)
!
!------------------------------------------------------------------------------

module kinds
  integer, parameter:: sp      = kind(1.0)
  integer, parameter:: dp      = kind(1.0d0)
end module kinds
!------------------------------------------------------------------------------

program readMoCaLaTA
  use kinds
  implicit none
  character(len=200)::    binfile
  real(dp)::              n_sofar,L_tot,z,d_L,R_obs,BW(2),obspar(4)
  integer::               pix,SpecRes1D,SpecRes2D,n
  real(dp), allocatable:: spec(:),CCD(:,:,:)

  call get_command_argument(1,binfile)

  open (1, file=trim(binfile),form='unformatted',status='old',action='read')
  read(1) n_sofar
  read(1) pix,SpecRes1D,SpecRes2D,n,L_tot,BW,z,d_L,R_obs
  read(1) obspar ! <NHI>, <EBV>, NHI_QSO, EBV_QSO
  allocate(spec(SpecRes1D))
  allocate(CCD(pix,pix,SpecRes2D))
  read(1) spec
  read(1) CCD
  close(1)

  print*, 'n_sofar     ', n_sofar
  print*, 'pix         ', pix
  print*, 'SpecRes1D   ', SpecRes1D
  print*, 'SpecRes2D   ', SpecRes2D
  print*, 'n           ', n
  print*, 'L_tot/erg/s ', L_tot
  print*, 'BW/AA       ', BW
  print*, 'z           ', z
  print*, 'd_L/kpc     ', d_L
  print*, 'R_obs/kpc   ', R_obs
  print*, 'obspar      ', obspar
  print*, 'spec(1:4)   ', spec(1:4)
  print*, 'CCD(1:4,1,1)', CCD(1:4,1,1)

end program readMoCaLaTA
