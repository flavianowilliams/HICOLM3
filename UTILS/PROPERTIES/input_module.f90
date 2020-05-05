module input_module

    implicit none

  integer nstp,nvar,nat,nopt,spctot
  real(8) dtime

contains

  subroutine entrada

  implicit none

  read(1,*)nstp,dtime,nvar,nat,spctot

  write(*,*)'RDF -> 1'
  write(*,*)'Medias -> 2'

!    nopt=2   
  read(*,*)nopt

 end subroutine entrada

end module input_module
