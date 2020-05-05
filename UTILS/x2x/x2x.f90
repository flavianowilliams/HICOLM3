program x2x

  implicit none

  integer ntmolecmax,atmax,natommax

  parameter(ntmolecmax=10,atmax=100)
  parameter(natommax=ntmolecmax*atmax*100)

  integer idna(natommax),fztp(natommax),atp(natommax)
  real(8) x(natommax),y(natommax),z(natommax)

  integer natom,i,j,k,ntmolecs,nmolec(ntmolecmax),namolec(ntmolecmax)
  integer tpmolec(ntmolecmax,atmax)
  real(8) v(6,6)
  character(2) at

  open(1,file='HICOLM.xyz',status='old')
  open(2,file='HICOLM.str',status='unknown')

  write(*,*)'Qde de tipos de moleculas (H2O,NH3,Ar,...):'
  read(*,*)ntmolecs

  do i=1,ntmolecs
     write(*,*)'Molecula:',i,':'
     write(*,*)'Quantidade de moleculas:'
     read(*,*)nmolec(i)
     write(*,*)'Quantidade de atomos por molecula:'
     read(*,*)namolec(i)
     write(*,*)'Classifique cada atomo por tipo:'
     read(*,*)(tpmolec(i,j),j=1,namolec(i))
  end do

  read(1,*)
  read(1,*)

  natom=1
  do i=1,ntmolecs
     do j=1,nmolec(i)
        do k=1,namolec(i)
           fztp(natom)=1
           atp(natom)=tpmolec(i,k)
           read(1,*)at,x(natom),y(natom),z(natom)
           call type(at,idna(natom))
           natom=natom+1
        end do
     end do
  end do

  natom=natom-1

  write(*,*)natom,'atomos.'

  call geometria(natom,idna,x,y,z,atp,fztp)

  return

end program x2x

subroutine type(at,atx)

  implicit none

  integer atx
  character(2) at

  select case(at)
  case('H')
     atx=1
  case('C')
     atx=6
  case('N')
     atx=7
  case('O')
     atx=8
  case('P')
     atx=15
  case('Ar')
     atx=18
  end select

  return

end subroutine type

 subroutine geometria(natom,idna,x,y,z,atp,fztp)

   implicit none

   integer i,natom,idna(natom),atp(natom),fztp(natom)
   real(8) v(6,6),x(natom),y(natom),z(natom)

   do i=1,natom
      write(2,'(i5,3f14.8,2x,2i7)')idna(i),x(i),y(i),z(i),atp(i),fztp(i)
   end do

    return

  end subroutine geometria
