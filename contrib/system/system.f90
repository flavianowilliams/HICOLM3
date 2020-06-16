program x2x

  implicit none

  integer ntmolecmax,atmax,natommax

  parameter(ntmolecmax=10,atmax=1000)
  parameter(natommax=ntmolecmax*atmax)

  real(8) x,y,z,qmolec(ntmolecmax,atmax)

  integer natom,i,j,k,ntmolecs,nmolec(ntmolecmax),namolec(ntmolecmax),zat
  character(2) at,tpmolec(ntmolecmax,atmax)
  character(10) nomemolec(10)

  open(1,file='HICOLM.xyz',status='old')
  open(2,file='HICOLM.sys',status='unknown')

  write(*,*)'Qde de tipos de moleculas (H2O,NH3,Ar,...):'
  read(*,*)ntmolecs

  do i=1,ntmolecs
     write(*,*)'Molecula:',i,':'
     write(*,*)
     write(*,*)'Nome:'
     read(*,*)nomemolec(i)
     write(*,*)'Quantidade de moleculas:'
     read(*,*)nmolec(i)
     write(*,*)'Quantidade de atomos por molecula:'
     read(*,*)namolec(i)
     write(*,*)'Classifique cada atomo por tipo:'
     read(*,*)(tpmolec(i,j),j=1,namolec(i))
     write(*,*)'Cargas parciais:'
     read(*,*)(qmolec(i,j),j=1,namolec(i))
  end do

  read(1,*)
  read(1,*)

  write(2,'(i3)')ntmolecs

  natom=1
  do i=1,ntmolecs
     write(2,'(a10,2i5)')nomemolec(i),nmolec(i),namolec(i)
     write(2,'(20a3)')(tpmolec(i,j),j=1,namolec(i))
     write(2,'(20f8.4)')(qmolec(i,j),j=1,namolec(i))
     do j=1,nmolec(i)
        do k=1,namolec(i)
           read(1,*)at,x,y,z
           call zatom(at,zat)
           write(2,'(i5,3f16.8)')zat,x,y,z
           natom=natom+1
        end do
     end do
  end do

  natom=natom-1

  write(2,*)
  write(2,*)natom,'atomos.'

  return

end program x2x

subroutine zatom(at,zat)

  implicit none

  integer zat
  character(2) at

  select case(at)
  case('H ')
     zat=1
  case('He')
     zat=2
  case('Li')
     zat=3
  case('Be')
     zat=4
  case('B ')
     zat=5
  case('C ')
     zat=6
  case('N ')
     zat=7
  case('O ')
     zat=8
  case('F ')
     zat=9
  case('Ne')
     zat=10
  case('Na')
     zat=11
  case('Mg')
     zat=12
  case('Al')
     zat=13
  case('Si')
     zat=14
  case('P ')
     zat=15
  case('S ')
     zat=16
  case('Cl')
     zat=17
  end select

  return

end subroutine zatom
