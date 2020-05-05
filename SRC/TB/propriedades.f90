!
! This file is part of the HICOLM distribution (https://github.com/flavianowilliams/HICOLM).
!
! Copyright (c) 2019 Flaviano Williams Fernandes.
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, version 3.
!
! This program is distributed in the hope that it will be useful, but
! WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
! General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program. If not, see <http://www.gnu.org/licenses/>.
!
module propriedades
  !******************************************************************************************
  !Modulo do calculo das propriedades realicionadas com a estrutura eletronica              *
  !                                                                                         *
  ! Flaviano Williams Fernandes, 05 de dezembro de 2018                                     *
  !******************************************************************************************
  use sistema
  use input
  use estrutura
  use brillouin
  use diagonalize
  use matriz
  use optimizing

contains

  subroutine solving
    !****************************************************************************************
    ! Flag para o calculo das propriedades eletronicas                                      *
    !****************************************************************************************
    implicit none

    !    integer i

    select case(prop)
    case('EBAND')
       !       write(6,*)('#',i=1,76)
       !       write(6,*)'$$$$ ',('PROPRIEDADES',i=1,5),' $$$$'
       !       write(6,*)('#',i=1,76)
       !       write(6,*)
       call ebands ! calculando da energia eletronica
    case('BANDS')
       !       write(6,*)('#',i=1,76)
       !       write(6,*)'$$$$ ',('PROPRIEDADES',i=1,5),' $$$$'
       !       write(6,*)('#',i=1,76)
       !       write(6,*)
       call overlap ! calculando elementos da matriz S
    case('BANDS3D')
       !       write(6,*)('#',i=1,76)
       !       write(6,*)'$$$$ ',('PROPRIEDADES',i=1,5),' $$$$'
       !       write(6,*)('#',i=1,76)
       !       write(6,*)
       call overlap_3d ! calculando elementos da matriz S
    case('SURFACE')
       !       write(6,*)('#',i=1,76)
       !       write(6,*)'$$$$ ',('PROPRIEDADES',i=1,5),' $$$$'
       !       write(6,*)('#',i=1,76)
       !       write(6,*)
       call fermi_surface ! calculando superficie de Fermi
    end select

    return

  end subroutine solving

  subroutine overlap
    !****************************************************************************************
    ! Diagonalizacao da matriz secular                                                      *
    !****************************************************************************************
    implicit none

    integer nkpt,i,j,k,nz,nz0,chv
    real(8) kpt(nkmax,3),mkpt(nbandmax),kptt(3)
    real(8) knorm,dsvres

    open(4,file='bandas.dat',status='unknown')

    !-definindo pontos K

    call high_symmetry(1,kpt,mkpt,nkpt)

    !-calculo das matrizes H e S

    chv=1
    nz0=0
    knorm=0.d0
    dsvres=0.d0
    do nz=1,nkpt
       do k=1,3
          kptt(k)=kpt(nz,k)
       end do

       !-calculando E(k)

       call bands(kptt)

       knorm=knorm+mkpt(chv)/dband(chv)

       write(4,'(100f8.3)')knorm,((bnd(j)-fermi)*econv,j=1,norb)

       if(nz.ge.(nz0+dband(chv)))then
          chv=chv+1
          nz0=nz
       end if

       !-interpolando pontos pelo método de Lagrange

       do i=1,nkrlx
          do j=1,dkrlx(i)
             if(nz.eq.kz(i,j))then
                dsvres=dsvres+abs(bndfit(nz,nrlx(i))-bnd(ntblx(i)))
             end if
          end do
       end do

    end do

    if(nopt.ne.0)write(6,100)'Erro(eV):',dsvres/econv

    write(6,*)

    return

100 format(5x,a9,f12.4)

  end subroutine overlap

  subroutine overlap_3d
    !**************************************************************************************
    ! Calculo tridimensional da energia da camada de eletronica                           *
    !**************************************************************************************
    implicit none

    integer nkpt,i,j,k,l,nz
    real(8) kpt(nkmax,3),kptt(3)

    open(3,file='bandas-3d.dat',status='unknown')

    !-definindo pontos K

    nkpt=1
    do i=1,mesh
       do j=1,mesh
          kpt(nkpt,1)=(i-0.5*mesh)*3.0d0*pi/(a*mesh)
          kpt(nkpt,2)=(j-0.5*mesh)*3.0d0*pi/(b*mesh)
          kpt(nkpt,3)=2.0d0*pi/b
          nkpt=nkpt+1
       end do
    end do

    nkpt=nkpt-1

    !-calculo das matrizes H e S

    do nz=1,nkpt
       do k=1,3
          kptt(k)=kpt(nz,k)+kptt0(k)
       end do

       !-calculando E(k)

       call bands(kptt)

       write(3,'(11f12.4)')(kptt(i)*kconv,i=1,3),(((bnd(l)-fermi)*econv),l=1,norb)
    end do

    return

  end subroutine overlap_3d

  subroutine fermi_surface
    !****************************************************************************************
    !Calculo da superficie de Fermi                                                         *
    !****************************************************************************************
    implicit none

    integer i,j,k,l,nkpt,ix,jx,kx
    real(8) kptt(3),kpt(nkmax,3),ek(100),dk1(3),dk2(3),dk3(3)

    open(2,file='HICOLM.BXSF',status='unknown')

    do i=1,3
       dk1(i)=vr(1,i)/mesh
       dk2(i)=vr(2,i)/mesh
       dk3(i)=vr(3,i)/mesh
    end do

    nkpt=1
    do i=0,mesh
       do j=0,mesh
          do k=0,mesh
             kpt(nkpt,1)=(i)*dk1(1)+(j)*dk2(1)+(k)*dk3(1)
             kpt(nkpt,2)=(i)*dk1(2)+(j)*dk2(2)+(k)*dk3(2)
             kpt(nkpt,3)=(i)*dk1(3)+(j)*dk2(3)+(k)*dk3(3)
             nkpt=nkpt+1
          end do
       end do
    end do

    nkpt=nkpt-1

    write(2,*)'BEGIN_INFO'
    write(2,*)'  #'
    write(2,*)'  # This is a Band_XCRYSDEN-Structure-File'
    write(2,*)'  # aimed for Visualization of Fermi Surface'
    write(2,*)'  #'
    write(2,*)'  # Case: Tight-binding calculation'
    write(2,*)'  #'
    write(2,*)'  # Launch as: xcrysden --bxsf TBMC.BXSF'
    write(2,*)'  #'
    write(2,*)'  Fermi Energy:',fermi/econv
    write(2,*)'END_INFO'
    write(2,*)
    write(2,*)'BEGIN_BLOCK_BANDGRID_3D'
    write(2,*)'Tigh-binding_1st_BZ'
    write(2,*)'BEGIN_BANDGRID_3D'
    write(2,'(4x,i3)')1
    write(2,'(4x,3i15)')nkpt,nkpt,nkpt
    write(2,'(4x,3f15.8)')(kptt0(i)/kconv,i=1,3)
    write(2,'(4x,3f15.8)')(vr(1,i)/kconv,i=1,3)
    write(2,'(4x,3f15.8)')(vr(2,i)/kconv,i=1,3)
    write(2,'(4x,3f15.8)')(vr(3,i)/kconv,i=1,3)

!    do l=1,norb
    l=lfermi
    write(4,'(4x,a5,i5)')'BAND:',l

    do ix=1,nkpt
       kptt(1)=kpt(ix,k)
       do jx=1,nkpt
          kptt(2)=kpt(jx,k)
          do kx=i,nkpt
             kptt(3)=kpt(kx,3)
             call bands(kptt)
             ek(k)=bnd(l)
             k=k+1
          end do
          write(2,'(4x,3f15.8)')(ek(i),i=1,nkpt)
       end do
       write(2,*)
    end do
 !end do

    write(2,*)'END_BANDGRID_3D'
    write(2,*)'END_BLOCK_BANDGRID_3D'

    return

  end subroutine fermi_surface

  subroutine ebands
    !***************************************************************************************
    ! Calculo da energia eletronica                                                        *
    !***************************************************************************************
    implicit none

    integer i,j,l,k
    real(8) kptt(3),dkx,dky,dkz
    real(8) kxshift,kyshift,kzshift,ebs

    real(8), allocatable :: sum(:)

    allocate(sum(norb))

    !-calculo da integral por simpson

    dkx=(kxmax-kxmin)/mpack(1)
    dky=(kymax-kymin)/mpack(2)
    dkz=(kzmax-kzmin)/mpack(3)

    kxshift=dkx*ksh(1)
    kyshift=dky*ksh(2)
    kzshift=dkz*ksh(3)

    do l=1,norb
       sum(l)=0.d0
    end do

    do l=1,norb
       do k=1,mpack(3)
          kptt(3)=kzshift+(k-0.5*mpack(3))*dkz
          do j=1,mpack(2)
             kptt(2)=kyshift+(j-0.5*mpack(2))*dky
             do i=1,mpack(1)
                kptt(1)=kxshift+(i-0.5*mpack(1))*dkx
                call bands(kptt)
                sum(l)=sum(l)+bnd(l)*fermiocp(l)
             end do
          end do
       end do
       sum(l)=dkx*dky*dkz*sum(l)
    end do

    ebs=0.d0
    do l=1,norb
       ebs=ebs+sum(l)
    end do

        write(6,'(5x,26a1)')('-',j=1,26)
        write(6,'(5x,a1,2(5x,a7))')'l','  Dens ','Energia'
        write(6,'(5x,26a1)')('-',j=1,26)
        do l=1,natom
           write(6,'(5x,i1,4(3x,f9.3))')l,frho(l)*rconv**(1.d0/3.d0),sum(l)*econv
        end do
        write(6,'(5x,26a1)')('-',j=1,26)
        write(6,'(5x,16x,f9.3)')ebs*econv
        write(6,*)

    return

  end subroutine ebands

  double precision function fermiocp(i)
    !****************************************************************************************
    ! Funcao antisimetrica de Fermi-Dirac                                                   *
    !****************************************************************************************
    implicit none

    integer i

    if(2*i.le.norb)then
       fermiocp=1.d0
    else
       fermiocp=0.d0
    end if

    return

  end function fermiocp

end module propriedades
