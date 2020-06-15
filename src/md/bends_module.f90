!
!MIT License
!
!Copyright (c) 2020 flavianowilliams
!
!Permission is hereby granted, free of charge, to any person obtaining a copy
!of this software and associated documentation files (the "Software"), to deal
!in the Software without restriction, including without limitation the rights
!to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
!copies of the Software, and to permit persons to whom the Software is
!furnished to do so, subject to the following conditions:
!
!THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
!IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
!FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
!AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
!LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
!OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
!SOFTWARE.
!
module bends_module
  !******************************************************************************************
  ! Contribuicao de deformacao para o campo de força:                                       *
  ! - Energia potencial;                                                                    *
  ! - Forças atômicas;                                                                      *
  ! - Virial;                                                                               *
  ! - Stress.                                                                               *
  !                                                                                         *
  ! Flaviano Williams Fernandes, 05 de dezembro de 2018                                     *
  !******************************************************************************************

  use input
  use utils
  use estrutura
  use alloc_arrays

  integer, allocatable :: bendij(:,:)
  integer, allocatable :: bendim(:),bendib(:)

  integer nbendstp

  save bendij,bendim,bendib,nbendstp

contains

  subroutine bends_alloc
    !****************************************************************************************
    ! Alocacao de variaveis globais                                                         *
    !****************************************************************************************

    implicit none

    integer ierr

    allocate(bendij(3,nbends),bendim(nbends),bendib(nbends),stat=ierr)

    if(ierr.ne.0)stop 'bends_alloc: allocation failed'

    return

  end subroutine bends_alloc

  subroutine bends_convert
    !****************************************************************************************
    ! Conversao de unidades de medida:                                                      *
    ! Unidades de entrada ---> a.u.                                                         *
    !****************************************************************************************

    implicit none

    integer i,j

    !-convertendo unidades de medida

    do i=1,nmolec
       do j=1,bendscnt(i)
          select case(bends(i,j))
          case(1)
             parbend(i,j,1)=parbend(i,j,1)/econv
             parbend(i,j,2)=parbend(i,j,2)/aconv
          case(2)
             parbend(i,j,1)=parbend(i,j,1)/econv
             parbend(i,j,2)=parbend(i,j,2)/aconv
          end select
       end do
    end do

    return

  end subroutine bends_convert

  subroutine bends_counts
    !****************************************************************************************
    ! Contagem do número de deformacoes                                                     *
    !****************************************************************************************

    implicit none

    integer, allocatable :: chk(:,:)

    integer nx,nxx,i,j,k,ni,nj,nk,nkb,np

    !-alocando array

    nx=0
    do i=1,nmolec
       nx=max(nx,bendscnt(i))
    end do

    allocate(chk(nmolec,nbends))

    !-checando viabilidade dos bonds

    nkb=1

    do i=1,nmolec
       do j=1,bendscnt(i)
          select case(bends(i,j))
          case(1)
             nkb=2
          case(2)
             nkb=2
          end select
          chk(i,j)=1
          do k=1,nkb
             if(parbend(i,j,k).eq.0.d0)chk(i,j)=0
          end do
       end do
    end do

    !-calculando bends intramoleculares

    np=0
    nx=1
    do i=1,nmolec
       nxx=0
       do j=1,ntmolec(i)
          do k=1,bendscnt(i)
             ni=np+molbend(i,k,2)
             nj=np+molbend(i,k,1)
             nk=np+molbend(i,k,3)
             if(chk(i,k).eq.1)then
                bendij(1,nx)=ni
                bendij(2,nx)=nj
                bendij(3,nx)=nk
                bendim(nx)=i
                bendib(nx)=k
                nxx=nxx+1
                nx=nx+1
             end if
          end do
          np=np+nxmolec(i)
       end do
       bendsmlc(i)=nxx
    end do

    nbendstp=nx-1

    !-limpando memoria

    deallocate(chk)

    return

  end subroutine bends_counts

  subroutine bends_calc(enbend,virbend)
    !****************************************************************************************
    ! - Energia potencial;                                                                  *
    ! - Contribuicao para o virial;                                                         *
    ! - Contribuicao para o stress;                                                         *
    ! - Forças atômicas.                                                                    *
    !****************************************************************************************

    implicit none

    integer i,ia,ib,ic,im,in
    real(8) pot,fa,drij(3),drik(3),theta,dr1,dr2
    real(8) enbend,virbend

    do i=1,nbendstp

       ia=bendij(1,i)
       ib=bendij(2,i)
       ic=bendij(3,i)
       im=bendim(i)
       in=bendib(i)

       call mic(ia,ib,drij(1),drij(2),drij(3))
       call mic(ia,ic,drik(1),drik(2),drik(3))

       dr1=max(1.d-4,sqrt(drij(1)**2+drij(2)**2+drij(3)**2))
       dr2=max(1.d-4,sqrt(drik(1)**2+drik(2)**2+drik(3)**2))

       theta=max(1.d-4,acos((drij(1)*drik(1)+drij(2)*drik(2)+drij(3)*drik(3))/(dr1*dr2)))

       call bends_flags(im,in,theta,pot,fa)
       call bends_force(ia,ib,ic,drij,drik,dr1,dr2,theta,fa,virbend)

       enbend=enbend+pot

    end do

    return

  end subroutine bends_calc

  subroutine bends_flags(im,in,theta,pot,fa)
    !****************************************************************************************
    ! Compontente angular                                                                   *
    !****************************************************************************************

    implicit none

    integer im,in
    real(8) theta,fa,pot

    !-atribuindo valores iniciais

    pot=0.d0
    fa=0.d0

    !-calculo da componente dU/dtheta

    select case(bends(im,in))
    case(1)
       pot=0.5d0*parbend(im,in,1)*(theta-parbend(im,in,2))**2
       fa=parbend(im,in,1)*(theta-parbend(im,in,2))
    case(2)
       pot=parbend(im,in,1)*(theta-parbend(im,in,2))**2
       fa=2.d0*parbend(im,in,1)*(theta-parbend(im,in,2))
    end select

    return

  end subroutine bends_flags

  subroutine bends_force(i1,i2,i3,drij,drik,dr1,dr2,theta,fa,virbend)
    !****************************************************************************************
    ! - Contribuicao para o virial                                                          *
    ! - Contribuicao para o stress                                                          *
    ! - Forcas dos atomos i, j e k                                                          *
    !****************************************************************************************

    implicit none

    integer i,j,i1,i2,i3,ix(3)
    real(8) drij(3),drik(3),derij(3,3),fbi(3),fbj(3),fbk(3),dr1,dr2,theta,virbend,fa

    !-calculo da componente dtheta/dr

    ix(1)=i1
    ix(2)=i2
    ix(3)=i3

    do j=1,3
       do i=1,3
          derij(i,j)=(kronij(ix(i),ix(2))-kronij(ix(i),ix(1)))*drik(j)/(dr1*dr2) &
               +(kronij(ix(i),ix(3))-kronij(ix(i),ix(1)))*drij(j)/(dr1*dr2) &
               -cos(theta)*((kronij(ix(i),ix(2))-kronij(ix(i),ix(1)))*drij(j)/dr1**2 &
               +(kronij(ix(i),ix(3))-kronij(ix(i),ix(1)))*drik(j)/dr2**2)
       end do
    end do

    !-calculo das forças atomicas  i, j e k

    fbi(1)=fa*derij(1,1)/sin(theta)
    fbi(2)=fa*derij(1,2)/sin(theta)
    fbi(3)=fa*derij(1,3)/sin(theta)

    fbj(1)=fa*derij(2,1)/sin(theta)
    fbj(2)=fa*derij(2,2)/sin(theta)
    fbj(3)=fa*derij(2,3)/sin(theta)

    fbk(1)=fa*derij(3,1)/sin(theta)
    fbk(2)=fa*derij(3,2)/sin(theta)
    fbk(3)=fa*derij(3,3)/sin(theta)

    fax(i1)=fax(i1)+fbi(1)
    fay(i1)=fay(i1)+fbi(2)
    faz(i1)=faz(i1)+fbi(3)

    fax(i2)=fax(i2)+fbj(1)
    fay(i2)=fay(i2)+fbj(2)
    faz(i2)=faz(i2)+fbj(3)

    fax(i3)=fax(i3)+fbk(1)
    fay(i3)=fay(i3)+fbk(2)
    faz(i3)=faz(i3)+fbk(3)

    !-contribuicao para o virial

    virbend=virbend+(fbj(1)*drij(1)+fbj(2)*drij(2)+fbj(3)*drij(3))&
         +(fbk(1)*drik(1)+fbk(2)*drik(2)+fbk(3)*drik(3))

    !-contribuicao para o stress

    str(1)=str(1)-(fbj(1)*drij(1)+fbk(1)*drik(1))
    str(2)=str(2)-(fbj(2)*drij(2)+fbk(2)*drik(2))
    str(3)=str(3)-(fbj(3)*drij(3)+fbk(3)*drik(3))
    str(4)=str(4)-(fbj(3)*drij(2)+fbk(3)*drik(2))
    str(5)=str(5)-(fbj(1)*drij(3)+fbk(1)*drik(3))
    str(6)=str(6)-(fbj(2)*drij(1)+fbk(2)*drik(1))

    return

  end subroutine bends_force

end module bends_module
