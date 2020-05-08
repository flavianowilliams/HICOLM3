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
module bonds_module
  !******************************************************************************************
  ! Contribuicoes de estiramento para o campo de força:                                     *
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

  integer, allocatable :: bondij(:,:)
  integer, allocatable :: bondim(:),bondib(:)

  integer nbondstp

  save bondij,bondim,bondib,nbondstp

contains

  subroutine bonds_alloc
    !****************************************************************************************
    ! Alocacao de variaveis globais                                                         *
    !****************************************************************************************

    implicit none

    integer ierr

    allocate(bondij(2,nbonds),bondim(nbonds),bondib(nbonds),stat=ierr)

    if(ierr.ne.0)stop 'bonds_alloc: allocation failed'

    return

  end subroutine bonds_alloc

  subroutine bonds_convert
    !****************************************************************************************
    ! Conversao de unidades de medida:                                                      *
    ! Unidades de entrada ---> a.u.                                                         *
    !****************************************************************************************

    implicit none

    integer i,j

    !-convertendo unidades de medida

    do i=1,nmolec
       do j=1,bondscnt(i)
          select case(bonds(i,j))
          case(1)
             parbnd(i,j,1)=parbnd(i,j,1)/econv
             parbnd(i,j,2)=parbnd(i,j,2)/kconv
             parbnd(i,j,3)=parbnd(i,j,3)/rconv
          case(2)
             parbnd(i,j,1)=parbnd(i,j,1)/(econv/rconv**2.d0)
             parbnd(i,j,2)=parbnd(i,j,2)/rconv
          case(3)
             parbnd(i,j,1)=parbnd(i,j,1)/(econv/rconv**2.d0)
             parbnd(i,j,2)=parbnd(i,j,2)/rconv
          end select
       end do
    end do

    return

  end subroutine bonds_convert

  subroutine bonds_counts
    !****************************************************************************************
    ! Contagem do número de diedros                                                         *
    !****************************************************************************************

    implicit none

    integer, allocatable :: chk(:,:)

    integer nx,nxx,i,j,k,ni,nj,np,nkb

    !-alocando array

    allocate(chk(nmolec,nbonds))

    !-checando viabilidade de bonds

    nkb=1

    do i=1,nmolec
       do j=1,bondscnt(i)
          select case(bonds(i,j))
          case(1)
             nkb=3
          case(2)
             nkb=2
          case(3)
             nkb=2
          end select
          chk(i,j)=1
          do k=1,nkb
             if(parbnd(i,j,k).eq.0.d0)chk(i,j)=0
          end do
       end do
    end do

    !-calculando bonds intramoleculares

    np=0
    nx=1
    do i=1,nmolec
       nxx=0
       do j=1,ntmolec(i)
          do k=1,bondscnt(i)
             ni=np+molbond(i,k,1)
             nj=np+molbond(i,k,2)
             if(chk(i,k).eq.1)then
                bondij(1,nx)=ni
                bondij(2,nx)=nj
                bondim(nx)=i
                bondib(nx)=k
                nxx=nxx+1
                nx=nx+1
             end if
          end do
          np=np+nxmolec(i)
       end do
       bondsmlc(i)=nxx
    end do

    nbondstp=nx-1

    !-limpando memoria

    deallocate(chk)

    return

  end subroutine bonds_counts

  subroutine bonds_calc(enbond,virbond)
    !****************************************************************************************
    ! - Energia potencial;                                                                  *
    ! - Contribuicao para o virial;                                                         *
    ! - Contribuicao para o stress;                                                         *
    ! - Forças atômicas.                                                                    *
    !****************************************************************************************

    implicit none

    integer i,ia,ib,im,in
    real(8) xvz,yvz,zvz,fr,dr,pot
    real(8) enbond,virbond

    do i=1,nbondstp

       ia=bondij(1,i)
       ib=bondij(2,i)
       im=bondim(i)
       in=bondib(i)

!       call ccpmm(ia,ib,xvz,yvz,zvz)
       call mic(ia,ib,xvz,yvz,zvz)

       !distancia interatomica

       dr=sqrt(xvz**2+yvz**2+zvz**2)

       call bonds_flags(im,in,dr,pot,fr)
       call bonds_force(ia,ib,xvz,yvz,zvz,fr,virbond)

       enbond=enbond+pot

    end do

    return

  end subroutine bonds_calc

  subroutine bonds_flags(im,in,dr,pot,fr)
    !****************************************************************************************
    ! Compontente angular                                                                   *
    !****************************************************************************************

    implicit none

    integer im,in,ii
    real(8) pot,fr,dr,prm(3)

    !-valores iniciais

    pot=0.d0
    fr=0.d0

    !-calculo do gradiente e energia potencial

    select case(bonds(im,in))
    case(1)
       do ii=1,3
          prm(ii)=parbnd(im,in,ii)
       end do
       pot=prm(1)*(exp(-2.d0*prm(2)*(dr-prm(3))) &
            -2.d0*exp(-prm(2)*(dr-prm(3))))

       fr=2.d0*prm(1)*prm(2)*(exp(-2.d0*prm(2)*(dr-prm(3))) &
            -exp(-prm(2)*(dr-prm(3))))/dr
    case(2)
       do ii=1,2
          prm(ii)=parbnd(im,in,ii)
       end do
       pot=0.5d0*prm(1)*(dr-prm(2))**2

       fr=-prm(1)*(dr-prm(2))/dr
    case(3)
       do ii=1,2
          prm(ii)=parbnd(im,in,ii)
       end do
       pot=0.5d0*prm(1)*(dr-prm(2))**2

       fr=-prm(1)*(dr-prm(2))/dr
    end select

    return

  end subroutine bonds_flags

  subroutine bonds_force(i,j,xvz,yvz,zvz,fr,virbond)
    !****************************************************************************************
    ! - Contribuicao para o virial                                                          *
    ! - Contribuicao para o stress                                                          *
    ! - Forcas dos atomos i e j                                                             *
    !****************************************************************************************

    implicit none

    integer i,j
    real(8) xvz,yvz,zvz,fr,virbond

    !-contribuicao para as forcas atomicas

    fax(i)=fax(i)-fr*xvz
    fay(i)=fay(i)-fr*yvz
    faz(i)=faz(i)-fr*zvz

    fax(j)=fax(j)+fr*xvz
    fay(j)=fay(j)+fr*yvz
    faz(j)=faz(j)+fr*zvz

    !-contribuicao para o virial

    virbond=virbond+fr*(xvz**2+yvz**2+zvz**2)

    !-contribuição para o stress

    str(1)=str(1)-(fr*xvz)*xvz
    str(2)=str(2)-(fr*yvz)*yvz
    str(3)=str(3)-(fr*zvz)*zvz
    str(4)=str(4)-(fr*zvz)*yvz
    str(5)=str(5)-(fr*xvz)*zvz
    str(6)=str(6)-(fr*yvz)*xvz

    return

  end subroutine bonds_force

end module bonds_module
