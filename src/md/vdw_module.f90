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
module vdw_module
  !******************************************************************************************
  ! Contribuicao de Van der Waals para o campo de força:                                    *
  ! - Energia potencial;                                                                    *
  ! - Forças atômicas;                                                                      *
  ! - Virial;                                                                               *
  ! - Stress.                                                                               *
  !                                                                                         *
  ! Flaviano Williams Fernandes, 05 de dezembro de 2018                                     *
  !******************************************************************************************

  use input
  use estrutura
  use alloc_arrays
  use dihedral_module

  integer nvdwstp
  real(8) envdw_corr,virvdw_corr

  ! integer nintp
!  real(8) interp_v(50,50,1000000),interp_f(50,50,1000000),interp_dr(1000000)
!  save nintp

  save nvdwstp,envdw_corr,virvdw_corr

contains

  subroutine vdw_prepare
    !****************************************************************************************
    ! Preparando parametros para o calculo de Van der Waals:                                *
    ! Conversao de unidades ---> a.u.                                                       *
    !****************************************************************************************
    implicit none

    integer, allocatable :: chk(:,:)

    integer i,j,k,ix,nx,ii,jj,ixx,nkb

!    real(8) rinit,drintp,pot,fr

    allocate(chk(spctot,spctot))

    !-convertendo parametros de Van der Waals

    do i=1,ntpmax
       do j=1,ntpmax
          select case(vdw(i,j))
          case(1)
             parvdw(i,j,1)=parvdw(i,j,1)/econv
             parvdw(i,j,2)=parvdw(i,j,2)/kconv
             parvdw(i,j,3)=parvdw(i,j,3)/rconv
          case(2)
             parvdw(i,j,1)=parvdw(i,j,1)/econv
             parvdw(i,j,2)=parvdw(i,j,2)/rconv
          case(3)
             parvdw(i,j,1)=parvdw(i,j,1)/econv
             parvdw(i,j,2)=parvdw(i,j,2)/rconv
          end select
       end do
    end do

    !-checando viabilidade de interacoes intermoleculares

    nkb=1

    do i=1,spctot
       do j=1,spctot
          select case(vdw(i,j))
          case(1)
             nkb=3
          case(2)
             nkb=2
          case(3)
             nkb=2
          end select
          chk(i,j)=1
          do k=1,nkb
             if(parvdw(i,j,k).eq.0.d0)chk(i,j)=0
          end do
       end do
    end do

    !-calculando qde de Van der Waals

    nvdw=0
    do i=1,spctot
       do j=i,spctot
          if(chk(i,j).eq.1)nvdw=nvdw+1
       end do
    end do

    ix=1
    nx=1
    nvdwstp=0
    do i=1,moltot-1
       do ii=1,nzmolec(i)
          ixx=nx+nzmolec(i)
          do j=i+1,moltot
             do jj=1,nzmolec(j)
                if(chk(atp(ix),atp(ixx)).eq.1)nvdwstp=nvdwstp+1
                ixx=ixx+1
             end do
          end do
          ix=ix+1
       end do
       nx=nx+nzmolec(i)
    end do

    deallocate(chk)

    return

  end subroutine vdw_prepare

  subroutine vdw_calc(envdw,virvdw,ni,nj,xvz,yvz,zvz)
    !****************************************************************************************
    ! - Energia potencial;                                                                  *
    ! - Contribuicao para o virial;                                                         *
    ! - Contribuicao para o stress;                                                         *
    ! - Forças atômicas.                                                                    *
    !****************************************************************************************

    implicit none

    integer ni,nj
    real(8) xvz,yvz,zvz,fr,pot
    real(8) envdw,virvdw

    call vdw_flags(atp(ni),atp(nj),xvz,yvz,zvz,pot,fr)
    call vdw_force(ni,nj,xvz,yvz,zvz,fr)

    virvdw=virvdw+fr*(xvz**2+yvz**2+zvz**2)
    envdw=envdw+pot

    return

  end subroutine vdw_calc

  subroutine vdw_14sf(envdw,virvdw)

    implicit none

    integer i,ni,nj,nk
    real(8) pot,fr,xvz,yvz,zvz
    real(8) envdw,virvdw

    do i=1,ntorsstp

       ni=torsijkn(1,i)
       nj=torsijkn(4,i)
       nk=torsim(i)

       call mic(ni,nj,xvz,yvz,zvz)

       call vdw_flags(atp(ni),atp(nj),xvz,yvz,zvz,pot,fr)

       pot=pot*sf_vdw(nk)
       fr=fr*sf_vdw(nk)

       call vdw_force(ni,nj,xvz,yvz,zvz,fr)

       virvdw=virvdw+fr*(xvz**2+yvz**2+zvz**2)
       envdw=envdw+pot

    end do

    return

  end subroutine vdw_14sf

  subroutine vdw_flags(i,j,xvz,yvz,zvz,pot,fr)
    !****************************************************************************************
    ! Compontente angular                                                                   *
    !****************************************************************************************

    implicit none

    integer i,j,ii,ptrm
    real(8) pot,fr,xvz,yvz,zvz,dr,prm(3)

    fr=0.d0
    pot=0.d0

    dr=sqrt(xvz**2+yvz**2+zvz**2)

    ptrm=vdw(i,j)

    select case(ptrm)
    case(1)

       do ii=1,3
          prm(ii)=parvdw(i,j,ii)
       end do

       pot=prm(1)*(exp(-2.d0*prm(2)*(dr-prm(3))) &
            -2.d0*exp(-prm(2)*(dr-prm(3))))

       fr=2.d0*prm(1)*prm(2)*(exp(-2.d0*prm(2)*(dr-prm(3))) &
            -exp(-prm(2)*(dr-prm(3))))/dr

    case(2)

       do ii=1,2
          prm(ii)=parvdw(i,j,ii)
       end do

       pot=4.d0*prm(1)*((prm(2)/dr)**12-(prm(2)/dr)**6)

       fr=24.d0*prm(1)*(2.d0*(prm(2)/dr)**12-(prm(2)/dr)**6)/dr**2

    case(3)

       do ii=1,2
          prm(ii)=parvdw(i,j,ii)
       end do

       !-versao LJ segundo AMBER

       pot=prm(1)*((prm(2)/dr)**12-2.d0*(prm(2)/dr)**6)

       fr=12.d0*prm(1)*((prm(2)/dr)**12-(prm(2)/dr)**6)/dr**2

       !-versão LJ normal

!       pot=4.d0*prm(1)*((prm(2)/dr)**12-(prm(2)/dr)**6)

!       fr=24.d0*prm(1)*(2.d0*(prm(2)/dr)**12-(prm(2)/dr)**6)/dr**2

    end select

    return

  end subroutine vdw_flags

  subroutine vdw_force(i,j,xvz,yvz,zvz,fr)
    !****************************************************************************************
    ! - Contribuicao para o virial                                                          *
    ! - Contribuicao para o stress                                                          *
    ! - Forcas dos atomos i, j, k e n                                                       *
    !****************************************************************************************
    implicit none

    integer i,j
    real(8) xvz,yvz,zvz,fr

    !-contribuicao para a forca atomica

    fax(i)=fax(i)-fr*xvz
    fay(i)=fay(i)-fr*yvz
    faz(i)=faz(i)-fr*zvz
    fax(j)=fax(j)+fr*xvz
    fay(j)=fay(j)+fr*yvz
    faz(j)=faz(j)+fr*zvz

    !-contribuição para o stress

    str(1)=str(1)+(fr*xvz)*xvz
    str(2)=str(2)+(fr*yvz)*yvz
    str(3)=str(3)+(fr*zvz)*zvz
    str(4)=str(4)+(fr*zvz)*yvz
    str(5)=str(5)+(fr*xvz)*zvz
    str(6)=str(6)+(fr*yvz)*xvz

    return

  end subroutine vdw_force

  subroutine vdw_corr()
    !****************************************************************************************
    ! - Correcao de Van der Waals para energia e virial                                     *
    !****************************************************************************************

    implicit none

    integer ii,i,j,ptrm
    real(8) prm(5),es,vs

    envdw_corr=0.d0
    virvdw_corr=0.d0

    do i=1,spctot
       do j=i,spctot
          ptrm=vdw(i,j)
          select case(ptrm)
          case(1)
             do ii=1,5
                prm(ii)=parvdw(i,j,ii)
             end do
             envdw_corr=0.d0
          case(2)
             do ii=1,2
                prm(ii)=parvdw(i,j,ii)
             end do
             es=4.d0*prm(1)*(prm(2)**12-3.d0*(rcutoff*prm(2))**6)/(9.d0*rcutoff**9)
             vs=24.d0*prm(1)*&
                  (2.d0*prm(2)**12-3.d0*(rcutoff*prm(2))**6)/(9.d0*rcutoff**9)
             envdw_corr=envdw_corr+es*natnp(i)*natnp(j)
             virvdw_corr=virvdw_corr+vs*natnp(i)*natnp(j)
          case(3)
             do ii=1,2
                prm(ii)=parvdw(i,j,ii)
             end do
             es=prm(1)*&
                  (prm(2)**12-6.d0*(rcutoff*prm(2))**6)/(9.d0*rcutoff**9)
             vs=12.d0*prm(1)*&
                  (prm(2)**12-3.d0*(rcutoff*prm(2))**6)/(9.d0*rcutoff**9)
             envdw_corr=envdw_corr+es*natnp(i)*natnp(j)
             virvdw_corr=virvdw_corr+vs*natnp(i)*natnp(j)
          end select
          envdw_corr=envdw_corr+es
          virvdw_corr=virvdw_corr+vs
       end do
    end do

    envdw_corr=envdw_corr*2*pi/volume
    virvdw_corr=virvdw_corr*2*pi/volume

    return

  end subroutine vdw_corr

end module vdw_module
