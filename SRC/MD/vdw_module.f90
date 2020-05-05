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
  use neighbour_list

  real(8) nvdwstp,envdw_corr,virvdw_corr

  save nvdwstp,envdw_corr,virvdw_corr

contains

  subroutine vdw_convert
    !****************************************************************************************
    ! Conversao de unidades de medida:                                                      *
    ! Unidades de entrada ---> a.u.                                                         *
    !****************************************************************************************

    implicit none

    integer nx,i,j,ii,jj,ni,nj

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
          end select
       end do
    end do

    !-calculando qde de Van der Waals e coulomb

    nx=0
    do i=1,moltot
       do ii=1,nzmolec(i)
          ni=(i-1)*nzmolec(i)+ii
          do j=i+1,moltot
             do jj=1,nzmolec(j)
                nj=(j-1)*nzmolec(j)+jj
                nvdwstp=nvdwstp+1
             end do
          end do
       end do
    end do

    return

  end subroutine vdw_convert

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

    do i=1,nvdw
       do j=i,nvdw
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
             vs=-8.d0*prm(1)*(2.d0*(prm(2)/rcutoff)**6/3.d0-1.d0)*prm(2)**6/rcutoff**3
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
