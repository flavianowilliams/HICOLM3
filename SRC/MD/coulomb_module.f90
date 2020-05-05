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
module coulomb_module
  !******************************************************************************************
  ! Contribuicao eletrostatica para o campo de força:                                       *
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
  use neighbour_list

contains

  subroutine coulomb_convert
    !****************************************************************************************
    ! Conversao de unidades de medida:                                                      *
    ! Unidades de entrada ---> a.u.                                                         *
    !****************************************************************************************

    implicit none

    integer i

    !-convertendo Van der Waals

    do i=1,ntpmax
       parcoul(i,1)=parcoul(i,1)/elconv
    end do

    return

  end subroutine coulomb_convert

  subroutine coulomb_calc(encoul,vircoul,ni,nj,xvz,yvz,zvz)
    !****************************************************************************************
    ! - Energia potencial;                                                                  *
    ! - Contribuicao para o virial;                                                         *
    ! - Contribuicao para o stress;                                                         *
    ! - Forças atômicas.                                                                    *
    !****************************************************************************************

    implicit none

    integer ni,nj
    real(8) xvz,yvz,zvz,fr,pot
    real(8) encoul,vircoul

    call coulomb_flags(atp(ni),atp(nj),xvz,yvz,zvz,pot,fr)
    call coulomb_force(ni,nj,xvz,yvz,zvz,fr)

    vircoul=vircoul+fr*(xvz**2+yvz**2+zvz**2)
    encoul=encoul+pot

    return

  end subroutine coulomb_calc

  subroutine coulomb_flags(i,j,xvz,yvz,zvz,pot,fr)
    !****************************************************************************************
    ! Compontente angular                                                                   *
    !****************************************************************************************

    implicit none

    integer i,j
    real(8) pot,fr,xvz,yvz,zvz,dr,prm(2),alcoul,xij,lambda

    !-atribuindo valores iniciais

    pot=0.d0
    fr=0.d0

    !-calculo do gradiente e potencial

    dr=sqrt(xvz**2+yvz**2+zvz**2)

    prm(1)=parcoul(i,1)
    prm(2)=parcoul(j,1)
    !
    select case(coulop)
    case('coul')

       alcoul=1.d-1/kconv

       pot=prm(1)*prm(2)*(erfc(alcoul*dr)/dr-erfc(alcoul*rcutoff)/rcutoff &
            +(erfc(alcoul*rcutoff)/rcutoff**2+(2.d0*alcoul) &
            *exp(-(alcoul*rcutoff)**2)/(sqrt(pi)*rcutoff))*(dr-rcutoff))

       fr=-prm(1)*prm(2)*(erfc(alcoul*dr)/dr**2+(2.d0*alcoul) &
            *exp(-(alcoul*dr)**2)/(sqrt(pi)*dr) &
            -(erfc(alcoul*rcutoff)/rcutoff**2 &
            +(2.d0*alcoul)*exp(-(alcoul*rcutoff)**2)/(sqrt(pi)*rcutoff)))

       fr=-fr/dr ! -(1/r)*dU/dr

    case('escl')

       lambda=lambdafi

       xij=sqrt(2.d0*(1.d0-lambda)**2+dr**2)

       pot=prm(1)*prm(2)*lambda*(1.d0/xij+xij/rcutoff**2-2.d0/rcutoff)

       fr=prm(1)*prm(2)*lambda*(1.d0/xij**2-1.d0/rcutoff**2)/xij

    end select
!
    return
!
  end subroutine coulomb_flags

  subroutine coulomb_force(i,j,xvz,yvz,zvz,fr)
    !****************************************************************************************
    ! - Contribuicao para o virial                                                          *
    ! - Contribuicao para o stress                                                          *
    ! - Forcas dos atomos i, j, k e n                                                       *
    !****************************************************************************************

    implicit none

    integer i,j
    real(8) xvz,yvz,zvz,fr

    !-calculo da contribuicao da forca atomica

    fax(i)=fax(i)-fr*xvz
    fay(i)=fay(i)-fr*yvz
    faz(i)=faz(i)-fr*zvz
    fax(j)=fax(j)+fr*xvz
    fay(j)=fay(j)+fr*yvz
    faz(j)=faz(j)+fr*zvz

    !-calculo da contribuição do stress

    str(1)=str(1)+(fr*xvz)*xvz
    str(2)=str(2)+(fr*yvz)*yvz
    str(3)=str(3)+(fr*zvz)*zvz
    str(4)=str(4)+(fr*zvz)*yvz
    str(5)=str(5)+(fr*xvz)*zvz
    str(6)=str(6)+(fr*yvz)*xvz

    return

   end subroutine coulomb_force

end module coulomb_module
