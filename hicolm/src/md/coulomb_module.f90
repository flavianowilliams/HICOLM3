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
  use dihedral_module

  integer ncoulstp

  save ncoulstp

contains

  subroutine coulomb_prepare

    implicit none

    integer i,j

    !-convertendo unidades de medida

    do i=1,nmolec
       do j=1,nxmolec(i)
          qatmolec(i,j)=qatmolec(i,j)/elconv
       end do
    end do

    ncoulstp=0
    do i=1,natom
       do j=i+1,natom
          if(abs(qat(i)*qat(j)).gt.1.d-8)ncoulstp=ncoulstp+1
       end do
    end do

    return

  end subroutine coulomb_prepare

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

    call coulomb_flags(ni,nj,xvz,yvz,zvz,pot,fr)
    call coulomb_force(ni,nj,xvz,yvz,zvz,fr)

    vircoul=vircoul+fr*(xvz**2+yvz**2+zvz**2)
    encoul=encoul+pot

    return

  end subroutine coulomb_calc

  subroutine coulomb_sf(ni,nj,nk,xvz,yvz,zvz,encoul,vircoul)

    implicit none

    integer i,ni,nj,nk
    real(8) pot,fr,xvz,yvz,zvz
    real(8) encoul,vircoul

       call coulomb_flags(ni,nj,xvz,yvz,zvz,pot,fr)

       pot=pot*sf_coul(nk)
       fr=fr*sf_coul(nk)

       call coulomb_force(ni,nj,xvz,yvz,zvz,fr)

       vircoul=vircoul+fr*(xvz**2+yvz**2+zvz**2)
       encoul=encoul+pot

    return

  end subroutine coulomb_sf

  subroutine coulomb_flags(i,j,xvz,yvz,zvz,pot,fr)
    !****************************************************************************************
    ! Compontente angular                                                                   *
    !****************************************************************************************

    implicit none

    integer i,j
    real(8) pot,fr,xvz,yvz,zvz,dr,alcoul,xij,lambda

    !-atribuindo valores iniciais

    pot=0.d0
    fr=0.d0

    !-calculo do gradiente e potencial

    dr=max(1.d-4,sqrt(xvz**2+yvz**2+zvz**2))

    select case(coulop)

    case('coul')

       pot=qat(i)*qat(j)/dr

       fr=-qat(i)*qat(j)/dr**2

       fr=-fr/dr ! -(1/r)*dU/dr

    case('fscs')

       alcoul=1.d-1/kconv

       pot=qat(i)*qat(j)*(erfc(alcoul*dr)/dr-erfc(alcoul*rcutoff)/rcutoff &
            +(erfc(alcoul*rcutoff)/rcutoff**2+(2.d0*alcoul) &
            *exp(-(alcoul*rcutoff)**2)/(sqrt(pi)*rcutoff))*(dr-rcutoff))

       fr=-qat(i)*qat(j)*(erfc(alcoul*dr)/dr**2+(2.d0*alcoul) &
            *exp(-(alcoul*dr)**2)/(sqrt(pi)*dr) &
            -(erfc(alcoul*rcutoff)/rcutoff**2 &
            +(2.d0*alcoul)*exp(-(alcoul*rcutoff)**2)/(sqrt(pi)*rcutoff)))

       fr=-fr/dr ! -(1/r)*dU/dr

    case('escl')

       lambda=lambdafi

       xij=sqrt(2.d0*(1.d0-lambda)**2+dr**2)

       pot=qat(i)*qat(j)*lambda*(1.d0/xij+xij/rcutoff**2-2.d0/rcutoff)

       fr=qat(i)*qat(j)*lambda*(1.d0/xij**2-1.d0/rcutoff**2)/xij

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
