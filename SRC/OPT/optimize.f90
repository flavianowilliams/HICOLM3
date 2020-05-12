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
module optimize

  use input
  use estrutura
  use force_field

contains

  subroutine opt

    implicit none

    integer i
    real(8) encoul,enbond,enbend,entors,envdw
    real(8) virvdw,virbond,virbend,virtors,vircoul

    !-valores iniciais

    !-energia

    encoul=0.d0   !coulombiano
    enbond=0.d0   !estiramento
    enbend=0.d0   !deformacao
    entors=0.d0   !torção
    envdw=0.d0    !Van der waals

    !-virial

    virvdw=0.d0   !Van der Waals
    virbond=0.d0  !estiramento
    virbend=0.d0  !deformacao
    virtors=0.d0  !torção
    vircoul=0.d0  !coulombiano

    !-forcas atomicas

    do i=1,natom
       fax(i)=0.d0
       fay(i)=0.d0
       faz(i)=0.d0
    end do

    !-stress

    do i=1,6
       str(i)=0.d0
    end do

    !-preparando Campo de Força

    call ff_prepare

    !-calculando contribuição intramolecular

    call ff_modules_intra&
         (enbond,enbend,entors,envdw,encoul,virbond,virbend,virtors,virvdw,vircoul)

    return

  end subroutine opt

end module optimize
