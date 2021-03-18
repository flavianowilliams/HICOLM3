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
module utils
  !****************************************************************************************
  ! Utilitários                                                                           *
  !                                                                                       *
  ! Flaviano Williams Fernandes, 05 de dezembro de 2018                                   *
  !****************************************************************************************

!  use input

contains

  double precision function heaviside(t)
    !************************************************************************************
    ! Funcao de truncamento Heaviside                                                   *
    !************************************************************************************

    implicit none

    real(8) t

    heaviside=0.5d0*(1.d0+t/max(abs(t),1.d-8))

    return

  end function heaviside

  integer function kronij(i,j)
    !************************************************************************************
    ! Delta de Kronecker                                                                *
    !************************************************************************************

    implicit none

    integer i,j

    kronij=int((float(i+j)-abs(i-j))/(float(i+j)+abs(i-j)))

    return

  end function kronij

end module utils
