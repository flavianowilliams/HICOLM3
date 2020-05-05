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
