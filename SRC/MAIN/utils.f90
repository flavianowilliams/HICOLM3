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

  use input

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

  subroutine spline_cubic(dh,func,funcl,ai,bi,ci,di)
    !************************************************************************************
    ! Interpolacao cubica - spline                                                      *
    !************************************************************************************

    implicit none

    integer i,INFO
    integer IPIV(3)
    real(8) dh(3),func(3),funcl(3),ai(3),bi(3),ci(3),di(3)
    real(8) spa(3,3),spb(3,1)

  !-montando matriz A do sistema Ax=B

    spa(1,1)=2.d0*dh(1)
    spa(1,2)=dh(1)
    spa(1,3)=0.d0

    spa(2,1)=dh(1)
    spa(2,2)=2.d0*(dh(1)+dh(2))
    spa(2,3)=dh(2)

    spa(3,1)=0.d0
    spa(3,2)=dh(2)
    spa(3,3)=2.d0*dh(2)

  !-montando matriz B do sistema Ax=B

    spb(1,1)=3.d0*(func(2)-func(1))/dh(1)-3.d0*funcl(1)
    spb(2,1)=3.d0*(func(3)-func(2))/dh(2)-3.d0*(func(2)-func(1))/dh(1)
    spb(3,1)=-3.d0*(func(3)-func(2))/dh(2)+3.d0*funcl(3)

  !-fatorizando matriz principal

    call dgetrf(3,3,spa,3,IPIV,INFO)

    if(INFO.eq.0)then
       call dgetrs('N',3,1,spa,3,IPIV,spb,3,INFO)
    else
       write(*,*)'spline_cubic: Error!!!!',INFO
       stop
    end if

    do i=1,3
       ci(i)=spb(i,1)
    end do

    do i=1,2
       ai(i)=func(i)
       bi(i)=(func(i+1)-func(i))/dh(i)-dh(i)*(2.d0*ci(i)+ci(i+1))/3.d0
       di(i)=(ci(i+1)-ci(i))/(3.d0*dh(i))
    end do

    return

  end subroutine spline_cubic

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
