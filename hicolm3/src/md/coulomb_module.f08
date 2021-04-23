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
  !*******************************************************************************************
  !*******************************************************************************************

  implicit none

!  integer i,j,k,l

  private
  public :: coulomb

  type :: coulomb
     real(8), private      :: encoul
     real(8), private      :: vircoul
     real(8), private      :: rcutoff
     real(8), private      :: kconv
     real(8), private      :: pi
     real(8), private      :: force
     character(4), private :: coulop
   contains
     procedure :: set_coulomb
     procedure :: coulomb_prepare
     procedure :: set_encoul
     procedure :: get_encoul
     procedure :: set_vircoul
     procedure :: get_vircoul
     procedure :: get_force
  end type coulomb

contains

  subroutine coulomb_prepare(this,coulop,kconv,rcutoff,pi)
    implicit none
    class(coulomb), intent(inout) :: this
    real(8), intent(in)           :: rcutoff,kconv,pi
    character(4), intent(in)      :: coulop
    this%coulop=coulop
    this%kconv=kconv
    this%rcutoff=rcutoff
    this%pi=pi
  end subroutine coulomb_prepare

  subroutine set_coulomb(this,dr,qi,qj)
    implicit none
    class(coulomb), intent(inout) :: this
    real(8), intent(in)           :: dr,qi,qj
    real(8)                       :: alcoul
    select case(this%coulop)
    case('coul')
       this%encoul=qi*qj/dr
       this%force=-qi*qj/dr**2
       this%force=-this%force/dr
    case('fscs')
       alcoul=1.d-1/this%kconv
       this%encoul=qi*qj*(erfc(alcoul*dr)/dr-erfc(alcoul*this%rcutoff)/this%rcutoff &
            +(erfc(alcoul*this%rcutoff)/this%rcutoff**2+(2.d0*alcoul) &
            *exp(-(alcoul*this%rcutoff)**2)/(sqrt(this%pi)*this%rcutoff))*(dr-this%rcutoff))
       this%force=-qi*qj*(erfc(alcoul*dr)/dr**2+(2.d0*alcoul) &
            *exp(-(alcoul*dr)**2)/(sqrt(this%pi)*dr) &
            -(erfc(alcoul*this%rcutoff)/this%rcutoff**2 &
            +(2.d0*alcoul)*exp(-(alcoul*this%rcutoff)**2)/(sqrt(this%pi)*this%rcutoff)))
       this%force=-this%force/dr
    end select
  end subroutine set_coulomb

  subroutine set_encoul(this,encoul)
    implicit none
    class(coulomb), intent(inout) :: this
    real(8), intent(in)           :: encoul
    this%encoul=encoul
  end subroutine set_encoul

  double precision function get_encoul(this)
    implicit none
    class(coulomb), intent(in) :: this
    get_encoul=this%encoul
  end function get_encoul

  subroutine set_vircoul(this,fr,dr)
    implicit none
    class(coulomb), intent(inout) :: this
    real(8), intent(in)           :: fr,dr
    this%vircoul=-fr*dr**2
  end subroutine set_vircoul

  double precision function get_vircoul(this)
    implicit none
    class(coulomb), intent(in) :: this
    get_vircoul=this%vircoul
  end function get_vircoul

  double precision function get_force(this)
    implicit none
    class(coulomb), intent(in) :: this
    get_force=this%force
  end function get_force

end module coulomb_module
