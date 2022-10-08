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
module idihedrals_module
  !*******************************************************************************************
  !*******************************************************************************************

  implicit none

  private
  public :: idihedrals

  type :: idihedrals
     real(8), private :: enitors
     real(8), private :: viritors
     real(8), private :: force
   contains
     procedure :: set_idihedrals
     procedure :: set_enitors
     procedure :: get_enitors
     procedure :: set_viritors
     procedure :: get_viritors
     procedure :: get_force
  end type idihedrals

contains

  subroutine set_idihedrals(this,phi,prm,ptrm)
    implicit none
    class(idihedrals), intent(inout) :: this
    character(7), intent(in)        :: ptrm
    real(8), intent(in)             :: phi,prm(3)
    select case(ptrm)
    case('charmm2')
       this%enitors=prm(1)*(1.0d0+cos(prm(2)*phi-prm(3)))
       this%force=-prm(1)*sin(prm(2)*phi-prm(3))
    case('charmm')
       this%enitors=prm(1)*(phi-prm(2))**2
       this%force=2.0d0*prm(1)*(phi-prm(2))
    case('harm')
       this%enitors=0.5d0*prm(1)*(phi-prm(2))**2
       this%force=prm(1)*(phi-prm(2))
    end select
  end subroutine set_idihedrals

  subroutine set_enitors(this,enitors)
    implicit none
    class(idihedrals), intent(inout) :: this
    real(8), intent(in)             :: enitors
    this%enitors=enitors
  end subroutine set_enitors

  double precision function get_enitors(this)
    implicit none
    class(idihedrals), intent(in) :: this
    get_enitors=this%enitors
  end function get_enitors

  subroutine set_viritors(this,fbi,fbj,fbk,fbl,ri,rj,rk,rl)
    implicit none
    class(idihedrals), intent(inout) :: this
    integer                         :: i
    real(8), intent(in)             :: fbi(3),fbj(3),fbk(3),fbl(3)
    real(8), intent(in)             :: ri(3),rj(3),rk(3),rl(3)
    real(8)                         :: viritors
    viritors=0.d0
    do i=1,3
       viritors=viritors-(fbi(i)*ri(i)+fbj(i)*rj(i)+fbk(i)*rk(i)+fbl(i)*rl(i))
    end do
    this%viritors=viritors
  end subroutine set_viritors

  double precision function get_viritors(this)
    implicit none
    class(idihedrals), intent(in) :: this
    get_viritors=this%viritors
  end function get_viritors

  double precision function get_force(this)
    implicit none
    class(idihedrals), intent(in) :: this
    get_force=this%force
  end function get_force

end module idihedrals_module
