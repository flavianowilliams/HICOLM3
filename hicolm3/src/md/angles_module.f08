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
module angles_module
  !*******************************************************************************************
  !*******************************************************************************************

  implicit none

  private
  public :: angles

  type :: angles
     real(8), private :: enbend
     real(8), private :: virbend
     real(8), private :: force
   contains
     procedure :: set_angles
     procedure :: set_enbend
     procedure :: get_enbend
     procedure :: set_virbend
     procedure :: get_virbend
     procedure :: get_force
  end type angles

contains

  subroutine set_angles(this,theta,prm,ptrm)
    implicit none
    class(angles), intent(inout) :: this
    character(6), intent(in)     :: ptrm
    real(8), intent(in)          :: theta,prm(2)
    select case(ptrm)
    case('charmm')
       this%enbend=prm(1)*(theta-prm(2))**2
       this%force=2.d0*prm(1)*(theta-prm(2))
    case('harm')
       this%enbend=0.5d0*prm(1)*(theta-prm(2))**2
       this%force=prm(1)*(theta-prm(2))
    end select
  end subroutine set_angles

  subroutine set_enbend(this,enbend)
    implicit none
    class(angles), intent(inout) :: this
    real(8), intent(in)          :: enbend
    this%enbend=enbend
  end subroutine set_enbend

  double precision function get_enbend(this)
    implicit none
    class(angles), intent(in) :: this
    get_enbend=this%enbend
  end function get_enbend

  subroutine set_virbend(this,fbj,fbk,drij,drik)
    implicit none
    class(angles), intent(inout) :: this
    real(8), intent(in)          :: fbj(3),fbk(3),drij(3),drik(3)
    this%virbend=-((fbj(1)*drij(1)+fbj(2)*drij(2)+fbj(3)*drij(3))&
         +(fbk(1)*drik(1)+fbk(2)*drik(2)+fbk(3)*drik(3)))
  end subroutine set_virbend

  double precision function get_virbend(this)
    implicit none
    class(angles), intent(in) :: this
    get_virbend=this%virbend
  end function get_virbend

  double precision function get_force(this)
    implicit none
    class(angles), intent(in) :: this
    get_force=this%force
  end function get_force

end module angles_module
