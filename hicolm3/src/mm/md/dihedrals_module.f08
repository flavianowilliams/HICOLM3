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
module dihedrals_module
  !*******************************************************************************************
  !*******************************************************************************************

  implicit none

  private
  public :: dihedrals

  type :: dihedrals
     real(8), private :: entors
     real(8), private :: virtors
     real(8), private :: force
   contains
     procedure :: set_dihedrals
     procedure :: set_entors
     procedure :: get_entors
     procedure :: set_virtors
     procedure :: get_virtors
     procedure :: get_force
  end type dihedrals

contains

  subroutine set_dihedrals(this,phi,prm,ptrm)
    implicit none
    class(dihedrals), intent(inout) :: this
    character(7), intent(in)        :: ptrm
    real(8), intent(in)             :: phi,prm(3)
    select case(ptrm)
    case('charmm')
       this%entors=prm(1)*(1.0d0+cos(prm(2)*phi-prm(3)))
       this%force=-prm(1)*sin(prm(2)*phi-prm(3))
    case('charmm2')
       this%entors=prm(1)*(phi-prm(2))**2
       this%force=2.0d0*prm(1)*(phi-prm(2))
    case('harm')
       this%entors=0.5d0*prm(1)*(phi-prm(2))**2
       this%force=prm(1)*(phi-prm(2))
    case('opls')
       this%entors=0.5d0*(prm(1)*(1.d0+cos(1.d0*phi))&
            +prm(2)*(1.d0-cos(2.d0*phi))&
            +prm(3)*(1.d0+cos(3.d0*phi)))
       this%force=0.5d0*(-prm(1)*sin(1.d0*phi)&
            +prm(2)*sin(2.d0*phi)&
            -prm(3)*sin(3.d0*phi))
    end select
  end subroutine set_dihedrals

  subroutine set_entors(this,entors)
    implicit none
    class(dihedrals), intent(inout) :: this
    real(8), intent(in)             :: entors
    this%entors=entors
  end subroutine set_entors

  double precision function get_entors(this)
    implicit none
    class(dihedrals), intent(in) :: this
    get_entors=this%entors
  end function get_entors

  subroutine set_virtors(this,fbi,fbj,fbk,fbl,ri,rj,rk,rl)
    implicit none
    class(dihedrals), intent(inout) :: this
    integer                         :: i
    real(8), intent(in)             :: fbi(3),fbj(3),fbk(3),fbl(3),ri(3),rj(3),rk(3),rl(3)
    real(8)                         :: virtors
    virtors=0.d0
    do i=1,3
       virtors=virtors-(fbi(i)*ri(i)+fbj(i)*rj(i)+fbk(i)*rk(i)+fbl(i)*rl(i))
    end do
    this%virtors=virtors
  end subroutine set_virtors

  double precision function get_virtors(this)
    implicit none
    class(dihedrals), intent(in) :: this
    get_virtors=this%virtors
  end function get_virtors

  double precision function get_force(this)
    implicit none
    class(dihedrals), intent(in) :: this
    get_force=this%force
  end function get_force

end module dihedrals_module
