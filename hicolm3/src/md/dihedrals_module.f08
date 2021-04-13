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
    character(6), intent(in)        :: ptrm
    real(8), intent(in)             :: phi,prm(3)
    select case(ptrm)
    case('charmm')
       this%entors=prm(1)*(1.0d0+cos(prm(2)*phi-prm(3)))
       this%force=-prm(1)*sin(prm(2)*phi-prm(3))
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
    class(dihedrals), intent(inout) :: this
    get_entors=this%entors
  end function get_entors

  subroutine set_virtors(this,fbj,fbk,drij,drik)
    implicit none
    class(dihedrals), intent(inout) :: this
    real(8), intent(in)          :: fbj(3),fbk(3),drij(3),drik(3)
    this%virtors=-((fbj(1)*drij(1)+fbj(2)*drij(2)+fbj(3)*drij(3))&
         +(fbk(1)*drik(1)+fbk(2)*drik(2)+fbk(3)*drik(3)))
  end subroutine set_virtors

  double precision function get_virtors(this)
    implicit none
    class(dihedrals), intent(inout) :: this
    get_virtors=this%virtors
  end function get_virtors

  double precision function get_force(this)
    implicit none
    class(dihedrals), intent(inout) :: this
    get_force=this%force
  end function get_force

end module dihedrals_module
