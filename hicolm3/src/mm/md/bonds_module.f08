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
module bonds_module
  !*******************************************************************************************
  !*******************************************************************************************

  implicit none

  private
  public :: bonds

  type :: bonds
     real(8), private :: enbond
     real(8), private :: virbond
     real(8), private :: force
   contains
     procedure :: set_bonds
     procedure :: set_enbond
     procedure :: get_enbond
     procedure :: set_virbond
     procedure :: get_virbond
     procedure :: get_force
  end type bonds

contains

  subroutine set_bonds(this,dr,prm,ptrm)
    implicit none
    class(bonds), intent(inout) :: this
    character(6), intent(in)    :: ptrm
    real(8), intent(in)         :: dr,prm(2)
    select case(ptrm)
    case('charmm')
       this%enbond=prm(1)*(dr-prm(2))**2
       this%force=-2.d0*prm(1)*(dr-prm(2))/dr
    case('harm')
       this%enbond=0.5d0*prm(1)*(dr-prm(2))**2
       this%force=-prm(1)*(dr-prm(2))/dr
    end select
  end subroutine set_bonds

  subroutine set_enbond(this,enbond)
    implicit none
    class(bonds), intent(inout) :: this
    real(8), intent(in)         :: enbond
    this%enbond=enbond
  end subroutine set_enbond

  double precision function get_enbond(this)
    implicit none
    class(bonds), intent(in) :: this
    get_enbond=this%enbond
  end function get_enbond

  subroutine set_virbond(this,fr,dr)
    implicit none
    class(bonds), intent(inout) :: this
    real(8), intent(in)         :: fr,dr
    this%virbond=-fr*dr**2
  end subroutine set_virbond

  double precision function get_virbond(this)
    implicit none
    class(bonds), intent(in) :: this
    get_virbond=this%virbond
  end function get_virbond

  double precision function get_force(this)
    implicit none
    class(bonds), intent(in) :: this
    get_force=this%force
  end function get_force

end module bonds_module
