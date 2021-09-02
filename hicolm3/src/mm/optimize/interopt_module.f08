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
module interopt_module
  !*******************************************************************************************
  !*******************************************************************************************

  use neighbourlist_module

  implicit none

  private
  public :: interopt

  type, extends(neighbourlist) :: interopt
     real(8), private :: d1bond
     real(8), private :: d2bond
     real(8), private :: d1bend
     real(8), private :: d2bend
   contains
     procedure :: set_bondopt
     procedure :: set_angleopt
     procedure :: set_vanderwaals
     procedure :: set_coulomb
     procedure :: get_d1bond
     procedure :: get_d2bond
     procedure :: get_d1bend
     procedure :: get_d2bend
  end type interopt

contains

  subroutine set_bondopt(this,dr,prm,ptrm,en)
    implicit none
    class(interopt), intent(inout) :: this
    character(6), intent(in)       :: ptrm
    real(8), intent(in)            :: dr,prm(3)
    real(8), intent(out)           :: en
    select case(ptrm)
    case('charmm')
       en=prm(1)*(dr-prm(2))**2
       this%d1bond=2.d0*prm(1)*(dr-prm(2))
       this%d2bond=2.d0*prm(1)
    case('harm')
       en=0.5d0*prm(1)*(dr-prm(2))**2
       this%d1bond=prm(1)*(dr-prm(2))
       this%d2bond=prm(1)
    end select
  end subroutine set_bondopt

  subroutine set_angleopt(this,theta,prm,ptrm,en)
    implicit none
    class(interopt), intent(inout) :: this
    character(6), intent(in)       :: ptrm
    real(8), intent(in)            :: theta,prm(3)
    real(8), intent(out)           :: en
    select case(ptrm)
    case('charmm')
       en=prm(1)*(theta-prm(2))**2
       this%d1bend=2.d0*prm(1)*(theta-prm(2))
       this%d2bend=2.d0*prm(1)
    case('harm')
       en=0.5d0*prm(1)*(theta-prm(2))**2
       this%d1bend=prm(1)*(theta-prm(2))
       this%d2bend=prm(1)
    end select
  end subroutine set_angleopt

  subroutine set_vanderwaals(this,dr,prm,ptrm,en)
    implicit none
    class(interopt), intent(inout) :: this
    character(6), intent(in)       :: ptrm
    real(8), intent(in)            :: dr,prm(2)
    real(8), intent(out)           :: en
    select case(ptrm)
    case('charmm')
       en=prm(1)*((prm(2)/dr)**12-2.0d0*(prm(2)/dr)**6)
       this%d1bond=-12.d0*prm(1)*((prm(2)/dr)**12-(prm(2)/dr)**6)/dr
       this%d2bond=12.d0*prm(1)*(13.0d0*(prm(2)/dr)**12-7.0d0*(prm(2)/dr)**6)/dr**2
    end select
  end subroutine set_vanderwaals

  subroutine set_coulomb(this,dr,qi,qj,en)
    implicit none
    class(interopt), intent(inout) :: this
    real(8), intent(in)            :: dr,qi,qj
    real(8), intent(out)           :: en
    en=qi*qj/dr
    this%d1bond=-qi*qj/dr**2
    this%d2bond=2.0d0*qi*qj/dr**3
  end subroutine set_coulomb

  double precision function get_d1bond(this)
    implicit none
    class(interopt), intent(in) :: this
    get_d1bond=this%d1bond
  end function get_d1bond

  double precision function get_d2bond(this)
    implicit none
    class(interopt), intent(in) :: this
    get_d2bond=this%d2bond
  end function get_d2bond

  double precision function get_d1bend(this)
    implicit none
    class(interopt), intent(in) :: this
    get_d1bend=this%d1bend
  end function get_d1bend

  double precision function get_d2bend(this)
    implicit none
    class(interopt), intent(in) :: this
    get_d2bend=this%d2bend
  end function get_d2bend

end module interopt_module
