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
module vanderwaals_module
  !*******************************************************************************************
  !*******************************************************************************************

  implicit none

!  integer i,j,k,l

  private
  public :: vanderwaals

  type :: vanderwaals
     real(8), private      :: envdw
     real(8), private      :: virvdw
     real(8), private      :: rcutoff
     real(8), private      :: force
     real(8), private      :: vdwcorr
   contains
     procedure :: set_vanderwaals
     procedure :: vanderwaals_prepare
     procedure :: set_envdw
     procedure :: get_envdw
     procedure :: set_virvdw
     procedure :: get_virvdw
     procedure :: get_force
     procedure :: set_vdwcorr
     procedure :: get_vdwcorr
  end type vanderwaals

contains

  subroutine vanderwaals_prepare(this,rcutoff)
    implicit none
    class(vanderwaals), intent(inout) :: this
    real(8), intent(in)               :: rcutoff
    this%rcutoff=rcutoff
  end subroutine vanderwaals_prepare

  subroutine set_vanderwaals(this,dr,prm,ptrm)
    implicit none
    class(vanderwaals), intent(inout) :: this
    character(4), intent(in)          :: ptrm
    real(8), intent(in)               :: dr,prm(2)
    select case(ptrm)
    case('amber')
       this%envdw=prm(1)*((prm(2)/dr)**12-2.d0*(prm(2)/dr)**6)
       this%force=12.d0*prm(1)*((prm(2)/dr)**12-(prm(2)/dr)**6)/dr**2
    case('lj')
       this%envdw=4.d0*prm(1)*((prm(2)/dr)**12-(prm(2)/dr)**6)
       this%force=24.d0*prm(1)*(2.d0*(prm(2)/dr)**12-(prm(2)/dr)**6)/dr**2
    end select
  end subroutine set_vanderwaals

  subroutine set_envdw(this,envdw)
    implicit none
    class(vanderwaals), intent(inout) :: this
    real(8), intent(in)           :: envdw
    this%envdw=envdw
  end subroutine set_envdw

  double precision function get_envdw(this)
    implicit none
    class(vanderwaals), intent(inout) :: this
    get_envdw=this%envdw
  end function get_envdw

  subroutine set_virvdw(this,virvdw)
    implicit none
    class(vanderwaals), intent(inout) :: this
    real(8), intent(in)               :: virvdw
    this%virvdw=virvdw
  end subroutine set_virvdw

  double precision function get_virvdw(this)
    implicit none
    class(vanderwaals), intent(inout) :: this
    get_virvdw=this%virvdw
  end function get_virvdw

  subroutine set_vdwcorr(this)
    implicit none
    class(vanderwaals), intent(inout) :: this
  end subroutine set_vdwcorr

  double precision function get_vdwcorr(this)
    implicit none
    class(vanderwaals), intent(inout) :: this
    get_vdwcorr=this%vdwcorr
  end function get_vdwcorr

  double precision function get_force(this)
    implicit none
    class(vanderwaals), intent(inout) :: this
    get_force=this%force
  end function get_force

end module vanderwaals_module
