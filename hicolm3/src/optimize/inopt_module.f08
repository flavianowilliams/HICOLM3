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
module inopt_module
  !*******************************************************************************************
  !*******************************************************************************************

  use forcefield_module

  implicit none

  private
  public :: inopt

  type, extends(forcefield) :: inopt
     integer, private :: nstep
     real(8), private :: eps
   contains
     procedure :: set_nstep
     procedure :: get_nstep
     procedure :: set_eps
     procedure :: get_eps
  end type inopt

contains

  subroutine set_nstep(this,nstep)
    implicit none
    class(inopt), intent(inout) :: this
    integer, intent(in)         :: nstep
    this%nstep=nstep
  end subroutine set_nstep

  integer function get_nstep(this)
    implicit none
    class(inopt), intent(inout) :: this
    get_nstep=this%nstep
  end function get_nstep

  subroutine set_eps(this,eps)
    implicit none
    class(inopt), intent(inout) :: this
    real(8), intent(in)         :: eps
    this%eps=eps
  end subroutine set_eps

  double precision function get_eps(this)
    implicit none
    class(inopt), intent(inout) :: this
    get_eps=this%eps
  end function get_eps

end module inopt_module
