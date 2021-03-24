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
module system_module
  !*******************************************************************************************
  !*******************************************************************************************

  implicit none

  private
  public :: system

  type :: system
     real(8), private :: qtotal
     real(8), private :: mtotal
   contains
     procedure :: set_qtotal
     procedure :: get_qtotal
     procedure :: set_mtotal
     procedure :: get_mtotal
  end type system

contains

  subroutine set_qtotal(this,qatom)
    class(system), intent(inout) :: this
    real(8), intent(in)          :: qatom(:,:)
    this%qtotal=sum(qatom)
  end subroutine set_qtotal

  double precision function get_qtotal(this)
    class(system), intent(in) :: this
    get_qtotal=this%qtotal
  end function get_qtotal

  subroutine set_mtotal(this,mtotal)
    class(system), intent(inout) :: this
    real(8), intent(in)          :: mtotal
    this%mtotal=mtotal
  end subroutine set_mtotal

  double precision function get_mtotal(this)
    class(system), intent(in) :: this
    get_mtotal=this%mtotal
  end function get_mtotal

end module system_module
