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
module ensemble_module
  !*******************************************************************************************
  !*******************************************************************************************

  use thermodynamics_module
  use nve_module

  implicit none

  private
  public :: ensemble

  type, extends(thermodynamics) :: ensemble
     type(nve) :: nve
   contains
     procedure :: set_ensemble
  end type ensemble

contains

  subroutine set_ensemble(this)
    implicit none
    class(ensemble), intent(inout) :: this
    call this%nve%teste()
    call this%set_ekinetic()
    call this%set_etotal()
    call this%set_temperature()
    call this%set_pressure()
  end subroutine set_ensemble

end module ensemble_module