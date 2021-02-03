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
module moleculardynamics_module
  !*******************************************************************************************
  !*******************************************************************************************

  use input_module

  implicit none

  private
  public :: moleculardynamics

  type, extends(input) :: moleculardynamics
   contains
  end type moleculardynamics

  interface moleculardynamics
     module procedure constructor
  end interface moleculardynamics

contains

  type(moleculardynamics) function constructor()
    implicit none
    call constructor%set_nstep(1000)
    call constructor%set_nrelax(1000)
    call constructor%set_nframes(200)
    call constructor%set_timestep(0.001d0)
    call constructor%set_press(1.d0)
    call constructor%set_temp(298.d0)
    call constructor%set_rcutoff(8.0d0)
    call constructor%set_drcutoff(0.1d0)
  end function constructor

end module moleculardynamics_module
