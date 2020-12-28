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
module amber_module
  !*******************************************************************************************
  !*******************************************************************************************

  implicit none

  integer i

  private
  public :: amber

  type :: amber
     integer      :: namber
     real(8)      :: prms(115,2)
     character(2) :: tpam(115,2)
   contains
     procedure :: amber_init
     procedure :: set_amber
  end type amber

  interface amber
     module procedure constructor
  end interface amber

contains

  type(amber) function constructor(namber)
    implicit none
    integer, intent(in) :: namber
    call constructor%amber_init(namber)
  end function constructor

  subroutine amber_init(this,namber)
    implicit none
    class(amber), intent(inout) :: this
    integer, intent(in)         :: namber
    this%namber=namber
  end subroutine amber_init

  subroutine set_amber(this)
    implicit none
    class(amber), intent(inout) :: this
    real(4) x1,x2
    character(2) pa,pb
    open(12,file='/tmp/amber/amber_bonds.prm',status='old')
    do i=1,115
       read(12,*,end=1)pa,pb,x1,x2
       this%tpam(i,1)=pa
       this%tpam(i,2)=pb
       this%prms(i,1)=dble(x1)
       this%prms(i,2)=dble(x2)
    end do
1  close(12)
  end subroutine set_amber

end module amber_module
