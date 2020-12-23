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
module prepare_module
  !*******************************************************************************************
  ! class to prepare the physical environment for simulation                                 *
  ! HICOLM.sys and HICOLM.top will be create                                                 *
  !*******************************************************************************************

  use system_module

  private
  public :: prepare

  type, extends(system) :: prepare
   contains
     procedure, private :: prepare_init
     procedure          :: check
     procedure          :: print_sys
  end type prepare

  interface prepare
     module procedure constructor
  end interface prepare

contains

  type(prepare) function constructor()
    implicit none
    call constructor%prepare_init()
  end function constructor

  subroutine prepare_init(this)
    implicit none
    class(prepare), intent(inout) :: this
    integer                       :: i,j
    do i=1,3
       do j=1,3
          this%v(i,j)=0.d0
       end do
    end do
  end subroutine prepare_init

  subroutine check(this)
    implicit none
    class(prepare), intent(inout) :: this
    integer                       :: i
    if(this%nmol.le.0)then
       write(6,*)'ERROR: The number of types of molecules does not be zero!'
       write(6,*)'Hint: Check the input in the HICOLM.in.'
       stop
    else
       do i=1,this%nmol
          if(this%ntmol(i).le.0)then
             write(6,*)'ERROR: The number of molecules does not be zero!'
             write(6,*)'Hint: Check the input in the HICOLM.in.'
             stop
          end if
          if(this%nxmol(i).le.0)then
             write(6,*)'ERROR: The number of sites at each molecule does not be zero!'
             write(6,*)'Hint: Check the input in the HICOLM.in.'
             stop
          end if
       end do
    end if
    if(this%a.le.1.e-4.or.this%b.le.1.e-4.or.this%c.le.1.e-4)then
       write(6,*)'ERROR: Lattice constant too small!'
       write(6,*)'Hint: Increase the value in the HICOLM.in.'
       stop
    end if
  end subroutine check

  subroutine print_sys(this)
    implicit none
    class(prepare), intent(inout) :: this
    integer                       :: i,j
    open(10,file='HICOLM.sys',status='unknown')
    do i=1,3
       write(10,'(5x,3f16.8)')(this%v(i,j),j=1,3)
    end do
    do i=1,this%get_natom()
       write(10,'(a5,3f16.8)')'#',this%xa(i),this%ya(i),this%za(i)
    end do
  end subroutine print_sys

end module prepare_module
