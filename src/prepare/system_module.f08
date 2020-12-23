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

  use constants_module

  private
  public :: system

  type, extends(constants) :: system
     integer                    :: nmol
     integer                    :: natom
     integer, allocatable       :: ntmol(:)
     integer, allocatable       :: nxmol(:)
     character(10), allocatable :: namemol(:)
     real(8)                    :: a
     real(8)                    :: b
     real(8)                    :: c
     real(8)                    :: v(3,3)
     real(8), allocatable       :: xa(:)
     real(8), allocatable       :: ya(:)
     real(8), allocatable       :: za(:)
   contains
     procedure :: sites
     procedure :: molecules
     procedure :: unit_cell
     procedure :: set_natom
     procedure :: get_natom
  end type system

contains

  subroutine molecules(this)
    implicit none
    class(system), intent(inout) :: this
    integer                      :: nx,i
    character(4)                 :: key
1   read(5,*,end=2)key
    if(key.ne.'&SYS')goto 1
    do while (key.ne.'&END')
       read(5,*)key
       if(key.eq.'nmol')then
          backspace(5)
          read(5,*)key,nx
          allocate(this%namemol(nx),this%ntmol(nx),this%nxmol(nx))
          do i=1,nx
             read(5,*)this%namemol(i),this%ntmol(i),this%nxmol(i)
          end do
       end if
    end do
    this%nmol=nx
    rewind(5)
2   return
  end subroutine molecules

  subroutine set_natom(this)
    implicit none
    class(system), intent(inout) :: this
    integer                      :: i,nx
    nx=0
    do i=1,this%nmol
       nx=nx+this%ntmol(i)*this%nxmol(i)
    end do
    this%natom=nx
  end subroutine set_natom

  function get_natom(this)
    implicit none
    class(system), intent(inout) :: this
    integer                      :: get_natom
    get_natom=this%natom
  end function get_natom

  subroutine sites(this)
    implicit none
    class(system), intent(inout) :: this
    integer                        :: i,nx
    character(2)                   :: at
    nx=this%get_natom()
    allocate(this%xa(nx),this%ya(nx),this%za(nx))
    read(9,*)nx
    read(9,*)
    if(nx.ne.this%get_natom())goto 1
    do i=1,this%get_natom()
       read(9,*)at,this%xa(i),this%ya(i),this%za(i)
    end do
    return
1   write(6,*)&
         'ERROR: The number of atoms in HICOLM.xyz does not match with defined in HICOLM.in!'
    write(6,*)'Hint: Change the HICOLM.xyz, or the inputs in HICOLM.in.'
    stop
  end subroutine sites

  subroutine unit_cell(this)
    implicit none
    class(system), intent(inout) :: this
    integer                      :: nx,i
    character(4)                 :: key
1   read(5,*,end=2)key
    if(key.ne.'&SYS')goto 1
    do while (key.ne.'&END')
       read(5,*)key
       if(key.eq.'cell')then
          backspace(5)
          read(5,*)key,this%v(1,1),this%v(2,2),this%v(3,3)
       end if
    end do
    this%nmol=nx
    this%a=sqrt(this%v(1,1)**2+this%v(1,2)**2+this%v(1,3)**2)
    this%b=sqrt(this%v(2,1)**2+this%v(2,2)**2+this%v(2,3)**2)
    this%c=sqrt(this%v(3,1)**2+this%v(3,2)**2+this%v(3,3)**2)
2   return
  end subroutine unit_cell

  subroutine system_check(this)
    implicit none
    class(system), intent(inout) :: this
  end subroutine system_check

end module system_module
