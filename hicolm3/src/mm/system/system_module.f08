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
  !********************************************************************************
  !********************************************************************************

  use zmatrix_module

  implicit none

  private
  public :: system

  type, extends(zmatrix) :: system
     real(8), private :: qtotal
     real(8), private :: mtotal
   contains
     procedure :: set_system
     procedure :: set_qtotal
     procedure :: get_qtotal
     procedure :: set_mtotal
     procedure :: get_mtotal
     procedure :: set_global
  end type system

contains

  subroutine set_system(this)
    implicit none
    class(system), intent(inout) :: this
    integer                      :: nmol,k,i,nxmolmax
    integer                      :: ntmol(10),nxmol(10),zatmol(10,1000)
    character(6)                 :: tpmol(10,1000)
    real(8)                      :: qatmol(10,1000)
    real(8)                      :: zmatrix_tol
    character(10)                :: namemol(10)
    character(11)                :: key
    logical                      :: check
    nmol=0
    nxmolmax=0
    check=.true.
    do while(check)
       read(5,*,end=1)key
       if(key.eq.'&SYSTEM'.or.key.eq.'&system')check=.false.
    end do
    check=.true.
    do while (check)
       read(5,*)key
       if(key.eq.'cell')then
          backspace(5)
          read(5,*)key,this%v(1,1),this%v(2,2),this%v(3,3)
       elseif(key.eq.'nmol')then
          backspace(5)
          read(5,*)key,nmol
          do i=1,nmol
             read(5,*)namemol(i),ntmol(i),nxmol(i)
             read(5,*)(zatmol(i,k),k=1,nxmol(i))
             read(5,*)(tpmol(i,k),k=1,nxmol(i))
             read(5,*)(qatmol(i,k),k=1,nxmol(i))
             nxmolmax=max(nxmolmax,nxmol(i))
          end do
       elseif(key.eq.'translate')then
          backspace(5)
          read(5,*)key,(this%sys_shift(i),i=1,3)
       elseif(key.eq.'zmatrix')then
          backspace(5)
          read(5,*)key,zmatrix_tol
          call this%set_zmatrix_tol(zmatrix_tol)
       elseif(key.eq.'&END_SYSTEM'.or.key.eq.'&end_system')then
          check=.false.
       end if
    end do
    call this%set_nmol(nmol)
    if(this%get_nmol().le.0)goto 2
    allocate(this%namemol(nmol),this%ntmol(nmol),this%nxmol(nmol))
    allocate(this%zatmol(nmol,nxmolmax))
    allocate(this%qatmol(nmol,nxmolmax))
    allocate(this%tpmol(nmol,nxmolmax))
    call this%set_namemol(namemol)
    call this%set_ntmol(ntmol)
    call this%set_nxmol(nxmol)
    call this%set_tpmol(tpmol)
    call this%set_zatmol(zatmol)
    call this%set_qatmol(qatmol)
1   rewind(5)
    return
2   write(6,*)'ERROR: The number of molecules does not be zero!'
    write(6,*)'Hint: Check the input in the &SYSTEM section.'
    stop
  end subroutine set_system

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

  subroutine set_global(this)
    class(system), intent(inout) :: this
    integer                        :: i,j
    real(8)                        :: mtotal
    call this%set_qtotal(this%qatmol)
    mtotal=0.d0
    do i=1,this%get_nmol()
       do j=1,this%nxmol(i)
          mtotal=mtotal+this%ntmol(i)*this%massmol(i,j)
       end do
    end do
    call this%set_mtotal(mtotal)
  end subroutine set_global

end module system_module
