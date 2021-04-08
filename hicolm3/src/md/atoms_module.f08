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
module atoms_module
  !*******************************************************************************************
  !*******************************************************************************************

  use input_module

  implicit none

  private
  public :: atoms

  type, extends(input) :: atoms
     integer, allocatable      :: zat(:)
     real(8), allocatable      :: mass(:)
     real(8), allocatable      :: qat(:)
     real(8), allocatable      :: vax(:)
     real(8), allocatable      :: vay(:)
     real(8), allocatable      :: vaz(:)
     real(8), allocatable      :: fax(:)
     real(8), allocatable      :: fay(:)
     real(8), allocatable      :: faz(:)
     character(6), allocatable :: tpa(:)
   contains
     procedure :: set_zat
     procedure :: set_mass
     procedure :: set_qat
     procedure :: set_tpa
     procedure :: set_velocity
  end type atoms

contains

  subroutine set_zat(this)
    implicit none
    class(atoms), intent(inout) :: this
    integer                     :: nx,i,j,k
    allocate(this%zat(this%get_natom()))
    nx=1
    do i=1,this%get_nmol()
       do j=1,this%ntmol(i)
          do k=1,this%nxmol(i)
             this%zat(nx)=this%zatmol(i,k)
             nx=nx+1
          end do
       end do
    end do
  end subroutine set_zat

  subroutine set_mass(this)
    implicit none
    class(atoms), intent(inout) :: this
    integer                     :: nx,i,j,k
    allocate(this%mass(this%get_natom()))
    nx=1
    do i=1,this%get_nmol()
       do j=1,this%ntmol(i)
          do k=1,this%nxmol(i)
             this%mass(nx)=this%massmol(i,k)
             nx=nx+1
          end do
       end do
    end do
  end subroutine set_mass

  subroutine set_qat(this)
    implicit none
    class(atoms), intent(inout) :: this
    integer                     :: nx,i,j,k
    allocate(this%qat(this%get_natom()))
    nx=1
    do i=1,this%get_nmol()
       do j=1,this%ntmol(i)
          do k=1,this%nxmol(i)
             this%qat(nx)=this%qatmol(i,k)
             nx=nx+1
          end do
       end do
    end do
  end subroutine set_qat

  subroutine set_tpa(this)
    implicit none
    class(atoms), intent(inout) :: this
    integer                     :: nx,i,j,k
    allocate(this%tpa(this%get_natom()))
    nx=1
    do i=1,this%get_nmol()
       do j=1,this%ntmol(i)
          do k=1,this%nxmol(i)
             this%tpa(nx)=this%tpmol(i,k)
             nx=nx+1
          end do
       end do
    end do
  end subroutine set_tpa

  subroutine set_velocity(this)
    implicit none
    class(atoms), intent(inout) :: this
    integer                     :: i,j,k,nx
    allocate(this%vax(this%get_natom()),this%vay(this%get_natom()),this%vaz(this%get_natom()))
    nx=1
    do i=1,this%get_nmol()
       do j=1,this%ntmol(i)
          do k=1,this%nxmol(i)
             this%vax(nx)=sqrt(this%get_temp()/this%massmol(i,k))
             this%vay(nx)=sqrt(this%get_temp()/this%massmol(i,k))
             this%vaz(nx)=sqrt(this%get_temp()/this%massmol(i,k))
             nx=nx+1
          end do
       end do
    end do
  end subroutine set_velocity

end module atoms_module
