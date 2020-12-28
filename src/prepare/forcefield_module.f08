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
module forcefield_module
  !*******************************************************************************************
  !*******************************************************************************************

  use zmatrix_module
  use amber_module

  implicit none

  integer i,j,k,l

  private
  public :: forcefield

  type, extends(zmatrix) :: forcefield
     type(amber)          :: amber
     real(8), allocatable :: parbnd(:,:,:)
   contains
     procedure :: set_potentials
     procedure :: set_extra_potentials
  end type forcefield

contains

  subroutine set_potentials(this)
    implicit none
    class(forcefield), intent(inout) :: this
    integer                          :: i1,i2
    allocate(this%parbnd(this%nmol,2,2))
    call this%amber%set_amber()
    do i=1,this%nmol
       do j=1,this%bondscnt(i)
          i1=this%molbond(i,j,1)
          i2=this%molbond(i,j,2)
          do k=1,115
             if(this%tpmol(i,i1).eq.this%amber%tpam(k,1).and.&
                  this%tpmol(i,i2).eq.this%amber%tpam(k,2))then
                do l=1,2
                   this%parbnd(i,j,l)=this%amber%prms(k,l)
                end do
             end if
          end do
       end do
    end do
  end subroutine set_potentials

  subroutine set_extra_potentials(this)
    implicit none
    class(forcefield), intent(inout) :: this
  end subroutine set_extra_potentials

end module forcefield_module
