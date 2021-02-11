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
module interaction_module
  !*******************************************************************************************
  !*******************************************************************************************

  use neighbourlist_module

  implicit none

  integer i,j

  private
  public :: interaction

  type, extends(neighbourlist) :: interaction
     real(8), private :: encoul
     real(8), private :: vircoul
   contains
     procedure :: interaction_prepare
     procedure :: set_forces
     procedure :: set_encoul
     procedure :: get_encoul
     procedure :: set_vircoul
     procedure :: get_vircoul
  end type interaction

contains

  subroutine set_forces(this)
    implicit none
    class(interaction), intent(inout) :: this
    integer                           :: ni,nj
    real(8)                           :: xvz,yvz,zvz
    do i=1,this%get_natom()
       this%fax(i)=0.d0
       this%fay(i)=0.d0
       this%faz(i)=0.d0
    end do
    call this%set_encoul(0.d0)
    call this%set_vircoul(0.d0)
    do i=1,this%get_natom()
       do j=1,this%nlist(i)
          ni=i
          nj=this%ilist(i,j)
          call this%mic(ni,nj,xvz,yvz,zvz)
          if(abs(this%qat(ni)*this%qat(nj)).gt.1.d-8)print*,'ok'
       end do
    end do
  end subroutine set_forces

  subroutine interaction_prepare(this)
    implicit none
    class(interaction), intent(inout) :: this
    allocate(this%fax(this%get_natom()),this%fay(this%get_natom()),this%faz(this%get_natom()))
  end subroutine interaction_prepare

  subroutine set_encoul(this,encoul)
    implicit none
    class(interaction), intent(inout) :: this
    real(8), intent(in)               :: encoul
    this%encoul=encoul
  end subroutine set_encoul

  double precision function get_encoul(this)
    implicit none
    class(interaction), intent(inout) :: this
    get_encoul=this%encoul
  end function get_encoul

  subroutine set_vircoul(this,vircoul)
    implicit none
    class(interaction), intent(inout) :: this
    real(8), intent(in)               :: vircoul
    this%vircoul=vircoul
  end subroutine set_vircoul

  double precision function get_vircoul(this)
    implicit none
    class(interaction), intent(inout) :: this
    get_vircoul=this%vircoul
  end function get_vircoul

end module interaction_module
