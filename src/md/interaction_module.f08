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
  use coulomb_module

  implicit none

  integer i,j

  private
  public :: interaction

  type, extends(neighbourlist) :: interaction
     type(coulomb)    :: coul
     real(8), private :: enpot
     real(8), private :: virtot
   contains
     procedure :: interaction_prepare
     procedure :: set_forces
     procedure :: set_enpot
     procedure :: get_enpot
     procedure :: set_virtot
     procedure :: get_virtot
  end type interaction

contains

  subroutine interaction_prepare(this)
    implicit none
    class(interaction), intent(inout) :: this
    allocate(this%fax(this%get_natom()),this%fay(this%get_natom()),this%faz(this%get_natom()))
  end subroutine interaction_prepare

  subroutine set_forces(this)
    implicit none
    class(interaction), intent(inout) :: this
    integer                           :: ni,nj
    real(8)                           :: xvz,yvz,zvz,dr,enpot
    do i=1,this%get_natom()
       this%fax(i)=0.d0
       this%fay(i)=0.d0
       this%faz(i)=0.d0
    end do
    call this%coul%coulomb_prepare&
         (this%get_coulop(),this%get_kconv(),this%get_rcutoff(),this%get_pi())
    call this%coul%set_encoul(0.d0)
    call this%coul%set_vircoul(0.d0)
    enpot=0.d0
    do i=1,this%get_natom()
       do j=1,this%nlist(i)
          ni=i
          nj=this%ilist(i,j)
          call this%mic(ni,nj,xvz,yvz,zvz)
          dr=sqrt(xvz**2+yvz**2+zvz**2)
          if(abs(this%qat(ni)*this%qat(nj)).gt.1.d-8)&
               call this%coul%set_coulomb(dr,this%qat(ni),this%qat(nj))
          enpot=enpot+this%coul%get_encoul()
          print*,enpot,i
       end do
    end do
    call this%set_enpot(enpot)
  end subroutine set_forces

  subroutine set_enpot(this,enpot)
    implicit none
    class(interaction), intent(inout) :: this
    real(8), intent(in)               :: enpot
    this%enpot=enpot
  end subroutine set_enpot

  double precision function get_enpot(this)
    implicit none
    class(interaction), intent(inout) :: this
    get_enpot=this%enpot
  end function get_enpot

  subroutine set_virtot(this,virtot)
    implicit none
    class(interaction), intent(inout) :: this
    real(8), intent(in)           :: virtot
    this%virtot=virtot
  end subroutine set_virtot

  double precision function get_virtot(this)
    implicit none
    class(interaction), intent(inout) :: this
    get_virtot=this%virtot
  end function get_virtot

end module interaction_module
