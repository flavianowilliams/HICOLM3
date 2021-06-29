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
module gradientdescent_module
  !*******************************************************************************************
  !*******************************************************************************************

  use atoms_module

  implicit none

  private
  public :: gradientdescent

  type, extends(atoms) :: gradientdescent
     integer, private     :: nmatrix
     real(8), private     :: h1
     real(8), private     :: h2
     real(8), allocatable :: res(:)
     real(8), allocatable :: hess(:,:)
   contains
     procedure :: gd_init
     procedure :: set_vanderwaals
     procedure :: set_residue
     procedure :: set_hessian
     procedure :: set_nmatrix
     procedure :: get_nmatrix
  end type gradientdescent

contains

  subroutine gd_init(this)
    implicit none
    class(gradientdescent), intent(inout) :: this
    call this%set_nmatrix(3*this%get_natom())
    allocate(this%res(this%nmatrix),this%hess(this%nmatrix,this%nmatrix))
  end subroutine gd_init

  subroutine set_nmatrix(this,nmatrix)
    implicit none
    class(gradientdescent), intent(inout) :: this
    integer, intent(in)                   :: nmatrix
    this%nmatrix=nmatrix
  end subroutine set_nmatrix

  integer function get_nmatrix(this)
    implicit none
    class(gradientdescent), intent(in) :: this
    get_nmatrix=this%nmatrix
  end function get_nmatrix

  subroutine set_vanderwaals(this,dr,prm,ptrm)
    implicit none
    class(gradientdescent), intent(inout) :: this
    character(6), intent(in)              :: ptrm
    real(8), intent(in)                   :: dr,prm(2)
    select case(ptrm)
    case('charmm')
       this%h1=-12.d0*prm(1)*((prm(2)/dr)**12-(prm(2)/dr)**6)/dr
       this%h2=12.d0*prm(1)*(13.0d0*(prm(2)/dr)**12-7.0d0*(prm(2)/dr)**6)/dr**2
    end select
  end subroutine set_vanderwaals

  subroutine set_residue(this,i1,i2,dr,xvz,yvz,zvz)
    implicit none
    class(gradientdescent), intent(inout) :: this
    integer, intent(in)                   :: i1,i2
    integer                               :: ix,ixx
    real(8), intent(in)                   :: xvz,yvz,zvz,dr
    real(8)                               :: fr
    fr=-this%h1/dr
    ix=3*i1-2
    ixx=3*i2-2
    this%res(ix)=this%res(ix)-fr*xvz
    this%res(ix+1)=this%res(ix+1)-fr*yvz
    this%res(ix+2)=this%res(ix+2)-fr*zvz
    this%res(ixx)=this%res(ixx)+fr*xvz
    this%res(ixx+1)=this%res(ixx+1)+fr*yvz
    this%res(ixx+2)=this%res(ixx+2)+fr*zvz
  end subroutine set_residue

  subroutine set_hessian(this,i1,i2,dr,xvz,yvz,zvz)
    implicit none
    class(gradientdescent), intent(inout) :: this
    integer                               :: i1,i2,ix,ixx,i,j
    real(8), intent(in)                   :: xvz,yvz,zvz,dr
    real(8)                               :: h1,h2,dx(3)
    dx(1)=xvz
    dx(2)=yvz
    dx(3)=zvz
    h1=this%h1
    h2=this%h2
    ix=3*i1-2
    ixx=3*i2-2
    do i=1,3 !xvz,yvz,zvz
       do j=i,3 !xvz,yvz,zvz
          this%hess(ix+i-1,ix+j-1)=(h2/dr**2-h1/dr**3)*dx(i)*dx(j)+h1*kronij(i,j)/dr
          this%hess(ix+i-1,ixx+j-1)=-((h2/dr**2-h1/dr**3)*dx(i)*dx(j)+h1*kronij(i,j)/dr)
          this%hess(ixx+i-1,ixx+j-1)=(h2/dr**2-h1/dr**3)*dx(i)*dx(j)+h1*kronij(i,j)/dr
       end do
    end do
  end subroutine set_hessian

  integer function kronij(i,j)
    integer, intent(in) :: i,j
    kronij=int((float(i+j)-abs(i-j))/(float(i+j)+abs(i-j)))
  end function kronij

end module gradientdescent_module
