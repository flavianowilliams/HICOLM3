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
module interopt_module
  !*******************************************************************************************
  !*******************************************************************************************

  use gradientdescent_module

  implicit none

  private
  public :: interopt

  type, extends(gradientdescent) :: interopt
     real(8), private :: maxforce
   contains
     procedure :: set_force
     procedure :: set_maxforce
     procedure :: get_maxforce
  end type interopt

contains

  subroutine set_maxforce(this)
    implicit none
    class(interopt), intent(inout) :: this
    integer                           :: i
    this%maxforce=0.d0
    do i=1,this%get_nmatrix()
       this%maxforce=max(this%maxforce,abs(this%res(i)))
    end do
  end subroutine set_maxforce

  double precision function get_maxforce(this)
    implicit none
    class(interopt), intent(in) :: this
    get_maxforce=this%maxforce
  end function get_maxforce

  subroutine set_force(this)
    implicit none
    class(interopt), intent(inout) :: this
    integer                               :: i,j,k,l
    real(8)                               :: xvz,yvz,zvz,dr
    real(8)                               :: prm(3)
    character(7)                          :: ptrm
    do i=1,this%get_nmatrix()
       this%res(i)=0.d0
       do j=1,this%get_nmatrix()
          this%hess(i,j)=0.d0
       end do
    end do
    do i=1,this%get_natom()
       do j=i+1,this%get_natom()
          call this%mic(i,j,xvz,yvz,zvz)
          dr=max(sqrt(xvz**2+yvz**2+zvz**2),1.d-8)
          do k=1,this%get_nvdw()
             if(this%tpa(i).eq.this%spcvdw(k,1).and.this%tpa(j).eq.this%spcvdw(k,2).or.&
                  this%tpa(i).eq.this%spcvdw(k,2).and.this%tpa(j).eq.this%spcvdw(k,1))then
                do l=1,2
                   prm(l)=this%parvdw(k,l)
                end do
                ptrm=this%tvdw(k)
                call this%set_vanderwaals(dr,prm,ptrm)
                call this%set_residue(i,j,dr,xvz,yvz,zvz)
                call this%set_hessian(i,j,dr,xvz,yvz,zvz)
             end if
          end do
       end do
    end do
    do i=1,this%get_nmatrix()
       do j=i+1,this%get_nmatrix()
          this%hess(j,i)=this%hess(i,j)
       end do
    end do
  end subroutine set_force

end module interopt_module
