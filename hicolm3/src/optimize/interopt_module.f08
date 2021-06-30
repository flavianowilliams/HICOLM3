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
     real(8), private :: enpot
   contains
     procedure :: set_loop
     procedure :: set_maxforce
     procedure :: get_maxforce
     procedure :: set_enpot
     procedure :: get_enpot
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

  subroutine set_loop(this)
    implicit none
    class(interopt), intent(inout) :: this
    integer                               :: i,j,k,l,ni,nj,nx
    real(8)                               :: xvz,yvz,zvz,dr
    real(8)                               :: prm(3)
    character(7)                          :: ptrm
    character(6)                          :: ptrm2
    call this%set_enpot(0.d0)
    do i=1,this%get_nmatrix()
       this%res(i)=0.d0
       do j=1,this%get_nmatrix()
          this%hess(i,j)=0.d0
       end do
    end do
    nx=0
    do i=1,this%get_nmol()
       do j=1,this%ntmol(i)
          do k=1,this%bondscnt(i)
             ni=nx+this%molbond(i,k,1)
             nj=nx+this%molbond(i,k,2)
             call this%mic(ni,nj,xvz,yvz,zvz)
             dr=sqrt(xvz**2+yvz**2+zvz**2)
             do l=1,2
                prm(l)=this%parbnd(i,k,l)
             end do
             ptrm2=this%tbonds(i,k)
             call this%set_bondopt(dr,prm,ptrm2)
             call this%set_residue(ni,nj,dr,xvz,yvz,zvz)
             call this%set_hessian(ni,nj,dr,xvz,yvz,zvz)
             call this%set_enpot(this%get_en())
          end do
          nx=nx+this%nxmol(i)
       end do
    end do
    do i=1,this%get_natom()
       do j=1,this%nlist(i)
          ni=i
          nj=this%ilist(i,j)
          call this%mic(ni,nj,xvz,yvz,zvz)
          dr=max(sqrt(xvz**2+yvz**2+zvz**2),1.d-8)
          if(abs(this%qat(ni)*this%qat(nj)).gt.1.d-8)then
             call this%set_coulomb(dr,this%qat(ni),this%qat(nj))
             call this%set_residue(ni,nj,dr,xvz,yvz,zvz)
             call this%set_hessian(ni,nj,dr,xvz,yvz,zvz)
             call this%set_enpot(this%get_en())
          end if
          do k=1,this%get_nvdw()
             if(this%tpa(ni).eq.this%spcvdw(k,1).and.this%tpa(nj).eq.this%spcvdw(k,2).or.&
                  this%tpa(ni).eq.this%spcvdw(k,2).and.this%tpa(nj).eq.this%spcvdw(k,1))then
                do l=1,2
                   prm(l)=this%parvdw(k,l)
                end do
                ptrm=this%tvdw(k)
                call this%set_vanderwaals(dr,prm,ptrm)
                call this%set_residue(ni,nj,dr,xvz,yvz,zvz)
                call this%set_hessian(ni,nj,dr,xvz,yvz,zvz)
                call this%set_enpot(this%get_en())
             end if
          end do
       end do
    end do
    do i=1,this%get_nmatrix()
       do j=i+1,this%get_nmatrix()
          this%hess(j,i)=this%hess(i,j)
       end do
    end do
  end subroutine set_loop

  subroutine set_enpot(this,enpot)
    implicit none
    class(interopt), intent(inout) :: this
    real(8), intent(in)            :: enpot
    this%enpot=this%enpot+enpot
  end subroutine set_enpot

  double precision function get_enpot(this)
    implicit none
    class(interopt), intent(in) :: this
    get_enpot=this%enpot
  end function get_enpot

end module interopt_module
