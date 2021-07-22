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

  use interopt_module

  implicit none

  private
  public :: gradientdescent

  type, extends(interopt) :: gradientdescent
     integer, private     :: nmatrix
     real(8), private     :: lsearch
     real(8), private     :: enpot
     real(8), private     :: maxforce
     real(8), allocatable :: res(:)
     real(8), allocatable :: hess(:,:)
   contains
     procedure :: gd_init
     procedure :: set_loop
     procedure :: set_residue2
     procedure :: set_residue3
     generic   :: set_residue => set_residue2, set_residue3
     procedure :: set_hessian
     procedure :: set_nmatrix
     procedure :: get_nmatrix
     procedure :: set_lsearch
     procedure :: get_lsearch
     procedure :: set_enpot
     procedure :: get_enpot
     procedure :: set_maxforce
     procedure :: get_maxforce
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

  subroutine set_loop(this)
    implicit none
    class(gradientdescent), intent(inout) :: this
    integer                               :: i,j,k,l,ni,nj,nk,nx
    real(8)                               :: xvz,yvz,zvz,dr,en,enpot,dr1,dr2,theta
    real(8)                               :: prm(3),drij(3),drik(3)
    character(7)                          :: ptrm
    character(6)                          :: ptrm2
    enpot=0.d0
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
             call this%set_bondopt(dr,prm,ptrm2,en)
             call this%set_residue(ni,nj,dr,xvz,yvz,zvz)
             call this%set_hessian(ni,nj,dr,xvz,yvz,zvz)
             enpot=enpot+en
          end do
          do k=1,this%bendscnt(i)
             ni=nx+this%molbend(i,k,2)
             nj=nx+this%molbend(i,k,1)
             nk=nx+this%molbend(i,k,3)
             call this%mic(ni,nj,drij(1),drij(2),drij(3))
             call this%mic(ni,nk,drik(1),drik(2),drik(3))
             dr1=sqrt(drij(1)**2+drij(2)**2+drij(3)**2)
             dr2=sqrt(drik(1)**2+drik(2)**2+drik(3)**2)
             theta=acos((drij(1)*drik(1)+drij(2)*drik(2)+drij(3)*drik(3))/(dr1*dr2))
             do l=1,2
                prm(l)=this%parbend(i,k,l)
             end do
             ptrm=this%tbends(i,k)
             call this%set_angleopt(theta,prm,ptrm,en)
             call this%set_residue(ni,nj,nk,drij,drik,dr1,dr2,theta)
             enpot=enpot+en
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
             call this%set_coulomb(dr,this%qat(ni),this%qat(nj),en)
             call this%set_residue(ni,nj,dr,xvz,yvz,zvz)
             call this%set_hessian(ni,nj,dr,xvz,yvz,zvz)
             enpot=enpot+en
          end if
          do k=1,this%get_nvdw()
             if(this%tpa(ni).eq.this%spcvdw(k,1).and.this%tpa(nj).eq.this%spcvdw(k,2).or.&
                  this%tpa(ni).eq.this%spcvdw(k,2).and.this%tpa(nj).eq.this%spcvdw(k,1))then
                do l=1,2
                   prm(l)=this%parvdw(k,l)
                end do
                ptrm=this%tvdw(k)
                call this%set_vanderwaals(dr,prm,ptrm,en)
                call this%set_residue(ni,nj,dr,xvz,yvz,zvz)
                call this%set_hessian(ni,nj,dr,xvz,yvz,zvz)
                enpot=enpot+en
             end if
          end do
       end do
    end do
    call this%set_enpot(enpot)
    do i=1,this%get_nmatrix()
       do j=i+1,this%get_nmatrix()
          this%hess(j,i)=this%hess(i,j)
       end do
    end do
  end subroutine set_loop

  subroutine set_enpot(this,enpot)
    implicit none
    class(gradientdescent), intent(inout) :: this
    real(8), intent(in)            :: enpot
    this%enpot=this%enpot+enpot
  end subroutine set_enpot

  double precision function get_enpot(this)
    implicit none
    class(gradientdescent), intent(in) :: this
    get_enpot=this%enpot
  end function get_enpot

  subroutine set_residue2(this,i1,i2,dr,xvz,yvz,zvz)
    implicit none
    class(gradientdescent), intent(inout) :: this
    integer, intent(in)                   :: i1,i2
    integer                               :: ix,ixx
    real(8), intent(in)                   :: xvz,yvz,zvz,dr
    real(8)                               :: fr
    fr=-this%get_d1bond()/dr
    ix=3*i1-2
    ixx=3*i2-2
    this%res(ix)=this%res(ix)-fr*xvz
    this%res(ix+1)=this%res(ix+1)-fr*yvz
    this%res(ix+2)=this%res(ix+2)-fr*zvz
    this%res(ixx)=this%res(ixx)+fr*xvz
    this%res(ixx+1)=this%res(ixx+1)+fr*yvz
    this%res(ixx+2)=this%res(ixx+2)+fr*zvz
  end subroutine set_residue2

  subroutine set_residue3(this,i1,i2,i3,drij,drik,dr1,dr2,theta)
    implicit none
    class(gradientdescent), intent(inout) :: this
    integer, intent(in)               :: i1,i2,i3
    integer                           :: ix(3),ix1,ix2,ix3,i,j
    real(8), intent(in)               :: dr1,dr2,theta
    real(8), intent(in)               :: drij(3),drik(3)
    real(8)                           :: derij(3,3),fa
    ix(1)=i1
    ix(2)=i2
    ix(3)=i3
    do j=1,3
       do i=1,3
          derij(i,j)=(kronij(ix(i),ix(2))-kronij(ix(i),ix(1)))*drik(j)/(dr1*dr2) &
               +(kronij(ix(i),ix(3))-kronij(ix(i),ix(1)))*drij(j)/(dr1*dr2) &
               -cos(theta)*((kronij(ix(i),ix(2))-kronij(ix(i),ix(1)))*drij(j)/dr1**2 &
               +(kronij(ix(i),ix(3))-kronij(ix(i),ix(1)))*drik(j)/dr2**2)
       end do
    end do
    ix1=3*ix(1)-2
    ix2=3*ix(2)-2
    ix3=3*ix(3)-2
    fa=this%get_d1bend()/sin(theta)
    this%res(ix1)=this%res(ix1)+fa*derij(1,1)
    this%res(ix1+1)=this%res(ix1+1)+fa*derij(1,2)
    this%res(ix1+2)=this%res(ix1+2)+fa*derij(1,3)
    this%res(ix2)=this%res(ix2)+fa*derij(2,1)
    this%res(ix2+1)=this%res(ix2+1)+fa*derij(2,2)
    this%res(ix2+2)=this%res(ix2+2)+fa*derij(2,3)
    this%res(ix3)=this%res(ix3)+fa*derij(3,1)
    this%res(ix3+1)=this%res(ix3+1)+fa*derij(3,2)
    this%res(ix3+2)=this%res(ix3+2)+fa*derij(3,3)
  end subroutine set_residue3

  subroutine set_hessian(this,i1,i2,dr,xvz,yvz,zvz)
    implicit none
    class(gradientdescent), intent(inout) :: this
    integer                               :: i1,i2,ix,ixx,i,j
    real(8), intent(in)                   :: xvz,yvz,zvz,dr
    real(8)                               :: h1,h2,dx(3)
    dx(1)=xvz
    dx(2)=yvz
    dx(3)=zvz
    h1=this%get_d1bond()
    h2=this%get_d2bond()
    ix=3*i1-2
    ixx=3*i2-2
    do i=1,3 !xvz,yvz,zvz
       do j=i,3 !xvz,yvz,zvz
          this%hess(ix+i-1,ix+j-1)=this%hess(ix+i-1,ix+j-1)&
               +(h2/dr**2-h1/dr**3)*dx(i)*dx(j)+h1*kronij(i,j)/dr
          this%hess(ix+i-1,ixx+j-1)=this%hess(ix+i-1,ixx+j-1)&
               -((h2/dr**2-h1/dr**3)*dx(i)*dx(j)+h1*kronij(i,j)/dr)
          this%hess(ixx+i-1,ixx+j-1)=this%hess(ixx+i-1,ixx+j-1)&
               +(h2/dr**2-h1/dr**3)*dx(i)*dx(j)+h1*kronij(i,j)/dr
       end do
    end do
  end subroutine set_hessian

  integer function kronij(i,j)
    integer, intent(in) :: i,j
    kronij=int((float(i+j)-abs(i-j))/(float(i+j)+abs(i-j)))
  end function kronij

  subroutine set_maxforce(this)
    implicit none
    class(gradientdescent), intent(inout) :: this
    integer                           :: i
    this%maxforce=0.d0
    do i=1,this%get_nmatrix()
       this%maxforce=max(this%maxforce,abs(this%res(i)))
    end do
  end subroutine set_maxforce

  double precision function get_maxforce(this)
    implicit none
    class(gradientdescent), intent(in) :: this
    get_maxforce=this%maxforce
  end function get_maxforce

  subroutine set_lsearch(this,lsearch)
    implicit none
    class(gradientdescent), intent(inout) :: this
    real(8), intent(in)                   :: lsearch
    this%lsearch=lsearch
  end subroutine set_lsearch

  double precision function get_lsearch(this)
    implicit none
    class(gradientdescent), intent(in) :: this
    get_lsearch=this%lsearch
  end function get_lsearch

end module gradientdescent_module
