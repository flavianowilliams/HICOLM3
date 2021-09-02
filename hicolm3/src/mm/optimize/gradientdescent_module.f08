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
  !********************************************************************************
  !********************************************************************************

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
     real(8)              :: dralpha
     real(8)              :: drbetaalpha
   contains
     procedure :: gd_init
     procedure :: set_loop
     procedure :: set_dralpha
     procedure :: get_dralpha
     procedure :: set_drbetaalpha
     procedure :: get_drbetaalpha
     procedure :: set_residue2
     procedure :: set_residue3
     generic   :: set_residue => set_residue2, set_residue3
     procedure :: set_hessian2
     procedure :: set_hessian3
     generic   :: set_hessian => set_hessian2, set_hessian3
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
    real(8)                               :: dr,en,enpot,dr1,dr2,theta
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
             call this%mic(ni,nj,drij(1),drij(2),drij(3))
             dr=sqrt(drij(1)**2+drij(2)**2+drij(3)**2)
             do l=1,2
                prm(l)=this%parbnd(i,k,l)
             end do
             ptrm2=this%tbonds(i,k)
             call this%set_bondopt(dr,prm,ptrm2,en)
             call this%set_residue(ni,nj,dr,drij)
             call this%set_hessian(ni,nj,dr,drij)
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
             theta=acos&
                  ((drij(1)*drik(1)+drij(2)*drik(2)+drij(3)*drik(3))/(dr1*dr2))
             do l=1,2
                prm(l)=this%parbend(i,k,l)
             end do
             ptrm=this%tbends(i,k)
             call this%set_angleopt(theta,prm,ptrm,en)
             call this%set_residue(ni,nj,nk,drij,drik,dr1,dr2,theta)
             call this%set_hessian(ni,nj,nk,drij,drik,dr1,dr2,theta)
             enpot=enpot+en
          end do
          nx=nx+this%nxmol(i)
       end do
    end do
    do i=1,this%get_natom()
       do j=1,this%nlist(i)
          ni=i
          nj=this%ilist(i,j)
          call this%mic(ni,nj,drij(1),drij(2),drij(3))
          dr=max(sqrt(drij(1)**2+drij(2)**2+drij(3)**2),1.d-8)
          if(abs(this%qat(ni)*this%qat(nj)).gt.1.d-8)then
             call this%set_coulomb(dr,this%qat(ni),this%qat(nj),en)
             call this%set_residue(ni,nj,dr,drij)
             call this%set_hessian(ni,nj,dr,drij)
             enpot=enpot+en
          end if
          do k=1,this%get_nvdw()
             if(this%tpa(ni).eq.this%spcvdw(k,1).and.this%tpa(nj)&
                  .eq.this%spcvdw(k,2).or.&
                  this%tpa(ni).eq.this%spcvdw(k,2).and.this%tpa(nj)&
                  .eq.this%spcvdw(k,1))then
                do l=1,2
                   prm(l)=this%parvdw(k,l)
                end do
                ptrm=this%tvdw(k)
                call this%set_vanderwaals(dr,prm,ptrm,en)
                call this%set_residue(ni,nj,dr,drij)
                call this%set_hessian(ni,nj,dr,drij)
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

  subroutine set_residue2(this,i1,i2,dr,drij)
    implicit none
    class(gradientdescent), intent(inout) :: this
    integer, intent(in)                   :: i1,i2
    integer                               :: ix(2),i,j
    real(8), intent(in)                   :: drij(3),dr
    real(8)                               :: fr
    fr=-this%get_d1bond()/dr
    ix(1)=3*i1-2
    ix(2)=3*i2-2
    do i=1,2 !i,j
       do j=1,3 !x,y,z
          this%res(ix(i)+j-1)=this%res(ix(i)+j-1)&
               +fr*drij(j)*(kronij(2,i)-kronij(1,i))
       end do
    end do
  end subroutine set_residue2

  subroutine set_residue3(this,i1,i2,i3,drij,drik,dr1,dr2,theta)
    implicit none
    class(gradientdescent), intent(inout) :: this
    integer, intent(in)               :: i1,i2,i3
    integer                           :: i,j,ix(3)
    real(8), intent(in)               :: theta,drij(3),drik(3),dr1,dr2
    real(8)                           :: fa
    ix(1)=i1 !i
    ix(2)=i2 !j
    ix(3)=i3 !k
    fa=this%get_d1bend()/sin(theta)
    do i=1,3 !i,j,k
       do j=1,3 !x,y,z
          call this%set_dralpha(i1,i2,i3,ix(i),drij(j),drik(j),dr1,dr2,theta)
          this%res(3*ix(i)+j-3)=this%res(3*ix(i)+j-3)+fa*this%dralpha
       end do
    end do
  end subroutine set_residue3

  subroutine set_hessian2(this,i1,i2,dr,drij)
    implicit none
    class(gradientdescent), intent(inout) :: this
    integer                               :: i1,i2,ix(2),i,j,k,l
    real(8), intent(in)                   :: drij(3),dr
    real(8)                               :: h1,h2,dx(3)
    dx(1)=drij(1)
    dx(2)=drij(2)
    dx(3)=drij(3)
    h1=this%get_d1bond()
    h2=this%get_d2bond()
    ix(1)=i1 !alpha
    ix(2)=i2 !beta
    do i=1,3 !x,y,z
       do j=i,3 !x,y,z
          do k=1,2 !alpha, beta
             do l=k,2 !alpha, beta
                this%hess(3*ix(k)+i-3,3*ix(l)+j-3)=&
                     this%hess(3*ix(k)+i-3,3*ix(l)+j-3)&
                     +((h2/dr**2-h1/dr**3)*dx(i)*dx(j)+h1*kronij(i,j)/dr)&
                     *(kronij(ix(k),ix(2))-kronij(ix(k),ix(1)))&
                     *(kronij(ix(l),ix(2))-kronij(ix(l),ix(1)))
             end do
          end do
       end do
    end do
  end subroutine set_hessian2

  subroutine set_hessian3(this,i1,i2,i3,drij,drik,dr1,dr2,theta)
    implicit none
    class(gradientdescent), intent(inout) :: this
    integer, intent(in)                   :: i1,i2,i3
    integer                               :: ix(3),i,j,k,l
    real(8), intent(in)                   :: dr1,dr2,theta
    real(8), intent(in)                   :: drij(3),drik(3)
    real(8)                               :: h1,h2,dalpha,dbeta,dbetaalpha
    h1=this%get_d1bend()
    h2=this%get_d2bend()
    ix(1)=i1
    ix(2)=i2
    ix(3)=i3
    do i=1,3 !x,y,z
       do j=i,3 !x,y,z
          do k=1,3 !i1, i2, i3
             do l=k,3 !i1, i2, i3
                call this%set_dralpha&
                     (ix(1),ix(2),ix(3),ix(k),drij(i),drik(i),dr1,dr2,theta)
                dalpha=this%dralpha
                call this%set_dralpha&
                     (ix(1),ix(2),ix(3),ix(l),drij(j),drik(j),dr1,dr2,theta)
                dbeta=this%dralpha
                call this%set_drbetaalpha(ix(1),ix(2),ix(3),ix(k),ix(l),i,j,&
                     drij(i),drik(i),drij(j),drik(j),dr1,dr2,theta)
                dbetaalpha=this%drbetaalpha
                this%hess(3*ix(k)+i-3,3*ix(l)+j-3)=&
                     this%hess(3*ix(k)+i-3,3*ix(l)+j-3)&
                     +(h2-h1*cos(theta)/sin(theta))*dbeta*dalpha/sin(theta)**2&
                     -h1*dbetaalpha/sin(theta)
             end do
          end do
       end do
    end do
  end subroutine set_hessian3

  subroutine set_dralpha(this,ai,aj,ak,aa,drij,drik,dr1,dr2,theta)
    implicit none
    class(gradientdescent), intent(inout) :: this
    integer, intent(in)                   :: ai,aj,ak,aa
    real(8), intent(in)                   :: dr1,dr2,theta,drij,drik
    this%dralpha=drik*(kronij(aa,aj)-kronij(aa,ai))/(dr1*dr2)&
         +drij*(kronij(aa,ak)-kronij(aa,ai))/(dr1*dr2)&
         -cos(theta)*(drij*(kronij(aa,aj)-kronij(aa,ai))/dr1**2 &
         +drik*(kronij(aa,ak)-kronij(aa,ai))/dr2**2)
  end subroutine set_dralpha

  double precision function get_dralpha(this)
    implicit none
    class(gradientdescent), intent(in) :: this
    get_dralpha=this%dralpha
  end function get_dralpha

  subroutine set_drbetaalpha &
       (this,ai,aj,ak,aa,ab,xx,xl,drija,drika,drijb,drikb,dr1,dr2,theta)
    implicit none
    class(gradientdescent), intent(inout) :: this
    integer, intent(in)                   :: ai,aj,ak,aa,ab,xx,xl
    real(8), intent(in)                   :: dr1,dr2,theta,drija,drika,drijb,drikb
    real(8)                               :: func
    func=((kronij(xx,xl)/(dr1*dr2)-drika*drikb/(dr1*dr2**3))&
         *(kronij(ab,ak)-kronij(ab,ai))-drika*drijb/(dr2*dr1**3)&
         *(kronij(ab,aj)-kronij(ab,ai)))*(kronij(aa,aj)-kronij(aa,ai))
    func=func+((kronij(xx,xl)/(dr1*dr2)-drija*drijb/(dr2*dr1**3))&
         *(kronij(ab,aj)-kronij(ab,ai))-drija*drikb/(dr1*dr2**3)&
         *(kronij(ab,ak)-kronij(ab,ai)))*(kronij(aa,ak)-kronij(aa,ai))
    func=func-((kronij(xx,xl)/dr1**2-2.d0*drija*drijb/dr1**4)&
         *(kronij(ab,aj)-kronij(ab,ai))*(kronij(aa,aj)-kronij(aa,ai))&
         +(kronij(xx,xl)/dr2**2-2.d0*drika*drikb/dr2**4)&
         *(kronij(ab,ak)-kronij(ab,ai))*(kronij(aa,ak)-kronij(aa,ai)))*cos(theta)
    call this%set_dralpha(ai,aj,ak,ab,drijb,drikb,dr1,dr2,theta)
    func=func-(drija/dr1**2*(kronij(aa,aj)-kronij(aa,ai))+drika/dr2**2&
         *(kronij(aa,ak)-kronij(aa,ai)))*this%get_dralpha()
    this%drbetaalpha=func
  end subroutine set_drbetaalpha

  double precision function get_drbetaalpha(this)
    implicit none
    class(gradientdescent), intent(in) :: this
    get_drbetaalpha=this%drbetaalpha
  end function get_drbetaalpha

  integer function kronij(i,j)
    implicit none
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
