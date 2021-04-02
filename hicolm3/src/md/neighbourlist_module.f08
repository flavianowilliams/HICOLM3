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
module neighbourlist_module
  !*******************************************************************************************
  !*******************************************************************************************

  use atoms_module

  implicit none

  private
  public :: neighbourlist

  type, extends(atoms) :: neighbourlist
     integer, private     :: verlchk
     integer, allocatable :: ilist(:,:)
     integer, allocatable :: nlist(:)
   contains
     procedure :: neighbour_prepare
     procedure :: verlet_list
  end type neighbourlist

contains

  subroutine neighbour_prepare(this)
    implicit none
    class(neighbourlist), intent(inout) :: this
    integer                             :: numb
    numb=this%get_natom()**2
    if(numb.le.0)goto 1
    allocate(this%nlist(this%get_natom()),this%ilist(this%get_natom(),numb))
    this%verlchk=1
    return
1   write(6,*)'ERROR: Fail to allocate the neighbour list array!'
    write(6,*)'Hint: Check if the number of atoms is correct.'
    stop
  end subroutine neighbour_prepare

  subroutine verlet_list(this)
    implicit none
    class(neighbourlist), intent(inout) :: this
    integer                             :: i,j,k,l,m,n,n1,n2,nx
    real(8)                             :: xvz,yvz,zvz,dr
    n1=1
    n2=1+this%nxmol(1)
    do i=1,this%get_nmol()
       do j=1,this%ntmol(i)
          do k=1,this%nxmol(i)
             nx=1
             do l=n2,this%get_natom()
                call this%mic(n1,n2,xvz,yvz,zvz)
                dr=sqrt(xvz**2+yvz**2+zvz**2)
                if(dr.le.this%get_rcutoff())then
                   this%ilist(n1,nx)=l
                   nx=nx+1
                end if
             end do
             this%nlist(n1)=nx-1
             n1=n1+1
          end do
          n2=n2+this%nxmol(i)
       end do
    end do
!    n1=1
!    do i=1,this%get_nmol()
!       do j=1,this%ntmol(i)
!          do k=1,this%nxmol(i)
!             nx=1
    !             n2=j*this%nxmol(i)+1
    !             do l=1,this%get_nmol()
!                do m=j+1,this%ntmol(l)
!                   do n=1,this%nxmol(l)
!                      call this%mic(n1,n2,xvz,yvz,zvz)
!                      dr=sqrt(xvz**2+yvz**2+zvz**2)
!                      if(dr.le.this%get_rcutoff())then
!                         this%ilist(n1,nx)=n2!
!                         nx=nx+1
!                      end if
!                      n2=n2+1
!                   end do
!                end do
!             end do
!             n2=1
!             do l=1,i
!                n2=n2+this%ntmol(l)*this%nxmol(l)
!             end do
!             do l=i+1,this%get_nmol()
!                do m=1,this%ntmol(l)
!                   do n=1,this%nxmol(l)
!                      call this%mic(n1,n2,xvz,yvz,zvz)
!                      dr=sqrt(xvz**2+yvz**2+zvz**2)
!                      if(dr.le.this%get_rcutoff())then
!                         this%ilist(n1,nx)=n2
!                         nx=nx+1
!                      end if
!                      n2=n2+1
!                   end do
!                end do
!             end do
    !             this%nlist(n1)=nx-1
    !             n1=n1+1
!          end do
!       end do
!    end do
  end subroutine verlet_list

end module neighbourlist_module
