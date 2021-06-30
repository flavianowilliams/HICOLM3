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
     procedure :: set_verlchk
     procedure :: get_verlchk
  end type neighbourlist

contains

  subroutine neighbour_prepare(this)
    implicit none
    class(neighbourlist), intent(inout) :: this
    integer                             :: numb
    numb=this%get_natom()**2
    if(numb.le.0)goto 1
    allocate(this%nlist(this%get_natom()),this%ilist(this%get_natom(),numb))
    return
1   write(6,*)'ERROR: Fail to allocate the neighbour list array!'
    write(6,*)'Hint: Check if the number of atoms is correct.'
    stop
  end subroutine neighbour_prepare

  subroutine verlet_list(this)
    implicit none
    class(neighbourlist), intent(inout) :: this
    integer                             :: i,j,k,l,n1,n2,nx
    real(8)                             :: xvz,yvz,zvz,dr
    n1=1
    n2=1+this%nxmol(1)
    do i=1,this%get_nmol()
       do j=1,this%ntmol(i)
          do k=1,this%nxmol(i)
             nx=1
             do l=n2,this%get_natom()
                call this%mic(n1,l,xvz,yvz,zvz)
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
  end subroutine verlet_list

  subroutine set_verlchk(this)
    implicit none
    class(neighbourlist), intent(inout) :: this
    integer                             :: i
    real(8)                             :: drmax,drx,dry,drz
    drmax=0.d0
    do i=1,this%get_natom()
       drx=this%vax(i)*this%get_timestep()&
            +0.5d0*this%fax(i)*this%get_timestep()**2/this%mass(i)
       dry=this%vay(i)*this%get_timestep()&
            +0.5d0*this%fay(i)*this%get_timestep()**2/this%mass(i)
       drz=this%vaz(i)*this%get_timestep()&
            +0.5d0*this%faz(i)*this%get_timestep()**2/this%mass(i)
       drmax=max(drmax,sqrt(drx**2+dry**2+drz**2))
    end do
    if(drmax.ge.1.d-8)this%verlchk=nint(this%get_drcutoff()/drmax)
  end subroutine set_verlchk

  integer function get_verlchk(this)
    implicit none
    class(neighbourlist), intent(in) :: this
    get_verlchk=this%verlchk
  end function get_verlchk

end module neighbourlist_module
