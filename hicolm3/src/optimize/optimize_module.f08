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
module optimize_module
  !*******************************************************************************************
  !*******************************************************************************************

  use input_module

  implicit none

  private
  public :: optimize

  type, extends(input) :: optimize
   contains
     procedure :: print
     procedure :: convert_units
  end type optimize

  interface optimize
     module procedure constructor
  end interface optimize

contains

  type(optimize) function constructor()
    implicit none
    call constructor%set_nstep(1000)
    call constructor%set_tolerance(0.001d0)
  end function constructor

  subroutine convert_units(this)
    implicit none
    class(optimize), intent(inout) :: this
    integer                        :: i,j
    do i=1,3
       do j=1,3
!          this%v(i,j)=this%v(i,j)/this%get_rconv()
       end do
    end do
!    call this%set_a(this%get_a()/this%get_rconv())
!    call this%set_b(this%get_b()/this%get_rconv())
!    call this%set_c(this%get_c()/this%get_rconv())
!    do i=1,this%get_natom()
!       this%xa(i)=this%xa(i)/this%get_rconv()
!       this%ya(i)=this%ya(i)/this%get_rconv()
!       this%za(i)=this%za(i)/this%get_rconv()
!    end do
!    do i=1,this%get_nmol()
!       do j=1,this%nxmol(i)
!          this%qatmol(i,j)=this%qatmol(i,j)/this%get_elconv()
!       end do
!    end do
!    do i=1,this%get_nmol()
!       do j=1,this%bondscnt(i)
!          select case(this%tbonds(i,j))
!          case('charmm')
!             this%parbnd(i,j,1)=this%parbnd(i,j,1)/(this%get_econv()/this%get_rconv()**2)
!             this%parbnd(i,j,2)=this%parbnd(i,j,2)/this%get_rconv()
!          case('harm')
!             this%parbnd(i,j,1)=this%parbnd(i,j,1)/(this%get_econv()/this%get_rconv()**2)
!             this%parbnd(i,j,2)=this%parbnd(i,j,2)/this%get_rconv()
!          end select
!       end do
!       do j=1,this%bendscnt(i)
!          select case(this%tbends(i,j))
!          case('charmm')
!             this%parbend(i,j,1)=this%parbend(i,j,1)/this%get_econv()
!             this%parbend(i,j,2)=this%parbend(i,j,2)/this%get_aconv()
!          case('harm')
!             this%parbend(i,j,1)=this%parbend(i,j,1)/this%get_econv()
!             this%parbend(i,j,2)=this%parbend(i,j,2)/this%get_aconv()
!          end select
!       end do
!       do j=1,this%torscnt(i)
!          select case(this%ttors(i,j))
!          case('charmm')
!             this%partors(i,j,1)=this%partors(i,j,1)/this%get_econv()
!             this%partors(i,j,3)=this%partors(i,j,3)/this%get_aconv()
!          case('icharmm')
!             this%partors(i,j,1)=this%partors(i,j,1)/this%get_econv()
!             this%partors(i,j,2)=this%partors(i,j,2)/this%get_aconv()
!          case('harm')
!             this%partors(i,j,1)=this%partors(i,j,1)/this%get_econv()
!             this%partors(i,j,2)=this%partors(i,j,2)/this%get_aconv()
!          end select
!       end do
!    end do
!    do i=1,this%get_nvdw()
!       select case(this%tvdw(i))
!       case('charmm')
!          this%parvdw(i,1)=this%parvdw(i,1)/this%get_econv()
!          this%parvdw(i,2)=this%parvdw(i,2)/this%get_rconv()
!       case('lj')
!          this%parvdw(i,1)=this%parvdw(i,1)/this%get_econv()
!          this%parvdw(i,2)=this%parvdw(i,2)/this%get_rconv()
!       end select
!    end do
    call this%set_tolerance(this%get_tolerance()/(this%get_econv()/this%get_rconv()))
  end subroutine convert_units

  subroutine print(this)
    implicit none
    class(optimize), intent(inout) :: this
    integer                        :: i,j
    write(6,*)('#',i=1,93)
    write(6,*)(' MINIMIZATION ',i=1,6)
    write(6,*)('#',i=1,93)
    write(6,*)
    write(6,'(28x,a38)')'Minimization of the energy information'
    write(6,'(28x,38a1)')('-',j=1,38)
    write(6,'(28x,a16,1x,i10)')'Number of steps:',this%get_nstep()
    write(6,'(28x,a16,1x,es10.3,1x,a10)')'      Tolerance:',&
         this%get_tolerance()*(this%get_econv()/this%get_rconv()),'kcal/mol*A'
    write(6,'(28x,38a1)')('-',j=1,38)
    write(6,*)
    write(6,*)('#',i=1,93)
    write(6,*)('SYSTEM ',i=1,13)
    write(6,*)('#',i=1,93)
    write(6,*)
    write(6,'(a18,i5)')'Total of atoms:',this%get_natom()
    write(6,*)
    write(6,*)'Real space:'
    write(6,*)
    write(6,'(a16,3f15.8)')'Lattice constts:',this%get_a()*this%get_rconv(),&
         this%get_b()*this%get_rconv(),this%get_c()*this%get_rconv()
    write(6,'(a16,3f15.8)')' Lattice angles:',this%get_alpha()*this%get_aconv(),&
         this%get_beta()*this%get_aconv(),this%get_gamma()*this%get_aconv()
    write(6,*)
    write(6,'(a16,3f15.8)')'Lattice vectors:',(this%v(1,i)*this%get_rconv(),i=1,3)
    write(6,'(16x,3f15.8)')(this%v(2,i)*this%get_rconv(),i=1,3)
    write(6,'(16x,3f15.8)')(this%v(3,i)*this%get_rconv(),i=1,3)
    write(6,*)
    write(6,'(a14,f14.4)')'       VOLUME:',this%get_volume()*this%get_rconv()**3
    write(6,*)
    write(6,'(a14,1x,a9)')'     Symmetry:',this%get_gsym()
    write(6,*)
  end subroutine print

end module optimize_module
