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
module zmatrix_module
  !*******************************************************************************************
  !*******************************************************************************************

  use molecule_module

  implicit none

  integer i,j

  private
  public :: zmatrix

  type, extends(molecule) :: zmatrix
     real(8)              :: zmatrix_tol
     integer, private     :: bondmax
     integer, allocatable :: bondscnt(:)
     integer, allocatable :: bendscnt(:)
     integer, allocatable :: torscnt(:)
     integer, allocatable :: itorscnt(:)
     integer, allocatable :: molbond(:,:,:)
   contains
     procedure, private :: zmatrix_init
     procedure          :: set_bonds
     procedure          :: set_bends
     procedure          :: set_torsion
     procedure          :: covalent_radius
     procedure          :: set_internal_coordinates
     procedure          :: set_bondmax
     procedure          :: get_bondmax
  end type zmatrix

  interface zmatrix
     module procedure constructor
  end interface zmatrix

contains

  type(zmatrix) function constructor()
    call constructor%zmatrix_init()
  end function constructor

  subroutine zmatrix_init(this)
    class(zmatrix), intent(inout) :: this
    call this%set_bondmax()
    allocate(this%bondscnt(this%get_nmol()))
    allocate(this%molbond(this%get_nmol(),this%bondmax,2))
  end subroutine zmatrix_init

  subroutine set_internal_coordinates(this)
    class(zmatrix), intent(inout) :: this
    call this%zmatrix_init()
    call this%set_bonds()
  end subroutine set_internal_coordinates

  subroutine set_bonds(this)
    class(zmatrix), intent(inout) :: this
    integer                       :: nx,nxx,imol,ia,ib
    real(8)                       :: dr,rca,rcb
    do imol=1,this%get_nmol()
       nx=0
       do i=1,imol-1
          nx=nx+this%nxmol(i)*this%ntmol(i)
       end do
       nxx=1
       do i=1,this%nxmol(imol)
          ia=i+nx
          do j=i+1,this%nxmol(imol)
             ib=j+nx
             dr=sqrt((this%xa(ia)-this%xa(ib))**2+(this%ya(ia)-this%ya(ib))**2+&
                  (this%za(ia)-this%za(ib))**2)
             call this%covalent_radius(imol,i,rca)
             call this%covalent_radius(imol,j,rcb)
             if(dr.gt.(rca+rcb-this%zmatrix_tol).and.dr.le.(rca+rcb+this%zmatrix_tol))then
                this%molbond(imol,nxx,1)=i
                this%molbond(imol,nxx,2)=j
                nxx=nxx+1
             end if
          end do
       end do
       this%bondscnt(imol)=nxx-1
    end do
  end subroutine set_bonds

  subroutine set_bends(this)
    class(zmatrix), intent(inout) :: this
    allocate(this%bendscnt(this%get_nmol()))
  end subroutine set_bends

  subroutine set_torsion(this)
    class(zmatrix), intent(inout) :: this
    allocate(this%torscnt(this%get_nmol()))
  end subroutine set_torsion

  subroutine set_zmatrix_tol(this)
    class(zmatrix), intent(inout) :: this
    

  subroutine set_bondmax(this)
    class(zmatrix), intent(inout) :: this
    integer                       :: nx,nxx
    nxx=1
    do i=1,this%get_nmol()
       nx=1
       do j=0,this%nxmol(i)-1
          nx=nx*(this%nxmol(i)-j)
       end do
       nx=int(0.5*nx)
       nxx=max(nx,nxx)
    end do
    this%bondmax=nxx
  end subroutine set_bondmax

  integer function get_bondmax(this)
    class(zmatrix), intent(in) :: this
    get_bondmax=this%bondmax
  end function get_bondmax

  subroutine covalent_radius(this,ix,jx,rc)
    class(zmatrix), intent(inout) :: this
    integer, intent(in)            :: ix,jx
    real(8), intent(out)           :: rc
    select case(this%zatmol(ix,jx))
    case(1)
       rc=0.37d0
    case(2)
       rc=0.32d0
    case(3)
       rc=1.34d0
    case(4)
       rc=0.90d0
    case(5)
       rc=0.82d0
    case(6)
       rc=0.77d0
    case(7)
       rc=0.75d0
    case(8)
       rc=0.73d0
    case(9)
       rc=0.71d0
    case(10)
       rc=0.69d0
    case(11)
       rc=1.54d0
    case(12)
       rc=1.30d0
    case(13)
       rc=1.18d0
    case(14)
       rc=1.11d0
    case(15)
       rc=1.06d0
    case(16)
       rc=1.02d0
    case(17)
       rc=0.99d0
    case(18)
       rc=0.97d0
    case(19)
       rc=1.96d0
    case(20)
       rc=1.74d0
    case(21)
       rc=1.44d0
    case(22)
       rc=1.36d0
    case(23)
       rc=1.25d0
    case(24)
       rc=1.27d0
    case(25)
       rc=1.39d0
    case(26)
       rc=1.25d0
    case(27)
       rc=1.26d0
    case(28)
       rc=1.21d0
    case(29)
       rc=1.38d0
    case(30)
       rc=1.31d0
    case(31)
       rc=1.26d0
    case(32)
       rc=1.22d0
    case(33)
       rc=1.19d0
    case(34)
       rc=1.16d0
    case(35)
       rc=1.14d0
    case(36)
       rc=1.10d0
    case(37)
       rc=2.11d0
    case(38)
       rc=1.92d0
    case(39)
       rc=1.62d0
    case(40)
       rc=1.48d0
    case(41)
       rc=1.37d0
    case(42)
       rc=1.45d0
    case(43)
       rc=1.56d0
    case(44)
       rc=1.26d0
    case(45)
       rc=1.35d0
    case(46)
       rc=1.31d0
    case(47)
       rc=1.53d0
    case(48)
       rc=1.48d0
    case(49)
       rc=1.44d0
    case(50)
       rc=1.41d0
    case(51)
       rc=1.38d0
    case(52)
       rc=1.35d0
    case(53)
       rc=1.33d0
    case(54)
       rc=1.30d0
    case(55)
       rc=2.25d0
    case(56)
       rc=1.98d0
    case(57)
       rc=1.69d0
    case(71)
       rc=1.60d0
    case(72)
       rc=1.50d0
    case(73)
       rc=1.38d0
    case(74)
       rc=1.46d0
    case(75)
       rc=1.59d0
    case(76)
       rc=1.28d0
    case(77)
       rc=1.37d0
    case(78)
       rc=1.28d0
    case(79)
       rc=1.44d0
    case(80)
       rc=1.49d0
    case(81)
       rc=1.48d0
    case(82)
       rc=1.47d0
    case(83)
       rc=1.46d0
    end select
  end subroutine covalent_radius

end module zmatrix_module
