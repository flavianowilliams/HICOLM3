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

  integer i,j,k,l

  private
  public :: zmatrix

  type, extends(molecule) :: zmatrix
     real(8), private     :: zmatrix_tol
     integer, private     :: bondmax
     integer, private     :: bendmax
     integer, private     :: torsmax
     integer, private     :: itorsmax
     integer, allocatable :: bondscnt(:)
     integer, allocatable :: bendscnt(:)
     integer, allocatable :: torscnt(:)
     integer, allocatable :: molbond(:,:,:)
     integer, allocatable :: molbend(:,:,:)
     integer, allocatable :: moltors(:,:,:)
   contains
     procedure, private :: zmatrix_init
     procedure          :: set_bonds
     procedure          :: set_bends
     procedure          :: set_torsion
     procedure          :: covalent_radius
     procedure          :: set_internal_coordinates
     procedure          :: set_bondmax
     procedure          :: get_bondmax
     procedure          :: set_bendmax
     procedure          :: get_bendmax
     procedure          :: set_torsmax
     procedure          :: get_torsmax
     procedure          :: set_itorsmax
     procedure          :: get_itorsmax
     procedure          :: set_zmatrix_tol
     procedure          :: get_zmatrix_tol
     procedure          :: set_zmatrixtol
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
    call this%set_zmatrix_tol()
  end subroutine zmatrix_init

  subroutine set_internal_coordinates(this)
    implicit none
    class(zmatrix), intent(inout) :: this
    call this%zmatrix_init()
    call this%set_bonds()
    call this%set_bends()
    call this%set_torsion()
  end subroutine set_internal_coordinates

  subroutine set_bonds(this)
    implicit none
    class(zmatrix), intent(inout) :: this
    integer                       :: nx,nxx,nxxx,imol,ia,ib
    real(8)                       :: dr,rca,rcb
    allocate(this%bondscnt(this%get_nmol()))
    allocate(this%molbond(this%get_nmol(),this%bondmax,2))
    nxxx=0
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
       nxxx=max(nxxx,this%bondscnt(imol))
    end do
    this%bondmax=nxxx
  end subroutine set_bonds

  subroutine set_bondmax(this,bondmax)
    implicit none
    class(zmatrix), intent(inout) :: this
    integer, intent(in)           :: bondmax
    this%bondmax=bondmax
  end subroutine set_bondmax

  integer function get_bondmax(this)
    implicit none
    class(zmatrix), intent(in) :: this
    get_bondmax=this%bondmax
  end function get_bondmax

  subroutine set_bends(this)
    implicit none
    class(zmatrix), intent(inout) :: this
    integer                       :: nx,nxx,imol
    allocate(this%bendscnt(this%get_nmol()))
    allocate(this%molbend(this%get_nmol(),this%bendmax,3))
    nxx=0
    do imol=1,this%get_nmol()
       nx=1
       do i=1,this%bondscnt(imol)
          do j=i+1,this%bondscnt(imol)
             do k=1,2
                do l=1,2
                   if(this%molbond(imol,i,k).eq.this%molbond(imol,j,l))then
                      this%molbend(imol,nx,1)=this%molbond(imol,i,3-k)
                      this%molbend(imol,nx,2)=this%molbond(imol,i,k)
                      this%molbend(imol,nx,3)=this%molbond(imol,j,3-l)
                      nx=nx+1
                   end if
                end do
             end do
          end do
       end do
       this%bendscnt(imol)=nx-1
       nxx=max(nxx,this%bendscnt(imol))
    end do
    this%bendmax=nxx
  end subroutine set_bends

  subroutine set_bendmax(this,bendmax)
    implicit none
    class(zmatrix), intent(inout) :: this
    integer, intent(in)           :: bendmax
    this%bendmax=bendmax
  end subroutine set_bendmax

  integer function get_bendmax(this)
    implicit none
    class(zmatrix), intent(in) :: this
    get_bendmax=this%bendmax
  end function get_bendmax

  subroutine set_torsion(this)
    implicit none
    class(zmatrix), intent(inout) :: this
    integer                       :: nx,nxx,imol,ii,kk
    allocate(this%torscnt(this%get_nmol()))
    allocate(this%moltors(this%get_nmol(),this%torsmax,4))
    nxx=0
    do imol=1,this%get_nmol()
       nx=1
       do i=1,this%bondscnt(imol)
          do j=1,i-1
             do k=1,j-1
                do ii=1,2
                   do kk=1,2
                      if(this%molbond(imol,j,1).eq.this%molbond(imol,i,3-ii))then
                         if(this%molbond(imol,j,2).eq.this%molbond(imol,k,kk))then
                            this%moltors(imol,nx,1)=this%molbond(imol,i,ii)
                            this%moltors(imol,nx,2)=this%molbond(imol,j,1)
                            this%moltors(imol,nx,3)=this%molbond(imol,j,2)
                            this%moltors(imol,nx,4)=this%molbond(imol,k,3-kk)
                            nx=nx+1
                         end if
                      end if
                   end do
                end do
             end do
             do k=j+1,this%bondscnt(imol)
                do ii=1,2
                   do kk=1,2
                      if(this%molbond(imol,j,1).eq.this%molbond(imol,i,3-ii))then
                         if(this%molbond(imol,j,2).eq.this%molbond(imol,k,kk))then
                            this%moltors(imol,nx,1)=this%molbond(imol,i,ii)
                            this%moltors(imol,nx,2)=this%molbond(imol,j,1)
                            this%moltors(imol,nx,3)=this%molbond(imol,j,2)
                            this%moltors(imol,nx,4)=this%molbond(imol,k,3-kk)
                            nx=nx+1
                         end if
                      end if
                   end do
                end do
             end do
          end do
          do j=i+1,this%bondscnt(imol)
             do k=1,j-1
                do ii=1,2
                   do kk=1,2
                      if(this%molbond(imol,j,1).eq.this%molbond(imol,i,3-ii))then
                         if(this%molbond(imol,j,2).eq.this%molbond(imol,k,kk))then
                            this%moltors(imol,nx,1)=this%molbond(imol,i,ii)
                            this%moltors(imol,nx,2)=this%molbond(imol,j,1)
                            this%moltors(imol,nx,3)=this%molbond(imol,j,2)
                            this%moltors(imol,nx,4)=this%molbond(imol,k,3-kk)
                            nx=nx+1
                         end if
                      end if
                   end do
                end do
             end do
             do k=j+1,this%bondscnt(imol)
                do ii=1,2
                   do kk=1,2
                      if(this%molbond(imol,j,1).eq.this%molbond(imol,i,3-ii))then
                         if(this%molbond(imol,j,2).eq.this%molbond(imol,k,kk))then
                            this%moltors(imol,nx,1)=this%molbond(imol,i,ii)
                            this%moltors(imol,nx,2)=this%molbond(imol,j,1)
                            this%moltors(imol,nx,3)=this%molbond(imol,j,2)
                            this%moltors(imol,nx,4)=this%molbond(imol,k,3-kk)
                            nx=nx+1
                         end if
                      end if
                   end do
                end do
             end do
          end do
       end do
       this%torscnt(imol)=nx-1
       nxx=max(nxx,this%torscnt(imol))
    end do
  end subroutine set_torsion

  subroutine set_torsmax(this,torsmax)
    implicit none
    class(zmatrix), intent(inout) :: this
    integer, intent(in)           :: torsmax
    this%torsmax=torsmax
  end subroutine set_torsmax

  integer function get_torsmax(this)
    implicit none
    class(zmatrix), intent(in) :: this
    get_torsmax=this%torsmax
  end function get_torsmax

  subroutine set_itorsmax(this,itorsmax)
    implicit none
    class(zmatrix), intent(inout) :: this
    integer, intent(in)           :: itorsmax
    this%itorsmax=itorsmax
  end subroutine set_itorsmax

  integer function get_itorsmax(this)
    implicit none
    class(zmatrix), intent(in) :: this
    get_itorsmax=this%itorsmax
  end function get_itorsmax

  subroutine set_zmatrixtol(this,zmatrix_tol)
    implicit none
    class(zmatrix), intent(inout) :: this
    real(8), intent(in)           :: zmatrix_tol
    this%zmatrix_tol=zmatrix_tol
  end subroutine set_zmatrixtol

  subroutine set_zmatrix_tol(this)
    implicit none
    class(zmatrix), intent(inout) :: this
    integer                       :: nx
    character(7)                  :: key
    nx=0
1   read(5,*,end=2)key
    if(key.ne.'&SYS')goto 1
    do while (key.ne.'&END')
       read(5,*)key
       if(key.eq.'zmatrix')then
          backspace(5)
          read(5,*)key,this%zmatrix_tol
       end if
    end do
2   rewind(5)
  end subroutine set_zmatrix_tol

  double precision function get_zmatrix_tol(this)
    implicit none
    class(zmatrix), intent(in) :: this
    get_zmatrix_tol=this%zmatrix_tol
  end function get_zmatrix_tol

  subroutine covalent_radius(this,ix,jx,rc)
    implicit none
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
