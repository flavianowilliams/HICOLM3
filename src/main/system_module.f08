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
module system_module
  !*******************************************************************************************
  !*******************************************************************************************

  use constants_module

  implicit none

  integer i

  private
  public :: system

  type, extends(constants) :: system
     integer                    :: nmol
     integer                    :: natom
     character(9)               :: gsym
     integer, allocatable       :: ntmol(:)
     integer, allocatable       :: nxmol(:)
     character(10), allocatable :: namemol(:)
     real(8)                    :: volume
     real(8)                    :: a
     real(8)                    :: b
     real(8)                    :: c
     real(8)                    :: alpha
     real(8)                    :: beta
     real(8)                    :: gamma
     real(8)                    :: v(3,3)
     real(8), allocatable       :: xa(:)
     real(8), allocatable       :: ya(:)
     real(8), allocatable       :: za(:)
   contains
     procedure, private :: system_init
     procedure          :: sites
     procedure          :: molecules
     procedure          :: set_lattice_constants
     procedure          :: set_lattice_angles
     procedure          :: set_symmetry
     procedure          :: set_nmol
     procedure          :: get_nmol
     procedure          :: set_natom
     procedure          :: get_natom
     procedure          :: set_namemol
     procedure          :: set_ntmol
     procedure          :: set_nxmol
     procedure          :: set_volume
     procedure          :: get_volume
     procedure          :: get_a
     procedure          :: get_b
     procedure          :: get_c
  end type system

  interface system
     module procedure constructor
  end interface system

contains

  type(system) function constructor()
    call constructor%system_init()
  end function constructor

  subroutine system_init(this)
    class(system), intent(inout) :: this
    integer                      :: nx
    character(4)                 :: key
    nx=0
1   read(5,*,end=2)key
    if(key.ne.'&SYS')goto 1
    do while (key.ne.'&END')
       read(5,*)key
       if(key.eq.'nmol')then
          backspace(5)
          read(5,*)key,nx
       end if
    end do
    allocate(this%namemol(nx),this%ntmol(nx),this%nxmol(nx))
    call this%set_nmol(nx)
2   rewind(5)
  end subroutine system_init

  subroutine set_nmol(this,nmol)
    class(system), intent(inout) :: this
    integer, intent(in)          :: nmol
    this%nmol=nmol
  end subroutine set_nmol

  integer function get_nmol(this)
    class(system), intent(in) :: this
    get_nmol=this%nmol
  end function get_nmol

  subroutine set_namemol(this,namemol,i)
    class(system), intent(inout) :: this
    integer, intent(in)          :: i
    character(10), intent(in)    :: namemol
    this%namemol(i)=namemol
  end subroutine set_namemol

  subroutine set_ntmol(this,ntmol,i)
    class(system), intent(inout) :: this
    integer, intent(in)          :: i
    integer, intent(in)          :: ntmol
    this%ntmol(i)=ntmol
  end subroutine set_ntmol

  subroutine set_nxmol(this,nxmol,i)
    class(system), intent(inout) :: this
    integer, intent(in)          :: i
    integer, intent(in)          :: nxmol
    this%nxmol(i)=nxmol
  end subroutine set_nxmol

  subroutine set_volume(this)
    implicit none
    class(system), intent(inout) :: this
    real(8)                      :: volume,vl(3)
    vl(1)=(this%v(1,2)*this%v(2,3)-this%v(1,3)*this%v(2,2))
    vl(2)=(this%v(1,3)*this%v(2,1)-this%v(1,1)*this%v(2,3))
    vl(3)=(this%v(1,1)*this%v(2,2)-this%v(1,2)*this%v(2,1))
    volume=0.d0
    do i=1,3
       volume=volume+this%v(3,i)*vl(i)
    end do
    this%volume=abs(volume)
  end subroutine set_volume

  double precision function get_volume(this)
    class(system), intent(in) :: this
    get_volume=this%volume
  end function get_volume

  double precision function get_a(this)
    class(system), intent(in) :: this
    get_a=this%a
  end function get_a

  double precision function get_b(this)
    class(system), intent(in) :: this
    get_b=this%b
  end function get_b

  double precision function get_c(this)
    class(system), intent(in) :: this
    get_c=this%c
  end function get_c

  subroutine molecules(this)
    class(system), intent(inout) :: this
    integer                      :: nx
    character(4)                 :: key
    call this%system_init()
    nx=0
1   read(5,*,end=2)key
    if(key.ne.'&SYS')goto 1
    do while (key.ne.'&END')
       read(5,*)key
       if(key.eq.'nmol')then
          backspace(5)
          read(5,*)key,nx
          do i=1,nx
             read(5,*)this%namemol(i),this%ntmol(i),this%nxmol(i)
          end do
       end if
    end do
2   rewind(5)
  end subroutine molecules

  subroutine set_natom(this)
    implicit none
    class(system), intent(inout) :: this
    integer                      :: nx
    nx=0
    do i=1,this%nmol
       nx=nx+this%ntmol(i)*this%nxmol(i)
    end do
    this%natom=nx
  end subroutine set_natom

  function get_natom(this)
    implicit none
    class(system), intent(inout) :: this
    integer                      :: get_natom
    get_natom=this%natom
  end function get_natom

  subroutine sites(this)
    implicit none
    class(system), intent(inout) :: this
    integer                        :: nx
    character(2)                   :: at
    nx=this%get_natom()
    allocate(this%xa(nx),this%ya(nx),this%za(nx))
    read(9,*)nx
    read(9,*)
    if(nx.ne.this%get_natom())goto 1
    do i=1,this%get_natom()
       read(9,*)at,this%xa(i),this%ya(i),this%za(i)
    end do
    return
1   write(6,*)&
         'ERROR: The number of atoms in HICOLM.xyz does not match with defined in HICOLM.in!'
    write(6,*)'Hint: Change the HICOLM.xyz, or the inputs in HICOLM.in.'
    stop
  end subroutine sites

  subroutine set_lattice_constants(this)
    implicit none
    class(system), intent(inout) :: this
    character(4)                 :: key
1   read(5,*,end=2)key
    if(key.ne.'&SYS')goto 1
    do while (key.ne.'&END')
       read(5,*)key
       if(key.eq.'cell')then
          backspace(5)
          read(5,*)key,this%v(1,1),this%v(2,2),this%v(3,3)
       end if
    end do
    this%a=sqrt(this%v(1,1)**2+this%v(1,2)**2+this%v(1,3)**2)
    this%b=sqrt(this%v(2,1)**2+this%v(2,2)**2+this%v(2,3)**2)
    this%c=sqrt(this%v(3,1)**2+this%v(3,2)**2+this%v(3,3)**2)
2   rewind(5)
  end subroutine set_lattice_constants

  subroutine set_lattice_angles(this)
    implicit none
    class(system), intent(inout) :: this
    real(8)                      :: sum
    sum=0.d0
    do i=1,3
       sum=sum+this%v(2,i)*this%v(3,i)
    end do
    this%alpha=acos(sum/(this%b*this%c))
    sum=0.d0
    do i=1,3
       sum=sum+this%v(1,i)*this%v(3,i)
    end do
    this%beta=acos(sum/(this%a*this%c))
    sum=0.d0
    do i=1,3
       sum=sum+this%v(1,i)*this%v(2,i)
    end do
    this%gamma=acos(sum/(this%a*this%b))
  end subroutine set_lattice_angles

  subroutine set_symmetry(this)
    implicit none
    class(system), intent(inout) :: this
    real(8)                      :: prec
    character(9)                 :: cvar
    prec=1.d-2
    !tryclinic
    cvar='Tryclinic'
    !-cubic
    if(abs(2.d0*this%alpha-this%pi).le.prec)then
       if(abs(2.d0*this%beta-this%pi).le.prec)then
          if(abs(2.d0*this%gamma-this%pi).le.prec)then
             if(abs(this%a-this%b).le.prec)then
                if(abs(this%a-this%c).le.prec)cvar='Cubic'
             end if
          end if
       end if
    end if
    !-Hexagonal
    if(abs(2.d0*this%alpha-acos(-1.d0)).le.prec)then
       if(abs(2.d0*this%beta-acos(-1.d0)).le.prec)then
          if(abs(this%a-this%b).le.prec)then
             if(abs(2.d0*this%gamma-this%pi/6.d0).le.prec.or.&
                  abs(2.d0*this%gamma-this%pi/3.d0).le.prec)cvar='Hexagonal'
          end if
       end if
    end if
    this%gsym=cvar
  end subroutine set_symmetry

end module system_module
