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
module structure_module
  !*******************************************************************************************
  !*******************************************************************************************

  use constants_module

  implicit none

  integer i

  private
  public :: structure

  type, extends(constants) :: structure
     integer, private           :: nmol
     integer, private           :: natom
     character(9), private      :: gsym
     integer, allocatable       :: ntmol(:)
     integer, allocatable       :: nxmol(:)
     character(10), allocatable :: namemol(:)
     real(8), private           :: volume
     real(8), private           :: a
     real(8), private           :: b
     real(8), private           :: c
     real(8), private           :: alpha
     real(8), private           :: beta
     real(8), private           :: gamma
     real(8), allocatable       :: xa(:)
     real(8), allocatable       :: ya(:)
     real(8), allocatable       :: za(:)
     real(8)                    :: v(3,3)
     real(8)                    :: sys_shift(3)
   contains
     procedure, private :: structure_init
     procedure          :: sites
     procedure          :: molecules
     procedure          :: translate
     procedure          :: set_latticevectors
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
     procedure          :: set_sys_shift
     procedure          :: set_volume
     procedure          :: get_volume
     procedure          :: get_a
     procedure          :: get_b
     procedure          :: get_c
     procedure          :: get_alpha
     procedure          :: get_beta
     procedure          :: get_gamma
     procedure          :: get_gsym
     procedure          :: ccp
     procedure          :: mic
  end type structure

  interface structure
     module procedure constructor
  end interface structure

contains

  type(structure) function constructor()
    call constructor%structure_init()
  end function constructor

  subroutine structure_init(this)
    class(structure), intent(inout) :: this
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
  end subroutine structure_init

  subroutine set_sys_shift(this)
    class(structure), intent(inout) :: this
    character(4)                 :: key
    character(9)                 :: key2
    do i=1,3
       this%sys_shift(i)=0.d0
    end do
1   read(5,*,end=2)key
    if(key.ne.'&SYS')goto 1
    do while (key2.ne.'&END     ')
       read(5,*)key2
       if(key2.eq.'translate')then
          backspace(5)
          read(5,*)key,(this%sys_shift(i),i=1,3)
       end if
    end do
2   rewind(5)
  end subroutine set_sys_shift

  subroutine translate(this)
    class(structure), intent(inout) :: this
    call this%set_sys_shift()
    do i=1,this%get_natom()
       this%xa(i)=this%xa(i)+1.0d0*this%sys_shift(1)*this%a
       this%ya(i)=this%ya(i)+1.0d0*this%sys_shift(2)*this%b
       this%za(i)=this%za(i)+1.0d0*this%sys_shift(3)*this%c
    end do
  end subroutine translate

  subroutine set_nmol(this,nmol)
    class(structure), intent(inout) :: this
    integer, intent(in)          :: nmol
    this%nmol=nmol
  end subroutine set_nmol

  integer function get_nmol(this)
    class(structure), intent(in) :: this
    get_nmol=this%nmol
  end function get_nmol

  subroutine set_namemol(this,namemol,i)
    class(structure), intent(inout) :: this
    integer, intent(in)          :: i
    character(10), intent(in)    :: namemol
    this%namemol(i)=namemol
  end subroutine set_namemol

  subroutine set_ntmol(this,ntmol,i)
    class(structure), intent(inout) :: this
    integer, intent(in)          :: i
    integer, intent(in)          :: ntmol
    this%ntmol(i)=ntmol
  end subroutine set_ntmol

  subroutine set_nxmol(this,nxmol,i)
    class(structure), intent(inout) :: this
    integer, intent(in)          :: i
    integer, intent(in)          :: nxmol
    this%nxmol(i)=nxmol
  end subroutine set_nxmol

  subroutine set_volume(this)
    implicit none
    class(structure), intent(inout) :: this
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
    class(structure), intent(in) :: this
    get_volume=this%volume
  end function get_volume

  double precision function get_a(this)
    class(structure), intent(in) :: this
    get_a=this%a
  end function get_a

  double precision function get_b(this)
    class(structure), intent(in) :: this
    get_b=this%b
  end function get_b

  double precision function get_c(this)
    class(structure), intent(in) :: this
    get_c=this%c
  end function get_c

  double precision function get_alpha(this)
    class(structure), intent(in) :: this
    get_alpha=this%alpha
  end function get_alpha

  double precision function get_beta(this)
    class(structure), intent(in) :: this
    get_beta=this%beta
  end function get_beta

  double precision function get_gamma(this)
    class(structure), intent(in) :: this
    get_gamma=this%gamma
  end function get_gamma

  character(9) function get_gsym(this)
    class(structure), intent(in) :: this
    get_gsym=this%gsym
  end function get_gsym

  subroutine molecules(this)
    class(structure), intent(inout) :: this
    integer                      :: nx
    character(4)                 :: key
    call this%structure_init()
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
    class(structure), intent(inout) :: this
    integer                      :: nx
    nx=0
    do i=1,this%nmol
       nx=nx+this%ntmol(i)*this%nxmol(i)
    end do
    this%natom=nx
  end subroutine set_natom

  function get_natom(this)
    class(structure), intent(inout) :: this
    integer                      :: get_natom
    get_natom=this%natom
  end function get_natom

  subroutine sites(this)
    class(structure), intent(inout) :: this
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
         'ERROR: The number of atoms in HICOLM.xyz does not match with defined in INPUT file!'
    write(6,*)'Hint: Change the HICOLM.xyz, or the inputs in INPUT file.'
    stop
  end subroutine sites

  subroutine set_latticevectors(this)
    class(structure), intent(inout) :: this
    character(4)                    :: key
1   read(5,*,end=2)key
    if(key.ne.'&SYS')goto 1
    do while (key.ne.'&END')
       read(5,*)key
       if(key.eq.'cell')then
          backspace(5)
          read(5,*)key,this%v(1,1),this%v(2,2),this%v(3,3)
       end if
    end do
2   rewind(5)
  end subroutine set_latticevectors

  subroutine set_lattice_constants(this)
    class(structure), intent(inout) :: this
    this%a=sqrt(this%v(1,1)**2+this%v(1,2)**2+this%v(1,3)**2)
    this%b=sqrt(this%v(2,1)**2+this%v(2,2)**2+this%v(2,3)**2)
    this%c=sqrt(this%v(3,1)**2+this%v(3,2)**2+this%v(3,3)**2)
  end subroutine set_lattice_constants

  subroutine set_lattice_angles(this)
    class(structure), intent(inout) :: this
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
    class(structure), intent(inout) :: this
    real(8)                      :: prec
    character(9)                 :: cvar
    prec=1.d-2
    !tryclinic
    cvar='Tryclinic'
    !-cubic
    if(abs(2.d0*this%alpha-this%get_pi()).le.prec)then
       if(abs(2.d0*this%beta-this%get_pi()).le.prec)then
          if(abs(2.d0*this%gamma-this%get_pi()).le.prec)then
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
             if(abs(2.d0*this%gamma-this%get_pi()/6.d0).le.prec.or.&
                  abs(2.d0*this%gamma-this%get_pi()/3.d0).le.prec)cvar='Hexagonal'
          end if
       end if
    end if
    this%gsym=cvar
  end subroutine set_symmetry

  subroutine ccp(this)
    implicit none
    class(structure), intent(inout) :: this
    integer                         :: i
    real(8)                         :: xvz,yvz,zvz
    do i=1,this%get_natom()
       xvz=this%v(1,1)+this%v(2,1)+this%v(3,1)
       yvz=this%v(1,2)+this%v(2,2)+this%v(3,2)
       zvz=this%v(1,3)+this%v(2,3)+this%v(3,3)
       this%xa(i)=this%xa(i)-xvz*nint(this%xa(i)/xvz)
       this%ya(i)=this%ya(i)-yvz*nint(this%ya(i)/yvz)
       this%za(i)=this%za(i)-zvz*nint(this%za(i)/zvz)
    end do
  end subroutine ccp

  subroutine mic(this,i1,i2,xvz,yvz,zvz)
    implicit none
    class(structure), intent(inout) :: this
    integer, intent(in)             :: i1,i2
    real(8), intent(out)            :: xvz,yvz,zvz
    real(8)                         :: xx,yy,zz
    xx=this%v(1,1)+this%v(2,1)+this%v(3,1)
    yy=this%v(1,2)+this%v(2,2)+this%v(3,2)
    zz=this%v(1,3)+this%v(2,3)+this%v(3,3)
    xvz=(this%xa(i2)-this%xa(i1))-xx*int(2.d0*(this%xa(i2)-this%xa(i1))/xx)
    yvz=(this%ya(i2)-this%ya(i1))-yy*int(2.d0*(this%ya(i2)-this%ya(i1))/yy)
    zvz=(this%za(i2)-this%za(i1))-zz*int(2.d0*(this%za(i2)-this%za(i1))/zz)
  end subroutine mic

end module structure_module
