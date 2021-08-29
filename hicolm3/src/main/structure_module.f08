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
     integer, private           :: nfree
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
     procedure          :: sites
     procedure          :: translate
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
     procedure          :: set_volume2
     procedure          :: set_nfree
     procedure          :: get_nfree
     procedure          :: set_a
     procedure          :: get_a
     procedure          :: set_b
     procedure          :: get_b
     procedure          :: set_c
     procedure          :: get_c
     procedure          :: get_alpha
     procedure          :: get_beta
     procedure          :: get_gamma
     procedure          :: get_gsym
     procedure          :: ccp
     procedure          :: mic
     procedure          :: random_coordinates
  end type structure

contains

  subroutine translate(this)
    class(structure), intent(inout) :: this
    do i=1,this%get_natom()
       this%xa(i)=this%xa(i)+1.0d0*this%sys_shift(1)*this%a
       this%ya(i)=this%ya(i)+1.0d0*this%sys_shift(2)*this%b
       this%za(i)=this%za(i)+1.0d0*this%sys_shift(3)*this%c
    end do
  end subroutine translate

  subroutine set_nmol(this,nmol)
    class(structure), intent(inout) :: this
    integer, intent(in)             :: nmol
    this%nmol=nmol
  end subroutine set_nmol

  integer function get_nmol(this)
    class(structure), intent(in) :: this
    get_nmol=this%nmol
  end function get_nmol

  subroutine set_namemol(this,namemol)
    class(structure), intent(inout) :: this
    integer                         :: i
    character(10), intent(in)       :: namemol(:)
    do i=1,this%get_nmol()
       this%namemol(i)=namemol(i)
    end do
  end subroutine set_namemol

  subroutine set_ntmol(this,ntmol)
    class(structure), intent(inout) :: this
    integer                         :: i
    integer, intent(in)             :: ntmol(:)
    do i=1,this%get_nmol()
       this%ntmol(i)=ntmol(i)
    end do
  end subroutine set_ntmol

  subroutine set_nxmol(this,nxmol)
    class(structure), intent(inout) :: this
    integer                         :: i
    integer, intent(in)             :: nxmol(:)
    do i=1,this%get_nmol()
       this%nxmol(i)=nxmol(i)
    end do
  end subroutine set_nxmol

  subroutine set_volume(this)
    implicit none
    class(structure), intent(inout) :: this
    real(8)                         :: volume,vl(3)
    vl(1)=(this%v(1,2)*this%v(2,3)-this%v(1,3)*this%v(2,2))
    vl(2)=(this%v(1,3)*this%v(2,1)-this%v(1,1)*this%v(2,3))
    vl(3)=(this%v(1,1)*this%v(2,2)-this%v(1,2)*this%v(2,1))
    volume=0.d0
    do i=1,3
       volume=volume+this%v(3,i)*vl(i)
    end do
    this%volume=abs(volume)
  end subroutine set_volume

  subroutine set_volume2(this,volume)
    implicit none
    class(structure), intent(inout) :: this
    real(8), intent(in)             :: volume
    this%volume=volume
  end subroutine set_volume2

  double precision function get_volume(this)
    class(structure), intent(in) :: this
    get_volume=this%volume
  end function get_volume

  subroutine set_a(this,a)
    implicit none
    class(structure), intent(inout) :: this
    real(8), intent(in)             :: a
    this%a=a
  end subroutine set_a

  double precision function get_a(this)
    class(structure), intent(in) :: this
    get_a=this%a
  end function get_a

  subroutine set_b(this,b)
    implicit none
    class(structure), intent(inout) :: this
    real(8), intent(in)             :: b
    this%b=b
  end subroutine set_b

  double precision function get_b(this)
    class(structure), intent(in) :: this
    get_b=this%b
  end function get_b

  subroutine set_c(this,c)
    implicit none
    class(structure), intent(inout) :: this
    real(8), intent(in)             :: c
    this%c=c
  end subroutine set_c

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

  subroutine set_nfree(this)
    class(structure), intent(inout) :: this
    this%nfree=3*this%natom-3
  end subroutine set_nfree

  integer function get_nfree(this)
    class(structure), intent(inout) :: this
    get_nfree=this%nfree
  end function get_nfree

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
         'ERROR: The number of atoms in HICOLM.xyz does not match with that defined in the INPUT file!'
    write(6,*)'Hint: Change the HICOLM.xyz or check the input in the &SYSTEM section.'
    stop
  end subroutine sites

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

  subroutine random_coordinates(this)
    implicit none
    class(structure), intent(inout) :: this
    integer                         :: i,j,k,nx
    real(8)                         :: srand,vx,vy,vz,vr,dr
    dr=0.1d0
    nx=1
    do i=1,this%get_nmol()
       do j=1,this%ntmol(i)
          call random_number(srand)
          vx=(2.0d0*srand-1.0d0)*dr
          call random_number(srand)
          vy=(2.0d0*srand-1.0d0)*dr
          call random_number(srand)
          vz=(2.0d0*srand-1.0d0)*dr
          vr=sqrt(vx**2+vy**2+vz**2)
          vx=vx*dr/vr
          vy=vy*dr/vr
          vz=vz*dr/vr
          do k=1,this%nxmol(i)
             this%xa(nx)=this%xa(nx)+vx
             this%ya(nx)=this%ya(nx)+vy
             this%za(nx)=this%za(nx)+vz
             nx=nx+1
          end do
       end do
    end do
  end subroutine random_coordinates

end module structure_module
