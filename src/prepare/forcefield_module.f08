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
module forcefield_module
  !*******************************************************************************************
  !*******************************************************************************************

  use zmatrix_module
  use amber_module

  implicit none

  integer i,j,k,l

  private
  public :: forcefield

  type, extends(zmatrix) :: forcefield
     type(amber)               :: amber
     integer, private          :: nspcs
     integer, private          :: nvdw
     real(8), allocatable      :: parbnd(:,:,:)
     real(8), allocatable      :: parbend(:,:,:)
     real(8), allocatable      :: parvdw(:,:,:)
     character(2), allocatable :: spcs(:)
     character(5), allocatable :: tbonds(:,:)
     character(5), allocatable :: tbends(:,:)
   contains
     procedure, private :: forcefield_init
     procedure          :: set_spcs
     procedure          :: set_nspcs
     procedure          :: get_nspcs
     procedure          :: set_parbnd
     procedure          :: set_parbend
     procedure          :: set_extra_parbnd
     procedure          :: set_extra_parbend
     procedure          :: set_parvdw
     procedure          :: get_nvdw
  end type forcefield

  interface forcefield
     module procedure constructor
  end interface forcefield

contains

  type(forcefield) function constructor()
    call constructor%forcefield_init()
  end function constructor

  subroutine forcefield_init(this)
    implicit none
    class(forcefield), intent(inout) :: this
    call this%set_nspcs()
    allocate(this%parbnd(this%get_nmol(),this%get_bondmax(),2))
    allocate(this%parbend(this%get_nmol(),this%get_bendmax(),2))
    allocate(this%parvdw(this%nspcs,this%nspcs,2))
    allocate(this%tbonds(this%get_nmol(),this%get_bondmax()))
    allocate(this%tbends(this%get_nmol(),this%get_bendmax()))
    do i=1,this%get_nmol()
       do j=1,this%get_bondmax()
          this%tbonds(i,j)='amber'
       end do
       do j=1,this%get_bendmax()
          this%tbends(i,j)='amber'
       end do
    end do
  end subroutine forcefield_init

  subroutine set_nspcs(this)
    class(forcefield), intent(inout) :: this
    integer                          :: nx
    logical                          :: check
    nx=0
    check=.FALSE.
    do i=1,this%get_nmol()
       do j=1,this%nxmol(i)
          do k=1,(i-1)
             do l=1,this%nxmol(k)
                if(this%tpmol(k,l).eq.this%tpmol(i,j))check=.TRUE.
             end do
          end do
          do l=1,(j-1)
             if(this%tpmol(i,l).eq.this%tpmol(i,j))check=.TRUE.
          end do
          if(check.neqv..TRUE.)nx=nx+1
          check=.FALSE.
       end do
    end do
    this%nspcs=nx
  end subroutine set_nspcs

  integer function get_nspcs(this)
    class(forcefield), intent(in) :: this
    get_nspcs=this%nspcs
  end function get_nspcs

  subroutine set_spcs(this)
    class(forcefield), intent(inout) :: this
    integer                          :: nx
    allocate(this%spcs(this%nspcs))
    this%spcs(1)=this%tpmol(1,1)
    nx=1
    do i=1,this%get_nmol()
       do j=1,this%nxmol(i)
          do k=1,nx
             if(this%tpmol(i,j).eq.this%spcs(k))goto 1
          end do
          this%spcs(nx+1)=this%tpmol(i,j)
          nx=nx+1
1         continue
       end do
    end do
  end subroutine set_spcs

  subroutine set_parvdw(this)
    class(forcefield), intent(inout) :: this
    integer                          :: nx
    real(8)                          :: e1,e2,s1,s2
    call this%set_spcs()
    nx=0
    do i=1,this%nspcs
       call this%amber%set_amber(this%spcs(i))
       e1=this%amber%prms_vdw(2)
       s1=this%amber%prms_vdw(1)
       do j=i,this%nspcs
          call this%amber%set_amber(this%spcs(j))
          e2=this%amber%prms_vdw(2)
          s2=this%amber%prms_vdw(1)
          this%parvdw(i,j,1)=sqrt(e1*e2)
          this%parvdw(i,j,2)=s1+s2
          if(this%parvdw(i,j,1).ge.1.d-8.and.this%parvdw(i,j,2).ge.1.d-2)nx=nx+1
       end do
    end do
    this%nvdw=nx
  end subroutine set_parvdw

  integer function get_nvdw(this)
    class(forcefield), intent(in) :: this
    get_nvdw=this%nvdw
  end function get_nvdw

  subroutine set_parbnd(this)
    class(forcefield), intent(inout) :: this
    integer                          :: i1,i2
    call this%forcefield_init()
    do i=1,this%get_nmol()
       do j=1,this%bondscnt(i)
          i1=this%molbond(i,j,1)
          i2=this%molbond(i,j,2)
          call this%amber%set_amber(this%tpmol(i,i1),this%tpmol(i,i2))
          do k=1,2
             this%parbnd(i,j,k)=this%amber%prms_bonds(k)
          end do
       end do
    end do
  end subroutine set_parbnd

  subroutine set_parbend(this)
    class(forcefield), intent(inout) :: this
    integer                          :: i1,i2,i3
    do i=1,this%get_nmol()
       do j=1,this%bendscnt(i)
          i1=this%molbend(i,j,1)
          i2=this%molbend(i,j,2)
          i3=this%molbend(i,j,3)
          call this%amber%set_amber(this%tpmol(i,i1),this%tpmol(i,i2),this%tpmol(i,i3))
          do k=1,2
             this%parbend(i,j,k)=this%amber%prms_bends(k)
          end do
       end do
    end do
  end subroutine set_parbend

  subroutine set_extra_parbnd(this)
    class(forcefield), intent(inout) :: this
    integer                          :: i1,i2,i3
    character(12)                    :: key
    character(10)                    :: cvar
    i3=0
1   read(5,*,end=3)key
    if(key.ne.'&FORCE_FIELD')goto 1
    do j=1,this%get_nmol()
       do while (key.ne.'&END')
          read(5,*)key
          if(key.eq.'molecule')then
             backspace(5)
             read(5,*)key,cvar
             do k=1,this%get_nmol()
                if(cvar.eq.this%namemol(k))then
                   i3=k
                end if
             end do
             read(5,*)
             read(5,*)
             read(5,*)
             do while (key.ne.'end_molecule')
                read(5,*)key
                if(key.eq.'end_molecule')goto 2
                if(key.eq.'bonds')then
                   backspace(5)
                   read(5,*)key,i1
                   do k=1,i1
                      read(5,*)i2
                      backspace(5)
                      read(5,*)i2,this%molbond(i3,i2,1),this%molbond(i3,i2,2),&
                           this%tbonds(i3,i2)
                      backspace(5)
                      select case(this%tbonds(i3,i2))
                      case('amber')
                         read(5,*)i2,this%molbond(i3,i2,1),this%molbond(i3,i2,2),&
                              this%tbonds(i3,i2),(this%parbnd(i3,i2,l),l=1,2)
                      case('harm')
                         read(5,*)i2,this%molbond(i3,i2,1),this%molbond(i3,i2,2),&
                              this%tbonds(i3,i2),(this%parbnd(i3,i2,l),l=1,2)
                      end select
                   end do
                end if
             end do
          end if
2         continue
       end do
    end do
3   rewind(5)
    return
  end subroutine set_extra_parbnd

  subroutine set_extra_parbend(this)
    class(forcefield), intent(inout) :: this
    integer                          :: i1,i2,i3
    character(12)                    :: key
    character(10)                    :: cvar
    i3=0
1   read(5,*,end=3)key
    if(key.ne.'&FORCE_FIELD')goto 1
    do j=1,this%get_nmol()
       do while (key.ne.'&END')
          read(5,*)key
          if(key.eq.'molecule')then
             backspace(5)
             read(5,*)key,cvar
             do k=1,this%get_nmol()
                if(cvar.eq.this%namemol(k))then
                   i3=k
                end if
             end do
             read(5,*)
             read(5,*)
             read(5,*)
             do while (key.ne.'end_molecule')
                read(5,*)key
                if(key.eq.'end_molecule')goto 2
                if(key.eq.'bends')then
                   backspace(5)
                   read(5,*)key,i1
                   do k=1,i1
                      read(5,*)i2
                      backspace(5)
                      read(5,*)i2,this%molbend(i3,i2,1),this%molbend(i3,i2,2),&
                           this%molbend(i3,i2,3),this%tbends(i3,i2)
                      backspace(5)
                      select case(this%tbends(i3,i2))
                      case('amber')
                         read(5,*)i2,this%molbend(i3,i2,1),this%molbend(i3,i2,2),&
                              this%molbend(i3,i2,3),this%tbonds(i3,i2),&
                              (this%parbend(i3,i2,l),l=1,2)
                      case('harm')
                         read(5,*)i2,this%molbend(i3,i2,1),this%molbend(i3,i2,2),&
                              this%molbend(i3,i2,2),this%tbonds(i3,i2),&
                              (this%parbend(i3,i2,l),l=1,2)
                      end select
                   end do
                end if
             end do
          end if
2         continue
       end do
    end do
3   rewind(5)
    return
  end subroutine set_extra_parbend

end module forcefield_module
