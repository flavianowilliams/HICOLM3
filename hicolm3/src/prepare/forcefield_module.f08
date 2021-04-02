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
     integer, private          :: itorsmax
     integer, allocatable      :: itorscnt(:)
     integer, allocatable      :: molitors(:,:,:)
     character(4), private     :: coulop
     real(8), allocatable      :: parbnd(:,:,:)
     real(8), allocatable      :: parbend(:,:,:)
     real(8), allocatable      :: partors(:,:,:)
     real(8), allocatable      :: paritors(:,:,:)
     real(8), allocatable      :: parvdw(:,:)
     character(2), allocatable :: spcs(:)
     character(2), allocatable :: spcvdw(:,:)
     character(5), allocatable :: tbonds(:,:)
     character(5), allocatable :: tbends(:,:)
     character(5), allocatable :: ttors(:,:)
     character(5), allocatable :: titors(:,:)
     character(5), allocatable :: tvdw(:)
   contains
     procedure, private :: forcefield_init
     procedure          :: set_spcs
     procedure          :: set_nspcs
     procedure          :: get_nspcs
     procedure          :: set_parbnd
     procedure          :: set_parbend
     procedure          :: set_partors
     procedure          :: set_paritors
     procedure          :: set_extra_parbnd
     procedure          :: set_extra_parbend
     procedure          :: set_extra_partors
     procedure          :: set_extra_parvdw
     procedure          :: set_parvdw
     procedure          :: set_nvdw
     procedure          :: get_nvdw
     procedure          :: set_coulop
     procedure          :: get_coulop
     procedure          :: set_itorsmax
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
    allocate(this%tbonds(this%get_nmol(),this%get_bondmax()))
    allocate(this%tbends(this%get_nmol(),this%get_bendmax()))
    allocate(this%ttors(this%get_nmol(),this%get_torsmax()))
    allocate(this%tvdw(this%nspcs))
    do i=1,this%get_nmol()
       do j=1,this%get_bondmax()
          this%tbonds(i,j)='amber'
       end do
       do j=1,this%get_bendmax()
          this%tbends(i,j)='amber'
       end do
       do j=1,this%get_torsmax()
          this%ttors(i,j)='amber'
       end do
    end do
    do i=1,this%nspcs
       this%tvdw(i)='amber'
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
    implicit none
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
    integer                          :: i,j,k,nx
    real(8)                          :: e1,e2,s1,s2
    call this%set_spcs()
    allocate(this%parvdw(this%nspcs,2),this%spcvdw(this%nspcs,2))
    nx=1
    do i=1,this%nspcs
       do k=1,this%amber%get_natp()
          if(this%spcs(i).eq.this%amber%atp(k))then
             e1=this%amber%prms_vdw(k,2)
             s1=this%amber%prms_vdw(k,1)
          end if
       end do
       if(e1.ge.1.d-4.and.s1.ge.1.d-1)then
          do j=i,this%nspcs
             do k=1,this%amber%get_natp()
                if(this%spcs(j).eq.this%amber%atp(k))then
                   e2=this%amber%prms_vdw(k,2)
                   s2=this%amber%prms_vdw(k,1)
                end if
             end do
             if(e2.ge.1.d-4.and.s2.ge.1.d-1)then
                this%parvdw(nx,1)=sqrt(e1*e2)
                this%parvdw(nx,2)=s1+s2
                this%spcvdw(nx,1)=this%spcs(i)
                this%spcvdw(nx,2)=this%spcs(j)
                nx=nx+1
             end if
          end do
       end if
    end do
    this%nvdw=nx-1
  end subroutine set_parvdw

  subroutine set_nvdw(this,nvdw)
    implicit none
    class(forcefield), intent(inout) :: this
    integer, intent(in)              :: nvdw
    this%nvdw=nvdw
  end subroutine set_nvdw

  integer function get_nvdw(this)
    class(forcefield), intent(in) :: this
    get_nvdw=this%nvdw
  end function get_nvdw

  subroutine set_parbnd(this)
    class(forcefield), intent(inout) :: this
    integer                          :: i,j,k,l,m,i1,i2
    allocate(this%parbnd(this%get_nmol(),this%get_bondmax(),2))
    call this%forcefield_init()
    do i=1,this%get_nmol()
       do j=1,this%bondscnt(i)
          i1=this%molbond(i,j,1)
          i2=this%molbond(i,j,2)
          do k=1,this%amber%get_natp()
             do l=1,this%amber%get_natp()
                if(this%amber%atp(k).eq.this%tpmol(i,i1).and.&
                     this%amber%atp(l).eq.this%tpmol(i,i2))then
                   do m=1,2
                      this%parbnd(i,j,m)=this%amber%prms_bonds(k,l,m)
                   end do
                elseif(this%amber%atp(k).eq.this%tpmol(i,i2).and.&
                     this%amber%atp(l).eq.this%tpmol(i,i1))then
                   do m=1,2
                      this%parbnd(i,j,m)=this%amber%prms_bonds(k,l,m)
                   end do
                end if
             end do
          end do
       end do
    end do
  end subroutine set_parbnd

  subroutine set_parbend(this)
    class(forcefield), intent(inout) :: this
    integer                          :: i,j,k,l,m,n,i1,i2,i3
    allocate(this%parbend(this%get_nmol(),this%get_bendmax(),2))
    do i=1,this%get_nmol()
       do j=1,this%bendscnt(i)
          i1=this%molbend(i,j,1)
          i2=this%molbend(i,j,2)
          i3=this%molbend(i,j,3)
          do k=1,this%amber%get_natp()
             do l=1,this%amber%get_natp()
                do m=1,this%amber%get_natp()
                   if(this%amber%atp(k).eq.this%tpmol(i,i1).and.&
                        this%amber%atp(l).eq.this%tpmol(i,i2).and.&
                        this%amber%atp(m).eq.this%tpmol(i,i3))then
                      do n=1,2
                         this%parbend(i,j,n)=this%amber%prms_angles(k,l,m,n)
                      end do
                   elseif(this%amber%atp(k).eq.this%tpmol(i,i3).and.&
                        this%amber%atp(l).eq.this%tpmol(i,i2).and.&
                        this%amber%atp(m).eq.this%tpmol(i,i1))then
                      do n=1,2
                         this%parbend(i,j,n)=this%amber%prms_angles(k,l,m,n)
                      end do
                   end if
                end do
             end do
          end do
       end do
    end do
  end subroutine set_parbend

  subroutine set_partors(this)
    class(forcefield), intent(inout) :: this
    integer                          :: i1,i2,i3,i4
    allocate(this%partors(this%get_nmol(),this%get_bendmax(),4))
    do i=1,this%get_nmol()
       do j=1,this%torscnt(i)
          i1=this%moltors(i,j,1)
          i2=this%moltors(i,j,2)
          i3=this%moltors(i,j,3)
          i4=this%moltors(i,j,4)
          call this%amber%set_amber&
               (this%tpmol(i,i1),this%tpmol(i,i2),this%tpmol(i,i3),this%tpmol(i,i4))
          do k=1,4
             this%partors(i,j,k)=this%amber%prms_tors(k)
          end do
       end do
    end do
  end subroutine set_partors

  subroutine set_itorsmax(this)
    class(forcefield), intent(inout) :: this
    integer                          :: i1,nxx
    character(12)                    :: key
    nxx=0
1   read(5,*,end=3)key
    if(key.ne.'&FORCE_FIELD')goto 1
    do j=1,this%get_nmol()
       do while (key.ne.'&END')
          read(5,*)key
          if(key.eq.'molecule')then
             read(5,*)
             read(5,*)
             read(5,*)
             do while (key.ne.'end_molecule')
                read(5,*)key
                if(key.eq.'end_molecule')goto 2
                if(key.eq.'idihedrals')then
                   backspace(5)
                   read(5,*)key,i1
                end if
             end do
          end if
2         continue
       end do
       nxx=max(nxx,i1)
    end do
    this%itorsmax=nxx
    call this%set_torsmax(this%get_torsmax()+this%itorsmax)
3   rewind(5)
  end subroutine set_itorsmax

  subroutine set_paritors(this)
    class(forcefield), intent(inout) :: this
    integer                          :: i2,i3
    character(12)                    :: key
    character(10)                    :: cvar
    i3=0
    call this%set_itorsmax()
    allocate(this%itorscnt(this%get_nmol()))
    allocate(this%molitors(this%get_nmol(),this%itorsmax,4))
    allocate(this%paritors(this%get_nmol(),this%itorsmax,4))
    allocate(this%titors(this%get_nmol(),this%itorsmax))
    if(this%itorsmax.eq.0)goto 3
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
                if(key.eq.'idihedrals')then
                   backspace(5)
                   read(5,*)key,this%itorscnt(i3)
                   do k=1,this%itorscnt(i3)
                      read(5,*)i2
                      backspace(5)
                      read(5,*)i2,this%molitors(i3,i2,1),this%molitors(i3,i2,2),&
                           this%molitors(i3,i2,3),this%molitors(i3,i2,4),this%titors(i3,i2)
                      backspace(5)
                      select case(this%titors(i3,i2))
                      case('amber')
                         read(5,*)i2,this%molitors(i3,i2,1),this%molitors(i3,i2,2),&
                              this%molitors(i3,i2,3),this%molitors(i3,i2,4),&
                              this%titors(i3,i2),(this%paritors(i3,i2,l),l=1,4)
                      case('harm')
                         read(5,*)i2,this%molitors(i3,i2,1),this%molitors(i3,i2,2),&
                              this%molitors(i3,i2,3),this%molitors(i3,i2,4),&
                              this%titors(i3,i2),(this%paritors(i3,i2,l),l=1,2)
                      end select
                   end do
                end if
             end do
          end if
2         continue
       end do
    end do
    rewind(5)
    return
3   do i=1,this%get_nmol()
       this%itorscnt(i)=0
    end do
    rewind(5)
  end subroutine set_paritors

  subroutine set_extra_parvdw(this)
    implicit none
    class(forcefield), intent(inout) :: this
    integer                          :: nvdw
    real(8)                          :: p1,p2
    character(2)                     :: spcs1,spcs2
    character(12)                    :: key
    character(5)                     :: tvdw
    nvdw=this%get_nvdw()
1   read(5,*,end=3)key
    if(key.ne.'&FORCE_FIELD')goto 1
    do while (key.ne.'&END')
       read(5,*)key
       if(key.eq.'vdw')then
          backspace(5)
          read(5,*)key,nvdw
          do i=1,nvdw
             read(5,*)spcs1,spcs2,tvdw,p1,p2
             do j=1,this%get_nvdw()
                if(spcs1.eq.this%spcvdw(j,1).and.spcs2.eq.this%spcvdw(j,2).or.&
                     spcs1.eq.this%spcvdw(j,2).and.spcs2.eq.this%spcvdw(j,1))then
                   this%parvdw(j,1)=p1
                   this%parvdw(j,2)=p2
                   this%tvdw(j)=tvdw
                   goto 2
                end if
             end do
             do j=1,this%get_nspcs()
                do k=1,this%get_nspcs()
                   if(spcs1.eq.this%spcs(j).and.spcs2.eq.this%spcs(k))then
                      this%parvdw(nvdw+1,1)=p1
                      this%parvdw(nvdw+1,2)=p2
                      this%tvdw(nvdw+1)=tvdw
                      this%spcvdw(nvdw+1,1)=spcs1
                      this%spcvdw(nvdw+1,2)=spcs2
                      nvdw=nvdw+1
                      goto 2
                   end if
                end do
             end do
             goto 4
2            continue
          end do
       end if
    end do
    call this%set_nvdw(nvdw)
3   rewind(5)
    return
4   write(6,*)'ERROR: The type does not match with that defined in the TOPOLOGY file!'
    write(6,*)'Hint: Check the input in the &FORCE_FIELD section.'
    stop
  end subroutine set_extra_parvdw

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
                              this%molbend(i3,i2,3),this%tbends(i3,i2),&
                              (this%parbend(i3,i2,l),l=1,2)
                      case('harm')
                         read(5,*)i2,this%molbend(i3,i2,1),this%molbend(i3,i2,2),&
                              this%molbend(i3,i2,3),this%tbends(i3,i2),&
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

  subroutine set_extra_partors(this)
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
                if(key.eq.'dihedrals')then
                   backspace(5)
                   read(5,*)key,i1
                   do k=1,i1
                      read(5,*)i2
                      backspace(5)
                      read(5,*)i2,this%moltors(i3,i2,1),this%moltors(i3,i2,2),&
                           this%moltors(i3,i2,3),this%moltors(i3,i2,4),this%ttors(i3,i2)
                      backspace(5)
                      select case(this%ttors(i3,i2))
                      case('amber')
                         read(5,*)i2,this%moltors(i3,i2,1),this%moltors(i3,i2,2),&
                              this%moltors(i3,i2,3),this%moltors(i3,i2,4),this%ttors(i3,i2),&
                              (this%partors(i3,i2,l),l=1,4)
                      case('harm')
                         read(5,*)i2,this%moltors(i3,i2,1),this%moltors(i3,i2,2),&
                              this%moltors(i3,i2,3),this%moltors(i3,i2,4),this%ttors(i3,i2),&
                              (this%partors(i3,i2,l),l=1,2)
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
  end subroutine set_extra_partors

  subroutine set_coulop(this,coulop)
    class(forcefield), intent(inout) :: this
    character(4)                     :: coulop
    character(13)                    :: key
    this%coulop=coulop
1   read(5,*,end=2)key
    if(key.ne.'&FORCE_FIELD')goto 1
    do while (key.ne.'&END')
       read(5,*)key
       if(key.eq.'electrostatic')then
          backspace(5)
          read(5,*)key,this%coulop
          goto 2
       end if
    end do
2   rewind(5)
  end subroutine set_coulop

  character(4) function get_coulop(this)
    class(forcefield), intent(inout) :: this
    get_coulop=this%coulop
  end function get_coulop

end module forcefield_module
