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
  !********************************************************************************
  !********************************************************************************

  use system_module
  use charmm_module
  use opls_module

  implicit none

  integer i,j,k,l

  private
  public :: forcefield

  type, extends(system) :: forcefield
     type(charmm)              :: charmm
     type(opls)                :: opls
     integer, private          :: nspcs
     integer, private          :: nvdw
     character(4), private     :: coulop
     character(6), private     :: ffmodel
     real(8), private          :: fscsalpha
     real(8), allocatable      :: parbnd(:,:,:)
     real(8), allocatable      :: parbend(:,:,:)
     real(8), allocatable      :: partors(:,:,:)
     real(8), allocatable      :: paritors(:,:,:)
     real(8), allocatable      :: parvdw(:,:)
     character(6), allocatable :: spcs(:)
     character(6), allocatable :: spcvdw(:,:)
     character(6), allocatable :: tbonds(:,:)
     character(6), allocatable :: tbends(:,:)
     character(7), allocatable :: ttors(:,:)
     character(7), allocatable :: titors(:,:)
     character(6), allocatable :: tvdw(:)
   contains
     procedure, private :: forcefield_init
     procedure          :: set_spcs
     procedure          :: set_nspcs
     procedure          :: get_nspcs
     procedure          :: set_parbnd
     procedure          :: set_parbend
     procedure          :: set_partors
     procedure          :: set_paritors
     procedure          :: set_parvdw
     procedure          :: check_vdw
     procedure          :: check_parbnd
     procedure          :: check_parbend
     procedure          :: check_partors
     procedure          :: check_paritors
     procedure          :: set_extra_bonds
     procedure          :: set_extra_angles
     procedure          :: set_extra_dihedrals
     procedure          :: set_extra_idihedrals
     procedure          :: set_extra_vdw
     procedure          :: set_extra_coulop
     procedure          :: set_nvdw
     procedure          :: get_nvdw
     procedure          :: set_fscsalpha
     procedure          :: get_fscsalpha
     procedure          :: set_coulop
     procedure          :: get_coulop
     procedure          :: set_ffmodel
     procedure          :: get_ffmodel
     procedure          :: set_topology
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
    integer                          :: i,j,nx
    call this%set_nspcs()
    nx=0
    do i=1,this%nspcs
       do j=i,this%nspcs
          nx=nx+1
       end do
    end do
    call this%set_nvdw(nx)
    allocate(this%tbonds(this%get_nmol(),this%get_bondmax()))
    allocate(this%tbends(this%get_nmol(),this%get_bendmax()))
    allocate(this%ttors(this%get_nmol(),this%get_torsmax()))
    allocate(this%titors(this%get_nmol(),this%get_itorsmax()))
    allocate(this%tvdw(this%get_nvdw()))
  end subroutine forcefield_init

  subroutine set_topology(this)
    implicit none
    class(forcefield), intent(inout) :: this
    integer                          :: i3,nx
    integer, allocatable             :: freeze(:)
    character(16)                    :: key
    character(10)                    :: cvar
    logical                          :: check,check2
    allocate(freeze(this%get_natom()))
    call this%forcefield_init()
    check=.true.
    check2=.false.
    do while(check)
       read(5,*,end=1)key
       if(key.eq.'&FORCE_FIELD'.or.key.eq.'&force_field')then
          check2=.true.
          do while(check2)
             read(5,*)key
             if(key.eq.'ffmodel')then
                backspace(5)
                read(5,*)key,cvar
                call this%set_ffmodel(cvar)
                check2=.false.
             elseif(key.eq.'&END_FORCE_FIELD'.or.key.eq.'&end_force_field')then
                check2=.false.
                check=.false.
             end if
          end do
          check=.false.
       end if
    end do
    call this%set_parbnd()
    call this%set_parbend()
    call this%set_partors()
    call this%set_paritors()
    call this%set_parvdw()
    rewind(5)
    check=.true.
    do while(check)
       read(5,*,end=1)key
       if(key.eq.'&FORCE_FIELD'.or.key.eq.'&force_field')check=.false.
    end do
    nx=0
    check=.true.
    do while(check)
       read(5,*)key
       if(key.eq.'molecule')then
          backspace(5)
          read(5,*)key,cvar
          i3=0
          do k=1,this%get_nmol()
             if(cvar.eq.this%namemol(k))i3=k
          end do
          if(i3.eq.0)goto 2
          check2=.true.
          do while(check2)
             read(5,*)key
             if(key.eq.'bonds')then
                backspace(5)
                call this%set_extra_bonds(i3)
             elseif(key.eq.'angles')then
                backspace(5)
                call this%set_extra_angles(i3)
             elseif(key.eq.'dihedrals')then
                backspace(5)
                call this%set_extra_dihedrals(i3)
             elseif(key.eq.'idihedrals')then
                backspace(5)
                call this%set_extra_idihedrals(i3)
             elseif(key.eq.'END_MOLECULE'.or.key.eq.'end_molecule')then
                check2=.false.
             end if
          end do
       elseif(key.eq.'vdw')then
          backspace(5)
          read(5,*)key,i3
          call this%set_extra_vdw(i3)
       elseif(key.eq.'electrostatic')then
          backspace(5)
          read(5,*)key,cvar
          call this%set_extra_coulop(cvar)
       elseif(key.eq.'freeze')then
          backspace(5)
          read(5,*)key,i3
          read(5,*)(freeze(i),i=1,i3)
          call this%set_extra_freeze(i3,freeze)
       elseif(key.eq.'&END_FORCE_FIELD'.or.key.eq.'&end_force_field')then
          check=.false.
       end if
    end do
1   rewind(5)
    return
2   write(6,*)&
         'ERROR: There is a molecule that does not belong to the physical system!'
    write(6,*)'Hint: Check the input in the &FORCE_FIELD section.'
    stop
  end subroutine set_topology

  subroutine set_nspcs(this)
    implicit none
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
    implicit none
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

  subroutine set_nvdw(this,nvdw)
    implicit none
    class(forcefield), intent(inout) :: this
    integer, intent(in)              :: nvdw
    this%nvdw=nvdw
  end subroutine set_nvdw

  integer function get_nvdw(this)
    implicit none
    class(forcefield), intent(in) :: this
    get_nvdw=this%nvdw
  end function get_nvdw

  subroutine set_parvdw(this)
    implicit none
    class(forcefield), intent(inout) :: this
    integer                          :: i,j,nx
    real(8)                          :: e1,e2,s1,s2
    call this%set_spcs()
    allocate(this%parvdw(this%get_nvdw(),2),this%spcvdw(this%get_nvdw(),2))
    if(this%get_ffmodel().eq.'charmm')then
       do i=1,this%get_nvdw()
          this%tvdw(i)='charmm'
       end do
       open(12,file='/tmp/hicolm3/charmm/charmm_vdw.prm',status='old')
       nx=1
       do i=1,this%get_nspcs()
          call this%charmm%set_charmmvdw(this%spcs(i))
          e1=this%charmm%prms_vdw(1)
          s1=this%charmm%prms_vdw(2)
          rewind(12)
          do j=i,this%get_nspcs()
             call this%charmm%set_charmmvdw(this%spcs(j))
             e2=this%charmm%prms_vdw(1)
             s2=this%charmm%prms_vdw(2)
             this%parvdw(nx,1)=sqrt(e1*e2)
             this%parvdw(nx,2)=s1+s2
             this%spcvdw(nx,1)=this%spcs(i)
             this%spcvdw(nx,2)=this%spcs(j)
             rewind(12)
             nx=nx+1
          end do
       end do
    elseif(this%get_ffmodel().eq.'opls')then
       do i=1,this%get_nvdw()
          this%tvdw(i)='opls'
       end do
       open(12,file='/tmp/hicolm3/opls/opls_vdw.prm',status='old')
       nx=1
       do i=1,this%get_nspcs()
          call this%opls%set_oplsvdw(this%spcs(i))
          e1=this%opls%prms_vdw(1)
          s1=this%opls%prms_vdw(2)
          rewind(12)
          do j=i,this%get_nspcs()
             call this%opls%set_oplsvdw(this%spcs(j))
             e2=this%opls%prms_vdw(1)
             s2=this%opls%prms_vdw(2)
             this%parvdw(nx,1)=sqrt(e1*e2)
             this%parvdw(nx,2)=0.5d0*(s1+s2)
             this%spcvdw(nx,1)=this%spcs(i)
             this%spcvdw(nx,2)=this%spcs(j)
             rewind(12)
             nx=nx+1
          end do
       end do
    end if
    call this%set_nvdw(nx-1)
    close(12)
  end subroutine set_parvdw

  subroutine check_vdw(this)
    implicit none
    class(forcefield), intent(inout) :: this
    integer                          :: i,nx
    real(8)                          :: ex,sx
    character(6)                     :: a1,a2
    nx=0
    do i=1,this%get_nvdw()
       ex=this%parvdw(i,1)
       sx=this%parvdw(i,2)
       a1=this%spcvdw(i,1)
       a2=this%spcvdw(i,2)
       if(ex.ge.1.d-4.and.sx.ge.5.d-4)then
          this%parvdw(nx+1,1)=ex
          this%parvdw(nx+1,2)=sx
          this%spcvdw(nx+1,1)=a1
          this%spcvdw(nx+1,2)=a2
          nx=nx+1
       end if
    end do
    call this%set_nvdw(nx)
  end subroutine check_vdw

  subroutine check_parbnd(this)
    implicit none
    class(forcefield), intent(inout) :: this
    integer                          :: i,nx
    real(8)                          :: p1,p2
    integer                          :: i1,i2
    logical                          :: check
    check=.true.
    do i=1,this%get_nmol()
       nx=0
       do j=1,this%bondscnt(i)
          p1=this%parbnd(i,j,1)
          p2=this%parbnd(i,j,2)
          i1=this%molbond(i,j,1)
          i2=this%molbond(i,j,2)
          if(p1.ge.1.d-4)then
             this%parbnd(i,nx+1,1)=p1
             this%parbnd(i,nx+1,2)=p2
             this%molbond(i,nx+1,1)=i1
             this%molbond(i,nx+1,2)=i2
             nx=nx+1
          end if
       end do
       if(nx.lt.this%bondscnt(i))check=.false.
       this%bondscnt(i)=nx
    end do
    if(check.eqv..false.)goto 1
    return
1   write(6,*)'Warning: There is unexpected unbound atoms in some molecules.'
    write(6,*)
  end subroutine check_parbnd

  subroutine check_parbend(this)
    implicit none
    class(forcefield), intent(inout) :: this
    integer                          :: i,nx
    real(8)                          :: p1,p2
    integer                          :: i1,i2,i3
    do i=1,this%get_nmol()
       nx=0
       do j=1,this%bendscnt(i)
          p1=this%parbend(i,j,1)
          p2=this%parbend(i,j,2)
          i1=this%molbend(i,j,1)
          i2=this%molbend(i,j,2)
          i3=this%molbend(i,j,3)
          if(p1.ge.1.d-4)then
             this%parbend(i,nx+1,1)=p1
             this%parbend(i,nx+1,2)=p2
             this%molbend(i,nx+1,1)=i1
             this%molbend(i,nx+1,2)=i2
             this%molbend(i,nx+1,3)=i3
             nx=nx+1
          end if
       end do
       this%bendscnt(i)=nx
    end do
  end subroutine check_parbend

  subroutine check_partors(this)
    implicit none
    class(forcefield), intent(inout) :: this
    integer                          :: i,nx
    real(8)                          :: p1,p2,p3
    integer                          :: i1,i2,i3,i4
    do i=1,this%get_nmol()
       nx=0
       do j=1,this%torscnt(i)
          p1=this%partors(i,j,1)
          p2=this%partors(i,j,2)
          p3=this%partors(i,j,3)
          i1=this%moltors(i,j,1)
          i2=this%moltors(i,j,2)
          i3=this%moltors(i,j,3)
          i4=this%moltors(i,j,4)
          if(this%get_ffmodel().eq.'charmm')then
             if(p1.ge.1.d-4)then
                this%partors(i,nx+1,1)=p1
                this%partors(i,nx+1,2)=p2
                this%partors(i,nx+1,3)=p3
                this%moltors(i,nx+1,1)=i1
                this%moltors(i,nx+1,2)=i2
                this%moltors(i,nx+1,3)=i3
                this%moltors(i,nx+1,4)=i4
                nx=nx+1
             end if
          elseif(this%get_ffmodel().eq.'opls')then
             if(p1.ge.1.d-4.or.p2.ge.1.d-4.or.p3.ge.1.d-4)then
                this%partors(i,nx+1,1)=p1
                this%partors(i,nx+1,2)=p2
                this%partors(i,nx+1,3)=p3
                this%moltors(i,nx+1,1)=i1
                this%moltors(i,nx+1,2)=i2
                this%moltors(i,nx+1,3)=i3
                this%moltors(i,nx+1,4)=i4
                nx=nx+1
             end if
          end if
       end do
       this%torscnt(i)=nx
    end do
  end subroutine check_partors

  subroutine check_paritors(this)
    implicit none
    class(forcefield), intent(inout) :: this
    integer                          :: i,nx
    real(8)                          :: p1,p2,p3
    integer                          :: i1,i2,i3,i4
    do i=1,this%get_nmol()
       nx=0
       do j=1,this%itorscnt(i)
          p1=this%paritors(i,j,1)
          p2=this%paritors(i,j,2)
          p3=this%paritors(i,j,3)
          i1=this%molitors(i,j,1)
          i2=this%molitors(i,j,2)
          i3=this%molitors(i,j,3)
          i4=this%molitors(i,j,4)
          if(p1.ge.1.d-4)then
             this%paritors(i,nx+1,1)=p1
             this%paritors(i,nx+1,2)=p2
             this%paritors(i,nx+1,3)=p3
             this%molitors(i,nx+1,1)=i1
             this%molitors(i,nx+1,2)=i2
             this%molitors(i,nx+1,3)=i3
             this%molitors(i,nx+1,4)=i4
             nx=nx+1
          end if
       end do
       this%itorscnt(i)=nx
    end do
  end subroutine check_paritors

  subroutine set_parbnd(this)
    implicit none
    class(forcefield), intent(inout) :: this
    integer                          :: i,j,m,i1,i2
    allocate(this%parbnd(this%get_nmol(),this%get_bondmax(),2))
    if(this%get_ffmodel().eq.'charmm')then
       do i=1,this%get_nmol()
          do j=1,this%get_bondmax()
             this%tbonds(i,j)='charmm'
          end do
       end do
       open(12,file='/tmp/hicolm3/charmm/charmm_bonds.prm',status='old')
       do i=1,this%get_nmol()
          do j=1,this%bondscnt(i)
             i1=this%molbond(i,j,1)
             i2=this%molbond(i,j,2)
             call this%charmm%set_charmmbonds(this%tpmol(i,i1),this%tpmol(i,i2))
             do m=1,2
                this%parbnd(i,j,m)=this%charmm%prms_bonds(m)
             end do
             rewind(12)
          end do
       end do
    elseif(this%get_ffmodel().eq.'opls')then
       do i=1,this%get_nmol()
          do j=1,this%get_bondmax()
             this%tbonds(i,j)='opls'
          end do
       end do
       open(12,file='/tmp/hicolm3/opls/opls_bonds.prm',status='old')
       do i=1,this%get_nmol()
          do j=1,this%bondscnt(i)
             i1=this%molbond(i,j,1)
             i2=this%molbond(i,j,2)
             call this%opls%set_oplsbonds(this%tpmol(i,i1),this%tpmol(i,i2))
             do m=1,2
                this%parbnd(i,j,m)=this%opls%prms_bonds(m)
             end do
             rewind(12)
          end do
       end do
    end if
    close(12)
  end subroutine set_parbnd

  subroutine set_parbend(this)
    implicit none
    class(forcefield), intent(inout) :: this
    integer                          :: i,j,m,i1,i2,i3
    allocate(this%parbend(this%get_nmol(),this%get_bendmax(),2))
    if(this%get_ffmodel().eq.'charmm')then
       do i=1,this%get_nmol()
          do j=1,this%get_bendmax()
             this%tbends(i,j)='charmm'
          end do
       end do
       open(12,file='/tmp/hicolm3/charmm/charmm_angles.prm',status='old')
       do i=1,this%get_nmol()
          do j=1,this%bendscnt(i)
             i1=this%molbend(i,j,1)
             i2=this%molbend(i,j,2)
             i3=this%molbend(i,j,3)
             call this%charmm%set_charmmangles&
                  (this%tpmol(i,i1),this%tpmol(i,i2),this%tpmol(i,i3))
             do m=1,2
                this%parbend(i,j,m)=this%charmm%prms_angles(m)
             end do
             rewind(12)
          end do
       end do
    elseif(this%get_ffmodel().eq.'opls')then
       do i=1,this%get_nmol()
          do j=1,this%get_bendmax()
             this%tbends(i,j)='opls'
          end do
       end do
       open(12,file='/tmp/hicolm3/opls/opls_angles.prm',status='old')
       do i=1,this%get_nmol()
          do j=1,this%bendscnt(i)
             i1=this%molbend(i,j,1)
             i2=this%molbend(i,j,2)
             i3=this%molbend(i,j,3)
             call this%opls%set_oplsangles&
                  (this%tpmol(i,i1),this%tpmol(i,i2),this%tpmol(i,i3))
             do m=1,2
                this%parbend(i,j,m)=this%opls%prms_angles(m)
             end do
             rewind(12)
          end do
       end do
    end if
    close(12)
  end subroutine set_parbend

  subroutine set_partors(this)
    implicit none
    class(forcefield), intent(inout) :: this
    integer                          :: i1,i2,i3,i4,i,j,m
    allocate(this%partors(this%get_nmol(),this%get_torsmax(),3))
    if(this%get_ffmodel().eq.'charmm')then
       do i=1,this%get_nmol()
          do j=1,this%get_torsmax()
             this%ttors(i,j)='charmm'
          end do
       end do
       open(12,file='/tmp/hicolm3/charmm/charmm_dihedrals.prm',status='old')
       do i=1,this%get_nmol()
          do j=1,this%torscnt(i)
             i1=this%moltors(i,j,1)
             i2=this%moltors(i,j,2)
             i3=this%moltors(i,j,3)
             i4=this%moltors(i,j,4)
             call this%charmm%set_charmmdihedrals(this%tpmol(i,i1),&
                  this%tpmol(i,i2),this%tpmol(i,i3),this%tpmol(i,i4))
             do m=1,3
                this%partors(i,j,m)=this%charmm%prms_tors(m)
             end do
             rewind(12)
          end do
       end do
    elseif(this%get_ffmodel().eq.'opls')then
       do i=1,this%get_nmol()
          do j=1,this%get_torsmax()
             this%ttors(i,j)='opls'
          end do
       end do
       open(12,file='/tmp/hicolm3/opls/opls_dihedrals.prm',status='old')
       do i=1,this%get_nmol()
          do j=1,this%torscnt(i)
             i1=this%moltors(i,j,1)
             i2=this%moltors(i,j,2)
             i3=this%moltors(i,j,3)
             i4=this%moltors(i,j,4)
             call this%opls%set_oplsdihedrals(this%tpmol(i,i1),&
                  this%tpmol(i,i2),this%tpmol(i,i3),this%tpmol(i,i4))
             do m=1,3
                this%partors(i,j,m)=this%opls%prms_tors(m)
             end do
             rewind(12)
          end do
       end do
    end if
    close(12)
  end subroutine set_partors

  subroutine set_paritors(this)
    implicit none
    class(forcefield), intent(inout) :: this
    integer                          :: i,j,m,i1,i2,i3,i4
    allocate(this%paritors(this%get_nmol(),this%get_itorsmax(),3))
    if(this%get_ffmodel().eq.'charmm')then
       do i=1,this%get_nmol()
          do j=1,this%get_itorsmax()
             this%titors(i,j)='charmm'
          end do
       end do
       open(12,file='/tmp/hicolm3/charmm/charmm_idihedrals.prm',status='old')
       do i=1,this%get_nmol()
          do j=1,this%itorscnt(i)
             i1=this%molitors(i,j,1)
             i2=this%molitors(i,j,2)
             i3=this%molitors(i,j,3)
             i4=this%molitors(i,j,4)
             call this%charmm%set_charmmidihedrals&
                  (this%tpmol(i,i1),this%tpmol(i,i2),&
                  this%tpmol(i,i3),this%tpmol(i,i4))
             do m=1,3
                this%paritors(i,j,m)=this%charmm%prms_itors(m)
             end do
             rewind(12)
          end do
       end do
    end if
    close(12)
  end subroutine set_paritors

  subroutine set_extra_vdw(this,nvdw)
    implicit none
    class(forcefield), intent(inout) :: this
    integer, intent(in)              :: nvdw
    integer                          :: i,j,k,nx
    real(8)                          :: p1,p2
    character(6)                     :: spcs1,spcs2
    character(6)                     :: tvdw
    nx=1
    nx=this%get_nvdw()
    do i=1,nvdw
       read(5,*)spcs1,spcs2,tvdw,p1,p2
       do j=1,this%get_nvdw()
          if(spcs1.eq.this%spcvdw(j,1).and.spcs2.eq.this%spcvdw(j,2).or.&
               spcs1.eq.this%spcvdw(j,2).and.spcs2.eq.this%spcvdw(j,1))then
             this%parvdw(j,1)=p1
             this%parvdw(j,2)=p2
             this%tvdw(j)=tvdw
             goto 1
          end if
       end do
       do j=1,this%get_nspcs()
          do k=1,this%get_nspcs()
             if(spcs1.eq.this%spcs(j).and.spcs2.eq.this%spcs(k))then
                this%parvdw(nx+1,1)=p1
                this%parvdw(nx+1,2)=p2
                this%tvdw(nx+1)=tvdw
                this%spcvdw(nx+1,1)=spcs1
                this%spcvdw(nx+1,2)=spcs2
                nx=nx+1
                goto 1
             end if
          end do
       end do
       goto 2
1      continue
    end do
    call this%set_nvdw(nx)
    return
2   write(6,*)'ERROR: The type does not match with that defined in the SYSTEM section!'
    write(6,*)'Hint: Check the input in the &FORCE_FIELD section.'
    stop
  end subroutine set_extra_vdw

  subroutine set_extra_bonds(this,i3)
    implicit none
    class(forcefield), intent(inout) :: this
    integer, intent(in)              :: i3
    integer                          :: i,j,i1,i2
    character(12)                    :: key
    read(5,*)key,i1
    do i=1,i1
       read(5,*)i2
       backspace(5)
       read(5,*)i2,this%molbond(i3,i2,1),this%molbond(i3,i2,2),this%tbonds(i3,i2)
       backspace(5)
       select case(this%tbonds(i3,i2))
       case('charmm')
          read(5,*)i2,this%molbond(i3,i2,1),this%molbond(i3,i2,2),&
               this%tbonds(i3,i2),(this%parbnd(i3,i2,j),j=1,2)
       case('harm')
          read(5,*)i2,this%molbond(i3,i2,1),this%molbond(i3,i2,2),&
               this%tbonds(i3,i2),(this%parbnd(i3,i2,j),j=1,2)
       case('opls')
          read(5,*)i2,this%molbond(i3,i2,1),this%molbond(i3,i2,2),&
               this%tbonds(i3,i2),(this%parbnd(i3,i2,j),j=1,2)
       end select
    end do
  end subroutine set_extra_bonds

  subroutine set_extra_angles(this,i3)
    implicit none
    class(forcefield), intent(inout) :: this
    integer, intent(in)              :: i3
    integer                          :: i1,i2
    character(12)                    :: key
    read(5,*)key,i1
    do k=1,i1
       read(5,*)i2
       backspace(5)
       read(5,*)i2,this%molbend(i3,i2,1),this%molbend(i3,i2,2),&
            this%molbend(i3,i2,3),this%tbends(i3,i2)
       backspace(5)
       select case(this%tbends(i3,i2))
       case('charmm')
          read(5,*)i2,this%molbend(i3,i2,1),this%molbend(i3,i2,2),&
               this%molbend(i3,i2,3),this%tbends(i3,i2),&
               (this%parbend(i3,i2,l),l=1,2)
       case('harm')
          read(5,*)i2,this%molbend(i3,i2,1),this%molbend(i3,i2,2),&
               this%molbend(i3,i2,3),this%tbends(i3,i2),&
               (this%parbend(i3,i2,l),l=1,2)
       case('opls')
          read(5,*)i2,this%molbend(i3,i2,1),this%molbend(i3,i2,2),&
               this%molbend(i3,i2,3),this%tbends(i3,i2),&
               (this%parbend(i3,i2,l),l=1,2)
       end select
    end do
  end subroutine set_extra_angles

  subroutine set_extra_dihedrals(this,i3)
    implicit none
    class(forcefield), intent(inout) :: this
    integer, intent(in)              :: i3
    integer                          :: i1,i2,k,l
    character(12)                    :: key
    read(5,*)key,i1
    do k=1,i1
       read(5,*)i2
       backspace(5)
       read(5,*)i2,this%moltors(i3,i2,1),this%moltors(i3,i2,2),&
            this%moltors(i3,i2,3),this%moltors(i3,i2,4),this%ttors(i3,i2)
       backspace(5)
       select case(this%ttors(i3,i2))
       case('charmm')
          read(5,*)i2,this%moltors(i3,i2,1),this%moltors(i3,i2,2),&
               this%moltors(i3,i2,3),this%moltors(i3,i2,4),this%ttors(i3,i2),&
               (this%partors(i3,i2,l),l=1,3)
       case('charmm2')
          read(5,*)i2,this%moltors(i3,i2,1),this%moltors(i3,i2,2),&
               this%moltors(i3,i2,3),this%moltors(i3,i2,4),&
               this%ttors(i3,i2),(this%partors(i3,i2,l),l=1,2)
       case('harm')
          read(5,*)i2,this%moltors(i3,i2,1),this%moltors(i3,i2,2),&
               this%moltors(i3,i2,3),this%moltors(i3,i2,4),this%ttors(i3,i2),&
               (this%partors(i3,i2,l),l=1,2)
       case('opls')
          read(5,*)i2,this%moltors(i3,i2,1),this%moltors(i3,i2,2),&
               this%moltors(i3,i2,3),this%moltors(i3,i2,4),this%ttors(i3,i2),&
               (this%partors(i3,i2,l),l=1,3)
       end select
    end do
  end subroutine set_extra_dihedrals

  subroutine set_extra_idihedrals(this,i3)
    implicit none
    class(forcefield), intent(inout) :: this
    integer, intent(in)              :: i3
    integer                          :: i1,i2,k,l
    character(12)                    :: key
    read(5,*)key,i1
    do k=1,i1
       read(5,*)i2
       backspace(5)
       read(5,*)i2,this%molitors(i3,i2,1),this%molitors(i3,i2,2),&
            this%molitors(i3,i2,3),this%molitors(i3,i2,4),this%titors(i3,i2)
       backspace(5)
       select case(this%titors(i3,i2))
       case('charmm')
          read(5,*)i2,this%molitors(i3,i2,1),this%molitors(i3,i2,2),&
               this%molitors(i3,i2,3),this%molitors(i3,i2,4),this%titors(i3,i2),&
               this%paritors(i3,i2,1),this%paritors(i3,i2,3)
          this%paritors(i3,i2,2)=0.d0
       case('charmm2')
          read(5,*)i2,this%molitors(i3,i2,1),this%molitors(i3,i2,2),&
               this%molitors(i3,i2,3),this%molitors(i3,i2,4),&
               this%titors(i3,i2),(this%paritors(i3,i2,l),l=1,3)
       case('harm')
          read(5,*)i2,this%molitors(i3,i2,1),this%molitors(i3,i2,2),&
               this%molitors(i3,i2,3),this%molitors(i3,i2,4),this%titors(i3,i2),&
               (this%paritors(i3,i2,l),l=1,2)
       end select
    end do
  end subroutine set_extra_idihedrals

  subroutine set_coulop(this,coulop)
    implicit none
    class(forcefield), intent(inout) :: this
    character(4), intent(in)         :: coulop
    this%coulop=coulop
  end subroutine set_coulop

  subroutine set_extra_coulop(this,coulop)
    implicit none
    class(forcefield), intent(inout) :: this
    character(4), intent(in)         :: coulop
    real(8)                          :: fscsalpha
    character(4)                     :: cvar
    character(13)                    :: key
    if(coulop.eq.'fscs')then
       backspace(5)
       read(5,*)key,cvar,fscsalpha
       this%fscsalpha=fscsalpha
       goto 1
    elseif(coulop.eq.'coul')then
       goto 1
    end if
    goto 2
1   this%coulop=coulop
    return
2   write(6,*)'ERROR: The electrostatic does not an option!'
    write(6,*)'Hint: Check the input in the &FORCE_FIELD section.'
    stop
  end subroutine set_extra_coulop

  character(4) function get_coulop(this)
    implicit none
    class(forcefield), intent(in) :: this
    get_coulop=this%coulop
  end function get_coulop

  subroutine set_fscsalpha(this,fscsalpha)
    implicit none
    class(forcefield), intent(inout) :: this
    real(8), intent(in)              :: fscsalpha
    this%fscsalpha=fscsalpha
  end subroutine set_fscsalpha

  double precision function get_fscsalpha(this)
    implicit none
    class(forcefield), intent(in) :: this
    get_fscsalpha=this%fscsalpha
  end function get_fscsalpha

  subroutine set_ffmodel(this,ffmodel)
    implicit none
    class(forcefield), intent(inout) :: this
    character(6), intent(in)         :: ffmodel
    this%ffmodel=ffmodel
  end subroutine set_ffmodel

  character(6) function get_ffmodel(this)
    implicit none
    class(forcefield), intent(in) :: this
    get_ffmodel=this%ffmodel
  end function get_ffmodel

end module forcefield_module
