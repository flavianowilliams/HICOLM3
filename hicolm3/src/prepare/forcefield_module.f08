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
  use charmm_module

  implicit none

  integer i,j,k,l

  private
  public :: forcefield

  type, extends(zmatrix) :: forcefield
     type(charmm)              :: charmm
     integer, private          :: nspcs
     integer, private          :: nvdw
     integer, allocatable      :: itorscnt(:)
     integer, allocatable      :: molitors(:,:,:)
     character(4), private     :: coulop
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
     procedure          :: set_extra_bonds
     procedure          :: set_extra_angles
     procedure          :: set_extra_dihedrals
     procedure          :: set_extra_parvdw
     procedure          :: set_parvdw
     procedure          :: set_nvdw
     procedure          :: get_nvdw
     procedure          :: set_fscsalpha
     procedure          :: get_fscsalpha
     procedure          :: set_coulop
     procedure          :: get_coulop
     procedure          :: set_coulop2
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
    do i=1,this%get_nmol()
       do j=1,this%get_bondmax()
          this%tbonds(i,j)='charmm'
       end do
       do j=1,this%get_bendmax()
          this%tbends(i,j)='charmm'
       end do
       do j=1,this%get_torsmax()
          this%ttors(i,j)='charmm'
       end do
       do j=1,this%get_itorsmax()
          this%titors(i,j)='icharmm'
       end do
    end do
    do i=1,this%get_nvdw()
       this%tvdw(i)='charmm'
    end do
  end subroutine forcefield_init

  subroutine set_topology(this)
    implicit none
    class(forcefield), intent(inout) :: this
    integer                          :: i3,nx
    character(16)                    :: key
    character(10)                    :: cvar
    logical                          :: check
!    allocate(this%zatmol(this%get_nmol(),this%get_natom()))
!    allocate(this%qatmol(this%get_nmol(),this%get_natom()))
!    allocate(this%tpmol(this%get_nmol(),this%get_natom()))
    call this%forcefield_init()
    call this%set_parbnd()
    call this%set_parbend()
    call this%set_partors()
    call this%set_paritors()
    call this%set_parvdw()
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
          read(5,*)key
          if(key.eq.'bonds')then
             backspace(5)
             call this%set_extra_bonds(i3)
          end if
          if(key.eq.'angles')then
             backspace(5)
             call this%set_extra_angles(i3)
          end if
          if(key.eq.'dihedrals')then
             backspace(5)
             call this%set_extra_dihedrals(i3)
          end if
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
    integer                          :: i,j,k,nx
    real(8)                          :: e1,e2,s1,s2
    call this%set_spcs()
    allocate(this%parvdw(this%get_nvdw(),2),this%spcvdw(this%get_nvdw(),2))
    nx=1
    do i=1,this%nspcs
       do k=1,this%charmm%get_natp()
          if(this%spcs(i).eq.this%charmm%atp(k))then
             e1=this%charmm%prms_vdw(k,1)
             s1=this%charmm%prms_vdw(k,2)
          end if
       end do
       if(e1.ge.1.d-4.and.s1.ge.5.d-4)then
          do j=i,this%nspcs
             do k=1,this%charmm%get_natp()
                if(this%spcs(j).eq.this%charmm%atp(k))then
                   e2=this%charmm%prms_vdw(k,1)
                   s2=this%charmm%prms_vdw(k,2)
                end if
             end do
             if(e2.ge.1.d-4.and.s2.ge.5.d-4)then
                this%parvdw(nx,1)=sqrt(e1*e2)
                this%parvdw(nx,2)=s1+s2
                this%spcvdw(nx,1)=this%spcs(i)
                this%spcvdw(nx,2)=this%spcs(j)
                nx=nx+1
             end if
          end do
       end if
    end do
    call this%set_nvdw(nx-1)
  end subroutine set_parvdw

  subroutine set_parbnd(this)
    implicit none
    class(forcefield), intent(inout) :: this
    integer                          :: i,j,k,m,i1,i2,n1,n2
    logical                          :: check
    allocate(this%parbnd(this%get_nmol(),this%get_bondmax(),2))
!    call this%forcefield_init()
    do i=1,this%get_nmol()
       do j=1,this%bondscnt(i)
          i1=this%molbond(i,j,1)
          i2=this%molbond(i,j,2)
          n1=1
          n2=1
          do k=1,this%charmm%get_natp()
             if(this%charmm%atp(k).eq.this%tpmol(i,i1))n1=k
             if(this%charmm%atp(k).eq.this%tpmol(i,i2))n2=k
          end do
          check=.true.
          do m=1,2
             this%parbnd(i,j,m)=this%charmm%prms_bonds(n1,n2,m)
             if(this%parbnd(i,j,m).lt.1.d-8)check=.false.
          end do
          if(check.eqv..false.)then
             do m=1,2
                this%parbnd(i,j,m)=this%charmm%prms_bonds(n2,n1,m)
             end do
          end if
       end do
    end do
  end subroutine set_parbnd

  subroutine set_parbend(this)
    implicit none
    class(forcefield), intent(inout) :: this
    integer                          :: i,j,k,m,i1,i2,i3,n1,n2,n3
    logical                          :: check
    allocate(this%parbend(this%get_nmol(),this%get_bendmax(),2))
    do i=1,this%get_nmol()
       do j=1,this%bendscnt(i)
          i1=this%molbend(i,j,1)
          i2=this%molbend(i,j,2)
          i3=this%molbend(i,j,3)
          n1=1
          n2=1
          n3=1
          do k=1,this%charmm%get_natp()
             if(this%charmm%atp(k).eq.this%tpmol(i,i1))n1=k
             if(this%charmm%atp(k).eq.this%tpmol(i,i2))n2=k
             if(this%charmm%atp(k).eq.this%tpmol(i,i3))n3=k
          end do
          check=.true.
          do m=1,2
             this%parbend(i,j,m)=this%charmm%prms_angles(n1,n2,n3,m)
             if(this%parbend(i,j,m).lt.1.d-8)check=.false.
          end do
          if(check.eqv..false.)then
             do m=1,2
                this%parbend(i,j,m)=this%charmm%prms_angles(n3,n2,n1,m)
             end do
          end if
       end do
    end do
  end subroutine set_parbend

  subroutine set_partors(this)
    implicit none
    class(forcefield), intent(inout) :: this
    integer                          :: i1,i2,i3,i4,i,j,k,m,n1,n2,n3,n4
    logical                          :: check
    allocate(this%partors(this%get_nmol(),this%get_torsmax(),3))
    do i=1,this%get_nmol()
       do j=1,this%torscnt(i)
          i1=this%moltors(i,j,1)
          i2=this%moltors(i,j,2)
          i3=this%moltors(i,j,3)
          i4=this%moltors(i,j,4)
          n1=1
          n2=1
          n3=1
          n4=1
          do k=1,this%charmm%get_natp()
             if(this%charmm%atp(k).eq.this%tpmol(i,i1))n1=k
             if(this%charmm%atp(k).eq.this%tpmol(i,i2))n2=k
             if(this%charmm%atp(k).eq.this%tpmol(i,i3))n3=k
             if(this%charmm%atp(k).eq.this%tpmol(i,i4))n4=k
          end do
          check=.true.
          do m=1,2
             this%partors(i,j,m)=this%charmm%prms_tors(n1,n2,n3,n4,m)
             if(this%partors(i,j,m).lt.1.d-8)check=.false.
          end do
          this%partors(i,j,3)=this%charmm%prms_tors(n1,n2,n3,n4,3)
          if(check.eqv..false.)then
             do m=1,3
                this%partors(i,j,m)=this%charmm%prms_tors(n4,n3,n2,n1,m)
             end do
          end if
       end do
    end do
  end subroutine set_partors

  subroutine set_paritors(this)
    implicit none
    class(forcefield), intent(inout) :: this
    integer                          :: i1,i2,i3,i4,i,j,nx,n1,n2,n3,n4
    allocate(this%paritors(this%get_nmol(),this%get_itorsmax(),4))
    allocate(this%itorscnt(this%get_nmol()))
    allocate(this%molitors(this%get_nmol(),this%get_itorsmax(),4))
    do i=1,this%get_nmol()
       nx=1
       do i1=1,this%nxmol(i)
          do j=1,this%charmm%get_natp()
             if(this%charmm%atp(j).eq.this%tpmol(i,i1))n1=j
          end do
          do i2=1,this%nxmol(i)
             if(i2.ne.i1)then
                do j=1,this%charmm%get_natp()
                   if(this%charmm%atp(j).eq.this%tpmol(i,i2))n2=j
                end do
                do i3=1,this%nxmol(i)
                   if(i3.ne.i2.and.i3.ne.i1)then
                      do j=1,this%charmm%get_natp()
                         if(this%charmm%atp(j).eq.this%tpmol(i,i3))n3=j
                      end do
                      do i4=1,this%nxmol(i)
                         if(i4.ne.i3.and.i4.ne.i2.and.i4.ne.i1)then
                            do j=1,this%charmm%get_natp()
                               if(this%charmm%atp(j).eq.this%tpmol(i,i4))n4=j
                            end do
                            if(this%charmm%prms_itors(n1,n2,n3,n4,1).gt.1.d-8)then
                               if(n1.eq.n4.and.n2.eq.n3)then
                                  if(i1.lt.i4)then
                                     do k=1,4
                                        this%paritors(i,nx,k)=&
                                             this%charmm%prms_itors(n1,n2,n3,n4,k)
                                     end do
                                     this%molitors(i,nx,1)=i1
                                     this%molitors(i,nx,2)=i2
                                     this%molitors(i,nx,3)=i3
                                     this%molitors(i,nx,4)=i4
                                     nx=nx+1
                                  end if
                               else
                                  do k=1,4
                                     this%paritors(i,nx,k)=&
                                          this%charmm%prms_itors(n1,n2,n3,n4,k)
                                  end do
                                  this%molitors(i,nx,1)=i1
                                  this%molitors(i,nx,2)=i2
                                  this%molitors(i,nx,3)=i3
                                  this%molitors(i,nx,4)=i4
                                  nx=nx+1
                               end if
                            end if
                         end if
                      end do
                   end if
                end do
             end if
          end do
       end do
       this%itorscnt(i)=nx-1
    end do
  end subroutine set_paritors

  subroutine set_extra_parvdw(this)
    implicit none
    class(forcefield), intent(inout) :: this
    integer                          :: nvdw
    real(8)                          :: p1,p2
    character(6)                     :: spcs1,spcs2
    character(16)                    :: key
    character(5)                     :: tvdw
    nvdw=this%get_nvdw()
!    check=.true.
!    do while(check)
!       read(5,*,end=2)key
!       if(key.eq.'&FORCE_FIELD'.or.key.eq.'&force_field')check=.false.
!    end do
!    check=.true.
!    do while (check)
!1   read(5,*,end=3)key
!    if(key.ne.'&FORCE_FIELD')goto 1
!    do while (key.ne.'&END')
!       read(5,*)key
!       if(key.eq.'vdw')then
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
             goto 1
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
                goto 1
             end if
          end do
       end do
       goto 3
1      continue
    end do
!       end if
!       if(key.eq.'&END_FORCE_FIELD'.or.key.eq.'&end_force_field')check=.false.
!    end do
       call this%set_nvdw(nvdw)
    rewind(5)
    return
3   write(6,*)'ERROR: The type does not match with that defined in the TOPOLOGY file!'
    write(6,*)'Hint: Check the input in the &FORCE_FIELD section.'
    stop
  end subroutine set_extra_parvdw

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
       case('icharmm')
          read(5,*)i2,this%moltors(i3,i2,1),this%moltors(i3,i2,2),&
               this%moltors(i3,i2,3),this%moltors(i3,i2,4),&
               this%ttors(i3,i2),(this%partors(i3,i2,l),l=1,2)
       case('harm')
          read(5,*)i2,this%moltors(i3,i2,1),this%moltors(i3,i2,2),&
               this%moltors(i3,i2,3),this%moltors(i3,i2,4),this%ttors(i3,i2),&
               (this%partors(i3,i2,l),l=1,2)
       end select
    end do
  end subroutine set_extra_dihedrals

  subroutine set_coulop(this,coulop)
    implicit none
    class(forcefield), intent(inout) :: this
    character(4), intent(in)         :: coulop
    this%coulop=coulop
  end subroutine set_coulop

  subroutine set_coulop2(this)
    implicit none
    class(forcefield), intent(inout) :: this
    real(8)                          :: fscsalpha
    character(4)                     :: coulop
    character(13)                    :: key
1   read(5,*,end=2)key
    if(key.ne.'&FORCE_FIELD')goto 1
    do while (key.ne.'&END')
       read(5,*)key
       if(key.eq.'electrostatic')then
          backspace(5)
          read(5,*)key,coulop
          this%coulop=coulop
          if(coulop.eq.'fscs')then
             backspace(5)
             read(5,*)key,coulop,fscsalpha
             this%fscsalpha=fscsalpha
          end if
          goto 2
       end if
    end do
2   rewind(5)
  end subroutine set_coulop2

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

end module forcefield_module
