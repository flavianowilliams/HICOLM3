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
module input_module
  !*******************************************************************************************
  !*******************************************************************************************

  use forcefield_module

  implicit none

  private
  public :: input

  type, extends(forcefield) :: input
     integer, private       :: nstep
     integer, private       :: nrelax
     integer, private       :: nframes
     integer, private       :: nhist
     real(8), private       :: timestep
     real(8), private       :: press
     real(8), private       :: temp
     real(8), private       :: rcutoff
     real(8), private       :: drcutoff
     real(8), private       :: pstat
     real(8), private       :: tstat
     real(8), private       :: bfactor
     real(8), private       :: checkenergy
     real(8), private       :: tolerance
     character(8), private  :: restart
     character(3), private  :: ensble
     character(9), private  :: ensble_mt
   contains
     procedure :: set_input
     procedure :: set_inopt
     procedure :: set_topology
     procedure :: set_molecules
     procedure :: set_latticevectors
     procedure :: set_restart
     procedure :: get_restart
     procedure :: set_nstep
     procedure :: get_nstep
     procedure :: set_nrelax
     procedure :: get_nrelax
     procedure :: set_nframes
     procedure :: get_nframes
     procedure :: set_nhist
     procedure :: get_nhist
     procedure :: set_timestep
     procedure :: get_timestep
     procedure :: set_press
     procedure :: get_press
     procedure :: set_temp
     procedure :: get_temp
     procedure :: set_rcutoff
     procedure :: get_rcutoff
     procedure :: set_drcutoff
     procedure :: get_drcutoff
     procedure :: set_pstat
     procedure :: get_pstat
     procedure :: set_tstat
     procedure :: get_tstat
     procedure :: set_ensble
     procedure :: get_ensble
     procedure :: set_ensble_mt
     procedure :: get_ensble_mt
     procedure :: set_bfactor
     procedure :: get_bfactor
     procedure :: set_checkenergy
     procedure :: get_checkenergy
     procedure :: set_tolerance
     procedure :: get_tolerance
  end type input

contains

  subroutine set_input(this)
    implicit none
    class(input), intent(inout) :: this
    real(8)                     :: tstat,pstat
    character(11)               :: key
    character(3)                :: ensble
    character(9)                :: ensble_mt
1   read(5,*,end=2)key
    if(key.ne.'&MD')goto 1
    do while (key.ne.'&END')
       read(5,*)key
       if(key.eq.'nstep')then
          backspace(5)
          read(5,*)key,this%nstep
       end if
       if(key.eq.'nrelax')then
          backspace(5)
          read(5,*)key,this%nrelax
       end if
       if(key.eq.'nframes')then
          backspace(5)
          read(5,*)key,this%nframes
       end if
       if(key.eq.'timestep')then
          backspace(5)
          read(5,*)key,this%timestep
       end if
       if(key.eq.'pressure')then
          backspace(5)
          read(5,*)key,this%press
       end if
       if(key.eq.'temperature')then
          backspace(5)
          read(5,*)key,this%temp
       end if
       if(key.eq.'rcutoff')then
          backspace(5)
          read(5,*)key,this%rcutoff,this%drcutoff
       end if
       if(key.eq.'restart')then
          backspace(5)
          read(5,*)key,this%restart
       end if
       if(key.eq.'ensemble')then
          backspace(5)
          read(5,*)key,ensble
          if(ensble.eq.'nve')then
             this%ensble=ensble
          elseif(ensble.eq.'nvt')then
             backspace(5)
             read(5,*)key,ensble,ensble_mt
             if(ensble_mt.eq.'berendsen'.or.ensble_mt.eq.'hoover')then
                backspace(5)
                read(5,*)key,ensble,ensble_mt,tstat
                this%ensble=ensble
                this%ensble_mt=ensble_mt
                this%tstat=tstat
             end if
          elseif(ensble.eq.'npt')then
             backspace(5)
             read(5,*)key,ensble,ensble_mt
             if(ensble_mt.eq.'berendsen'.or.ensble_mt.eq.'hoover')then
                backspace(5)
                read(5,*)key,ensble,ensble_mt,tstat,pstat
                this%ensble=ensble
                this%ensble_mt=ensble_mt
                this%tstat=tstat
                this%pstat=pstat
             end if
          end if
       end if
    end do
2   rewind(5)
  end subroutine set_input

  subroutine set_inopt(this)
    implicit none
    class(input), intent(inout) :: this
    integer                     :: nstep
    real(8)                     :: tolerance
    character(13)               :: key
1   read(5,*,end=2)key
    if(key.ne.'&OPTIMIZATION')goto 1
    do while (key.ne.'&END')
       read(5,*)key
       if(key.eq.'nstep')then
          backspace(5)
          read(5,*)key,nstep
          call this%set_nstep(nstep)
       end if
       if(key.eq.'tolerance')then
          backspace(5)
          read(5,*)key,tolerance
          call this%set_tolerance(tolerance)
       end if
    end do
2   rewind(5)
  end subroutine set_inopt

  subroutine set_topology(this)
    implicit none
    class(input), intent(inout) :: this
    integer                     :: nmol,bondmax,bendmax,torsmax,nspcs,nvdw,i1,i,j,k,nx
    real(8)                     :: f1,f2,fscsalpha
    character(2)                :: mtd
    character(4)                :: coulop
    character(6)                :: tvdw,spcvdw1,spcvdw2
    character(7)                :: ttors
    character(10), allocatable  :: namemol(:)
    logical                     :: chk
    open(11,file='TOPOLOGY',status='old')
    read(11,'(1x,a2)')mtd
    if(mtd.eq.'MM')then
       read(11,'(1x,a4)')coulop
       if(coulop.eq.'fscs')then
          backspace(11)
          read(11,'(1x,a4,1x,f5.3)')coulop,fscsalpha
          call this%set_fscsalpha(fscsalpha)
       end if
       call this%set_coulop(coulop)
       read(11,'(1x,i2,4(1x,i3))')nmol,bondmax,bendmax,torsmax,nspcs
       call this%set_nmol(nmol)
       allocate(this%zatmol(nmol,this%get_natom()),this%qatmol(nmol,this%get_natom()))
       allocate(this%tpmol(nmol,this%get_natom()))
       allocate(this%sf_coul(nmol),this%sf_vdw(nmol))
       allocate(this%massmol(this%get_nmol(),this%get_natom()))
       allocate(this%bondscnt(nmol),this%bendscnt(nmol),this%torscnt(nmol))
       allocate(this%tbonds(nmol,bondmax),this%tbends(nmol,bendmax),this%ttors(nmol,torsmax))
       allocate(this%molbond(nmol,bondmax,2),this%molbend(nmol,bendmax,3))
       allocate(this%moltors(nmol,torsmax,4))
       allocate(this%parbnd(nmol,bondmax,2),this%parbend(nmol,bendmax,2))
       allocate(this%partors(nmol,torsmax,4))
       allocate(namemol(this%get_nmol()))
       do i=1,this%get_nmol()
          do j=1,bondmax
             do k=1,2
                this%parbnd(i,j,k)=0.d0
             end do
          end do
          do j=1,bendmax
             do k=1,2
                this%parbend(i,j,k)=0.d0
             end do
          end do
          do j=1,torsmax
             do k=1,4
                this%partors(i,j,k)=0.d0
             end do
          end do
       end do
       do i=1,this%get_nmol()
          read(11,'(1x,a10,2(1x,f8.6))')namemol(i)
          do j=1,this%get_nmol()
             if(namemol(i).eq.this%namemol(j))then
                nx=j
                goto 1
             end if
          end do
          goto 2
1         backspace(11)
          read(11,'(1x,a10,2(1x,f8.6))')this%namemol(nx),this%sf_coul(nx),this%sf_vdw(nx)
          read(11,'(15(1x,i2))')(this%zatmol(nx,j),j=1,this%nxmol(nx))
          read(11,'(15(1x,a8))')(this%tpmol(nx,j),j=1,this%nxmol(nx))
          read(11,'(15(1x,f8.4))')(this%massmol(nx,j),j=1,this%nxmol(nx))
          read(11,'(15(1x,f8.4))')(this%qatmol(nx,j),j=1,this%nxmol(nx))
          read(11,'(7x,i3)')this%bondscnt(nx)
          do j=1,this%bondscnt(nx)
             read(11,'(2(1x,i3),1x,a6,2(1x,f9.4))')this%molbond(nx,j,1),this%molbond(nx,j,2),&
                  this%tbonds(nx,j),(this%parbnd(nx,j,k),k=1,2)
          end do
          read(11,'(7x,i3)')this%bendscnt(nx)
          do j=1,this%bendscnt(nx)
             read(11,'(3(1x,i3),1x,a6,2(1x,f9.4))')(this%molbend(nx,j,k),k=1,3),&
                  this%tbends(nx,j),(this%parbend(nx,j,k),k=1,2)
          end do
          read(11,'(11x,i3)')this%torscnt(nx)
          do j=1,this%torscnt(nx)
             read(11,'(17x,a7)')ttors
             select case(ttors)
             case('charmm')
                backspace(11)
                read(11,'(4(1x,i3),1x,a7,1x,f8.4,1x,i2,1x,f8.4)')&
                     (this%moltors(nx,j,k),k=1,4),this%ttors(nx,j),f1,i1,f2
                this%partors(nx,j,1)=f1
                this%partors(nx,j,2)=i1
                this%partors(nx,j,3)=f2
             case('icharmm')
                backspace(11)
                read(11,'(4(1x,i3),1x,a7,1x,f8.4,1x,f8.4)')&
                     (this%moltors(nx,j,k),k=1,4),this%ttors(nx,j),f1,f2
                this%partors(nx,j,1)=f1
                this%partors(nx,j,2)=f2
             case('harm')
                backspace(11)
                read(11,'(4(1x,i3),1x,a7,2(1x,f9.4))')(this%moltors(nx,j,k),k=1,4),&
                     this%ttors(nx,j),this%partors(nx,j,1),this%partors(nx,j,2)
             end select
          end do
       end do
       read(11,'(4x,2(1x,i3))')nvdw
       allocate(this%tvdw(nvdw))
       call this%set_nspcs()
       call this%set_spcs()
       allocate(this%spcvdw(nvdw,2),this%parvdw(nvdw,2))
       do i=1,nvdw
          read(11,'(2(1x,a6),1x,a6)')spcvdw1,spcvdw2,tvdw
          backspace(11)
          select case(tvdw)
          case('charmm')
             read(11,'(2(1x,a6),1x,a6,2(1x,f9.4))')&
                  spcvdw1,spcvdw2,tvdw,(this%parvdw(i,j),j=1,2)
          case('lj')
             read(11,'(2(1x,a6),1x,a6,2(1x,f9.4))')&
                  spcvdw1,spcvdw2,tvdw,(this%parvdw(i,j),j=1,2)
          end select
          this%spcvdw(i,1)=spcvdw1
          this%spcvdw(i,2)=spcvdw2
          this%tvdw(i)=tvdw
       end do
       call this%set_nvdw(nvdw)
    end if
    do i=1,this%get_nmol()
       chk=.false.
       do j=1,this%get_nmol()
          if(this%namemol(i).eq.namemol(j))chk=.true.
       end do
    end do
    if(chk.eqv..false.)goto 3
    return
2   write(6,*)'ERROR: There is a molecule that does not belong to the physical system!'
    write(6,*)'Hint: Check the TOPOLOGY file.'
    stop
3   write(6,*)'ERROR: There is a lack in the description of molecules in the force field!'
    write(6,*)'Hint: Check the TOPOLOGY file.'
    stop
  end subroutine set_topology

  subroutine set_molecules(this)
    implicit none
    class(input), intent(inout) :: this
    integer                     :: nmol,natom,i,j,k,nx
    read(10,'(1x,2i5)')nmol,natom
    call this%set_nmol(nmol)
    allocate(this%namemol(nmol),this%ntmol(nmol),this%nxmol(nmol))
    allocate(this%xa(natom),this%ya(natom),this%za(natom))
    nx=1
    do i=1,this%get_nmol()
       read(10,'(1x,a10,2(1x,i5))')this%namemol(i),this%ntmol(i),this%nxmol(i)
       do j=1,this%ntmol(i)
          do k=1,this%nxmol(i)
             read(10,'(3f16.8)')this%xa(nx),this%ya(nx),this%za(nx)
             nx=nx+1
          end do
       end do
    end do
  end subroutine set_molecules

  subroutine set_latticevectors(this)
    implicit none
    class(input), intent(inout) :: this
    integer                     :: i,j
    open(10,file='SYSTEM',status='old')
    do i=1,3
       read(10,'(3f16.8)')(this%v(i,j),j=1,3)
    end do
  end subroutine set_latticevectors

  subroutine set_nstep(this,nstep)
    implicit none
    class(input), intent(inout) :: this
    integer, intent(in)         :: nstep
    this%nstep=nstep
  end subroutine set_nstep

  integer function get_nstep(this)
    implicit none
    class(input), intent(in) :: this
    get_nstep=this%nstep
  end function get_nstep

  subroutine set_nrelax(this,nrelax)
    implicit none
    class(input), intent(inout) :: this
    integer, intent(in)         :: nrelax
    this%nrelax=nrelax
  end subroutine set_nrelax

  integer function get_nrelax(this)
    implicit none
    class(input), intent(in) :: this
    get_nrelax=this%nrelax
  end function get_nrelax

  subroutine set_nframes(this,nframes)
    implicit none
    class(input), intent(inout) :: this
    integer, intent(in)         :: nframes
    this%nframes=nframes
  end subroutine set_nframes

  integer function get_nframes(this)
    implicit none
    class(input), intent(in) :: this
    get_nframes=this%nframes
  end function get_nframes

  subroutine set_nhist(this)
    implicit none
    class(input), intent(inout) :: this
    this%nhist=nint(real(this%get_nstep()-this%get_nrelax())/this%get_nframes())
  end subroutine set_nhist

  integer function get_nhist(this)
    implicit none
    class(input), intent(in) :: this
    get_nhist=this%nhist
  end function get_nhist

  subroutine set_timestep(this,timestep)
    implicit none
    class(input), intent(inout) :: this
    real(8), intent(in)         :: timestep
    this%timestep=timestep
  end subroutine set_timestep

  double precision function get_timestep(this)
    implicit none
    class(input), intent(in) :: this
    get_timestep=this%timestep
  end function get_timestep

  subroutine set_press(this,press)
    implicit none
    class(input), intent(inout) :: this
    real(8), intent(in)         :: press
    this%press=press
  end subroutine set_press

  double precision function get_press(this)
    implicit none
    class(input), intent(in) :: this
    get_press=this%press
  end function get_press

  subroutine set_temp(this,temp)
    implicit none
    class(input), intent(inout) :: this
    real(8), intent(in)         :: temp
    this%temp=temp
  end subroutine set_temp

  double precision function get_temp(this)
    implicit none
    class(input), intent(in) :: this
    get_temp=this%temp
  end function get_temp

  subroutine set_rcutoff(this,rcutoff)
    implicit none
    class(input), intent(inout) :: this
    real(8), intent(in)         :: rcutoff
    this%rcutoff=rcutoff
  end subroutine set_rcutoff

  double precision function get_rcutoff(this)
    implicit none
    class(input), intent(in) :: this
    get_rcutoff=this%rcutoff
  end function get_rcutoff

  subroutine set_drcutoff(this,drcutoff)
    implicit none
    class(input), intent(inout) :: this
    real(8), intent(in)         :: drcutoff
    this%drcutoff=drcutoff
  end subroutine set_drcutoff

  double precision function get_drcutoff(this)
    implicit none
    class(input), intent(in) :: this
    get_drcutoff=this%drcutoff
  end function get_drcutoff

  subroutine set_pstat(this,pstat)
    implicit none
    class(input), intent(inout) :: this
    real(8), intent(in)         :: pstat
    this%pstat=pstat
  end subroutine set_pstat

  double precision function get_pstat(this)
    implicit none
    class(input), intent(in) :: this
    get_pstat=this%pstat
  end function get_pstat

  subroutine set_tstat(this,tstat)
    implicit none
    class(input), intent(inout) :: this
    real(8), intent(in)         :: tstat
    this%tstat=tstat
  end subroutine set_tstat

  double precision function get_tstat(this)
    implicit none
    class(input), intent(in) :: this
    get_tstat=this%tstat
  end function get_tstat

  subroutine set_ensble(this,ensble)
    implicit none
    class(input), intent(inout) :: this
    character(3), intent(in)    :: ensble
    this%ensble=ensble
  end subroutine set_ensble

  character(3) function get_ensble(this)
    implicit none
    class(input), intent(in) :: this
    get_ensble=this%ensble
  end function get_ensble

  subroutine set_ensble_mt(this,ensble_mt)
    implicit none
    class(input), intent(inout) :: this
    character(9), intent(in)    :: ensble_mt
    this%ensble_mt=ensble_mt
  end subroutine set_ensble_mt

  character(9) function get_ensble_mt(this)
    implicit none
    class(input), intent(in) :: this
    get_ensble_mt=this%ensble_mt
  end function get_ensble_mt

  subroutine set_restart(this,restart)
    implicit none
    class(input), intent(inout) :: this
    character(8), intent(in)    :: restart
    this%restart=restart
  end subroutine set_restart

  character(8) function get_restart(this)
    implicit none
    class(input), intent(in) :: this
    get_restart=this%restart
  end function get_restart

  subroutine set_bfactor(this,bfactor)
    implicit none
    class(input), intent(inout) :: this
    real(8), intent(in)         :: bfactor
    this%bfactor=bfactor
  end subroutine set_bfactor

  double precision function get_bfactor(this)
    implicit none
    class(input), intent(in) :: this
    get_bfactor=this%bfactor
  end function get_bfactor

  subroutine set_checkenergy(this,checkenergy)
    implicit none
    class(input), intent(inout) :: this
    real(8), intent(in)         :: checkenergy
    this%checkenergy=checkenergy
  end subroutine set_checkenergy

  double precision function get_checkenergy(this)
    implicit none
    class(input), intent(in) :: this
    get_checkenergy=this%checkenergy
  end function get_checkenergy

  subroutine set_tolerance(this,tolerance)
    implicit none
    class(input), intent(inout) :: this
    real(8), intent(in)         :: tolerance
    this%tolerance=tolerance
  end subroutine set_tolerance

  double precision function get_tolerance(this)
    implicit none
    class(input), intent(inout) :: this
    get_tolerance=this%tolerance
  end function get_tolerance

end module input_module
