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
     integer, private       :: restart
     real(8), private       :: timestep
     real(8), private       :: press
     real(8), private       :: temp
     real(8), private       :: rcutoff
     real(8), private       :: drcutoff
     real(8), private       :: pstat
     real(8), private       :: tstat
     real(8), private       :: bfactor
     character(3), private  :: ensble
     character(9), private  :: ensble_mt
   contains
     procedure :: set_input
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
     procedure :: convert_units
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

  subroutine set_topology(this)
    implicit none
    class(input), intent(inout) :: this
    integer                     :: nmol,bondmax,bendmax,torsmax,nspcs,nvdw,i1,i2,i,j,k
    real(8)                     :: f1,f2
    character(2)                :: mtd,spcvdw1,spcvdw2
    character(4)                :: coulop
    character(5)                :: ttors,tvdw
    open(11,file='TOPOLOGY',status='old')
    read(11,'(1x,a2)')mtd
    if(mtd.eq.'MM')then
       read(11,'(1x,a4)')coulop
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
       do i=1,this%get_nmol()
          read(11,'(1x,a10,2(1x,f8.6))')this%namemol(i),this%sf_coul(i),this%sf_vdw(i)
          read(11,'(15(1x,i2))')(this%zatmol(i,j),j=1,this%nxmol(i))
          read(11,'(15(1x,a2))')(this%tpmol(i,j),j=1,this%nxmol(i))
          read(11,'(15(1x,f8.4))')(this%massmol(i,j),j=1,this%nxmol(i))
          read(11,'(15(1x,f8.4))')(this%qatmol(i,j),j=1,this%nxmol(i))
          read(11,'(7x,i3)')this%bondscnt(i)
          do j=1,this%bondscnt(i)
             read(11,'(2(1x,i3),1x,a5,2(1x,f9.4))')this%molbond(i,j,1),this%molbond(i,j,2),&
                  this%tbonds(i,j),(this%parbnd(i,j,k),k=1,2)
          end do
          read(11,'(7x,i3)')this%bendscnt(i)
          do j=1,this%bendscnt(i)
             read(11,'(3(1x,i3),1x,a5,2(1x,f9.4))')(this%molbend(i,j,k),k=1,3),&
                  this%tbends(i,j),(this%parbend(i,j,k),k=1,2)
          end do
          read(11,'(11x,i3)')this%torscnt(i)
          do j=1,this%torscnt(i)
             read(11,'(17x,a5)')ttors
             select case(ttors)
             case('amber')
                backspace(11)
                read(11,'(4(1x,i3),1x,a5,2x,i2,f8.2,f8.1,1x,i2)')&
                     (this%moltors(i,j,k),k=1,4),this%ttors(i,j),i1,f1,f2,i2
                this%partors(i,j,1)=i1
                this%partors(i,j,2)=f1
                this%partors(i,j,3)=f2
                this%partors(i,j,4)=i2
             case('harm')
                backspace(11)
                read(11,'(4(1x,i3),1x,a5,2(1x,f9.4))')(this%moltors(i,j,k),k=1,4),&
                     this%ttors(i,j),this%partors(i,j,1),this%partors(i,j,2)
             end select
          end do
       end do
       read(11,'(4x,2(1x,i3))')nvdw
       allocate(this%tvdw(nvdw))
       call this%set_nspcs()
       call this%set_spcs()
       allocate(this%spcvdw(nvdw,2),this%parvdw(nvdw,2))
       do i=1,nvdw
          read(11,'(2(1x,a2),1x,a5)')spcvdw1,spcvdw2,tvdw
          backspace(11)
          select case(tvdw)
          case('amber')
             read(11,'(2(1x,a2),1x,a5,2(1x,f9.4))')&
                  spcvdw1,spcvdw2,tvdw,(this%parvdw(i,j),j=1,2)
          case('lj')
             read(11,'(2(1x,a2),1x,a5,2(1x,f9.4))')&
                  spcvdw1,spcvdw2,tvdw,(this%parvdw(i,j),j=1,2)
          end select
          this%spcvdw(i,1)=spcvdw1
          this%spcvdw(i,2)=spcvdw2
          this%tvdw(i)=tvdw
       end do
       call this%set_nvdw(nvdw)
    end if
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
    integer, intent(in)         :: restart
    this%restart=restart
  end subroutine set_restart

  integer function get_restart(this)
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

  subroutine convert_units(this)
    implicit none
    class(input), intent(inout) :: this
    integer                     :: i,j
    do i=1,3
       do j=1,3
          this%v(i,j)=this%v(i,j)/this%get_rconv()
       end do
    end do
    call this%set_a(this%get_a()/this%get_rconv())
    call this%set_b(this%get_b()/this%get_rconv())
    call this%set_c(this%get_c()/this%get_rconv())
    do i=1,this%get_natom()
       this%xa(i)=this%xa(i)/this%get_rconv()
       this%ya(i)=this%ya(i)/this%get_rconv()
       this%za(i)=this%za(i)/this%get_rconv()
    end do
    do i=1,this%get_nmol()
       do j=1,this%nxmol(i)
          this%qatmol(i,j)=this%qatmol(i,j)/this%get_elconv()
          this%massmol(i,j)=this%massmol(i,j)/this%get_mconv()
       end do
    end do
    do i=1,this%get_nmol()
       do j=1,this%bondscnt(i)
          select case(this%tbonds(i,j))
          case('amber')
             this%parbnd(i,j,1)=this%parbnd(i,j,1)/(this%get_econv()/this%get_rconv()**2)
             this%parbnd(i,j,2)=this%parbnd(i,j,2)/this%get_rconv()
          case('harm')
             this%parbnd(i,j,1)=this%parbnd(i,j,1)/(this%get_econv()/this%get_rconv()**2)
             this%parbnd(i,j,2)=this%parbnd(i,j,2)/this%get_rconv()
          end select
       end do
       do j=1,this%bendscnt(i)
          select case(this%tbends(i,j))
          case('amber')
             this%parbend(i,j,1)=this%parbend(i,j,1)/(this%get_econv()/this%get_aconv()**2)
             this%parbend(i,j,2)=this%parbend(i,j,2)/this%get_aconv()
          case('harm')
             this%parbend(i,j,1)=this%parbend(i,j,1)/(this%get_econv()/this%get_aconv()**2)
             this%parbend(i,j,2)=this%parbend(i,j,2)/this%get_aconv()
          end select
       end do
       do j=1,this%torscnt(i)
          select case(this%ttors(i,j))
          case('amber')
             this%partors(i,j,2)=this%partors(i,j,2)/this%get_econv()
             this%partors(i,j,3)=this%partors(i,j,3)/this%get_aconv()
          case('harm')
             this%partors(i,j,1)=this%partors(i,j,1)/(this%get_econv()/this%get_aconv()**2)
             this%partors(i,j,2)=this%partors(i,j,2)/this%get_aconv()
          end select
       end do
    end do
    do i=1,this%get_nvdw()
       select case(this%tvdw(i))
       case('amber')
          this%parvdw(i,1)=this%parvdw(i,1)/this%get_econv()
          this%parvdw(i,2)=this%parvdw(i,2)/this%get_rconv()
       case('lj')
          this%parvdw(i,1)=this%parvdw(i,1)/this%get_econv()
          this%parvdw(i,2)=this%parvdw(i,2)/this%get_rconv()
       end select
    end do
    call this%set_timestep(this%get_timestep()/this%get_tconv())
    call this%set_temp(this%get_temp()/this%get_teconv())
    call this%set_press(this%get_press()/this%get_pconv())
    call this%set_tstat(this%get_tstat()/this%get_tconv())
    call this%set_pstat(this%get_pstat()/this%get_tconv())
    call this%set_rcutoff(this%get_rcutoff()/this%get_rconv())
    call this%set_drcutoff(this%get_drcutoff()/this%get_rconv())
    call this%set_bfactor(this%get_bfactor()*this%get_pconv())
  end subroutine convert_units

end module input_module
