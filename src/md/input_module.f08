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

  use molecule_module

  implicit none

!  integer i,j,k,l

  private
  public :: input

  type, extends(molecule) :: input
     integer, private      :: nstep
     integer, private      :: nrelax
     integer, private      :: nframes
     real(8), private      :: timestep
     real(8), private      :: press
     real(8), private      :: temp
     real(8), private      :: rcutoff
     real(8), private      :: drcutoff
     real(8), private      :: pstat
     real(8), private      :: tstat
     character(3), private :: ensble
     character(9), private :: ensble_mt
   contains
     procedure :: set_input
     procedure :: set_nstep
     procedure :: get_nstep
     procedure :: set_nrelax
     procedure :: get_nrelax
     procedure :: set_nframes
     procedure :: get_nframes
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
       if(key.eq.'ensemble')then
          backspace(5)
          read(5,*)key,ensble
          if(ensble.eq.'nve')then
             this%ensble=ensble
             print*,this%ensble
             stop
          elseif(.eq.'nvt')then
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


       if(key.eq.'ensemble')then
          backspace(5)
          read(5,*)key,this%
       end if
    end do
2   rewind(5)
  end subroutine set_input

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

end module input_module
