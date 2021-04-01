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
module ensemble_module
  !*******************************************************************************************
  !*******************************************************************************************

  use thermodynamics_module

  implicit none

  private
  public :: ensemble

  type, extends(thermodynamics) :: ensemble
   contains
     procedure :: set_nve
     procedure :: set_nvt_berendsen
     procedure :: set_npt_berendsen
     procedure :: set_nvt_nosehoover
     procedure :: set_npt_nosehoover
     procedure :: check_lattice
  end type ensemble

contains

  subroutine set_nve(this)
    implicit none
    class(ensemble), intent(inout) :: this
    integer                        :: i
    do i=1,this%get_natom()
       this%vax(i)=this%vax(i)+this%fax(i)*(0.5d0*this%get_timestep())/this%mass(i)
       this%vay(i)=this%vay(i)+this%fay(i)*(0.5d0*this%get_timestep())/this%mass(i)
       this%vaz(i)=this%vaz(i)+this%faz(i)*(0.5d0*this%get_timestep())/this%mass(i)
       this%xa(i)=this%xa(i)+this%vax(i)*this%get_timestep()
       this%ya(i)=this%ya(i)+this%vay(i)*this%get_timestep()
       this%za(i)=this%za(i)+this%vaz(i)*this%get_timestep()
    end do
    call this%ccp()
    call this%set_forcefield()
    do i=1,this%get_natom()
       this%vax(i)=this%vax(i)+this%fax(i)*(0.5d0*this%get_timestep())/this%mass(i)
       this%vay(i)=this%vay(i)+this%fay(i)*(0.5d0*this%get_timestep())/this%mass(i)
       this%vaz(i)=this%vaz(i)+this%faz(i)*(0.5d0*this%get_timestep())/this%mass(i)
    end do
    call this%set_ekinetic()
    call this%set_etotal()
    call this%set_temperature()
    call this%set_pressure()
  end subroutine set_nve

  subroutine set_nvt_berendsen(this)
    implicit none
    class(ensemble), intent(inout) :: this
    integer                        :: i
    real(8)                        :: qui
    do i=1,this%get_natom()
       this%vax(i)=this%vax(i)+this%fax(i)*(0.5d0*this%get_timestep())/this%mass(i)
       this%vay(i)=this%vay(i)+this%fay(i)*(0.5d0*this%get_timestep())/this%mass(i)
       this%vaz(i)=this%vaz(i)+this%faz(i)*(0.5d0*this%get_timestep())/this%mass(i)
       this%xa(i)=this%xa(i)+this%vax(i)*this%get_timestep()
       this%ya(i)=this%ya(i)+this%vay(i)*this%get_timestep()
       this%za(i)=this%za(i)+this%vaz(i)*this%get_timestep()
    end do
    call this%ccp()
    call this%set_forcefield()
    do i=1,this%get_natom()
       this%vax(i)=this%vax(i)+this%fax(i)*(0.5d0*this%get_timestep())/this%mass(i)
       this%vay(i)=this%vay(i)+this%fay(i)*(0.5d0*this%get_timestep())/this%mass(i)
       this%vaz(i)=this%vaz(i)+this%faz(i)*(0.5d0*this%get_timestep())/this%mass(i)
    end do
    call this%set_ekinetic()
    qui=sqrt(1.d0+this%get_timestep()*(this%get_sigma()/&
         this%get_ekinetic()-1.d0)/this%get_tstat())
    do i=1,this%get_natom()
       this%vax(i)=this%vax(i)*qui
       this%vay(i)=this%vay(i)*qui
       this%vaz(i)=this%vaz(i)*qui
    end do
    call this%set_ekinetic()
    call this%set_etotal()
    call this%set_temperature()
    call this%set_pressure()
  end subroutine set_nvt_berendsen

  subroutine set_nvt_nosehoover(this)
    implicit none
    class(ensemble), intent(inout) :: this
    integer                        :: i
    real(8)                        :: qui
  end subroutine set_nvt_nosehoover

  subroutine set_npt_berendsen(this)
    implicit none
    class(ensemble), intent(inout) :: this
    integer                        :: i,j
    real(8)                        :: qui,eta
    eta=(1.d0-this%get_bfactor()*this%get_timestep()*&
         (this%get_press()-this%get_pressure())/this%get_pstat())**(1.d0/3.d0)
    do i=1,this%get_natom()
       this%vax(i)=this%vax(i)+this%fax(i)*0.5d0*this%get_timestep()/this%mass(i)
       this%vay(i)=this%vay(i)+this%fay(i)*0.5d0*this%get_timestep()/this%mass(i)
       this%vaz(i)=this%vaz(i)+this%faz(i)*0.5d0*this%get_timestep()/this%mass(i)
       this%xa(i)=this%xa(i)*eta+this%vax(i)*this%get_timestep()
       this%ya(i)=this%ya(i)*eta+this%vay(i)*this%get_timestep()
       this%za(i)=this%za(i)*eta+this%vaz(i)*this%get_timestep()
    end do
    do i=1,3
       do j=1,3
          this%v(i,j)=this%v(i,j)*eta
       end do
    end do
    call this%set_lattice_constants()
    call this%set_volume2(this%get_volume()*eta**3)
    call this%check_lattice()
    call this%ccp()
    call this%set_forcefield()
    do i=1,this%get_natom()
       this%vax(i)=this%vax(i)+this%fax(i)*0.5d0*this%get_timestep()/this%mass(i)
       this%vay(i)=this%vay(i)+this%fay(i)*0.5d0*this%get_timestep()/this%mass(i)
       this%vaz(i)=this%vaz(i)+this%faz(i)*0.5d0*this%get_timestep()/this%mass(i)
    end do
    call this%set_ekinetic()
    qui=sqrt(1.d0+this%get_timestep()*&
         (this%get_sigma()/this%get_ekinetic()-1.d0)/this%get_tstat())
    do i=1,this%get_natom()
       this%vax(i)=this%vax(i)*qui
       this%vay(i)=this%vay(i)*qui
       this%vaz(i)=this%vaz(i)*qui
    end do
    call this%set_ekinetic()
    call this%set_etotal()
    call this%set_temperature()
    call this%set_pressure()
  end subroutine set_npt_berendsen

  subroutine set_npt_nosehoover(this)
    implicit none
    class(ensemble), intent(inout) :: this
    integer                        :: i
  end subroutine set_npt_nosehoover

  subroutine check_lattice(this)
    implicit none
    class(ensemble), intent(inout) :: this
    if(this%get_a().le.2.d0*this%get_rcutoff().or.this%get_b().le.2.d0*this%get_rcutoff()&
         .or.this%get_c().le.2.d0*this%get_rcutoff())then
       write(6,*)'ERROR: The half-box size exceeds the cutoff radius!'
       stop
    end if
  end subroutine check_lattice

end module ensemble_module
