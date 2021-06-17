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
     real(8), private :: tfcnvt
     real(8), private :: tfcnpt
     real(8), private :: bfc
     real(8), private :: rcm(3)
     real(8), private :: vcm(3)
   contains
     procedure :: ensemble_init
     procedure :: set_nve
     procedure :: set_nvt_berendsen
     procedure :: set_npt_berendsen
     procedure :: set_nvt_nosehoover
     procedure :: set_npt_nosehoover
     procedure :: set_tfcnvt
     procedure :: get_tfcnvt
     procedure :: set_tfcnpt
     procedure :: get_tfcnpt
     procedure :: set_bfc
     procedure :: get_bfc
     procedure :: set_bfc2
     procedure :: check_lattice
     procedure :: check_energy
     procedure :: set_rcm
     procedure :: set_vcm
  end type ensemble

contains

  subroutine ensemble_init(this)
    implicit none
    class(ensemble), intent(inout) :: this
    this%tfcnvt=0.d0
    this%tfcnpt=0.d0
    this%bfc=0.d0
  end subroutine ensemble_init

  subroutine set_tfcnvt(this)
    implicit none
    class(ensemble), intent(inout) :: this
    this%tfcnvt=this%tfcnvt+0.5d0*this%get_timestep()*&
         (this%get_ekinetic()-this%get_sigma())/this%get_qmass()
  end subroutine set_tfcnvt

  double precision function get_tfcnvt(this)
    implicit none
    class(ensemble), intent(in) :: this
    get_tfcnvt=this%tfcnvt
  end function get_tfcnvt

  subroutine set_tfcnpt(this)
    implicit none
    class(ensemble), intent(inout) :: this
    this%tfcnpt=this%tfcnpt+0.125d0*this%get_timestep()*(2.0d0*this%get_ekinetic()&
         +this%get_pmass()*this%get_bfc()**2-2.0d0*this%get_sigma()-this%get_temp())&
         /this%get_qmass()
  end subroutine set_tfcnpt

  double precision function get_tfcnpt(this)
    implicit none
    class(ensemble), intent(in) :: this
    get_tfcnpt= this%tfcnpt
  end function get_tfcnpt

  subroutine set_bfc(this)
    implicit none
    class(ensemble), intent(inout) :: this
    this%bfc=this%bfc+0.75d0*this%get_timestep()&
         *(this%get_pressure()-this%get_press())*this%get_volume()/this%get_pmass()
  end subroutine set_bfc

  double precision function get_bfc(this)
    implicit none
    class(ensemble), intent(in) :: this
    get_bfc=this%bfc
  end function get_bfc

  subroutine set_bfc2(this)
    implicit none
    class(ensemble), intent(inout) :: this
    this%bfc=this%bfc*exp(-0.125d0*this%get_tfcnpt()*this%get_timestep())
  end subroutine set_bfc2

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
    call this%check_energy()
    call this%set_temperature()
    call this%set_pressure()
  end subroutine set_nve

  subroutine set_nvt_berendsen(this)
    implicit none
    class(ensemble), intent(inout) :: this
    integer                        :: i
    real(8)                        :: qui
    do i=1,this%get_natom()
       this%vax(i)=this%vax(i)+this%fax(i)*0.5d0*this%get_timestep()/this%mass(i)
       this%vay(i)=this%vay(i)+this%fay(i)*0.5d0*this%get_timestep()/this%mass(i)
       this%vaz(i)=this%vaz(i)+this%faz(i)*0.5d0*this%get_timestep()/this%mass(i)
       this%xa(i)=this%xa(i)+this%vax(i)*this%get_timestep()
       this%ya(i)=this%ya(i)+this%vay(i)*this%get_timestep()
       this%za(i)=this%za(i)+this%vaz(i)*this%get_timestep()
    end do
    call this%ccp()
    call this%set_forcefield()
    do i=1,this%get_natom()
       this%vax(i)=this%vax(i)+this%fax(i)*0.5d0*this%get_timestep()/this%mass(i)
       this%vay(i)=this%vay(i)+this%fay(i)*0.5d0*this%get_timestep()/this%mass(i)
       this%vaz(i)=this%vaz(i)+this%faz(i)*0.5d0*this%get_timestep()/this%mass(i)
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
    call this%check_energy()
    call this%set_temperature()
    call this%set_pressure()
  end subroutine set_nvt_berendsen

  subroutine set_nvt_nosehoover(this)
    implicit none
    class(ensemble), intent(inout) :: this
    integer                        :: i
    call this%set_tfcnvt()
    do i=1,this%get_natom()
       this%vax(i)=this%vax(i)*exp(-0.5d0*this%get_tfcnvt()*this%get_timestep())
       this%vay(i)=this%vay(i)*exp(-0.5d0*this%get_tfcnvt()*this%get_timestep())
       this%vaz(i)=this%vaz(i)*exp(-0.5d0*this%get_tfcnvt()*this%get_timestep())
    end do
    call this%set_ekinetic()
    call this%set_tfcnvt()
    do i=1,this%get_natom()
       this%vax(i)=this%vax(i)+0.5d0*this%fax(i)*this%get_timestep()/this%mass(i)
       this%vay(i)=this%vay(i)+0.5d0*this%fay(i)*this%get_timestep()/this%mass(i)
       this%vaz(i)=this%vaz(i)+0.5d0*this%faz(i)*this%get_timestep()/this%mass(i)
       this%xa(i)=this%xa(i)+this%get_timestep()*this%vax(i)
       this%ya(i)=this%ya(i)+this%get_timestep()*this%vay(i)
       this%za(i)=this%za(i)+this%get_timestep()*this%vaz(i)
    end do
    call this%ccp()
    call this%set_forcefield()
    do i=1,this%get_natom()
       this%vax(i)=this%vax(i)+0.5d0*this%fax(i)*this%get_timestep()/this%mass(i)
       this%vay(i)=this%vay(i)+0.5d0*this%fay(i)*this%get_timestep()/this%mass(i)
       this%vaz(i)=this%vaz(i)+0.5d0*this%faz(i)*this%get_timestep()/this%mass(i)
    end do
    call this%set_ekinetic()
    call this%set_tfcnvt()
    do i=1,this%get_natom()
       this%vax(i)=this%vax(i)*exp(-0.5d0*this%get_tfcnvt()*this%get_timestep())
       this%vay(i)=this%vay(i)*exp(-0.5d0*this%get_tfcnvt()*this%get_timestep())
       this%vaz(i)=this%vaz(i)*exp(-0.5d0*this%get_tfcnvt()*this%get_timestep())
    end do
    call this%set_ekinetic()
    call this%set_etotal()
    call this%check_energy()
    call this%set_temperature()
    call this%set_pressure()
    call this%set_tfcnvt()
  end subroutine set_nvt_nosehoover

  subroutine set_npt_berendsen(this)
    implicit none
    class(ensemble), intent(inout) :: this
    integer                        :: i,j
    real(8)                        :: qui,eta
    eta=(1.d0-this%get_bfactor()*this%get_timestep()*&
         (this%get_press()-this%get_pressure())/this%get_pstat())**(1.d0/3.d0)
    do i=1,this%get_natom()
       this%vax(i)=this%vax(i)+0.5d0*this%fax(i)*this%get_timestep()/this%mass(i)
       this%vay(i)=this%vay(i)+0.5d0*this%fay(i)*this%get_timestep()/this%mass(i)
       this%vaz(i)=this%vaz(i)+0.5d0*this%faz(i)*this%get_timestep()/this%mass(i)
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
       this%vax(i)=this%vax(i)+0.5d0*this%fax(i)*this%get_timestep()/this%mass(i)
       this%vay(i)=this%vay(i)+0.5d0*this%fay(i)*this%get_timestep()/this%mass(i)
       this%vaz(i)=this%vaz(i)+0.5d0*this%faz(i)*this%get_timestep()/this%mass(i)
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
    call this%check_energy()
    call this%set_temperature()
    call this%set_pressure()
  end subroutine set_npt_berendsen

  subroutine set_npt_nosehoover(this)
    implicit none
    class(ensemble), intent(inout) :: this
    integer                        :: i,j
    real(8)                        :: eta
    call this%set_rcm()
    call this%set_tfcnpt() !thermostat
    do i=1,this%get_natom()
       this%vax(i)=this%vax(i)*exp(-0.25d0*this%get_tfcnpt()*this%get_timestep())
       this%vay(i)=this%vay(i)*exp(-0.25d0*this%get_tfcnpt()*this%get_timestep())
       this%vaz(i)=this%vaz(i)*exp(-0.25d0*this%get_tfcnpt()*this%get_timestep())
    end do
    call this%set_ekinetic()
    call this%set_tfcnpt()
    call this%set_bfc2() !barostat
    call this%set_bfc()
    call this%set_bfc2()
    do i=1,this%get_natom()
       this%vax(i)=this%vax(i)*exp(-0.5d0*this%get_bfc()*this%get_timestep())
       this%vay(i)=this%vay(i)*exp(-0.5d0*this%get_bfc()*this%get_timestep())
       this%vaz(i)=this%vaz(i)*exp(-0.5d0*this%get_bfc()*this%get_timestep())
    end do
    call this%set_bfc2()
    call this%set_bfc()
    call this%set_bfc2()
    call this%set_ekinetic() !thermostat
    call this%set_tfcnpt()
    do i=1,this%get_natom()
       this%vax(i)=this%vax(i)*exp(-0.25d0*this%get_tfcnpt()*this%get_timestep())
       this%vay(i)=this%vay(i)*exp(-0.25d0*this%get_tfcnpt()*this%get_timestep())
       this%vaz(i)=this%vaz(i)*exp(-0.25d0*this%get_tfcnpt()*this%get_timestep())
    end do
    call this%set_ekinetic()
    call this%set_tfcnpt()
    do i=1,this%get_natom()
       this%vax(i)=this%vax(i)+0.5d0*this%fax(i)*this%get_timestep()/this%mass(i)
       this%vay(i)=this%vay(i)+0.5d0*this%fay(i)*this%get_timestep()/this%mass(i)
       this%vaz(i)=this%vaz(i)+0.5d0*this%faz(i)*this%get_timestep()/this%mass(i)
    end do
    eta=exp(this%get_bfc()*this%get_timestep())
    do i=1,3
       do j=1,3
          this%v(i,j)=this%v(i,j)*eta
       end do
    end do
    call this%set_lattice_constants()
    call this%set_volume2(this%get_volume()*eta**3)
    call this%check_lattice()
    do i=1,this%get_natom()
       this%xa(i)=(this%xa(i)-this%rcm(1))*eta+this%get_timestep()*this%vax(i)+this%rcm(1)
       this%ya(i)=(this%ya(i)-this%rcm(2))*eta+this%get_timestep()*this%vay(i)+this%rcm(2)
       this%za(i)=(this%za(i)-this%rcm(3))*eta+this%get_timestep()*this%vaz(i)+this%rcm(3)
    end do
    call this%ccp()
    call this%set_forcefield()
    do i=1,this%get_natom()
       this%vax(i)=this%vax(i)+0.5d0*this%fax(i)*this%get_timestep()/this%mass(i)
       this%vay(i)=this%vay(i)+0.5d0*this%fay(i)*this%get_timestep()/this%mass(i)
       this%vaz(i)=this%vaz(i)+0.5d0*this%faz(i)*this%get_timestep()/this%mass(i)
    end do
    call this%set_ekinetic() !thermostat
    call this%set_tfcnpt()
    do i=1,this%get_natom()
       this%vax(i)=this%vax(i)*exp(-0.25d0*this%get_tfcnpt()*this%get_timestep())
       this%vay(i)=this%vay(i)*exp(-0.25d0*this%get_tfcnpt()*this%get_timestep())
       this%vaz(i)=this%vaz(i)*exp(-0.25d0*this%get_tfcnpt()*this%get_timestep())
    end do
    call this%set_ekinetic()
    call this%set_tfcnpt()
    call this%set_pressure() !barostat
    call this%set_bfc2()
    call this%set_bfc()
    call this%set_bfc2()
    do i=1,this%get_natom()
       this%vax(i)=this%vax(i)*exp(-0.5d0*this%get_bfc()*this%get_timestep())
       this%vay(i)=this%vay(i)*exp(-0.5d0*this%get_bfc()*this%get_timestep())
       this%vaz(i)=this%vaz(i)*exp(-0.5d0*this%get_bfc()*this%get_timestep())
    end do
    call this%set_bfc2()
    call this%set_bfc()
    call this%set_bfc2()
    call this%set_ekinetic() !thermostat
    call this%set_tfcnpt()
    do i=1,this%get_natom()
       this%vax(i)=this%vax(i)*exp(-0.25d0*this%get_tfcnpt()*this%get_timestep())
       this%vay(i)=this%vay(i)*exp(-0.25d0*this%get_tfcnpt()*this%get_timestep())
       this%vaz(i)=this%vaz(i)*exp(-0.25d0*this%get_tfcnpt()*this%get_timestep())
    end do
    call this%set_ekinetic()
    call this%set_tfcnpt()
    call this%set_vcm()
    do i=1,this%get_natom()
       this%vax(i)=this%vax(i)-this%vcm(1)
       this%vay(i)=this%vay(i)-this%vcm(2)
       this%vaz(i)=this%vaz(i)-this%vcm(3)
    end do
    call this%set_ekinetic()
    call this%set_etotal()
    call this%check_energy()
    call this%set_temperature()
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

  subroutine check_energy(this)
    implicit none
    class(ensemble), intent(inout) :: this
    if(this%get_etotal()/this%get_natom().ge.this%get_checkenergy())then
       write(6,*)
       write(6,'(a42,es7.1,a61)')'ERROR: The energy per atom is higher than ',&
            this%get_checkenergy()*this%get_econv(),&
            ' kcal/mol! The simulation was interrupted to avoid explosion.'
       stop
    end if
  end subroutine check_energy

  subroutine set_rcm(this)
    implicit none
    class(ensemble), intent(inout) :: this
    integer                        :: i
    real(8)                        :: rcm(3)
    rcm(1)=0.d0
    rcm(2)=0.d0
    rcm(3)=0.d0
    do i=1,this%get_natom()
       rcm(1)=rcm(1)+this%mass(i)*this%xa(i)
       rcm(2)=rcm(2)+this%mass(i)*this%ya(i)
       rcm(3)=rcm(3)+this%mass(i)*this%za(i)
    end do
    this%rcm(1)=rcm(1)/this%sys%get_mtotal()
    this%rcm(2)=rcm(2)/this%sys%get_mtotal()
    this%rcm(3)=rcm(3)/this%sys%get_mtotal()
  end subroutine set_rcm

  subroutine set_vcm(this)
    implicit none
    class(ensemble), intent(inout) :: this
    integer                        :: i
    real(8)                        :: vcm(3)
    vcm(1)=0.d0
    vcm(2)=0.d0
    vcm(3)=0.d0
    do i=1,this%get_natom()
       vcm(1)=vcm(1)+this%mass(i)*this%vax(i)
       vcm(2)=vcm(2)+this%mass(i)*this%vay(i)
       vcm(3)=vcm(3)+this%mass(i)*this%vaz(i)
    end do
    this%vcm(1)=vcm(1)/this%sys%get_mtotal()
    this%vcm(2)=vcm(2)/this%sys%get_mtotal()
    this%vcm(3)=vcm(3)/this%sys%get_mtotal()
  end subroutine set_vcm

end module ensemble_module
