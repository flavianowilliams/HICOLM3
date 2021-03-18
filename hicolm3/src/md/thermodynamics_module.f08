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
module thermodynamics_module
  !*******************************************************************************************
  !*******************************************************************************************

  use interaction_module

  implicit none

  private
  public :: thermodynamics

  type, extends(interaction) :: thermodynamics
     real(8)   :: temperature
     real(8)   :: pressure
     real(8)   :: ekinetic
     real(8)   :: etotal
   contains
     procedure :: set_temperature
     procedure :: get_temperature
     procedure :: set_pressure
     procedure :: get_pressure
     procedure :: set_ekinetic
     procedure :: get_ekinetic
     procedure :: set_etotal
     procedure :: get_etotal
  end type thermodynamics

contains

  subroutine set_temperature(this)
    implicit none
    class(thermodynamics), intent(inout) :: this
    this%temperature=2.d0*this%ekinetic/this%get_nfree()
  end subroutine set_temperature

  double precision function get_temperature(this)
    implicit none
    class(thermodynamics), intent(inout) :: this
    get_temperature=this%temperature
  end function get_temperature

  subroutine set_pressure(this)
    implicit none
    class(thermodynamics), intent(inout) :: this
    this%pressure=(2.d0*this%get_ekinetic()-this%get_virtot())/(3.d0*this%get_volume())
  end subroutine set_pressure

  double precision function get_pressure(this)
    implicit none
    class(thermodynamics), intent(inout) :: this
    get_pressure = this%pressure
  end function get_pressure

  subroutine set_ekinetic(this)
    implicit none
    class(thermodynamics), intent(inout) :: this
    integer                        :: i
    real(8)                        :: ekinetic
    ekinetic=0.d0
    do i=1,this%get_natom()
       ekinetic=ekinetic+this%mass(i)*(this%vax(i)**2+this%vay(i)**2+this%vaz(i)**2)
    end do
    this%ekinetic=0.5d0*ekinetic
  end subroutine set_ekinetic

  double precision function get_ekinetic(this)
    implicit none
    class(thermodynamics), intent(inout) :: this
    get_ekinetic = this%ekinetic
  end function get_ekinetic

  subroutine set_etotal(this)
    implicit none
    class(thermodynamics), intent(inout) :: this
    this%etotal = this%get_ekinetic()+this%get_enpot()
  end subroutine set_etotal

  double precision function get_etotal(this)
    implicit none
    class(thermodynamics), intent(inout) :: this
    get_etotal = this%etotal
  end function get_etotal

end module thermodynamics_module
