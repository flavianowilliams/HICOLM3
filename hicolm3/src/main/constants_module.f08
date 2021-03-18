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
module constants_module
  !*******************************************************************************************
  !*******************************************************************************************

  implicit none

  private
  public :: constants

  type :: constants
     real(8), private :: pi
     real(8), private :: kb
     real(8), private :: n0
     real(8), private :: kelect
     real(8), private :: mconv
     real(8), private :: aconv
     real(8), private :: hcconv
     real(8), private :: elconv
     real(8), private :: keconv
     real(8), private :: econv
     real(8), private :: rconv
     real(8), private :: pconv
     real(8), private :: tconv
     real(8), private :: teconv
     real(8), private :: kconv
   contains
     procedure :: constants_prepare
     procedure :: get_pi
     procedure :: get_kb
     procedure :: get_n0
     procedure :: get_kelect
     procedure :: get_mconv
     procedure :: get_aconv
     procedure :: get_hcconv
     procedure :: get_elconv
     procedure :: get_keconv
     procedure :: get_econv
     procedure :: get_rconv
     procedure :: get_pconv
     procedure :: get_tconv
     procedure :: get_teconv
     procedure :: get_kconv
  end type constants

contains

  subroutine constants_prepare(this)
    implicit none
    class(constants), intent(inout) :: this

    !-constantes físicas

    this%pi=dble(acos(-1.d0))
    this%kelect=1.d0
    this%n0=6.022d+23
    this%kb=8.6173324d-5

    !-unidades fundamentais (SI)

    this%mconv=9.10938291d-31                      ! massa (kg)
    this%aconv=180.d0/acos(-1.d0)                  ! ângulo (radianos)
    this%hcconv=1.054571726d-34                    ! Constante de Planck (J*s).
    this%elconv=1.602176565d-19                    ! carga eletrica (C)
    this%keconv=8.9875517873681d+9                 ! constante eletrostática (N*m^2/C^2)

    !-unidades derivadas (SI)

    this%econv=this%mconv*this%elconv**4*this%keconv**2/this%hcconv**2 ! energia (J)
    this%rconv=this%hcconv**2/(this%mconv*this%keconv*this%elconv**2)  ! comprimento (m)
    this%pconv=this%econv/this%rconv**3                                ! pressao (N/m^2)
    this%tconv=this%hcconv/this%econv                                  ! tempo (s)
    this%teconv=this%econv/this%kb                                     ! temperatura (K)

    !-convertendo unidades SI para unidades de entrada

    this%mconv=this%mconv*(1.d+3*this%n0)                              !
    this%rconv=this%rconv*(1.d+10)                                     !
    this%kconv=1.d0/this%rconv                                         !
    this%econv=this%econv*(1.438689d+20)                               ! Joule --> kcal/mol
    this%pconv=this%pconv*(0.9872d-5)                                  !
    this%tconv=this%tconv*(1.d+12)                                     ! second --> picosecond
    this%elconv=1.d0                                                   !
    this%teconv=this%teconv*(6.241506363094d+18)                       !

  end subroutine constants_prepare

  double precision function get_pi(this)
    class(constants), intent(in) :: this
    get_pi=this%pi
  end function get_pi

  double precision function get_kb(this)
    class(constants), intent(in) :: this
    get_kb=this%kb
  end function get_kb

  double precision function get_n0(this)
    class(constants), intent(in) :: this
    get_n0=this%n0
  end function get_n0

  double precision function get_kelect(this)
    class(constants), intent(in) :: this
    get_kelect=this%kelect
  end function get_kelect

  double precision function get_mconv(this)
    class(constants), intent(in) :: this
    get_mconv=this%mconv
  end function get_mconv

  double precision function get_aconv(this)
    class(constants), intent(in) :: this
    get_aconv=this%aconv
  end function get_aconv

  double precision function get_hcconv(this)
    class(constants), intent(in) :: this
    get_hcconv=this%hcconv
  end function get_hcconv

  double precision function get_elconv(this)
    class(constants), intent(in) :: this
    get_elconv=this%elconv
  end function get_elconv

  double precision function get_keconv(this)
    class(constants), intent(in) :: this
    get_keconv=this%keconv
  end function get_keconv

  double precision function get_econv(this)
    class(constants), intent(in) :: this
    get_econv=this%econv
  end function get_econv

  double precision function get_rconv(this)
    class(constants), intent(in) :: this
    get_rconv=this%rconv
  end function get_rconv

  double precision function get_pconv(this)
    class(constants), intent(in) :: this
    get_pconv=this%pconv
  end function get_pconv

  double precision function get_tconv(this)
    class(constants), intent(in) :: this
    get_tconv=this%tconv
  end function get_tconv

  double precision function get_teconv(this)
    class(constants), intent(in) :: this
    get_teconv=this%teconv
  end function get_teconv

  double precision function get_kconv(this)
    class(constants), intent(in) :: this
    get_kconv=this%kconv
  end function get_kconv

end module constants_module
