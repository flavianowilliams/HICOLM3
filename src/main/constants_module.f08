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

  private
  public :: constants

  type :: constants
     real(8) :: pi
     real(8) :: kb
     real(8) :: n0
     real(8) :: kelect
     real(8) :: mconv
     real(8) :: aconv
     real(8) :: hcconv
     real(8) :: elconv
     real(8) :: keconv
     real(8) :: econv
     real(8) :: rconv
     real(8) :: pconv
     real(8) :: tconv
     real(8) :: teconv
     real(8) :: kconv
   contains
     procedure :: constants_prepare
  end type constants

  interface constants
     module procedure constructor
  end interface constants

contains

  type(constants) function constructor()
    implicit none
    call constructor%constants_prepare
  end function constructor

  subroutine constants_prepare(this)
    implicit none
    class(constants), intent(inout) :: this

    !-constantes físicas

    this%pi=dble(acos(-1.d0))
    this%kb=8.6173324d-5
    this%kelect=1.d0
    this%n0=6.022d+23

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

end module constants_module
