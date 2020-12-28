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
module elements_module
  !*******************************************************************************************
  !*******************************************************************************************

  use system_module

  private
  public :: elements

  type, extends(system) :: elements
     integer, allocatable  :: zat(:)
     real(8), allocatable  :: mass(:)
   contains
     procedure :: covalent_radius
     procedure :: set_mass
  end type elements

contains

  subroutine set_zat(this)
    implicit none
    class(elements), intent(inout) :: this
    integer                        :: i
    allocate(this%zat(this%get_natom()))
    do i=1,this%get_natom()
       select case(this%at(i))
       case('H ')
          zat=1
       case('He')
          zat=2
       case('Li')
          zat=3
       case('Be')
          zat=4
       case('B ')
          zat=5
       case('C ')
          zat=6
       case('N ')
          zat=7
       case('O ')
          zat=8
       case('F ')
          zat=9
       case('Ne')
          zat=10
       case('Na')
          zat=11
       case('Mg')
          zat=12
       case('Al')
          zat=13
       case('Si')
          zat=14
       case('P ')
          zat=15
       case('S ')
          zat=16
       case('Cl')
          zat=17
       end select
    end do
  end subroutine set_zat

  subroutine covalent_radius(this,i,rc)
    implicit none
    class(elements), intent(inout) :: this
    integer, intent(in)            :: i
    real(8), intent(out)           :: rc
    select case(this%zat(i))
    case(1)
       rc=0.37d0
    case(2)
       rc=0.32d0
    case(3)
       rc=1.34d0
    case(4)
       rc=0.90d0
    case(5)
       rc=0.82d0
    case(6)
       rc=0.77d0
    case(7)
       rc=0.75d0
    case(8)
       rc=0.73d0
    case(9)
       rc=0.71d0
    case(10)
       rc=0.69d0
    case(11)
       rc=1.54d0
    case(12)
       rc=1.30d0
    case(13)
       rc=1.18d0
    case(14)
       rc=1.11d0
    case(15)
       rc=1.06d0
    case(16)
       rc=1.02d0
    case(17)
       rc=0.99d0
    case(18)
       rc=0.97d0
    case(19)
       rc=1.96d0
    case(20)
       rc=1.74d0
    case(21)
       rc=1.44d0
    case(22)
       rc=1.36d0
    case(23)
       rc=1.25d0
    case(24)
       rc=1.27d0
    case(25)
       rc=1.39d0
    case(26)
       rc=1.25d0
    case(27)
       rc=1.26d0
    case(28)
       rc=1.21d0
    case(29)
       rc=1.38d0
    case(30)
       rc=1.31d0
    case(31)
       rc=1.26d0
    case(32)
       rc=1.22d0
    case(33)
       rc=1.19d0
    case(34)
       rc=1.16d0
    case(35)
       rc=1.14d0
    case(36)
       rc=1.10d0
    case(37)
       rc=2.11d0
    case(38)
       rc=1.92d0
    case(39)
       rc=1.62d0
    case(40)
       rc=1.48d0
    case(41)
       rc=1.37d0
    case(42)
       rc=1.45d0
    case(43)
       rc=1.56d0
    case(44)
       rc=1.26d0
    case(45)
       rc=1.35d0
    case(46)
       rc=1.31d0
    case(47)
       rc=1.53d0
    case(48)
       rc=1.48d0
    case(49)
       rc=1.44d0
    case(50)
       rc=1.41d0
    case(51)
       rc=1.38d0
    case(52)
       rc=1.35d0
    case(53)
       rc=1.33d0
    case(54)
       rc=1.30d0
    case(55)
       rc=2.25d0
    case(56)
       rc=1.98d0
    case(57)
       rc=1.69d0
    case(71)
       rc=1.60d0
    case(72)
       rc=1.50d0
    case(73)
       rc=1.38d0
    case(74)
       rc=1.46d0
    case(75)
       rc=1.59d0
    case(76)
       rc=1.28d0
    case(77)
       rc=1.37d0
    case(78)
       rc=1.28d0
    case(79)
       rc=1.44d0
    case(80)
       rc=1.49d0
    case(81)
       rc=1.48d0
    case(82)
       rc=1.47d0
    case(83)
       rc=1.46d0
    end select
  end subroutine covalent_radius

  subroutine set_mass(this,nx)
    implicit none
    class(elements), intent (inout) :: this
    integer, intent(in)             :: nx
    integer                         :: i
    allocate(this%mass(nx))
    do i=1,nx
       select case(this%zat(i))
       case(1)
          this%mass(i)=1.007940000d0
       case(2)
          this%mass(i)=4.002602000d0
       case(3)
          this%mass(i)=6.941000000d0
       case(4)
          this%mass(i)=9.012182000d0
       case(5)
          this%mass(i)=10.81100000d0
       case(6)
          this%mass(i)=12.01070000d0
       case(7)
          this%mass(i)=14.00670000d0
       case(8)
          this%mass(i)=15.99940000d0
       case(9)
          this%mass(i)=18.99840300d0
       case(10)
          this%mass(i)=20.17970000d0
       case(11)
          this%mass(i)=22.98976928d0
       case(12)
          this%mass(i)=24.30500000d0
       case(13)
          this%mass(i)=26.98153860d0
       case(14)
          this%mass(i)=28.08550000d0
       case(15)
          this%mass(i)=30.97376200d0
       case(16)
          this%mass(i)=32.06500000d0
       case(17)
          this%mass(i)=35.45300000d0
       case(18)
          this%mass(i)=39.94800000d0
       case(19)
          this%mass(i)=39.09830000d0
       case(20)
          this%mass(i)=40.07800000d0
       case(21)
          this%mass(i)=44.95591200d0
       case(22)
          this%mass(i)=47.86700000d0
       case(23)
          this%mass(i)=50.94150000d0
       case(24)
          this%mass(i)=51.99610000d0
       case(25)
          this%mass(i)=54.93804500d0
       case(26)
          this%mass(i)=55.84500000d0
       case(27)
          this%mass(i)=58.93319500d0
       case(28)
          this%mass(i)=58.69340000d0
       case(29)
          this%mass(i)=58.69340000d0
       case(30)
          this%mass(i)=65.40900000d0
       case(31)
          this%mass(i)=69.72300000d0
       case(32)
          this%mass(i)=72.64000000d0
       case(33)
          this%mass(i)=74.92160000d0
       case(34)
          this%mass(i)=78.96000000d0
       case(35)
          this%mass(i)=79.90400000d0
       case(36)
          this%mass(i)=83.79800000d0
       case(37)
          this%mass(i)=85.46780000d0
       case(38)
          this%mass(i)=87.62000000d0
       case(39)
          this%mass(i)=88.90585000d0
       case(40)
          this%mass(i)=91.22400000d0
       case(41)
          this%mass(i)=92.90638000d0
       case(42)
          this%mass(i)=95.94000000d0
       case(43)
          this%mass(i)=98.00000000d0
       case(44)
          this%mass(i)=101.07000000d0
       case(45)
          this%mass(i)=102.90550000d0
       case(46)
          this%mass(i)=106.42000000d0
       case(47)
          this%mass(i)=107.86820000d0
       case(48)
          this%mass(i)=112.41100000d0
       case(49)
          this%mass(i)=114.81800000d0
       case(50)
          this%mass(i)=118.71000000d0
       case(51)
          this%mass(i)=121.76000000d0
       case(52)
          this%mass(i)=127.60000000d0
       case(53)
          this%mass(i)=126.90447000d0
       case(54)
          this%mass(i)=131.29300000d0
       case(55)
          this%mass(i)=132.90545190d0
       case(56)
          this%mass(i)=137.32700000d0
       case(57)
          this%mass(i)=138.90547000d0
       case(58)
          this%mass(i)=140.11600000d0
       case(59)
          this%mass(i)=140.90765000d0
       case(60)
          this%mass(i)=144.24200000d0
       case(61)
          this%mass(i)=145.00000000d0
       case(62)
          this%mass(i)=150.36000000d0
       case(63)
          this%mass(i)=151.96400000d0
       case(64)
          this%mass(i)=157.25000000d0
       case(65)
          this%mass(i)=158.92535000d0
       case(66)
          this%mass(i)=162.50000000d0
       case(67)
          this%mass(i)=164.93032000d0
       case(68)
          this%mass(i)=167.25900000d0
       case(69)
          this%mass(i)=168.93421000d0
       case(70)
          this%mass(i)=173.04000000d0
       case(71)
          this%mass(i)=174.96700000d0
       case(72)
          this%mass(i)=178.49000000d0
       case(73)
          this%mass(i)=180.94788000d0
       case(74)
          this%mass(i)=183.84000000d0
       case(75)
          this%mass(i)=186.20700000d0
       case(76)
          this%mass(i)=190.23000000d0
       case(77)
          this%mass(i)=192.21700000d0
       case(78)
          this%mass(i)=195.08400000d0
       case(79)
          this%mass(i)=196.96656900d0
       case(80)
          this%mass(i)=200.59000000d0
       case(81)
          this%mass(i)=204.38330000d0
       case(82)
          this%mass(i)=207.20000000d0
       case(83)
          this%mass(i)=208.98040000d0
       case(84)
          this%mass(i)=209.00000000d0
       case(85)
          this%mass(i)=210.00000000d0
       case(86)
          this%mass(i)=222.00000000d0
       case(87)
          this%mass(i)=223.00000000d0
       case(88)
          this%mass(i)=226.00000000d0
       case(89)
          this%mass(i)=227.00000000d0
       case(90)
          this%mass(i)=232.03806000d0
       case(91)
          this%mass(i)=231.03588000d0
       case(92)
          this%mass(i)=238.02891000d0
       case(93)
          this%mass(i)=237.00000000d0
       case(94)
          this%mass(i)=244.00000000d0
       case(95)
          this%mass(i)=243.00000000d0
       case(96)
          this%mass(i)=247.00000000d0
       case(97)
          this%mass(i)=247.00000000d0
       case(98)
          this%mass(i)=251.00000000d0
       case(99)
          this%mass(i)=251.00000000d0
       case(100)
          this%mass(i)=257.00000000d0
       case(101)
          this%mass(i)=258.00000000d0
       case(102)
          this%mass(i)=259.00000000d0
       case(103)
          this%mass(i)=262.00000000d0
       case(104)
          this%mass(i)=261.00000000d0
       case(105)
          this%mass(i)=262.00000000d0
       case(106)
          this%mass(i)=266.00000000d0
       case(107)
          this%mass(i)=264.00000000d0
       case(108)
          this%mass(i)=277.00000000d0
       case(109)
          this%mass(i)=268.00000000d0
       case(110)
          this%mass(i)=281.00000000d0
       case(111)
          this%mass(i)=272.00000000d0
       case(112)
          this%mass(i)=285.00000000d0
       case(113)
          this%mass(i)=284.00000000d0
       case(114)
          this%mass(i)=289.00000000d0
       case(115)
          this%mass(i)=288.00000000d0
       case(116)
          this%mass(i)=292.00000000d0
       case(118)
          this%mass(i)=294.00000000d0
       end select
    end do
  end subroutine set_mass

end module elements_module
