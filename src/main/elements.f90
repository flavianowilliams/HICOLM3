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
module elements

  use sistema

  real(8) mass(natmax)

  save mass

contains

  subroutine covalent_radius(zat,rc)

    implicit none

    integer zat
    real(8) rc

    select case(zat)
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

    return

  end subroutine covalent_radius

  subroutine atomic_mass(i,zat)
    !***************************************************************************************
    ! Massas atomicas                                                                      *
    !***************************************************************************************
    implicit none

    integer zat,i

    !-definindo massa atomica

    select case(zat)
    case(1)
       mass(i)=1.007940000d0
    case(2)
       mass(i)=4.002602000d0
    case(3)
       mass(i)=6.941000000d0
    case(4)
       mass(i)=9.012182000d0
    case(5)
       mass(i)=10.81100000d0
    case(6)
       mass(i)=12.01070000d0
    case(7)
       mass(i)=14.00670000d0
    case(8)
       mass(i)=15.99940000d0
    case(9)
       mass(i)=18.99840300d0
    case(10)
       mass(i)=20.17970000d0
    case(11)
       mass(i)=22.98976928d0
    case(12)
       mass(i)=24.30500000d0
    case(13)
       mass(i)=26.98153860d0
    case(14)
       mass(i)=28.08550000d0
    case(15)
       mass(i)=30.97376200d0
    case(16)
       mass(i)=32.06500000d0
    case(17)
       mass(i)=35.45300000d0
    case(18)
       mass(i)=39.94800000d0
    case(19)
       mass(i)=39.09830000d0
    case(20)
       mass(i)=40.07800000d0
    case(21)
       mass(i)=44.95591200d0
    case(22)
       mass(i)=47.86700000d0
    case(23)
       mass(i)=50.94150000d0
    case(24)
       mass(i)=51.99610000d0
    case(25)
       mass(i)=54.93804500d0
    case(26)
       mass(i)=55.84500000d0
    case(27)
       mass(i)=58.93319500d0
    case(28)
       mass(i)=58.69340000d0
    case(29)
       mass(i)=58.69340000d0
    case(30)
       mass(i)=65.40900000d0
    case(31)
       mass(i)=69.72300000d0
    case(32)
       mass(i)=72.64000000d0
    case(33)
       mass(i)=74.92160000d0
    case(34)
       mass(i)=78.96000000d0
    case(35)
       mass(i)=79.90400000d0
    case(36)
       mass(i)=83.79800000d0
    case(37)
       mass(i)=85.46780000d0
    case(38)
       mass(i)=87.62000000d0
    case(39)
       mass(i)=88.90585000d0
    case(40)
       mass(i)=91.22400000d0
    case(41)
       mass(i)=92.90638000d0
    case(42)
       mass(i)=95.94000000d0
    case(43)
       mass(i)=98.00000000d0
    case(44)
       mass(i)=101.07000000d0
    case(45)
       mass(i)=102.90550000d0
    case(46)
       mass(i)=106.42000000d0
    case(47)
       mass(i)=107.86820000d0
    case(48)
       mass(i)=112.41100000d0
    case(49)
       mass(i)=114.81800000d0
    case(50)
       mass(i)=118.71000000d0
    case(51)
       mass(i)=121.76000000d0
    case(52)
       mass(i)=127.60000000d0
    case(53)
       mass(i)=126.90447000d0
    case(54)
       mass(i)=131.29300000d0
    case(55)
       mass(i)=132.90545190d0
    case(56)
       mass(i)=137.32700000d0
    case(57)
       mass(i)=138.90547000d0
    case(58)
       mass(i)=140.11600000d0
    case(59)
       mass(i)=140.90765000d0
    case(60)
       mass(i)=144.24200000d0
    case(61)
       mass(i)=145.00000000d0
    case(62)
       mass(i)=150.36000000d0
    case(63)
       mass(i)=151.96400000d0
    case(64)
       mass(i)=157.25000000d0
    case(65)
       mass(i)=158.92535000d0
    case(66)
       mass(i)=162.50000000d0
    case(67)
       mass(i)=164.93032000d0
    case(68)
       mass(i)=167.25900000d0
    case(69)
       mass(i)=168.93421000d0
    case(70)
       mass(i)=173.04000000d0
    case(71)
       mass(i)=174.96700000d0
    case(72)
       mass(i)=178.49000000d0
    case(73)
       mass(i)=180.94788000d0
    case(74)
       mass(i)=183.84000000d0
    case(75)
       mass(i)=186.20700000d0
    case(76)
       mass(i)=190.23000000d0
    case(77)
       mass(i)=192.21700000d0
    case(78)
       mass(i)=195.08400000d0
    case(79)
       mass(i)=196.96656900d0
    case(80)
       mass(i)=200.59000000d0
    case(81)
       mass(i)=204.38330000d0
    case(82)
       mass(i)=207.20000000d0
    case(83)
       mass(i)=208.98040000d0
    case(84)
       mass(i)=209.00000000d0
    case(85)
       mass(i)=210.00000000d0
    case(86)
       mass(i)=222.00000000d0
    case(87)
       mass(i)=223.00000000d0
    case(88)
       mass(i)=226.00000000d0
    case(89)
       mass(i)=227.00000000d0
    case(90)
       mass(i)=232.03806000d0
    case(91)
       mass(i)=231.03588000d0
    case(92)
       mass(i)=238.02891000d0
    case(93)
       mass(i)=237.00000000d0
    case(94)
       mass(i)=244.00000000d0
    case(95)
       mass(i)=243.00000000d0
    case(96)
       mass(i)=247.00000000d0
    case(97)
       mass(i)=247.00000000d0
    case(98)
       mass(i)=251.00000000d0
    case(99)
       mass(i)=251.00000000d0
    case(100)
       mass(i)=257.00000000d0
    case(101)
       mass(i)=258.00000000d0
    case(102)
       mass(i)=259.00000000d0
    case(103)
       mass(i)=262.00000000d0
    case(104)
       mass(i)=261.00000000d0
    case(105)
       mass(i)=262.00000000d0
    case(106)
       mass(i)=266.00000000d0
    case(107)
       mass(i)=264.00000000d0
    case(108)
       mass(i)=277.00000000d0
    case(109)
       mass(i)=268.00000000d0
    case(110)
       mass(i)=281.00000000d0
    case(111)
       mass(i)=272.00000000d0
    case(112)
       mass(i)=285.00000000d0
    case(113)
       mass(i)=284.00000000d0
    case(114)
       mass(i)=289.00000000d0
    case(115)
       mass(i)=288.00000000d0
    case(116)
       mass(i)=292.00000000d0
    case(118)
       mass(i)=294.00000000d0
    end select

    return

  end subroutine atomic_mass

end module elements
